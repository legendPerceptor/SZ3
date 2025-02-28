//
// Created by yuanjian on 4/1/24.
//

#ifndef SZ3_COMPRESSIONMPIMANAGER_H
#define SZ3_COMPRESSIONMPIMANAGER_H

#include "DebugStream.hpp"
#include "SZ3/api/sz.hpp"
#include "split_common.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

namespace sz3_split {

template <typename T> class CompressionMPIManager {

    std::string input_file;
    std::string output_file;
    std::vector<size_t> dimension;
    size_t num_io_processes;
    T eb;
    size_t depth;
    bool is_compression_mode;
    size_t total_ranks;

  public:
    CompressionMPIManager(std::string input_file, std::string output_file,
                          std::vector<size_t> dimension, T eb, size_t depth, bool is_compression,
                          size_t num_io_processes)
        : input_file(std::move(input_file)), output_file(std::move(output_file)),
          dimension(std::move(dimension)), eb(eb), depth(depth),
          is_compression_mode(is_compression), num_io_processes(num_io_processes) {}

    void startMPI() {
        int is_mpi_initialized = 0;
        MPI_Initialized(&is_mpi_initialized);
        if (!is_mpi_initialized) {
            MPI_Init(NULL, NULL);
        }
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        double start_time, end_time;
        start_time = MPI_Wtime();
        this->total_ranks = size;
        if (size < 3) {
            std::cerr << "Error: At least 3 processes are needed!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (num_io_processes < 2) {
            std::cerr << "Error: At least 2 I/O processes are needed!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);
        printf("Process %d of %d on node %s\n", rank, size, processor_name);

        if (is_compression_mode) {
            if (rank >= 0 && rank < num_io_processes - 1) {
                debugStream << "read rank: " << rank << ", size:" << size << std::endl;
                readProcess(rank, num_io_processes - 1);
            } else if (rank == num_io_processes - 1) {
                debugStream << "write rank: " << rank << ", size:" << size << std::endl;
                writerProcess(rank);
            } else {
                debugStream << "worker rank: " << rank << ", size:" << size << std::endl;
                workerProcess(rank);
            }
        } else {
            if (rank == 0) {
                readProcess(rank, 1);
            } else if (rank >= 1 && rank < num_io_processes) {
                writerProcess(rank);
            } else {
                workerProcess(rank);
            }
        }
        end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        debugStream << "[stats] rank<" << rank << "> total run time: " << elapsed_time << "seconds"
                    << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }

  private:
    void readProcess(size_t rank, size_t number_of_read_processes) {
        assert(dimension.size() == 3);
        SZ3::Config conf;
        if (depth == 1) {
            conf.setDims(dimension.begin(), dimension.end() - 1);
        } else {
            std::vector<size_t> chunk_dimension = {depth, dimension[1], dimension[0]};
            conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
        }
        size_t chunk_size = sizeof(T) * conf.num;
        conf.absErrorBound = eb;

        MPI_File fh;
        MPI_Offset file_size;
        int status;

        status =
            MPI_File_open(MPI_COMM_SELF, input_file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_get_size(fh, &file_size);

        MPI_Offset dp_offset = 0;

        double total_read_time = 0;
        size_t num_iterations = dimension[2] / depth;
        size_t leftover = dimension[2] % depth;
        if (leftover > 0) {
            num_iterations += 1;
        }
        SZ3::Timer temp(false);
        size_t org_chunk_size = chunk_size;
        debugStream << "read process with rank " << rank << " start working!" << std::endl;
        std::vector<size_t> managed_worker;
        for (size_t worker_rank = num_io_processes; worker_rank < total_ranks; worker_rank++) {
            if (map_worker_to_reader(worker_rank) == rank) {
                managed_worker.push_back(worker_rank);
            }
        }
        size_t next_worker_index = 0;
        for (size_t i = rank; i < num_iterations; i += number_of_read_processes) {
            if (i == num_iterations - 1 && leftover > 0) {
                std::vector<size_t> chunk_dimension = {leftover, dimension[1], dimension[0]};
                conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
                chunk_size = sizeof(T) * conf.num;
            }
            DataChunk<T> chunk;
            chunk.id = i;
            chunk.sequenceNumber = i;
            chunk.conf = conf;

            if (is_compression_mode) {
                chunk.dataBuffer = std::vector<T>(conf.num);
                // Calculate start position for this thread, the last chunk may have different
                // chunk_size
                temp.start();
                // Calculate the offset for reading
                MPI_Offset offset = i * org_chunk_size;
                // Set the file view for each process
                MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
                MPI_File_read_all(fh, chunk.dataBuffer.data(), conf.num * sizeof(T), MPI_BYTE,
                                  MPI_STATUS_IGNORE);
                // Read data using MPI I/O
                debugStream << "[read] first 10 values in chunk [" << chunk.sequenceNumber << "]: ";
                for (int _index_test = 0; _index_test < 10; _index_test++) {
                    debugStream << chunk.dataBuffer[_index_test] << " ";
                }
                debugStream << std::endl;
                debugStream << "[read] last 10 values in chunk [" << chunk.sequenceNumber
                            << "] with offset " << offset << ": ";
                for (int _index_test = chunk.dataBuffer.size() - 10;
                     _index_test < chunk.dataBuffer.size(); _index_test++) {
                    debugStream << chunk.dataBuffer[_index_test] << " ";
                }
                debugStream << std::endl;
                total_read_time += temp.stop();
            } else {
                temp.start();
                int64_t compressed_chunk_size;
                MPI_File_set_view(fh, dp_offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
                MPI_File_read_all(fh, &compressed_chunk_size, sizeof(int64_t), MPI_BYTE,
                                  MPI_STATUS_IGNORE);
                chunk.cpdataBuffer = std::vector<char>(compressed_chunk_size);
                MPI_File_read_all(fh, chunk.cpdataBuffer.data(), compressed_chunk_size, MPI_BYTE,
                                  MPI_STATUS_IGNORE);
                dp_offset += sizeof(int64_t) + compressed_chunk_size;
                total_read_time += temp.stop();
            }

            size_t dest = managed_worker[next_worker_index];
            next_worker_index = (next_worker_index + 1) % managed_worker.size();

            debugStream << "read process with rank " << rank << " sending chunk "
                        << chunk.sequenceNumber << " to worker " << dest << std::endl;
            if (i == num_iterations - 1 && leftover > 0) {
                sendDataToWorker(dest, chunk, rank, leftover);
            } else {
                sendDataToWorker(dest, chunk, rank, depth);
            }
        }

        debugStream << "[read stats]read process with rank " << rank
                    << " finished. Total read time:" << total_read_time << std::endl;

        for (auto dest : managed_worker) {
            debugStream << "[signal worker] read process with rank " << rank
                        << " signaling worker rank " << dest << " to exit!" << std::endl;
            signalReadCompleteToWorkers(rank, dest);
        }

        MPI_File_close(&fh);
        debugStream << "[exit]read process with rank " << rank << "exited!" << std::endl;
    }

    void sendDataToWorker(int dest, DataChunk<T>& chunk, size_t read_rank, size_t layer_depth) {
        int tag = 0;
        // Serialize chunk
        if (is_compression_mode) {
            std::vector<char> buffer(chunk.dataBuffer.size() * sizeof(T) + 3 * sizeof(size_t));
            char* pointer = buffer.data();
            size_t chunk_buffer_size = chunk.dataBuffer.size();
            memcpy(pointer, &chunk_buffer_size, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, &chunk.sequenceNumber, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, &layer_depth, sizeof(size_t));
            pointer += sizeof(size_t);

            memcpy(pointer, chunk.dataBuffer.data(), chunk.dataBuffer.size() * sizeof(T));
            //                debugStream << "read rank: " << read_rank << ", about to MPI_Send,
            //                buffer size:" << buffer.size() << std::endl;
            // Send chunk to worker
            MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
            //                debugStream << "read rank: " << read_rank << ", finished MPI_Send,
            //                buffer size:" << buffer.size() << std::endl;
        } else {
            std::vector<char> buffer(chunk.cpdataBuffer.size() + 3 * sizeof(size_t));
            char* pointer = buffer.data();
            size_t chunk_buffer_size = chunk.cpdataBuffer.size();
            memcpy(pointer, &chunk_buffer_size, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, &chunk.sequenceNumber, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, &layer_depth, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, chunk.cpdataBuffer.data(), chunk.cpdataBuffer.size());
            // Send chunk to worker
            MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
        }
    }

    void sendDataToWriter(int dest, DataChunk<T>& chunk, size_t rank) {
        int tag = 1;
        if (is_compression_mode) {
            std::vector<char> buffer(chunk.cpdataBuffer.size() + 2 * sizeof(size_t));
            char* pointer = buffer.data();
            size_t chunk_buffer_size = chunk.cpdataBuffer.size();
            memcpy(pointer, &chunk_buffer_size, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, &chunk.sequenceNumber, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, chunk.cpdataBuffer.data(), chunk_buffer_size);
            MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
        } else {
            std::vector<char> buffer(chunk.dataBuffer.size() * sizeof(T) + 2 * sizeof(size_t));
            char* pointer = buffer.data();
            size_t chunk_buffer_size = chunk.dataBuffer.size();
            memcpy(pointer, &chunk_buffer_size, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, &chunk.sequenceNumber, sizeof(size_t));
            pointer += sizeof(size_t);
            memcpy(pointer, chunk.dataBuffer.data(), chunk.dataBuffer.size() * sizeof(T));
            MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
        }
    }

    void signalReadCompleteToWorkers(size_t read_rank, size_t dest) {
        size_t buffer_size = 0;
        MPI_Send(&buffer_size, sizeof(size_t), MPI_BYTE, dest, 0, MPI_COMM_WORLD);
    }

    void signalWorkerCompleteToWriter(size_t worker_rank, size_t dest) {
        size_t buffer_size = 0;
        MPI_Send(&buffer_size, sizeof(size_t), MPI_BYTE, dest, 1, MPI_COMM_WORLD);
    }

    int map_worker_to_reader(size_t rank) {
        if (is_compression_mode) {
            return (rank - num_io_processes) % (num_io_processes - 1);
        } else {
            return 0;
        }
    }

    int map_worker_to_writer(size_t rank) {
        if (is_compression_mode) {
            return num_io_processes - 1;
        } else {
            return 1 + (rank - num_io_processes) % (num_io_processes - 1);
        }
    }

    void workerProcess(size_t rank) {
        double total_read_time;
        size_t chunk_size = dimension[0] * dimension[1] * depth * sizeof(T) + sizeof(size_t) * 3;
        std::vector<char> buffer(chunk_size);

        while (true) {
            MPI_Recv(buffer.data(), chunk_size, MPI_BYTE, map_worker_to_reader(rank), 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char* pointer = buffer.data();
            size_t buffer_size = *reinterpret_cast<size_t*>(pointer);
            debugStream << "worker with rank " << rank << " received size: " << buffer_size
                        << std::endl;
            pointer += sizeof(size_t);
            if (buffer_size == 0) {
                printf("worker process %ld has completed!\n", rank);
                signalWorkerCompleteToWriter(rank, map_worker_to_writer(rank));
                return;
            }
            size_t sequence_number = *reinterpret_cast<size_t*>(pointer);
            pointer += sizeof(size_t);
            size_t layer_depth = *reinterpret_cast<size_t*>(pointer);
            pointer += sizeof(size_t);
            DataChunk<T> chunk;
            chunk.sequenceNumber = sequence_number;
            chunk.conf.absErrorBound = eb;
            std::vector<size_t> dims = {layer_depth, dimension[1], dimension[0]};

            chunk.conf.setDims(dims.begin(), dims.end());
            int dest_writer = map_worker_to_writer(rank);
            if (is_compression_mode) {
                chunk.dataBuffer = std::vector<T>(buffer_size);
                memcpy(chunk.dataBuffer.data(), pointer, buffer_size * sizeof(T));
                size_t outSize = 0;
                char* compressedData = SZ_compress<T>(chunk.conf, chunk.dataBuffer.data(), outSize);
                chunk.cpdataBuffer = std::vector<char>(outSize);
                memcpy(chunk.cpdataBuffer.data(), compressedData, outSize);
                delete[] compressedData;
                sendDataToWriter(dest_writer, chunk, rank);
            } else {
                chunk.cpdataBuffer = std::vector<char>(buffer_size);
                memcpy(chunk.cpdataBuffer.data(), pointer, buffer_size);
                T* decData = SZ_decompress<T>(chunk.conf, chunk.cpdataBuffer.data(), buffer_size);
                chunk.dataBuffer = std::vector<T>(chunk.conf.num);
                memcpy(chunk.dataBuffer.data(), decData, chunk.conf.num * sizeof(T));
                debugStream << "worker with rank " << rank << " finished decompressing chunk "
                            << chunk.sequenceNumber << std::endl;
                int64_t offset = sequence_number * dimension[0] * dimension[1] * depth * sizeof(T);
                debugStream << "[worker] first 10 values in chunk [" << chunk.sequenceNumber
                            << "] with offset " << offset << ": ";
                for (int _index_test = 0; _index_test < 10; _index_test++) {
                    debugStream << chunk.dataBuffer[_index_test] << " ";
                }
                debugStream << std::endl;
                debugStream << "[worker] last 10 values in chunk [" << chunk.sequenceNumber
                            << "] with offset " << offset << ": ";
                for (int _index_test = chunk.dataBuffer.size() - 10;
                     _index_test < chunk.dataBuffer.size(); _index_test++) {
                    debugStream << chunk.dataBuffer[_index_test] << " ";
                }
                debugStream << std::endl;
                sendDataToWriter(dest_writer, chunk, rank);
                delete[] decData;
            }
        }
        debugStream << "[exit] worker with rank " << rank << "finished and exited!" << std::endl;
    }

    void writerProcess(size_t rank) {
        MPI_File fh;
        MPI_Status status;

        MPI_File_open(MPI_COMM_SELF, output_file.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &fh);
        size_t chunk_cp_size = dimension[0] * dimension[1] * depth * sizeof(T) + sizeof(size_t) * 3;

        std::vector<char> buffer(chunk_cp_size);
        DataChunk<T> chunk;

        size_t end_counter = 0;
        MPI_Offset offset = 0;
        debugStream << "writer with rank " << rank << " start running!" << std::endl;
        SZ3::Timer writer_timer(false);
        double total_write_time = 0;

        std::vector<int> managed_workers;
        for (int worker_rank = num_io_processes; worker_rank < total_ranks; worker_rank++) {
            if (map_worker_to_writer(worker_rank) == rank) {
                managed_workers.push_back(worker_rank);
            }
        }
        int total_managed_workers = managed_workers.size();
        std::vector<bool> worker_finished(total_managed_workers, false);
        if (is_compression_mode) {
            while (true) {
                for (int i = 0; i < managed_workers.size(); i++) {
                    int src = managed_workers[i];
                    if (worker_finished[i]) {
                        continue;
                    }
                    // Receive a chunk
                    MPI_Recv(buffer.data(), chunk_cp_size, MPI_BYTE, src, 1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    char* pointer = buffer.data();
                    size_t buffer_size = *reinterpret_cast<size_t*>(pointer);
                    pointer += sizeof(size_t);
                    if (buffer_size == 0) {
                        worker_finished[i] = true;
                        end_counter++;
                        continue;
                    }
                    size_t sequence_number = *reinterpret_cast<size_t*>(pointer);
                    pointer += sizeof(size_t);
                    chunk.sequenceNumber = sequence_number;
                    chunk.cpdataBuffer = std::vector<char>(buffer_size);
                    memcpy(chunk.cpdataBuffer.data(), pointer, buffer_size);

                    // Calculate the offset for writing

                    writer_timer.start();

                    int64_t final_chunk_cp_size = buffer_size;
                    MPI_File_write_at_all(fh, offset, &final_chunk_cp_size, sizeof(int64_t),
                                          MPI_BYTE, &status);
                    offset += sizeof(int64_t);
                    // Write data using MPI collective I/O
                    debugStream << "[compress]writer with rank " << rank << " writing chunk "
                                << chunk.sequenceNumber << std::endl;
                    MPI_File_write_at_all(fh, offset, chunk.cpdataBuffer.data(),
                                          chunk.cpdataBuffer.size(), MPI_BYTE, &status);
                    offset += chunk.cpdataBuffer.size();

                    total_write_time += writer_timer.stop();
                }

                if (end_counter == managed_workers.size()) {
                    break;
                }
            }
        } else {

            while (true) {
                for (int i = 0; i < total_managed_workers; i++) {
                    int src = managed_workers[i];
                    if (worker_finished[i]) {
                        continue;
                    }
                    MPI_Recv(buffer.data(), chunk_cp_size, MPI_BYTE, src, 1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    char* pointer = buffer.data();
                    size_t data_buffer_num_elements = *reinterpret_cast<size_t*>(pointer);
                    pointer += sizeof(size_t);
                    if (data_buffer_num_elements == 0) {
                        worker_finished[i] = true;
                        end_counter++;
                        continue;
                    }
                    size_t sequence_number = *reinterpret_cast<size_t*>(pointer);
                    pointer += sizeof(size_t);
                    chunk.sequenceNumber = sequence_number;
                    chunk.dataBuffer = std::vector<T>(data_buffer_num_elements);
                    memcpy(chunk.dataBuffer.data(), pointer, data_buffer_num_elements * sizeof(T));
                    offset = sequence_number * dimension[0] * dimension[1] * depth * sizeof(T);
                    debugStream << "[decompress] first 10 values in chunk [" << chunk.sequenceNumber
                                << "] with offset " << offset << ": ";
                    for (int _index_test = 0; _index_test < 10; _index_test++) {
                        debugStream << chunk.dataBuffer[_index_test] << " ";
                    }
                    debugStream << std::endl;
                    debugStream << "[decompress] last 10 values in chunk [" << chunk.sequenceNumber
                                << "] with offset " << offset << ": ";
                    for (int _index_test = chunk.dataBuffer.size() - 10;
                         _index_test < chunk.dataBuffer.size(); _index_test++) {
                        debugStream << chunk.dataBuffer[_index_test] << " ";
                    }
                    debugStream << std::endl;

                    debugStream << "[decompress]writer with rank " << rank << " writing chunk "
                                << chunk.sequenceNumber << std::endl;
                    writer_timer.start();
                    MPI_File_write_at_all(fh, offset, chunk.dataBuffer.data(),
                                          chunk.dataBuffer.size() * sizeof(T), MPI_BYTE, &status);
                    total_write_time += writer_timer.stop();
                }
                if (end_counter == total_managed_workers) {
                    break;
                }
            }
        }
        // Close the file
        debugStream << "[write stats] writer with rank " << rank
                    << " finished, total write time: " << total_write_time << std::endl;
        MPI_File_close(&fh);
        debugStream << "[exit] writer with rank " << rank << " exited!" << std::endl;
    }
};

} // namespace sz3_split

#endif // SZ3_COMPRESSIONMPIMANAGER_H
