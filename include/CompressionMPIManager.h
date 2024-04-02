//
// Created by yuanjian on 4/1/24.
//

#ifndef SZ3_COMPRESSIONMPIMANAGER_H
#define SZ3_COMPRESSIONMPIMANAGER_H


#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <mpi.h>
#include "SZ3/api/sz.hpp"
#include "split_common.h"

namespace sz3_split {

    template<typename T>
    class CompressionMPIManager {

        std::string input_file;
        std::string output_file;
        std::vector<size_t> dimension;
        size_t num_io_processes;
        T eb;
        size_t depth;
        bool is_compression_mode;
        size_t total_ranks;

    public:
        CompressionMPIManager(std::string input_file, std::string output_file, std::vector<size_t> dimension, T eb,
                              size_t depth, bool is_compression, size_t num_io_processes)
                : input_file(std::move(input_file)), output_file(std::move(output_file)),
                  dimension(std::move(dimension)), eb(eb), depth(depth), is_compression_mode(is_compression), num_io_processes(num_io_processes) {
        }

        void startMPI() {
            MPI_Init(NULL, NULL);
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


            if (is_compression_mode) {
                if (rank >= 0 && rank < num_io_processes - 1) {
                    std::cout << "read rank: " << rank << ", size:" << size << std::endl;
                    readProcess(rank, num_io_processes - 1);
                } else if(rank == num_io_processes - 1) {
                    std::cout << "write rank: " << rank << ", size:" << size << std::endl;
                    writerProcess(rank);
                } else {
                    std::cout << "worker rank: " << rank << ", size:" << size << std::endl;
                    workerProcess(rank);
                }
            } else {
                if(rank == 0) {
                    readProcess(rank, 1);
                } else if (rank >= 1 && rank < num_io_processes) {
                    writerProcess(rank);
                } else {
                    workerProcess(rank);
                }
            }
            end_time = MPI_Wtime();
            double elapsed_time = end_time - start_time;
            std::cout << "[stats] rank<" << rank << "> total run time: " << elapsed_time << "seconds" << std::endl;
            MPI_Finalize();
        }

    private:

        void readProcess(size_t rank, size_t number_of_read_processes) {
            assert(dimension.size() == 3);
            SZ3::Config conf;
            if (depth == 1) {
                conf.setDims(dimension.begin(), dimension.end() - 1);
            } else {
                std::vector<size_t> chunk_dimension = {depth, dimension[0], dimension[1]};
                conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
            }
            size_t chunk_size = sizeof(T) * conf.num;
            conf.absErrorBound = eb;

            MPI_File fh;
            MPI_Offset file_size;
            int status;

            status = MPI_File_open(MPI_COMM_SELF, input_file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
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
            std::unordered_set<int> worker_ranks;
            std::cout << "read process with rank " << rank << " start working!" << std::endl;
            for (size_t i = 0; i < num_iterations; ++i) {
                if (i == num_iterations - 1 && leftover > 0) {
                    std::vector<size_t> chunk_dimension = {leftover, dimension[0], dimension[1]};
                    conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
                    chunk_size = sizeof(T) * conf.num;
                }
                DataChunk<T> chunk;
                chunk.id = rank + i * number_of_read_processes;
                chunk.sequenceNumber = rank + i * number_of_read_processes;
                chunk.conf = conf;

                if (is_compression_mode) {
                    chunk.dataBuffer = std::vector<T>(conf.num);
                    // Calculate start position for this thread, the last chunk may have different chunk_size
                    temp.start();
                    // Calculate the offset for reading
                    MPI_Offset offset = (rank + i * number_of_read_processes) * org_chunk_size;
                    // Set the file view for each process
                    MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
                    MPI_File_read_all(fh, chunk.dataBuffer.data(), conf.num * sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
                    // Read data using MPI I/O
                    total_read_time += temp.stop();
                } else {
                    temp.start();
                    int64_t compressed_chunk_size;
                    MPI_File_set_view(fh, dp_offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
                    MPI_File_read_all(fh, &compressed_chunk_size, sizeof(int64_t), MPI_BYTE, MPI_STATUS_IGNORE);
                    chunk.cpdataBuffer = std::vector<char>(compressed_chunk_size);
                    MPI_File_read_all(fh, chunk.cpdataBuffer.data(), compressed_chunk_size, MPI_BYTE, MPI_STATUS_IGNORE);
                    dp_offset += sizeof(int64_t) + compressed_chunk_size;
                    total_read_time += temp.stop();
                }

                size_t num_workers = total_ranks - num_io_processes;
                int dest = num_io_processes + (chunk.id % num_workers);
                worker_ranks.insert(dest);

                std::cout << "read process with rank " << rank << " sending chunk "
                          << chunk.sequenceNumber << " to worker " << dest << std::endl;
                if (i == num_iterations - 1 && leftover > 0) {
                    sendDataToWorker(chunk, rank, leftover);
                } else {
                    sendDataToWorker(chunk, rank, depth);
                }
            }

            std::cout << "[read stats]read process with rank " << rank << " finished. Total read time:" << total_read_time << std::endl;

            for(auto dest : worker_ranks) {
                signalReadCompleteToWorkers(rank, dest);
            }
            
            MPI_File_close(&fh);
        }

        void sendDataToWorker(DataChunk<T>& chunk, size_t read_rank, size_t layer_depth) {
            size_t num_workers = total_ranks - num_io_processes;
            int dest = num_io_processes + (chunk.id % num_workers); // read 0 1 write 2;   0-0->3  1-1->4  0-3->5
            int tag = 0;
            // Serialize chunk
            if(is_compression_mode) {
                std::vector<char> buffer(chunk.dataBuffer.size() * sizeof(T) + 3 * sizeof(size_t));
                char *pointer = buffer.data();
                size_t chunk_buffer_size = chunk.dataBuffer.size();
                memcpy(pointer, &chunk_buffer_size, sizeof(size_t));
                pointer += sizeof(size_t);
                memcpy(pointer, &chunk.sequenceNumber, sizeof(size_t));
                pointer += sizeof(size_t);
                memcpy(pointer, &layer_depth, sizeof(size_t));
                pointer += sizeof(size_t);

                memcpy(pointer, chunk.dataBuffer.data(), chunk.dataBuffer.size() * sizeof(T));
//                std::cout << "read rank: " << read_rank << ", about to MPI_Send, buffer size:" << buffer.size() << std::endl;
                // Send chunk to worker
                MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
//                std::cout << "read rank: " << read_rank << ", finished MPI_Send, buffer size:" << buffer.size() << std::endl;
            } else {
                std::vector<char> buffer(chunk.cpdataBuffer.size() + 3 * sizeof(size_t));
                char *pointer = buffer.data();
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

        void sendDataToWriter(DataChunk<T>& chunk, size_t rank) {
            size_t num_workers = total_ranks - num_io_processes;
            int dest;
            int tag = 1;
            if (is_compression_mode) {
                dest = map_worker_to_writer(rank);
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
                dest = map_worker_to_writer(rank);
                std::vector<char> buffer(chunk.dataBuffer.size() * sizeof(T) + 2 * sizeof(size_t));
                char* pointer = buffer.data();
                size_t chunk_buffer_size = chunk.cpdataBuffer.size();
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
            size_t num_workers = total_ranks - num_io_processes;
            if(is_compression_mode) {
                return (rank - num_io_processes) % (num_io_processes - 1);
            } else {
                return 0;
            }
        }

        int map_worker_to_writer(size_t rank) {
            if(is_compression_mode) {
                return num_io_processes - 1;
            } else {
                return 1 + (rank - num_io_processes) % (num_io_processes - 1);
            }
        }

        void workerProcess(size_t rank) {
            double total_read_time;
            size_t chunk_size = dimension[0] * dimension[1] * depth * sizeof(T) + sizeof(size_t) * 3;
            std::vector<char> buffer(chunk_size);

            while(true) {
//                std::cout << "worker with rank " << rank << ", src:" << map_worker_to_reader(rank) << ", chunk_size to receive: " << chunk_size << std::endl;
                MPI_Recv(buffer.data(), chunk_size, MPI_BYTE, map_worker_to_reader(rank), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* pointer = buffer.data();
                size_t buffer_size = *reinterpret_cast<size_t*>(pointer);
//                std::cout << "worker with rank " << rank << " recevied size: " << buffer_size << std::endl;
                pointer += sizeof (size_t);
                if (buffer_size == 0) {
                    printf("worker process %ld has completed!\n", rank);
                    signalWorkerCompleteToWriter(rank, map_worker_to_writer(rank));
                    return;
                }
                size_t sequence_number = *reinterpret_cast<size_t*>(pointer);
                pointer += sizeof (size_t);
                size_t layer_depth = *reinterpret_cast<size_t*>(pointer);
                pointer += sizeof (size_t);
                DataChunk<T> chunk;
                chunk.sequenceNumber = sequence_number;
                chunk.conf.absErrorBound = eb;
                std::vector<size_t> dims = {layer_depth, dimension[0], dimension[1]};

                chunk.conf.setDims(dims.begin(), dims.end());
//                std::cout << "decompression dims size: " << dims.size() << ", chunk.conf.N=" << (int)chunk.conf.N << ", sequence_number: " << sequence_number << ", layer_depth:" << layer_depth << std::endl;

                if (is_compression_mode) {
                    chunk.dataBuffer = std::vector<T>(buffer_size);
                    memcpy(chunk.dataBuffer.data(), pointer, buffer_size * sizeof(T));
                    size_t outSize = 0;
                    char *compressedData = SZ_compress<T>(chunk.conf, chunk.dataBuffer.data(), outSize);
                    chunk.cpdataBuffer = std::vector<char>(outSize);
//                    std::cout << "work with rank " << rank << " has outSide: " << outSize << std::endl;
//                    std::copy(compressedData, compressedData + outSize, chunk.cpdataBuffer.begin());
                    memcpy(chunk.cpdataBuffer.data(), compressedData, outSize);
                    delete[] compressedData;
//                    std::cout << "worker with rank " << rank << " finished compressing chunk " << chunk.sequenceNumber << std::endl;
                    sendDataToWriter(chunk, rank);
                } else {
                    chunk.cpdataBuffer = std::vector<char>(buffer_size);
                    memcpy(chunk.cpdataBuffer.data(), pointer, buffer_size);
//                    std::cout << "worker with rank " << rank << " start decompressing, with compressed data size: " << buffer_size  << ", conf.num=" << chunk.conf.num << ", conf.N=" << (int)chunk.conf.N << std::endl;
                    T *decData = SZ_decompress<T>(chunk.conf, chunk.cpdataBuffer.data(), buffer_size);
                    chunk.dataBuffer = std::vector<T>(chunk.conf.num);
//                    std::copy(decData, decData + chunk.conf.num, chunk.dataBuffer.begin());
                    memcpy(chunk.dataBuffer.data(), decData, chunk.conf.num * sizeof(T));
                    std::cout << "worker with rank " << rank << " finished decompressing chunk " << chunk.sequenceNumber << std::endl;
                    sendDataToWriter(chunk, rank);
                    delete[] decData;
                }

            }

        }

        void writerProcess(size_t rank){
            MPI_File fh;
            MPI_Status status;

            MPI_File_open(MPI_COMM_SELF, output_file.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
            size_t chunk_cp_size = dimension[0] * dimension[1] * depth * sizeof(T) + sizeof(size_t) * 3;

            std::vector<char> buffer(chunk_cp_size);
            DataChunk<T> chunk;

            size_t end_counter = 0;
            size_t num_workers = total_ranks - num_io_processes;
            std::vector<bool> worker_finished(total_ranks - num_io_processes, false);
            MPI_Offset offset = 0;
            std::cout << "writer with rank " << rank << " start running!" << std::endl;
            SZ3::Timer writer_timer(false);
            double total_write_time = 0;
            if(is_compression_mode) {
                while (true) {
                    for (int src = num_io_processes; src < total_ranks; src++) {
                        if (worker_finished[src - num_io_processes]) {
                            continue;
                        }
//                        std::cout << "writer with rank " << rank << " receiving from " << src << " with maximum chunk size: " << chunk_cp_size << std::endl;
                        // Receive a chunk
                        MPI_Recv(buffer.data(), chunk_cp_size, MPI_BYTE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                        std::cout << "writer with rank " << rank << " receiving from " << src << "finished MPI_Recv"<< std::endl;
                        char *pointer = buffer.data();
                        size_t buffer_size = *reinterpret_cast<size_t *>(pointer);
//                        std::cout << "writer with rank " << rank << " receiving from " << src << " read buffer size: " << buffer_size << std::endl;
                        pointer += sizeof(size_t);
                        if (buffer_size == 0) {
                            worker_finished[src - num_io_processes] = true;
                            end_counter++;
                        }
                        size_t sequence_number = *reinterpret_cast<size_t *>(pointer);
                        pointer += sizeof(size_t);
                        chunk.sequenceNumber = sequence_number;
                        chunk.cpdataBuffer = std::vector<char>(buffer_size);
                        memcpy(chunk.cpdataBuffer.data(), pointer, buffer_size);

                        // Calculate the offset for writing

                        writer_timer.start();

                        int64_t final_chunk_cp_size = buffer_size;
                        MPI_File_write_at_all(fh, offset, &final_chunk_cp_size, sizeof(int64_t), MPI_BYTE,
                                              &status);
                        offset += sizeof(int64_t);
                        // Write data using MPI collective I/O
                        std::cout << "[compress]writer with rank " << rank << " writing chunk " << chunk.sequenceNumber << std::endl;
                        MPI_File_write_at_all(fh, offset, chunk.cpdataBuffer.data(), chunk.cpdataBuffer.size(),
                                              MPI_BYTE,
                                              &status);
                        offset += chunk.cpdataBuffer.size();

                        total_write_time += writer_timer.stop();
                    }

                    if (end_counter == num_workers) {
                        break;
                    }
                }
            } else {
                int num_of_writers = num_io_processes - 1; //
                int total_worker_for_this_writer = (total_ranks - (num_io_processes + rank - 1)) / num_of_writers;
                while(true) {
                    for(int src = num_io_processes + rank - 1; src < total_ranks; src += num_of_writers) {
                        if (worker_finished[src - num_io_processes]) {
                            continue;
                        }
                        MPI_Recv(buffer.data(), chunk_cp_size, MPI_BYTE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        char *pointer = buffer.data();
                        size_t buffer_size = *reinterpret_cast<size_t *>(pointer);
                        pointer += sizeof(size_t);
                        if (buffer_size == 0) {
                            worker_finished[src - num_io_processes] = true;
                            end_counter++;
                        }
                        size_t sequence_number = *reinterpret_cast<size_t *>(pointer);
                        pointer += sizeof(size_t);
                        chunk.sequenceNumber = sequence_number;
                        chunk.dataBuffer = std::vector<T>(buffer_size);
                        memcpy(chunk.dataBuffer.data(), pointer, buffer_size);
                        offset = sequence_number * dimension[0] * dimension[1] * depth;
                        std::cout << "[decompress]writer with rank " << rank << " writing chunk " << chunk.sequenceNumber << std::endl;
                        writer_timer.start();
                        MPI_File_write_at_all(fh, offset, chunk.dataBuffer.data(), chunk.dataBuffer.size() * sizeof(T),
                                              MPI_BYTE, &status);
                        total_write_time += writer_timer.stop();
                    }
                    if(end_counter == total_worker_for_this_writer) {
                        break;
                    }
                }
            }
            // Close the file
            std::cout << "[write stats]writer with rank " << rank << " finished, total write time: " << total_write_time << std::endl;
            MPI_File_close(&fh);
        }
    };

}

#endif //SZ3_COMPRESSIONMPIMANAGER_H
