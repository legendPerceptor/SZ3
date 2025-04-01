//
// Created by yuanjian on 3/29/24.
//

#ifndef SZ3_COMPRESSIONTHREADMANAGER_H
#define SZ3_COMPRESSIONTHREADMANAGER_H

#include "DebugStream.hpp"
#include "SZ3/api/sz.hpp"
#include "split_common.h"
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

using namespace std::chrono_literals;

namespace sz3_split {

template <typename T> class CompressionThreadManager {
    std::queue<DataChunk<T>> readQueue;
    std::queue<DataChunk<T>> writeQueue;
    std::mutex readMutex;
    std::mutex writeMutex;
    std::condition_variable reader_can_read, worker_can_get;
    std::condition_variable writer_can_write, worker_can_put;
    bool read_done, all_done;
    size_t expectedSequenceNumber;

    std::string input_file;
    std::string output_file;
    std::vector<size_t> dimension;
    double eb;
    size_t depth;
    size_t num_of_threads;
    size_t num_of_readers;
    size_t read_queue_max_size;
    size_t write_queue_max_size;
    bool is_compression_mode;
    bool use_logscale;
    int skip_header_size;

    void readThread(size_t threadId, size_t totalReadThreads) {
        assert(dimension.size() == 3);
        SZ3::Config conf = defaultConfig();
        if (depth == 1) {
            conf.setDims(dimension.begin(), dimension.end() - 1);
        } else {
            std::vector<size_t> chunk_dimension = {depth, dimension[1], dimension[0]};
            conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
        }
        size_t chunk_size = sizeof(T) * conf.num;
        conf.absErrorBound = eb;
        std::ifstream fin(input_file.c_str(), std::ios::binary);
        if (!fin.is_open()) {
            std::cerr << "Error opening the file: " << input_file.c_str() << std::endl;
            exit(EXIT_FAILURE);
        }
        if (skip_header_size > 0) {
            fin.seekg(skip_header_size, std::ios::beg);
        }

        double total_read_time = 0;
        size_t num_iterations = dimension[2] / depth;
        size_t leftover = dimension[2] % depth;
        if (leftover > 0) {
            num_iterations += 1;
        }
        SZ3::Timer temp(false);
        size_t org_chunk_size = chunk_size;
        for (size_t i = threadId; i < num_iterations; i += totalReadThreads) {
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
                size_t start_position = i * org_chunk_size * totalReadThreads;
                // Set file pointer to the calculated start position
                fin.seekg(skip_header_size + start_position, std::ios::beg);
                temp.start();
                fin.read(reinterpret_cast<char*>(chunk.dataBuffer.data()), chunk_size);
                total_read_time += temp.stop();
            } else {
                temp.start();
                int64_t compressed_chunk_size;
                fin.read(reinterpret_cast<char*>(&compressed_chunk_size), sizeof(int64_t));
                chunk.cpdataBuffer = std::vector<char>(compressed_chunk_size);
                fin.read(chunk.cpdataBuffer.data(), compressed_chunk_size);
                total_read_time += temp.stop();
            }

            {
                std::unique_lock<std::mutex> lock(readMutex);
                //                    printf("chunk %u, read queue size: %u\n", i,
                //                    readQueue.size());
                reader_can_read.wait(lock,
                                     [this] { return readQueue.size() < read_queue_max_size; });
                //                    printf("putting chunk %u into the read queue\n", i);
                readQueue.push(chunk);
            }
            worker_can_get.notify_all();
        }
        // indicate that reading is done
        {
            std::lock_guard<std::mutex> lock(readMutex);
            debugStream << "finished all reading, total read time is " << total_read_time
                        << std::endl;
            read_done = true;
        }
        worker_can_put.notify_all();
    }

    void workerThread() {
        while (true) {
            // Wait for data to be available
            std::unique_lock<std::mutex> readLock(readMutex);
            while (readQueue.empty()) {
                if (worker_can_get.wait_for(readLock, 100ms) == std::cv_status::timeout) {
                    if (read_done) {
                        return;
                    }
                }
            }
            //                worker_can_get.wait(readLock, [this]{return !readQueue.empty() ||
            //                read_done;});

            // Check if reading is done and no more data is available
            if (read_done && readQueue.empty())
                break;

            // Retrieve a chunk of data from the read queue
            DataChunk<T> chunk = readQueue.front();
            readQueue.pop();

            readLock.unlock();
            reader_can_read.notify_one();

            // Compress the data chunk
            if (is_compression_mode) {
                size_t outSize = 0;
                debugStream << "start to compress chunk [" << chunk.sequenceNumber
                            << "], chunk_size: " << sizeof(T) * chunk.dataBuffer.size()
                            << std::endl;
                char* compressedData;
                if (use_logscale) {
                    if constexpr (std::is_integral_v<T>) {
                        std::vector<float> floatBuffer(chunk.dataBuffer.size());
                        std::transform(chunk.dataBuffer.begin(), chunk.dataBuffer.end(),
                                       floatBuffer.begin(),
                                       [](double val) { return std::log(val); });
                        compressedData =
                            SZ_compress<float>(chunk.conf, floatBuffer.data(), outSize);
                    } else {
                        std::transform(chunk.dataBuffer.begin(), chunk.dataBuffer.end(),
                                       chunk.dataBuffer.begin(),
                                       [](T val) { return static_cast<T>(std::log(val)); });
                        compressedData =
                            SZ_compress<T>(chunk.conf, chunk.dataBuffer.data(), outSize);
                    }
                } else {
                    compressedData = SZ_compress<T>(chunk.conf, chunk.dataBuffer.data(), outSize);
                }
                chunk.cpdataBuffer = std::vector<char>(outSize);
                std::copy(compressedData, compressedData + outSize, chunk.cpdataBuffer.begin());
                delete[] compressedData;
            } else {
                debugStream << "start to decompress chunk [" << chunk.sequenceNumber
                            << "], compressed chunk size: " << chunk.cpdataBuffer.size()
                            << std::endl;

                chunk.dataBuffer = std::vector<T>(chunk.conf.num);
                if (use_logscale) {
                    if constexpr (std::is_integral_v<T>) {
                        auto* decData = SZ_decompress<float>(chunk.conf, chunk.cpdataBuffer.data(),
                                                             chunk.cpdataBuffer.size());
                        std::vector<float> floatDPBuffer(chunk.conf.num);
                        std::copy_n(decData, chunk.conf.num, floatDPBuffer.begin());
                        delete[] decData;
                        std::transform(floatDPBuffer.begin(), floatDPBuffer.end(), chunk.dataBuffer.begin(),
                                       [](float val) { return static_cast<T>(std::exp(val)); });
                    } else {
                        auto* decData = SZ_decompress<T>(chunk.conf, chunk.cpdataBuffer.data(),
                                                         chunk.cpdataBuffer.size());
                        std::copy(decData, decData + chunk.conf.num, chunk.dataBuffer.begin());
                        delete[] decData;
                        std::transform(chunk.dataBuffer.begin(), chunk.dataBuffer.end(), chunk.dataBuffer.begin(),
                                       [](T val) { return static_cast<T>(std::exp(val)); });
                    }
                } else {
                    auto* decData = SZ_decompress<T>(chunk.conf, chunk.cpdataBuffer.data(),
                                                     chunk.cpdataBuffer.size());
                    std::copy(decData, decData + chunk.conf.num, chunk.dataBuffer.begin());
                    delete[] decData;
                }

                size_t decompressed_size = chunk.conf.num * sizeof(T);
                debugStream << "finished decompressing chunk [" << chunk.sequenceNumber
                            << "], decompressed_chunk_size: " << decompressed_size << std::endl;
            }
            //                printf("finished compressing chunk %u\n", chunk.sequenceNumber);
            size_t sequence_number = chunk.sequenceNumber;
            // Add the compressed data to the write queue
            {
                std::unique_lock<std::mutex> lock(writeMutex);
                worker_can_put.wait(lock, [this, &sequence_number] {
                    return writeQueue.size() < write_queue_max_size &&
                           expectedSequenceNumber == sequence_number;
                });
                //                    printf("putting chunk %u into the write queue\n",
                //                    chunk.sequenceNumber);
                writeQueue.push(chunk);
                ++expectedSequenceNumber;
                //                    printf("The next expected sequence number is %u\n",
                //                    expectedSequenceNumber);
            }

            // Signal the write thread that new data is available
            writer_can_write.notify_one();
        }
    }

    void writeThread() {
        std::ofstream fout(output_file.c_str(), std::ios::binary | std::ios::out);
        if (!fout.is_open()) {
            std::cerr << "Error opening the file to write: " << output_file.c_str() << std::endl;
            exit(EXIT_FAILURE);
        }
        if (skip_header_size > 0) {
            std::ifstream fin(input_file.c_str(), std::ios::binary | std::ios::in);
            if (!fin.is_open()) {
                std::cerr << "Error opening the input file in the write thread: "
                          << input_file.c_str() << std::endl;
            }
            std::vector<char> header(skip_header_size);
            fin.read(header.data(), skip_header_size);
            fout.write(header.data(), skip_header_size);
        }
        while (true) {
            std::unique_lock<std::mutex> writeLock(writeMutex);
            writer_can_write.wait(writeLock, [this] { return !writeQueue.empty() || all_done; });

            if (all_done && writeQueue.empty()) {
                break;
            }

            while (!writeQueue.empty()) {
                DataChunk<T> chunk = writeQueue.front();
                debugStream << "Writing chunk with sequence number " << chunk.sequenceNumber
                            << std::endl;

                debugStream << is_compression_mode << "[writer] first 10 values in chunk ["
                            << chunk.sequenceNumber << "]" << ": ";
                for (int _index_test = 0; _index_test < 10; _index_test++) {
                    debugStream << chunk.dataBuffer[_index_test] << " ";
                }
                debugStream << std::endl;
                debugStream << is_compression_mode << "[writer] last 10 values in chunk ["
                            << chunk.sequenceNumber << "]" << ": ";
                for (int _index_test = chunk.dataBuffer.size() - 10;
                     _index_test < chunk.dataBuffer.size(); _index_test++) {
                    debugStream << chunk.dataBuffer[_index_test] << " ";
                }
                debugStream << std::endl;

                if (is_compression_mode) {
                    int64_t compressed_chunk_size = chunk.cpdataBuffer.size();
                    fout.write(reinterpret_cast<const char*>(&compressed_chunk_size),
                               sizeof(int64_t));
                    fout.write(reinterpret_cast<const char*>(chunk.cpdataBuffer.data()),
                               compressed_chunk_size);
                } else {
                    fout.write(reinterpret_cast<const char*>(chunk.dataBuffer.data()),
                               chunk.conf.num * sizeof(T));
                }
                writeQueue.pop();
            }
            worker_can_put.notify_all();
            writeLock.unlock();
        }
    }

  public:
    CompressionThreadManager(std::string input_file, std::string output_file,
                             std::vector<size_t> dimension, double eb, size_t depth,
                             size_t num_of_threads, size_t num_of_readers, bool is_compression,
                             bool use_logscale, int skip_header_size)
        : input_file(std::move(input_file)), output_file(std::move(output_file)),
          dimension(std::move(dimension)), eb(eb), depth(depth), num_of_threads(num_of_threads),
          num_of_readers(num_of_readers) {
        read_queue_max_size = num_of_threads;
        write_queue_max_size = num_of_threads;
        read_done = false;
        all_done = false;
        expectedSequenceNumber = 0;
        is_compression_mode = is_compression;
        this->use_logscale = use_logscale;
        this->skip_header_size = skip_header_size;
    }

    void startThreads() {

        const int numReaders = num_of_readers;
        const int numWorkers = num_of_threads;
        std::thread writer(&CompressionThreadManager::writeThread, this);
        std::thread workers[numWorkers];
        std::thread readers[numReaders];
        for (int i = 0; i < num_of_readers; ++i) {
            readers[i] = std::thread(&CompressionThreadManager::readThread, this, i, numReaders);
        }
        for (int i = 0; i < numWorkers; ++i) {
            workers[i] = std::thread(&CompressionThreadManager::workerThread, this);
        }

        // Join threads
        for (int i = 0; i < numReaders; ++i) {
            readers[i].join();
            debugStream << "reader[" << i << "] has joined!" << std::endl;
        }
        debugStream << "Reader has joined!" << std::endl;
        for (int i = 0; i < numWorkers; ++i) {
            workers[i].join();
            debugStream << "workers[" << i << "] has joined" << std::endl;
        }
        all_done = true;
        writer_can_write.notify_one();
        debugStream << "Exit all threads and (de)compression finished" << std::endl;
        writer.join();

        auto outSize = std::filesystem::file_size(output_file);
        auto inputSize = std::filesystem::file_size(input_file);
        std::cout << "output file size: " << outSize << ", input file size: " << inputSize
                  << std::endl;
        if (is_compression_mode) {
            std::cout << "compression ratio: " << (double)inputSize / outSize << std::endl;
        } else {
            std::cout << "compression ratio: " << (double)outSize / inputSize << std::endl;
        }
    }
};

} // namespace sz3_split
#endif // SZ3_COMPRESSIONTHREADMANAGER_H
