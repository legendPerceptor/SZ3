//
// Created by yuanjian on 3/29/24.
//

#ifndef SZ3_COMPRESSIONTHREADMANAGER_H
#define SZ3_COMPRESSIONTHREADMANAGER_H

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <chrono>
#include <utility>
#include <vector>
#include "SZ3/api/sz.hpp"

using namespace std::chrono_literals;

namespace sz3_split {

    SZ3::Config defaultConfig() {
        SZ3::Config conf;
        conf.cmprAlgo = SZ3::ALGO_LORENZO_REG;
        conf.lorenzo = true; // only use 1st order lorenzo
        conf.lorenzo2 = false;
        conf.regression = false;
        conf.regression2 = false;
        conf.errorBoundMode = SZ3::EB_ABS; // refer to def.hpp for all supported error bound mode
        conf.absErrorBound = 1E-3; // absolute error bound 1e-3
        return conf;
    }


    template<typename T>
    struct DataChunk {
        size_t id;
        size_t sequenceNumber;
        std::vector<T> dataBuffer;
        std::vector<char> cpdataBuffer;
        SZ3::Config conf;
    };

    template<typename T>
    class CompressionThreadManager {
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
        T eb;
        size_t depth;
        size_t num_of_threads;
        size_t read_queue_max_size;
        size_t write_queue_max_size;
        bool is_compression_mode;

        void readThread(){
            assert(dimension.size() == 3);
            SZ3::Config conf = defaultConfig();
            if (depth == 1) {
                conf.setDims(dimension.begin(), dimension.end() - 1);
            } else {
                std::vector<size_t> chunk_dimension = {dimension[0], dimension[1], depth};
                conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
            }
            size_t chunk_size = sizeof(T) * conf.num;
            conf.absErrorBound = eb;
            std::ifstream fin(input_file.c_str(), std::ios::binary);
            if(!fin.is_open()) {
                std::cerr << "Error opening the file: " << input_file.c_str() << std::endl;
                exit(EXIT_FAILURE);
            }

            double total_read_time = 0;
            size_t num_iterations = dimension[2] / depth;
            size_t leftover = dimension[2] % depth;
            if(leftover > 0) {
                num_iterations += 1;
            }
            SZ3::Timer temp(false);
            for(size_t i = 0;i<num_iterations;i++) {
                if (i == num_iterations - 1 && leftover > 0) {
                    std::vector<size_t> chunk_dimension = {dimension[0], dimension[1], leftover};
                    conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
                    chunk_size = sizeof(T) * conf.num;
                }
                DataChunk<T> chunk;
                chunk.id = i;
                chunk.sequenceNumber = i;
                chunk.conf = conf;


                if (is_compression_mode) {
                    chunk.dataBuffer = std::vector<T>(conf.num);
                    temp.start();
                    fin.read(reinterpret_cast<char *>(chunk.dataBuffer.data()), chunk_size);
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
//                    printf("chunk %u, read queue size: %u\n", i, readQueue.size());
                    reader_can_read.wait(lock, [this]{
                        return readQueue.size() < read_queue_max_size;
                    });
//                    printf("putting chunk %u into the read queue\n", i);
                    readQueue.push(chunk);
                }
                worker_can_get.notify_all();
            }
            // indicate that reading is done
            {
                std::lock_guard<std::mutex> lock(readMutex);
                printf("finished all reading, total read time is %f\n", total_read_time);
                read_done = true;
            }
            worker_can_put.notify_all();
        }


        void workerThread() {
            while (true) {
                // Wait for data to be available
                std::unique_lock<std::mutex> readLock(readMutex);
                while(readQueue.empty()) {
                    if(worker_can_get.wait_for(readLock, 100ms) == std::cv_status::timeout) {
                        if(read_done){
                            return;
                        }
                    }
                }
//                worker_can_get.wait(readLock, [this]{return !readQueue.empty() || read_done;});

                // Check if reading is done and no more data is available
                if(read_done && readQueue.empty())
                    break;

                // Retrieve a chunk of data from the read queue
                DataChunk<T> chunk = readQueue.front();
                readQueue.pop();

                readLock.unlock();
                reader_can_read.notify_one();

                //Compress the data chunk
                if(is_compression_mode) {
                    size_t outSize = 0;
                    printf("start to compress chunk %lu, chunk_size: %lu\n", chunk.sequenceNumber, sizeof(T) * chunk.dataBuffer.size());
                    char *compressedData = SZ_compress<T>(chunk.conf, chunk.dataBuffer.data(), outSize);
                    printf("finished compressing chunk %lu, compressed_chunk_size: %lu\n", chunk.sequenceNumber, outSize);
                    chunk.cpdataBuffer = std::vector<char>(outSize);
                    std::copy(compressedData, compressedData + outSize, chunk.cpdataBuffer.begin());
                    delete[] compressedData;
                } else {
                    printf("start to decompress chunk %lu, compressed chunk size: %lu\n", chunk.sequenceNumber, chunk.cpdataBuffer.size());
                    T* decData = SZ_decompress<T>(chunk.conf, chunk.cpdataBuffer.data(), chunk.cpdataBuffer.size());
                    size_t decompressed_size = chunk.conf.num * sizeof(T);
                    printf("finished running SZ_decompress with decompressed_size: %lu\n", decompressed_size);
                    chunk.dataBuffer = std::vector<T>(chunk.conf.num);
                    std::copy(decData, decData + chunk.conf.num, chunk.dataBuffer.begin());
                    printf("finished decompressing chunk %lu, decompressed_chunk_size: %lu\n", chunk.sequenceNumber, decompressed_size);
                    delete[] decData;
                }
//                printf("finished compressing chunk %u\n", chunk.sequenceNumber);
                size_t sequence_number = chunk.sequenceNumber;
                // Add the compressed data to the write queue
                {
                    std::unique_lock<std::mutex> lock(writeMutex);
                    worker_can_put.wait(lock, [this, &sequence_number] {
                        return writeQueue.size() < write_queue_max_size && expectedSequenceNumber == sequence_number;
                    });
//                    printf("putting chunk %u into the write queue\n", chunk.sequenceNumber);
                    writeQueue.push(chunk);
                    ++expectedSequenceNumber;
//                    printf("The next expected sequence number is %u\n", expectedSequenceNumber);
                }

                // Signal the write thread that new data is available
                writer_can_write.notify_one();
            }
        }


        void writeThread() {
            std::ofstream fout(output_file.c_str(), std::ios::binary | std::ios::out);
            if(!fout.is_open()) {
                std::cerr << "Error opening the file to write: " << output_file.c_str() << std::endl;
                exit(EXIT_FAILURE);
            }
            while(true) {
                std::unique_lock<std::mutex> writeLock(writeMutex);
                writer_can_write.wait(writeLock, [this] { return !writeQueue.empty() || all_done; });

                if (all_done && writeQueue.empty()) {
                    break;
                }

                while(!writeQueue.empty()) {
                    DataChunk<T> chunk = writeQueue.front();
                    std::cout << "Writing chunk with sequence number " << chunk.sequenceNumber << std::endl;
                    if(is_compression_mode) {
                        int64_t compressed_chunk_size = chunk.cpdataBuffer.size();
                        fout.write(reinterpret_cast<const char *>(&compressed_chunk_size), sizeof(int64_t));
                        fout.write(reinterpret_cast<const char *>(chunk.cpdataBuffer.data()), compressed_chunk_size);
                    } else {
                        fout.write(reinterpret_cast<const char*>(chunk.dataBuffer.data()), chunk.conf.num * sizeof(T));
                    }
                    writeQueue.pop();
                }
                worker_can_put.notify_all();
                writeLock.unlock();
            }

        }

    public:
        CompressionThreadManager(std::string input_file, std::string output_file, std::vector<size_t> dimension, T eb,
                            size_t depth, size_t num_of_threads, bool is_compression)
                : input_file(std::move(input_file)), output_file(std::move(output_file)),
                  dimension(std::move(dimension)), eb(eb), depth(depth), num_of_threads(num_of_threads){
            read_queue_max_size = num_of_threads;
            write_queue_max_size = num_of_threads;
            read_done = false;
            all_done = false;
            expectedSequenceNumber = 0;
            is_compression_mode = is_compression;
        }

        void startThreads() {
            std::thread reader(&CompressionThreadManager::readThread, this);
            std::thread writer(&CompressionThreadManager::writeThread, this);
            const int numWorkers = num_of_threads;
            std::thread workers[numWorkers];
            for (int i = 0; i < numWorkers; ++i) {
                workers[i] = std::thread(&CompressionThreadManager::workerThread, this);
            }

            // Join threads
            reader.join();
            std::cout << "Reader has joined!" << std::endl;
            for (int i = 0; i < numWorkers; ++i) {
                workers[i].join();
                std::cout << "workers[" << i << "] has joined" << std::endl;
            }
            all_done = true;
            writer_can_write.notify_one();
            std::cout << "Exit all threads and (de)compression finished" << std::endl;
            writer.join();
        }
    };

}
#endif //SZ3_COMPRESSIONTHREADMANAGER_H
