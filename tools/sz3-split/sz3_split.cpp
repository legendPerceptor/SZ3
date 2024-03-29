// Created by yuanjian on 11/23/23.
// This tool is dedicated for extremely large file compression.
// The file cannot fit in the memory, and we need to split it to compress with multiple nodes or threads
#include <getopt.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include "SZ3/api/sz.hpp"

namespace sz3_split {

    void
    parseCompressOptions(int argc, char **argv, int &threads, std::string &raw_file, std::string &output_file,
                         std::vector<size_t> &data_dimension, float& eb, bool& is_float64, std::string& mode, size_t& depth) {
        optind = 1;
        const char *opt_index = "ht:i:d:e:o:";
        const int FLOAT_64 = 1008;
        const int MODE = 1009;
        const int DEPTH = 1010;
        struct option opts[] = {
                {"threads",    required_argument, nullptr, 't'},
                {"help",       no_argument,       nullptr, 'h'},
                {"input",      required_argument, nullptr, 'i'},
                {"output", required_argument, nullptr, 'o'},
                {"errorbound", required_argument, nullptr, 'e'},
                {"dimension",  required_argument, nullptr, 'd'},
                {"float64", no_argument, nullptr, FLOAT_64},
                {"mode", required_argument, nullptr, MODE},
                {"depth", required_argument, nullptr, DEPTH}
        };
        depth = 1;
        std::string compress_helper_info = "Usage: sz3_split (de)compress [options]\n"
                                           "options:  --threads/-t     INT   number of threads, default is 1\n"
                                           "          --input/-i       STR   the RAW file/compressed file\n"
                                           "          --output/-o      STR   the compressed file/decompressed file location\n"
                                           "          --help/-h              print this help information\n"
                                           "          --dimension/-d   STR   the data dimension of the file, e.g., 256 256 512\n"
                                           "          --errorbound/-e  FLOAT the error bound to use in compression\n"
                                           "          --float64              the default is float32 for each datapoint, this param changes it to float64\n"
                                           "          --mode           STR   select 'layer' for layer-by-layer compression, 'direct' for direct compression\n"
                                           "          --depth          INT   select the layer depth in layer-by-layer compression\n";
        is_float64 = false; // default is float32
        int c;
        while ((c = getopt_long(argc, argv, opt_index, opts, nullptr)) != -1) {
            switch (c) {
                case 'e':
                    eb = std::stof(optarg);
                    break;
                case 'h':
                    std::cout << compress_helper_info << std::endl;
                    exit(EXIT_SUCCESS);
                    break;
                case 't':
                    threads = std::stoi(optarg);
                    break;
                case 'i':
                    raw_file = optarg;
                    break;
                case 'd':
                    optind--;
                    for (; optind < argc && *argv[optind] != '-'; optind++) {
                        data_dimension.push_back(std::stoi(argv[optind]));
                    }
                    break;
                case 'o':
                    output_file = optarg;
                    break;
                case FLOAT_64:
                    is_float64 = true;
                    break;
                case MODE:
                    mode = optarg;
                    break;
                case DEPTH:
                    depth = std::stoi(optarg);
                    break;
                default:
                    std::cerr << "Usage: compress/decompress/test\n";
                    std::cout << compress_helper_info << std::endl;
                    exit(EXIT_FAILURE);
            }
        }
        if (raw_file.empty() || output_file.empty()) {
            std::cerr << "Please specify a proper input/output file path" << std::endl;
            std::cout << compress_helper_info << std::endl;
            exit(EXIT_FAILURE);
        }
        if (data_dimension.empty()) {
            std::cerr << "Please input the correct data dimension" << std::endl;
            std::cout << compress_helper_info << std::endl;
            exit(EXIT_FAILURE);
        }

    }

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
    template<typename TYPE>
    int compress_impl(const std::string& input_file, const std::string& output_file, std::vector<size_t> dimension, TYPE eb, const std::string& mode, size_t depth) {
        SZ3::Config conf = defaultConfig();

        if(dimension.size() == 3 && mode == "layer") {
            if (depth == 1) {
                conf.setDims(dimension.begin(), dimension.end() - 1);
            } else {
                std::vector<size_t> chunk_dimension = {dimension[0], dimension[1], depth};
                conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
            }
            conf.absErrorBound = eb;
            std::ifstream fin(input_file.c_str(), std::ios::binary);
            if(!fin.is_open()) {
                std::cerr << "Error opening the file: " << input_file.c_str() << std::endl;
                return -1;
            }
//            std::streampos read_start = 0;
            size_t chunk_size = sizeof(TYPE) * conf.num;
            std::ofstream fout(output_file.c_str(), std::ios::binary | std::ios::out);
            if(!fout.is_open()) {
                std::cerr << "Error opening the file: " << output_file.c_str() << std::endl;
                return -1;
            }
            SZ3::Timer total_timer(true);
            SZ3::Timer temp(false);
            double total_read_time = 0, total_write_time = 0, total_compress_time = 0;
            size_t num_iterations = dimension[2] / depth;
            size_t leftover = dimension[2] % depth;
            if(leftover > 0) {
                num_iterations += 1;
            }
            for(size_t i = 0;i<num_iterations;i++) {
                if(i == num_iterations - 1 && leftover > 0) {
                    std::vector<size_t> chunk_dimension = {dimension[0], dimension[1], leftover};
                    conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
                }
                std::vector<TYPE> buffer(conf.num);
//                fin.seekg(read_start);
                temp.start();
                fin.read(reinterpret_cast<char *>(buffer.data()), chunk_size);
                total_read_time += temp.stop();
//                read_start += chunk_size;
                size_t outSize;
                SZ3::Timer timer(true);
                char *compressedData = SZ_compress<TYPE>(conf, buffer.data(), outSize);
                double compress_time = timer.stop();
                total_compress_time += compress_time;
                auto compresed_chunk_size = static_cast<int64_t>(outSize);
                temp.start();
                fout.write(reinterpret_cast<const char*>(&compresed_chunk_size), sizeof(int64_t));
                fout.write(reinterpret_cast<const char*>(compressedData), compresed_chunk_size);
                fout.flush();
                total_write_time += temp.stop();
                std::cout << "Chunk " << i << " compression completed! compressed_chunk_size: " << compresed_chunk_size << "; Time elasped: " << compress_time
                          << " seconds, total time elapsed: " << total_timer.stop() << " seconds." << "fout.tellp(): " << fout.tellp() <<std::endl;
            }
            size_t file_size = (size_t)dimension[0] * (size_t)dimension[1] * (size_t)dimension[2] * sizeof(TYPE);
            std::cout << "Congratulations! Compression completed! Total file size compressed: " << file_size / 1024 / 1024 << " MB;\n"
                      << "Total compressed file size generated: " << fout.tellp() / 1024 /1024 << "MB;" << std::endl
                      << "total time elapsed: " << total_timer.stop() << "seconds" << std::endl;
            double cr = (double)file_size / (double)fout.tellp();
            std::cout << "compression ratio:" << std::fixed << std::setprecision(2) << cr <<", compression_time: " << total_compress_time << ", read time: " << total_read_time
                      << ", write time: " << total_write_time << std::endl;
        } else {
            // 300 is the fastest dimension
            conf.setDims(dimension.begin(), dimension.end());
            conf.absErrorBound = eb;
            std::vector<TYPE> buffer(conf.num);
            SZ3::readfile<TYPE>(input_file.c_str(), conf.num, buffer.data());
            size_t outSize;
            SZ3::Timer timer(true);
            char *compressedData = SZ_compress<TYPE>(conf, buffer.data(), outSize);
            double compress_time = timer.stop();
            std::cout << "Compression completed! Time elasped: " << compress_time << std::endl;
            SZ3::writefile(output_file.c_str(), compressedData, outSize);
            printf("compression ratio = %.2f \n", (double) conf.num * 1.0 * sizeof(TYPE) / (double) outSize);
            printf("compression time = %f\n", compress_time);
            printf("compressed data file = %s\n", output_file.c_str());
        }
        return 0;
    }

    int compress(int argc, char **argv) {
        int threads;
        std::vector<size_t> dimension;
        std::string input_file, output_file;
        float eb;
        bool isfloat64;
        std::string mode;
        size_t depth;
        parseCompressOptions(argc, argv, threads, input_file, output_file, dimension, eb, isfloat64, mode, depth);
        if(isfloat64){
            return compress_impl<double>(input_file, output_file, dimension, eb, mode, depth);
        }else{
            return compress_impl<float>(input_file, output_file, dimension, eb, mode, depth);
        }
    }

    template<typename TYPE>
    int decompress_impl(const std::string& input_file, const std::string& output_file, std::vector<size_t> dimension, TYPE eb, const std::string& mode, size_t depth) {
        SZ3::Config conf = defaultConfig(); // 300 is the fastest dimension
        if (dimension.size() == 3 && mode == "layer") {
            if (depth == 1) {
                conf.setDims(dimension.begin(), dimension.end() - 1);
            } else {
                std::vector<size_t> chunk_dimension = {dimension[0], dimension[1], depth};
                conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
            }
            conf.absErrorBound = eb;
            std::ifstream fin(input_file.c_str(), std::ios::binary);
            if(!fin.is_open()) {
                std::cerr << "Error opening the file: " << input_file.c_str() << std::endl;
                return -1;
            }
            int64_t compressed_chunk_size;
            double total_read_time = 0, total_write_time = 0, total_decompress_time = 0;
            std::ofstream fout(output_file.c_str(), std::ios::binary | std::ios::out);
            SZ3::Timer total_timer(true);
            SZ3::Timer timer(false);
            size_t num_iterations = dimension[2] / depth;
            size_t leftover = dimension[2] % depth;
            if(leftover > 0) {
                num_iterations += 1;
            }
            for(int i=0;i<num_iterations;i++) {
                if(i == num_iterations - 1 && leftover > 0) {
                    std::vector<size_t> chunk_dimension = {dimension[0], dimension[1], leftover};
                    conf.setDims(chunk_dimension.begin(), chunk_dimension.end());
                }
                timer.start();
                fin.read(reinterpret_cast<char*>(&compressed_chunk_size), sizeof(int64_t));
                std::vector<char> buffer(compressed_chunk_size);
                fin.read(buffer.data(), compressed_chunk_size);
                total_read_time += timer.stop();
                timer.start();
                auto* decData = SZ_decompress<TYPE>(conf, buffer.data(), compressed_chunk_size);
                double decompress_time = timer.stop();
                total_decompress_time += decompress_time;
                timer.start();
                fout.write(reinterpret_cast<const char*>(decData), conf.num * sizeof(TYPE));
                delete[] decData;
                total_write_time += timer.stop();
                std::cout << "Chunk " << i << " compression completed! Time elasped: " << decompress_time
                          << " seconds, total time elapsed: " << total_timer.stop() << " seconds." << std::endl;
            }
            size_t file_size = (size_t)dimension[0] * (size_t)dimension[1] * (size_t)dimension[2] * sizeof(TYPE);
            size_t compressed_size = fin.tellg();
            std::cout << "Congratulations! Deompression completed! Total file size decompressed: " << file_size / 1024 / 1024 << " MB;\n"
                      << "total time elapsed: " << total_timer.stop() << "seconds" << std::endl;
            double cr = (double)file_size / (double) compressed_size;
            std::cout << "compression ratio:" << std::fixed << std::setprecision(2) << cr << ", compression_time: " << total_decompress_time << ", read time: " << total_read_time
                      << ", write time: " << total_write_time << std::endl;
        } else {
            conf.setDims(dimension.begin(), dimension.end());
            conf.absErrorBound = eb;
            size_t cmpSize;
            auto cmpData = SZ3::readfile<char>(input_file.c_str(), cmpSize);

            SZ3::Timer timer(true);
            auto *decData = SZ_decompress<TYPE>(conf, cmpData.get(), cmpSize);
            double decompress_time = timer.stop();
            SZ3::writefile<TYPE>(output_file.c_str(), decData, conf.num);
            delete[] decData;
            printf("compression ratio = %f\n", (double)conf.num * (double) sizeof(TYPE) * 1.0 / (double) cmpSize);
            printf("decompression time = %f seconds.\n", decompress_time);
            printf("decompressed file = %s\n", output_file.c_str());
        }
        return 0;
    }

    int decompress(int argc, char**argv) {
        int threads;
        std::vector<size_t> dimension;
        std::string input_file, output_file;
        float eb;
        bool isfloat64;
        std::string mode;
        size_t depth;
        parseCompressOptions(argc, argv, threads, input_file, output_file, dimension, eb, isfloat64, mode, depth);
        if(isfloat64) {
            return decompress_impl<double>(input_file, output_file, dimension, eb, mode, depth);
        }else{
            return decompress_impl<float>(input_file, output_file, dimension, eb, mode, depth);
        }
    }

    int test(int argc, char**argv) {
        std::cout << "The test functionality has not yet been implemented!" << std::endl;
        exit(EXIT_SUCCESS);
    }
}

int main(int argc, char** argv) {
    std::string general_helper_info = "Usage: sz3_split compress/decompress/test [options]\n"
                                      "use --help or -h to get the usage for each command\n";
    if(argc <= 1) {
        std::cout << general_helper_info << std::endl;
        return 0;
    }
    std::string command = argv[1];
    if(command == "compress") {
        sz3_split::compress(argc-1,argv+1);
    }else if(command == "decompress") {
        sz3_split::decompress(argc-1, argv+1);
    }else if(command == "test") {
        sz3_split::test(argc-1, argv+1);
    }else{
        std::cout << general_helper_info << std::endl;
    }
    return 0;
}