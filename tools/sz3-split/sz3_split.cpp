// Created by yuanjian on 11/23/23.
// This tool is dedicated for extremely large file compression.
// The file cannot fit in the memory, and we need to split it to compress with multiple nodes or threads
#include <getopt.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "SZ3/api/sz.hpp"

namespace sz3_split {

    void
    parseCompressOptions(int argc, char **argv, int &threads, std::string &raw_file, std::string &output_file,
                         std::vector<int> &data_dimension, float& eb) {
        optind = 1;
        const char *opt_index = "ht:i:d:e:o:";
        struct option opts[] = {
                {"threads",    required_argument, nullptr, 't'},
                {"help",       no_argument,       nullptr, 'h'},
                {"input",      required_argument, nullptr, 'i'},
                {"output", required_argument, nullptr, 'o'},
                {"errorbound", required_argument, nullptr, 'e'},
                {"dimension",  required_argument, nullptr, 'd'}
        };
        std::string compress_helper_info = "Usage: sz3_split compress [options]\n"
                                           "options:  --threads/-t     INT   number of threads, default is 1\n"
                                           "          --input/-i       STR   the RAW file/compressed file\n"
                                           "          --output/-o      STR   the compressed file/decompressed file location\n"
                                           "          --help/-h              print this help information\n"
                                           "          --dimension/-d   STR   the data dimension of the file, e.g., 256 256 512\n"
                                           "          --errorbound/-e  FLOAT the error bound to use in compression\n";
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

    int compress(int argc, char **argv) {
        int threads;
        std::vector<int> dimension;
        std::string input_file, output_file;
        float eb;
        parseCompressOptions(argc, argv, threads, input_file, output_file, dimension, eb);
        SZ3::Config conf = defaultConfig(); // 300 is the fastest dimension
        conf.setDims(dimension.begin(), dimension.end());
        conf.absErrorBound = eb;
        auto *data = new float[conf.num];
        SZ3::readfile<float>(input_file.c_str(), conf.num, data);
        size_t outSize;
        SZ3::Timer timer(true);
        char *compressedData = SZ_compress(conf, data, outSize);
        double compress_time = timer.stop();
        std::cout << "Compression completed! Time elasped: " << compress_time << std::endl;
        SZ3::writefile(output_file.c_str(), compressedData, outSize);
        printf("compression ratio = %.2f \n", (double)conf.num * 1.0 * sizeof(float) / (double)outSize);
        printf("compression time = %f\n", compress_time);
        printf("compressed data file = %s\n", output_file.c_str());
    }

    int decompress(int argc, char**argv) {
        int threads;
        std::vector<int> dimension;
        std::string input_file, output_file;
        float eb;
        parseCompressOptions(argc, argv, threads, input_file, output_file, dimension, eb);
        SZ3::Config conf = defaultConfig(); // 300 is the fastest dimension
        conf.setDims(dimension.begin(), dimension.end());
        conf.absErrorBound = eb;
        size_t cmpSize;
        auto cmpData = SZ3::readfile<char>(input_file.c_str(), cmpSize);

        SZ3::Timer timer(true);
        auto *decData = SZ_decompress<float>(conf, cmpData.get(), cmpSize);
        double decompress_time = timer.stop();
        SZ3::writefile<float>(output_file.c_str(), decData, conf.num);
        printf("compression ratio = %f\n", conf.num * (double)sizeof(float) * 1.0 / (double)cmpSize);
        printf("decompression time = %f seconds.\n", decompress_time);
        printf("decompressed file = %s\n", output_file.c_str());
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