// Created by yuanjian on 11/23/23.
// This tool is dedicated for extremely large file compression.
// The file cannot fit in the memory, and we need to split it to compress with multiple nodes or threads

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include "SZ3/api/sz.hpp"
#include "sz3_split.h"


int main(int argc, char** argv) {
    std::string general_helper_info = "Usage: sz3_split compress/decompress/test [options]\n"
                                      "use --help or -h to get the usage for each command\n";
    std::cout << "sz3_split program started!" << std::endl;
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