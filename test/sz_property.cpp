//
// Created by apple on 2021/9/14.
//

#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/SZ_General_Compressor.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/LorenzoPredictor.hpp>
#include <predictor/RegressionPredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <utils/Iterator.hpp>
#include <utils/FileUtil.h>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <memory>
#include <chrono>
#include <tclap/CmdLine.h>
#include "qcat_dataAnalysis.h"
#include "csv.hpp"
#include <string>

int main(int argc, char**argv) {
    if(argc < 2)
    {
        printf("Usage: printDataProperty [dataType] tgtFilePath]\n");
        printf("Example: printDataProperty -f testfloat_8_8_128.dat\n");
        exit(0);
    }
    TCLAP::CmdLine cmd1("SZ3 Data Property", ' ', "0.1");
    TCLAP::ValueArg<std::string> inputFilePath("f","file", "The input data source file path",false,"","string");
    TCLAP::SwitchArg bigEndian("e", "bigendian", "Whether it's big endian", cmd1, false);
    TCLAP::SwitchArg debugArg("d", "debug", "Print additional information", cmd1, false);
    cmd1.add(inputFilePath);
    cmd1.parse(argc, argv);
    bool debug = debugArg.getValue();
    size_t num = 0;
    auto data = SZ::readfile<float>(inputFilePath.getValue().c_str(), num);
//    unsigned char* data  = NULL;

    int i = 0;
    if(debug) {
        printf("The first 10 values are: \n");
        for (i = 0; i < 10; i++)
            printf("%f ", data[i]);

        printf("....\n------------------------\n");
    }
    QCAT_DataProperty* property = computeProperty(QCAT_FLOAT, data.get(), num);
    std::stringstream ss;
    auto writer = csv::make_csv_writer(ss);
    writer << std::vector<std::string>({"size", "num", "min", "max", "valueRange","avgValue", "entropy", "zeromean_variance"});
    writer << std::vector<std::string>({std::to_string(property->totalByteSize),
                                        std::to_string(property->numOfElem),
                                        std::to_string(property->minValue),
                                        std::to_string(property->maxValue),
                                        std::to_string(property->valueRange),
                                        std::to_string(property->avgValue),
                                        std::to_string(property->entropy),
                                        std::to_string(property->zeromean_variance)});
    if(debug) {
        printProperty(property);
    }
    free(property);
    std::cout << ss.str() << std::endl;
    return 0;
}