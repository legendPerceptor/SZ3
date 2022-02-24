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


template<uint N>
float lorenzo_test(float* data, std::array<size_t, N> dims, float eb=0.1){
    auto P_l = std::make_shared<SZ::LorenzoPredictor<float, N, 1>>(eb);
    int stride = 6;
    auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<float, N>>(data, std::begin(dims),std::end(dims), 1,0);
    auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<float, N>>(data, std::begin(dims),std::end(dims), stride,0);
    double avg_err = 0;
    int nble = 1;
    for(int i=0;i<N;i++){
        nble *= dims[i];
    }
    auto inter_begin = inter_block_range->begin();
    auto inter_end = inter_block_range->end();
    std::array<size_t, N> intra_block_dims;
//    double sum = 0;
    for (auto block = inter_begin; block != inter_end; ++block) {
        for (int i = 0; i < intra_block_dims.size(); i++) {
            size_t cur_index = block.get_local_index(i);
            size_t _dims = inter_block_range->get_dimensions(i);
            intra_block_dims[i] = (cur_index == _dims - 1 &&
                                   dims[i] - cur_index * stride < 6) ?
                                  dims[i] - cur_index * stride : 6;
        }

        intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
        intra_block_range->set_offsets(block.get_offset());
        intra_block_range->set_starting_position(block.get_local_index());
        auto intra_begin = intra_block_range->begin();
        auto intra_end = intra_block_range->end();
        for(auto iter = intra_begin; iter != intra_end; iter++) {
            float cur_err = P_l->estimate_error(iter);
//            sum += cur_err;
            avg_err += cur_err / (double) nble;
        }
    }
//    sum = sum / nble;
//    std::cout<<"AVG SUM:" << sum << "; avgerr" << avg_err <<std::endl;
    return avg_err;
}

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
    TCLAP::SwitchArg debugArg("b", "debug", "Print additional information", cmd1, false);
    TCLAP::ValueArg<std::string> dimensionArg("d", "dimension", "the dimention of data", false, "", "string");
    TCLAP::SwitchArg logcalculation("l", "log", "Whether use the log before anything", cmd1, false);

    cmd1.add(inputFilePath);
    cmd1.add(dimensionArg);
    cmd1.parse(argc, argv);
    bool debug = debugArg.getValue();
    size_t num = 0;
    bool use_log = logcalculation.getValue();
    auto data = SZ::readfile<float>(inputFilePath.getValue().c_str(), num);
//    std::cout<<"num: " << num <<std::endl;
    if(use_log) {
        for(int i=0;i<num;i++) {
            if(data[i] < 0) {
                data[i] = -log10(-data[i]);
            } else if(data[i]>0) {
                data[i] = log10(data[i]);
            }
        }
    }
//    unsigned char* data  = NULL;

    std::vector<size_t> dims;
    float avg_err;
    if(dimensionArg.isSet()) {
        std::string dimsString = dimensionArg.getValue();
        {
            std::stringstream ss;
            ss << dimsString;
            int tmp_dim;
            while(ss >> tmp_dim) {
                dims.push_back(tmp_dim);
            }
        }
        if(dims.size()==3){
            std::array<size_t, 3> _dims = {{dims[0], dims[1], dims[2]}};
            avg_err = lorenzo_test<3>(data.get(), _dims);
        }else if(dims.size()==2) {
            std::array<size_t, 2> _dims = {{dims[0], dims[1]}};
            avg_err = lorenzo_test<2>(data.get(), _dims);
        }else if(dims.size()==1) {
            std::array<size_t, 1> _dims = {{dims[0]}};
            avg_err = lorenzo_test<1>(data.get(), _dims);
        }
        if(debug) {
            std::cout << "Avg Error in Lorenzo: " << avg_err << std::endl;
        }
    }

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
    writer << std::vector<std::string>({"size", "num", "min", "max", "valueRange","avgValue", "entropy", "zeromean_variance", "avg_lorenzo"});
    writer << std::vector<std::string>({std::to_string(property->totalByteSize),
                                        std::to_string(property->numOfElem),
                                        std::to_string(property->minValue),
                                        std::to_string(property->maxValue),
                                        std::to_string(property->valueRange),
                                        std::to_string(property->avgValue),
                                        std::to_string(property->entropy),
                                        std::to_string(property->zeromean_variance),
                                        std::to_string(avg_err)});
    if(debug) {
        printProperty(property);
    }
    free(property);
    std::cout << ss.str() << std::endl;
    return 0;
}