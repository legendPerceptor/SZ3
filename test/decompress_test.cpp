//
// Created by 刘远见 on 2021/1/14.
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
#include <type_traits>
#include <chrono>
#include "def.hpp"
//void myfunc(){
//    size_t my_size;
//    auto compressed_file = SZ::readfile<unsigned char>("/Users/hython/Development/globus/output/test.dat", my_size);
//    std::cout<<"read in size: "<<my_size<<std::endl;
//}
int main(int argc, char** argv) {
//    myfunc();
//    size_t num = 0;
    float eb = 0.1;
    // use Hurricane for testing 6113246
//    auto compressed = SZ::readfile<unsigned char>("/Users/hython/Development/globus/output/test.dat", num);
//    std::cout<<"read size: "<<num <<std::endl;

    size_t my_size;
    auto compressed_file = SZ::readfile<unsigned char>(argv[1], my_size);
    std::cout<<"read in size: "<<my_size<<std::endl;
//    uint dim = std::stoi(argv[2]);
//    std::vector<int> dims_;
//    for (int i=3;i<=2+dim;i++){
//        dims_.push_back(std::stoi(argv[i]));
//    }

    auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 2, 1>>(eb);
    auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 2>>(6, 0.1 * eb);
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 2>>> predictors_;
    predictors_.push_back(P_l);
    predictors_.push_back(P_reg);
//    auto cp = std::make_shared<SZ::ComposedPredictor<float, 3>>(predictors_);
    SZ::Config<float, 2> conf(eb, std::array<size_t, 2>{1800,3600});
    auto sz = SZ::SZ_General_Compressor<float, 2, SZ::ComposedPredictor<float, 2>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf,
            SZ::ComposedPredictor<float, 2>(predictors_),
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd()
    );
    std::unique_ptr<float[]> dec_data;
    dec_data.reset(sz.decompress(compressed_file.get(), my_size));
    SZ::writefile(argv[2], dec_data.get(), sz.getNumElements());
    return 0;
}
