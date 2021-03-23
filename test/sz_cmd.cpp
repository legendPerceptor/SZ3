//
// Created by apple on 2021/3/23.
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
#include <filesystem>
#include <tclap/CmdLine.h>

namespace fs = std::filesystem;
static void convert(float* data, int num) {
    float t;
    char *bytes;
    for(int i=0;i<num;i++) {
        bytes = (char*)&data[i];
        t = bytes[3];
        bytes[3] = bytes[0];
        bytes[0] = t;
        t = bytes[2];
        bytes[2] = bytes[1];
        bytes[1] = t;
    }
}
static void printData(float* data) {
    int z =1;
    printf("special %g\n", data[9098138]);
    for (int y=0; y< 100; y++) {
        for(int x=0;x<500;x++){
            printf("%g  ",data[x + 500*(y+500*z)]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {

    TCLAP::CmdLine cmd("SZ3 Compress/Decompress tests", ' ', "0.1");
    TCLAP::ValueArg<std::string> inputFilePath("i","input", "The input data source file path",true,"","string");
    TCLAP::ValueArg<std::string> outputFilePath("o", "output", "The compressed data output file path", true, "", "string");
//    TCLAP::MultiArg<int> dimension("d", "dimension", "the dimension of data",true,"multiInt");
    TCLAP::ValueArg<std::string> valueRange("r", "range", "The ranges with low, high, and error bound; '20 50 0.1; 50 80 0.01;'", true,"", "string");
    TCLAP::SwitchArg bigEndian("e", "bigEndian", "Whether it's big endian", cmd, false);
    TCLAP::SwitchArg use_bitmapArg("m","bitmap","Whether to use the bitmap", cmd, false);
    TCLAP::SwitchArg preserve_signArg("s", "preserveSign","Whether to preserve sign", cmd, false);
    TCLAP::SwitchArg hasBackgroundData("b", "backgroundData", "Whether there is background data", cmd, false);
    cmd.add(inputFilePath);
    cmd.add(outputFilePath);
//    cmd.add(dimension);
    cmd.add(valueRange);
    try {
        cmd.parse(argc, argv);
    }catch (TCLAP::ArgException &e) {
        std::cerr<< "error: "<<e.error() <<" for arg " << e.argId() << std::endl;
        exit(10);
    }
    std::string ranges = valueRange.getValue();
    auto ebs = std::vector<std::tuple<float, float, float>>();
    int start = 0;
    int end = ranges.find(';', start);
    float eb_min = 10000;
    while(end!=-1){
        std::string cur_rang = ranges.substr(start, end-start);
        std::stringstream ss(cur_rang);
        float low, high, eb;
        ss>>low>>high>>eb;
        ebs.emplace_back(low, high, eb);
        if(eb<eb_min){eb_min=eb;}
        start = end+1;
        end = ranges.find(';', start);
    }
    std::string inputFileStr = inputFilePath.getValue();
    size_t num = 0;
    auto data = SZ::readfile<float>(inputFileStr.c_str(), num);
    std::cout << "Read " << num << " elements\n";
    std::cout << "Original Size: " << num*sizeof(float) << std::endl;
    if(bigEndian.getValue()) { // convert big endian data
        convert(data.get(), num);
    }
//    auto dimensions = dimension.getValue();
//    size_t num_d=1;
//    for(auto iter=dimensions.begin();iter!=dimensions.end();iter++){
//        num_d*=(*iter);
//    }
//    if(num_d!=num) {
//        std::cerr<<"error: "<<"dimension doesn't match! "<<"Required Dimension: "
//                 <<dimensions.size()<<"; Actual parsed dimension: "<<std::endl;
//        exit(10);
//    }

    float eb =eb_min;
    float low_range = std::get<0>(*ebs.begin()), high_range = std::get<1>(*ebs.end());
    float bg = 1.0000000e+35;
    bool has_bg = hasBackgroundData.getValue();
    bool preserve_sign = preserve_signArg.getValue();
    bool use_bitmap = use_bitmapArg.getValue();
    const size_t DIM = 3;
    auto P_l = std::make_shared<SZ::LorenzoPredictor<float, DIM, 1>>(eb);
    auto P_reg = std::make_shared<SZ::RegressionPredictor<float, DIM>>(6, 0.1* eb);
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, DIM>>> predictors_;
    predictors_.push_back(P_l);
    predictors_.push_back(P_reg);

    SZ::Config<float, DIM> conf(eb, std::array<size_t, DIM>{512, 512, 512});
    auto sz = SZ::SZ_General_Compressor<float, DIM, SZ::ComposedPredictor<float, DIM>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf,
            SZ::ComposedPredictor<float, DIM>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd()
    );
    size_t compressed_size = 0;
    // Change to use Cpp-style timer
    // struct timespec start, end;
    std::chrono::time_point<std::chrono::system_clock> startTime, endTime;


    startTime = std::chrono::system_clock::now();
    std::unique_ptr<unsigned char[]> compressed;
    compressed.reset(sz.compress_withBG(data.get(), compressed_size, bg, low_range, high_range, use_bitmap, preserve_sign, has_bg));
//    compressed.reset(sz.compress(data.get(), compressed_size));
    endTime = std::chrono::system_clock::now();
    std::cout << "Compression time: "
              //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime).count()) / 1000000000 << "s"
              << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression Ratio = " << num * sizeof (float)/ (float)compressed_size << std::endl;
    SZ::writefile(outputFilePath.getValue().c_str(), compressed.get(), compressed_size);

    startTime = std::chrono::system_clock::now();
//    err = clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<float[]> dec_data;
    dec_data.reset(sz.decompress_withBG(compressed.get(), compressed_size, bg, low_range, high_range, use_bitmap, preserve_sign, has_bg));
//    dec_data.reset(sz.decompress(compressed.get(), compressed_size));
//    err = clock_gettime(CLOCK_REALTIME, &end);
    endTime = std::chrono::system_clock::now();
    std::cout << "Decompression time: "
              //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime).count()) / 1000000000 << "s"
              << std::endl;
//    printData(dec_data.get());
    auto dataV = SZ::readfile<float>(inputFileStr.c_str(), num);
    if(bigEndian.getValue()) {
        convert(dataV.get(), num);
    }
    float max_err = 0;
    for (int i = 0; i < num; i++) {
//        max_err = std::max(max_err, std::abs(data[i] - dec_data[i]));
        if(dataV[i] - dec_data[i] > max_err || dataV[i] - dec_data[i] < -max_err) {
            max_err = (dataV[i] > dec_data[i]) ? dataV[i] - dec_data[i] : dec_data[i] - dataV[i];
//            printf("NOTG data: %g, dec_data: %g, i:%d\n", dataV[i], dec_data[i], i);
        }
    }
    std::cout << "Max error = " << max_err << std::endl;
    return 0;
}
