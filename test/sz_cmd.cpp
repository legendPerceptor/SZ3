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
#include <chrono>
#include <tclap/CmdLine.h>
#include <fstream>

//namespace fs = std::filesystem;
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
    bool FromFile = false;
    std::vector<std::string> myargv;
    if(argc<=2) {
        FromFile = true;
        std::fstream f;
        if(argc==2) {
            f.open(argv[1], std::ios::in);
        }else if(argc==1) {
            f.open("sz3.config", std::ios::in);
        }
        if(f.fail()){
            std::cerr<< "Please make sure the config file exists!\n "
            << "Use commandline arguments or config file\n"
            << "The default config file is sz3.config"<<std::endl;
            exit(5);
        }
        std::stringstream strStream;
        strStream << f.rdbuf();
        std::string tmp;
        while(strStream >> tmp){
            if(tmp[0]=='\"'){
                std::string t=tmp.substr(1, tmp.size()-1);
                do {
                    strStream >> tmp;
                    size_t a = tmp.find('\"');
                    if(a !=std::string::npos){
                        t += ' ' + tmp.substr(0, a);
                        break;
                    }
                    t+=' ' + tmp;
                } while(strStream);
                myargv.push_back(t);
            }else {
                myargv.push_back(tmp);
            }
        }
    }
    TCLAP::CmdLine cmd("SZ3 Compress/Decompress tests", ' ', "0.1");
    TCLAP::ValueArg<std::string> inputFilePath("i","input", "The input data source file path",false,"","string");
    TCLAP::ValueArg<std::string> outputFilePath("c", "compress", "The compressed data output file path", true, "", "string");
//    TCLAP::MultiArg<int> dimension("d", "dimension", "the dimension of data",true,"multiInt");
    TCLAP::ValueArg<std::string> dimensionArg("d", "dimension", "the dimention of data", true, "", "string");
    TCLAP::ValueArg<std::string> valueRange("r", "range", "The ranges with low, high, and error bound; '20 50 0.1; 50 80 0.01;'", true,"", "string");
    TCLAP::SwitchArg bigEndian("e", "bigEndian", "Whether it's big endian", cmd, false);
    TCLAP::SwitchArg use_bitmapArg("p","bitmap","Whether to use the bitmap", cmd, false);
    TCLAP::SwitchArg preserve_signArg("s", "preserveSign","Whether to preserve sign", cmd, false);
    TCLAP::SwitchArg hasBackgroundData("b", "backgroundData", "Whether there is background data", cmd, false);
    TCLAP::ValueArg<std::string> decFilePath("q", "decFile", "The decompressed data file", true, "", "string");
    TCLAP::ValueArg<std::string> logFilePath("l","logFile","The log file path", true, "", "string");
    TCLAP::SwitchArg fall_back("f", "fallback", "Whether to use old SZ3 compressor", cmd, false);
    TCLAP::ValueArg<std::string> modeArg("m", "mode", "The mode of the program (test, compress, decompress)", false, "test", "string");
    cmd.add(inputFilePath);
    cmd.add(outputFilePath);
    cmd.add(dimensionArg);
    cmd.add(valueRange);
    cmd.add(decFilePath);
    cmd.add(logFilePath);
    cmd.add(modeArg);
    try {
        if(FromFile){
            char** arr = new char*[myargv.size()+1];
            arr[0] = new char[strlen(argv[0])+1];
            strcpy(arr[0], argv[0]);
            for(size_t i=0; i < myargv.size(); i++) {
                arr[i+1] = new char[myargv[i].size()+1];
                strcpy(arr[i+1], myargv[i].c_str());
            }
            cmd.parse(myargv.size()+1, arr);
        }else {
            cmd.parse(argc, argv);
        }
    }catch (TCLAP::ArgException &e) {
        std::cerr<< "error: "<<e.error() <<" for arg " << e.argId() << std::endl;
        exit(10);
    }
    std::string mode = modeArg.getValue();
    if(mode=="test" || mode=="compress"){
        if(inputFilePath.getValue()=="") {
            std::cerr << "ERROR: compress/test mode; must provide the original data file" << std::endl;
            exit(10);
        }
    }
    std::string ranges = valueRange.getValue();
    auto ebs = std::vector<SZ::RangeTuple<float>>();
    int start = 0;
    int end = ranges.find(';', start);
    float eb_min = 10000;
    bool fallback = fall_back.getValue();
    while(end!=-1){
        std::string cur_rang = ranges.substr(start, end-start);
        std::stringstream ss(cur_rang);
        float low, high, eb;
        ss>>low>>high>>eb;
        ebs.push_back(SZ::RangeTuple(low, high, eb));
        if(eb<eb_min){eb_min=eb;}
        start = end+1;
        end = ranges.find(';', start);
    }

    std::string dimsString = dimensionArg.getValue();
    std::vector<size_t> dims;
    {
        std::stringstream ss;
        ss << dimsString;
        int tmp_dim;
        while(ss >> tmp_dim) {
            dims.push_back(tmp_dim);
        }
    }
    float eb =eb_min;
    float low_range = (*ebs.begin()).low, high_range = ebs[ebs.size()-1].high;
    float bg = 1.0000000e+35;
    bool has_bg = hasBackgroundData.getValue();
    bool preserve_sign = preserve_signArg.getValue();
    bool use_bitmap = use_bitmapArg.getValue();
    const size_t DIM = 3;



    SZ::Compressor<float> *sz, *sz_old;
    if(dims.size()==3) {
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 3, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 3>>(6, 0.0001);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, DIM>>> predictors_;
        predictors_.push_back(P_l);
        predictors_.push_back(P_reg);
        SZ::Config<float, 3> conf(eb, std::array<size_t, 3>{  dims[0], dims[1], dims[2]});
        sz = new SZ::SZ_General_Compressor<float, 3, SZ::ComposedPredictor<float, 3>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf,
                SZ::ComposedPredictor<float, 3>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
                SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd()
        );
        sz_old = new SZ::SZ_General_Compressor<float, 3, SZ::ComposedPredictor<float, 3>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf,
                SZ::ComposedPredictor<float, 3>(predictors_),
                SZ::LinearQuantizer<float>(eb),
//            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd()
        );
    } else if(dims.size()==2) {
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 2, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 2>>(6, 0.0001);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 2>>> predictors_;
        predictors_.push_back(P_l);
        predictors_.push_back(P_reg);
        SZ::Config<float, 2> conf(eb, std::array<size_t, 2>{  dims[0], dims[1]});
        sz = new SZ::SZ_General_Compressor<float, 2, SZ::ComposedPredictor<float, 2>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf,
                SZ::ComposedPredictor<float, 2>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
                SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd()
        );
        sz_old = new SZ::SZ_General_Compressor<float, 2, SZ::ComposedPredictor<float, 2>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf,
                SZ::ComposedPredictor<float, 2>(predictors_),
                SZ::LinearQuantizer<float>(eb),
//            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd()
        );
    } else if (dims.size()==1){
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 1, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 1>>(6, 0.0001);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 1>>> predictors_;
        predictors_.push_back(P_l);
        predictors_.push_back(P_reg);
        SZ::Config<float, 1> conf(eb, std::array<size_t, 1>{  dims[0]});
        sz = new SZ::SZ_General_Compressor<float, 1, SZ::ComposedPredictor<float, 1>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf,
                SZ::ComposedPredictor<float, 1>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
                SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd()
        );
        sz_old = new SZ::SZ_General_Compressor<float, 1, SZ::ComposedPredictor<float, 1>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf,
                SZ::ComposedPredictor<float, 1>(predictors_),
                SZ::LinearQuantizer<float>(eb),
//            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd()
        );
    }

    size_t compressed_size = 0;
    std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
    std::ofstream fs(logFilePath.getValue().c_str(), std::ios_base::app);
    size_t num = 0;
    std::string inputFileStr = inputFilePath.getValue();
    if(mode == "test" || mode == "compress") {
        auto data = SZ::readfile<float>(inputFileStr.c_str(), num);
        std::cout << "Read " << num << " elements\n";
        std::cout << "Original Size: " << num * sizeof(float) << std::endl;
        if (bigEndian.getValue()) { // convert big endian data
            convert(data.get(), num);
        }
//        std::cout<<"special: "<< data[5324929]<<std::endl;
//        auto quantizer = SZ::MultipleErrorBoundsQuantizer<float>(ebs);
//        float dp= data[5324929];
//        int tmp_quant = quantizer.quantize_and_overwrite(dp, 0.7888);
//        float dp2 = quantizer.recover(0.7888, tmp_quant);
//        printf("tests!!!!!!");
        startTime = std::chrono::system_clock::now();
        std::unique_ptr<unsigned char[]> compressed;
        if (fallback) {
            compressed.reset(sz_old->compress(data.get(), compressed_size));
        } else if (!has_bg) {
            compressed.reset(sz->compress(data.get(), compressed_size));
        } else {
            compressed.reset(
                    sz->compress_withBG(data.get(), compressed_size, bg, low_range, high_range, use_bitmap,
                                       preserve_sign,
                                       has_bg));
        }
        endTime = std::chrono::system_clock::now();
        std::cout << "Compression time: "
                  //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
                  << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) /
                     1000000000 << "s"
                  << std::endl;
        std::cout << "Compressed size = " << compressed_size << std::endl;
        std::cout << "Compression Ratio = " << num * sizeof(float) / (float) compressed_size << std::endl;
        SZ::writefile(outputFilePath.getValue().c_str(), compressed.get(), compressed_size);
        fs << outputFilePath.getValue() << " Compression Time: "
           << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) / 1000000000
           << "s;"
           << "Compression size: " << compressed_size << "; Compression ratio: "
           << num * sizeof(float) / (float) compressed_size << "; ";
    }
    if(mode=="test" || mode=="decompress") {
        auto compressed = SZ::readfile<unsigned char>(outputFilePath.getValue().c_str(), compressed_size);
        startTime = std::chrono::system_clock::now();
        std::unique_ptr<float[]> dec_data;
        if (fallback) {
            dec_data.reset(sz_old->decompress(compressed.get(), compressed_size));
        } else if (!has_bg) {
            dec_data.reset(sz->decompress(compressed.get(), compressed_size));
        } else {
            dec_data.reset(
                    sz->decompress_withBG(compressed.get(), compressed_size, bg, low_range, high_range, use_bitmap,
                                         preserve_sign, has_bg));
        }
        SZ::writefile(decFilePath.getValue().c_str(), dec_data.get(), num);
//    err = clock_gettime(CLOCK_REALTIME, &end);
        endTime = std::chrono::system_clock::now();
        std::cout << "Decompression time: "
                  //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
                  << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) /
                     1000000000 << "s"
                  << std::endl;
        fs << "Decompression Time: "
           << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) / 1000000000
           << "s; ";

        if(mode=="test") {
            auto dataV = SZ::readfile<float>(inputFileStr.c_str(), num);
            if (bigEndian.getValue()) {
                convert(dataV.get(), num);
            }
            float max_err = 0;
            std::cout << "Low: " << low_range << ", high: " << high_range << std::endl;
            for (int i = 0; i < num; i++) {
                if (!has_bg) {
                    dataV[i] = fmin(fmax(low_range, dataV[i]), high_range);
                }
                if (dataV[i] - dec_data[i] > max_err || dataV[i] - dec_data[i] < -max_err) {
                    max_err = (dataV[i] > dec_data[i]) ? dataV[i] - dec_data[i] : dec_data[i] - dataV[i];
                }
            }
            std::cout << "Max error = " << max_err << std::endl;
            fs << "Max error = " << max_err << std::endl;
        }
    }
    fs.close();
    return 0;
}
