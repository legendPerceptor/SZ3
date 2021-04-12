//
// Created by apple on 2021/4/11.
//
#include "H5Z_SZ3.h"
#include <iostream>
#include <cmath>
#include <memory>
#include <chrono>
#include <tclap/CmdLine.h>
#include <fstream>
#include <string>

#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/SZ_General_Compressor.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/LorenzoPredictor.hpp>
#include <predictor/RegressionPredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <utils/Iterator.hpp>
#include <utils/FileUtil.h>

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    int dimension;
    if(r1==0)
    {
        dimension = 0;
    }
    else if(r2==0)
    {
        dimension = 1;
    }
    else if(r3==0)
    {
        dimension = 2;
    }
    else if(r4==0)
    {
        dimension = 3;
    }
    else if(r5==0)
    {
        dimension = 4;
    }
    else
    {
        dimension = 5;
    }
    return dimension;
}

 void longToBytes_bigEndian(unsigned char *b, unsigned long num)
{
    b[0] = (unsigned char)(num>>56);
    b[1] = (unsigned char)(num>>48);
    b[2] = (unsigned char)(num>>40);
    b[3] = (unsigned char)(num>>32);
    b[4] = (unsigned char)(num>>24);
    b[5] = (unsigned char)(num>>16);
    b[6] = (unsigned char)(num>>8);
    b[7] = (unsigned char)(num);
//	if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_8bytes(*b);
}

 int bytesToInt_bigEndian(unsigned char* bytes)
{
    int temp = 0;
    int res = 0;

    res <<= 8;
    temp = bytes[0] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[1] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[2] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[3] & 0xff;
    res |= temp;

    return res;
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    size_t dataLength;
    if(r1==0)
    {
        dataLength = 0;
    }
    else if(r2==0)
    {
        dataLength = r1;
    }
    else if(r3==0)
    {
        dataLength = r1*r2;
    }
    else if(r4==0)
    {
        dataLength = r1*r2*r3;
    }
    else if(r5==0)
    {
        dataLength = r1*r2*r3*r4;
    }
    else
    {
        dataLength = r1*r2*r3*r4*r5;
    }
    return dataLength;
}

 void intToBytes_bigEndian(unsigned char *b, unsigned int num)
{
    b[0] = (unsigned char)(num >> 24);
    b[1] = (unsigned char)(num >> 16);
    b[2] = (unsigned char)(num >> 8);
    b[3] = (unsigned char)(num);

    //note: num >> xxx already considered endian_type...
//if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_4bytes(*b); //change to BIG_ENDIAN_DATA
}

 long bytesToLong_bigEndian(unsigned char* b) {
    long temp = 0;
    long res = 0;

    res <<= 8;
    temp = b[0] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[1] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[2] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[3] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[4] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[5] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[6] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[7] & 0xff;
    res |= temp;

    return res;
}


void H5Z_SZ3_Init(SZ::Compressor<float> *& sz, SZ::Compressor<float>*& sz_old, SZ3_config_params & sz3conf, int r3, int r2, int r1) {
    std::vector<std::string> myargv;
    bool FromFile;
    FromFile = true;
    std::fstream f;
    f.open("sz3.config", std::ios::in);
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

    TCLAP::CmdLine cmd("SZ3 Compress/Decompress tests", ' ', "0.1");
//    TCLAP::ValueArg<std::string> inputFilePath("i","input", "The input data source file path",false,"","string");
//    TCLAP::ValueArg<std::string> outputFilePath("c", "compress", "The compressed data output file path", true, "", "string");
//    TCLAP::MultiArg<int> dimension("d", "dimension", "the dimension of data",true,"multiInt");
//    TCLAP::ValueArg<std::string> dimensionArg("d", "dimension", "the dimention of data", true, "", "string");
    TCLAP::ValueArg<std::string> valueRange("r", "range", "The ranges with low, high, and error bound; '20 50 0.1; 50 80 0.01;'", true,"", "string");
    TCLAP::SwitchArg bigEndian("e", "bigEndian", "Whether it's big endian", cmd, false);
    TCLAP::SwitchArg use_bitmapArg("p","bitmap","Whether to use the bitmap", cmd, false);
    TCLAP::SwitchArg preserve_signArg("s", "preserveSign","Whether to preserve sign", cmd, false);
    TCLAP::ValueArg<float> hasBackgroundData("b", "backgroundData", "Whether there is background data", false,0, "float");
//    TCLAP::ValueArg<std::string> decFilePath("q", "decFile", "The decompressed data file", true, "", "string");
    TCLAP::ValueArg<std::string> logFilePath("l","logFile","The log file path", true, "", "string");
    TCLAP::SwitchArg fall_back("f", "fallback", "Whether to use old SZ3 compressor", cmd, false);
//    TCLAP::ValueArg<std::string> modeArg("m", "mode", "The mode of the program (test, compress, decompress)", false, "test", "string");
//    cmd.add(inputFilePath);
//    cmd.add(outputFilePath);
//    cmd.add(dimensionArg);
    cmd.add(valueRange);
    cmd.add(hasBackgroundData);
//    cmd.add(decFilePath);
    cmd.add(logFilePath);
//    cmd.add(modeArg);
    try {
        char** arr = new char*[myargv.size()+1];
        arr[0] = new char[4];
        strcpy(arr[0], "n\0");
        for(size_t i=0; i < myargv.size(); i++) {
            arr[i+1] = new char[myargv[i].size()+1];
            strcpy(arr[i+1], myargv[i].c_str());
        }
        cmd.parse(myargv.size()+1, arr);
        for(size_t i=0; i<myargv.size()+1; i++) {
            delete arr[i];
        }
        delete[]arr;
    }catch (TCLAP::ArgException &e) {
        std::cerr<< "error: "<<e.error() <<" for arg " << e.argId() << std::endl;
        exit(10);
    }
//    std::string mode = modeArg.getValue();
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

//    std::string dimsString = dimensionArg.getValue();
    std::vector<size_t> dims;
//    {
//        std::stringstream ss;
//        ss << dimsString;
//        int tmp_dim;
//        while(ss >> tmp_dim) {
//            dims.push_back(tmp_dim);
//        }
//    }
    printf("The dimensions are r1=%d, r2=%d, r3=%d;\n", r1, r2, r3);
    if(r1!=0){dims.push_back(r1);}
    if(r2!=0){dims.push_back(r2);}
    if(r3!=0){dims.push_back(r3);}
    float eb =eb_min;
    float low_range = (*ebs.begin()).low, high_range = ebs[ebs.size()-1].high;
    float bg = hasBackgroundData.getValue();
//    bool has_bg = hasBackgroundData.getValue();
    bool preserve_sign = preserve_signArg.getValue();
    bool use_bitmap = use_bitmapArg.getValue();
    sz3conf.low = low_range;
    sz3conf.high = high_range;
    sz3conf.has_bg = bg!=0;
    sz3conf.preserve_sign = preserve_sign;
    sz3conf.usebitmap = use_bitmap;
    sz3conf.bg_value = bg;
    sz3conf.fallback = fallback;
    if(dims.size()==3) {
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 3, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 3>>(6, 0.0001);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 3>>> predictors_;
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
}