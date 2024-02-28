//
// Created by apple on 2021/7/28.
//

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
#include "csv.hpp"

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

const std::string help_info = R"(
Additional Help Information for sz_region

Both range and region based compression are included in the sz_region executable.

Examples:

Set error bound 0.2 for range [-1, 0.98), error bound 0.001 for range [0.98, 1.02),
and error bound 0.01 for range [1.02, 10) for the Katrina dataset.

sz_region -i katrina.bin -c katrina.bin.sz -q katrina.bin.out \
-r "-1 0.98 0.2; 0.98 1.02 0.001; 1.02 10 0.01;" -d "417642 162;" -m test

For the QMCPACK dataset of dimension (33120, 69, 69), we set error bound 0.01 for range [-17,-8),
error bound 0.1 for range [-8, -5), and error bound 1 for range [-5, 17).

sz_region -i QMCPACK.f32 -c QMCPACK.sz -q QMCPACK.dp \
-r "-17 -8 0.01; -8 -5 0.1; -5 17 1;" -d "33120 69 69;" -m test

For region based compression, the format is as follows: default_error_bound>x_lower y_lower z_lower:x_length y_length z_length error_bound.

sz_region -i CLDHGH_1_1800_3600.dat -c CLDHGH.sz -q CLDHGH.dp \
--region "0.001>5 5 5: 8 10 12 0.1; 20 20 20: 90 90 90 0.01;" -d "1800 3600;" -l log.txt -m test
)";

int main(int argc, char **argv) {
    bool FromFile = false;
    std::vector<std::string> myargv;

    TCLAP::CmdLine cmd(help_info, ' ', "0.1");
    TCLAP::ValueArg<std::string> inputFilePath("i","input", "The input data source file path",false,"","string");
    TCLAP::ValueArg<std::string> outputFilePath("c", "compress", "The compressed data output file path", true, "", "string");
//    TCLAP::MultiArg<int> dimension("d", "dimension", "the dimension of data",true,"multiInt");
    TCLAP::ValueArg<std::string> dimensionArg("d", "dimension", "the dimention of data", true, "", "string");
    TCLAP::ValueArg<std::string> valueRange("r", "range", "The ranges with low, high, and error bound; '20 50 0.1; 50 80 0.01;'", false,"", "string");
    TCLAP::ValueArg<std::string> regions("w", "region", "The regions; '0.1> 5 5 5: 8 10 12 0.1; 200 200 200: 60 60 60 0.01;'", false, "", "string");
    TCLAP::SwitchArg bigEndian("e", "bigEndian", "Whether it's big endian", cmd, false);
    TCLAP::SwitchArg use_bitmapArg("p","bitmap","Whether to use the bitmap", cmd, false);
    TCLAP::SwitchArg preserve_signArg("s", "preserveSign","Whether to preserve sign", cmd, false);
    TCLAP::ValueArg<std::string> hasBackgroundData("b", "background", "set the background data to be this value", false, "", "string");
    TCLAP::ValueArg<std::string> decFilePath("q", "decFile", "The decompressed data file", true, "", "string");
    TCLAP::SwitchArg logcalculation("l", "log", "Whether use the log before anything", cmd, false);
    TCLAP::SwitchArg fall_back("f", "fallback", "Whether to use old SZ3 compressor", cmd, false);
    TCLAP::ValueArg<std::string> from_file_arg("", "cmd_from_file", "The commands should be read from the specified file", false, "", "string");
    TCLAP::ValueArg<std::string> modeArg("m", "mode", "The mode of the program (test, compress, decompress)", false, "test", "string");
    TCLAP::SwitchArg random_accessArg("a", "random", "Independent access between regions", cmd, false);
    cmd.add(inputFilePath);
    cmd.add(outputFilePath);
    cmd.add(dimensionArg);
    cmd.add(valueRange);
    cmd.add(decFilePath);
//    cmd.add(logcalculation);
    cmd.add(hasBackgroundData);
    cmd.add(modeArg);
    cmd.add(regions);
    try {
        if(from_file_arg.isSet()){

            std::fstream f;
            f.open(from_file_arg.getValue(), std::ios::in);
            if(f.fail()){
                std::cerr<< "Please make sure the command config file exists!\n "
                        << "Use commandline arguments or config file\n"<<std::endl;
                std::cout << help_info << std::endl;
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

    bool use_log = logcalculation.getValue();
    std::string mode = modeArg.getValue();
    if(mode=="test" || mode=="compress"){
        if(inputFilePath.getValue()=="") {
            std::cerr << "ERROR: compress/test mode; must provide the original data file" << std::endl;
            exit(10);
        }
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

    auto ebs = std::vector<SZ::RangeTuple<float>>();
    auto region_ebs = std::vector<SZ::RegionTuple<float>>();
    float eb_min = 10000;
    float default_eb;
    bool fallback = fall_back.getValue();
    bool random_access = random_accessArg.getValue();
    float low_range, high_range;
    if(valueRange.isSet() && !regions.isSet()) {
        std::string ranges = valueRange.getValue();
        int start = 0;
        int end = ranges.find(';', start);

        while (end != -1) {
            std::string cur_rang = ranges.substr(start, end - start);
            std::stringstream ss(cur_rang);
            float low, high, eb;
            ss >> low >> high >> eb;
            ebs.push_back(SZ::RangeTuple(low, high, eb));
            if (eb < eb_min) { eb_min = eb; }
            start = end + 1;
            end = ranges.find(';', start);
        }
        low_range = (*ebs.begin()).low;
        high_range = ebs[ebs.size()-1].high;
    } else if(regions.isSet() && ! valueRange.isSet()) {
        std::string regionStr = regions.getValue();
        int start = 0;
        int end = regionStr.find('>', start);
        if(end==std::string::npos) {
            std::cerr<<"The region requires a default error bound!" << std::endl;
            exit(1);
        }
        default_eb = std::stof(regionStr.substr(start, end - start));
        start = end + 1;
        while(end!=std::string::npos) {
            int *start_position = new int[dims.size()];
            int *region_length = new int[dims.size()];
            end = regionStr.find(':', start);
            if(end == std::string::npos) {
                break;
            }
            std::stringstream ss(regionStr.substr(start, end - start));
            for (int i = 0; i < dims.size(); i++) {
                ss >> start_position[i];
            }
            start = end + 1;
            end = regionStr.find(';', start);
            std::stringstream ss2(regionStr.substr(start, end - start));
            for (int i = 0; i < dims.size(); i++) {
                ss2 >> region_length[i];
            }
            float eb;
            ss2 >> eb;
            if(eb < eb_min) {
                eb_min = eb;
            }
            region_ebs.emplace_back(start_position, region_length, eb);
            start = end + 1;
        }
    }


    float eb;
    if(regions.isSet()) {
        eb = default_eb;
    } else {
        eb = eb_min;
    }

    float bg = -99999;
    bool has_bg = hasBackgroundData.isSet() && mode != "bg_pre";
    if(hasBackgroundData.isSet()) {
        bg = std::stof(hasBackgroundData.getValue());
    }
    bool preserve_sign = preserve_signArg.getValue();
    bool use_bitmap = use_bitmapArg.getValue();
//    const size_t DIM = 3;

    size_t compressed_size = 0;
    std::chrono::time_point<std::chrono::system_clock> startTime, endTime;

    size_t num = 0;
    std::string inputFileStr = inputFilePath.getValue();

    if (mode == "bg_pre") {
        auto data = SZ::readfile<float>(inputFileStr.c_str(), num);
        startTime = std::chrono::system_clock::now();
        std::cout << "Read " << num << " elements\n";
        std::cout << "Original Size: " << num * sizeof(float) << std::endl;
        if (bigEndian.getValue()) { // convert big endian data
            std::cout<<"converting data from bigEndian..."<<std::endl;
            convert(data.get(), num);
        }
        for(int i=0;i<num;i++) {
            if(data[i] - bg < std::numeric_limits<float>::epsilon()){
                data[i]=0;
            }
        }
        SZ::writefile("tmp.bg", data.get(), num);
        inputFileStr = "tmp.bg";
        endTime = std::chrono::system_clock::now();
        std::cout << "Preprocessing Time: " << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) / 1000000000 << std::endl;
        mode = "test";
    }


    SZ::Compressor<float> *sz, *sz_old, *sz_region;
    if(dims.size()==3) {
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 3, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 3>>(6, 0.1*eb);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 3>>> predictors_;
        predictors_.push_back(P_l);
        predictors_.push_back(P_reg);
        SZ::Config<float, 3> conf(eb, std::array<size_t, 3>{  dims[0], dims[1], dims[2]});
        if(fallback) {
            sz_old = new SZ::SZ_General_Compressor<float, 3, SZ::ComposedPredictor<float, 3>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 3>(predictors_),
                    SZ::LinearQuantizer<float>(eb),
//            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        } else if(valueRange.isSet()) {
            sz = new SZ::SZ_General_Compressor<float, 3, SZ::ComposedPredictor<float, 3>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 3>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
                    SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        }else if(regions.isSet()) {
            sz_region = new SZ::SZ_General_Compressor<float, 3, SZ::ComposedPredictor<float, 3>, SZ::RegionBasedQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 3>(predictors_),
                    SZ::RegionBasedQuantizer<float>(region_ebs, 3, default_eb),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        }
    } else if(dims.size()==2) {
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 2, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 2>>(6, 0.1*eb);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 2>>> predictors_;
        predictors_.push_back(P_l);
        predictors_.push_back(P_reg);
        SZ::Config<float, 2> conf(eb, std::array<size_t, 2>{  dims[0], dims[1]});
        if(fallback) {
            sz_old = new SZ::SZ_General_Compressor<float, 2, SZ::ComposedPredictor<float, 2>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 2>(predictors_),
                    SZ::LinearQuantizer<float>(eb),
//            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        } else if(valueRange.isSet()) {
            sz = new SZ::SZ_General_Compressor<float, 2, SZ::ComposedPredictor<float, 2>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 2>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
                    SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        } else if (regions.isSet()) {
            sz_region = new SZ::SZ_General_Compressor<float, 2, SZ::ComposedPredictor<float, 2>, SZ::RegionBasedQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 2>(predictors_),
                    SZ::RegionBasedQuantizer<float>(region_ebs, 2, default_eb),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        }
    } else if (dims.size()==1){
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 1, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 1>>(6, 0.1*eb);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 1>>> predictors_;
        predictors_.push_back(P_l);
        predictors_.push_back(P_reg);
        SZ::Config<float, 1> conf(eb, std::array<size_t, 1>{  dims[0]});
        if(fallback) {
            sz_old = new SZ::SZ_General_Compressor<float, 1, SZ::ComposedPredictor<float, 1>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 1>(predictors_),
                    SZ::LinearQuantizer<float>(eb),
//            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        } else if(valueRange.isSet()) {
            sz = new SZ::SZ_General_Compressor<float, 1, SZ::ComposedPredictor<float, 1>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 1>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
                    SZ::MultipleErrorBoundsQuantizer<float>(ebs),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        } else if (regions.isSet()) {
            sz_region = new SZ::SZ_General_Compressor<float, 1, SZ::ComposedPredictor<float, 1>, SZ::RegionBasedQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                    conf,
                    SZ::ComposedPredictor<float, 1>(predictors_),
                    SZ::RegionBasedQuantizer<float>(region_ebs, 1, default_eb),
                    SZ::HuffmanEncoder<int>(),
                    SZ::Lossless_zstd()
            );
        }
    }


    std::vector<std::string> csv_result;
    std::shared_ptr<bool[]> negative;
    if(mode == "test" || mode == "compress") {
        auto data = SZ::readfile<float>(inputFileStr.c_str(), num);

        std::cout << "Read " << num << " elements\n";
        std::cout << "Original Size: " << num * sizeof(float) << std::endl;
        if (bigEndian.getValue()) { // convert big endian data
            convert(data.get(), num);
        }
        if(use_log) {
            negative = std::shared_ptr<bool[]>(new bool[num]);
            for(int i=0;i<num;i++) {
                if(data[i] < 0) {
                    negative[i] = true;
                    data[i] = log10(-data[i]);
                } else if(data[i]>0) {
                    negative[i] = false;
                    data[i] = log10(data[i]);
                }
            }
        }
//        std::cout<<"special before global range: "<<data[13692013]<<std::endl;
//        for(int i=0;i<num;i++){
//            data[i] = fmin(fmax(low_range, data[i]), high_range);
//        }
//        std::cout<<"special after global range: "<<data[13692013]<<std::endl;
//        std::cout<<"special: "<< data[13692013]<<std::endl;
//        auto quantizer = SZ::MultipleErrorBoundsQuantizer<float>(ebs);
//        float dp= data[13692013];
//        int tmp_quant = quantizer.quantize_and_overwrite(dp, -4.921669006);
//        float dp2 = quantizer.recover(-4.921669006, tmp_quant);
//        printf("tests!!!!!!");
        startTime = std::chrono::system_clock::now();
        std::unique_ptr<unsigned char[]> compressed;
        if (fallback) {
            compressed.reset(sz_old->compress(data.get(), compressed_size));
        } else if (valueRange.isSet() && !regions.isSet() && !has_bg) {
            compressed.reset(sz->compress(data.get(), compressed_size));
        } else if (regions.isSet() && !valueRange.isSet() && !has_bg) {
            compressed.reset(sz_region->compress_region(data.get(), compressed_size,random_access));
        } else if (valueRange.isSet() && has_bg) {
            compressed.reset(
                    sz->compress_withBG(data.get(), compressed_size, bg, low_range, high_range, use_bitmap,
                                        preserve_sign,
                                        has_bg));
        } else {
            std::cerr << "The configuration must be wrong; You can use fallback, regions, multi-range, or background!"<< std::endl;
            exit(-1);
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
//        fs << outputFilePath.getValue() << " Compression Time: "
//           << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) / 1000000000
//           << "s;"
//           << "Compression size: " << compressed_size << "; Compression ratio: "
//           << num * sizeof(float) / (float) compressed_size << "; ";
        csv_result.push_back(std::to_string(num * sizeof(float) / (float) compressed_size));
        csv_result.push_back(std::to_string(double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) /1000000000));
    }
    if(mode=="test" || mode=="decompress") {
        auto compressed = SZ::readfile<unsigned char>(outputFilePath.getValue().c_str(), compressed_size);
        startTime = std::chrono::system_clock::now();
        std::unique_ptr<float[]> dec_data;
        if (fallback) {
            dec_data.reset(sz_old->decompress(compressed.get(), compressed_size));
        } else if (valueRange.isSet() && !regions.isSet() && !has_bg) {
            dec_data.reset(sz->decompress(compressed.get(), compressed_size));
        } else if(regions.isSet() && !valueRange.isSet() && !has_bg){
            dec_data.reset(sz_region->decompress_region(compressed.get(), compressed_size, random_access));
        } else if(valueRange.isSet() && has_bg){
            dec_data.reset(
                    sz->decompress_withBG(compressed.get(), compressed_size, bg, low_range, high_range, use_bitmap,
                                          preserve_sign, has_bg));
        } else {
            std::cerr << "The configuration must be wrong; You can use fallback, regions, multi-range, or background!"<< std::endl;
            exit(-1);
        }
        if(use_log) {
            for(int i=0;i<num;i++) {
                if(negative[i]) {
                    dec_data[i] = -pow(10.0, dec_data[i]);
                } else {
                    dec_data[i] = pow(10.0, dec_data[i]);
                }
            }
        }

        SZ::writefile(decFilePath.getValue().c_str(), dec_data.get(), num);
//    err = clock_gettime(CLOCK_REALTIME, &end);
        endTime = std::chrono::system_clock::now();
        std::cout << "Decompression time: "
                  //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
                  << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) /
                     1000000000 << "s"
                  << std::endl;
//        fs << "Decompression Time: "
//           << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) / 1000000000
//           << "s; ";
        csv_result.push_back(std::to_string(double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()) /1000000000));
        if(mode=="test") {
            auto dataV = SZ::readfile<float>(inputFileStr.c_str(), num);
            if (bigEndian.getValue()) {
                convert(dataV.get(), num);
            }
            if(use_log) {
                for(int i=0;i<num;i++) {
                    if(dataV[i] < 0) {
                        dataV[i] = log10(-dataV[i]);
                        dec_data[i] = log10(-dec_data[i]);
                    } else if(dataV[i]>0) {
                        dataV[i] = log10(dataV[i]);
                        dec_data[i] = log10(dec_data[i]);
                    }
                }
            }
            float max_err = 0;
            if(valueRange.isSet()) {
                std::cout << "Low: " << low_range << ", high: " << high_range << std::endl;
                float *rmse = new float[ebs.size()];
                memset(rmse, 0, sizeof(float) * ebs.size());
                float *psnr = new float[ebs.size()];
                memset(psnr, 0, sizeof(float) * ebs.size());
                float dmin = 10000, dmax = -10000;
                int *nums = new int[ebs.size()];
                memset(nums, 0, sizeof(int) * ebs.size());
                for (int i = 0; i < num; i++) {
//                if (!has_bg) {
//                    dataV[i] = fmin(fmax(low_range, dataV[i]), high_range);
//                }
                    if (dataV[i] > dmax) {
                        dmax = dataV[i];
                    }
                    if (dataV[i] < dmin) {
                        dmin = dataV[i];
                    }
                    for (int j = 0; j < ebs.size(); j++) {
                        if (dataV[i] >= ebs[j].low && dataV[i] < ebs[j].high) {
                            nums[j]++;
                            break;
                        }
                    }
                }
                int sum = 0;
                for (int j = 0; j < ebs.size(); j++) {
                    sum += nums[j];
                }
                std::cout << "total num: " << sum << std::endl;
                std::cout<<"dmax: "<< dmax << "; dmin: " << dmin << std::endl;
                float rmse_overall = 0, psnr_overall = 0;
                for (int i = 0; i < num; i++) {
                    for (int j = 0; j < ebs.size(); j++) {
                        float diff = dataV[i] - dec_data[i];
                        if (dataV[i] >= ebs[j].low && dataV[i] < ebs[j].high) {

                            float t = diff * diff / nums[j];
                            rmse[j] += t;
                        }
                        rmse_overall += diff * diff / num;
                    }
                    if (dataV[i] - dec_data[i] > max_err || dataV[i] - dec_data[i] < -max_err) {
                        max_err = (dataV[i] > dec_data[i]) ? dataV[i] - dec_data[i] : dec_data[i] - dataV[i];
                    }
                }
                for (int j = 0; j < ebs.size(); j++) {
                    psnr[j] = 20 * log10(dmax - dmin) - 10 * log10(rmse[j]);
                    rmse[j] = sqrt(rmse[j]);
                    std::cout << "[" << ebs[j].low << ", " << ebs[j].high << "] RMSE: " << rmse[j] << ", PSNR: "
                              << psnr[j] << std::endl;
                }
                psnr_overall = 20 * log10(dmax - dmin) - 10 * log10(rmse_overall);
                delete[] rmse;
                delete[] psnr;
                std::cout << "Max error = " << max_err << std::endl;
                std::cout<<" RMSE_OVERALL: " <<sqrt(rmse_overall)<< "; PSNR_OVERALL: "<< psnr_overall << std::endl;
//                fs << "Max error = " << max_err << std::endl;
            } else {
                std::cout << "Region Based - Num of elements: " << num << std::endl;
                max_err = 0;
                for (int i = 0; i < num; i++) {
                    if (dataV[i] - dec_data[i] > max_err || dataV[i] - dec_data[i] < -max_err) {
                        max_err = (dataV[i] > dec_data[i]) ? dataV[i] - dec_data[i] : dec_data[i] - dataV[i];
//                        if(max_err > 0.1) {
//                            std::cout << "dataV[" << i <<"] = "<<dataV[i] << " | dec_data[" << i <<"] = "<< dec_data[i]<<std::endl;
//                        }
                    }
                }
                std::cout << "Max error = " << max_err << std::endl;
            }
        }
    }
//    std::ofstream fs(logFilePath.getValue().c_str(), std::ios_base::trunc);
//    auto writer = csv::make_csv_writer(fs);
//    writer<< std::vector<std::string>({"CR", "CPTime", "DPTime"});
//    writer << csv_result;
//    fs.close();
    return 0;
}
