//
// Created by apple on 2021/6/24.
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <tclap/CmdLine.h>
#include <sstream>

#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/SZ_General_Compressor.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/LorenzoPredictor.hpp>
#include <predictor/RegressionPredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <utils/Iterator.hpp>
#include <utils/FileUtil.h>
#include "mpi.h"


static void convert(float* data, int num) {
    float t;
    char *bytes;
    for(int i=0;i<num;i++) {
        bytes = (char*)&data[i];
//        t = bytes[3] | (bytes[2] << 24) | (bytes[1] << 16) | (bytes[0] << 8);
        t = bytes[3];
        bytes[3] = bytes[0];
        bytes[0] = t;
        t = bytes[2];
        bytes[2] = bytes[1];
        bytes[1] = t;
    }
}

int main(int argc, char** argv) {
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
    TCLAP::ValueArg<std::string> outputFilePath("c", "compress", "The temporary compress file folder", true, "", "string");
//    TCLAP::MultiArg<int> dimension("d", "dimension", "the dimension of data",true,"multiInt");
    TCLAP::ValueArg<std::string> dimensionArg("d", "dimension", "the dimention of data", true, "", "string");
    TCLAP::ValueArg<std::string> valueRange("r", "range", "The ranges with low, high, and error bound; '20 50 0.1; 50 80 0.01;'", true,"", "string");
    TCLAP::SwitchArg bigEndian("e", "bigEndian", "Whether it's big endian", cmd, false);
    TCLAP::SwitchArg use_bitmapArg("p","bitmap","Whether to use the bitmap", cmd, false);
    TCLAP::SwitchArg preserve_signArg("s", "preserveSign","Whether to preserve sign", cmd, false);
    TCLAP::SwitchArg hasBackgroundData("b", "backgroundData", "Whether there is background data", cmd, false);
    TCLAP::ValueArg<std::string> decFilePath("q", "inpath", "The input raw data file", true, "", "string");
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

    // MPI Initialization
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) printf ("Start parallel compressing ... \n");
    if (world_rank == 0) printf("size: %d\n", world_size);
    double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    SZ::Compressor<float> *sz, *sz_old;
    if(dims.size()==3) {
        auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 3, 1>>(eb);
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 3>>(6, 0.1*eb_min);
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
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 2>>(6, 0.1*eb_min);
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
        auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 1>>(6, 0.1*eb_min);
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
//    std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
    double startTime, endTime;
    std::ofstream fs(logFilePath.getValue().c_str(), std::ios_base::app);
    std::string inputFileStr = inputFilePath.getValue();
    float *dataIn;
    size_t num = 0;
    char zip_filename[100];
    const char *folder=decFilePath.getValue().c_str();
    const char *tmp_folder= outputFilePath.getValue().c_str();
    char filename[100];
    sprintf(zip_filename,"%s",outputFilePath.getValue().c_str());
    int folder_index = world_rank;
    std::unique_ptr<float[]> data;
    if(mode == "test" || mode == "compress") {
        sprintf(filename, "%s/%s", folder, inputFileStr.c_str());
        if (world_rank == 0) {
            startTime = MPI_Wtime();
            data = SZ::readfile<float>(filename, num);
            if (bigEndian.getValue()) { // convert big endian data
                convert(data.get(), num);
            }
            endTime = MPI_Wtime();
            std::cout << "Read " << num << " elements\n";
            std::cout << "Original Size: " << num * sizeof(float) << std::endl;
            std::cout << "RAW File read time: "<< endTime-startTime << std::endl;
            startTime = MPI_Wtime();
            dataIn = data.get();
            MPI_Bcast(&num, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(data.get(), num, MPI_FLOAT, 0, MPI_COMM_WORLD);
            endTime = MPI_Wtime();
            std::cout << "RAW file broadcast time: " << endTime - startTime << std::endl;
        } else {
            MPI_Bcast(&num, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            dataIn = (float *) malloc(num * sizeof(float));
            MPI_Bcast(dataIn, num, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank == 0){
            endTime = MPI_Wtime();
            costReadOri += endTime - startTime;
        }
        //Compress Input Data
        size_t out_size;
        if (world_rank == 0) printf ("Compressing %s\n", filename);
        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank == 0) startTime = MPI_Wtime();
        std::unique_ptr<unsigned char[]> compressed;
        if (fallback) {
            compressed.reset(sz_old->compress(dataIn, compressed_size));
        } else if (!has_bg) {
            compressed.reset(sz->compress(dataIn, compressed_size));
        } else {
            compressed.reset(
                    sz->compress_withBG(dataIn, compressed_size, bg, low_range, high_range, use_bitmap,
                                        preserve_sign,
                                        has_bg));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank == 0) {
            endTime = MPI_Wtime();
            costComp += endTime - startTime;
        }
        // No need to free dataIn because it is a unique_ptr
//        struct stat st = {0};
//        if (stat("/lcrc/globalscratch/yuanjian", &st) == -1) {
//            mkdir("/lcrc/globalscratch/yuanjian", 0777);
//        }

        int folder_index = world_rank;
        sprintf(zip_filename, "%s/yuan_%d_%d.out", tmp_folder, folder_index, rand());
        //Write compressed data
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) printf("write compressed file to disk %s \n", zip_filename);
        if(world_rank == 0) startTime = MPI_Wtime();
        std::cout<<"compressed size: " << compressed_size << std::endl;
        SZ::writefile(zip_filename, compressed.get(), compressed_size);
        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank == 0){
            endTime = MPI_Wtime();
            costWriteZip += endTime - startTime;
        }
    }
    if(mode=="test" || mode=="decompress") {
        // Read Compressed Data
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) printf("read compressed file from disk %s \n", zip_filename);
        if(world_rank == 0) startTime = MPI_Wtime();
        auto compressed = SZ::readfile<unsigned char>(zip_filename, compressed_size);
        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank == 0){
            endTime = MPI_Wtime();
            costReadZip += endTime - startTime;
        }
        // delete the compressed file
        remove(zip_filename);

        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) printf("decompress field\n");
        if(world_rank == 0) startTime = MPI_Wtime();
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
        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank == 0){
            endTime = MPI_Wtime();
            costDecomp += endTime - startTime;
        }
//         SZ::writefile(decFilePath.getValue().c_str(), dec_data.get(), num);
        if(world_rank == 0) {
            printf ("Yuan Finish parallel compressing, total compression ratio %.4g.\n", (double)(num*sizeof(float))/(double)compressed_size);
            printf("\n");
            printf ("Timecost of reading original files = %.2f seconds\n", costReadOri);
            printf ("Timecost of reading compressed files = %.2f seconds\n", costReadZip);
            printf ("Timecost of writing compressed files = %.2f seconds\n", costWriteZip);
            printf ("Timecost of writing decompressed files = %.2f seconds\n", costWriteOut);
            printf ("Timecost of compressing using %d processes = %.2f seconds\n", world_size, costComp);
            printf ("Timecost of decompressing using %d processes = %.2f seconds\n\n", world_size, costDecomp);
        }
    }
    MPI_Finalize();
    return 0;
}