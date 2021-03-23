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
namespace fs = std::filesystem;
//void myfunc(){
//    size_t my_size;
//    auto compressed_file = SZ::readfile<unsigned char>("/Users/hython/Development/globus/output/test.dat", my_size);
//    std::cout<<"read in size: "<<my_size<<std::endl;
//}
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
static void printData(float* data) {
    int z =1;
    printf("special %g\n", data[9098138]);
//    exit(1);
    for (int y=0; y< 100; y++) {
        for(int x=0;x<500;x++){
            printf("%g  ",data[x + 500*(y+500*z)]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
//    myfunc();
    size_t num = 0;
    // use Hurricane for testing
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";
    std::cout << "Original Size: " << num*sizeof(float) << std::endl;
    convert(data.get(), num);
//    data.get()[0] = 0;
    printData(data.get());

    float eb =0.1;
//    float low_range=-83.00402, high_range=31.51576;
    float low_range = -6000, high_range = 5000;
    float bg = 1.0000000e+35;
    bool preserve_sign = false;
    bool use_bitmap = true;
    int count = 0;
    for(int i =0; i<num;i++) {
        if(data[i]==0){
            count++;
        }
    }
    std::cout<<"count: "<< count<<std::endl;



    const size_t DIM = 3;

    auto P_l = std::make_shared<SZ::LorenzoPredictor<float, DIM, 1>>(eb);
    auto P_reg = std::make_shared<SZ::RegressionPredictor<float, DIM>>(6, 0.1* eb);
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, DIM>>> predictors_;
    predictors_.push_back(P_l);
    predictors_.push_back(P_reg);
//    auto cp = std::make_shared<SZ::ComposedPredictor<float, 3>>(predictors_);
    auto ebs = std::vector<std::tuple<float, float, float>>();
//    ebs.push_back(std::tuple<float,float,float>(low_range, high_range, eb));

    ebs.emplace_back(low_range, -30, eb);
    ebs.emplace_back(-30, 0, eb);
    ebs.emplace_back(0, 40, eb);
    ebs.emplace_back(40, high_range, eb);
//    ebs.emplace_back(low_range, high_range, eb);
    SZ::Config<float, DIM> conf(eb, std::array<size_t, DIM>{500, 500, 100});
    auto sz = SZ::SZ_General_Compressor<float, DIM, SZ::ComposedPredictor<float, DIM>, SZ::MultipleErrorBoundsQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf,
            SZ::ComposedPredictor<float, DIM>(predictors_),
//            SZ::LinearQuantizer<float>(eb),
            SZ::MultipleErrorBoundsQuantizer<float>(ebs),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd()
    );

//    auto quantizer = SZ::MultipleErrorBoundsQuantizer<float>(ebs);
//    float dp= data[9098138];
//    quantizer.quantize_and_overwrite(dp, 39.5);
//    printf("tests!!!!!!");
//    exit(1);
    size_t compressed_size = 0;
    // Change to use Cpp-style timer
    // struct timespec start, end;
    std::chrono::time_point<std::chrono::system_clock> startTime, endTime;

//    int err = 0;
//    err = clock_gettime(CLOCK_REALTIME, &start);
    startTime = std::chrono::system_clock::now();
    std::unique_ptr<unsigned char[]> compressed;
    compressed.reset(sz.compress_withBG(data.get(), compressed_size, bg, low_range, high_range, use_bitmap, preserve_sign,true));
//    compressed.reset(sz.compress(data.get(), compressed_size));
//    err = clock_gettime(CLOCK_REALTIME, &end);
    endTime = std::chrono::system_clock::now();
    std::cout << "Compression time: "
              //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime).count()) / 1000000000 << "s"
              << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression Ratio = " << num * sizeof (float)/ (float)compressed_size << std::endl;
    SZ::writefile(argv[2], compressed.get(), compressed_size);



//    size_t my_size;
//    auto compressed = SZ::readfile<unsigned char>("/Users/hython/Development/globus/output/test.dat", my_size);
//    std::cout<<"read in size: "<<my_size<<std::endl;



    startTime = std::chrono::system_clock::now();
//    err = clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<float[]> dec_data;
    dec_data.reset(sz.decompress_withBG(compressed.get(), compressed_size, bg, low_range, high_range, use_bitmap, preserve_sign, true));
//    dec_data.reset(sz.decompress(compressed.get(), compressed_size));
//    err = clock_gettime(CLOCK_REALTIME, &end);
    endTime = std::chrono::system_clock::now();
    std::cout << "Decompression time: "
              //              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime).count()) / 1000000000 << "s"
              << std::endl;
//    printData(dec_data.get());
    auto dataV = SZ::readfile<float>(argv[1], num);
    convert(dataV.get(), num);
    float max_err = 0;
    for (int i = 0; i < num; i++) {
//        max_err = std::max(max_err, std::abs(data[i] - dec_data[i]));
        if(dataV[i] - dec_data[i] > max_err || dataV[i] - dec_data[i] < -max_err) {
            max_err = (dataV[i] > dec_data[i]) ? dataV[i] - dec_data[i] : dec_data[i] - dataV[i];
//            printf("NOTG data: %g, dec_data: %g, i:%d\n", dataV[i], dec_data[i], i);
        }
//        if(abs(dataV[i] - dec_data[i]) > 0.1) {
//            printf("data: %g, dec_data: %g, i:%d\n", dataV[i], dec_data[i], i);
//            break;
//        }
    }
    std::cout << "Max error = " << max_err << std::endl;
    return 0;
}
