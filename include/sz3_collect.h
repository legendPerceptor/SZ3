#ifndef SZ3_COLLECT_H
#define SZ3_COLLECT_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

#include "SZ3/api/sz.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "qcat_dataAnalysis.h"
#include <getopt.h>

namespace sz3_collect {

struct LorenzoResult {
    double avg_err;
    double predict_cr;
    double predict_bitrate;
    double quant_entropy;
    double overhead_time;
    double p0;
    double P0;
};

double calculateQuantizationEntroy(const std::vector<int>& quant_inds, int binNumber, int nble) {
    std::vector<int> bucket(binNumber, 0);
    for (auto iter = quant_inds.begin(); iter != quant_inds.end(); iter++) {
        ++bucket[*iter];
    }
    double entVal = 0;
    for (auto iter = bucket.begin(); iter != bucket.end(); iter++) {
        if (*iter != 0) {
            double prob = double(*iter) / nble;
            entVal -= prob * log(prob) / log(2);
        }
    }
    return entVal;
}

template <uint N>
LorenzoResult lorenzo_test(float* data, std::array<size_t, N> dims, float eb = 1e-6, int r = 512) {
    auto P_l = std::make_shared<SZ3::LorenzoPredictor<float, N, 1>>(eb);
    auto quantizer = std::make_shared<SZ3::LinearQuantizer<float>>(eb, r);
    int block_size = 6;
    auto element_range = std::make_shared<SZ3::multi_dimensional_range<float, N>>(
        data, std::begin(dims), std::end(dims), 1, 0);
    auto block_range = std::make_shared<SZ3::multi_dimensional_range<float, N>>(
        data, std::begin(dims), std::end(dims), block_size, 0);
    double avg_err = 0;
    int nble = 1;
    for (int i = 0; i < N; i++) {
        nble *= dims[i];
    }

    ska::unordered_map<int, size_t> pre_freq;
    int ii = 0;
    int pre_num = 0;

    std::vector<int> quant_inds(nble);
    int quant_count = 0;

    for (auto block = block_range->begin(); block != block_range->end(); ++block) {
        element_range->update_block_range(block, block_size);
        for (auto element = element_range->begin(); element != element_range->end(); ++element) {
            float org = *element;
            ii++;
            quant_inds[quant_count++] =
                quantizer->quantize_and_overwrite(*element, P_l->predict(element));

            if (ii % 100 == 0) {
                pre_num++;
                pre_freq[quant_inds[quant_count - 1]]++;
            }

            float cur_err = fabs(*element - org);
            avg_err += cur_err / (double)nble;
        }
    }

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    double prediction = 0;
    float temp_bit = 0;
    float p_0 = (float)pre_freq[r] / pre_num;
    float P_0;
    float C_1 = 1.0;
    float pre_lossless = 1.0;
    for (int i = 1; i < r * 2 - 1; i++) {
        if (pre_freq[i] != 0) {
            temp_bit = -log2((float)pre_freq[i] / pre_num);
            if (temp_bit < 32) {
                if (temp_bit < 1) {
                    if (i == r)
                        prediction += ((float)pre_freq[i] / pre_num) * 1;
                    else if (i == r - 1)
                        prediction += ((float)pre_freq[i] / pre_num) * 2.5;
                    else if (i == r + 1)
                        prediction += ((float)pre_freq[i] / pre_num) * 2.5;
                    else
                        prediction += ((float)pre_freq[i] / pre_num) * 4;
                } else
                    prediction += ((float)pre_freq[i] / pre_num) * temp_bit;
            }
        }
    }
    if (pre_freq[0] != 0)
        prediction += ((float)pre_freq[0] / pre_num) * 32;
    if (pre_freq[r] > pre_num / 2)
        P_0 = (((float)pre_freq[r] / pre_num) * 1.0) / prediction;
    else
        P_0 = -(((float)pre_freq[r] / pre_num) * log2((float)pre_freq[r] / pre_num)) / prediction;

    pre_lossless = 1 / (C_1 * (1 - p_0) * P_0 + (1 - P_0));
    if (pre_lossless < 1)
        pre_lossless = 1;
    prediction = prediction / pre_lossless;

    double quant_entropy = calculateQuantizationEntroy(quant_inds, r * 2, nble);
    clock_gettime(CLOCK_REALTIME, &end);
    double overhead_time = (double)(end.tv_sec - start.tv_sec) +
                           (double)(end.tv_nsec - start.tv_nsec) / (double)1000000000;

    LorenzoResult result = {.avg_err = avg_err,
                            .predict_cr = 32 / prediction,
                            .predict_bitrate = prediction,
                            .quant_entropy = quant_entropy,
                            .overhead_time = overhead_time,
                            .p0 = p_0,
                            .P0 = P_0};

    return result;
}

struct CompressionResult {
    double CR;
    double CPTime;
    std::unique_ptr<char> cpdata;
    size_t outsize;
};

struct DecompressionResult {
    double PSNR;
    double RMSE;
    double DPTime;
    double WriteTime;
};

template <class T> CompressionResult compress(T* data, char* cmpPath, SZ3::Config conf) {

    size_t outSize;
    SZ3::Timer timer(true);
    std::unique_ptr<char> bytes(SZ_compress<T>(conf, data, outSize));
    double compress_time = timer.stop();

    char outputFilePath[1024];
    if (cmpPath == nullptr) {
        sprintf(outputFilePath, "tmp.sz3");
    } else {
        strcpy(outputFilePath, cmpPath);
    }
    SZ3::writefile(outputFilePath, bytes.get(), outSize);

    CompressionResult result;
    result.CR = conf.num * 1.0 * sizeof(T) / outSize;
    result.CPTime = compress_time;
    result.cpdata = std::move(bytes);
    result.outsize = outSize;

    printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(T) / outSize);
    printf("compression time = %f\n", compress_time);
    printf("compressed data file = %s\n", outputFilePath);

    return result;
}

template <class T>
DecompressionResult decompress(float* ori_data, char* cmpData, size_t cmpSize, char* decPath,
                               SZ3::Config conf, int binaryOutput) {

    SZ3::Timer timer(true);
    std::unique_ptr<T> decData(SZ_decompress<T>(conf, cmpData, cmpSize));
    double compress_time = timer.stop();
    std::cout << "Decompression finished successfully" << std::endl;
    char outputFilePath[1024];
    if (decPath == nullptr) {
        sprintf(outputFilePath, "tmp.sz3.dp");
    } else {
        strcpy(outputFilePath, decPath);
    }
    SZ3::Timer timer2(true);
    if (binaryOutput == 1) {
        SZ3::writefile<T>(outputFilePath, decData.get(), conf.num);
    } else {
        SZ3::writeTextFile<T>(outputFilePath, decData.get(), conf.num);
    }
    double write_time = timer2.stop();

    DecompressionResult dpresult;
    dpresult.DPTime = compress_time;
    dpresult.WriteTime = write_time;
    SZ3::verify<T>(ori_data, decData.get(), conf.num, dpresult.PSNR, dpresult.RMSE);
    printf("compression ratio = %f\n", conf.num * sizeof(T) * 1.0 / cmpSize);
    printf("decompression time = %f seconds.\n", compress_time);
    printf("decompressed file = %s\n", outputFilePath);
    return dpresult;
}

std::pair<double, double> verify_decompressed_file_in_memory(std::string oriFilePath,
                                                             std::string decFilePath) {
    size_t num = 0, ori_num;
    auto dp_data = SZ3::readfile<float>(decFilePath.c_str(), num);
    auto ori_data = SZ3::readfile<float>(oriFilePath.c_str(), ori_num);
    double PSNR, RMSE;
    SZ3::verify<float>(ori_data.get(), dp_data.get(), num, PSNR, RMSE);
    return std::make_pair(PSNR, RMSE);
}

std::vector<std::vector<std::string>>
collect_prediction_data(std::string inputFilePath, std::vector<size_t> dims, float error_bound,
                        std::string config_path = "", bool use_relative_error_bound = false,
                        bool use_log = false) {
    std::vector<std::vector<std::string>> result;
    std::vector<std::string> header = std::vector<std::string>({"filename",
                                                                "size",
                                                                "num",
                                                                "min",
                                                                "max",
                                                                "valueRange",
                                                                "avgValue",
                                                                "entropy",
                                                                "zeromean_variance",
                                                                "total_overhead_time",
                                                                "total_overhead_percentage",
                                                                "prediction_overhead_time",
                                                                "prediction_overhead_percentage",
                                                                "ABS Error Bound",
                                                                "Set Error Bound",
                                                                "avg_lorenzo",
                                                                "quant_entropy",
                                                                "predicted CR",
                                                                "predicted bitrate",
                                                                "p0",
                                                                "P0",
                                                                "CPTime",
                                                                "CR",
                                                                "DPTime",
                                                                "WriteTime",
                                                                "PSNR",
                                                                "RMSE"});
    result.push_back(header);
    size_t num = 0;
    auto data = SZ3::readfile<float>(inputFilePath.c_str(), num);
    if (use_log) {
        // convert data to log scale
        for (int i = 0; i < num; i++) {
            if (data[i] < 0) {
                data[i] = -log10(-data[i]);
            } else if (data[i] > 0) {
                data[i] = log10(data[i]);
            }
        }
    }

    LorenzoResult lorenzoResult;

    QCAT_DataProperty* property = computeProperty(QCAT_FLOAT, data.get(), num);
    float eb = error_bound;
    if (use_relative_error_bound) {
        eb = (property->maxValue - property->minValue) * eb;
    }

    // run lorenzo test on the data
    SZ3::Timer timer(true);
    if (dims.size() == 3) {
        std::array<size_t, 3> _dims = {{dims[0], dims[1], dims[2]}};
        lorenzoResult = lorenzo_test<3>(data.get(), _dims, eb);
    } else if (dims.size() == 2) {
        std::array<size_t, 2> _dims = {{dims[0], dims[1]}};
        lorenzoResult = lorenzo_test<2>(data.get(), _dims, eb);
    } else if (dims.size() == 1) {
        std::array<size_t, 1> _dims = {{dims[0]}};
        lorenzoResult = lorenzo_test<1>(data.get(), _dims, eb);
    }
    double overhead_time = timer.stop();

    // Do the actual compression and decompression
    SZ3::Config conf;
    if (dims.size() == 1) {
        conf = SZ3::Config(dims[0]);
    } else if (dims.size() == 2) {
        conf = SZ3::Config(dims[1], dims[0]);
    } else if (dims.size() == 3) {
        conf = SZ3::Config(dims[2], dims[1], dims[0]);
    } else {
        conf = SZ3::Config(dims[3], dims[2], dims[1], dims[0]);
    }
    if (config_path != "") {
        conf.loadcfg(config_path);
    }
    conf.absErrorBound = eb;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.num = property->numOfElem;

    auto cp_result = compress<float>(data.get(), NULL, conf);
    auto ori_data = SZ3::readfile<float>(inputFilePath.c_str(), num);
    DecompressionResult dp_result;
    dp_result =
        decompress<float>(ori_data.get(), cp_result.cpdata.get(), cp_result.outsize, NULL, conf, 1);
    std::string filename = inputFilePath.substr(inputFilePath.rfind('/') + 1);
    std::string ebArg = "ABS";
    if (use_relative_error_bound) {
        ebArg = "REL";
    }
    printf("prediction data has been collcected for %s\n", filename.c_str());
    std::vector<std::string> row =
        std::vector<std::string>({filename,
                                  std::to_string(property->totalByteSize),
                                  std::to_string(property->numOfElem),
                                  std::to_string(property->minValue),
                                  std::to_string(property->maxValue),
                                  std::to_string(property->valueRange),
                                  std::to_string(property->avgValue),
                                  std::to_string(property->entropy),
                                  std::to_string(property->zeromean_variance),
                                  std::to_string(overhead_time),
                                  std::to_string(overhead_time / cp_result.CPTime),
                                  std::to_string(lorenzoResult.overhead_time),
                                  std::to_string(lorenzoResult.overhead_time / cp_result.CPTime),
                                  std::to_string(eb),
                                  ebArg,
                                  std::to_string(lorenzoResult.avg_err),
                                  std::to_string(lorenzoResult.quant_entropy),
                                  std::to_string(lorenzoResult.predict_cr),
                                  std::to_string(lorenzoResult.predict_bitrate),
                                  std::to_string(lorenzoResult.p0),
                                  std::to_string(lorenzoResult.P0),
                                  std::to_string(cp_result.CPTime),
                                  std::to_string(cp_result.CR),
                                  std::to_string(dp_result.DPTime),
                                  std::to_string(dp_result.WriteTime),
                                  std::to_string(dp_result.PSNR),
                                  std::to_string(dp_result.RMSE)});
    result.push_back(row);
    return result;
}

void parseCollectOptions(int argc, char** argv) {}

} // namespace sz3_collect
#endif