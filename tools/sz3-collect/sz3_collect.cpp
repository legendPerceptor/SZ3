#include "sz3_collect.h"

// We plan to let users use the python binding to collect prediction data.
// This main function is just for testing the correctness of the collect_prediction_data function.
int main() {
    std::string inputFilePath = "/media/briannas/Research/compression-data/nyx/temperature.f32";
    std::vector<size_t> dims = {512, 512, 512};
    float error_bound = 0.01;
    std::string config_path = "";
    bool use_relative_error_bound = false;
    bool use_log = false;

    auto result = sz3_collect::collect_prediction_data(
        inputFilePath, dims, error_bound, config_path, use_relative_error_bound, use_log);
    for (auto& r : result) {
        for (auto& v : r) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}