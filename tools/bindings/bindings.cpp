#include "SZ3/api/sz.hpp"
#include "sz3_collect.h"
#include "sz3_split.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(sz3py, m) {
    m.doc() = "Python bindings for the SZ3 compression library";

    m.def("safeMPIFinalize", []() { sz3_split::safe_call_MPI_finalize(); }, "end the MPI program");

    m.def(
        "compress",
        [](std::vector<std::string> args) {
            std::vector<char*> argv(args.size());
            for (size_t i = 0; i < args.size(); i++) {
                argv[i] = const_cast<char*>(args[i].c_str());
            }
            return sz3_split::compress(args.size(), argv.data());
        },
        "compress a file");

    m.def(
        "decompress",
        [](std::vector<std::string> args) {
            std::vector<char*> argv(args.size());
            for (size_t i = 0; i < args.size(); i++) {
                argv[i] = const_cast<char*>(args[i].c_str());
            }
            return sz3_split::decompress(args.size(), argv.data());
        },
        "decompress a file");

    m.def(
        "test",
        [](std::vector<std::string> args) {
            std::vector<char*> argv(args.size());
            for (size_t i = 0; i < args.size(); i++) {
                argv[i] = const_cast<char*>(args[i].c_str());
            }
            return sz3_split::test(args.size(), argv.data());
        },
        "Test function for SZ3");

    m.def(
        "verify_decompressed_file_in_memory",
        [](std::string original_file, std::string decompressed_file) {
            return sz3_collect::verify_decompressed_file_in_memory(original_file,
                                                                   decompressed_file);
        },
        py::arg("original_file"), py::arg("decompressed_file"),
        "Compare the original data and decompressed file in memory.");

    m.def(
        "collect_prediction_data",
        [](std::string inputFilePath, std::vector<size_t> dims, float error_bound,
           std::string config_path = "", bool use_relative_error_bound = false,
           bool use_log = false) {
            return sz3_collect::collect_prediction_data(
                inputFilePath, dims, error_bound, config_path, use_relative_error_bound, use_log);
        },
        py::arg("inputFilePath"), py::arg("dims"), py::arg("error_bound"),
        py::arg("config_path") = "", py::arg("use_relative_error_bound") = false,
        py::arg("use_log") = false, "Collect data for predicting compression performance.");
}