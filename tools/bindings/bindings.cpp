#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SZ3/api/sz.hpp"
#include "sz3_split.h"

namespace py = pybind11;

PYBIND11_MODULE(sz3py, m) {
    m.doc() = "Python bindings for the SZ3 compression library";

    m.def("compress", [](std::vector<std::string> args) {
        std::vector<char*> argv(args.size());
        for (size_t i = 0; i < args.size(); i++) {
            argv[i] = const_cast<char*>(args[i].c_str());
        }
        return sz3_split::compress(args.size(), argv.data());
    }, "compress a file");

    m.def("decompress", [](std::vector<std::string> args) {
        std::vector<char*> argv(args.size());
        for (size_t i = 0; i < args.size(); i++) {
            argv[i] = const_cast<char*>(args[i].c_str());
        }
        return sz3_split::decompress(args.size(), argv.data());
    }, "decompress a file");

    m.def("test", [](std::vector<std::string> args) {
        std::vector<char*> argv(args.size());
        for (size_t i = 0; i < args.size(); i++) {
            argv[i] = const_cast<char*>(args[i].c_str());
        }
        return sz3_split::test(args.size(), argv.data());
    }, "Test function for SZ3");
}