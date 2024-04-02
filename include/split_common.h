//
// Created by yuanjian on 4/1/24.
//

#ifndef SZ3_SPLIT_COMMON_H
#define SZ3_SPLIT_COMMON_H
#include "SZ3/api/sz.hpp"

namespace sz3_split {
    SZ3::Config defaultConfig() {
        SZ3::Config conf;
        conf.errorBoundMode = SZ3::EB_ABS; // refer to def.hpp for all supported error bound mode
        conf.absErrorBound = 1E-3; // absolute error bound 1e-3
        return conf;
    }

    template<typename T>
    struct DataChunk {
        size_t id;
        size_t sequenceNumber;
        std::vector<T> dataBuffer;
        std::vector<char> cpdataBuffer;
        SZ3::Config conf;
    };
}

#endif //SZ3_SPLIT_COMMON_H
