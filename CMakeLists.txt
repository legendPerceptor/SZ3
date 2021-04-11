cmake_minimum_required(VERSION 3.0)
project(SZ)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")
add_definitions(-DPOLY_REGRESSION_COEF_FOLDER="${CMAKE_CURRENT_SOURCE_DIR}/data/")

message("The Project Source DIR:" ${PROJECT_SOURCE_DIR})



find_library(ZSTD_LIBS zstd HINTS "${PROJECT_SOURCE_DIR}/tools/zstd-1.4.5/lib")
set(ZSTD_INCLUDES "${PROJECT_SOURCE_DIR}/tools/zstd-1.4.5/lib/install/usr/local/include")
message("ZSTD_LIBS:" ${ZSTD_LIBS})
message("ZSTD_INCLUDES:" ${ZSTD_INCLUDES})

#find_library(HDF5_LIBS hdf5 HINTS "${PROJECT_SOURCE_DIR}/vendors/hdf5-1.10.1/install/lib")

include_directories(${ZSTD_INCLUDES}  ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/vendors/tclap-1.4.0-rc1/include)
link_directories(${PROJECT_SOURCE_DIR}/tools/zstd-1.4.5/lib)

set(HEADER_FILES
        include/compressor/Compressor.hpp
        include/compressor/SZ_General_Compressor.hpp
        include/encoder/Encoder.hpp
        include/encoder/HuffmanEncoder.hpp
        include/lossless/Lossless_zstd.hpp
        include/lossless/Lossless_bypass.hpp
        include/lossless/Lossless.hpp
        include/predictor/ComposedPredictor.hpp
        include/predictor/LorenzoPredictor.hpp
        include/predictor/PolyRegressionPredictor.hpp
        include/predictor/Predictor.hpp
        include/predictor/RegressionPredictor.hpp
        include/quantizer/IntegerQuantizer.hpp
        include/quantizer/Quantizer.hpp
        include/utils/Concepts.hpp
        include/utils/Config.hpp
        include/utils/FileUtil.h
        include/utils/Iterator.hpp
        include/utils/MemoryOps.hpp
        include/utils/Verification.hpp
        include/def.hpp
        include/SZ.hpp
        )

add_library(SZ3 SHARED
        src/dummy.cpp ${HEADER_FILES}
        )
target_link_libraries(SZ3 ${ZSTD_LIBS})

if (lzbench_LIBS)
    target_link_libraries(SZ3 ${lzbench_LIBS})
    target_include_directories(SZ3 PUBLIC ${lzbench_INCLUDES})
endif ()

add_subdirectory(src)
add_subdirectory(test)

add_subdirectory(tools/hdf5-sz3)

