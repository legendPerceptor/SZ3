add_executable(sz_huffman_test sz_huffman_test.cpp)
target_include_directories(sz_huffman_test PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZSTD_INCLUDES})
target_link_libraries(sz_huffman_test SZ3)

add_executable(sz_test sz_compress_test.cpp)
target_include_directories(sz_test PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZSTD_INCLUDES})
target_link_libraries(sz_test SZ3 ${ZSTD_LIBS})

add_executable(sz_demo sz_demo.cpp)
target_include_directories(sz_demo PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZSTD_INCLUDES})
target_link_libraries(sz_demo SZ3)



add_executable(sz_decompress_test decompress_test.cpp)
target_include_directories(sz_decompress_test PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZSTD_INCLUDES})
target_link_libraries(sz_decompress_test SZ3 ${ZSTD_LIBS})

add_executable(sz_cmd sz_cmd.cpp)
target_include_directories(sz_cmd PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZSTD_INCLUDES})
target_link_libraries(sz_cmd SZ3 ${ZSTD_LIBS})