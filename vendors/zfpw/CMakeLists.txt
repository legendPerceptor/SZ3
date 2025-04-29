
add_subdirectory(zfp)

add_library(zfpw
    src/zfp_rw.c src/zfp_toolkit.c src/zfpw.c
        include/zfp_rw.h include/zfp_toolkit.h include/zfpw.h
)
target_link_libraries(zfpw zfp)

target_include_directories(zfpw PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/zfp/include
        ${CMAKE_CURRENT_SOURCE_DIR}/zfp/include/zfp
)

add_executable(test_zfp_fixtol examples/test_zfp_fixtol.c)
target_link_libraries(test_zfp_fixtol zfpw)

add_executable(test_zfp_fixratio examples/test_zfp_fixratio.c)
target_link_libraries(test_zfp_fixratio zfpw)
