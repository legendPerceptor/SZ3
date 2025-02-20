# cmake_policy(SET CMP0148 NEW)
find_package(pybind11 REQUIRED)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})

pybind11_add_module(sz3py bindings.cpp)
target_link_libraries(sz3py PRIVATE ${PROJECT_NAME} ${MPI_C_LIBRARIES})
message("PROJECT_ROOT: " ${PROJECT_ROOT})
install(TARGETS sz3py
    LIBRARY DESTINATION ${PROJECT_ROOT}/tools/bindings/sz3py
)
