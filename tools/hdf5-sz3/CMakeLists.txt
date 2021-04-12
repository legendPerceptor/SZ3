#find_package(HDF5 REQUIRED HINTS )
add_library(
        hdf5sz SHARED
        H5Z_SZ3.cpp
        H5Z_SZ3.h
        sz_functions.cpp)

#set(HDF5_DIR ${PROJECT_SOURCE_DIR}/vendors/hdf5-1.10.1/install)
#find_package(HDF5)
set(HDF5_LIBRARIES ${PROJECT_SOURCE_DIR}/vendors/hdf5-1.10.1/install/lib/libhdf5.dylib)
set(HDF5_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/vendors/hdf5-1.10.1/install/include)
message("HDF5LIB PATH:" ${HDF5_LIBRARIES})
message("HDF5INCLUDE PATH:" ${HDF5_INCLUDE_DIRS})
#message("CMAKE MODULEPATH:" ${CMAKE_MODULE_PATH})
#include_directories(${ZSTD_INCLUDES}  ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/vendors/tclap-1.4.0-rc1/include)
include_directories(${HDF5_INCLUDE_DIRS})

target_link_libraries(
        hdf5sz # STATIC
        SZ3
        ${HDF5_LIBRARIES}
)
#target_include_directories(
#        hdf5sz
#        PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#        $<INSTALL_INTERFACE:hdf5_sz>
#)
#
install(TARGETS hdf5sz EXPORT HDF5SZ
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz)
export(TARGETS hdf5sz SZ3 FILE HDF5SZ.cmake)


