
add_executable(sz3_split sz3_split.cpp
        ../../include/CompressionThreadManager.h)
target_link_libraries(sz3_split ${PROJECT_NAME})