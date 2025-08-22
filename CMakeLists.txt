cmake_minimum_required(VERSION 3.10)
project(pbgz_main VERSION 2.0.0 LANGUAGES C CXX)

# Set C++ standard (add these lines right after project declaration)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(code)
add_subdirectory(test)
add_subdirectory(3rd_party/googletest)