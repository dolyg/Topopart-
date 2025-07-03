cmake_minimum_required(VERSION 3.5.0)
project(partitioner VERSION 0.1.0 LANGUAGES C CXX)
set(CMAKE_BUILD_TYPE Debug)
# add_compile_options(-Wall -ggdb)
add_compile_options(-O2)
add_executable(partitioner db.cpp graph.cpp log.cpp main.cpp net.cpp node.cpp parse.cpp tensor.cpp)

