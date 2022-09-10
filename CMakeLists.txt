cmake_minimum_required(VERSION 3.22)
project(hw2)

set(CMAKE_CXX_STANDARD 14)

add_executable(hw2 src/Shading.h src/slIO.h src/slVector.cpp src/slVector.h src/trace.cpp src/trace.h src/values.h)

add_executable(hw2Prof src/Shading.h src/slIO.h src/slVector.cpp src/slVector.h src/traceProf.cpp src/trace.h src/values.h)