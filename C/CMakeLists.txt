cmake_minimum_required(VERSION 2.6)
project(polymer)

add_executable(integrator integrator.cpp helpers.cpp)
target_link_libraries(integrator -lfftw3 -lm)

set(CMAKE_CXX_FLAGS "-g -O2 -Wall -Wextra -W -pedantic -std=c++0x")
