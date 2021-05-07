# cmake -DCMAKE_BUILD_TYPE=Release -Dbamboo_DIR=$(python -c 'from pkg_resources import resource_filename; print(resource_filename("bamboo", "cmake"))') .
cmake_minimum_required(VERSION 3.11)

project(BambooHMEEvaluator LANGUAGES CXX)

find_package(ROOT REQUIRED COMPONENTS Physics)
find_package(bamboo REQUIRED COMPONENTS BambooRandom)

add_library(BambooHMEEvaluator SHARED ../src/HME.cc)
target_include_directories(BambooHMEEvaluator PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BambooHMEEvaluator ROOT::Physics BambooRandom)
