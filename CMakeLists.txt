project(lebesgue-averaging)
cmake_minimum_required(VERSION 3.0)
add_definitions("--std=c++17 -DASIO_STANDALONE")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -pthread")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/run)

file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(main ${SOURCES})

# add_subdirectory(libs/yaml-cpp/build)
# add_executable(yaml_config_parser "${CMAKE_SOURCE_DIR}/src/yaml_config_parser.cpp")
# target_link_libraries(yaml_config_parser PRIVATE yaml-cpp)
