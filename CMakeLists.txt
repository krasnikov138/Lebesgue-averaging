project(lebesgue-averaging)
cmake_minimum_required(VERSION 3.0)
add_definitions("--std=c++17 -DASIO_STANDALONE")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -pthread")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/run)

file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(main ${SOURCES})

add_library(yaml-cpp STATIC IMPORTED)
set_target_properties(yaml-cpp PROPERTIES 
	IMPORTED_LOCATION libs/yaml-cpp/lib64/libyaml-cpp.a
	INCLUDE_DIRECTORIES libs/yaml-cpp/include
)
target_link_libraries(main PRIVATE yaml-cpp)
