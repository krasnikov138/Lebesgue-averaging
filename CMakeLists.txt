project(lebesgue-averaging)
cmake_minimum_required(VERSION 3.0)
add_definitions("--std=c++17 -DASIO_STANDALONE")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -pthread")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/run)

if (NOT PYTHON_EXECUTABLE)
    set(SOURCES 
        src/cpp/main.cpp
        src/cpp/endf.cpp
        src/cpp/measure.cpp
        src/cpp/resonances.cpp
        src/cpp/spectrum.cpp
        src/cpp/yaml_config_parser.cpp
    )
    add_executable(main ${SOURCES})

    add_library(yaml-cpp STATIC IMPORTED)
    set_target_properties(yaml-cpp PROPERTIES 
        IMPORTED_LOCATION libs/yaml-cpp/lib64/libyaml-cpp.a
        INCLUDE_DIRECTORIES libs/yaml-cpp/include
    )
    target_link_libraries(main PRIVATE yaml-cpp)
else()
    add_subdirectory(libs/pybind11)
    set(SOURCES
        src/cpp/bindings.cpp
        src/cpp/endf.cpp
        src/cpp/measure.cpp
        src/cpp/resonances.cpp
        src/cpp/spectrum.cpp
    )
    pybind11_add_module(bindings ${SOURCES})
    target_link_libraries(bindings PRIVATE pybind11::module)
endif()
