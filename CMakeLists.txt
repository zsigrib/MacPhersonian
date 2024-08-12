cmake_minimum_required(VERSION 3.5.0)
project(MacPhersonian VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(MacPhersonian src/main.cpp)
add_subdirectory(src/OMtools)
#target_link_libraries(MacPhersonian PUBLIC OMtools)
target_include_directories(
    MacPhersonian PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}/src/OMtools"
)

include(CTest)
enable_testing()

