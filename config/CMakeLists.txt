cmake_minimum_required(VERSION 3.21)
project(KyneDraw)

set(CMAKE_CXX_STANDARD 20)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(../lib/fftw-3.3.10)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib/emsdk/upstream/emscripten/system/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib/emsdk/upstream/emscripten/cache/ports/boost_headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)