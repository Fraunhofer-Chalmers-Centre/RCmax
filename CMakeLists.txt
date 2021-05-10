cmake_minimum_required (VERSION 2.8.11)
project(RCMAX)
set(CMAKE_CXX_STANDARD 17)

set(CXX_COMPILE_FLAGS "-fp:precise -GL -O2 -Oi -Ot")

include_directories(include)
file(GLOB SOURCES "src/*.cpp")
add_executable(RCmax ${SOURCES})
