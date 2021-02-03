cmake_minimum_required(VERSION 3.17)
project(hesse_curve_new)

set(CMAKE_CXX_STANDARD 20)

find_package(PkgConfig REQUIRED)
pkg_check_modules(tommath REQUIRED IMPORTED_TARGET libtommath)

add_executable(hesse_curve_new main.cpp)

target_link_libraries(hesse_curve_new PRIVATE PkgConfig::tommath)