cmake_minimum_required(VERSION 3.12)
project(TwistedHessian)
 
include(FindPkgConfig)
 
pkg_check_modules(LTM REQUIRED IMPORTED_TARGET libtommath)
 
add_executable(${PROJECT_NAME} "main.c")
set(CMAKE_C_COMPILER clang)
target_link_libraries(${PROJECT_NAME} PRIVATE PkgConfig::LTM)
