cmake_minimum_required(VERSION 3.9.1)
project(FactorizationOfLargeSemiprimes)
set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(IntegerFactorizer src/testing.cpp src/BigInt.cpp src/BasicAlgorithms.cpp src/QuadraticSieve.cpp Dependencies/blanczos/src/blanczos.cpp)
find_path(GMP_INCLUDE_DIR gmp.h)
include_directories(${GMP_INCLUDE_DIR} Dependencies/blanczos/src)
add_library(gmp STATIC IMPORTED)
add_library(gmpxx STATIC IMPORTED)
find_library(GMP_LIBRARY_DIR gmp)
set_target_properties(gmp PROPERTIES IMPORTED_LOCATION ${GMP_LIBRARY_DIR}) 
find_library(GMPXX_LIBRARY_DIR gmpxx)
set_target_properties(gmpxx PROPERTIES IMPORTED_LOCATION ${GMPXX_LIBRARY_DIR})
target_link_libraries(IntegerFactorizer gmp gmpxx)
