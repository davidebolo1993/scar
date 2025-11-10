cmake_minimum_required(VERSION 3.14)
project(scar)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
include(ExternalProject)

# Find required compression libraries
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

# Download htslib
FetchContent_Declare(
  htslib
  GIT_REPOSITORY https://github.com/samtools/htslib.git
  GIT_TAG 1.20
)

FetchContent_GetProperties(htslib)
if(NOT htslib_POPULATED)
  FetchContent_Populate(htslib)
  
  set(HTSLIB_INSTALL_DIR ${htslib_BINARY_DIR}/install)
  set(HTSLIB_LIB ${HTSLIB_INSTALL_DIR}/lib/libhts.a)
  set(HTSLIB_INCLUDE ${HTSLIB_INSTALL_DIR}/include)
  
  ExternalProject_Add(htslib_build
    SOURCE_DIR ${htslib_SOURCE_DIR}
    BINARY_DIR ${htslib_SOURCE_DIR}
    CONFIGURE_COMMAND bash -c "cd ${htslib_SOURCE_DIR} && ./configure --prefix=${HTSLIB_INSTALL_DIR} || (autoreconf -i && ./configure --prefix=${HTSLIB_INSTALL_DIR})"
    BUILD_COMMAND make -j 4
    INSTALL_COMMAND make install
  )
endif()

add_executable(scar
    main.cpp
    scar.cpp
)

target_include_directories(scar PRIVATE ${HTSLIB_INCLUDE})
target_link_directories(scar PRIVATE ${HTSLIB_INSTALL_DIR}/lib)

target_link_libraries(scar 
    hts
    stdc++fs
    pthread
    ZLIB::ZLIB
    BZip2::BZip2
    LibLZMA::LibLZMA
)

target_compile_options(scar PRIVATE 
    -O3 
    -Wall 
    -Wextra 
    -Wno-sign-compare
    -std=c++17
)

add_dependencies(scar htslib_build)

