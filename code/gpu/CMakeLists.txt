#use cmake -DCMAKE_CUDA_FLAGS=”-arch=sm_30” to compile for specific architectures
#for release build:
#cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_FLAGS="-arch=sm_30" ..


if(WIN32)
    cmake_minimum_required(VERSION 3.9)
else()
    cmake_minimum_required(VERSION 3.8)
endif()

project(2018_pagerank LANGUAGES CXX C CUDA)

set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED true)
set(CMAKE_CUDA_STANDARD 11)

set(SOURCES
    ../mmio.c
    )

set(GPU_SOURCES
    main.cu
    )


set(GPU_EPS_SOURCES
    main_epsilon.cu
    )

set(GPU_BENCHMARK_SOURCES
    runbenchmark.cu
    )

set(HEADERS
    ../constants.hpp
    ../precision.hpp
    ../macros.hpp
    ../split_pointer.hpp
    ../types.hpp
    ../ell_matrix.hpp
    ../coo_matrix.hpp
    ../csr_matrix.hpp
    ../mmio.h
    )

set (GPU_HEADERS
    split_pointer_helper.cuh
    gpu_kernels.cuh
    ell_page_rank.cuh
    csr_page_rank.cuh
    csr_benchmark.cuh
    ell_benchmark.cuh
    )

find_library(CUBLAS cublas
    HINT ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
find_library(CUSPARSE cusparse
    HINT ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})


include_directories(. ../)

add_executable(gpu_2018_pagerank_csr ${SOURCES} ${GPU_SOURCES} ${HEADERS} ${GPU_HEADERS})
target_compile_definitions(gpu_2018_pagerank_csr PRIVATE USE_ELL=false)
add_executable(gpu_2018_pagerank_ell ${SOURCES} ${GPU_SOURCES} ${HEADERS} ${GPU_HEADERS})
target_compile_definitions(gpu_2018_pagerank_ell PRIVATE USE_ELL=true)

add_executable(gpu_benchmark ${SOURCES} ${GPU_BENCHMARK_SOURCES} ${HEADERS} ${GPU_HEADERS})
add_executable(gpu_create_heatmap ${SOURCES} ${GPU_EPS_SOURCES} ${HEADERS} ${GPU_HEADERS})


target_link_libraries(gpu_benchmark PRIVATE ${CUBLAS} ${CUSPARSE})
