#use cmake -DCMAKE_CUDA_FLAGS=”-arch=sm_30” to compile for specific architectures
#for release build:
#cmake -DCMAKE_BUILD_TYPE=Release ..

#To look for AVX and SSE support, write:
# gcc -march=native -dM -E - < /dev/null | egrep "SSE|AVX" | sort
# or
# icpc -march=native -dM -E - < /dev/null | egrep "SSE|AVX" | sort

if(WIN32)
    cmake_minimum_required(VERSION 3.9)
else()
    cmake_minimum_required(VERSION 3.8)
endif()

project(2018_pagerank LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  # -mavx2 -fno-memsafety
    set(CMAKE_CXX_FLAGS "-march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  # -mavx2 -fno-stack-protector
  #-Wall -Wno-reorder 
    set(CMAKE_CXX_FLAGS "-march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -ftree-vectorize")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # -march=core-avx2 -fno-stack-protector
  # To list optimization report: -qopt-report=5
    set(CMAKE_CXX_FLAGS "-march=native")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()


set(SOURCES
    ../mmio.c
    )

set(CPU_SOURCES
    main.cpp
    )


set(CPU_EPS_SOURCES
    main_epsilon.cpp
    )

set(CPU_BENCHMARK_SOURCES
    runbenchmark.cpp
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


set (CPU_HEADERS
    avx2_helper.hpp
    cpu_kernels.hpp
    cpu_kernel_selector.hpp
    csr_page_rank.hpp
    csr_benchmark.hpp
    split_pointer_helper.hpp
    omp_helper.hpp
    ell_benchmark.hpp
    ell_page_rank.hpp
    thread_buffer.hpp
    )


include_directories( . .. )

add_executable(cpu_2018_pagerank_csr ${SOURCES} ${CPU_SOURCES} ${HEADERS} ${CPU_HEADERS})
target_link_libraries(cpu_2018_pagerank_csr PRIVATE OpenMP::OpenMP_CXX)
target_compile_definitions(cpu_2018_pagerank_csr PRIVATE USE_ELL=false)

add_executable(cpu_2018_pagerank_ell ${SOURCES} ${CPU_SOURCES} ${HEADERS} ${CPU_HEADERS})
target_link_libraries(cpu_2018_pagerank_ell PRIVATE OpenMP::OpenMP_CXX)
target_compile_definitions(cpu_2018_pagerank_ell PRIVATE USE_ELL=true)


add_executable(cpu_benchmark ${SOURCES} ${CPU_BENCHMARK_SOURCES} ${HEADERS} ${CPU_HEADERS})
target_link_libraries(cpu_benchmark  PRIVATE OpenMP::OpenMP_CXX)
add_executable(cpu_create_heatmap ${SOURCES} ${CPU_EPS_SOURCES} ${HEADERS} ${CPU_HEADERS})
target_link_libraries(cpu_create_heatmap  PRIVATE OpenMP::OpenMP_CXX)

