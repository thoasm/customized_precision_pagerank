#!/bin/bash

LOCAL_TIME=$(date +"%y%m%d-%H%M")

if [[ $0 == "."* ]] || [[ $0 == "/"* ]]; then


cd $(dirname $0)

#DEVICE="gpu"
DEVICE="cpu"


PROJECT_LOCATION=".."
BINARY_LOCATION="${PROJECT_LOCATION}/code/${DEVICE}/build"
MATRIX_FOLDER="${PROJECT_LOCATION}/scripts/mm_matrices/"
OUTPUT_FOLDER="${PROJECT_LOCATION}/benchmarks/${DEVICE}/"

BINARY_CSR="${BINARY_LOCATION}/${DEVICE}_2018_pagerank_csr"
BINARY_ELL="${BINARY_LOCATION}/${DEVICE}_2018_pagerank_ell"

EPSILON="1e-10"
DAMPING_FACTOR="0.85"
MAX_ITERS=1000


mkdir -p "${PROJECT_LOCATION}/benchmarks"
mkdir -p ${OUTPUT_FOLDER}

# build on target system, so the "-march=native" argument for the compiler will use the proper system!


for mtxFile in $(ls ${MATRIX_FOLDER}*.mtx); do
        echo "Running matrix ${mtxFile}"
        ${BINARY_CSR} ${mtxFile} ${EPSILON} ${DAMPING_FACTOR} ${MAX_ITERS} >> ${OUTPUT_FOLDER}${LOCAL_TIME}_${DEVICE}_result_csr.m 2>> ${OUTPUT_FOLDER}${LOCAL_TIME}_${DEVICE}_error_csr.txt
        ${BINARY_ELL} ${mtxFile} ${EPSILON} ${DAMPING_FACTOR} ${MAX_ITERS} >> ${OUTPUT_FOLDER}${LOCAL_TIME}_${DEVICE}_result_ell.m 2>> ${OUTPUT_FOLDER}${LOCAL_TIME}_${DEVICE}_error_ell.txt
done


else
echo "Do not source this file. It needs to be executed!"
fi
