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

BINARY="${BINARY_LOCATION}/${DEVICE}_create_heatmap"

DAMPING_FACTOR="0.85"
MAX_ITERS=1000


mkdir -p "${PROJECT_LOCATION}/benchmarks"
mkdir -p ${OUTPUT_FOLDER}


for mtxFile in $(ls ${MATRIX_FOLDER}*.mtx); do
        echo "Running matrix ${mtxFile}"
        ${BINARY} ${mtxFile} ${DAMPING_FACTOR} ${MAX_ITERS} >> ${OUTPUT_FOLDER}${LOCAL_TIME}_${DEVICE}_heatmap_result.m 2>> ${OUTPUT_FOLDER}${LOCAL_TIME}_${DEVICE}_heatmap_error.txt
done

else
echo "Do not source this file. It needs to be executed!"
fi
