#!/bin/bash

if [[ $0 == "."* ]] || [[ $0 == "/"* ]]; then

cd $(dirname $0)


#DEVICE="gpu"
DEVICE="cpu"

PROJECT_LOCATION=".."

CODE_LOCATION="${PROJECT_LOCATION}/code/${DEVICE}"
FOLDER_BUILD="build"
BUILD_LOCATION="${CODE_LOCATION}/${FOLDER_BUILD}"

cd ${CODE_LOCATION}

mkdir -p ${FOLDER_BUILD}

cd ${FOLDER_BUILD}

if [ ${DEVICE} == "cpu" ]; then

echo "Currently at: $(pwd)"
cmake -DCMAKE_BUILD_TYPE=Release ..

make -j 10

elif [ ${DEVICE} == "gpu" ]; then

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_FLAGS="-arch=sm_30" ..

make -j 10

fi


else
echo "Do not source this file. It needs to be executed!"
fi
