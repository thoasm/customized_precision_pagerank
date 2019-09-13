#!/bin/bash

MTX_FOLDER="mm_matrices"
TEMP_FOLDER="temp"

if [[ $0 == "."* ]] || [[ $0 == "/"* ]]; then


cd $(dirname $0)

mkdir -p ${MTX_FOLDER}
cd ${MTX_FOLDER}
mkdir -p ${TEMP_FOLDER}
cd ${TEMP_FOLDER}

wget https://sparse.tamu.edu/MM/DIMACS10/adaptive.tar.gz
wget https://sparse.tamu.edu/MM/DIMACS10/delaunay_n22.tar.gz
wget https://sparse.tamu.edu/MM/DIMACS10/europe_osm.tar.gz
wget https://sparse.tamu.edu/MM/DIMACS10/hugebubbles-00020.tar.gz
wget https://sparse.tamu.edu/MM/DIMACS10/rgg_n_2_24_s0.tar.gz
wget https://sparse.tamu.edu/MM/DIMACS10/road_usa.tar.gz
wget https://sparse.tamu.edu/MM/Kamvar/Stanford.tar.gz
wget https://sparse.tamu.edu/MM/Gleich/wb-edu.tar.gz
wget https://sparse.tamu.edu/MM/SNAP/web-BerkStan.tar.gz
wget https://sparse.tamu.edu/MM/SNAP/web-Google.tar.gz

for f in $(ls *.tar.gz); do
    tar -xzf ${f}
done

for f in */; do
    if [ ${f} != ${TEMP_FOLDER}/ ]
    then
        cd ${f}
        mv *.mtx ../../
        cd ..
    fi
done

cd ../

#rm -r ${TEMP_FOLDER}


else
echo "Do not source this file. It needs to be executed!"
fi
