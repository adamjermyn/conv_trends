#!/bin/bash

export MESASDK_ROOT="/mnt/home/ajermyn/Software/MESA/mesasdk"
source $MESASDK_ROOT/bin/mesasdk_init.sh
export MESA_DIR="/mnt/home/ajermyn/Software/MESA/mesa-r21.12.1"
export OMP_NUM_THREADS=2


echo "Compiling..."
cd template
./mk &> ../mk.out
cd ..
echo "Done compiling."

echo "Starting runs."
mkdir runs

for j in 1.0 1.025 1.05 1.075 1.1 1.125 1.15 1.175 1.2 1.225 1.25 1.275 1.3 1.325 1.35 1.375 1.4 1.425 1.45 1.475 1.5 1.525 1.55 1.575 1.6 1.625 1.65 1.675 1.7 1.725 1.75 1.775 1.8 1.825 1.85 1.875 1.9 1.925 1.95 1.975 1.2
do
    rm -r runs/$j
    cp -R template runs/$j    
    cd runs/$j
    sed -i "s/MMM/$j/" inlist_project
    echo "Running Model..." $j
    ./rn &> rn.out &
    cd ../../
done

wait
