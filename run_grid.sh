#!/bin/bash

export MESASDK_ROOT="/mnt/home/ajermyn/Software/MESA/mesasdk"
source $MESASDK_ROOT/bin/mesasdk_init.sh
export MESA_DIR="/mnt/home/ajermyn/Software/MESA/mesa-15140"
export OMP_NUM_THREADS=1


echo "Compiling..."
cd template
./mk > ../mk.out
cd ..
echo "Done compiling."

echo "Starting runs."
mkdir runs

for j in 1.1 3.0 5.0 9.0 12 20 30 40 50 60 1.2 2.0 4.0 7.0 11 25 35 45 55 1.3 10 14 18 22 27 32 37 42 47 52 57 1.4 2.2 3.2 4.2 1.5 1.6 1.7 1.8 1.9 21 23 24 26 28 29 16 8.5 9.5 10.5 11.5 2.4 2.6 2.8 3.4 3.6 3.8 4.4 4.6 4.8 5.2 5.4 5.6 5.8 6.0 6.5 7.5 8.0

do
    rm -r runs/$j
    cp -R template runs/$j    
    cd runs/$j
    sed -i "s/MMM/$j/" inlist_project
    echo "Running Model..." $j
    ./rn > rn.out &
    cd ../../
done

wait
