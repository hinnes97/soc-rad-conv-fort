#!/bin/bash

if [ -f "omp_tests/318_bands_timings" ] ; then
    rm "omp_tests/318_bands_timings"
fi
touch "omp_tests/318_bands_timings"

for i in 1 2 4 10 20 40 ; do
    echo "OMP_NUM_THREADS = ${i}"
    export OMP_NUM_THREADS=${i}
    time_average=0.0
    python run.py
    cat "omp_tests/timing.txt"
    time=$(cat omp_tests/timing.txt)
    mv "omp_tests/timing.txt" "omp_tests/timing_40bands_$i.txt"
    echo "${OMP_NUM_THREADS} threads:    ${time} s" >> 40_bands_timings
done
	 
