# run cloud parcel model with OpenMP

export OMP_NUM_THREADS=5

nohup mpirun -host n0002 ./exemain > main.log 2>&1 &

#print the time and date again
date

