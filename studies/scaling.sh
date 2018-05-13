#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --output=out/serial.out
#SBATCH --error=out/serial.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

module load mvapich2

# weak scaling study
for i in 1,600 2,850 4,1200 8,1700 16,2400 32,3400 64,4800
do
    IFS=',' read ranks N <<< "${i}"

    mpirun -n $ranks ./cg $N parallel
    echo
done

# strong scaling study
for ranks in 1 2 4 8 16 32 64
do
    N=1600
    mpirun -n $ranks ./cg $N parallel
    echo ''
done
