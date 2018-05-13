#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --output=out/serial.out
#SBATCH --error=out/serial.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

module load mvapich2
./cg 75 serial_dense
./cg 75 serial_sparse
timeout 5 ./cg 10000 serial_dense
timeout 5 ./cg 10000 serial_sparse
