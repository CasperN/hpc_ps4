#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --output=out/scaling1.out
#SBATCH --error=out/scaling1.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

# module load mvapich2

# weak scaling study
#ws=1,600 2,850 4,1200 8,1700 16,2400 32,3400 64,4800

ws='32,3400 64,4800'

for i in echo $(echo $ws)
do
    IFS=',' read ranks N <<< "${i}"

    mpirun -n $ranks ./cg $N parallel
    echo
done

# strong scaling study
ss='16 32 64'
for ranks in $(echo $ss)
do
    N=1600
    mpirun -n $ranks ./cg $N parallel
    echo ''
done
