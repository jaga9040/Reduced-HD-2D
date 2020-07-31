#!/bin/bash
#PBS -N RHD2D
#PBS -q serialq
#PBS -l select=1:ncpus=1:mpiprocs=1
##PBS -l walltime=00:10:00
#PBS -o ./pbsOUT/
#PBS -e ./pbsERR/
#PBS -V

cd $PBS_O_WORKDIR

./a.out #-d OUTPUT_DIR #restart append 

