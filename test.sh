#!/bin/sh

#$ -cwd
#$ -V
#$ -q all.q
#$ -pe openmpi1 72

mpirun -np 72  ./a.out 
