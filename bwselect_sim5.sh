#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M jlstiles@berkeley.edu 
#$ -m beas

export OMP_NUM_THREADS=1

R --vanilla < bwselect_sim5.R > bwselect_sim5.Rout

