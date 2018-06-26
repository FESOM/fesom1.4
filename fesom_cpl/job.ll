#!/bin/bash
#PBS -j oe
#PBS -N coupled_ucr
#PBS -l walltime=00:15:00
#PBS -l nodes=32:ppn=24
#PBS -q testq ##mpp1q
#PBS -A hbk00032
#PBS -V

echo $NCPUS
export OMP_WAIT_POLICY=PASSIVE
export CRAY_OMP_CHECK_AFFINITY=TRUE
export OMP_NUM_THREADS=1

cd /home/h/hbkdsido/fesom_echam6_oasis3-mct/fesom_cpl 
exit
date
aprun -n 768 fesom_ini.x
date

qstat -f $PBS_JOBID

#qsub job.ll
