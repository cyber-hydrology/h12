#PBS -q workq
#PBS -S /bin/bash 
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -l walltime=24:00:00 
#PBS -V



cd $PBS_O_WORKDIR

./run
