#PBS -q workq
#PBS -S /bin/bash 
#PBS -l select=4:ncpus=6:mpiprocs=1:ompthreads=6
#PBS -l walltime=24:00:00 
#PBS -V

NCORE=`cat $PBS_NODEFILE | wc -l`

cd $PBS_O_WORKDIR

mpirun -machinefile $PBS_NODEFILE ./run
