#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
cd $PBS_O_WORKDIR
#source /data/profile.d/mpi_intelmpi-2017.4.239.sh 
mpirun -n 16 ../build_mpi/mitgcmuv
