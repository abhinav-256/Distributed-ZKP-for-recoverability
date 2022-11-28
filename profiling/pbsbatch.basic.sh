#PBS -q standard
#PBS -N zkproofshpc
### Set the project name, your department code by default
#PBS -P zkp
### Request email when job begins and ends, don't change anything on the below line
###PBS -m bea
### Specify email address to use for notification, don't change anything on the below line
#PBS -M $USER@iitd.ac.in
#### Request your resources, just change the numbers
#PBS -l select=1:ncpus=4
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=03:00:00

# After job starts, must goto working directory.
# $PBS_O_WORKDIR is the directory from where the job is fired.
echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo "==============================="
cd $HOME/voting/profiling

#export PATH=/opt/am/python/bin/:$PATH
#module load compiler/python/3.6.0/ucs4/gnu/447
which python3

#mpirun -n 36 python3 test_mpi.py
mpirun -n 4 python3 zkproofs_mpi.py 1000000 basic
