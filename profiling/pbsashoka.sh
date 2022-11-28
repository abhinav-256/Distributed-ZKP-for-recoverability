#! /bin/bash
#PBS -N diszkp
#PBS -o eurospreport_100_cores_10000_votes_4_mixers.out
#PBS -e eurospreport_100_cores_10000_votes_4_mixers.err
#PBS -l ncpus=100
#PBS -q cpu

cd /home/prashant.agrawal_iitd/voting/profiling
source /home/prashant.agrawal_iitd/.bashrc
mpirun -np 100 python3 dis_zkproofs.py 100 4
