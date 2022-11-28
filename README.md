# Recoverability-Distributed-ZKP
This is an implementation of a novel distributed ZKP that makes E2E-V voting recoverable in case of tally mismatch.

# Running the benchmarks #
==========================

0. Ensure that you have python3.7
1. Follow install instructions for installing Charm crypto library: https://jhuisi.github.io/charm/install_source.html 
   1. Install GMP 5.x
   2. Install PBC
   3. Install OpenSSL
   4. Install Charm
      - Download the source code from https://github.com/JHUISI/charm, unzip the contents and cd to that directory.
      - ./configure.sh
      - make
      - make install
2. Install MPICC and MPI4PY:
   - sudo apt install mpich
   - python3.7 -m pip install mpi4py
3. Run on c number of cores with n votes and alpha authorities (with optional argument verbosity=0/1 to print more detailed timings).:
   - mpirun -np c python3.7 dis_zkproofs.py n alpha [verbosity]
   IMPORTANT: n must be divisible by c.

By Prashant Agrawal and Abhinav Nakarmi
