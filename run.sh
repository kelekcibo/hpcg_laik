# This script is to compile and run all in one, if changes in the code were made.
# compile LAIK
PROCS=8
# compile HPCG
make arch=Linux_MPI
# clean up
cd ./bin
rm -f hpcg*.txt
clear
# run HPCG
mpirun -np ${PROCS} ./xhpcg
