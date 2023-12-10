# This script is to compile and run all in one, if changes in the code were made.
# compile LAIK
PROCS=2
export LAIK_BACKEND=tcp2
export LAIK_SIZE=${PROCS}
cd laik
make
# compile HPCG
cd ..
make arch=Linux_LAIK
# clean up
cd ./bin
rm -f hpcg*.txt
clear
# run HPCG
mpirun -np ${PROCS} ./xhpcg
