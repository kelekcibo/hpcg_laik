# This script is to compile and run all in one, if changes in the code were made.
# compile LAIK
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
mpirun -np 2 ./xhpcg
