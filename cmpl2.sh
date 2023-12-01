# compile LAIK
cd laik
make
# compile HPCG
cd ..
make arch=Linux_LAIK
# clean up
# cd ./bin
rm hpcg*.txt
clear
# run HPCG
# mpirun -np 2 ./xhpcg
 ./resize_feature/launcher/tcp2run -n 3 ./bin/xhpcg 