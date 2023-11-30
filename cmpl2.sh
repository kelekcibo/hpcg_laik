# compile LAIK
cd laik
make
# compile HPCG
cd ..
make arch=Linux_LAIK
# clean up
# cd ./bin
# rm *.txt
clear
# run HPCG
# mpirun -np 2 ./xhpcg
 ./resize_feature/launcher/tcp2run -n 2 ./bin/xhpcg 