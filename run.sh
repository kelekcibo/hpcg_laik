# compile LAIK
cd laik
make
# compile HPCG
cd ..
make arch=Linux_LAIK
# clean up
# cd ./bin
clear
# run HPCG
# mpirun -np 2 ./xhpcg
cd ./bin
rm -f hpcg*.txt

./../resize_feature/launcher/tcp2run -n 1 -s 1 ./xhpcg 
