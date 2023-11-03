clear

make arch=Linux_LAIK

cd ./bin

rm *.txt

clear

mpirun -np 2  ./xhpcg