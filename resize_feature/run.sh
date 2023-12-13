#!/bin/bash

# script to run and debug app

# compile hpcg app
cd ..
make arch=Linux_LAIK
cd resize_feature

# delete old result files
rm -f *.txt
clear

./launcher/tcp2run -n 4 -s 4 ./../bin/xhpcg
