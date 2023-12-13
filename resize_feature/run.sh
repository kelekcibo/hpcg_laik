#!/bin/bash

# script to run and debug app

# compile hpcg app
cd ..
make
cd resize_feature

# delete old result files
rm -f *.txt
clear

./launcher/tcp2run -n 2 -s 2 ./../bin/xhpcg