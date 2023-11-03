#!/bin/sh

dir="./test-cases/expansion"

# TEST CASE: 1 initial process, 1 new joining process after calling laik_allow_world_resize()
./launcher/tcp2run -n 1 -s 1 ./../bin/xhpcg

# Move result-file into .out file to be able to compare the results
mv ./hpcg*.txt ${dir}/test-hpcg-resize-1-2.out 

# Compare output with expected output
if cmp ${dir}/test-hpcg-resize-1-2.out "${dir}/test-hpcg-resize-1-2.expected"; 
    then echo "Test Case: Result is as expected";
else
    echo "Test Case: Result is not as expected";
fi