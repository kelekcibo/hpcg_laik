#!/bin/sh

dir="./test-cases/shrinking"

# TEST CASE: 2 initial process, 1 process removed after calling laik_allow_world_resize()
./launcher/tcp2run -n 2 -r 1 ./../bin/xhpcg

# Move result-file into .out file to be able to compare the results
mv ./hpcg*.txt ${dir}/test-hpcg-resize-2-1.out 

# Compare output with expected output
if cmp ${dir}/test-hpcg-resize-2-1.out "${dir}/test-hpcg-resize-2-1.expected"; 
    then echo "Test Case: Result is as expected";
else
    echo "Test Case: Result is not as expected";
fi