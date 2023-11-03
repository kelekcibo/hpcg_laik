#!/bin/sh

dir="./test-cases/expansion"

# TEST CASE: 2 initial processes, 2 new joining processes after calling laik_allow_world_resize()
./launcher/tcp2run -n 2 -s 2 ./../bin/xhpcg

# Move result-file into .out file to be able to compare the results
mv ./hpcg*.txt ${dir}/test-hpcg-resize-2-2.out 

# Compare output with expected output
if cmp ${dir}/test-hpcg-resize-2-2.out "${dir}/test-hpcg-resize-2-2.expected"; 
    then echo "Test Case: Result is as expected";
else
    echo "Test Case: Result is not as expected";
fi