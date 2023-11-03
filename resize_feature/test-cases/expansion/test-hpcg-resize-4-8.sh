#!/bin/sh

dir="./test-cases/expansion"

# TEST CASE: 4 initial processes, 4 new joining processes after calling laik_allow_world_resize()
./launcher/tcp2run -n 4 -s 4 ./../bin/xhpcg

# Move result-file into .out file to be able to compare the results
mv ./hpcg*.txt ${dir}/test-hpcg-resize-4-8.out 

# Compare output with expected output
if cmp ${dir}/test-hpcg-resize-4-8.out "${dir}/test-hpcg-resize-4-8.expected"; 
    then echo "Test Case: Result is as expected";
else
    echo "Test Case: Result is not as expected";
fi