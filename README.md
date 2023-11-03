########################################################
# High Performance Conjugate Gradient Benchmark (HPCG) #
########################################################

This is the README for porting LULESH to LAIK

Building
========

==

1. build laik

Clone with "git clone --recurse ..." to also clone LAIK as submodule in subdirectory "laik/".
By default, this version will be linked with HPCG.
First, run "configure". MPI needs to be installed and detected for HPCG to work.

    ```bash
    cd laik
    ./configure
    make liblaik.so
    ```
==
2. build laik-hpcg

To link HPCG with an existing LAIK installation, set variables (see *./setup/Make.Linux_LAIK, line 94*) to the path you installed
LAIK into (default: base directory with the compiled LAIK sources).

To build HPCG, just run

    ```bash
    make arch=Linux_LAIK
    ```

3. run laik-hpcg

To run hpcg with LAIK enabled

    ```bash
    mpirun -np 2  ./bin/xhpcg
    ```

4. view results

hpcg will generate a result-file *hpcg\*.txt* and a report of the benchmark *HPCG-\*.txt*
