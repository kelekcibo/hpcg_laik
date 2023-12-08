# High Performance Conjugate Gradient Benchmark (HPCG)

This is the README is for porting HPCG to LAIK

## Build laik

Clone with "git clone --recurse-submodules  ..." to also clone LAIK as submodule in subdirectory "laik/".
By default, this version will be linked with HPCG.
First, run "configure". MPI needs to be installed and detected for HPCG to work.

```bash
    cd laik
    ./configure
    make liblaik.so
```

## Build hpcg with laik enabled

To link HPCG with an existing LAIK installation, set variables (see *./setup/Make.Linux_LAIK, line 94*) to the path you installed
LAIK into (default: base directory with the compiled LAIK sources).

To build HPCG, just run

```bash
    make arch=Linux_LAIK
```

## Run hpcg-laik

To run hpcg with LAIK enabled

```bash
    mpirun -np 2  ./bin/xhpcg
```

## View results

hpcg will generate a result-file *hpcg\*.txt* and a report of the benchmark *HPCG-\*.txt*

## En-/Disabling REPARTITIONING

Go to file ./src/laik/hpcg_laik.hpp and (un)comment the definition *#define REPARTITION*
Then run

```bash
    # -n 2: two inital processes; -s 2: add two processes
    ./resize_feature/launcher/tcp2run -n 2 -s 2 ./../bin/xhpcg
```

For more information, see `./resize_feature/README.md`
