# RESIZE FEATURE

This feature is part of the LAIK library. It is about a resize of the world. This means that the world size (the number of processes taking part in the computation) can change dynamically with a call into following function:

```C
Laik_Group *newworld = laik_allow_world_resize(hpcg_instance, k + 1);
```

## EXAMPLE

Start with two processes "-n 2" and add two new processes (-s 2) after a call into *laik_allow_world_resize*

```bash
# ./../bin/xhpcg is the path to the executable
./launcher/tcp2run -n 2 -s 2 ./../bin/xhpcg
```

## EXPLANATION

After starting the launcher *tcp2run* will start 4 processes in total in regard to the previous example. Now, all processes jump into the function *laik_init()* and here, LAIK will allow only two processes to return and go on with computation. The important struct *Laik_Group* will return *2*, if we call *laik_size()*. This means, that we have only two processes sharing data in this group.
After a call into *laik_allow_world_resize()*, the function *laik_init()* will let the other two prcoesses "free". If we call *laik_size()* now, it will return *4*. And ressources can be shared among all four processes now.

## RUN TEST CASES

There are different test cases in this directory to illustrate a resize of the world.
Let us assume, we want to test a resize from 2 to 4 (previous example)

```bash
# Assumption: Current directory is resize_feature
./test-cases/expansion/test-hpcg-resize-2-2.sh
```

## NOTE

* Certain configurations for resizing are not supported yet. For instance, expanding from one to 3 processes.
* The scripts may echo *Result is not as expected* although the results are correct
  * *numberOfCGsets* is calculated dynamically
  