# Current state

## 0 Getting to know to LAIK and HPCG Code Base

* Concepts of LAIK
* Jump into the HPCG Code Base

## 1 Replacing all MPI calls with LAIK calls

* Implementation of one function  with buffers of variable size and data type, which is used for all reductions and broadcast operations
  * Preventing reduntant code

## 2 Implementation of the partitioner

* First try: Partitioned the LAIK container as it were an local partitiong. As we are partitioning a global problem (, and LAIK does the local part under the hood), the partitioner was adjusted.  
* Second try: Identification of which information calculated in SetupHalo is needed for our partitioner
  * First approach of the global partitioning resulted in a deadlock, because, processes did not calculate the same partitioning.
    * The reason was that we only specified what proc i needs
  * Second approach of global partitioning Specification of what other procs need from proc i, so the partitiong was the same for every proc
    * No deadlock anymore

## 3 Implementation of the mapping (due to the lex_layout)

* Mapping from local indices to global indices, and then to allocation indices
  * Identifying which information is needed. New Struct called L2A_map.
* Tested the partitioner by copying all values from Vector x to the container, calling a switchto, performing the SPMV operation, and copying the result back to the Vector x
  * The result file was the same => Partitioner and thus, exchange via LAIK was succesfull

## 4 Making use of LAIK Containers during the whole application

* First try: Adding parameters to the already existing functions where exchange is done.
  * It lead to confusing code as backwards compatibality should be supported.
  * Due to the lex_layout, it was not possible to reuse the access pattern of the hpcg application (=> implement hpcg_layout (?))
  * In some functions, it would lead to adding reduntant parameters as it takes vectors with exchange and no exchange
    * also would lead to confusing code within the functions and when calling the functions (passing NULL)
      * e.g.

```c
ComputeWAXPBY(nrow, 1.0, x, alpha, p, x, A.isWaxpbyOptimized, NULL, p_blob, NULL);
ComputeWAXPBY(nrow, 1.0, r, -alpha, Ap, r, A.isWaxpbyOptimized, NULL, NULL, NULL); 
ComputeDotProduct(nrow, r, r, normr, t4, A.isDotProductOptimized, NULL, NULL);
```

* Second try: Implementation of prototypes of the operations/functions with LAIK only containers.
  * Simplicity: Even vectors with no exchange are implemented with LAIK containers.
  * => cleaner and backwards compatibilty is guruanteed in a nicer way.

## 5 Partitioning the LAIK containers

* First try: In SetupHalo, all data structures needed for the partitioner algorithm were intitialized.
   Then, for every Vector, an Struct called Laik_Blob was introduced.
   In the first try, after initializing each laik_blob resulted in creating the partionings in every call to init_blob and allocating new buffer, copying the data for the partitioner algorithm and mapping, and saving the pointer in Laik_Blob.
   This resulted in very bad performance => Reduntant bytes and delay (creating same partitiong n times)
   In this try, after the first calls into ComputeSPMV, an Bus error occured and the program execution was stopped.
* Second try:
  In SetupHalo, all data structures needed for the partitioner algorithm were intitialized.
  This time, the partitiongs were initialized once and saved only in the SparseMatrix struct.
  When calling init_blob, only the new Laik_data is created and partitioned with the partitioning stored in the SparseMatrix
  With this try, the performance was optimized. After trying to compile with this new approach of using laik containers, building issues occured.

## 6 Building issues after implementation of laik functions

* Some forward decl. were needed
* Linking errors:

```c
  // Function declaration in ComputeMG_ref.hpp was
   int ComputeMG_laik_ref(const SparseMatrix &A, Laik_Blob *r, Laik_Blob *x);
  // instead of
   int ComputeMG_laik_ref(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
  /* "const" was missing. It resolved the building issues. */
```

## 7 Discrepancies in the results

* In the first try, bus errors occured. But with the new approach, there were no bus errors anymore.
* There were at first some segmentation faults:
  Did not use the mapping from local to allocation indices at the correct place. (In functions like Prolongation_ref.cpp).
* Afterwards, the program run without any errors, but result was wrong.
* Approach to find the errors was to call the laik version and original version, and then compare both vector results.
* SPMV operation gave the same results.
* The next function in the main function of the HPCG code was MG_REF. And the discrepancy came from there.
* Inside of MG_ref, the first operation was SYMGS_REF. After comparing both vectors after a call into symgs_ref, an discrepancy occured.
* After comparing the original version of the function with the laik version, there was no copy mistake in sight.
  But then, I realized that at each iteration the variable sum is initialized with an entry of the r vector.
  And the problem here was, that the local index was not mapped to the allocation index, which resulted in big discrepancies in the result file.
* After this fix, the tests of equality of the vectores passed.
* The result file contained the same results as if the original application produced the output.
* v1.0 of HCPG LAIK

## 8 Enabling expand/shrink features

* Getting to know this feature => examples in LAIK repo
* Re-run of setup functions with new world size
  * Just re-calling setup functions with updated config parameters does not work
    * Need to dig into initialization logic
  * Need to adjust nx, ny, nz as npx, npy, npz could change according to the new size and global size is calculated using npx and nx, (y,z resp.)
    * dynamically calculating new nx, creating new geom for every proc and re-generating problem
      * except for removed procs, they only calculate needed data to do the last switch
* Re-Setup and Switch worked after some assertion errors and seg faults
* Need to update local lengths for the already existing laik containers
  * Forgot Afx_blob, but now all lengths are updated and switched to
* Disceprancies in the results for the test case -np 2 and shrinking to np 1
  * non deterministic behaviour
  * Analyzing the discrepancy was done in step 13

## 9 Refactoring of code

* Maintainability was bad
  * Refactoring of code
  * Splitting the one big file into multiple smaller files

## 11 Partitioner algorithm for the SparseMatrix

* In 8, setup functions were recalled
  * Values of matrix are deleted and re-init'ed
* New Approach
  * Partitioner algorithm for matrix A, so there is no need to call setup functions again
    * Analysing which data of A needs to be distributed
  * Code for repartitioning within CG_REF.cpp adjusted
    * But code for new joining procs not done yet. See step 15

## 12 Adjusting setup and other functions

* Usage of LAIK Containers in the matrix A results in changes of code which uses data of the matrix, as we need to get base ptr
* Due to the lex_layout, we need a mapping here as well
  * Same approach as in step 3
* Discussion of how to implement the member matrixValues and mtxIndG
  * Original HPCG code uses 2D array
  * LAIK 1D container instead of 2D container (flattening the 2D arrays)

## 13 Disceprancies in the result

* After implementing step 11 and 12, there were disceprancies in the result
* In SPMV_ref with repartitioning feature
  * "Off by one" mistake during initialization of matrixValues, matrixDiagonal and mtxIndG
* After call into TestCG, again discrepancies
  * Functions replacing matrix diagonal with exxaggerated diagonal was the origin of this discrepancy
    * The data in matrixValues was different than in the original application
    * This means, we were accessing wrong values in the LAIK version
    * Pointers are stored in the original application, but in LAIK we do not store pointers in matrixDiagonal
    * Assigning new values in matrixValues as well fixed the issue

## 14 Disceprancies in the result with shrinking

* Test case: 2 procs shrinking to 1 proc
* Seg fault
  * Issue with f2c container
  * Using wrong data type (laik_UInt64) was the reason for the seg fault
    * changing to to laik_Int32 fixed the problem
* Almost the same results
  * Reason: f2c container is local to every own proc. So this is not a global variable.
  * After reinitializing f2cOperator during repartitioning fixed the discrepancy and we get the same results
  * Test Case
    * mpirun -np 2; 16 16 16
    * Shrinking from two to one process

## 15 Code for new joining processes

* Need to consider new size
  * For certain configurations, it is not possible to get the same problem size
  * Example
    * -np 2, adding one new process; with problem size: 8192 totalRows
      * not possible to get 8192 rows with 3 procs
* Need to add code for new joining processes as they need to do only certain things (not everything in setup functions)
