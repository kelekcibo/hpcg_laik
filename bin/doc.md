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
  * Handling the case for new joining procs within setup functions and init functions
  * Need to call every laik_switchto_ in the correct order as old procs do
* Discrepancies in the result
  * Issue was the variable rtz, which also needs to be sent from old proc to new pro

## 16 Test Cases and README

* Guide for how to use compile and run HPCG with LAIK
  * Changed README and created new file in setup to enable LAIK
* Scripts to reproduce test cases with shrinking/expanding

## 17 Implementing the sparse vector layout

* Overview of the interface
  * Analyzed lex_layout and generic layout interface
  * Asking advisor to understand it better
* Need to refactor code
  * current code in data.c is adjusted only to use of lex_layout
    * need a generic way
    * using a flag to know which layout is being used
* Our custom layout needs extra data from the application
  * Quick solution
    * Add pointer in Laik_Data to layout_data
    * add function to set layout_Data in Laik_Data
* Implement function to get required range for the vector layout
  * instead of expanding in both directions, we only add the size of ranges to the upper bound of the required range
* Approach for calculating the mapping from owned global to local indices by LAIK
  * Implementation of a struct storing all intervals a process contains
  * with that, we can calculate in which interval the global index is
    * during calculation, we add the interval sizes, to get the offset into allocation buffer
  * little optimization
    * if we know that an index is between two intervals, than we can immediatley break, because idx is external
* Approach for calculating the mapping from external global to local indices by LAIK
  * Analyzing in which order HPCG copies the external values to the vector after the local values
    * starting in ascending order of rank id's and starting in ascending order with global indexes to be sent
  * Just using a member in the layout object, which is increased by 1 for every calculated offset
    * if it reaches the numberOfexternalvalues stored in the layout as well, then communication should be finished and we reset the counter

## 18 Discrepancies within sparse vector layout
  
* Comparing again every function call with original application to find the discrepancies
  * First function ComputeSPMV_ref.cpp is correct
  * But I realized, that I needed to delete the old mappings due to the lex layout at some points
    * deleting all mappings gave us identical result files
  
## 19 Repartitioning with sparse vector layout

* Now, need to test it with resizing enabled
* segfaults
  * GOing iteratively through code with "exit" function to detect the error
  * did not delete mappings here as well
* Testing shrinking
  * deadlock when redistributing data according to new world size
    * Reason was when redistributing data of matrix members
    * a little if (n>0) was omitted and thus it resulted that strange things happend
    * this was the case for procs which were removed, and thus have n = 0 value, so no layout
  * FIxing this it worked, but issue with redistributing data of LAIK vectors
    * Again, removed proc does strange things as mentioned
    * But reason was, that not removed procs switched to new external partitoning. deleted that,
    * then both procs switched from old layout again. and everything worked
  * Running it further, segfault came.
    * Due to not setting layout data,little details
  * Implementing reuse function was important, because errors due to correct reuse function
  * Another error of free() and malloc()  unaligned tcache chunk detected
    * incremental search for the place where error occurs.
    * in reuse_function
      * if vector size is the same, then this means we switched from a local to external partitioning.
        // this means, that we need the Map from the local partitioning as it is not calculated for external partitioning
        // when we repartition, we are not able to make the optimisation of switching to the external partitioning first
    * After fixing this everything worked
* Testing expanding
  * Expanding from 2 to 4 procs worked immediately.
  * expanding from 1 to 2 gave us an assertion error, when we calculate the Mapping for the sparse vector layout.
    * Need incrementally detect bug.
    * bug was that we needed to initialie last chunk as well, if the last two ranges are not neighbours.
  * Compared results with orig app
    * correct results

## 20 Refactoring

* Some refactoring in the code
* Documentation of code
* Code for measuring performance of repartitioning
