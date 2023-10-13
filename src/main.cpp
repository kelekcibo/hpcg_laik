
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file main.cpp

 HPCG routine
 */

// Main routine of a program that calls the HPCG conjugate gradient
// solver to solve the problem, and then prints results.

// #####################

// #define HPCG_DETAILED_DEBUG
// #define HPCG_DEBUG

#ifndef HPCG_NO_MPI
#include <mpi.h>

#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik/hpcg_laik.hpp"
#endif

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#ifdef HPCG_DETAILED_DEBUG
using std::cin;
#endif
using std::endl;

#include <vector>

#include "hpcg.hpp"

#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "GenerateCoarseProblem.hpp"
#include "SetupHalo.hpp"
#include "CheckProblem.hpp"
#include "ExchangeHalo.hpp"
#include "OptimizeProblem.hpp"
#include "WriteProblem.hpp"
#include "ReportResults.hpp"
#include "mytimer.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeMG_ref.hpp"
#include "ComputeResidual.hpp"
#include "CG.hpp"
#include "CG_ref.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "CGData.hpp"
#include "TestCG.hpp"
#include "TestSymmetry.hpp"
#include "TestNorms.hpp"

/*!
  Main driver program: Construct synthetic problem, run V&V tests, compute benchmark parameters, run benchmark, report results.

  @param[in]  argc Standard argument count.  Should equal 1 (no arguments passed in) or 4 (nx, ny, nz passed in)
  @param[in]  argv Standard argument array.  If argc==1, argv is unused.  If argc==4, argv[1], argv[2], argv[3] will be interpreted as nx, ny, nz, resp.

  @return Returns zero on success and a non-zero value otherwise.

*/
int main(int argc, char *argv[])
{

  HPCG_Params params;

#ifndef HPCG_NO_MPI
  hpcg_instance = laik_init(&argc, &argv);
  world = laik_world(hpcg_instance);
#endif

  HPCG_Init(&argc, &argv, params);

  // Check if QuickPath option is enabled.
  // If the running time is set to zero, we minimize all paths through the program
  bool quickPath = (params.runningTime == 0);

  int size = params.comm_size, rank = params.comm_rank; // Number of MPI processes, My process ID
  // bool doIO = rank == 0;
  bool doIO = false;

#ifdef HPCG_DETAILED_DEBUG
  if (size < 100 && rank == 0)
    HPCG_fout << "Process " << rank << " of " << size << " is alive with " << params.numThreads << " threads." << endl;

  if (rank == 0)
  {
    printf("LAIK %d\tStart Application. %d arguments specified\n\n", rank, argc);

    char c;
    std::cout << "Press key to continue" << std::endl;
    std::cin.get(c);
  }
#ifndef HPCG_NO_MPI
  MPI_Barrier(MPI_COMM_WORLD); // Is there a way to barrier with LAIK?
#endif
#endif

  local_int_t nx, ny, nz;
  nx = (local_int_t)params.nx;
  ny = (local_int_t)params.ny;
  nz = (local_int_t)params.nz;
  int ierr = 0; // Used to check return codes on function calls

  ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank == 0);
  if (ierr)
    return ierr;

  /////////////////////////
  // Problem setup Phase //
  /////////////////////////

#ifdef HPCG_DEBUG
  double t1 = mytimer();
#endif

  if (doIO)
    printf("######## HPCG LAIK v1.1 ########\n#\n# New Features\n#\t-Expanding/Shrinking the group of processes\n#\n\n");

  if (doIO)
    printf("Start Setup Phase\n");

  // Construct the geometry and linear system
 
  Geometry *geom = new Geometry;
  GenerateGeometry(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz, geom);

  // print_HPCG_PARAMS(params, rank == 0);
  // exit_hpcg_run("PARAMS");

  ierr = CheckAspectRatio(0.125, geom->npx, geom->npy, geom->npz, "process grid", rank == 0);
  if (ierr)
    return ierr;

  // Use this array for collecting timing information
  std::vector<double> times(10, 0.0);

  double setup_time = mytimer();
  SparseMatrix A;
  InitializeSparseMatrix(A, geom);

  #ifdef REPARTITION
  std::memcpy(&hpcg_params, &params, sizeof(params));  
  // currently hpcg_params is an global variable, fix this
  // assert(hpcg_params.npx== params.npx);
  // assert(hpcg_params.npy == params.npy);
  // assert(hpcg_params.npz == params.npz);
  // assert(hpcg_params.numThreads == params.numThreads);
  // assert(hpcg_params.nx == params.nx);
  // assert(hpcg_params.ny == params.ny);
  // assert(hpcg_params.nz == params.nz);
  // assert(hpcg_params.pz == params.pz);
  // assert(hpcg_params.runningTime == params.runningTime);
  // assert(hpcg_params.zl == params.zl);
  // assert(hpcg_params.zu == params.zu);
  // exit_hpcg_run("PARAMSARETHE SAME!");
  #endif

  Vector b, x, xexact;
  GenerateProblem(A, &b, &x, &xexact);

  SetupHalo(A);

#ifdef USE_LAIK

  Laik_Blob *b_l = init_blob(A, false);
  b_l->name = "b_l";
  Laik_Blob *x_l = init_blob(A, false);
  x_l->name = "x_l";
  Laik_Blob *xexact_l = init_blob(A, false);
  xexact_l->name = "xexact_l";

  CopyVectorToLaikVector(b, b_l, A.mapping);
  CopyVectorToLaikVector(x, x_l, A.mapping);
  CopyVectorToLaikVector(xexact, xexact_l, A.mapping);

#ifdef REPARTITION
  x_l->xexact_l_ptr = xexact_l; /* See @Laik_Blob */
#endif

#endif

      int numberOfMgLevels = 4; // Number of levels including first
  SparseMatrix *curLevelMatrix = &A;
  for (int level = 1; level < numberOfMgLevels; ++level)
  {
    // HPCG_fout << "\nCoarse Problem level " << level << std::endl;
    // std::cout << "\nCoarse Problem level " << level << std::endl;
    GenerateCoarseProblem(*curLevelMatrix);
    curLevelMatrix = curLevelMatrix->Ac; // Make the just-constructed coarse grid the next level

    // #### Debug
    // printSPM(curLevelMatrix, level);
    // #### Debug
  }


  setup_time = mytimer() - setup_time; // Capture total time of setup
  times[9] = setup_time;               // Save it for reporting

  curLevelMatrix = &A;
  Vector *curb = &b;
  Vector *curx = &x;
  Vector *curxexact = &xexact;
  for (int level = 0; level < numberOfMgLevels; ++level)
  {
    CheckProblem(*curLevelMatrix, curb, curx, curxexact);
    curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    curb = 0;                            // No vectors after the top level
    curx = 0;
    curxexact = 0;
  }

  CGData data;
  InitializeSparseCGData(A, data);

  if (doIO)
    printf("End Setup Phase\n");

  ////////////////////////////////////
  // Reference SpMV+MG Timing Phase //
  ////////////////////////////////////

  if (doIO)
    printf("Start Reference SpMV+MG Timing Phase\n");

  // Call Reference SpMV and MG. Compute Optimization time as ratio of times in these routines

// #ifdef USE_LAIK
  Laik_Blob *x_overlap_l = init_blob(A, true);
  Laik_Blob *b_computed_l = init_blob(A, true);

  // fillRandomLaikVector(x_overlap_l, A.mapping);
// #else

  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;
  Vector x_overlap, b_computed;
  InitializeVector(x_overlap, ncol);  // Overlapped copy of x vector
  InitializeVector(b_computed, nrow); // Computed RHS vector

  // Record execution time of reference SpMV and MG kernels for reporting times
  // First load vector with random values
  FillRandomVector(x_overlap);
  // #endif

  CopyVectorToLaikVector(x_overlap, x_overlap_l, A.mapping);
  
  int numberOfCalls = 10;

  if (quickPath)
    numberOfCalls = 1; // QuickPath means we do on one call of each block of repetitive code
  double t_begin = mytimer();
  for (int i = 0; i < numberOfCalls; ++i)
  {
    #ifdef USE_LAIK
    ierr = ComputeSPMV_laik_ref(A, x_overlap_l, b_computed_l); // b_computed = A*x_overlap
    #else
    ierr = ComputeSPMV_ref(A, x_overlap, b_computed); // b_computed = A*x_overlap
    #endif

    if (ierr)
      HPCG_fout << "Error in call to SpMV: " << ierr << ".\n"
                << endl;

    #ifdef USE_LAIK
    ierr = ComputeMG_laik_ref(A, b_computed_l, x_overlap_l); // b_computed = Minv*y_overlap
    #else
    ierr = ComputeMG_ref(A, b_computed, x_overlap); // b_computed = Minv*y_overlap
    #endif

    if (ierr)
      HPCG_fout << "Error in call to MG: " << ierr << ".\n"
                << endl;
  }
  times[8] = (mytimer() - t_begin) / ((double)numberOfCalls); // Total time divided by number of calls.

  // printf("Reference SpMV+MG Timing Phase: %.5f seconds\n", times[8]);

#ifdef HPCG_DEBUG
  if (rank == 0)
    HPCG_fout << "Total SpMV+MG timing phase execution time in main (sec) = " << mytimer() - t1 << endl;
#endif

  if (doIO)
    printf("End Reference SpMV+MG Timing Phase\n");

  ///////////////////////////////
  // Reference CG Timing Phase //
  ///////////////////////////////

  if (doIO)
    printf("Start Reference CG Timing Phase\n");

#ifdef HPCG_DEBUG
      t1 = mytimer();
#endif
  int global_failure = 0; // assume all is well: no failures

  int niters = 0;
  int totalNiters_ref = 0;
  double normr = 0.0;
  double normr0 = 0.0;
  int refMaxIters = 50;
  numberOfCalls = 1; // Only need to run the residual reduction analysis once

  // Compute the residual reduction for the natural ordering and reference kernels
  std::vector<double> ref_times(9, 0.0);
  double tolerance = 0.0; // Set tolerance to zero to make all runs do maxIters iterations
  int err_count = 0;
  for (int i = 0; i < numberOfCalls; ++i)
  {
    #ifdef USE_LAIK
      ZeroLaikVector(x_l, A.mapping);
      
    #ifdef REPARTITION
      if (i == 0)
        A.repartition_me = true; /* Repartitioning is only done in Reference CG Timing Phase */
    #endif // REPARTITION
  
      ierr = CG_laik_ref(A, data, b_l, x_l, refMaxIters, tolerance, niters, normr, normr0, &ref_times[0], true);
    
    #ifdef REPARTITION
      if (i == 0)
        A.repartition_me = false;  /* Repartitioning is only done in Reference CG Timing Phase */
    #endif // REPARTITION

    #else
      ZeroVector(x);
      ierr = CG_ref(A, data, b, x, refMaxIters, tolerance, niters, normr, normr0, &ref_times[0], true);
    #endif // USE_LAIK

    if (ierr)
      ++err_count; // count the number of errors in CG
    totalNiters_ref += niters;
  }

  if (rank == 0 && err_count)
    HPCG_fout << err_count << " error(s) in call(s) to reference CG." << endl;
  double refTolerance = normr / normr0;

  // Call user-tunable set up function.
  double t7 = mytimer();
  OptimizeProblem(A, data, b, x, xexact);
  t7 = mytimer() - t7;
  times[7] = t7;
#ifdef HPCG_DEBUG
  if (rank == 0)
    HPCG_fout << "Total problem setup time in main (sec) = " << mytimer() - t1 << endl;
#endif

#ifdef HPCG_DETAILED_DEBUG
// TODO. Copy Laik_blobs to b, x, xexact, if LAIK with one proc is used
  if (geom->size == 1)
    WriteProblem(*geom, A, b, x, xexact);
#endif

  if (doIO)
    printf("End Reference CG Timing Phase\n");



  //////////////////////////////
  // Validation Testing Phase //
  //////////////////////////////

  if (doIO)
    printf("Start Validation Testing Phase\n");

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif
  // TestCGData testcg_data_laik;
  TestCGData testcg_data;

  testcg_data.count_pass = testcg_data.count_fail = 0;
  // testcg_data_laik.count_pass = testcg_data_laik.count_fail = 0;

#ifdef USE_LAIK
  TestCG_laik(A, data, b_l, x_l, testcg_data);
#else
  TestCG(A, data, b, x, testcg_data);
#endif

  // printf("Test TestCGData\n");

  // assert(testcg_data.count_fail == testcg_data_laik.count_fail);
  // assert(testcg_data.count_pass == testcg_data_laik.count_pass);
  // assert(testcg_data.expected_niters_no_prec == testcg_data_laik.expected_niters_no_prec);
  // assert(testcg_data.expected_niters_prec == testcg_data_laik.expected_niters_prec);
  // assert(testcg_data.niters_max_no_prec == testcg_data_laik.niters_max_no_prec);
  // assert(testcg_data.niters_max_prec == testcg_data_laik.niters_max_prec);
  // assert(testcg_data.normr == testcg_data_laik.normr);


  TestSymmetryData testsymmetry_data;
  // TestSymmetryData testsymmetry_data_laik;

#ifdef USE_LAIK
  TestSymmetry_laik(A, b_l, xexact_l, testsymmetry_data);
#else
    TestSymmetry(A, b, xexact, testsymmetry_data);
  #endif

    // printf("Test TestSymmetryData\n");

    // assert(testsymmetry_data.count_fail == testsymmetry_data_laik.count_fail);
    // // will fail since both calls TestSymmetry fill x vector with random values
    // compare2(testsymmetry_data.depsym_mg, testsymmetry_data_laik.depsym_mg, false, 0);
    // compare2(testsymmetry_data.depsym_spmv, testsymmetry_data_laik.depsym_spmv, false, 0);

    // printf("normal; %.20f\t laik: %.20f\n", testsymmetry_data.depsym_mg, testsymmetry_data_laik.depsym_mg);
    // printf("normal; %.20f\t laik: %.20f\n", testsymmetry_data.depsym_spmv, testsymmetry_data_laik.depsym_spmv);

#ifdef HPCG_DEBUG
  if (rank == 0)
    HPCG_fout << "Total validation (TestCG and TestSymmetry) execution time in main (sec) = " << mytimer() - t1 << endl;
#endif

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif

  if (doIO)
    printf("End Validation Testing Phase\n");

  //////////////////////////////
  // Optimized CG Setup Phase //
  //////////////////////////////

  if (doIO)
    printf("Start Optimized CG Setup Phase \n");

  niters = 0;
  normr = 0.0;
  normr0 = 0.0;
  err_count = 0;
  int tolerance_failures = 0;

  int optMaxIters = 10 * refMaxIters;
  int optNiters = refMaxIters;
  double opt_worst_time = 0.0;

  std::vector<double> opt_times(9, 0.0);

  // Compute the residual reduction and residual count for the user ordering and optimized kernels.
  for (int i = 0; i < numberOfCalls; ++i)
  {
    double last_cummulative_time = opt_times[0];

    #ifdef USE_LAIK
      ZeroLaikVector(x_l, A.mapping); // start x at all zeros
      ierr = CG_laik(A, data, b_l, x_l, optMaxIters, refTolerance, niters, normr, normr0, &opt_times[0], true);
    #else
      ZeroVector(x); // start x at all zeros
      ierr = CG(A, data, b, x, optMaxIters, refTolerance, niters, normr, normr0, &opt_times[0], true);
    #endif
    
    if (ierr)
      ++err_count; // count the number of errors in CG
    // Convergence check accepts an error of no more than 6 significant digits of relTolerance
    if (normr / normr0 > refTolerance * (1.0 + 1.0e-6))
      ++tolerance_failures; // the number of failures to reduce residual

    // pick the largest number of iterations to guarantee convergence
    if (niters > optNiters)
      optNiters = niters;

    double current_time = opt_times[0] - last_cummulative_time;
    if (current_time > opt_worst_time)
      opt_worst_time = current_time;
  }

#ifndef HPCG_NO_MPI
  // Get the absolute worst time across all MPI ranks (time in CG can be different)
  double local_opt_worst_time = opt_worst_time;
  laik_allreduce(&local_opt_worst_time, &opt_worst_time, 1, laik_Double, LAIK_RO_Max);
#endif

  if (rank == 0 && err_count)
    HPCG_fout << err_count << " error(s) in call(s) to optimized CG." << endl;
  if (tolerance_failures)
  {
    global_failure = 1;
    if (rank == 0)
      HPCG_fout << "Failed to reduce the residual " << tolerance_failures << " times." << endl;
  }

  if (doIO)
    printf("End Optimized CG Setup Phase \n");

  ///////////////////////////////
  // Optimized CG Timing Phase //
  ///////////////////////////////

  if (doIO)
    printf("Start Optimized CG Timing Phase \n");

  // Here we finally run the benchmark phase
  // The variable total_runtime is the target benchmark execution time in seconds

  double total_runtime = params.runningTime;
  int numberOfCgSets = int(total_runtime / opt_worst_time) + 1; // Run at least once, account for rounding

#ifdef HPCG_DEBUG
  if (rank == 0)
  {
    HPCG_fout << "Projected running time: " << total_runtime << " seconds" << endl;
    HPCG_fout << "Number of CG sets: " << numberOfCgSets << endl;
  }
#endif

  /* This is the timed run for a specified amount of time. */

  optMaxIters = optNiters;
  double optTolerance = 0.0; // Force optMaxIters iterations
  TestNormsData testnorms_data;
  testnorms_data.samples = numberOfCgSets;
  testnorms_data.values = new double[numberOfCgSets];

  for (int i = 0; i < numberOfCgSets; ++i)
  {
    #ifdef USE_LAIK
      ZeroLaikVector(x_l, A.mapping); // Zero out x
      ierr = CG_laik(A, data, b_l, x_l, optMaxIters, optTolerance, niters, normr, normr0, &times[0], true);
    #else
      ZeroVector(x); // Zero out x
      ierr = CG(A, data, b, x, optMaxIters, optTolerance, niters, normr, normr0, &times[0], true);
    #endif

    if (ierr)
      HPCG_fout << "Error in call to CG: " << ierr << ".\n"
                << endl;
    if (rank == 0)
      HPCG_fout << "Call [" << i << "] Scaled Residual [" << normr / normr0 << "]" << endl;
    testnorms_data.values[i] = normr / normr0; // Record scaled residual from this run
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.
#ifdef HPCG_DEBUG
  double residual = 0;
  ierr = ComputeResidual(A.localNumberOfRows, x, xexact, residual);
  if (ierr)
    HPCG_fout << "Error in call to compute_residual: " << ierr << ".\n"
              << endl;
  if (rank == 0)
    HPCG_fout << "Difference between computed and exact  = " << residual << ".\n"
              << endl;
#endif

  // Test Norm Results
  ierr = TestNorms(testnorms_data);
  if (ierr)
    HPCG_fout << "Error in call to TestNorms: " << ierr << ".\n"
              << endl;
  if (doIO)
    printf("End  Optimized CG Timing Phase \n");

  ////////////////////
  // Report Results //
  ////////////////////

  // Report results to YAML file
  ReportResults(A, numberOfMgLevels, numberOfCgSets, refMaxIters, optMaxIters, &times[0], testcg_data, testsymmetry_data, testnorms_data, global_failure, quickPath);

  // Clean up
  // DeleteMatrix(A); // This delete will recursively delete all coarse grid data
  // DeleteCGData(data);
  // DeleteVector(x);
  // DeleteVector(b);
  // DeleteVector(xexact);

  #ifdef USE_LAIK
    // delete laik containers

  #else
    DeleteVector(x_overlap);
    DeleteVector(b_computed);
  #endif

      delete[] testnorms_data.values;
  
  HPCG_Finalize();
  
  if (doIO)
    printf("Done\n");

  // Finish up
#ifndef HPCG_NO_MPI
  laik_finalize(hpcg_instance);
#endif

  return 0;
}
