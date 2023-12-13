
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


#ifndef HPCG_NO_MPI
#ifndef HPCG_NO_LAIK
#include "laik/hpcg_laik.hpp"
#else
#include <mpi.h>
#endif
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
  #ifndef HPCG_NO_LAIK
  hpcg_instance = laik_init(&argc, &argv);
  world = laik_world(hpcg_instance);
  #else
    MPI_Init(&argc, &argv);
  #endif // HPCG_NO_LAIK
#endif // HPCG_NO_MPI

  HPCG_Init(&argc, &argv, params); 
  printf("LAIK %d\t HI\n", laik_myid(world));

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
  std::memcpy(&hpcg_params, &params, sizeof(params));
#endif
#endif

  // Check if QuickPath option is enabled.
  // If the running time is set to zero, we minimize all paths through the program
  bool quickPath = (params.runningTime == 0);

  int size = params.comm_size, rank = params.comm_rank; // Number of MPI processes, My process ID

  bool doIO = rank == 0;

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
  #ifndef HPCG_NO_LAIK
    laik_barrier();
  #else
    MPI_Barrier(MPI_COMM_WORLD);
  #endif // HPCG_NO_LAIK
#endif // HPCG_NO_MPI
#endif // HPCG_DETAILED_DEBUG

  local_int_t nx, ny, nz;
  nx = (local_int_t)params.nx; ny = (local_int_t)params.ny; nz = (local_int_t)params.nz;
  int ierr = 0; // Used to check return codes on function calls
  ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank == 0);
  if (ierr) return ierr;

  /////////////////////////
  // Problem setup Phase //
  /////////////////////////

#ifdef HPCG_DEBUG
  double t1 = mytimer();
#endif

#ifndef HPCG_NO_LAIK
  if (doIO) HPCG_fout << "######## HPCG LAIK v1.2 ########\n#\n# New Features\n#\t-Custom layout: sparse vector\n#\n\n";
#endif // HPCG_NO_LAIK

  // Use this array for collecting timing information
  std::vector<double> times(10, 0.0);

  // Construct the geometry and linear system
  Geometry *geom = new Geometry;
  GenerateGeometry(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz, geom);

  ierr = CheckAspectRatio(0.125, geom->npx, geom->npy, geom->npz, "process grid", rank == 0);
  if (ierr) return ierr;

  int iter = 0;

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
    // Need this to know, if this proc is a new joining or old initial process.
    iter = laik_phase(hpcg_instance);
#endif
#endif

  double setup_time = mytimer();

  SparseMatrix A;
  InitializeSparseMatrix(A, geom);

  Vector b, x, xexact;
  GenerateProblem(A, &b, &x, &xexact);   
  SetupHalo(A);

#ifndef HPCG_NO_LAIK
  std::string name{""};
  name = "b_l";
  Laik_Blob *b_l = init_blob(A, false, name.data());
  name = "x_l";
  Laik_Blob *x_l = init_blob(A, false, name.data());
  name = "xexact_l";
  Laik_Blob *xexact_l = init_blob(A, false, name.data());

  // Only initial processes will do this copy
  if (iter == 0)
  {
    CopyVectorToLaikVector(b, b_l);
    CopyVectorToLaikVector(x, x_l);
    CopyVectorToLaikVector(xexact, xexact_l);
  }

#ifdef REPARTITION
  A.ptr_to_xexact = xexact_l; /* See @Laik_Blob */
#endif
#endif // HPCG_NO_LAIK

 
  int numberOfMgLevels = 4; // Number of levels including first
  SparseMatrix *curLevelMatrix = &A;
  for (int level = 1; level < numberOfMgLevels; ++level)
  {
    GenerateCoarseProblem(*curLevelMatrix);
    curLevelMatrix = curLevelMatrix->Ac; // Make the just-constructed coarse grid the next level
  }

  setup_time = mytimer() - setup_time; // Capture total time of setup
  times[9] = setup_time;               // Save it for reporting

  curLevelMatrix = &A;
  Vector *curb = &b;
  Vector *curx = &x;
  Vector *curxexact = &xexact;

  // Only initial processes will call CheckProblem
  if(iter == 0)
  {
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
      CheckProblem(*curLevelMatrix, curb, curx, curxexact);
      curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
      curb = 0;                            // No vectors after the top level
      curx = 0;
      curxexact = 0;
    }
  }
 
  CGData data;
  InitializeSparseCGData(A, data);

#ifdef REPARTITION
  // Distribute LAIK vectors now, if I am a new process
  // initial processes are calling the function to exchange values within the CG iteration, after repartitioning was done 
  if (iter > 0)
  {
    std::vector<Laik_Blob *> list{};
    /* Vectors in MG_data will be recursively handled in re_switch_LaikVectors */
    list.push_back(b_l); list.push_back(x_l); list.push_back(A.ptr_to_xexact); 
    list.push_back(data.r_blob); list.push_back(data.z_blob); list.push_back(data.p_blob);
    list.push_back(data.Ap_blob); re_switch_LaikVectors(A, list);
  }
#endif

  // Measure memory consumption
  double total_doubles = 0;

  double *ptr;
  uint64_t count;
  laik_get_map_1d(b_l->values, 0, (void **)&ptr, &count);
  total_doubles += count;
  laik_get_map_1d(x_l->values, 0, (void **)&ptr, &count);
  total_doubles += count;
  laik_get_map_1d(data.r_blob->values, 0, (void **)&ptr, &count);
  total_doubles += count;
  laik_get_map_1d(data.z_blob->values, 0, (void **)&ptr, &count);
  total_doubles += count;
  laik_get_map_1d(data.p_blob->values, 0, (void **)&ptr, &count);
  total_doubles += count;
  laik_get_map_1d(data.Ap_blob->values, 0, (void **)&ptr, &count);
  total_doubles += count;

  curLevelMatrix = &A;
  MGData *curMGData;
  for (int level = 1; level < numberOfMgLevels; ++level)
  {
    curMGData = curLevelMatrix->mgData;
    laik_get_map_1d(curMGData->Axf_blob->values, 0, (void **)&ptr, &count);
    total_doubles += count;
    laik_get_map_1d(curMGData->xc_blob->values, 0, (void **)&ptr, &count);
    total_doubles += count;
    laik_get_map_1d(curMGData->rc_blob->values, 0, (void **)&ptr, &count);
    total_doubles += count;

    curLevelMatrix = curLevelMatrix->Ac; // Next level
  }

  printf("%.1f total doubles (mem %.6f MB)\n",
         (double)total_doubles, 8 * total_doubles * 0.000001); /// 1000000);

  exit(1);


  printf("\x1B[34m ROWS: %lld \x1B[0m\n", A.totalNumberOfRows);
  ////////////////////////////////////
  // Reference SpMV+MG Timing Phase //
  ////////////////////////////////////

  // Call Reference SpMV and MG. Compute Optimization time as ratio of times in these routines

#ifndef HPCG_NO_LAIK
  name = "xexact_l";
  Laik_Blob *x_overlap = init_blob(A, true, name.data());
  name = "xexact_l";
  Laik_Blob *b_computed = init_blob(A, false, name.data());
  // Only initial processes will do this copy
  if(iter == 0)
    fillRandomLaikVector(x_overlap);
#else
  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;
  Vector x_overlap, b_computed;
  InitializeVector(x_overlap, ncol);  // Overlapped copy of x vector
  InitializeVector(b_computed, nrow); // Computed RHS vector
  // Record execution time of reference SpMV and MG kernels for reporting times
  // First load vector with random values
  FillRandomVector(x_overlap);
#endif // HPCG_NO_LAIK

  int numberOfCalls = 10;

  if (quickPath) numberOfCalls = 1; // QuickPath means we do on one call of each block of repetitive code

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
  // New processes need to "jump" to the current iteration in the CG computation
  if(iter > 0)
    numberOfCalls = 0;
#endif
#endif

  double t_begin = mytimer();
  for (int i = 0; i < numberOfCalls; ++i)
  {
#ifndef HPCG_NO_LAIK
    ierr = ComputeSPMV_laik_ref(A, x_overlap, b_computed); // b_computed = A*x_overlap
#else
    ierr = ComputeSPMV_ref(A, x_overlap, b_computed); // b_computed = A*x_overlap
#endif
    if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
#ifndef HPCG_NO_LAIK
    ierr = ComputeMG_laik_ref(A, b_computed, x_overlap); // b_computed = Minv*y_overlap
#else
    ierr = ComputeMG_ref(A, b_computed, x_overlap); // b_computed = Minv*y_overlap
#endif
    if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
  }

    if (iter > 0) times[8] = 0; // new processes skipped this part, so store 0
    else times[8] = (mytimer() - t_begin) / ((double)numberOfCalls); // Total time divided by number of calls.

#ifdef HPCG_DEBUG
  if (rank == 0)
    HPCG_fout << "Total SpMV+MG timing phase execution time in main (sec) = " << mytimer() - t1 << endl;
#endif

  printf("\x1B[31m LAIK %d \t Checkpoint 0 \x1B[0m\n", laik_myid(world));

  ///////////////////////////////
  // Reference CG Timing Phase //
  ///////////////////////////////

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
#ifndef HPCG_NO_LAIK
    ZeroLaikVector(x_l);
#ifdef REPARTITION
    if (i == 0) A.repartition_me = true; /* Repartitioning is only done in Reference CG Timing Phase */
#endif // REPARTITION
    ierr = CG_laik_ref(A, data, b_l, x_l, refMaxIters, tolerance, niters, normr, normr0, &ref_times[0], true);
#ifdef REPARTITION
    if (rank == 0) HPCG_fout << "REPARTITIONG: Call [" << i << "] Scaled Residual [" << normr / normr0 << "]" << endl;
    if (i == 0) A.repartition_me = false; /* Repartitioning is only done in Reference CG Timing Phase */
#endif // REPARTITION
#else
    ZeroVector(x);
    ierr = CG_ref(A, data, b, x, refMaxIters, tolerance, niters, normr, normr0, &ref_times[0], true);
#endif // USE_LAIK

    if (ierr) ++err_count; // count the number of errors in CG
    totalNiters_ref += niters;
  }

  if (rank == 0 && err_count) HPCG_fout << err_count << " error(s) in call(s) to reference CG." << endl;
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

  printf("\x1B[34m ROWS: %lld; n \x1B[0m\n", A.totalNumberOfRows);

  printf("\x1B[34m LAIK %d \t  Checkpoint 1 \x1B[0m\n", laik_myid(world));

  //////////////////////////////
  // Validation Testing Phase //
  //////////////////////////////

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif

  TestCGData testcg_data;
  testcg_data.count_pass = testcg_data.count_fail = 0;
#ifndef HPCG_NO_LAIK
    TestCG_laik(A, data, b_l, x_l, testcg_data);
#else
  TestCG(A, data, b, x, testcg_data);
#endif

  TestSymmetryData testsymmetry_data;
#ifndef HPCG_NO_LAIK
  TestSymmetry_laik(A, b_l, xexact_l, testsymmetry_data);
#else
  TestSymmetry(A, b, xexact, testsymmetry_data);
#endif

#ifdef HPCG_DEBUG
  if (rank == 0)
    HPCG_fout << "Total validation (TestCG and TestSymmetry) execution time in main (sec) = " << mytimer() - t1 << endl;
#endif

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif

  printf("\x1B[33m LAIK %d \t Checkpoint 2 \x1B[0m\n", laik_myid(world));

  //////////////////////////////
  // Optimized CG Setup Phase //
  //////////////////////////////

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

#ifndef HPCG_NO_LAIK
    ZeroLaikVector(x_l); // start x at all zeros
    ierr = CG_laik(A, data, b_l, x_l, optMaxIters, refTolerance, niters, normr, normr0, &opt_times[0], true);
#else
    ZeroVector(x); // start x at all zeros
    ierr = CG(A, data, b, x, optMaxIters, refTolerance, niters, normr, normr0, &opt_times[0], true);
#endif
    
    if (ierr) ++err_count; // count the number of errors in CG
    // Convergence check accepts an error of no more than 6 significant digits of relTolerance
    if (normr / normr0 > refTolerance * (1.0 + 1.0e-6)) ++tolerance_failures; // the number of failures to reduce residual

    // pick the largest number of iterations to guarantee convergence
    if (niters > optNiters) optNiters = niters;

    double current_time = opt_times[0] - last_cummulative_time;
    if (current_time > opt_worst_time) opt_worst_time = current_time;
  }

#ifndef HPCG_NO_MPI
  // Get the absolute worst time across all MPI ranks (time in CG can be different)
  double local_opt_worst_time = opt_worst_time;
#ifndef HPCG_NO_LAIK
  laik_allreduce(&local_opt_worst_time, &opt_worst_time, 1, laik_Double, LAIK_RO_Max);
#else
  MPI_Allreduce(&local_opt_worst_time, &opt_worst_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif // HPCG_NO_LAIK
#endif // HPCG_NO_MPI


  if (rank == 0 && err_count) HPCG_fout << err_count << " error(s) in call(s) to optimized CG." << endl;
  if (tolerance_failures)
  {
    global_failure = 1;
    if (rank == 0) HPCG_fout << "Failed to reduce the residual " << tolerance_failures << " times." << endl;
  }

  printf("\x1B[36m LAIK %d \t Checkpoint 3 \x1B[0m\n", laik_myid(world));

  ///////////////////////////////
  // Optimized CG Timing Phase //
  ///////////////////////////////

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
  HPCG_fout << "Number of CG sets: " << numberOfCgSets << "\n";
  for (int i = 0; i < numberOfCgSets; ++i)
  {
#ifndef HPCG_NO_LAIK
    ZeroLaikVector(x_l); // Zero out x
    ierr = CG_laik(A, data, b_l, x_l, optMaxIters, optTolerance, niters, normr, normr0, &times[0], true);
#else
    ZeroVector(x); // Zero out x
    ierr = CG(A, data, b, x, optMaxIters, optTolerance, niters, normr, normr0, &times[0], true);
#endif

    if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
    if (rank == 0) HPCG_fout << "Call [" << i << "] Scaled Residual [" << normr / normr0 << "]" << endl;
    testnorms_data.values[i] = normr / normr0; // Record scaled residual from this run
  }
  // Compute difference between known exact solution and computed solution
  // All processors are needed here.
#ifdef HPCG_DEBUG
  double residual = 0;
  ierr = ComputeResidual(A.localNumberOfRows, x, xexact, residual);
  if (ierr) HPCG_fout << "Error in call to compute_residual: " << ierr << ".\n" << endl;
  if (rank == 0) HPCG_fout << "Difference between computed and exact  = " << residual << ".\n" << endl;
#endif

  // Test Norm Results
  ierr = TestNorms(testnorms_data);
  if (ierr) HPCG_fout << "Error in call to TestNorms: " << ierr << ".\n" << endl;

  printf("\x1B[32m LAIK %d \t Checkpoint 4 \x1B[0m\n", laik_myid(world));

  ////////////////////
  // Report Results //
  ////////////////////

  // Report results to YAML file
  ReportResults(A, numberOfMgLevels, numberOfCgSets, refMaxIters, optMaxIters, &times[0], testcg_data, testsymmetry_data, testnorms_data, global_failure, quickPath);

  // Clean up
  // DeleteMatrix(A); // This delete will recursively delete all coarse grid data
  DeleteCGData(data);
  // DeleteVector(x); // TODO NEW joining procs do not init them. Handle that
  // DeleteVector(b);
  // DeleteVector(xexact);

#ifndef HPCG_NO_LAIK
  DeleteLaikVector(x_overlap);
  DeleteLaikVector(b_computed);
  DeleteLaikVector(x_l);
  DeleteLaikVector(b_l);
  DeleteLaikVector(xexact_l);

#else
  DeleteVector(x_overlap);
  DeleteVector(b_computed);
#endif

  delete[] testnorms_data.values;

  HPCG_Finalize();

  // Finish up
#ifndef HPCG_NO_MPI
  laik_finalize(hpcg_instance); // Gives double free error, if new joining procs call this functions
#else
  MPI_Finalize();
#endif
  printf("LAIK %d\tEnding program\n", laik_myid(world));
  exit_hpcg_run("Ending program", false);
  return 0;
}
