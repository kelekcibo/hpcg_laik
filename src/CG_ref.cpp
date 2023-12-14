
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
 @file CG_ref.cpp

 HPCG routine
 */

#include <fstream>
#include <cmath>

#include "laik/hpcg_laik.hpp"
#include "hpcg.hpp"
#include "CG_ref.hpp"
#include "mytimer.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeMG_ref.hpp"
#include "ComputeDotProduct_ref.hpp"
#include "ComputeWAXPBY_ref.hpp"

// Use TICK and TOCK to time a code section in MATLAB-like fashion
#define TICK()  t0 = mytimer() //!< record current time in 't0'
#define TOCK(t) t += mytimer() - t0 //!< store time difference in 't' using time in 't0'

#ifndef HPCG_NO_LAIK
/*!
  Reference routine to compute an approximate solution to Ax = b

  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess; on exit: the new approximate solution
  @param[in]    max_iter  The maximum number of iterations to perform, even if tolerance is not met.
  @param[in]    tolerance The stopping criterion to assert convergence: if norm of residual is <= to tolerance.
  @param[out]   niters    The number of iterations actually performed.
  @param[out]   normr     The 2-norm of the residual vector after the last iteration.
  @param[out]   normr0    The 2-norm of the residual vector before the first iteration.
  @param[out]   times     The 7-element vector of the timing information accumulated during all of the iterations.
  @param[in]    doPreconditioning The flag to indicate whether the preconditioner should be invoked at each iteration.

  @return Returns zero on success and a non-zero value otherwise.

  @see CG()
*/
int CG_laik_ref(SparseMatrix &A, CGData &data, Laik_Blob *b, Laik_Blob *x,
           const int max_iter, const double tolerance, int &niters, double &normr, double &normr0,
           double *times, bool doPreconditioning)
{

  double t_begin = mytimer(); // Start timing right away
  normr = 0.0;
  double rtz = 0.0, oldrtz = 0.0, alpha = 0.0, beta = 0.0, pAp = 0.0;

  double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;
  // #ifndef HPCG_NO_MPI
  //   double t6 = 0.0;
  // #endif
  double t_rep = 0;
  double t_it_before = 0;
  double t_it_after = 0;

  local_int_t nrow = A.localNumberOfRows;

  Laik_Blob * r = data.r_blob; // Residual vector
  Laik_Blob * z = data.z_blob;  // Preconditioned residual vector
  Laik_Blob * p = data.p_blob;   // Direction vector (in MPI mode ncol>=nrow)
  Laik_Blob * Ap = data.Ap_blob;

  if (!doPreconditioning && A.geom->rank == 0)
    HPCG_fout << "WARNING: PERFORMING UNPRECONDITIONED ITERATIONS" << std::endl;

#ifdef HPCG_DEBUG
  int print_freq = 1;
  if (print_freq > 50)
    print_freq = 50;
  if (print_freq < 1)
    print_freq = 1;
#endif

#ifdef REPARTITION
  // Need this to know, if this proc is a new joining or old initial process.
  int iter = laik_phase(hpcg_instance);
#endif


  if(A.repartition_me)
  {
    if (iter == 0)
    {
      // copy x to p for sparse MV operation
      CopyLaikVectorToLaikVector(x, p, A.mapping);
      TICK();
      ComputeSPMV_laik_ref(A, p, Ap);
      TOCK(t3); // Ap = A*p
      TICK();
      ComputeWAXPBY_laik_ref(nrow, 1.0, b, -1.0, Ap, r, A.mapping);
      TOCK(t2); // r = b - Ax (x stored in p)
      TICK();
      ComputeDotProduct_laik_ref(nrow, r, r, normr, t4, A.mapping);
      TOCK(t1);
      normr = sqrt(normr);

      // Record initial residual for convergence testing
      normr0 = normr;
    }
    else
    {
      // Joining procs need to get this from other proc
      laik_broadcast((void *)&normr0, (void *)&normr0, 1, laik_Double);
      laik_broadcast((void *)&normr, (void *)&normr, 1, laik_Double);
      laik_broadcast((void *)&rtz, (void *)&rtz, 1, laik_Double);
      // Values of the vector were recieved with a switchto before
    }
  }
  else
  {
    // Normal CG_ref call without repartitioning
    // copy x to p for sparse MV operation
    CopyLaikVectorToLaikVector(x, p, A.mapping);
    TICK();
    ComputeSPMV_laik_ref(A, p, Ap);
    TOCK(t3); // Ap = A*p
    TICK();
    ComputeWAXPBY_laik_ref(nrow, 1.0, b, -1.0, Ap, r, A.mapping);
    TOCK(t2); // r = b - Ax (x stored in p)
    TICK();
    ComputeDotProduct_laik_ref(nrow, r, r, normr, t4, A.mapping);
    TOCK(t1);
    normr = sqrt(normr);

    // Record initial residual for convergence testing
    normr0 = normr;
  }
  

#ifdef HPCG_DEBUG
  if (A.geom->rank == 0)
    HPCG_fout << "Initial Residual = " << normr << std::endl;
#endif

  int k = 1;

#ifdef REPARTITION
  if(A.repartition_me && laik_phase(hpcg_instance) > 0)
  {
    k = laik_phase(hpcg_instance); // should equal 11
    // assert(k == 11);
  }
#endif

  double t_before_start = 0;
  double t_after_start = 0;
  int k_before = 0;
  int k_after = 0;
  double local_time = 0;
  int max_iter2 = 100;
  for (; k <= max_iter2 && normr / normr0 > tolerance; k++)
  {
    printf("%dth iteration\n", k);
    if (k <= 50)
    {
      t_before_start = mytimer();
      k_before++;
    }
    else if (k > 50 && k <= 100)
    {
      t_after_start = mytimer();
      k_after++;
    }

    TICK();
    if (doPreconditioning)
      ComputeMG_laik_ref(A, r, z, k); // Apply preconditioner
    else
      ComputeWAXPBY_laik_ref(nrow, 1.0, r, 0.0, r, z, A.mapping); // copy r to z (no preconditioning)
    TOCK(t5); // Preconditioner apply time

   

    if (k == 1)
    {
      CopyLaikVectorToLaikVector(z, p, A.mapping);
      TOCK(t2); // Copy Mr to p
      TICK();
      ComputeDotProduct_laik_ref(nrow, r, z, rtz, t4, A.mapping);
      TOCK(t1); // rtz = r'*z
    }
    else
    {
      oldrtz = rtz;
      TICK();
      ComputeDotProduct_laik_ref(nrow, r, z, rtz, t4, A.mapping);
      TOCK(t1); // rtz = r'*z
      beta = rtz / oldrtz;
      TICK();
      ComputeWAXPBY_laik_ref(nrow, 1.0, z, beta, p, p, A.mapping);
      TOCK(t2); // p = beta*p + z
    }

    TICK();
    ComputeSPMV_laik_ref(A, p, Ap);
    TOCK(t3); // Ap = A*p
    TICK();
    ComputeDotProduct_laik_ref(nrow, p, Ap, pAp, t4, A.mapping);
    TOCK(t1); // alpha = p'*Ap
    alpha = rtz / pAp;
    TICK();
    ComputeWAXPBY_laik_ref(nrow, 1.0, x, alpha, p, x, A.mapping); // x = x + alpha*p
    ComputeWAXPBY_laik_ref(nrow, 1.0, r, -alpha, Ap, r, A.mapping);
    TOCK(t2); // r = r - alpha*Ap
    TICK();
    ComputeDotProduct_laik_ref(nrow, r, r, normr, t4, A.mapping);
    TOCK(t1);
    normr = sqrt(normr);
#ifdef HPCG_DEBUG
    if (A.geom->rank == 0 && (k % print_freq == 0 || k == max_iter))
      HPCG_fout << "Iteration = " << k << "   Scaled Residual = " << normr / normr0 << std::endl;
#endif
    niters = k;
    if (k <= 50)
    {
      local_time = mytimer() - t_before_start;
      t_it_before += local_time; // avg time of all iterations before repartitioning
    }
    else if (k > 50 && k <= 100)
    {
      local_time = mytimer() - t_after_start;
      t_it_after += local_time; // avg time of all iterations before repartitioning
    }
#ifdef REPARTITION
    // Repartitioning / Resizing of current world (group of proccesses) in the 10th iteration
    // For now, repartitioning is only done once
    if (k == 50 && A.repartition_me && !A.repartitioned)
    {
      TICK();
      A.repartition_me = false;             // do not repartition the matrix anymore
      uint32_t old_size = laik_size(world); // Old world size for output

      // allow resize of world and get new world
      Laik_Group *newworld = laik_allow_world_resize(hpcg_instance, k + 1);
      uint32_t new_size = laik_size(newworld); // Old world size for output

      // laik_finish_world_resize(hpcg_instance);
      if (newworld != world)
      {
        // Old partitionings for the x vector, which will be freed after re_switch_LaikVectors()
        Laik_Partitioning *old_local = A.local;
        Laik_Partitioning *old_ext = A.ext;

        // Assign new world and release old world
        laik_release_group(world);
        world = newworld;

        // Re-run setup functions and run partitioners for the new group
        repartition_SparseMatrix(A);
        nrow = A.localNumberOfRows; /* update local value after repartitioning */

        // Switch to the new partitioning on all Laik_data containers
        std::vector<Laik_Blob *> list{};
        /* Vectors in MG_data will be recursively handled in re_switch_LaikVectors */
        list.push_back(b);
        list.push_back(x);
        list.push_back(A.ptr_to_xexact); /* x_exact is out of scope but we need to switch this vector as well*/
        list.push_back(r);
        list.push_back(z);
        list.push_back(p);
        list.push_back(Ap);

        re_switch_LaikVectors(A, list);
        laik_free_partitioning(old_local);
        laik_free_partitioning(old_ext);

        // Need to send normr to new joining procs
        if(new_size > old_size)
        {
          laik_broadcast((void *)&normr0, (void *) &normr0, 1, laik_Double);
          laik_broadcast((void *)&normr, (void *) &normr, 1, laik_Double);
          laik_broadcast((void *)&rtz, (void *)&rtz, 1, laik_Double);
        }

        // Releasing old, not needed ressources in re_setup_problem(A);
        // Exit, if we got removed from the world
        if (laik_myid(world) < 0)
        {
          // DeleteMatrix(A); TODO fix seg fault
          DeleteCGData(data);
          DeleteLaikVector(A.ptr_to_xexact);
          DeleteLaikVector(x);
          DeleteLaikVector(b);

          HPCG_Finalize();
          laik_finalize(hpcg_instance);
          exit_hpcg_run("", false);
        }

        A.repartitioned = true;
        HPCG_fout << "REPARTIONING: Old world size [" << old_size << "] New world size [" << new_size << "]" << std::endl;
      }
      TOCK(t_rep);
    }
#endif // REPARTITION
  }

  // change according to test size TODO do it automatically
  int old_size = 12 - 1;
  int new_size = 16 - 1;

  if (laik_myid(world) <= old_size)
    t_it_before = (t_it_before / (double)k_before);

  t_it_after = (t_it_after / (double)k_after);

  // if(A.geom->rank == 0)
  {
    std::string result{"LAIK \t"};
    result += to_string(laik_myid(world)) + "\n";
    if (laik_myid(world) <= old_size)
    {
      result += "repartitioning time: " + to_string(t_rep) + " seconds\n";
      result += "avg before repartitioning: " + to_string(t_it_before) + " seconds\n";
    }

    result += "avg after repartitioning: " + to_string(t_it_after) + " seconds\n";
    std::cout << result;
  }

  // #ifdef REPARTITION
  //   laik_set_phase(hpcg_instance, 0, 0, 0);
  // #endif

  // Store times
  times[1] += t1; // dot product time
  times[2] += t2; // WAXPBY time
  times[3] += t3; // SPMV time
  times[4] += t4; // AllReduce time
  times[5] += t5; // preconditioner apply time
                  // #ifndef HPCG_NO_MPI
  //   times[6] += t6; // exchange halo time
  // #endif
  times[0] += mytimer() - t_begin; // Total time. All done...
  exit(1);
  // printf("Reference CG Timing Phase: %.5f seconds\n", times[0]);

  return 0;
}
#else
/*!
  Reference routine to compute an approximate solution to Ax = b

  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess; on exit: the new approximate solution
  @param[in]    max_iter  The maximum number of iterations to perform, even if tolerance is not met.
  @param[in]    tolerance The stopping criterion to assert convergence: if norm of residual is <= to tolerance.
  @param[out]   niters    The number of iterations actually performed.
  @param[out]   normr     The 2-norm of the residual vector after the last iteration.
  @param[out]   normr0    The 2-norm of the residual vector before the first iteration.
  @param[out]   times     The 7-element vector of the timing information accumulated during all of the iterations.
  @param[in]    doPreconditioning The flag to indicate whether the preconditioner should be invoked at each iteration.

  @return Returns zero on success and a non-zero value otherwise.

  @see CG()
*/
int CG_ref(const SparseMatrix &A, CGData &data, const Vector &b, Vector &x,
           const int max_iter, const double tolerance, int &niters, double &normr, double &normr0,
           double *times, bool doPreconditioning)
{

  double t_begin = mytimer(); // Start timing right away
  normr = 0.0;
  double rtz = 0.0, oldrtz = 0.0, alpha = 0.0, beta = 0.0, pAp = 0.0;

  double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;
  // #ifndef HPCG_NO_MPI
  //   double t6 = 0.0;
  // #endif

  local_int_t nrow = A.localNumberOfRows;

  Vector &r = data.r; // Residual vector
  Vector &z = data.z; // Preconditioned residual vector
  Vector &p = data.p; // Direction vector (in MPI mode ncol>=nrow)
  Vector &Ap = data.Ap;

  if (!doPreconditioning && A.geom->rank == 0)
    HPCG_fout << "WARNING: PERFORMING UNPRECONDITIONED ITERATIONS" << std::endl;

#ifdef HPCG_DEBUG
  int print_freq = 1;
  if (print_freq > 50)
    print_freq = 50;
  if (print_freq < 1)
    print_freq = 1;
#endif
  // p is of length ncols, copy x to p for sparse MV operation
  CopyVector(x, p);
  TICK();
  ComputeSPMV_ref(A, p, Ap);
  TOCK(t3); // Ap = A*p
  TICK();
  ComputeWAXPBY_ref(nrow, 1.0, b, -1.0, Ap, r);
  TOCK(t2); // r = b - Ax (x stored in p)
  TICK();
  ComputeDotProduct_ref(nrow, r, r, normr, t4);
  TOCK(t1);
  normr = sqrt(normr);
#ifdef HPCG_DEBUG
  if (A.geom->rank == 0)
    HPCG_fout << "Initial Residual = " << normr << std::endl;
#endif

  // Record initial residual for convergence testing
  normr0 = normr;

  // Start iterations

  for (int k = 1; k <= max_iter && normr / normr0 > tolerance; k++)
  {
    TICK();
    if (doPreconditioning)
      ComputeMG_ref(A, r, z); // Apply preconditioner
    else
      ComputeWAXPBY_ref(nrow, 1.0, r, 0.0, r, z); // copy r to z (no preconditioning)
    TOCK(t5);                                     // Preconditioner apply time

    if (k == 1)
    {
      CopyVector(z, p);
      TOCK(t2); // Copy Mr to p
      TICK();
      ComputeDotProduct_ref(nrow, r, z, rtz, t4);
      TOCK(t1); // rtz = r'*z
    }
    else
    {
      oldrtz = rtz;
      TICK();
      ComputeDotProduct_ref(nrow, r, z, rtz, t4);
      TOCK(t1); // rtz = r'*z
      beta = rtz / oldrtz;
      TICK();
      ComputeWAXPBY_ref(nrow, 1.0, z, beta, p, p);
      TOCK(t2); // p = beta*p + z
    }

    TICK();
    ComputeSPMV_ref(A, p, Ap);
    TOCK(t3); // Ap = A*p
    TICK();
    ComputeDotProduct_ref(nrow, p, Ap, pAp, t4);
    TOCK(t1); // alpha = p'*Ap
    alpha = rtz / pAp;
    TICK();
    ComputeWAXPBY_ref(nrow, 1.0, x, alpha, p, x); // x = x + alpha*p
    ComputeWAXPBY_ref(nrow, 1.0, r, -alpha, Ap, r);
    TOCK(t2); // r = r - alpha*Ap
    TICK();
    ComputeDotProduct_ref(nrow, r, r, normr, t4);
    TOCK(t1);
    normr = sqrt(normr);
#ifdef HPCG_DEBUG
    if (A.geom->rank == 0 && (k % print_freq == 0 || k == max_iter))
      HPCG_fout << "Iteration = " << k << "   Scaled Residual = " << normr / normr0 << std::endl;
#endif
    niters = k;
  }

  // Store times
  times[1] += t1; // dot product time
  times[2] += t2; // WAXPBY time
  times[3] += t3; // SPMV time
  times[4] += t4; // AllReduce time
  times[5] += t5; // preconditioner apply time
                  // #ifndef HPCG_NO_MPI
  //   times[6] += t6; // exchange halo time
  // #endif
  times[0] += mytimer() - t_begin; // Total time. All done...
  return 0;
}
#endif