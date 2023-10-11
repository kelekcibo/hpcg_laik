
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

#include "laik_instance.hpp"
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

  #ifdef REPARTITION
  int laik_iter = laik_get_iteration(hpcg_instance);
  #endif

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
// Old processes are in the 2nd iteration. Skip this part
  if(laik_iter < 2)
  {
#endif
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
#ifdef HPCG_DEBUG
  if (A.geom->rank == 0)
    HPCG_fout << "Initial Residual = " << normr << std::endl;
#endif

  // Record initial residual for convergence testing
  normr0 = normr; /* TODO. If new joining processes join and skip first part, this means they need to get normr0 from old procs */
  #ifdef REPARTITION
  }
  else
  {
    laik_broadcast(&normr, &normr0, 1, laik_Double); /* Get normr0 by proc 0 */
  }
  #endif


  std::string debug_str{""};

  // Start iterations
  int k = 1;
  // For new joining processes
  if (laik_iter > 1)
    k = laik_iter;

  for (; k <= max_iter && normr / normr0 > tolerance; k++)
  {

#ifdef REPARTITION

  if(k == 1)
  {
    // Repartitioning / Resizing of current world (group of proccesses)
    laik_set_iteration(hpcg_instance, k); /* Current iteration */

    // allow resize of world and get new world
    Laik_Group *newworld = laik_allow_world_resize(hpcg_instance, k);

    int shrink_count = 1;
    int plist[1];
    plist[0] = 1; // remove proc 1 as test with config: size = 2

    debug_str += "Before shrinking: size " + std::to_string(laik_size(newworld)) + " (id " + std::to_string(laik_myid(newworld)) + ")\n";

    newworld = laik_new_shrinked_group(world, shrink_count, plist);

    debug_str += "After shrinking: size " + std::to_string(laik_size(newworld)) + " (id " + std::to_string(laik_myid(newworld)) + ")";

    // exit_hpcg_run(debug_str.c_str());

    laik_finish_world_resize(hpcg_instance);

    if (newworld != world)
    {

      // Assign new world and release old world
      laik_release_group(world);
      world = newworld;

      // Re-run setup functions and run partitioners for the new group
      re_setup_problem(A);
      // exit_hpcg_run("RE-SETUP PROBLEM WORKS (without crashing)!");

      // TODO. Switch to the new partitioning on all Laik_data containers
      std::vector<Laik_Blob *> list{};

      list.push_back(b);
      list.push_back(x);
      list.push_back(r);
      list.push_back(z);
      list.push_back(p);
      list.push_back(Ap);

      /* Send normr0 to new procs, old procs have already this value */
      laik_broadcast(&normr0, &normr0, 1, laik_Double); /* Get normr0 by proc 0 */

      /* Vectors in MG_data will be recursively handled in re_switch_LaikVectors */
      re_switch_LaikVectors(A, list);


      exit_hpcg_run("RE_SWITCH_LAIKVECTORS PROBLEM WORKS (without crashing)!");

      // Releasing old, not needed ressources in re_setup_problem(A);
      // Exit, if we got removed from the world
      if (laik_myid(world) < 0)
      {
        DeleteMatrix_repartition(A, true);
        DeleteCGData(data);
        DeleteLaikVector(x);
        DeleteLaikVector(b);

        laik_finalize(hpcg_instance);
        exit(0);
        // return 0;
      }

      debug_str += "IN IF ";
      exit_hpcg_run(debug_str.c_str());
    }
  }

#endif

    TICK();
    if (doPreconditioning)
      ComputeMG_laik_ref(A, r, z); // Apply preconditioner
    else
      ComputeWAXPBY_laik_ref(nrow, 1.0, r, 0.0, r, z, A.mapping); // copy r to z (no preconditioning)
    TOCK(t5);                                     // Preconditioner apply time

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
  // printf("Reference CG Timing Phase: %.5f seconds\n", times[0]);

  return 0;
}

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
