
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
 @file ComputeSYMGS_ref.cpp

 HPCG routine
 */

#include <cassert>
#include <iostream>

#include "laik/hpcg_laik.hpp"
#include "ComputeMG_ref.hpp"
#include "ComputeSYMGS_ref.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeRestriction_ref.hpp"
#include "ComputeProlongation_ref.hpp"

#ifndef HPCG_NO_LAIK
/*!

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG
*/
int ComputeMG_laik_ref(const SparseMatrix &A, const Laik_Blob * r, Laik_Blob * x, int k)
{
  assert(x->localLength == A.localNumberOfRows);
  assert(x->localLength == r->localLength);

  ZeroLaikVector(x); // initialize x to zero

  int ierr = 0;
  if (A.mgData != 0)
  { // Go to next coarse level if defined
    int numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps;
    for (int i = 0; i < numberOfPresmootherSteps; ++i) ierr += ComputeSYMGS_laik_ref(A, r, x);
    if (k == 11)
    {
      std::string debug{"\x1B[33mCheckpoint END\x1B[0m"};
      exit_hpcg_run(debug.data(), false);
    }
    if (ierr != 0)
      return ierr;
    ierr = ComputeSPMV_laik_ref(A, x, A.mgData->Axf_blob);
    if (ierr != 0)
      return ierr;

    // Perform restriction operation using simple injection
    ierr = ComputeRestriction_laik_ref(A, r, k);
    if (ierr != 0)
      return ierr;

    ierr = ComputeMG_laik_ref(*A.Ac, A.mgData->rc_blob, A.mgData->xc_blob, k);
    if (ierr != 0)
      return ierr;
    ierr = ComputeProlongation_laik_ref(A, x);
    if (ierr != 0)
      return ierr;
    int numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps;
    for (int i = 0; i < numberOfPostsmootherSteps; ++i)
      ierr += ComputeSYMGS_laik_ref(A, r, x);
    if (ierr != 0)
      return ierr;
  }
  else
  {
    ierr = ComputeSYMGS_laik_ref(A, r, x);
    if (ierr != 0)
      return ierr;
  }
  return 0;
}
#else
/*!

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG
*/
int ComputeMG_ref(const SparseMatrix &A, const Vector &r, Vector &x)
{
  assert(x.localLength == A.localNumberOfColumns); // Make sure x contain space for halo values

  ZeroVector(x); // initialize x to zero

  int ierr = 0;
  if (A.mgData != 0)
  { // Go to next coarse level if defined
    int numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps;
    for (int i = 0; i < numberOfPresmootherSteps; ++i)
      ierr += ComputeSYMGS_ref(A, r, x);
    if (ierr != 0)
      return ierr;
    ierr = ComputeSPMV_ref(A, x, *A.mgData->Axf);
    if (ierr != 0)
      return ierr;
    // Perform restriction operation using simple injection
    ierr = ComputeRestriction_ref(A, r);
    if (ierr != 0)
      return ierr;
    ierr = ComputeMG_ref(*A.Ac, *A.mgData->rc, *A.mgData->xc);
    if (ierr != 0)
      return ierr;
    ierr = ComputeProlongation_ref(A, x);
    if (ierr != 0)
      return ierr;
    int numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps;
    for (int i = 0; i < numberOfPostsmootherSteps; ++i)
      ierr += ComputeSYMGS_ref(A, r, x);
    if (ierr != 0)
      return ierr;
  }
  else
  {
    ierr = ComputeSYMGS_ref(A, r, x);
    if (ierr != 0)
      return ierr;
  }
  return 0;
}
#endif