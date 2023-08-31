
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

#include "ComputeMG_ref.hpp"
#include "ComputeSYMGS_ref.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeRestriction_ref.hpp"
#include "ComputeProlongation_ref.hpp"
#include <cassert>
#include <iostream>

/*!

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG
*/
int ComputeMG_ref(const SparseMatrix & A, const Vector & r, Vector & x, Laik_Blob * x_blob)
{
  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

  bool use_laik = x_blob != 0;

  if (use_laik)
    assert(x_blob->localLength == A.localNumberOfRows);


  ZeroVector(x); // initialize x to zero

  if(use_laik)
    ZeroLaikVector(x_blob);

  int ierr = 0;
  if (A.mgData!=0) { // Go to next coarse level if defined
    int numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps;
    for (int i=0; i< numberOfPresmootherSteps; ++i)
    {
      if(use_laik)
        ierr += ComputeSYMGS_ref(A, r, x, x_blob);
      else
        ierr += ComputeSYMGS_ref(A, r, x, NULL);
    }

    if (ierr!=0) return ierr;

    if(use_laik)
      ierr = ComputeSPMV_ref(A, x, *A.mgData->Axf, x_blob); if (ierr!=0) return ierr;
    else
      ierr = ComputeSPMV_ref(A, x, *A.mgData->Axf, NULL); if (ierr!=0) return ierr;

    exit(1); // DELETE.

    // Perform restriction operation using simple injection
    ierr = ComputeRestriction_ref(A, r);  if (ierr!=0) return ierr;

    // printf("Matrix layer (%d)\tA locCol (%d), xc locLeng (%d)\n", A.Ac->level, A.Ac->localNumberOfColumns, A.mgData->xc->localLength);
    
    // didnt implement xc laik blob yet TODO
    // if(use_laik)
    //   ierr = ComputeMG_ref(*A.Ac,*A.mgData->rc, *A.mgData->xc, A.mgData->xc_blob);  if (ierr!=0) return ierr;
    // else

    ierr = ComputeMG_ref(*A.Ac,*A.mgData->rc, *A.mgData->xc, NULL);  if (ierr!=0) return ierr;


    if(use_laik)
      ierr = ComputeProlongation_ref(A, x, x_blob);  if (ierr!=0) return ierr;
    else
      ierr = ComputeProlongation_ref(A, x, NULL);  if (ierr!=0) return ierr;
  
    int numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps;
    for (int i=0; i< numberOfPostsmootherSteps; ++i)
    {
      if(use_laik)
        ierr += ComputeSYMGS_ref(A, r, x, x_blob);
      else
        ierr += ComputeSYMGS_ref(A, r, x, NULL);
    } 
    if (ierr!=0) return ierr;
  }
  else {
    if(use_laik)
      ierr = ComputeSYMGS_ref(A, r, x, x_blob);
    else
      ierr = ComputeSYMGS_ref(A, r, x, NULL);

    if (ierr!=0) return ierr;
  }
  return 0;
}
