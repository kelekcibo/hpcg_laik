
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
 @file ComputeProlongation_ref.cpp

 HPCG routine
 */

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "ComputeProlongation_ref.hpp"

/*!
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeProlongation_ref(const SparseMatrix & Af, Vector & xf, Laik_Blob * xf_blob) {

  bool use_laik = xf_blob != 0;
  double * xfv = xf.values;

  if(use_laik)
  {
    double *base;
    uint64_t count;
    laik_get_map_1d(xf_blob->values, 0, (void **)&base, &count);

    xfv = base;
  }

  double * xcv = Af.mgData->xc->values;
  local_int_t * f2c = Af.mgData->f2cOperator;
  local_int_t nc = Af.mgData->rc->localLength;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
// TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for (local_int_t i=0; i<nc; ++i)
  {
    if(use_laik)
      xfv[map_l2a(xf_blob->mapping, f2c[i], false)] += xcv[i]; // This loop is safe to vectorize
    else
      xfv[f2c[i]] += xcv[i]; // This loop is safe to vectorize
  
  }
  return 0;
}
