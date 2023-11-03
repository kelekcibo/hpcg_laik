
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
 @file ComputeRestriction_ref.cpp

 HPCG routine
 */


#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "laik/hpcg_laik.hpp"
#include "ComputeRestriction_ref.hpp"

#ifndef HPCG_NO_MPI
/*!
  Routine to compute the coarse residual vector.

  @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf, the fine grid matrix-vector product and mgData->rc the coarse residual vector.
  @param[in]    rf - Fine grid RHS.


  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeRestriction_laik_ref(const SparseMatrix &A, const Laik_Blob *rf, int k)
{
  double *rfv;
  double *rcv;
  double *Axfv;

  laik_get_map_1d(rf->values, 0, (void **)&rfv, 0);
  laik_get_map_1d(A.mgData->rc_blob->values, 0, (void **)&rcv, 0);
  laik_get_map_1d(A.mgData->Axf_blob->values, 0, (void **)&Axfv, 0);

  local_int_t *f2c = A.mgData->f2cOperator;
  local_int_t nc = A.mgData->rc_blob->localLength;

  // rc vector is for next layer, thus need mapping from next level matrix
  assert(A.Ac != NULL); 
  L2A_map * mapping_rc_blob = A.Ac->mapping;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for (local_int_t i = 0; i < nc; ++i)
  {
    local_int_t j = map_l2a_x(A.mapping, f2c[i], false);
    rcv[map_l2a_x(mapping_rc_blob, i, false)] = rfv[j] - Axfv[j]; 
  }

  return 0;
}
#else
/*!
  Routine to compute the coarse residual vector.

  @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf, the fine grid matrix-vector product and mgData->rc the coarse residual vector.
  @param[in]    rf - Fine grid RHS.


  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeRestriction_ref(const SparseMatrix & A, const Vector & rf) {

  double * Axfv = A.mgData->Axf->values;
  double * rfv = rf.values;
  double * rcv = A.mgData->rc->values;
  local_int_t * f2c = A.mgData->f2cOperator;
  local_int_t nc = A.mgData->rc->localLength;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for (local_int_t i=0; i<nc; ++i) rcv[i] = rfv[f2c[i]] - Axfv[f2c[i]];

  return 0;
}
#endif