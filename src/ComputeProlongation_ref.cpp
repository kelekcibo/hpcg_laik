
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

#include "laik/hpcg_laik.hpp"
#include "ComputeProlongation_ref.hpp"

#ifndef HPCG_NO_MPI

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
int ComputeProlongation_laik_repartition_ref(const SparseMatrix &Af, Laik_Blob *xf)
{

  double *xfv;
  double *xcv;

  laik_get_map_1d(xf->values, 0, (void **)&xfv, 0);
  laik_get_map_1d(Af.mgData->xc_blob->values, 0, (void **)&xcv, 0);

  local_int_t *f2c;
  laik_get_map_1d(Af.mgData->f2cOperator_d, 0, (void **)&f2c, 0);

  local_int_t nc = Af.mgData->rc_blob->localLength;

  // xc vector is for next layer, thus need mapping from next level matrix
  assert(Af.Ac != NULL);
  L2A_map *mapping_xc_blob = Af.Ac->mapping;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  // TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for (local_int_t i = 0; i < nc; ++i)
    xfv[map_l2a_x(Af.mapping, f2c[map_l2a_A(Af, i)], false)] += xcv[map_l2a_x(mapping_xc_blob, i, false)]; // This loop is safe to vectorize

  return 0;
}
#endif
#endif // HPCG_NO_LAIK

/*!
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeProlongation_laik_ref(const SparseMatrix & Af, Laik_Blob * xf) {

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
  return ComputeProlongation_laik_repartition_ref(Af, xf);
#endif
#endif

  double * xfv;
  double *xcv;
  
  laik_get_map_1d(xf->values, 0, (void **)&xfv, 0);
  laik_get_map_1d(Af.mgData->xc_blob->values, 0, (void **)&xcv, 0);

  local_int_t * f2c = Af.mgData->f2cOperator;
  local_int_t nc = Af.mgData->rc_blob->localLength;

  // xc vector is for next layer, thus need mapping from next level matrix
  assert(Af.Ac != NULL);
  L2A_map *mapping_xc_blob = Af.Ac->mapping;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
// TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for (local_int_t i=0; i<nc; ++i)
    xfv[map_l2a_x(Af.mapping, f2c[i], false)] += xcv[map_l2a_x(mapping_xc_blob, i, false)]; // This loop is safe to vectorize

  return 0;
}

#else
/*!
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeProlongation_ref(const SparseMatrix &Af, Vector &xf)
{

  double *xfv = xf.values;
  double *xcv = Af.mgData->xc->values;
  local_int_t *f2c = Af.mgData->f2cOperator;
  local_int_t nc = Af.mgData->rc->localLength;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  // TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for (local_int_t i = 0; i < nc; ++i)
    xfv[f2c[i]] += xcv[i]; // This loop is safe to vectorize

  return 0;
}
#endif