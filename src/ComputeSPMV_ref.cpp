
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
 @file ComputeSPMV_ref.cpp

 HPCG routine
 */

#include "ComputeSPMV_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#include "laik_instance.hpp"
#include <iostream>
#include <cstdlib>
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference SPMV implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV
*/
int ComputeSPMV_ref(const SparseMatrix &A, Vector &x, Vector &y, Laik_Blob *x_blob)
{

  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);
  assert(x_blob->localLength >= A.localNumberOfRows);

#ifndef HPCG_NO_MPI

  bool use_laik = x_blob != 0;
  
  double *base;
  uint64_t count;
  if (use_laik)
    laik_switchto_partitioning(x_blob->values, x_blob->x_ext, LAIK_DF_Preserve, LAIK_RO_None);
  else
    ExchangeHalo(A, x);
#endif

  // const double * const xv = x.values;
  double * xv = x.values;
  if (use_laik)
  {
    laik_get_map_1d(x_blob->values, 0, (void **)&base, &count);
    xv = base;
  }

  double * const yv = y.values;
  const local_int_t nrow = A.localNumberOfRows;
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< nrow; i++)  {
    double sum = 0.0;
    const double * const cur_vals = A.matrixValues[i];
    const local_int_t * const cur_inds = A.mtxIndL[i];
    const int cur_nnz = A.nonzerosInRow[i];

    for (int j=0; j< cur_nnz; j++)
    {
      if (use_laik)
        sum += cur_vals[j] * base[map_l2a(x_blob->mapping, cur_inds[j], true)];
      else
        sum += cur_vals[j] * xv[cur_inds[j]];
    }

    yv[i] = sum;
  }

#ifndef HPCG_NO_MPI
  if (use_laik)
    laik_switchto_partitioning(x_blob->values, x_blob->x_local, LAIK_DF_None, LAIK_RO_None);
#endif

  return 0;
}
