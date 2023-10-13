
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
 @file ComputeWAXPBY_ref.cpp

 HPCG routine
 */

#include "ComputeWAXPBY_ref.hpp"
#include "laik/hpcg_laik.hpp"
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This is the reference WAXPBY impmentation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY
*/
int ComputeWAXPBY_laik_ref(const local_int_t n, const double alpha, const Laik_Blob *x,
                      const double beta, const Laik_Blob *y, const Laik_Blob *w, L2A_map * mapping)
{

  assert(x->localLength == n); // Test vector lengths
  assert(y->localLength == n);
  assert(w->localLength == n);
  assert(mapping->localNumberOfRows == n);

  const double * xv;
  const double * yv;
  double * wv;
  laik_get_map_1d(x->values, 0, (void **)&xv, 0);
  laik_get_map_1d(y->values, 0, (void **)&yv, 0);
  laik_get_map_1d(w->values, 0, (void **)&wv, 0);

  if (alpha == 1.0)
  {
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
    for (local_int_t i = 0; i < n; i++)
    {
      allocation_int_t j = map_l2a(mapping, i, false);
      wv[j] = xv[j] + beta * yv[j];
    }
  }
  else if (beta == 1.0)
  {
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
    for (local_int_t i = 0; i < n; i++)
    {
      allocation_int_t j = map_l2a(mapping, i, false);
      wv[j] = alpha * xv[j] + yv[j];
    }
  }
  else
  {
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
    for (local_int_t i = 0; i < n; i++)
    {
      allocation_int_t j = map_l2a(mapping, i, false);
      wv[j] = alpha * xv[j] + beta * yv[j];
    }
  }

  return 0;
}

/*!
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This is the reference WAXPBY impmentation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY
*/
int ComputeWAXPBY_ref(const local_int_t n, const double alpha, const Vector &x,
                      const double beta, const Vector &y, Vector &w)
{

  assert(x.localLength >= n); // Test vector lengths
  assert(y.localLength >= n);

  const double *const xv = x.values;
  const double *const yv = y.values;
  double *const wv = w.values;

  if (alpha == 1.0)
  {
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
    for (local_int_t i = 0; i < n; i++)
      wv[i] = xv[i] + beta * yv[i];
  }
  else if (beta == 1.0)
  {
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
    for (local_int_t i = 0; i < n; i++)
      wv[i] = alpha * xv[i] + yv[i];
  }
  else
  {
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
    for (local_int_t i = 0; i < n; i++)
      wv[i] = alpha * xv[i] + beta * yv[i];
  }

  return 0;
}