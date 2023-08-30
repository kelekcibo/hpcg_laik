
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
#include "laik_instance.hpp"
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
int ComputeWAXPBY_ref(const local_int_t n, const double alpha, const Vector & x,
    const double beta, const Vector & y, Vector & w, Laik_Blob * x_blob, Laik_Blob * y_blob, Laik_Blob * w_blob) {

  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);

  bool x_blob_active = false;
  bool y_blob_active = false;
  bool w_blob_active = false;

  if(x_blob != 0)
    x_blob_active = true;
  if (y_blob != 0)
    y_blob_active = true;
  if (w_blob != 0)
    w_blob_active = true;

  const double * const xv = x.values;
  const double * const yv = y.values;
  double * const wv = w.values;

  if (alpha==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = xv[i] + beta * yv[i];
  } else if (beta==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + yv[i];
  } else  {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + beta * yv[i];
  }

  return 0;
}
