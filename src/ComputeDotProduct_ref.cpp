
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
 @file ComputeDotProduct_ref.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include "laik_instance.hpp"
#include "mytimer.hpp"
#endif
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
#include "ComputeDotProduct_ref.hpp"

/*!
  Routine to compute the dot product of two vectors where:

  This is the reference dot-product implementation.  It _CANNOT_ be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] x, y the input vectors
  @param[in] result a pointer to scalar value, on exit will contain result.
  @param[out] time_allreduce the time it took to perform the communication between processes

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct
*/
int ComputeDotProduct_ref(const local_int_t n, const Vector & x, const Vector & y,
    double & result, double & time_allreduce, Laik_Blob * x_blob, Laik_Blob * y_blob) {
  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);

  double local_result = 0.0;
  double * xv = x.values;
  double * yv = y.values;

  bool x_blob_active = false;
  bool y_blob_active = false;

  double *base_x;
  uint64_t count_x;
  double *base_y;
  uint64_t count_y;

  if (x_blob != 0)
  {
    x_blob_active = true;
    laik_get_map_1d(x_blob->values, 0, (void **)&base_x, &count_x);
    xv = base_x;
  }
  if (y_blob != 0)
  {
    y_blob_active = true;
    laik_get_map_1d(y_blob->values, 0, (void **)&base_y, &count_y);
    yv = base_y;
  }

  if (yv==xv) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for reduction (+:local_result)
#endif
    for (local_int_t i=0; i<n; i++) 
    {
      // As yv == xv, both blobs will be active in case
      if(x_blob_active && y_blob_active)
        local_result += xv[map_l2a(x_blob->mapping, i, false)] * xv[map_l2a(x_blob->mapping, i, false)];
      else
        local_result += xv[i] * xv[i];
    }
  } else {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for reduction (+:local_result)
#endif
    for (local_int_t i=0; i<n; i++)
    {
      if (x_blob_active && y_blob_active)
        local_result += xv[map_l2a(x_blob->mapping, i, false)] * yv[map_l2a(y_blob->mapping, i, false)];
      else if(x_blob_active && !y_blob_active)
        local_result += xv[map_l2a(x_blob->mapping, i, false)] * yv[i];
      else if(!x_blob_active && y_blob_active)
        local_result += xv[i] * yv[map_l2a(y_blob->mapping, i, false)];
      else
        local_result += xv[i] * yv[i];
    }
  }

#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to collect all partial sums
  double t0 = mytimer();
  double global_result = 0.0;
  laik_allreduce(&local_result, &global_result, 1, laik_Double, LAIK_RO_Sum);
  result = global_result;
  time_allreduce += mytimer() - t0;
#else
  time_allreduce += 0.0;
  result = local_result;
#endif

  return 0;
}
