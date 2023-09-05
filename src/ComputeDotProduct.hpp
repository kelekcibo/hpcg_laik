
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

#ifndef COMPUTEDOTPRODUCT_HPP
#define COMPUTEDOTPRODUCT_HPP

#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik_instance.hpp"
#include "Vector.hpp"

#ifdef USE_LAIK
int ComputeDotProduct(const local_int_t n, const Laik_Blob *x, const Laik_Blob *y,
                      double &result, double &time_allreduce, bool &isOptimized, L2A_map *mapping);
#else
int ComputeDotProduct(const local_int_t n, const Vector &x, const Vector &y,
                      double &result, double &time_allreduce, bool &isOptimized);
#endif
#endif // COMPUTEDOTPRODUCT_HPP
