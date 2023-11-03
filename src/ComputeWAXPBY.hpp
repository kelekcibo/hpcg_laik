
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

#ifndef COMPUTEWAXPBY_HPP
#define COMPUTEWAXPBY_HPP

#include "laik/hpcg_laik.hpp"
#include "Vector.hpp"
int ComputeWAXPBY_laik(const local_int_t n, const double alpha, const Laik_Blob *x,
                       const double beta, const Laik_Blob *y, Laik_Blob *w, bool &isOptimized, L2A_map *mapping);
int ComputeWAXPBY(const local_int_t n, const double alpha, const Vector &x,
                  const double beta, const Vector &y, Vector &w, bool &isOptimized);
#endif // COMPUTEWAXPBY_HPP
