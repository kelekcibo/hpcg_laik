
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

#ifndef COMPUTERESIDUAL_HPP
#define COMPUTERESIDUAL_HPP
#include "Vector.hpp"
#ifdef USE_LAIK
int ComputeResidual(const local_int_t n, const Laik_Blob *v1, const Laik_Blob *v2, double &residual, L2A_map *mapping);
#else
int ComputeResidual(const local_int_t n, const Vector &v1, const Vector &v2, double &residual);
#endif
#endif // COMPUTERESIDUAL_HPP
