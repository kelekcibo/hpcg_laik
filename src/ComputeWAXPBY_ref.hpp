
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

#ifndef COMPUTEWAXPBY_REF_HPP
#define COMPUTEWAXPBY_REF_HPP

#include "laik/hpcg_laik.hpp"
#include "Vector.hpp"

int ComputeWAXPBY_laik_ref(const local_int_t n, const double alpha, const Laik_Blob *x,
                           const double beta, const Laik_Blob *y, const Laik_Blob *w);
int ComputeWAXPBY_ref(const local_int_t n, const double alpha, const Vector & x,
    const double beta, const Vector & y, Vector & w);
#endif // COMPUTEWAXPBY_REF_HPP
