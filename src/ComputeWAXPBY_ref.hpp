
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
#include "Vector.hpp"
#include "laik_instance.hpp"
int ComputeWAXPBY_ref(const local_int_t n, const double alpha, const Vector & x,
    const double beta, const Vector & y, Vector & w, Laik_Blob * x_blob, Laik_Blob * y_blob, Laik_Blob * w_blob);
#endif // COMPUTEWAXPBY_REF_HPP
