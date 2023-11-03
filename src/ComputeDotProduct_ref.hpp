
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

#ifndef COMPUTEDOTPRODUCT_REF_HPP
#define COMPUTEDOTPRODUCT_REF_HPP

#include "laik/hpcg_laik.hpp"
#include "Vector.hpp"

int ComputeDotProduct_laik_ref(const local_int_t n, const Laik_Blob *x, const Laik_Blob *y,
                          double &result, double &time_allreduce, L2A_map * mapping);

int ComputeDotProduct_ref(const local_int_t n, const Vector & x, const Vector & y,
    double & result, double & time_allreduce);

#endif // COMPUTEDOTPRODUCT_REF_HPP
