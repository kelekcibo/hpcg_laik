
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
#include "Vector.hpp"
#include "laik_instance.hpp"
int ComputeDotProduct(const local_int_t n, const Vector & x, const Vector & y,
    double & result, double & time_allreduce, bool & isOptimized, Laik_Blob * x_blob, Laik_Blob * y_blob);

#endif // COMPUTEDOTPRODUCT_HPP
