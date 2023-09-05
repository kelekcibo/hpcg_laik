
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
#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "Vector.hpp"
#include "laik_instance.hpp"

#ifdef USE_LAIK
int ComputeDotProduct_ref(const local_int_t n, const Laik_Blob *x, const Laik_Blob *y,
                          double &result, double &time_allreduce, L2A_map * mapping);
#else
int ComputeDotProduct_ref(const local_int_t n, const Vector & x, const Vector & y,
    double & result, double & time_allreduce);
#endif

#endif // COMPUTEDOTPRODUCT_REF_HPP
