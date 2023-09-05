
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

#ifndef COMPUTEPROLONGATION_REF_HPP
#define COMPUTEPROLONGATION_REF_HPP

#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik_instance.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#ifdef USE_LAIK
int ComputeProlongation_ref(const SparseMatrix &Af, Laik_Blob *xf_blob);
#else
int ComputeProlongation_ref(const SparseMatrix & Af, Vector & xf);
#endif
#endif // COMPUTEPROLONGATION_REF_HPP
