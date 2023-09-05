
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

#ifndef COMPUTESPMV_REF_HPP
#define COMPUTESPMV_REF_HPP

#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik_instance.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"


#ifdef USE_LAIK
int ComputeSPMV_ref(const SparseMatrix &A, Laik_Blob *x, Laik_Blob *y);
#else
int ComputeSPMV_ref(const SparseMatrix &A, Vector &x, Vector &y);
#endif

#endif  // COMPUTESPMV_REF_HPP
