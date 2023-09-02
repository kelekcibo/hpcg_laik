
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

#ifndef COMPUTESPMV_HPP
#define COMPUTESPMV_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "laik_instance.hpp"

#ifdef USE_LAIK
int ComputeSPMV(const SparseMatrix &A, Laik_Blob *x_blob, Laik_Blob *y_blob);
#else
int ComputeSPMV(const SparseMatrix &A, Vector &x, Vector &y);
#endif

#endif  // COMPUTESPMV_HPP
