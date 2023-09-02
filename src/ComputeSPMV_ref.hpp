
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
#include "Vector.hpp"
#include "SparseMatrix.hpp"

// #ifndef HPCG_NO_MPI
#include "laik_instance.hpp"
// #endif

#ifdef USE_LAIK
int ComputeSPMV_ref(const SparseMatrix &A, Laik_Blob *x, Laik_Blob *y);
#else
int ComputeSPMV_ref(const SparseMatrix &A, Vector &x, Vector &y);
#endif

#endif  // COMPUTESPMV_REF_HPP
