
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

#ifndef COMPUTESYMGS_REF_HPP
#define COMPUTESYMGS_REF_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

#ifdef USE_LAIK
int ComputeSYMGS_ref(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
#else
int ComputeSYMGS_ref(const SparseMatrix &A, const Vector &r, Vector &x);
#endif

#endif // COMPUTESYMGS_REF_HPP
