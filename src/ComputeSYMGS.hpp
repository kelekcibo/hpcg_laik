
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

#ifndef COMPUTESYMGS_HPP
#define COMPUTESYMGS_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

#ifdef USE_LAIK
int ComputeSYMGS(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
#else
int ComputeSYMGS(const SparseMatrix &A, const Vector &r, Vector &x);
#endif
#endif // COMPUTESYMGS_HPP
