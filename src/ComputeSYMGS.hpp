
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

#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

int ComputeSYMGS_laik(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
int ComputeSYMGS(const SparseMatrix &A, const Vector &r, Vector &x);
#endif // COMPUTESYMGS_HPP
