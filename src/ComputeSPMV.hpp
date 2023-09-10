
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

#include "laik_instance.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int ComputeSPMV_laik(const SparseMatrix &A, Laik_Blob *x_blob, Laik_Blob *y_blob);
int ComputeSPMV(const SparseMatrix &A, Vector &x, Vector &y);

#endif  // COMPUTESPMV_HPP
