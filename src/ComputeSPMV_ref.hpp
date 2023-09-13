
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

#include "laik_instance.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int ComputeSPMV_laik_ref(const SparseMatrix &A, Laik_Blob *x, Laik_Blob *y);
int ComputeSPMV_ref(const SparseMatrix &A, Vector &x, Vector &y);

#endif  // COMPUTESPMV_REF_HPP
