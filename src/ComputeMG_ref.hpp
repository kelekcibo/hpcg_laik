
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

#ifndef COMPUTEMG_REF_HPP
#define COMPUTEMG_REF_HPP

// #ifndef HPCG_NO_MPI
#include "laik/hpcg_laik.hpp"
// #endif
#include "SparseMatrix.hpp"
#include "Vector.hpp"


int ComputeMG_laik_ref(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
int ComputeMG_ref(const SparseMatrix &A, const Vector &r, Vector &x);

#endif // COMPUTEMG_REF_HPP
