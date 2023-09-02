
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
#include "SparseMatrix.hpp"
#include "Vector.hpp"

// #ifndef HPCG_NO_MPI
#include "laik_instance.hpp"
// #endif

#ifdef USE_LAIK
int ComputeMG_ref(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
#else
int ComputeMG_ref(const SparseMatrix &A, const Vector &r, Vector &x);
#endif

#endif // COMPUTEMG_REF_HPP
