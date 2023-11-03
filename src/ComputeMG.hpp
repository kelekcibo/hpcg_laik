
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

#ifndef COMPUTEMG_HPP
#define COMPUTEMG_HPP

#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

#ifndef HPCG_NO_MPI
int ComputeMG_laik(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
#else
int ComputeMG(const SparseMatrix  & A, const Vector & r, Vector & x);
#endif

#endif // COMPUTEMG_HPP
