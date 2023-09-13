
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
 
#include "laik_instance.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

int ComputeMG_laik(const SparseMatrix &A, const Laik_Blob *r, Laik_Blob *x);
int ComputeMG(const SparseMatrix  & A, const Vector & r, Vector & x);

#endif // COMPUTEMG_HPP
