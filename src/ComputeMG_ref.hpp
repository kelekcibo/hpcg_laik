
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

int ComputeMG_ref(const SparseMatrix  & A, const Vector & r, Vector & x, Laik_Blob * x_blob);

#endif // COMPUTEMG_REF_HPP
