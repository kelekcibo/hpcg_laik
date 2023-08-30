
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
#include "Vector.hpp"
#include "SparseMatrix.hpp"

// #ifndef HPCG_NO_MPI
#include "laik_instance.hpp"
// #endif

int ComputeSPMV_ref( const SparseMatrix & A, Vector  & x, Vector & y, Laik_Blob * x_blob);

#endif  // COMPUTESPMV_REF_HPP
