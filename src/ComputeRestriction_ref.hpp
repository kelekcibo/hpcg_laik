
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

#ifndef COMPUTERESTRICTION_REF_HPP
#define COMPUTERESTRICTION_REF_HPP
#include "laik/hpcg_laik.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"
int ComputeRestriction_laik_ref(const SparseMatrix &A, const Laik_Blob *rf, int k);
int ComputeRestriction_ref(const SparseMatrix & A, const Vector & rf);
#endif // COMPUTERESTRICTION_REF_HPP
