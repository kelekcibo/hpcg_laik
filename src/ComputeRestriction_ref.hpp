
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
#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik_instance.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#ifdef USE_LAIK
int ComputeRestriction_ref(const SparseMatrix &A, const Laik_Blob *rf);
#else
int ComputeRestriction_ref(const SparseMatrix & A, const Vector & rf);
#endif
#endif // COMPUTERESTRICTION_REF_HPP
