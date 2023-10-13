
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

#ifndef COMPUTEPROLONGATION_REF_HPP
#define COMPUTEPROLONGATION_REF_HPP

#include "laik/hpcg_laik.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int ComputeProlongation_laik_ref(const SparseMatrix &Af, Laik_Blob *xf_blob);
int ComputeProlongation_ref(const SparseMatrix & Af, Vector & xf);
#endif // COMPUTEPROLONGATION_REF_HPP
