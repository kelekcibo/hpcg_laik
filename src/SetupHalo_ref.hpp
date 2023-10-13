
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

#ifndef SETUPHALO_REF_HPP
#define SETUPHALO_REF_HPP
#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"

void SetupHalo_ref(SparseMatrix & A);

#endif // SETUPHALO_REF_HPP
