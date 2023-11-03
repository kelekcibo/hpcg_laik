
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

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
void SetupHalo_repartition_ref(SparseMatrix &A);
#endif
#endif

#endif // SETUPHALO_REF_HPP
