
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

#ifndef SETUPHALO_HPP
#define SETUPHALO_HPP
#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik_instance.hpp"
#include "SparseMatrix.hpp"

void SetupHalo(SparseMatrix & A);

#endif // SETUPHALO_HPP
