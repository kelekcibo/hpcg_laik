
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

#ifndef GENERATECOARSEPROBLEM_HPP
#define GENERATECOARSEPROBLEM_HPP

#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "laik_instance.hpp"
#include "SparseMatrix.hpp"

void GenerateCoarseProblem(const SparseMatrix & A);
#endif // GENERATECOARSEPROBLEM_HPP
