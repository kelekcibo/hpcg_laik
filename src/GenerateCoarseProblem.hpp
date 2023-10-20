
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

#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"

void GenerateCoarseProblem(const SparseMatrix & A);

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
void GenerateCoarseProblem_repartition_ref(const SparseMatrix &A);
#endif // HPCG_NO_LAIK
#endif // REPARTITION

#endif // GENERATECOARSEPROBLEM_HPP
