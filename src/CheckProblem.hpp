
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

#ifndef CHECKPROBLEM_HPP
#define CHECKPROBLEM_HPP

#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void CheckProblem(SparseMatrix & A, Vector * b, Vector * x, Vector * xexact);

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
void CheckProblem_repartition(SparseMatrix &A, Vector *b, Vector *x, Vector *xexact);
#endif // HPCG_NO_LAIK
#endif // REPARTITION

#endif // CHECKPROBLEM_HPP
