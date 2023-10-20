
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

#ifndef GENERATEPROBLEM_REF_HPP
#define GENERATEPROBLEM_REF_HPP

#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void GenerateProblem_ref(SparseMatrix & A, Vector * b, Vector * x, Vector * xexact);

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
void GenerateProblem_repartition_ref(SparseMatrix &A, Vector *b, Vector *x, Vector *xexact);
#endif // HPCG_NO_LAIK
#endif // REPARTITION

#endif // GENERATEPROBLEM_REF_HPP
