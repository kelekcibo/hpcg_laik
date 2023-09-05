
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

#ifndef USE_LAIK
#define USE_LAIK
#endif
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void CheckProblem(SparseMatrix & A, Vector * b, Vector * x, Vector * xexact);
#endif // CHECKPROBLEM_HPP
