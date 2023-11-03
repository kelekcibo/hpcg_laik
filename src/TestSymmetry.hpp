
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

/*!
 @file TestSymmetry.hpp

 HPCG data structures for symmetry testing
 */

#ifndef TESTSYMMETRY_HPP
#define TESTSYMMETRY_HPP

#include "laik/hpcg_laik.hpp"
#include "hpcg.hpp"
#include "SparseMatrix.hpp"
#include "CGData.hpp"

struct TestSymmetryData_STRUCT {
  double depsym_spmv;  //!< departure from symmetry for the SPMV kernel
  double depsym_mg; //!< departure from symmetry for the MG kernel
  int    count_fail;   //!< number of failures in the symmetry tests
};
typedef struct TestSymmetryData_STRUCT TestSymmetryData;

#ifndef HPCG_NO_MPI
extern int TestSymmetry_laik(SparseMatrix &A, Laik_Blob *b, Laik_Blob *xexact, TestSymmetryData &testsymmetry_data);
#else
extern int TestSymmetry(SparseMatrix &A, Vector &b, Vector &xexact, TestSymmetryData &testsymmetry_data);
#endif

#endif  // TESTSYMMETRY_HPP
