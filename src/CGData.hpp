
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
 @file CGData.hpp

 HPCG data structure
 */

#ifndef CGDATA_HPP
#define CGDATA_HPP

#ifndef USE_LAIK
#define USE_LAIK
#endif

#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "laik_instance.hpp"

struct CGData_STRUCT {
#ifdef USE_LAIK
  Laik_Blob * r;  //!< pointer to residual vector
  Laik_Blob * z;  //!< pointer to preconditioned residual vector
  Laik_Blob * p;  //!< pointer to direction vector
  Laik_Blob * Ap; //!< pointer to Krylov vector
#else
  Vector r;  //!< pointer to residual vector
  Vector z;  //!< pointer to preconditioned residual vector
  Vector p;  //!< pointer to direction vector
  Vector Ap; //!< pointer to Krylov vector
#endif
};
typedef struct CGData_STRUCT CGData;

/*!
 Constructor for the data structure of CG vectors.

 @param[in]  A    the data structure that describes the problem matrix and its structure
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 */
inline void InitializeSparseCGData(SparseMatrix & A, CGData & data) {
  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;

#ifdef USE_LAIK
  data.r = init_blob(A, false);
  data.z = init_blob(A, true);
  data.p = init_blob(A, true);
  data.Ap = init_blob(A, false);
#else
  InitializeVector(data.r, nrow);
  InitializeVector(data.z, ncol);
  InitializeVector(data.p, ncol);
  InitializeVector(data.Ap, nrow);
#endif

  return;
}

/*!
 Destructor for the CG vectors data.

 @param[inout] data the CG vectors data structure whose storage is deallocated
 */
inline void DeleteCGData(CGData & data) {

#ifdef USE_LAIK
  // TODO delete Laik vectors
#else
  DeleteVector (data.r);
  DeleteVector (data.z);
  DeleteVector (data.p);
  DeleteVector (data.Ap);
#endif
  return;
}

#endif // CGDATA_HPP

