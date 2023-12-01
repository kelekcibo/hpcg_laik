
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

#include "laik/hpcg_laik.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

struct CGData_STRUCT {

#ifndef HPCG_NO_LAIK
  Laik_Blob * r_blob;  //!< pointer to residual vector
  Laik_Blob * z_blob;  //!< pointer to preconditioned residual vector
  Laik_Blob * p_blob;  //!< pointer to direction vector
  Laik_Blob * Ap_blob; //!< pointer to Krylov vector
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

#ifndef HPCG_NO_LAIK
  data.r_blob = init_blob(A, false);
  data.r_blob->name = "CG_Data r";
  data.z_blob = init_blob(A, true);
  data.z_blob->name = "CG_Data z";
  data.p_blob = init_blob(A, true);
  data.p_blob->name = "CG_Data p";
  data.Ap_blob = init_blob(A, false);
  data.Ap_blob->name = "CG_Data Ap";
#else
  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;

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

#ifndef HPCG_NO_LAIK
  DeleteLaikVector(data.r_blob);
  DeleteLaikVector(data.Ap_blob);
  DeleteLaikVector(data.z_blob);
  DeleteLaikVector(data.p_blob);

#else
  DeleteVector (data.r);
  DeleteVector (data.z);
  DeleteVector (data.p);
  DeleteVector (data.Ap);
#endif
  return;
}

#endif // CGDATA_HPP

