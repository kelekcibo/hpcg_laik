
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
 @file TestSymmetry.cpp

 HPCG routine
 */

// The MPI include must be first for Windows platforms
#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif
#include <fstream>
#include <iostream>
#include <cfloat>
using std::endl;
#include <vector>
#include <cmath>

#include "hpcg.hpp"

#include "ComputeSPMV.hpp"
#include "ComputeMG.hpp"
#include "ComputeDotProduct.hpp"
#include "ComputeResidual.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "TestSymmetry.hpp"

/*!
  Tests symmetry-preserving properties of the sparse matrix vector multiply and multi-grid routines.

  @param[in]    geom   The description of the problem's geometry.
  @param[in]    A      The known system matrix
  @param[in]    b      The known right hand side vector
  @param[in]    xexact The exact solution vector
  @param[inout] testsymmetry_data The data structure with the results of the CG symmetry test including pass/fail information

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct
  @see ComputeDotProduct_ref
  @see ComputeSPMV
  @see ComputeSPMV_ref
  @see ComputeMG
  @see ComputeMG_ref
*/
int TestSymmetry(SparseMatrix & A, Vector & b, Vector & xexact, TestSymmetryData & testsymmetry_data) {

 local_int_t nrow = A.localNumberOfRows;
 local_int_t ncol = A.localNumberOfColumns;

 Vector x_ncol, y_ncol, z_ncol;
 InitializeVector(x_ncol, ncol);
 InitializeVector(y_ncol, ncol);
 InitializeVector(z_ncol, ncol);

 Laik_Blob * x_ncol_blob;
 Laik_Blob * y_ncol_blob;
 Laik_Blob * z_ncol_blob;
 x_ncol_blob = init_blob(A.totalNumberOfRows, nrow, A.A_map_data, A.A_local, A.A_ext);
 y_ncol_blob = init_blob(A.totalNumberOfRows, nrow, A.A_map_data, A.A_local, A.A_ext);
 z_ncol_blob = init_blob(A.totalNumberOfRows, nrow, A.A_map_data, A.A_local, A.A_ext);

 double t4 = 0.0; // Needed for dot-product call, otherwise unused
 testsymmetry_data.count_fail = 0;

 // Test symmetry of matrix

 // First load vectors with random values
 FillRandomVector(x_ncol);
 FillRandomVector(y_ncol);
 CopyVectorToLaikVector(x_ncol, x_ncol_blob);
 CopyVectorToLaikVector(y_ncol, y_ncol_blob);

 //  fillRandomLaikVector(x_ncol_blob);
 //  fillRandomLaikVector(y_ncol_blob);
 //  fillRandomLaikVector(z_ncol_blob);

 double xNorm2, yNorm2;
 double ANorm = 2 * 26.0;

 // Next, compute x'*A*y
 ComputeDotProduct(nrow, y_ncol, y_ncol, yNorm2, t4, A.isDotProductOptimized, y_ncol_blob, y_ncol_blob);
 int ierr = ComputeSPMV(A, y_ncol, z_ncol, y_ncol_blob); // z_nrow = A*y_overlap
 CopyVectorToLaikVector(z_ncol, z_ncol_blob);            // FIXME. use another parameter in computeMG for const Vector& r.
 if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
 double xtAy = 0.0;
 ierr = ComputeDotProduct(nrow, x_ncol, z_ncol, xtAy, t4, A.isDotProductOptimized, x_ncol_blob, z_ncol_blob); // x'*A*y
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 // Next, compute y'*A*x
 ComputeDotProduct(nrow, x_ncol, x_ncol, xNorm2, t4, A.isDotProductOptimized, x_ncol_blob, x_ncol_blob);
 ierr = ComputeSPMV(A, x_ncol, z_ncol, x_ncol_blob); // b_computed = A*x_overlap
 CopyVectorToLaikVector(z_ncol, z_ncol_blob);        // FIXME. use another parameter in computeMG for const Vector& r.
 if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
 double ytAx = 0.0;
 ierr = ComputeDotProduct(nrow, y_ncol, z_ncol, ytAx, t4, A.isDotProductOptimized, y_ncol_blob, z_ncol_blob); // y'*A*x
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 testsymmetry_data.depsym_spmv = std::fabs((long double) (xtAy - ytAx))/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON));
 if (testsymmetry_data.depsym_spmv > 1.0) ++testsymmetry_data.count_fail;  // If the difference is > 1, count it wrong
 if (A.geom->rank==0) HPCG_fout << "Departure from symmetry (scaled) for SpMV abs(x'*A*y - y'*A*x) = " << testsymmetry_data.depsym_spmv << endl;

 // Test symmetry of multi-grid

 // Compute x'*Minv*y
 ierr = ComputeMG(A, y_ncol, z_ncol, z_ncol_blob); // z_ncol = Minv*y_ncol
 CopyVectorToLaikVector(z_ncol, z_ncol_blob);      // FIXME. use another parameter in computeMG for const Vector& r.

 if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
 double xtMinvy = 0.0;
 ierr = ComputeDotProduct(nrow, x_ncol, z_ncol, xtMinvy, t4, A.isDotProductOptimized, x_ncol_blob, z_ncol_blob); // x'*Minv*y
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 // Next, compute z'*Minv*x
 CopyLaikVectorToVector(x_ncol_blob, x_ncol); // FIXME. use another parameter in computeMG for const Vector& r.
 ierr = ComputeMG(A, x_ncol, z_ncol, z_ncol_blob); // z_ncol = Minv*x_ncol
 CopyVectorToLaikVector(z_ncol, z_ncol_blob);      // FIXME. use another parameter in computeMG for const Vector& r.
 if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
 double ytMinvx = 0.0;
 ierr = ComputeDotProduct(nrow, y_ncol, z_ncol, ytMinvx, t4, A.isDotProductOptimized, y_ncol_blob, z_ncol_blob); // y'*Minv*x
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 testsymmetry_data.depsym_mg = std::fabs((long double) (xtMinvy - ytMinvx))/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON));
 if (testsymmetry_data.depsym_mg > 1.0) ++testsymmetry_data.count_fail;  // If the difference is > 1, count it wrong
 if (A.geom->rank==0) HPCG_fout << "Departure from symmetry (scaled) for MG abs(x'*Minv*y - y'*Minv*x) = " << testsymmetry_data.depsym_mg << endl;

 CopyVector(xexact, x_ncol); // Copy exact answer into overlap vector
 CopyVectorToLaikVector(x_ncol, x_ncol_blob);

 int numberOfCalls = 2;
 double residual = 0.0;
 for (int i=0; i< numberOfCalls; ++i) {
   ierr = ComputeSPMV(A, x_ncol, z_ncol, x_ncol_blob); // b_computed = A*x_overlap
   if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
   if ((ierr = ComputeResidual(A.localNumberOfRows, b, z_ncol, residual)))
     HPCG_fout << "Error in call to compute_residual: " << ierr << ".\n" << endl;
   if (A.geom->rank==0) HPCG_fout << "SpMV call [" << i << "] Residual [" << residual << "]" << endl;
 }
 DeleteVector(x_ncol);
 DeleteVector(y_ncol);
 DeleteVector(z_ncol);

 return 0;
}

