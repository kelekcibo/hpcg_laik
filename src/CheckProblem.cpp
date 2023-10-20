
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
 @file CheckProblem.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include "laik/hpcg_laik.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#if defined(HPCG_DEBUG) || defined(HPCG_DETAILED_DEBUG)
#include <fstream>
using std::endl;
#include "hpcg.hpp"
#endif
#include <cassert>
#include <iostream>

#include "CheckProblem.hpp"


/*!
  Check the contents of the generated sparse matrix to see if values match expected contents.

  @param[in]  A      The known system matrix
  @param[inout] b      The newly allocated and generated right hand side vector (if b!=0 on entry)
  @param[inout] x      The newly allocated solution vector with entries set to 0.0 (if x!=0 on entry)
  @param[inout] xexact The newly allocated solution vector with entries set to the exact solution (if the xexact!=0 non-zero on entry)

  @see GenerateGeometry
*/

void CheckProblem(SparseMatrix & A, Vector * b, Vector * x, Vector * xexact) {

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
  CheckProblem_repartition(A, b, x, xexact);
  return;
#endif
#endif

  // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  // below may result in global range values.
  global_int_t nx = A.geom->nx;
  global_int_t ny = A.geom->ny;
  global_int_t nz = A.geom->nz;
  global_int_t gnx = A.geom->gnx;
  global_int_t gny = A.geom->gny;
  global_int_t gnz = A.geom->gnz;
  global_int_t gix0 = A.geom->gix0;
  global_int_t giy0 = A.geom->giy0;
  global_int_t giz0 = A.geom->giz0;

  local_int_t localNumberOfRows = nx*ny*nz; // This is the size of our subblock
  global_int_t totalNumberOfRows = gnx*gny*gnz; // Total number of grid points in mesh

  double * bv = 0;
  double * xv = 0;
  double * xexactv = 0;
  if (b!=0) bv = b->values; // Only compute exact solution if requested
  if (x!=0) xv = x->values; // Only compute exact solution if requested
  if (xexact!=0) xexactv = xexact->values; // Only compute exact solution if requested

  local_int_t localNumberOfNonzeros = 0;
  // TODO:  This triply nested loop could be flattened or use nested parallelism
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t iz=0; iz<nz; iz++) {
    global_int_t giz = giz0+iz;
    for (local_int_t iy=0; iy<ny; iy++) {
      global_int_t giy = giy0+iy;
      for (local_int_t ix=0; ix<nx; ix++) {
        global_int_t gix = gix0+ix;
        local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
        global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
        assert(A.localToGlobalMap[currentLocalRow] == currentGlobalRow);
#ifdef HPCG_DETAILED_DEBUG
        HPCG_fout << " rank, globalRow, localRow = " << A.geom->rank << " " << currentGlobalRow << " " << A.globalToLocalMap.find(currentGlobalRow)->second << endl;
#endif
        char numberOfNonzerosInRow = 0;
        double * currentValuePointer = A.matrixValues[currentLocalRow]; // Pointer to current value in current row
        global_int_t * currentIndexPointerG = A.mtxIndG[currentLocalRow]; // Pointer to current index in current row
        for (int sz=-1; sz<=1; sz++) {
          if (giz+sz>-1 && giz+sz<gnz) {
            for (int sy=-1; sy<=1; sy++) {
              if (giy+sy>-1 && giy+sy<gny) {
                for (int sx=-1; sx<=1; sx++) {
                  if (gix+sx>-1 && gix+sx<gnx) {
                    global_int_t curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx;
                    if (curcol==currentGlobalRow) {
                      assert(A.matrixDiagonal[currentLocalRow] == currentValuePointer);
                      assert(*currentValuePointer++ == 26.0);
                    } else {
                      assert(*currentValuePointer++ == -1.0);
                    }
                    assert(*currentIndexPointerG++ == curcol);
                    numberOfNonzerosInRow++;
                  } // end x bounds test
                } // end sx loop
              } // end y bounds test
            } // end sy loop
          } // end z bounds test
        } // end sz loop
        assert(A.nonzerosInRow[currentLocalRow] == numberOfNonzerosInRow);
#ifndef HPCG_NO_OPENMP
        #pragma omp critical
#endif
        localNumberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
        if (b!=0)      assert(bv[currentLocalRow] == 26.0 - ((double) (numberOfNonzerosInRow-1)));
        if (x!=0)      assert(xv[currentLocalRow] == 0.0);
        if (xexact!=0) assert(xexactv[currentLocalRow] == 1.0);
      } // end ix loop
    } // end iy loop
  } // end iz loop
#ifdef HPCG_DETAILED_DEBUG
  HPCG_fout     << "Process " << A.geom->rank << " of " << A.geom->size <<" has " << localNumberOfRows    << " rows."     << endl
      << "Process " << A.geom->rank << " of " << A.geom->size <<" has " << localNumberOfNonzeros<< " nonzeros." <<endl;
#endif

  global_int_t totalNumberOfNonzeros = 0;
#ifndef HPCG_NO_MPI
  // Use reduce function to sum all nonzeros
#ifdef HPCG_NO_LONG_LONG
  laik_allreduce((void *)&localNumberOfNonzeros, (void *)&totalNumberOfNonzeros, 1, laik_Int64, LAIK_RO_Sum);
#else
  long long lnnz = localNumberOfNonzeros, gnnz = 0; // convert to 64 bit for MPI call
  laik_allreduce((void *)&lnnz, (void *)&gnnz, 1, laik_UInt64, LAIK_RO_Sum); // Int64 == long_long_int ?
  totalNumberOfNonzeros = gnnz; // Copy back
#endif
#else
  totalNumberOfNonzeros = localNumberOfNonzeros;
#endif

  assert(A.totalNumberOfRows == totalNumberOfRows);
  assert(A.totalNumberOfNonzeros == totalNumberOfNonzeros);
  assert(A.localNumberOfRows == localNumberOfRows);
  assert(A.localNumberOfNonzeros == localNumberOfNonzeros);

  return;
}

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
void CheckProblem_repartition(SparseMatrix &A, Vector *b, Vector *x, Vector *xexact)
{
  // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  // below may result in global range values.
  global_int_t nx = A.geom->nx;
  global_int_t ny = A.geom->ny;
  global_int_t nz = A.geom->nz;
  global_int_t gnx = A.geom->gnx;
  global_int_t gny = A.geom->gny;
  global_int_t gnz = A.geom->gnz;
  global_int_t gix0 = A.geom->gix0;
  global_int_t giy0 = A.geom->giy0;
  global_int_t giz0 = A.geom->giz0;

  local_int_t localNumberOfRows = nx * ny * nz;     // This is the size of our subblock
  global_int_t totalNumberOfRows = gnx * gny * gnz; // Total number of grid points in mesh

  char *nonzerosInRow;
  laik_get_map_1d(A.nonzerosInRow_d, 0, (void **)&nonzerosInRow, 0);
  double *matrixDiagonal;
  laik_get_map_1d(A.matrixDiagonal_d, 0, (void **)&matrixDiagonal, 0);
  uint64_t *mtxIndG;
  laik_get_map_1d(A.mtxIndG_d, 0, (void **)&mtxIndG, 0);
  double *matrixValues;
  laik_get_map_1d(A.matrixValues_d, 0, (void **)&matrixValues, 0);

  double *bv = 0;
  double *xv = 0;
  double *xexactv = 0;
  if (b != 0)
    bv = b->values; // Only compute exact solution if requested
  if (x != 0)
    xv = x->values; // Only compute exact solution if requested
  if (xexact != 0)
    xexactv = xexact->values; // Only compute exact solution if requested

    local_int_t localNumberOfNonzeros = 0;
    // TODO:  This triply nested loop could be flattened or use nested parallelism
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for (local_int_t iz = 0; iz < nz; iz++)
  {
    global_int_t giz = giz0 + iz;
    for (local_int_t iy = 0; iy < ny; iy++)
    {
      global_int_t giy = giy0 + iy;
      for (local_int_t ix = 0; ix < nx; ix++)
      {
        global_int_t gix = gix0 + ix;
        local_int_t currentLocalRow = iz * nx * ny + iy * nx + ix;
        global_int_t currentGlobalRow = giz * gnx * gny + giy * gnx + gix;
        assert(A.localToGlobalMap[currentLocalRow] == currentGlobalRow);
#ifdef HPCG_DETAILED_DEBUG
        HPCG_fout << " rank, globalRow, localRow = " << A.geom->rank << " " << currentGlobalRow << " " << A.globalToLocalMap.find(currentGlobalRow)->second << endl;
#endif
        char numberOfNonzerosInRow = 0;
        uint64_t currentValuePointer_index = 0; // Index to current value in current row
        global_int_t currentIndexPointerG_index = 0;  // Index to current index in current row
        for (int sz = -1; sz <= 1; sz++)
        {
          if (giz + sz > -1 && giz + sz < gnz)
          {
            for (int sy = -1; sy <= 1; sy++)
            {
              if (giy + sy > -1 && giy + sy < gny)
              {
                for (int sx = -1; sx <= 1; sx++)
                {
                  if (gix + sx > -1 && gix + sx < gnx)
                  {
                    global_int_t curcol = currentGlobalRow + sz * gnx * gny + sy * gnx + sx;
                    if (curcol == currentGlobalRow)
                    {
                      assert(matrixDiagonal[map_l2a_A(A, currentLocalRow)] == matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + currentValuePointer_index]);
                      assert(matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + ++currentValuePointer_index] == 26.0);
                    }
                    else
                    {
                      assert(matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + ++currentValuePointer_index] == -1.0);
                    }
                    assert(mtxIndG[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + ++currentIndexPointerG_index] == curcol);
                    numberOfNonzerosInRow++;
                  } // end x bounds test
                }   // end sx loop
              }     // end y bounds test
            }       // end sy loop
          }         // end z bounds test
        }           // end sz loop
        assert(nonzerosInRow[map_l2a_A(A, currentLocalRow)] == numberOfNonzerosInRow);
#ifndef HPCG_NO_OPENMP
#pragma omp critical
#endif
        localNumberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
        if (b != 0)
          assert(bv[currentLocalRow] == 26.0 - ((double)(numberOfNonzerosInRow - 1)));
        if (x != 0)
          assert(xv[currentLocalRow] == 0.0);
        if (xexact != 0)
          assert(xexactv[currentLocalRow] == 1.0);
      } // end ix loop
    }   // end iy loop
  }     // end iz loop
#ifdef HPCG_DETAILED_DEBUG
  HPCG_fout << "Process " << A.geom->rank << " of " << A.geom->size << " has " << localNumberOfRows << " rows." << endl
            << "Process " << A.geom->rank << " of " << A.geom->size << " has " << localNumberOfNonzeros << " nonzeros." << endl;
#endif

  global_int_t totalNumberOfNonzeros = 0;
#ifndef HPCG_NO_MPI
  // Use reduce function to sum all nonzeros
#ifdef HPCG_NO_LONG_LONG
  laik_allreduce((void *)&localNumberOfNonzeros, (void *)&totalNumberOfNonzeros, 1, laik_Int64, LAIK_RO_Sum);
#else
  long long lnnz = localNumberOfNonzeros, gnnz = 0;                          // convert to 64 bit for MPI call
  laik_allreduce((void *)&lnnz, (void *)&gnnz, 1, laik_UInt64, LAIK_RO_Sum); // Int64 == long_long_int ?
  totalNumberOfNonzeros = gnnz;                                              // Copy back
#endif
#else
  totalNumberOfNonzeros = localNumberOfNonzeros;
#endif

  assert(A.totalNumberOfRows == totalNumberOfRows);
  assert(A.totalNumberOfNonzeros == totalNumberOfNonzeros);
  assert(A.localNumberOfRows == localNumberOfRows);
  assert(A.localNumberOfNonzeros == localNumberOfNonzeros);

  return;
}
#endif // HPCG_NO_LAIK
#endif // REPARTITION