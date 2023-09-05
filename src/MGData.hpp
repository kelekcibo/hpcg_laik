
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
 @file MGData.hpp

 HPCG data structure
 */

#ifndef MGDATA_HPP
#define MGDATA_HPP

#include <cassert>

struct Laik_Blob;

#include "Vector.hpp"
#include "laik_instance.hpp"


struct MGData_STRUCT {
  int numberOfPresmootherSteps; // Call ComputeSYMGS this many times prior to coarsening
  int numberOfPostsmootherSteps; // Call ComputeSYMGS this many times after coarsening
  local_int_t * f2cOperator; //!< 1D array containing the fine operator local IDs that will be injected into coarse space.

#ifdef USE_LAIK
  Laik_Blob *rc;  // coarse grid residual vector
  Laik_Blob *xc;  // coarse grid solution vector
  Laik_Blob *Axf; // fine grid residual vector
#else
  Vector *rc;  // coarse grid residual vector
  Vector *xc;  // coarse grid solution vector
  Vector *Axf; // fine grid residual vector
#endif
 

  /*!
   This is for storing optimized data structres created in OptimizeProblem and
   used inside optimized ComputeSPMV().
   */
  void * optimizationData;
};
typedef struct MGData_STRUCT MGData;

#ifdef USE_LAIK
/*!
 Constructor for the data structure of CG vectors.

 @param[in] Ac - Fully-formed coarse matrix
 @param[in] f2cOperator -
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 */
inline void InitializeMGData(local_int_t *f2cOperator, Laik_Blob *rc, Laik_Blob *xc, Laik_Blob *Axf, MGData &data)
{
  data.numberOfPresmootherSteps = 1;
  data.numberOfPostsmootherSteps = 1;
  data.f2cOperator = f2cOperator; // Space for injection operator
  data.rc = rc;
  data.xc = xc;
  data.Axf = Axf;
  return;
}
#else
/*!
 Constructor for the data structure of CG vectors.

 @param[in] Ac - Fully-formed coarse matrix
 @param[in] f2cOperator -
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 */
inline void InitializeMGData(local_int_t * f2cOperator, Vector * rc, Vector * xc, Vector * Axf, MGData & data) {
  data.numberOfPresmootherSteps = 1;
  data.numberOfPostsmootherSteps = 1;
  data.f2cOperator = f2cOperator; // Space for injection operator
  data.rc = rc;
  data.xc = xc;
  data.Axf = Axf;
  return;
}
#endif

/*!
 Destructor for the CG vectors data.

 @param[inout] data the MG data structure whose storage is deallocated
 */
inline void DeleteMGData(MGData & data) {

  delete [] data.f2cOperator;

#ifdef USE_LAIK
  // delete laik vectors
#else
  DeleteVector(*data.Axf);
  DeleteVector(*data.rc);
  DeleteVector(*data.xc);
  delete data.Axf;
  delete data.rc;
  delete data.xc;
#endif
  
  return;
}

#endif // MGDATA_HPP

