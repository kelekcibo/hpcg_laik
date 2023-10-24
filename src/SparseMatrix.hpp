
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
 @file SparseMatrix.hpp

 HPCG data structures for the sparse matrix
 */

#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <vector>
#include <cassert>
#include <map>

#ifndef HPCG_NO_LAIK
#include "laik/hpcg_laik.hpp"

// forw. decl.
struct Local2Allocation_map_x;
typedef Local2Allocation_map_x L2A_map;
void free_L2A_map(L2A_map *mapping);
#ifdef REPARTITION
allocation_int_t map_l2a_A(const SparseMatrix &A, local_int_t localIndex);
void replaceMatrixValues(SparseMatrix &A);
#endif
// forw. decl.
#endif

#include "Geometry.hpp"
#include "Vector.hpp"
#include "MGData.hpp"
#if __cplusplus < 201103L
// for C++03
#include <map>
    typedef std::map<global_int_t, local_int_t> GlobalToLocalMap;
#else
// for C++11 or greater
#include <unordered_map>
using GlobalToLocalMap = std::unordered_map< global_int_t, local_int_t >;
#endif

struct SparseMatrix_STRUCT {
  char  * title; //!< name of the sparse matrix
  Geometry * geom; //!< geometry associated with this matrix
  global_int_t totalNumberOfRows; //!< total number of matrix rows across all processes
  global_int_t totalNumberOfNonzeros; //!< total number of matrix nonzeros across all processes
  local_int_t localNumberOfRows; //!< number of rows local to this process
  local_int_t localNumberOfColumns;  //!< number of columns local to this process
  local_int_t localNumberOfNonzeros;  //!< number of nonzeros local to this process
  char  * nonzerosInRow;  //!< The number of nonzeros in a row will always be 27 or fewer
  global_int_t ** mtxIndG; //!< matrix indices as global values
  local_int_t ** mtxIndL; //!< matrix indices as local values
  double ** matrixValues; //!< values of matrix entries
  double ** matrixDiagonal; //!< values of matrix diagonal entries
  GlobalToLocalMap globalToLocalMap; //!< global-to-local mapping
  std::vector< global_int_t > localToGlobalMap; //!< local-to-global mapping
  mutable bool isDotProductOptimized;
  mutable bool isSpmvOptimized;
  mutable bool isMgOptimized;
  mutable bool isWaxpbyOptimized;
  /*!
   This is for storing optimized data structres created in OptimizeProblem and
   used inside optimized ComputeSPMV().
   */
  mutable struct SparseMatrix_STRUCT * Ac; // Coarse grid matrix
  mutable MGData * mgData; // Pointer to the coarse level data for this fine matrix
  void * optimizationData;  // pointer that can be used to store implementation-specific data

#ifndef HPCG_NO_MPI
  local_int_t numberOfExternalValues; //!< number of entries that are external to this process
  int numberOfSendNeighbors; //!< number of neighboring processes that will be send local data
  local_int_t totalToBeSent; //!< total number of entries to be sent
  local_int_t * elementsToSend; //!< elements to send to neighboring processes
  int * neighbors; //!< neighboring processes
  local_int_t * receiveLength; //!< lenghts of messages received from neighboring processes
  local_int_t * sendLength; //!< lenghts of messages sent to neighboring processes
  double * sendBuffer; //!< send buffer for non-blocking sends

#ifndef HPCG_NO_LAIK
  // ############### Data needed to create partitionings and Laik_Data container
  std::map<local_int_t, global_int_t> localToExternalMap; /* Needed for LAIK (@see L2A_map)*/
  L2A_map * mapping;
  Laik_Space * space;
  Laik_Partitioning * ext;
  Laik_Partitioning * local;

#ifdef REPARTITION

    bool repartition_me; /* Tell the app, that a reseize should happen. We want to test it during the call to CG_REFin CG Reference Timing Phase */

    uint64_t * mapping_; // @see L2A_map. Same applies for allocation buffers of A
    int offset_;
    
    // Special space for 2D arrays implemented as 1D array
    Laik_Space *space2d; 

    // Partitionings for ressources below
    Laik_Partitioning * partitioning_1d;
    Laik_Partitioning * partitioning_2d;

    // Ressources of this matrix, which will be partitioned
    Laik_Data * nonzerosInRow_d;    //!< The number of nonzeros in a row will always be 27 or fewer
    Laik_Data * mtxIndG_d;          //!< matrix indices as global values
    Laik_Data * matrixValues_d;     //!< values of matrix entries
    Laik_Data * matrixDiagonal_d;   //!< values of matrix diagonal entries

    /*
      This variable is only for x_l
      He will store the pointer to xexact_l
      We need to re-switch this vector as well
      But when reseizing, xexact_l is out of scope
      Quick solution is this here
    */
    Laik_Blob * ptr_to_xexact = 0;

#endif // REPARTITION
#endif // HPCG_NO_LAIK
#endif // HPCG_NO_MPI
};
typedef struct SparseMatrix_STRUCT SparseMatrix;

/*!
  Initializes the known system matrix data structure members to 0.

  @param[in] A the known system matrix
 */
inline void InitializeSparseMatrix(SparseMatrix & A, Geometry * geom) {
  A.title = 0;
  A.geom = geom;
  A.totalNumberOfRows = 0;
  A.totalNumberOfNonzeros = 0;
  A.localNumberOfRows = 0;
  A.localNumberOfColumns = 0;
  A.localNumberOfNonzeros = 0;
  A.nonzerosInRow = 0;
  A.mtxIndG = 0;
  A.mtxIndL = 0;
  A.matrixValues = 0;
  A.matrixDiagonal = 0;

  // Optimization is ON by default. The code that switches it OFF is in the
  // functions that are meant to be optimized.
  A.isDotProductOptimized = true;
  A.isSpmvOptimized       = true;
  A.isMgOptimized      = true;
  A.isWaxpbyOptimized     = true;

#ifndef HPCG_NO_MPI
  A.numberOfExternalValues = 0;
  A.numberOfSendNeighbors = 0;
  A.totalToBeSent = 0;
  A.elementsToSend = 0;
  A.neighbors = 0;
  A.receiveLength = 0;
  A.sendLength = 0;
  A.sendBuffer = 0;
  
  #ifndef HPCG_NO_LAIK
  A.mapping = 0;
  A.space = 0;
  A.ext = 0;
  A.local = 0;
  
    #ifdef REPARTITION
      A.repartition_me = false;
      A.mapping_ = 0;
      A.offset_ = 0;
      A.space2d = 0;
      A.partitioning_1d = 0;
      A.partitioning_2d = 0;
      A.nonzerosInRow_d = 0;
      A.mtxIndG_d = 0;
      A.matrixValues_d = 0;
      A.matrixDiagonal_d = 0;
      A.ptr_to_xexact = 0;
    #endif // REPARTITION
  #endif // HPCG_NO_LAIK
#endif // HPCG_NO_MPI
  A.mgData = 0; // Fine-to-coarse grid transfer initially not defined.
  A.Ac =0;
  return;
}

/*!
  Copy values from matrix diagonal into user-provided vector.

  @param[in] A the known system matrix.
  @param[inout] diagonal  Vector of diagonal values (must be allocated before call to this function).
 */
inline void CopyMatrixDiagonal(SparseMatrix & A, Vector & diagonal) {

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
    double *matrixDiagonal;
    laik_get_map_1d(A.matrixDiagonal_d, 0, (void **)&matrixDiagonal, 0);
    double *dia_v = diagonal.values;
    assert(A.localNumberOfRows == diagonal.localLength);
    for (local_int_t i=0; i<A.localNumberOfRows; ++i) dia_v[i] = matrixDiagonal[map_l2a_A(A, i)];
  return;
#endif
#endif

    double ** curDiagA = A.matrixDiagonal;
    double * dv = diagonal.values;
    assert(A.localNumberOfRows==diagonal.localLength);
    for (local_int_t i=0; i<A.localNumberOfRows; ++i) dv[i] = *(curDiagA[i]);
  return;
}
/*!
  Replace specified matrix diagonal value.

  @param[inout] A The system matrix.
  @param[in] diagonal  Vector of diagonal values that will replace existing matrix diagonal values.
 */
inline void ReplaceMatrixDiagonal(SparseMatrix & A, Vector & diagonal) {

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
    double *matrixDiagonal;
    laik_get_map_1d(A.matrixDiagonal_d, 0, (void **)&matrixDiagonal, 0);
    double *dia_v = diagonal.values;
    assert(A.localNumberOfRows == diagonal.localLength);
    for (local_int_t i=0; i<A.localNumberOfRows; ++i) matrixDiagonal[map_l2a_A(A, i)] = dia_v[i];

    replaceMatrixValues(A);
    
    return;
#endif
#endif

    double ** curDiagA = A.matrixDiagonal;
    double * dv = diagonal.values;
    assert(A.localNumberOfRows==diagonal.localLength);
    for (local_int_t i=0; i<A.localNumberOfRows; ++i) *(curDiagA[i]) = dv[i];
  return;
}

/*!
  Deallocates the members of the data structure of the known system matrix provided they are not 0.

  @param[in] A the known system matrix
 */
inline void DeleteMatrix(SparseMatrix & A) {

#ifndef HPCG_CONTIGUOUS_ARRAYS
  for (local_int_t i = 0; i< A.localNumberOfRows; ++i) {
    delete [] A.matrixValues[i];
    delete [] A.mtxIndG[i];
    delete [] A.mtxIndL[i];
  }
#else
  delete [] A.matrixValues[0];
  delete [] A.mtxIndG[0];
  delete [] A.mtxIndL[0];
#endif
  if (A.title)                  delete [] A.title;
  if (A.nonzerosInRow)             delete [] A.nonzerosInRow;
  if (A.mtxIndG) delete [] A.mtxIndG;
  if (A.mtxIndL) delete [] A.mtxIndL;
  if (A.matrixValues) delete [] A.matrixValues;
  if (A.matrixDiagonal)           delete [] A.matrixDiagonal;

#ifndef HPCG_NO_MPI
  if (A.elementsToSend)       delete [] A.elementsToSend;
  if (A.neighbors)              delete [] A.neighbors;
  if (A.receiveLength)            delete [] A.receiveLength;
  if (A.sendLength)            delete [] A.sendLength;
  if (A.sendBuffer)            delete [] A.sendBuffer;

  #ifndef HPCG_NO_LAIK
    // Delete LAIK specific data
    A.globalToLocalMap.clear();
    if (A.mapping) free_L2A_map(A.mapping);
    if (A.space) { laik_free_space(A.space); A.space = 0; }
    if (A.local) { laik_free_partitioning(A.local); A.local = 0; };
    if (A.ext) { laik_free_partitioning(A.ext); A.ext = 0; };

    #ifdef REPARTITION
      if (A.mapping_) {  delete [] A.mapping_; A.mapping_ = 0;};
      if (A.space2d) { laik_free_space(A.space2d); A.space2d = 0; }
      if (A.partitioning_1d) { laik_free_partitioning(A.partitioning_1d); A.partitioning_1d = 0; };
      if (A.partitioning_2d) { laik_free_partitioning(A.partitioning_2d); A.partitioning_2d = 0; };
      if (A.mtxIndG_d) { laik_free(A.mtxIndG_d); };
      if (A.matrixValues_d) { laik_free(A.matrixValues_d); };
      if (A.nonzerosInRow_d) { laik_free(A.nonzerosInRow_d); };
      if (A.matrixDiagonal_d) { laik_free(A.matrixDiagonal_d); };
    #endif // REPARTITION
  #endif // HPCG_NO_LAIK
#endif // HPCG_NO_MPI

  if (A.geom!=0) { DeleteGeometry(*A.geom); delete A.geom; A.geom = 0;}
  if (A.Ac!=0) { DeleteMatrix(*A.Ac); delete A.Ac; A.Ac = 0;} // Delete coarse matrix
  if (A.mgData!=0) { DeleteMGData(*A.mgData); delete A.mgData; A.mgData = 0;} // Delete MG data

  return;
}

  

#endif // SPARSEMATRIX_HPP
