
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
 @file SetupHalo_ref.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include <map>
#include <set>
#include <cstring>
#include <iostream>
#include <cstdlib>

#include "laik/hpcg_laik.hpp"

#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

// #define HPCG_DETAILED_DEBUG
#ifdef HPCG_DETAILED_DEBUG
#include <fstream>
using std::endl;
#include "hpcg.hpp"
#include <cassert>
#endif

#include "SetupHalo_ref.hpp"
#include "mytimer.hpp"

/*!
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/
void SetupHalo_ref(SparseMatrix & A) {

  // Extract Matrix pieces
#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
  SetupHalo_repartition_ref(A);
  return ;
#endif
#endif


  local_int_t localNumberOfRows = A.localNumberOfRows;
  char  * nonzerosInRow = A.nonzerosInRow;
  global_int_t ** mtxIndG = A.mtxIndG;
  local_int_t ** mtxIndL = A.mtxIndL;

#ifdef HPCG_NO_MPI  // In the non-MPI case we simply copy global indices to local index storage
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< localNumberOfRows; i++) {
    int cur_nnz = nonzerosInRow[i];
    for (int j=0; j<cur_nnz; j++) mtxIndL[i][j] = mtxIndG[i][j];
  }

#else // Run this section if compiling for MPI

  // Scan global IDs of the nonzeros in the matrix.  Determine if the column ID matches a row ID.  If not:
  // 1) We call the ComputeRankOfMatrixRow function, which tells us the rank of the processor owning the row ID.
  //  We need to receive this value of the x vector during the halo exchange.
  // 2) We record our row ID since we know that the other processor will need this value from us, due to symmetry.

  std::map< int, std::set< global_int_t> > sendList, receiveList;
  typedef std::map< int, std::set< global_int_t> >::iterator map_iter;
  typedef std::set<global_int_t>::iterator set_iter;
  std::map< global_int_t, local_int_t > externalToLocalMap;

  // TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time
  for (local_int_t i=0; i< localNumberOfRows; i++) {
    global_int_t currentGlobalRow = A.localToGlobalMap[i];
    for (int j=0; j<nonzerosInRow[i]; j++) {
      global_int_t curIndex = mtxIndG[i][j];
      int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
#ifdef HPCG_DETAILED_DEBUG
      // HPCG_fout << "rank, row , col, globalToLocalMap[col] = " << A.geom->rank << " " << currentGlobalRow << " "
      //     << curIndex << " " << A.globalToLocalMap[curIndex] << endl;
#endif  
      if (A.geom->rank!=rankIdOfColumnEntry) {// If column index is not a row index, then it comes from another processor
        receiveList[rankIdOfColumnEntry].insert(curIndex);
        sendList[rankIdOfColumnEntry].insert(currentGlobalRow); // Matrix symmetry means we know the neighbor process wants my value
      }
    }
  }

  // Count number of matrix entries to send and receive
  local_int_t totalToBeSent = 0;
  for (map_iter curNeighbor = sendList.begin(); curNeighbor != sendList.end(); ++curNeighbor) {
    totalToBeSent += (curNeighbor->second).size();
  }
  local_int_t totalToBeReceived = 0;
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor) {
    totalToBeReceived += (curNeighbor->second).size();
  }
#ifdef HPCG_DETAILED_DEBUG
  // These are all attributes that should be true, due to symmetry
  HPCG_fout << "totalToBeSent = " << totalToBeSent << " totalToBeReceived = " << totalToBeReceived << endl;
  assert(totalToBeSent==totalToBeReceived); // Number of sent entry should equal number of received
  assert(sendList.size()==receiveList.size()); // Number of send-to neighbors should equal number of receive-from
  // Each receive-from neighbor should be a send-to neighbor, and send the same number of entries
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor) {
    assert(sendList.find(curNeighbor->first)!=sendList.end());
    assert(sendList[curNeighbor->first].size()==receiveList[curNeighbor->first].size());
  }
#endif

  // Build the arrays and lists needed by the ExchangeHalo function.
  double * sendBuffer = new double[totalToBeSent];
  local_int_t * elementsToSend = new local_int_t[totalToBeSent];
  int * neighbors = new int[sendList.size()];
  local_int_t * receiveLength = new local_int_t[receiveList.size()];
  local_int_t * sendLength = new local_int_t[sendList.size()];
  int neighborCount = 0;
  local_int_t receiveEntryCount = 0;
  local_int_t sendEntryCount = 0;
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor, ++neighborCount) {
    int neighborId = curNeighbor->first; // rank of current neighbor we are processing
    neighbors[neighborCount] = neighborId; // store rank ID of current neighbor
    receiveLength[neighborCount] = receiveList[neighborId].size();
    sendLength[neighborCount] = sendList[neighborId].size(); // Get count if sends/receives
    for (set_iter i = receiveList[neighborId].begin(); i != receiveList[neighborId].end(); ++i, ++receiveEntryCount) {
      externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount; // The remote columns are indexed at end of internals
      A.localToExternalMap[localNumberOfRows + receiveEntryCount] = *i;
    }
    for (set_iter i = sendList[neighborId].begin(); i != sendList[neighborId].end(); ++i, ++sendEntryCount) {
      //if (geom.rank==1) HPCG_fout << "*i, globalToLocalMap[*i], sendEntryCount = " << *i << " " << A.globalToLocalMap[*i] << " " << sendEntryCount << endl;
      elementsToSend[sendEntryCount] = A.globalToLocalMap[*i]; // store local ids of entry to send
    }
  }


  // Convert matrix indices to local IDs
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< localNumberOfRows; i++) {
    for (int j=0; j<nonzerosInRow[i]; j++) {
      global_int_t curIndex = mtxIndG[i][j];
      int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
      if (A.geom->rank==rankIdOfColumnEntry) { // My column index, so convert to local index
        mtxIndL[i][j] = A.globalToLocalMap[curIndex];
      } else { // If column index is not a row index, then it comes from another processor
        mtxIndL[i][j] = externalToLocalMap[curIndex];
      }
    }
  }

  // Store contents in our matrix struct
  A.numberOfExternalValues = externalToLocalMap.size();
  A.localNumberOfColumns = A.localNumberOfRows + A.numberOfExternalValues;
  A.numberOfSendNeighbors = sendList.size();
  A.totalToBeSent = totalToBeSent;
  A.elementsToSend = elementsToSend;
  A.neighbors = neighbors;
  A.receiveLength = receiveLength;
  A.sendLength = sendLength;
  A.sendBuffer = sendBuffer;

#ifndef HPCG_NO_LAIK
  // ########## Data for partitioning algorithm
  partition_d *pt_data_local = (partition_d *)malloc(sizeof(partition_d));
  partition_d *pt_data_ext = (partition_d *)malloc(sizeof(partition_d));

  std::memcpy((void *)&pt_data_ext->receiveList, (void *)&receiveList, sizeof(receiveList));

  init_partition_data(A, pt_data_local, pt_data_ext);
  // ########## Data for partitioning algorithm

  // ########## Data to calculate mapping
  L2A_map * map_data = (L2A_map *) malloc(sizeof(L2A_map));

  map_data->localNumberOfRows = A.localNumberOfRows;
  map_data->offset = -1;
  map_data->offset_ext = -1;

  std::memcpy((void *)&map_data->localToGlobalMap, (void *)&A.localToGlobalMap, sizeof(A.localToGlobalMap));
  std::memcpy((void *)&map_data->localToExternalMap, (void *)&A.localToExternalMap, sizeof(A.localToExternalMap));

  A.mapping = map_data;
  // ########## Data to calculate mapping

  A.space = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows);

  init_partitionings(A, pt_data_local, pt_data_ext);
#endif // HPCG_NO_LAIK

#ifdef HPCG_DETAILED_DEBUG
  HPCG_fout << " For rank " << A.geom->rank << " of " << A.geom->size << ", number of neighbors = " << A.numberOfSendNeighbors << endl;
  for (int i = 0; i < A.numberOfSendNeighbors; i++) {
    HPCG_fout << "     rank " << A.geom->rank << " neighbor " << neighbors[i] << " send/recv length = " << sendLength[i] << "/" << receiveLength[i] << endl;
    for (local_int_t j = 0; j<sendLength[i]; ++j)
      HPCG_fout << "       rank " << A.geom->rank << " elementsToSend[" << j << "] = " << elementsToSend[j] << endl;
  }
#endif

#endif
// ifdef HPCG_NO_MPI

  return;
}

#ifndef HPCG_NO_LAIK
#ifdef REPARTITION
/*!
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/
void SetupHalo_repartition_ref(SparseMatrix &A)
{

  // Extract Matrix pieces
  local_int_t localNumberOfRows = A.localNumberOfRows;
  local_int_t **mtxIndL = A.mtxIndL;

    
  char * nonzerosInRow;
  laik_get_map_1d(A.nonzerosInRow_d, 0, (void **) &nonzerosInRow, 0);

  global_int_t * mtxIndG;
  laik_get_map_1d(A.mtxIndG_d, 0, (void **)&mtxIndG, 0);

  // Scan global IDs of the nonzeros in the matrix.  Determine if the column ID matches a row ID.  If not:
  // 1) We call the ComputeRankOfMatrixRow function, which tells us the rank of the processor owning the row ID.
  //  We need to receive this value of the x vector during the halo exchange.
  // 2) We record our row ID since we know that the other processor will need this value from us, due to symmetry.

  std::map<int, std::set<global_int_t>> sendList, receiveList;
  typedef std::map<int, std::set<global_int_t>>::iterator map_iter;
  typedef std::set<global_int_t>::iterator set_iter;
  std::map<global_int_t, local_int_t> externalToLocalMap;

  // TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time

  for (local_int_t i = 0; i < localNumberOfRows; i++)
  {
    global_int_t currentGlobalRow = A.localToGlobalMap[i];
    for (int j = 0; j < nonzerosInRow[map_l2a_A(A, i)]; j++)
    {
      global_int_t curIndex = mtxIndG[map_l2a_A(A, i) * numberOfNonzerosPerRow + j];
      int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
#ifdef HPCG_DETAILED_DEBUG
      // HPCG_fout << "rank, row , col, globalToLocalMap[col] = " << A.geom->rank << " " << currentGlobalRow << " "
      //     << curIndex << " " << A.globalToLocalMap[curIndex] << endl;
#endif
      if (A.geom->rank != rankIdOfColumnEntry)
      { // If column index is not a row index, then it comes from another processor
        receiveList[rankIdOfColumnEntry].insert(curIndex);
        sendList[rankIdOfColumnEntry].insert(currentGlobalRow); // Matrix symmetry means we know the neighbor process wants my value
      }
    }
  }


  // HPCG_fout << structure;

  // Count number of matrix entries to send and receive
  local_int_t totalToBeSent = 0;
  for (map_iter curNeighbor = sendList.begin(); curNeighbor != sendList.end(); ++curNeighbor)
  {
    totalToBeSent += (curNeighbor->second).size();
  }
  local_int_t totalToBeReceived = 0;
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor)
  {
    totalToBeReceived += (curNeighbor->second).size();
  }
#ifdef HPCG_DETAILED_DEBUG
  // These are all attributes that should be true, due to symmetry
  HPCG_fout << "totalToBeSent = " << totalToBeSent << " totalToBeReceived = " << totalToBeReceived << endl;
  assert(totalToBeSent == totalToBeReceived);    // Number of sent entry should equal number of received
  assert(sendList.size() == receiveList.size()); // Number of send-to neighbors should equal number of receive-from
  // Each receive-from neighbor should be a send-to neighbor, and send the same number of entries
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor)
  {
    assert(sendList.find(curNeighbor->first) != sendList.end());
    assert(sendList[curNeighbor->first].size() == receiveList[curNeighbor->first].size());
  }
#endif
  if (laik_phase(hpcg_instance) > 0 && laik_myid(world) < 2)
  {
    printf("OLD PROCS ARE EXITING!\n");
    exit_hpcg_run("re_switch_LaikVectors", false);
  }

  // Build the arrays and lists needed by the ExchangeHalo function.
  double *sendBuffer = new double[totalToBeSent];
  local_int_t *elementsToSend = new local_int_t[totalToBeSent];
  int *neighbors = new int[sendList.size()];
  local_int_t *receiveLength = new local_int_t[receiveList.size()];
  local_int_t *sendLength = new local_int_t[sendList.size()];
  int neighborCount = 0;
  local_int_t receiveEntryCount = 0;
  local_int_t sendEntryCount = 0;
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor, ++neighborCount)
  {
    int neighborId = curNeighbor->first;   // rank of current neighbor we are processing
    neighbors[neighborCount] = neighborId; // store rank ID of current neighbor
    receiveLength[neighborCount] = receiveList[neighborId].size();
    sendLength[neighborCount] = sendList[neighborId].size(); // Get count if sends/receives
    for (set_iter i = receiveList[neighborId].begin(); i != receiveList[neighborId].end(); ++i, ++receiveEntryCount)
    {
      externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount; // The remote columns are indexed at end of internals
      A.localToExternalMap[localNumberOfRows + receiveEntryCount] = *i;
    }
    for (set_iter i = sendList[neighborId].begin(); i != sendList[neighborId].end(); ++i, ++sendEntryCount)
    {
      // if (geom.rank==1) HPCG_fout << "*i, globalToLocalMap[*i], sendEntryCount = " << *i << " " << A.globalToLocalMap[*i] << " " << sendEntryCount << endl;
      elementsToSend[sendEntryCount] = A.globalToLocalMap[*i]; // store local ids of entry to send
    }
  }

  // Convert matrix indices to local IDs
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for (local_int_t i = 0; i < localNumberOfRows; i++)
  {
    for (int j = 0; j < nonzerosInRow[map_l2a_A(A, i)]; j++)
    {
      global_int_t curIndex = mtxIndG[map_l2a_A(A, i) * numberOfNonzerosPerRow + j];
      int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
       
  
      if (A.geom->rank == rankIdOfColumnEntry)
      { // My column index, so convert to local index
        mtxIndL[i][j] = A.globalToLocalMap[curIndex];
      }
      else
      { // If column index is not a row index, then it comes from another processor
        mtxIndL[i][j] = externalToLocalMap[curIndex];
      }
    }
  }
  
  // Store contents in our matrix struct
  A.numberOfExternalValues = externalToLocalMap.size();
  A.localNumberOfColumns = A.localNumberOfRows + A.numberOfExternalValues;
  A.numberOfSendNeighbors = sendList.size();
  A.totalToBeSent = totalToBeSent;
  A.elementsToSend = elementsToSend;
  A.neighbors = neighbors;
  A.receiveLength = receiveLength;
  A.sendLength = sendLength;
  A.sendBuffer = sendBuffer;

  // ########## Data for partitioning algorithm
  partition_d *pt_data_local = (partition_d *)malloc(sizeof(partition_d));
  partition_d *pt_data_ext = (partition_d *)malloc(sizeof(partition_d));

  std::memcpy((void *)&pt_data_ext->receiveList, (void *)&receiveList, sizeof(receiveList));

  init_partition_data(A, pt_data_local, pt_data_ext);
  // ########## Data for partitioning algorithm

  // ########## Data to calculate mapping
  L2A_map *map_data = (L2A_map *)malloc(sizeof(L2A_map));

  map_data->localNumberOfRows = A.localNumberOfRows;
  map_data->offset = -1;
  map_data->offset_ext = -1;

  std::memcpy((void *)&map_data->localToGlobalMap, (void *)&A.localToGlobalMap, sizeof(A.localToGlobalMap));
  std::memcpy((void *)&map_data->localToExternalMap, (void *)&A.localToExternalMap, sizeof(A.localToExternalMap));

  A.mapping = map_data;
  // ########## Data to calculate mapping

  init_partitionings(A, pt_data_local, pt_data_ext);

#ifdef HPCG_DETAILED_DEBUG
  HPCG_fout << " For rank " << A.geom->rank << " of " << A.geom->size << ", number of neighbors = " << A.numberOfSendNeighbors << endl;
  for (int i = 0; i < A.numberOfSendNeighbors; i++)
  {
    HPCG_fout << "     rank " << A.geom->rank << " neighbor " << neighbors[i] << " send/recv length = " << sendLength[i] << "/" << receiveLength[i] << endl;
    for (local_int_t j = 0; j < sendLength[i]; ++j)
      HPCG_fout << "       rank " << A.geom->rank << " elementsToSend[" << j << "] = " << elementsToSend[j] << endl;
  }
#endif
  
  return;
}

#endif // HPCG_NO_LAIK
#endif // REPARTITION