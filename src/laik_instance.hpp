/**
 * @file laik_instance.hpp
 * @brief Header file for leveraging LAIK API in the HPCG Application
 * @version 0.1
 * @date 2023-08-13
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef LAIK_INSTANCE_HPP
#define LAIK_INSTANCE_HPP


extern "C" {
    #include <laik.h>
}

#ifndef USE_LAIK
#define USE_LAIK
#endif

#include <vector>
#include <set>
#include <map>

struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;

#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "Geometry.hpp"

/**
 * @brief Struct (data) needed for the partitioner algorithms
 *
 * @param size                  size of vector x
 * @param elementsToSend        local index to vector x of elements to be sent (converting to global later)
 * @param geom                  Geometry of the problem
 * @param numberOfNeighbours    count neighbors we need to exchange data with
 * @param neighbors             Process ID's of neighbours
 * @param receiveLength         number of elements to be received/sent by/to neighbours
 * @param receiveList           Global Indices proc i needs to receive
 * @param localToGlobalMap      local-to-global mapping
 * @param halo                  x will be partioned such that indices to external values are owned by other procs as well
 * @param offset                offset to allocation buffer of partitioning with local values
 * @param offset_ext            offset to allocation buffer of partitioning with local+external values
 *
 */
typedef struct partition_data
{
    local_int_t size;                       /* size of vector x */
    local_int_t * elementsToSend;           /* local index to vector x of elements to be sent (converting to global later) */
    Geometry * geom;                        /* Geometry of the problem */
    int numberOfNeighbours;                 /* count neighbors we need to exchange data with*/
    int * neighbors;                        /* Process ID's of neighbours */
    local_int_t *receiveLength;             /* number of elements to be received/sent by/to neighbours */
    std::map<int, std::set<global_int_t>> receiveList;  /* Global Indices proc i needs to receive */
    std::vector<global_int_t> *localToGlobalMap;        /* local-to-global mapping */
    bool halo; /* x will be partioned such that indices to external values are owned by other procs as well */
    int offset; /* offset to allocation buffer of partitioning */
} pt_data;

typedef long long allocation_int_t;

/**
 * @brief Making use of the lex_layout, we need following mapping of indices:
 * 
 * Local indices -> Global indices -> Allocation indices
 * 
 * As LAIK allocates a buffer under the hood according to the lex layout, we need to specify a mapping for the corresponding
 * indices within that buffer.
 * 
 */
struct L2A_map
{
    // std::map<global_int_t, allocation_int_t> globalToAllocationMap; /* Mapping from Global to Allocation Indices */
    long long offset_ext; /* Offset into allocation buffer with external values */
    long long offset;      /* Offset into allocation buffer without external values */

    /* Mapping from Local to Global Indices */
    std::map<local_int_t, global_int_t>  localToExternalMap;        /* External global indices */
    std::vector<global_int_t> localToGlobalMap;                     /* Owned global indices */
    local_int_t localNumberOfRows;                                  /* Border between owned and external indices */

};

struct Laik_Blob
{
    Laik_Data * values;
    uint64_t localLength;

    bool exchange; /* This blob will exchange values if true */
};


extern Laik_Instance *hpcg_instance;
extern Laik_Group * world;

extern void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type);
extern void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type);
extern void laik_barrier(void);

extern void fillRandomLaikVector(Laik_Blob *x, L2A_map *mapping);
extern void ZeroLaikVector(Laik_Blob *x, L2A_map *mapping);
extern void CopyLaikVectorToLaikVector(Laik_Blob *x, Laik_Blob *y, L2A_map *mapping);
extern void CopyVectorToLaikVector(Vector &v, Laik_Blob *x_blob, L2A_map *mapping);
extern void CopyLaikVectorToVector(Laik_Blob *x_blob, Vector &v, L2A_map *mapping);
extern void ScaleLaikVectorValue(Laik_Blob *v, local_int_t index, double value);

extern void init_partitionings(SparseMatrix &A, pt_data * local, pt_data * ext);
extern Laik_Blob *init_blob(const SparseMatrix &A, bool exchangeHalo);

extern allocation_int_t map_l2a(L2A_map *mapping, local_int_t local_index, bool halo);

#endif