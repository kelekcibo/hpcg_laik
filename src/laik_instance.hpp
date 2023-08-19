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

#include <vector>
#include <set>
#include <map>

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
typedef struct Local2Allocation_map
{
    // std::map<global_int_t, allocation_int_t> globalToAllocationMap; /* Mapping from Global to Allocation Indices */
    long long offset_halo; /* Offset into allocation buffer with external values */
    long long offset;      /* Offset into allocation buffer without external values */

    /* Mapping from Local to Global Indices */
    std::map<local_int_t, global_int_t>  localToExternalMap;        /* External global indices */
    std::vector<global_int_t> localToGlobalMap;                     /* Owned global indices */
    local_int_t localNumberOfRows;                                  /* Border between owned and external indices */

} L2A_map;

extern Laik_Data * x_vector;

extern Laik_Instance *hpcg_instance;
extern Laik_Group * world;

extern void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type);
extern void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type);
extern void laik_barrier(void);
extern void init_partitionings(pt_data *data, pt_data *data2);
extern void exchangeValues(bool halo);
extern void init_map_data(L2A_map * map_data);
extern allocation_int_t map_l2a(local_int_t local_index, bool halo);

#endif