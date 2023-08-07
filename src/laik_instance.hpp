/**
 * @brief Helper Functions + global variable instance + world (reason: using laik would result in more ifdef's)
 */
#ifndef LAIK_INSTANCE_HPP
#define LAIK_INSTANCE_HPP

extern "C" {
    #include <laik.h>
}

#include "Geometry.hpp" /* For local_int_t*/
#include <set>
#include <map>

/**
 * @brief Struct (data) needed for the partitioner algorithms
 * 
 * @param size size of vector x
 */
typedef struct partition_data
{
    local_int_t size;                       /* size of vector x */
    local_int_t local_portion;              /* Local portion of each proc in vector x */
    local_int_t numberOfExternalValues;     /* number of elements which need to be received */
    local_int_t * elementsToSend;           /* local index to vector x of elements to be sent */
    int numberOfNeighbours;                 /* count neighbors we need to exchange data with*/
    int * neighbors;                         /* Process ID's of neighbours */
    local_int_t * receiveLength;            /* number of elements to be received/sent by/to neighbours */
    std::map<int, std::set<global_int_t>> receiveList; /* Global Indices proc i needs to receive */
    bool halo; /* x will be partioned such that indices to external values are owned by other procs as well */
} pt_data;


extern Laik_Space * x_space;
extern Laik_Data * x_vector;
extern Laik_Data *x_vector_halo;
extern Laik_Instance * hpcg_instance;
extern Laik_Group * world;
extern Laik_Partitioning * x_pt;
extern Laik_Partitioning * x_halo_pt;

extern void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type);
extern void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type);
extern void laik_barrier(void);
extern void init_partitionings(pt_data *data);
extern void exchangeValues(bool halo);


#endif