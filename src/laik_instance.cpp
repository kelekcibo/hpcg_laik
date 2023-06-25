/**
 * @brief Definiton for helper functions
 * 
 */

#include "laik_instance.hpp"

#include <iostream>
#include <cstring>

// should be initialized at the very beginning of the program. No use without init. 
Laik_Instance *hpcg_instance;
Laik_Group *world;


/**
 * @brief Helper function implements Allreduce / Broadcast
 * 
 * @param partitioner needed to specify Allreduce or Broadcast
 */
void laik_helper(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type, Laik_ReductionOperation ro_type, Laik_Partitioner * partitioner)
{

    // reducing some buffer occurs very often in hpcg code base
    // implement it here once and call it from the hpcg code base

    // Definition and initialization of laik-specific data
    Laik_Space * space = laik_new_space_1d(hpcg_instance, n);
    Laik_Data * data = laik_new_data(space, data_type);

    // If broadcast: Root process has only access to container
    // else: everyone
    laik_switchto_new_partitioning(data, world, partitioner, LAIK_DF_None, LAIK_RO_None);

    // fill sendBuf into container
    uint64_t count;
    if(data_type == laik_Int32)
    {
        int * sendBuf_tmp = (int *) sendBuf;
        int * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            base[i] = sendBuf_tmp[i];

    }else if(data_type == laik_Double)
    { 
        double * sendBuf_tmp = (double *)sendBuf;
        double * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            base[i] = sendBuf_tmp[i];
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t * sendBuf_tmp = (uint64_t *)sendBuf;
        uint64_t * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            base[i] = sendBuf_tmp[i];

    }

    // Broadcast/Allreduce data with ro_type as reduction operation
    laik_switchto_new_partitioning(data, world, laik_All, LAIK_DF_Preserve, ro_type);

    // Store received result in recvBuf
    if (data_type == laik_Int32)
    {
        int * recvBuf_tmp = (int *)sendBuf;
        int * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            recvBuf_tmp[i] = base[i];

        std::memcpy(recvBuf, recvBuf_tmp, count * sizeof(int32_t));
    }
    else if (data_type == laik_Double)
    {
        double * recvBuf_tmp = (double *)sendBuf;
        double * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            recvBuf_tmp[i] = base[i];

        std::memcpy(recvBuf, recvBuf_tmp, count * sizeof(double)); // count should be equal to n
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t * recvBuf_tmp = (uint64_t *)sendBuf;
        uint64_t * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            recvBuf_tmp[i] = base[i];

        std::memcpy(recvBuf, (void *)recvBuf_tmp, count * sizeof(uint64_t));
    }

    // Free ressources
    laik_free(data);
    laik_free_space(space);
    
    return;
}

/**
 * @brief "Allreduce" a buffer based on the ro_type.
 *
 * @param sendBuf send buffer
 * @param recvBuf receive buffer
 * @param n size of send/recv buffer
 * @param data_type of the send buffer to be sent
 * @param ro_type type of reduction to be applied
 * @param allreduce if true, use allreduce
 */
void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type){
    laik_helper(sendBuf, recvBuf, n, data_type, ro_type, laik_All);
    return;
}

/**
 * @brief Broadcast a buffer to all processes in the world (Broadcast done by root process (rank==0)
 *
 * @param sendBuf send buffer
 * @param recvBuf receive buffer
 * @param n size of buffer
 * @param data_type of the buffer to be broadcasted
 */
void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type)
{
    laik_helper(sendBuf, recvBuf, n, data_type, LAIK_RO_None, laik_Master);
    return;
}
