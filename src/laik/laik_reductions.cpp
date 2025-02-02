/**
 * @file laik_reductions.cpp
 * @brief Using LAIK for Allreduction-/Broadcastoperations
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

/*
    Includes
*/
#include <cstring>

#include "laik_reductions.hpp"
/*
    Includes -END
*/

/*
    Functions
*/

/* LAIK has a hard limit of “Laik_Data” objects. Reuse “Laik_Data” objects with same size and Laik_Type */
std::map<std::pair<int, Laik_Type *>, Laik_Data *> data_objects;

/**
 * @brief Helper function to avoid redundant code for operations like Allreduce / Broadcast
 *
 * @param[in] sendBuf buffer to be sent
 * @param[out] recvBuf received bytes will be stored in this buffer
 * @param[in] n size of send-/recvbuffer
 * @param[in] data_type of bytes in sendBuf
 * @param[in] ro_type reduction operation after Allreduce / Broadcast
 * @param[in] partitioner partitioner algorithm used for the first partitioning
 * @param[in] partitioner2 partitioner algorithm used for the second partitioning
 */
void laik_helper(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type, Laik_Partitioner * partitioner, Laik_Partitioner * partitioner2)
{
    // Definition and initialization of laik-specific data (Reuse of laik-specific data, if used more than once)
    Laik_Data *data;
    Laik_Space *space;
    if (data_objects.find({n, data_type}) != data_objects.end())
    {
        data = data_objects[{n, data_type}];
        // space = laik_data_get_space(data);
    }
    else
    {
        space = laik_new_space_1d(hpcg_instance, n);
        data = laik_new_data(space, data_type);
        data_objects.insert({{n, data_type} , data});
    }
 
    // Switch to partitioning created by running "partitioner"
    laik_switchto_new_partitioning(data, world, partitioner, LAIK_DF_None, LAIK_RO_None);

    // fill sendBuf into container
    uint64_t count;
    if (data_type == laik_Int32)
    {
        int *base; laik_get_map_1d(data, 0, (void **)&base, &count);
        std::memcpy((void *)base, sendBuf, count * sizeof(int));
    }
    else if (data_type == laik_Double)
    {
        double *base; laik_get_map_1d(data, 0, (void **)&base, &count);
        std::memcpy((void *)base, sendBuf, count * sizeof(double));
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t *base; laik_get_map_1d(data, 0, (void **)&base, &count);
        std::memcpy((void *)base, sendBuf, count * sizeof(uint64_t));
    }

    // Switch to partitioning created by running "partitioner2", exchange values, and apply reduction "ro_type"
    laik_switchto_new_partitioning(data, world, partitioner2, LAIK_DF_Preserve, ro_type);

    // Store received result in recvBuf
    if (data_type == laik_Int32)
    {
        int *base; laik_get_map_1d(data, 0, (void **)&base, &count);
        std::memcpy(recvBuf, (void *) base, count * sizeof(int32_t));
    }
    else if (data_type == laik_Double)
    {
        double *base; laik_get_map_1d(data, 0, (void **)&base, &count);
        std::memcpy(recvBuf, (void *) base, count * sizeof(double));
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t *base; laik_get_map_1d(data, 0, (void **)&base, &count);
        std::memcpy(recvBuf, (void *)base, count * sizeof(uint64_t));
    }

    return;
}

/**
 * @brief "Allreduce" a buffer based on the ro_type.
 *
 * @param[in] sendBuf send buffer
 * @param[out] recvBuf receive buffer
 * @param[in] n size of send/recv buffer
 * @param[in] data_type of the send buffer to be sent
 * @param[in] ro_type type of reduction to be applied
 * @param[in] allreduce if true, use allreduce
 */
void laik_allreduce(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type, Laik_ReductionOperation ro_type)
{
    laik_helper(sendBuf, recvBuf, n, data_type, ro_type, laik_All, laik_All);
    return;
}

/**
 * @brief Broadcast a buffer to all processes in the world (Broadcast done by root process (rank==0)
 *
 * @param[in] sendBuf send buffer
 * @param[out] recvBuf receive buffer
 * @param[in] n size of buffer
 * @param[in] data_type of the buffer to be broadcasted
 */
void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type)
{
    laik_helper(sendBuf, recvBuf, n, data_type, LAIK_RO_None, laik_Master, laik_All);
    return;
}

/**
 * @brief Synchronize all processes
 */
void laik_barrier()
{
    // arbitrary data
    int32_t data = 71;

    // Synchronize all processes by making use of an All-to-All-Reduction
    laik_helper((void *)&data, (void *)&data, 1, laik_Int32, LAIK_RO_None, laik_All, laik_All);
    return;
}

/**
 * @brief  Broadcast a buffer to certain processes in the world (Broadcast done by root process (rank==0)
 *
 * Not in use yet. Using broadcast instead.
 *
 * @param[in] sendBuf send buffer
 * @param[out] recvBuf receive buffer
 * @param[in] n size of buffer
 * @param[in] data_type of the buffer to be broadcasted
 * @param[in] old_size of the world
 */
void laik_partial_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type, int old_size)
{
    printf("WARNING: This is only enabled, if repartitioning should be done!\n");    
#ifdef REPARTITION
    Laik_Partitioner *partitioner = laik_new_partitioner("Partial_broadcast", new_joining_procs, (void *)&old_size, LAIK_PF_None);
    laik_helper(sendBuf, recvBuf, n, data_type, LAIK_RO_None, laik_Master, partitioner);
#endif
    return;
}

/*
    Functions -END
*/
