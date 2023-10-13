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
 * @brief Helper function implements Allreduce / Broadcast
 *
 * @param sendBuf buffer to be sent
 * @param recvBuf received bytes will be stored in this buffer
 * @param n size of send-/recvbuffer
 * @param data_type of bytes in sendBuf
 * @param ro_type reduction operation after Allreduce / Broadcast
 * @param partitioner partitioner algorithm (to implement Broadcast and Allreduce)
 */
void laik_helper(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type, Laik_ReductionOperation ro_type, Laik_Partitioner *partitioner)
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

    // If broadcast: Root process has only access to container
    // else: everyone
    laik_switchto_new_partitioning(data, world, partitioner, LAIK_DF_None, LAIK_RO_None);

    // fill sendBuf into container
    uint64_t count;
    if (data_type == laik_Int32)
    {
        int *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        std::memcpy((void *)base, sendBuf, count * sizeof(int));
    }
    else if (data_type == laik_Double)
    {
        double *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        std::memcpy((void *)base, sendBuf, count * sizeof(double));
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        std::memcpy((void *)base, sendBuf, count * sizeof(uint64_t));
    }


    // Broadcast/Allreduce data with ro_type as reduction operation
    laik_switchto_new_partitioning(data, world, laik_All, LAIK_DF_Preserve, ro_type);

    // Store received result in recvBuf
    if (data_type == laik_Int32)
    {
        int *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        std::memcpy(recvBuf, (void *) base, count * sizeof(int32_t));
    }
    else if (data_type == laik_Double)
    {
        double *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        std::memcpy(recvBuf, (void *) base, count * sizeof(double));
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        std::memcpy(recvBuf, (void *)base, count * sizeof(uint64_t));
    }

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
void laik_allreduce(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type, Laik_ReductionOperation ro_type)
{
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

/**
 * @brief Synchronize all processes
 */
void laik_barrier()
{
    // arbitrary data
    int32_t data = 71;

    // Synchronize all processes by making use of an All-to-All-Reduction
    laik_helper((void *)&data, (void *)&data, 1, laik_Int32, LAIK_RO_None, laik_All);
    return;
}

/*
    Functions -END
*/
