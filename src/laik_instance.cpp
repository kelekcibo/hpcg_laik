/**
 * @file laik_instance.cpp
 * @brief Implementation of the definitions in th header file for leveraging LAIK API in the HPCG Application
 * @version 0.1
 * @date 2023-08-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <iostream>
#include <cstring>
#include <cassert>

#include "laik_instance.hpp"
#include "Geometry.hpp"

// should be initialized at the very beginning of the program. No use without init.
Laik_Instance *hpcg_instance; /* Laik instance during HPCG run */
Laik_Group *world;

/* Space and containers */
Laik_Space *x_space; /* space of vector x */
Laik_Data *x_vector; /* container / x vector */

/* Partitioning to switch between */
Laik_Partitioning *x_pt;
Laik_Partitioning *x_halo_pt;

/**
 * @brief Algorithm to partition the x vector.
 *
 * x_pt (Partitioning 1): Every process has exclusive access to its indices
 *
 * x_halo_pt (Partitioning 2): Every process has now access to needed external values as well
 *
 */
void partitioner_x(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    pt_data *data = (pt_data *)laik_partitioner_data(p->partitioner);
    

    Laik_Range range;

    // std::string a{""};

    // a += "LAIK " + std::to_string(laik_myid(world)) + "\t[";

    for (uint32_t i = 0; i < data->size; i++)
    {
        // assign every process its global part
        laik_range_init_1d(&range, x_space, i, i + 1);
        laik_append_range(r, ComputeRankOfMatrixRow(*data->geom, i), &range, 1, 0);

        // a += std::to_string(ComputeRankOfMatrixRow(*data->geom, i)) + ", ";
    }

    // a += "]\n";

    // printf("%s", a.data());
    
    // Specifying ranges which need to be exchanged
    if (data->halo)
    {
        // process i needs access to external values
        int rank = data->geom->rank;
        typedef std::set<global_int_t>::iterator set_iter;

        for (int nb = 0; nb < data->numberOfNeighbours; nb++)
        {
            for (set_iter i = data->receiveList[data->neighbors[nb]].begin(); i != data->receiveList[data->neighbors[nb]].end(); ++i)
            {
                int globalIndex = *i;

                laik_range_init_1d(&range, x_space, globalIndex, globalIndex + 1);
                laik_append_range(r, rank, &range, 1, 0);

                // printf("I (%d) need to have access to global index %d of x vector (updated by proc %d)\n", rank, globalIndex, data->neighbors[nb]);
            }
        }


        // proc i needs to know which processes need access to its own ranges
        local_int_t index_elementsToSend = 0;

        // Process i needs to give access to indices to neighbour procs
        for (int nb = 0; nb < data->numberOfNeighbours; nb++)
        {
            local_int_t numbOfElementsToSend = data->receiveLength[nb]; // due to symmetry, this is equal to sendLength[nb]

            for (local_int_t i = 0; i < numbOfElementsToSend; i++)
            {

                global_int_t j = (*data->localToGlobalMap)[data->elementsToSend[index_elementsToSend++]]; // neighbour n needs this value from proc i
             
                assert(ComputeRankOfMatrixRow(*(data->geom), j) == data->geom->rank);
                // printf("I (%d) need to give access to proc %d at global offset %llu\n", data->geom->rank, data->neighbors[nb], j);

                laik_range_init_1d(&range, x_space, j, j + 1);
                laik_append_range(r, data->neighbors[nb], &range, 1, 0);
            }
        }
    }
}

/**
 * @brief Initialize the two partitiongs to switch between them (exchange data between processes)
 *
 * @param data,data2 necessary information needed for the partitioner algorithm
 */
void init_partitionings(pt_data *data, pt_data *data2)
{
    x_space = laik_new_space_1d(hpcg_instance, data->size);
    x_vector = laik_new_data(x_space, laik_Double);

    // Initialize partitioners to partition data accordingly
    Laik_Partitioner *x_halo_pr = laik_new_partitioner("x_pt_halo", partitioner_x, (void *)data, LAIK_PF_None);
    Laik_Partitioner *x_pr = laik_new_partitioner("x_pt", partitioner_x, (void *)data2, LAIK_PF_None);

    // Initialize partitionings
    x_halo_pt = laik_new_partitioning(x_halo_pr, world, x_space, NULL);
    x_pt = laik_new_partitioning(x_pr, world, x_space, NULL);

    // Start with / Switch to partitioning x_pt
    exchangeValues(false);
}

/**
 * @brief Switch between partitionings to exchange data
 *
 * @param halo allow access to external values, if true
 */
void exchangeValues(bool halo)
{

    if (halo)
    {
        laik_switchto_partitioning(x_vector, x_halo_pt, LAIK_DF_Preserve, LAIK_RO_None);
        return;
    }

    laik_switchto_partitioning(x_vector, x_pt, LAIK_DF_Preserve, LAIK_RO_None);
}

/* LAIK has a hard limit of “Laik_Data” objects. Reuse “Laik_Data” objects with same size and Laik_Type */
std::map<std::pair<int, Laik_Type *>, Laik_Data *> data_objects;

/**
 * @brief Helper function implements Allreduce / Broadcast
 *
 * @param partitioner needed to specify Allreduce or Broadcast
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

        assert(count == n); // elements in buffer; should be equal
        std::memcpy(recvBuf, (void *) base, count * sizeof(int32_t));
    }
    else if (data_type == laik_Double)
    {
        double *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        assert(count == n); // elements in buffer; should be equal
        std::memcpy(recvBuf, (void *) base, count * sizeof(double));
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t *base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        assert(count == n); // elements in buffer; should be equal
        std::memcpy(recvBuf, (void *)base, count * sizeof(uint64_t));
    }

    // Free ressources
    // laik_free(data);
    // laik_free_space(space);

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
 *
 */
void laik_barrier()
{
    // arbitrary data
    int32_t data = 71;

    // Synchronize all processes by making use of an All-to-All-Reduction
    laik_helper((void *)&data, (void *)&data, 1, laik_Int32, LAIK_RO_None, laik_All);
    return;
}