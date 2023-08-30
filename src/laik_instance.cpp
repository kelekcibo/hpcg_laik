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
#include "Vector.hpp"

// should be initialized at the very beginning of the program. No use without init.
Laik_Instance *hpcg_instance; /* Laik instance during HPCG run */
Laik_Group *world;

// ############### Layer 0 Matrix Data needed to create partitionings and Laik_Data container
L2A_map *A0_map_data;
pt_data *A0_ext;
pt_data *A0_local;
// ############### Layer 1 Matrix Data needed to create partitionings and Laik_Data container
L2A_map *A1_map_data;
pt_data *A1_ext;
pt_data *A1_local;
// ############### Layer 2 Matrix Data needed to create partitionings and Laik_Data container
L2A_map *A2_map_data;
pt_data *A2_ext;
pt_data *A2_local;
// ############### Layer 3 Matrix Data needed to create partitionings and Laik_Data container
L2A_map *A3_map_data;
pt_data *A3_ext;
pt_data *A3_local;

/**
 * @brief Map the local index (mtxIndL) to the corresponding global index. Then map this global index to the allocation index
 *
 * @param[in] mapping used to calculate the allocation index
 * @param[in] localIndex which will be mapped to the corresponding allocation index
 * @param[in] halo determines the allocation buffer to be used
 * 
 * @return index to the allocation buffer
 * 
 * @see L2A_map in laik_instance.hpp
 */
allocation_int_t map_l2a(L2A_map * mapping, local_int_t local_index, bool halo)
{

    std::string a{""};
    a += std::to_string(local_index) + "   (LI) ->\t";

    global_int_t global_index;

    // Map local index to global index
    if (local_index < mapping->localNumberOfRows)
    {
        global_index = mapping->localToGlobalMap[local_index]; /* No need to handle errors, since we map values from mtxIndL only */
        a += std::to_string(global_index) + "   (GI) ->\t";
    }
    else
    {
        global_index = mapping->localToExternalMap[local_index];
        a += std::to_string(global_index) + "   (GI) ->\t";
    }

    // Map global index to allocation index
    // We have different allocation buffers depending on the current active partitioning
    allocation_int_t allocation_index;
    if (halo)
        allocation_index = global_index - mapping->offset_ext;
    else
        allocation_index = global_index - mapping->offset;

    a += std::to_string(allocation_index) + "   (AI)\n";

    // if (laik_myid(world) == 0)
        // printf("%s", a.data());

    return allocation_index;
}

/**
 * @brief Algorithm to partition an vector.
 *
 * local (Partitioning 1): Every process has exclusive access to its indices
 *
 * external (Partitioning 2): Every process has now access to needed external values as well
 *
 */
void partitioner_x(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    pt_data *data = (pt_data *)laik_partitioner_data(p->partitioner);
    Laik_Space * x_space = p->space;

    assert(data->size == laik_space_size(x_space));

    int rank = data->geom->rank;
    Laik_Range range;

    std::string a{""};
    a += "LAIK " + std::to_string(laik_myid(world)) + "\t[";

    for (long long i = 0; i < data->size; i++)
    {
        // assign every process its global part
        int proc = ComputeRankOfMatrixRow(*data->geom, i);

        laik_range_init_1d(&range, x_space, i, i + 1);
        laik_append_range(r, proc, &range, 1, 0);

        // Prepare data for mapping from global index to allocation index (lex_layout)

        if (data->offset == -1 && rank == proc)
        {
            data->offset = i;
        }
        a += std::to_string(ComputeRankOfMatrixRow(*data->geom, i)) + ", ";
    }

    a += "]\n";

    // printf("%s", a.data());

    // Specifying ranges which need to be exchanged
    if (data->halo)
    {
        data->offset_ext = data->offset;

        // process i needs access to external values
        typedef std::set<global_int_t>::iterator set_iter;

        for (int nb = 0; nb < data->numberOfNeighbours; nb++)
        {
            for (set_iter i = data->receiveList[data->neighbors[nb]].begin(); i != data->receiveList[data->neighbors[nb]].end(); ++i)
            {
                global_int_t globalIndex = *i;

                laik_range_init_1d(&range, x_space, globalIndex, globalIndex + 1);
                laik_append_range(r, rank, &range, 1, 0);

                if (globalIndex < data->offset_ext)
                {
                    data->offset_ext = globalIndex;
                }

                printf("I (%d) need to have access to global index %lld of x vector (updated by proc %d)\n", rank, globalIndex, data->neighbors[nb]);
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
 * @brief Fill the input vector with pseudo-random values.
 *
 * @param[inout] x_blob contains input vector
 *
 * @see FillRandomVector in Vector.hpp
 */
void fillRandomLaikVector(Laik_Blob *x_blob)
{
    double *base;
    uint64_t count;

    laik_get_map_1d(x_blob->values, 0, (void **)&base, &count);

    for (uint64_t i = 0; i < x_blob->localLength; i++)
    {
        base[map_l2a(x_blob->mapping, i, false)] = rand() / (double)(RAND_MAX) + 1.0;
    }
    return;
}

/**
 * @brief Fill the input vector with zero values.
 *
 * @param[inout] x_blob contains input vector
 *
 * @see ZeroVector in Vector.hpp
 */
void ZeroLaikVector(Laik_Blob * x_blob)
{
    double *base;
    uint64_t count;
    
    laik_get_map_1d(x_blob->values, 0, (void **)&base, &count);

    for (uint64_t i = 0; i < x_blob->localLength; i++)
    {
        base[map_l2a(x_blob->mapping, i, false)] = 0;
    }
    return;
}

/**
 * @brief Copy input vector to output Laik-vector.
 *
 * @param[in] v input vector
 * @param[in] x_blob contains output vector
 *
 * @see CopyVector in Vector.hpp
 */
void CopyVectorToLaikVector(Vector & v, Laik_Blob * x_blob)
{
    double *base;
    uint64_t count;
    laik_get_map_1d(x_blob->values, 0, (void **)&base, &count);

    const double *vv = v.values;

    for (uint64_t i = 0; i < x_blob->localLength; i++)
    {
        base[map_l2a(x_blob->mapping, i, false)] =  vv[i];
    }

    return;
}

/**
 * @brief Copy input Laik-vector to output vector.
 *
 * @param[in] x_blob contains input vector
 * @param[in] v output vector
 */
void CopyLaikVectorToVector(Laik_Blob *x_blob, Vector &v)
{
    double *base;
    uint64_t count;
    laik_get_map_1d(x_blob->values, 0, (void **)&base, &count);

    double *vv = v.values;

    for (uint64_t i = 0; i < x_blob->localLength; i++)
    {
        vv[i] = base[map_l2a(x_blob->mapping, i, false)];
    }

    return;
}

/**
 * @brief Populate an Laik_Blob and start with the partitioning such that every process has only access to local values.
 *
 * @param global_size of vector
 * @param local_size of vector for this process
 * @param map_data to calculate mapping from local to allocation indices (see L2A_map)
 * @param data_local needed for the partitioner algorithm (Access to local values)
 * @param data_ext needed for the partitioner algorithm (Access to local and external values)
 * @return populated blob
 */
Laik_Blob * init_blob(int64_t global_size, uint64_t local_size, L2A_map *map_data, pt_data *data_local, pt_data *data_ext)
{
    assert(global_size >= 0);
    assert(map_data != NULL);
    assert(data_local != NULL);
    assert(data_ext != NULL);
    assert(map_data->localNumberOfRows == local_size);

    Laik_Blob * blob = (Laik_Blob *) malloc(sizeof(Laik_Blob));

    Laik_Space * x_space = laik_new_space_1d(hpcg_instance, global_size);
    Laik_Data  * x_data = laik_new_data(x_space, laik_Double);

    // Initialize partitioners to partition data accordingly

    /* TODO. Do I need a seperate data buffer for every vector? */
    // pt_data * pt_data_ext_ = (pt_data *)malloc(sizeof(*data_ext));
    // pt_data * pt_data_local_ = (pt_data *)malloc(sizeof(*data_local));
    // std::memcpy((void *) pt_data_ext_, (void *) data_ext, sizeof(*data_ext));
    // std::memcpy((void *)pt_data_local_, (void *)data_local, sizeof(*data_local));

    /* Is it possible to reuse the same map_data object for each vector? Should be the case for vectors with same exchange data */

    Laik_Partitioner *x_extPR = laik_new_partitioner("x_extPR", partitioner_x, (void *)data_ext, LAIK_PF_None);
    Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR", partitioner_x, (void *)data_local, LAIK_PF_None);

    Laik_Partitioning * x_ext = laik_new_partitioning(x_extPR, world, x_space, NULL);
    Laik_Partitioning *x_local = laik_new_partitioning(x_localPR, world, x_space, NULL);


    blob->values = x_data;
    blob->x_ext = x_ext;
    blob->x_local = x_local;
    blob->mapping = map_data;
    blob->mapping->offset = data_local->offset;
    blob->mapping->offset_ext = data_ext->offset_ext;
    blob->localLength = local_size;

    // Start with x_local partitioning
    laik_switchto_partitioning(blob->values, blob->x_local, LAIK_DF_None, LAIK_RO_None);

    return blob;
}

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