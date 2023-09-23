/**
 * @file laik_instance.cpp
 * @brief Implementation of the definitions in th header file for leveraging LAIK API in the HPCG Application
 * @version 0.1
 * @date 2023-08-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cassert>

#include "laik_instance.hpp"
#include "Geometry.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"

// should be initialized at the very beginning of the program. No use without init.
Laik_Instance *hpcg_instance; /* Laik instance during HPCG run */
Laik_Group *world;

// Debug testsymmetry
Vector x_ncol_test;
Vector y_ncol_test;

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
    //     printf("%s", a.data());

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
    // a += "LAIK " + std::to_string(laik_myid(world)) + "\t[";

    for (long long i = 0; i < data->size; i++)
    {
        // assign every process its global part
        int proc = ComputeRankOfMatrixRow(*data->geom, i);

        laik_range_init_1d(&range, x_space, i, i + 1);
        laik_append_range(r, proc, &range, 1, 0);

        // Prepare data for mapping from global index to allocation index (lex_layout)

        if (data->offset == -1 && rank == proc)
            data->offset = i;

        // if(i== 0)
        //     printf("Global row 0 belongs to proc %d\n", ComputeRankOfMatrixRow(*data->geom, 0));

        // a += std::to_string(ComputeRankOfMatrixRow(*data->geom, i)) + ", ";
    }

    // a += "]\n";

    // printf("%s", a.data());

    // Specifying ranges which need to be exchanged
    if (data->halo)
    {
        // process i needs access to external values
        typedef std::set<global_int_t>::iterator set_iter;

        for (int nb = 0; nb < data->numberOfNeighbours; nb++)
        {
            for (set_iter i = data->receiveList[data->neighbors[nb]].begin(); i != data->receiveList[data->neighbors[nb]].end(); ++i)
            {
                global_int_t globalIndex = *i;

                laik_range_init_1d(&range, x_space, globalIndex, globalIndex + 1);
                laik_append_range(r, rank, &range, 1, 0);

                if (globalIndex < data->offset)
                {
                    data->offset = globalIndex;
                }

                // printf("I (%d) need to have access to global index %lld of x vector (updated by proc %d)\n", rank, globalIndex, data->neighbors[nb]);
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
 * @brief Fill the input vector with zero values.
 *
 * @param[inout] x_blob contains input vector
 *
 * @see ZeroVector in Vector.hpp
 */
void ZeroLaikVector(Laik_Blob *x, L2A_map *mapping)
{
    assert(x->localLength == mapping->localNumberOfRows);

    double *base;
    uint64_t count;
    
    laik_get_map_1d(x->values, 0, (void **)&base, &count);

    for (uint64_t i = 0; i < x->localLength; i++)
        base[map_l2a(mapping, i, false)] = 0;
    
    return;
}

/*!
  Multiply (scale) a specific vector entry by a given value.

  @param[inout] v Vector to be modified
  @param[in] index Local index of entry to scale
  @param[in] value Value to scale by
 */
void ScaleLaikVectorValue(Laik_Blob *v, local_int_t index, double value, L2A_map * mapping)
{
    assert(index >= 0 && index < v->localLength);
    assert(v->localLength == mapping->localNumberOfRows);

    double *vv;
    laik_get_map_1d(v->values, 0, (void **)&vv, 0);
    vv[map_l2a(mapping, index, false)] *= value;
    return;
}

/**
 * @brief Copy input Laik-vector to output Laik-vector.
 *
 * @param[in] x input vector
 * @param[in] y output vector
 */
void CopyLaikVectorToLaikVector(Laik_Blob *x, Laik_Blob *y, L2A_map *mapping)
{
    assert(x->localLength == y->localLength); /* they use the same mapping */
    assert(x->localLength == mapping->localNumberOfRows);

    double * xv;
    double * yv;

    laik_get_map_1d(x->values, 0, (void **)&xv, 0);
    laik_get_map_1d(y->values, 0, (void **)&yv, 0);

    for (uint64_t i = 0; i < x->localLength; i++)
        yv[map_l2a(mapping, i, false)] = xv[map_l2a(mapping, i, false)];

    return;
}

/**
 * @brief Fill the input vector with pseudo-random values.
 *
 * @param[inout] x contains input vector
 *
 * @see FillRandomVector in Vector.hpp
 */
void fillRandomLaikVector(Laik_Blob *x, L2A_map *mapping)
{
    assert(x != NULL);
    assert(mapping != NULL);
    assert(x->localLength == mapping->localNumberOfRows);

    double *xv;
    uint64_t count;

    laik_get_map_1d(x->values, 0, (void **)&xv, &count);

    for (uint64_t i = 0; i < x->localLength; i++)
        xv[map_l2a(mapping, i, false)] = rand() / (double)(RAND_MAX) + 1.0;
    
    return;
}

/**
 * @brief Copy input vector to output Laik-vector.
 *
 * @param[in] v input vector
 * @param[in] x contains output vector
 *
 * @see CopyVector in Vector.hpp
 */
void CopyVectorToLaikVector(Vector &v, Laik_Blob *x, L2A_map *mapping)
{
    assert(v.localLength >= x->localLength);
    assert(x->localLength == mapping->localNumberOfRows);

    double *xv;
    uint64_t count;
    laik_get_map_1d(x->values, 0, (void **)&xv, &count);

    const double *vv = v.values;

    for (uint64_t i = 0; i < x->localLength; i++)
        xv[map_l2a(mapping, i, false)] = vv[i];

    return;
}

/**
 * @brief Copy input Laik-vector to output vector.
 *
 * @param[in] x contains input vector
 * @param[in] v output vector
 */
void CopyLaikVectorToVector(const Laik_Blob *x, Vector &v, L2A_map *mapping)
{
    assert(x->localLength == v.localLength);
    assert(x->localLength == mapping->localNumberOfRows);

    double *xv;
    uint64_t count;
    laik_get_map_1d(x->values, 0, (void **)&xv, &count);

    double *vv = v.values;

    for (uint64_t i = 0; i < x->localLength; i++)
        vv[i] = xv[map_l2a(mapping, i, false)];

    return;
}
void CopyLaikVectorToVector(Laik_Blob *x, Vector &v, L2A_map *mapping)
{
    const Laik_Blob * x_const = x;
    CopyLaikVectorToVector(x_const, v, mapping);
}

/**
 * @brief Populate an Laik_Blob and start with the partitioning such that every process has only access to local values.
 *
 * @return populated blob
 */
Laik_Blob * init_blob(const SparseMatrix &A, bool exchangeHalo)
{
    Laik_Blob * blob = (Laik_Blob *) malloc(sizeof(Laik_Blob));

    blob->values = laik_new_data(A.space, laik_Double);
    blob->localLength = A.localNumberOfRows;
    blob->exchange = exchangeHalo;

    // Start with partitioning containing only access to local elements
    laik_switchto_partitioning(blob->values, A.local, LAIK_DF_None, LAIK_RO_None);

    return blob;
}

void init_partitionings(SparseMatrix &A, pt_data *local, pt_data *ext)
{
    A.space = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows);

    Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR", partitioner_x, (void *)local, LAIK_PF_None);
    Laik_Partitioner *x_extPR = laik_new_partitioner("x_extPR", partitioner_x, (void *)ext, LAIK_PF_None);


    A.local = laik_new_partitioning(x_localPR, world, A.space, NULL);
    A.ext = laik_new_partitioning(x_extPR, world, A.space, NULL);
    

    A.mapping->offset = local->offset;
    A.mapping->offset_ext = ext->offset;

    return;
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

// #######################################################################################
// Clean-up functions

void free_L2A_map(L2A_map *mapping)
{
    mapping->localNumberOfRows = 0;
    mapping->offset = 0;
    mapping->offset_ext = 0;

    mapping->localToExternalMap.clear();
    mapping->localToGlobalMap.clear();
    free((void*)mapping);
}

// #######################################################################################
// Debug functions

void compare2(double x, double y, bool doIO, allocation_int_t curIndex)
{
    double delta = std::abs(x - y);

    if(doIO)
        printf("map_l2a index: %lld\t| %.10f - %.10f | = %.10f\n", curIndex, x, y, delta);

    if (delta != 0.0)
    {
        if(doIO)
            printf("Difference is not tolerated: %.20f\nindex of map_l2a: %lld\n", delta, curIndex);
        assert(delta == 0);
    }
}

void compareResult(Vector &x, Laik_Blob *y, L2A_map *mapping, bool doIO)
{
    assert(x.localLength >= y->localLength); // Test vector lengths
    assert(y->localLength == mapping->localNumberOfRows);

    double *xv = x.values;
    double *yv;
    laik_get_map_1d(y->values, 0, (void **)&yv, 0);

    size_t length = y->localLength;

    for (size_t i = 0; i < length; i++)
    {
        double delta = std::abs(xv[i] - yv[map_l2a(mapping, i, false)]);
        if (doIO)
        // printf("Index %lld: Delta %.10f\n", i, delta);
        if (doIO)
            printf("xv[%ld]=%.10f\tyv_blob[%lld]=%.10f\n", i, xv[i], map_l2a(mapping, i, false), yv[map_l2a(mapping, i, false)]);
        if (delta != 0)
        {
            if (doIO)
                printf("Difference is not tolerated: %.20f\n", delta);
            assert(false);
        }
    }

    if (doIO)
        printf("Compare done\n");
}

void printResultVector(Vector &x)
{
    if (laik_myid(world) == 0)
    {
        printf("Print result of vector\n");
        double *xv = x.values;
        size_t length = x.localLength;

        for (size_t i = 0; i < length; i++)
            printf("xv[%ld]=%.10f\n", i, xv[i]);
    }
     
}

void printResultLaikVector(Laik_Blob *x, L2A_map *mapping)
{
    if (laik_myid(world) == 0)
        printf("Print result of Laik-vector\n");
    double *xv;
    laik_get_map_1d(x->values, 0, (void **)&xv, 0);

    size_t length = x->localLength;

    if (laik_myid(world) == 0)
        for (size_t i = 0; i < length; i++)
            printf("xv[%ld]=%.10f\n", i, xv[map_l2a(mapping, i, false)]);
}

void compareAfterExchange(Vector &x, Laik_Blob * y, SparseMatrix &A)
{
    const local_int_t nrow = A.localNumberOfRows;
    double **matrixDiagonal = A.matrixDiagonal; // An array of pointers to the diagonal entries A.matrixValues
}
