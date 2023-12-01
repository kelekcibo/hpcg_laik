/**
 * @file laik_x_vector.cpp
 * @brief Implementation - Everything needed for using Laik Vectors
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

/*
    Includes
*/
    #include "laik_x_vector.hpp"
// forw. decl.
struct Local2Allocation_map_x;
typedef Local2Allocation_map_x L2A_map;
// forw. decl.

/*
    Includes -END
*/

/*
    Functions needed to exchange values via LAIK
*/
/**
 * @brief Algorithm to partition an Laik Vector.
 *
 * local (Partitioning 1): Every process has exclusive access to its owned indices
 *
 * external (Partitioning 2): Every process has now access to needed external values as well
 *
 */
void partitioner_alg_for_x_vector(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    partition_d *data = (partition_d *)laik_partitioner_data(p->partitioner);
    Laik_Space *x_space = p->space;

    // printf("data->size=%d\tspace->size=%ld\n", data->size, laik_space_size(x_space));
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
 * @brief Map the local index (mtxIndL) to the corresponding global index. Then map this global index to the allocation index. 
 * 
 * This extra step has to be done since Laik makes use of the lex_layout under the hood.  
 * v1.2 will hopefully implement the custom hpcg_layout, so that there is no need for this anymore
 *
 * @param[in] mapping used to calculate the allocation index
 * @param[in] localIndex which will be mapped to the corresponding allocation index
 * @param[in] halo determines the allocation buffer to be used
 *
 * @return index to the allocation buffer
 *
 * @see L2A_map 
 */
allocation_int_t map_l2a_x(L2A_map *mapping, local_int_t local_index, bool halo)
{
    global_int_t global_index;

    // Map local index to global index
    if (local_index < mapping->localNumberOfRows)
        global_index = mapping->localToGlobalMap[local_index]; /* No need to handle errors, since we map values from mtxIndL only */
    else
        global_index = mapping->localToExternalMap[local_index];

    // Map global index to allocation index
    // We have different allocation buffers depending on the current active partitioning
    allocation_int_t allocation_index;
    if (halo)
        allocation_index = global_index - mapping->offset_ext;
    else
        allocation_index = global_index - mapping->offset;

    return allocation_index;
}

/**
 * @brief Populate an Laik_Blob and start with the partitioning such that every process has only access to local values.
 *
 * @param A SparseMatrix
 * @return populated blob
 */
Laik_Blob *init_blob(const SparseMatrix &A, bool exchangesValues)
{
    Laik_Blob *blob = (Laik_Blob *)malloc(sizeof(Laik_Blob));

    blob->values = laik_new_data(A.space, laik_Double);
    blob->localLength = A.localNumberOfRows;
    blob->exchangesValues = exchangesValues;

    // debug
    std::string b{"local "};
    b += to_string(laik_myid(world));
    std::string c{"ext "};
    c += to_string(laik_myid(world));

    laik_partitioning_set_name(A.local, b.data());
    laik_partitioning_set_name(A.ext, c.data());
    std::string a{"Blob "};
    a += to_string(laik_myid(world)) + "\t";
    a.append("b_l");
    laik_data_set_name(blob->values, a.data());

    // Preparation to use the vector layout
    laik_data_set_layout_factory(blob->values, laik_new_layout_vector);
    laik_data_set_layout_flag(blob->values, LAIK_Vector_Layout);

    Laik_vector_data *layout_data = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));
    layout_data->localLength = blob->localLength;
    layout_data->numberOfExternalValues = A.numberOfExternalValues;
    layout_data->id = laik_myid(world);
    laik_data_set_layout_data(blob->values, layout_data);

    // New joining procs will switch later
    if (laik_phase(hpcg_instance) == 0)
    {
        // First, switchto A.ext then switch to A.local immediately
        // Optimization, as LAIK will create a buffer under the hood and reuse the buffer, if we switch to A.ext first and then to A.local
        // Not all vectors have external values, handle that case,so we only switch to local pt
        if(blob->exchangesValues)
            laik_switchto_partitioning(blob->values, A.ext, LAIK_DF_None, LAIK_RO_None);

        // Start with partitioning containing only access to local elements
        laik_switchto_partitioning(blob->values, A.local, LAIK_DF_None, LAIK_RO_None);
    }

    return blob;
}

/**
 * @brief Initialize partition_d data needed to partition the Laik vectors.
 *
 * Needs to be called before init_partitionings()
 *
 * @param A SparseMatrix
 * @param local partitioning data
 * @param ext partitioning data
 */
void init_partition_data(SparseMatrix &A, partition_d *local, partition_d *ext)
{
    ext->size = A.totalNumberOfRows;
    ext->geom = A.geom;
    ext->numberOfNeighbours = A.numberOfSendNeighbors;
    ext->neighbors = A.neighbors;
    ext->localToGlobalMap = &A.localToGlobalMap;
    ext->elementsToSend = A.elementsToSend;
    ext->receiveLength = A.receiveLength;
    ext->halo = true;
    ext->offset = -1;

    local->halo = false;
    local->geom = ext->geom;
    local->size = ext->size;
    local->offset = -1;
    /* These values are not needed for the 2nd partitioning */
    local->neighbors = NULL;
    local->localToGlobalMap = NULL;
    local->elementsToSend = NULL;
    local->receiveLength = NULL;
    local->numberOfNeighbours = -1;

    return;
}

/**
 * @brief Initialize partitionings needed to partition the Laik vectors.
 *
 * This partitionings are needed to partition the x vectors and to be able to exchange values via LAIK
 *
 * @param A SparseMatrix
 * @param local partitioning
 * @param ext partitioning
 */
void init_partitionings(SparseMatrix &A, partition_d *local, partition_d *ext)
{
    Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR", partitioner_alg_for_x_vector, (void *)local, LAIK_PF_None);
    Laik_Partitioner *x_extPR = laik_new_partitioner("x_extPR", partitioner_alg_for_x_vector, (void *)ext, LAIK_PF_None);

    A.local = laik_new_partitioning(x_localPR, world, A.space, NULL);
    A.ext = laik_new_partitioning(x_extPR, world, A.space, NULL);

    A.mapping->offset = local->offset;
    A.mapping->offset_ext = ext->offset;

    return;
}

/*
    Functions needed to exchange values via LAIK -END
*/

/*
    Operations on laik vectors
*/

/**
 * @brief Fill the input vector with zero values.
 *
 * @param[inout] x_blob contains input vector
 *
 * @see ZeroVector in Vector.hpp
 */
void ZeroLaikVector(Laik_Blob *x)
{
    assert(x);

    double *base;
    uint64_t count;

    laik_get_map_1d(x->values, 0, (void **)&base, &count);

    for (uint64_t i = 0; i < x->localLength; i++)
        base[i] = 0;

    return;
}

/*!
  Multiply (scale) a specific vector entry by a given value.

  @param[inout] v Vector to be modified
  @param[in] index Local index of entry to scale
  @param[in] value Value to scale by
 */
void ScaleLaikVectorValue(Laik_Blob *v, local_int_t index, double value)
{
    assert(v);
    assert(index >= 0 && index < v->localLength);

    double *vv;
    laik_get_map_1d(v->values, 0, (void **)&vv, 0);
    vv[index] *= value;
    return;
}

/**
 * @brief Copy input Laik-vector to output Laik-vector.
 *
 * @param[in] x input vector
 * @param[in] y output vector
 */
void CopyLaikVectorToLaikVector(Laik_Blob *x, Laik_Blob *y)
{
    assert(x->localLength == y->localLength); /* they use the same mapping */

    double *xv;
    double *yv;

    laik_get_map_1d(x->values, 0, (void **)&xv, 0);
    laik_get_map_1d(y->values, 0, (void **)&yv, 0);

    for (uint64_t i = 0; i < x->localLength; i++)
        yv[i] = xv[i];

    return;
}

/**
 * @brief Fill the input vector with pseudo-random values.
 * 
 * @param[inout] x contains input vector
 *
 * @see FillRandomVector in Vector.hpp
 */
void fillRandomLaikVector(Laik_Blob *x)
{
    assert(x != NULL);

    double *xv;
    uint64_t count;
    laik_get_map_1d(x->values, 0, (void **)&xv, &count);
    for (uint64_t i = 0; i < x->localLength; i++)
        // xv[i] = rand() / (double)(RAND_MAX) + 1.0;
        xv[i] = i + 1.0;
        // xv[i] = laik_myid(world);
        // TODO. Hardcoded values for now to test it with the original application
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
void CopyVectorToLaikVector(Vector &v, Laik_Blob *x)
{
    assert(v.localLength >= x->localLength);

    double *xv;
    uint64_t count;
    laik_get_map_1d(x->values, 0, (void **)&xv, &count);

    const double *vv = v.values;

    for (uint64_t i = 0; i < x->localLength; i++)
        xv[i] = vv[i];

    return;
}

/**
 * @brief Copy input Laik-vector to output vector.
 *
 * @param[in] x contains input vector
 * @param[in] v output vector
 */
void CopyLaikVectorToVector(const Laik_Blob *x, Vector &v)
{
    assert(x->localLength == v.localLength);

    double *xv;
    uint64_t count;
    laik_get_map_1d(x->values, 0, (void **)&xv, &count);

    double *vv = v.values;

    for (uint64_t i = 0; i < x->localLength; i++)
        vv[i] = xv[i];

    return;
}

/**
 * @brief Copy input Laik-vector to output vector.
 * @see CopyLaikVectorToVector
 */
void CopyLaikVectorToVector(Laik_Blob *x, Vector &v)
{
    const Laik_Blob *x_const = x;
    CopyLaikVectorToVector(x_const, v);
    return;
}
/*
    Operations on laik vectors -END
*/

/*
    Clean-up functions
*/
/**
 * @brief Deallocate L2A_map Struct
 *
 * @param mapping to be deallocated
 */
void free_L2A_map(L2A_map *mapping)
{
    if (mapping != NULL)
    {
        mapping->localNumberOfRows = 0;
        mapping->offset = 0;
        mapping->offset_ext = 0;

        if (!mapping->localToExternalMap.empty())
        {
            mapping->localToExternalMap.clear();
            // Fixed segfault for now see also in laik_reaprtition:cpp in update_partitionings_x
            mapping->localToExternalMap[-1] = 144;
        }
        if (!mapping->localToGlobalMap.empty())
            mapping->localToGlobalMap.clear();

        free((void *)mapping);
    }
    return;
}

/**
 * @brief Deallocate a Laik vector
 *
 * @param x to be deallocated
 */
void DeleteLaikVector(Laik_Blob *x)
{
    x->localLength = 0;
    laik_free(x->values);
    x->values = NULL;
}

/*
    Clean-up functions -END
*/