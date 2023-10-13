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

#ifdef REPARTITION
#include "hpcg.hpp"
#include "Geometry.hpp"
#include "GenerateProblem.hpp"
#include "SetupHalo.hpp"
#include "GenerateCoarseProblem.hpp"
#include "GenerateGeometry.hpp"
#include "CheckProblem.hpp"
#endif

// debug forw. delc
void printSPM(SparseMatrix *spm, int coarseLevel);
// debug forw. delc

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

        /* DEBUG */
        if(rank == -1)
        {
            assert(proc == 0);
            assert((proc >= 0) && (proc < (int)r->list->tid_count));
            // a += std::to_string(ComputeRankOfMatrixRow(*data->geom, i));
            // printf("%s", a.data());
        }
        /* DEBUG */

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
    // For the case, if repartitioning is enabled. We do not need a new space
    if(A.space == NULL)
    {
        A.space = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows);
        if(A.localNumberOfRows==8192)
            printf("new_space should not be done");
    }

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

#ifdef REPARTITION
// #######################################################################################
// Needed functions/variables for shrink/expand feature

HPCG_Params hpcg_params; /* Global copy of params in main function */

/**
 * @brief Re-run the setup phase to update parameters and data for partitioner with new config (size of world / num of procs).
 *
 * @param A SparseMatrix
 */
void re_setup_problem(SparseMatrix &A)
{

    // Need this variables for generateGeometry, before we delete them
    global_int_t old_gnx = A.geom->gnx;
    global_int_t old_gny = A.geom->gny;
    global_int_t old_gnz = A.geom->gnz;

    // Delete old setup; Everything except next layer matrix, space and mgdata
    DeleteGeometry(*(A.geom));
    DeleteMatrix_repartition(A, false);

    // Update parameters of params struct (hpcg_params is a copy)
    hpcg_params.comm_rank = laik_myid(world);
    hpcg_params.comm_size = laik_size(world);

    // TODO. Send params_struct to new procs; old procs already have them
    // Broadcast updated values. New joining processes will need them, if broadcast was done in init
    // laik_broadcast(iparams, iparams, nparams, laik_Int32);

    // Construct the new geometry and linear system
    Geometry * new_geom = new Geometry;

    // need old gnx to calculate new nx, ny, nz. As npx, npy, npz may change due to new world size
    new_geom->gnx = old_gnx; 
    new_geom->gny = old_gny;
    new_geom->gnz = old_gnz;

    int dynamicCalculation = hpcg_params.numThreads;
    // use numThreads param to jump into the if to calculate new nx,ny,nz for now
    GenerateGeometry(hpcg_params.comm_size, hpcg_params.comm_rank, -1, hpcg_params.pz, hpcg_params.zl, hpcg_params.zu, hpcg_params.nx, hpcg_params.ny, hpcg_params.nz, hpcg_params.npx, hpcg_params.npy, hpcg_params.npz, new_geom);
    new_geom->numThreads = dynamicCalculation;

    // should be the same problem size
    assert(old_gnx == new_geom->gnx);
    assert(old_gny == new_geom->gny);
    assert(old_gnz == new_geom->gnz);

    // Only new geom will be set, other variables are assigned in generateProblem and SetupHalo
    A.geom = new_geom;


    if(hpcg_params.comm_rank >= 0)
    {
        // Only active procs need to generate problem again.
        GenerateProblem(A, 0, 0, 0);

        // save pointer to old partitionings, since they need to be freed after call to re_switch_LaikVectors
        assert(A.old_local == NULL);
        assert(A.old_ext == NULL);
        A.old_local = A.local;
        A.old_ext = A.ext;

        SetupHalo(A);

        int numberOfMgLevels = 4; // Number of levels including first
        SparseMatrix *curLevelMatrix = &A;
        for (int level = 1; level < numberOfMgLevels; ++level)
        {
            // std::cout << "\nCoarse Problem level " << level << std::endl;

            /* if all rows are mapped on proc 0 */
            for (int i = 0; i < curLevelMatrix->totalNumberOfRows; i++)
                assert(ComputeRankOfMatrixRow(*(curLevelMatrix->geom), i) == 0); // TEST CASE. -np 2 => shrinking to world size 1

            GenerateCoarseProblem(*curLevelMatrix);
            curLevelMatrix = curLevelMatrix->Ac; // Make the just-constructed coarse grid the next level
            // #### Debug
            // printSPM(curLevelMatrix, level);
            // #### Debug
        }

        curLevelMatrix = &A;
        Vector *curb = 0;
        Vector *curx = 0;
        Vector *curxexact = 0;
        for (int level = 0; level < numberOfMgLevels; ++level)
        {
            // std::cout << "\nCoarse Problem level " << level << std::endl;
            CheckProblem(*curLevelMatrix, curb, curx, curxexact);
            curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
            curb = 0;                            // No vectors after the top level
            curx = 0;
            curxexact = 0;
        }
    } 
    else
    {
        // Processes, which are removed from the world do not need to re-generate the problem
        // They need to create the local partitiniong, so they send their values to the active procs
        // But before that, create the geometry for the next three layers of the matrix

        int numberOfMgLevels = 4; // Number of levels including first
        SparseMatrix *curLevelMatrix = &A;

        // Need to update the TotalNumberOfRows value since we do not call generateProblem on kicked out procs
        A.totalNumberOfRows = A.geom->gnx * A.geom->gny * A.geom->gnz; /* Not needed as problem size should stay the same */

        for (int level = 1; level < numberOfMgLevels; ++level)
        {
    
            /* Next three layers do not have the new geometry yet */
            Geometry *geomc = new Geometry;
            Geometry *old_geomc = curLevelMatrix->Ac->geom;
            assert(old_geomc != NULL);

            local_int_t zlc = 0; // Coarsen nz for the lower block in the z processor dimension
            local_int_t zuc = 0; // Coarsen nz for the upper block in the z processor dimension
            int pz = curLevelMatrix->geom->pz;
            if (pz > 0)
            {
                zlc = curLevelMatrix->geom->partz_nz[0] / 2; // Coarsen nz for the lower block in the z processor dimension
                zuc = curLevelMatrix->geom->partz_nz[1] / 2; // Coarsen nz for the upper block in the z processor dimension
            }

            global_int_t nxf = curLevelMatrix->geom->nx;
            global_int_t nyf = curLevelMatrix->geom->ny;
            global_int_t nzf = curLevelMatrix->geom->nz;

            local_int_t nxc, nyc, nzc; // Coarse nx, ny, nz
            assert(nxf % 2 == 0);
            assert(nyf % 2 == 0);
            assert(nzf % 2 == 0); // Need fine grid dimensions to be divisible by 2
            nxc = nxf / 2;
            nyc = nyf / 2;
            nzc = nzf / 2;

            GenerateGeometry(curLevelMatrix->geom->size, curLevelMatrix->geom->rank, curLevelMatrix->geom->numThreads, curLevelMatrix->geom->pz, zlc, zuc, nxc, nyc, nzc, curLevelMatrix->geom->npx, curLevelMatrix->geom->npy, curLevelMatrix->geom->npz, geomc);
            
            // should be the same problem size
            assert(old_geomc->gnx == geomc->gnx);
            assert(old_geomc->gny == geomc->gny);
            assert(old_geomc->gnz == geomc->gnz);

            DeleteGeometry(*old_geomc);

            curLevelMatrix->Ac->geom = geomc;

            // Need to update the TotalNumberOfRows value since we do not call generateProblem on kicked out procs
            /* Not needed as problem size should be the same */
            curLevelMatrix->Ac->totalNumberOfRows = geomc->gnx * geomc->gny * geomc->gnz;

            curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
        }

        // Create the local partitioning for all Matrix levels
        curLevelMatrix = &A;
        for (int level = 0; level < numberOfMgLevels; ++level)
        {
            // std::cout << "\nCoarse Problem level " << level << std::endl;
         
            /* DEBUG */
            /* if all rows are mapped on proc 0 */
            for (int i = 0; i < curLevelMatrix->totalNumberOfRows; i++)
            {
                int rank = ComputeRankOfMatrixRow(*(curLevelMatrix->geom), i);
                // TEST CASE. -np 2 => shrinking to world size 1
                if (rank != 0)
                {
                    std::string str{"Rank not equal to 0: " + std::to_string(rank)};
                    std::cout << str;
                    while (1) ;
                }
            }
            /* DEBUG */

            pt_data *pt_local = (pt_data *) malloc(sizeof(pt_data));

            pt_local->halo = false;
            pt_local->geom = curLevelMatrix->geom;
            pt_local->size = curLevelMatrix->totalNumberOfRows;
            /* pt_ext and the rest of the information in pt_data is not needed. */

            Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR_temporary", partitioner_x, (void *)pt_local, LAIK_PF_None);

            // save pointer to old partitionings, since they need to be freed after call to re_switch_LaikVectors
            assert(curLevelMatrix->old_local == NULL);
            assert(curLevelMatrix->old_ext == NULL);
            curLevelMatrix->old_local = curLevelMatrix->local;

            curLevelMatrix->local = laik_new_partitioning(x_localPR, world, curLevelMatrix->space, NULL);

            curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
        }
    }

}

/**
 * @brief Switch to the new local partioning on every Laik_Blob 
 * 
 * @param A SparseMatrix
 * @param list of Laik_Blobs
 */
void re_switch_LaikVectors(SparseMatrix &A, std::vector<Laik_Blob *> list)
{
    for (uint32_t i = 0; i < list.size(); i++)
    {
        Laik_Blob * elem = list[i];
        elem->localLength = A.localNumberOfRows; /* Need to update local length, since it could have changed */

        laik_switchto_partitioning(elem->values, A.local, LAIK_DF_Preserve, LAIK_RO_None);
    }

    // MGData has two container with partitioning of next layer matrix and one of the current layer
    int numberOfMgLevels = 4; // Number of levels excluding last layer matrix (has no MGData)
    SparseMatrix *curLevelMatrix = &A;
    MGData * curMGData;
    for (int level = 1; level < numberOfMgLevels; ++level)
    {
        // std::cout << "\nCoarse Problem level " << level << std::endl;
        curMGData = curLevelMatrix->mgData;

        laik_switchto_partitioning(curMGData->Axf_blob->values, curLevelMatrix->local, LAIK_DF_Preserve, LAIK_RO_None);
        laik_switchto_partitioning(curMGData->rc_blob->values, curLevelMatrix->Ac->local, LAIK_DF_Preserve, LAIK_RO_None);
        laik_switchto_partitioning(curMGData->xc_blob->values, curLevelMatrix->Ac->local, LAIK_DF_Preserve, LAIK_RO_None);
        // Update local lengths as well
        curMGData->Axf_blob->localLength = curLevelMatrix->localNumberOfRows;
        curMGData->xc_blob->localLength = curLevelMatrix->Ac->localNumberOfRows;
        curMGData->rc_blob->localLength = curLevelMatrix->Ac->localNumberOfRows;

        curLevelMatrix = curLevelMatrix->Ac; // Next level
    }
}

#endif

// #######################################################################################
// Clean-up functions

/**
 * @brief Deallocate L2A_map Struct
 *
 * @param mapping to be deallocated
 */
void free_L2A_map(L2A_map * mapping)
{
    if(mapping != NULL)
    {
        mapping->localNumberOfRows = 0;
        mapping->offset = 0;
        mapping->offset_ext = 0;

        if(!mapping->localToExternalMap.empty())
            mapping->localToExternalMap.clear();
        if (!mapping->localToGlobalMap.empty())
            mapping->localToGlobalMap.clear();

        free((void*)mapping);

        mapping = NULL;
    }
}

/**
 * @brief Deallocate a Laik vector
 *
 * @param x to be deallocated
 */
void DeleteLaikVector(Laik_Blob *x)
{
    x->localLength = 0;
    x->exchange = 0;

    // TODO how to free Laik_data ?
}

// #######################################################################################
// Debug functions

/**
 * @brief Compare two double values x and y
 * 
 * @param x 
 * @param y 
 * @param doIO 
 * @param curIndex 
 */
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

/**
 * @brief Compare the two vectors x and y.
 * 
 * @param x 
 * @param y 
 * @param mapping 
 * @param doIO 
 */
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

/**
 * @brief Print the vector
 * 
 * @param x 
 */
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

/**
 * @brief Print the Laik Vector
 * 
 * @param x 
 * @param mapping 
 */
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

/**
 * @brief Print information about the SparseMatrix
 * 
 * @param spm 
 * @param coarseLevel 
 */
void printSPM(SparseMatrix *spm, int coarseLevel)
{
    
    // Global data
    HPCG_fout << "\n##################### Global stats #####################\n\n";

    HPCG_fout << "\nTotal # of rows " << spm->totalNumberOfRows
              << std::endl
              << "\nTotal # of nonzeros " << spm->totalNumberOfNonzeros
              << std::endl;


    HPCG_fout << "\n##################### Local stats #####################\n\n";

    // Local
    HPCG_fout << "\nLocal # of rows " << spm->localNumberOfRows
              << std::endl
              << "\nLocal # of nonzeros " << spm->localNumberOfNonzeros
              << std::endl;

    // HPCG_fout << "NumberOfExternalValues: " << (spm->localNumberOfColumns - spm->localNumberOfRows)
    //           << std::endl;

    // HPCG_fout << "\n##################### Mapping of rows #####################\n\n";

    // // Global to local mapping:
    // HPCG_fout << "\nLocal-to-global Map\n";
    // HPCG_fout << "Local\tGlobal\n";
    // for (int c = 0; c < spm->localToGlobalMap.size(); c++)
    // {
    //     HPCG_fout << c << "\t\t" << spm->localToGlobalMap[c] << std::endl;
    // }

    // HPCG_fout << "\nGlobal-to-local Map\n";
    // HPCG_fout << "Global\tlocal\n";
    // for (int c = 0; c < spm->globalToLocalMap.size(); c++)
    // {
    //   HPCG_fout << c << "\t\t" << spm->globalToLocalMap[c] << std::endl;
    // }

    // // Non zero indexes
    // HPCG_fout << "\n\n##################### Local subblock in matrix #####################\n\n";
    // for (uint64_t row_i = 0; row_i < spm->localNumberOfRows && row_i < 8; row_i++)
    // {
    //   HPCG_fout << "Row " << row_i << " (" << (int)spm->nonzerosInRow[row_i] << " non zeros) mtxIndL" << std::endl
    //             << std::endl;

    //   for (uint64_t nz_column_j = 0; nz_column_j < spm->nonzerosInRow[row_i]; nz_column_j++)
    //   {
    //     HPCG_fout << "Index (" << row_i << "," << spm->mtxIndL[row_i][nz_column_j] << ") = " << spm->matrixValues[row_i][nz_column_j] << std::endl;
    //   }
    //   HPCG_fout << std::endl;
    // }

    if(spm->geom->rank != 0)
    {
        // Global data
        std::cout << "\n##################### Global stats #####################\n\n";

        std::cout << "\nTotal # of rows " << spm->totalNumberOfRows
                  << std::endl
                  << "\nTotal # of nonzeros " << spm->totalNumberOfNonzeros
                  << std::endl;

        std::cout << "\n##################### Local stats #####################\n\n";

        // Local
        std::cout << "\nLocal # of rows " << spm->localNumberOfRows
                  << std::endl
                  << "\nLocal # of nonzeros " << spm->localNumberOfNonzeros
                  << std::endl;

        // std::cout << "NumberOfExternalValues: " << (spm->localNumberOfColumns - spm->localNumberOfRows)
        //           << std::endl;

        // std::cout << "\n##################### Mapping of rows #####################\n\n";

        // // Global to local mapping:
        // std::cout << "\nLocal-to-global Map\n";
        // std::cout << "Local\tGlobal\n";
        // for (int c = 0; c < spm->localToGlobalMap.size(); c++)
        // {
        //     std::cout << c << "\t\t" << spm->localToGlobalMap[c] << std::endl;
        // }
    }
}

/**
 * @brief Print members of HPCG_Params Struct
 * 
 * @param params 
 * @param doIO 
 */
void print_HPCG_PARAMS(HPCG_Params params, bool doIO)
{
    if(doIO)
    {
        std::string a{"####### PARAM values\n"};
        a += "npx: " + std::to_string(params.npx) + "\n";
        a += "npy: " + std::to_string(params.npy) + "\n";
        a += "npz: " + std::to_string(params.npz) + "\n";
        a += "numThreads: " + std::to_string(params.numThreads) + "\n";
        a += "nx: " + std::to_string(params.nx) + "\n";
        a += "ny: " + std::to_string(params.ny) + "\n";
        a += "nz: " + std::to_string(params.nz) + "\n";
        a += "pz: " + std::to_string(params.pz) + "\n";
        a += "runningTime: " + std::to_string(params.runningTime) + "\n";
        a += "zl: " + std::to_string(params.zl) + "\n";
        a += "zu: " + std::to_string(params.zu) + "\n";

        std::cout << a;
    }
    
    return ;
}

/**
 * @brief Print members of geometry struct
 * 
 * @param geom 
 * @param doIO 
 */
void print_GEOMETRY(Geometry * geom, bool doIO)
{
    if(doIO)
    {
        std::string a{"####### Geometry values\n"};

        a += "gix0: " + std::to_string(geom->gix0) + "\n";
        a += "giy0: " + std::to_string(geom->giy0) + "\n";
        a += "giz0: " + std::to_string(geom->giz0) + "\n";
        a += "gnx: " + std::to_string(geom->gnx) + "\n";
        a += "gny: " + std::to_string(geom->gny) + "\n";
        a += "gnz: " + std::to_string(geom->gnz) + "\n";
        a += "ipx: " + std::to_string(geom->ipx) + "\n";
        a += "ipy: " + std::to_string(geom->ipy) + "\n";
        a += "ipz: " + std::to_string(geom->ipz) + "\n";
        a += "npartz: " + std::to_string(geom->npartz) + "\n";
        a += "npx: " + std::to_string(geom->npx) + "\n";
        a += "npy: " + std::to_string(geom->npy) + "\n";
        a += "npz: " + std::to_string(geom->npz) + "\n";
        a += "numThreads: " + std::to_string(geom->numThreads) + "\n";
        a += "nx: " + std::to_string(geom->nx) + "\n";
        a += "ny: " + std::to_string(geom->ny) + "\n";
        a += "nz: " + std::to_string(geom->nz) + "\n";
        a += "pz: " + std::to_string(geom->pz) + "\n";
       
        a += "partz_ids: ";
        for (int i = 0; i < geom->npartz; ++i)
        {
            a += std::to_string(geom->partz_ids[i]);

            if(i != geom->npartz - 1)
                a += ", ";
        }
        a += "\n";

        a += "partz_nz: ";

        for (int i = 0; i < geom->npartz; ++i)
        {
            a += std::to_string(geom->partz_nz[i]);

            if (i != geom->npartz - 1)
                a += ", ";
        }
        a += "\n";

        a += "rank: " + std::to_string(geom->rank) + "\n";
        a += "size: " + std::to_string(geom->size) + "\n";



        std::cout << a;
    }
}

void print_LaikBlob(Laik_Blob * x)
{
    if(x != NULL)
    {
        std::string str{"####### Laik_Blob\n"};

        std::string str2{x->name};
        // str += "Exchange: " + std::to_string(x->exchange) + "\n";
        str += "LAIK " + std::to_string(laik_myid(world)) + "\n";
        str += "Vector name: " + str2 + "\n";
        str += "localLength: " + std::to_string(x->localLength) + "\n#######\n";

        std::cout << str;
    }
}

/**
 * @brief Debug function. 
 * 
 * Exit the program on all processes
 * 
 * @param msg will be printed (No need to add "\\n")
 */
void exit_hpcg_run(const char * msg)
{
    if(laik_myid(world) == 0)
        printf("\n####### %s\n####### Debug DONE -> Exiting #######\n", msg);

    exit(0);

    return ;
}