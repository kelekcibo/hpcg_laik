/**
 * @file laik_repartition.hpp
 * @brief Everything needed for enabling repartitioning
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

/*
    Includes
*/
#include <iostream>
#include <cstring>
#include <cassert>

#include "laik/hpcg_laik.hpp"
#include "hpcg.hpp"
#include "Geometry.hpp"
#include "GenerateProblem.hpp"
#include "SetupHalo.hpp"
#include "GenerateCoarseProblem.hpp"
#include "GenerateGeometry.hpp"
#include "CheckProblem.hpp"
/*
    Includes -END
*/

/*
    Partitioner algorithm for SparseMatrix
*/
#ifdef REPARTITION
/**
 * @brief Partitioner algorithm for 1d members of a SparseMatrix.
 *
 * The data which needs to be distributed/partitioned will be "splitted" here,
 *
 * For now, following data will be distributed:
 *  - nonzerosInRow (1d)
 *  - MGData.f2cOperator (1d)
 *  - matrixDiagonal (1d)
 *  - matrixValues (2d)
 *  - mtxIndG (2d)
 *
 */
void partitioner_1d_members_of_A(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    SparseMatrix * A = (SparseMatrix *)laik_partitioner_data(p->partitioner);
    Laik_Space * A_space = p->space;
    int rank = A->geom->rank;

    assert(A->totalNumberOfRows == laik_space_size(A_space));

    /* Calculate mapping. Needed to map local to allocation indices */
    uint64_t * mapping_ = A->mapping_;
    uint64_t localIndex = 0;

    Laik_Range range;
    for (long long i = 0; i < A->totalNumberOfRows; i++)
    {
        // assign every process its global part
        int proc = ComputeRankOfMatrixRow(*(A->geom), i);

        if(A->repartitioned)
            assert(proc == 0);

        laik_range_init_1d(&range, A_space, i, i + 1);
        laik_append_range(r, proc, &range, 1, 0);

        /* Calculating the offset into the allocation buffer due to lex_layout needs to be done once for all LAIK containers in A */
        if (A->offset_ == -1 && proc == rank)
            A->offset_ = i;

        if (A->offset_ != -1 && proc == rank)
            mapping_[localIndex++] = i - A->offset_; // GlobalIndex - Offset = Allocation Index;
    }
}

/**
 * @brief Partitioner algorithm for 2d members of a SparseMatrix. 
 *
 * 2d members of SparseMatrix A are realized by 1d array. 
 *
 * The data which needs to be distributed/partitioned will be "splitted" here,
 *
 * For now, following data will be distributed:
 *  - nonzerosInRow (1d)
 *  - MGData.f2cOperator (1d)
 *  - matrixDiagonal (1d)
 *  - matrixValues (2d)
 *  - mtxIndG (2d)
 *
 */
void partitioner_2d_members_of_A(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    SparseMatrix *A = (SparseMatrix *)laik_partitioner_data(p->partitioner);
    Laik_Space * space2d = p->space;

    Laik_Range range;
    for (long long i = 0; i < A->totalNumberOfRows; i++)
    {
        // assign every process its global part
        int proc = ComputeRankOfMatrixRow(*(A->geom), i);

        laik_range_init_1d(&range, space2d, i * numberOfNonzerosPerRow, i * numberOfNonzerosPerRow + numberOfNonzerosPerRow); // we are flattening the 2D array
        laik_append_range(r, proc, &range, 1, 0);
    }
}
/*
    Partitioner algorithm for SparseMatrix -END
*/

/*
    Needed functions/variables for shrink/expand feature
*/

/* Global copy of "HPCG_Params params" in main function */
HPCG_Params hpcg_params; 

/*
    Forw. Decl. of helper functions for repartition_SparseMatrix
*/
void update_Maps(SparseMatrix &A);
void update_Geometry(SparseMatrix &A);
void update_coarse_Geometry(SparseMatrix &Af);
void update_partitionings_x(SparseMatrix &A);
void update_partitionings_A(SparseMatrix &A);
void update_Local_Values(SparseMatrix &A);
void delete_mtxIndL(SparseMatrix &A);

/*
    Forw. Decl. of helper functions for repartition_SparseMatrix -END
*/

/**
 * @brief Repartition the SparseMatrix according to the new world size.
 *
 * Following has to be updated:
 *  - Geometry
 *  - Data which is distributed
 *  - Local values of A 
 *
 * This has to be done for all "matrix-layers"
 * 
 * @param A
 */
void repartition_SparseMatrix(SparseMatrix &A)
{



    int numberOfMgLevels = 4; // Number of levels including first
 
    // Delete old mtxIndL of the SparseMatrices (This is done here, because we need the old localNumberOfRows)
    SparseMatrix * curLevelMatrix = &A;
    // for (int level = 0; level < numberOfMgLevels; ++level)
    // {
    //     delete_mtxIndL(A); // double free fix that ? or garbage collecting
    //     curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    // }

    // Procs which will be removed need the old localNumberOfRows again as when we free old ressources
    local_int_t old_localNumberOfRows[numberOfMgLevels];
    if(laik_myid(world) < 0)
    {
        for (int level = 0; level < numberOfMgLevels; ++level)
        {
            old_localNumberOfRows[level] = curLevelMatrix->localNumberOfRows;
            curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
        }
    }


    /* Debug */
//     printf("DEBUGING START\n");
//     if (A.totalNumberOfRows == 8192 && A.geom->rank != 0)
//     {
//         std::string debug{"LAIK "};

//         debug += std::to_string(A.geom->rank) + "\tPrinting f2c values\n";
//         local_int_t *f2c;
//         laik_get_map_1d(A.mgData->f2cOperator_d, 0, (void **)&f2c, 0);
//         uint64_t count = 0;

//         for (size_t currentCoarseRow = 0; currentCoarseRow < 512; currentCoarseRow++)
//         {
//             count++;
//             // debug += "f2c[" + std::to_string(map_l2a_A(A, currentCoarseRow)) + "] = " + std::to_string(f2c[map_l2a_A(A, currentCoarseRow)]) + "\n";
//             // printf("f2c[%ld (AI:%lld) (GI: %lld)] = %d\n", currentCoarseRow, map_l2a_A(A, currentCoarseRow), map_l2a_A(A, currentCoarseRow) + A.offset_, f2c[map_l2a_A(A, currentCoarseRow)]);
//             printf("f2c[%ld (AI:%lld) (GI: %lld)] = %d\n", currentCoarseRow, map_l2a_A(A, currentCoarseRow), map_l2a_A(A, currentCoarseRow) + A.offset_, f2c[map_l2a_A(A, currentCoarseRow)]);

//             // printf("AI:%lld\n",map_l2a_A(A, currentCoarseRow));
//         }

//         debug += "\n I was " + std::to_string(count) + " times in the loop\n";

//         debug += "DONE DEBUGGING #########\n";
//         std::cout << debug;
//     }

//     /* Debug */


// while (1)
//     ;


update_Geometry(A);
// Update update_Geometry of the SparseMatrices
curLevelMatrix = &A;
for (int level = 1; level < numberOfMgLevels; ++level)
{
    update_coarse_Geometry(*curLevelMatrix);
    curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    }

    // Update maps of the SparseMatrices
    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
        update_Maps(*curLevelMatrix);
        curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    }

    // Update partitionings of the SparseMatrices
    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
        update_partitionings_A(*curLevelMatrix);
        curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    }
    
    // Update partitionings of the x vector
    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
        update_partitionings_x(*curLevelMatrix);
        curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    }

    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
        update_Local_Values(*curLevelMatrix);
        curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    }


    if(A.geom->rank < 0)
    {
        // Procs which will be removed need the old localNumberOfRows again as when we free old ressources
    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
        curLevelMatrix->localNumberOfRows = old_localNumberOfRows[level];
        curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
        }
    }

    return ;
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
        Laik_Blob *elem = list[i];
        elem->localLength = A.localNumberOfRows; /* Need to update local length, since it could have changed */

        laik_switchto_partitioning(elem->values, A.local, LAIK_DF_Preserve, LAIK_RO_None);
    }

    // MGData has two container with partitioning of next layer matrix and one of the current layer
    int numberOfMgLevels = 4; // Number of levels excluding last layer matrix (has no MGData)
    SparseMatrix *curLevelMatrix = &A;
    MGData *curMGData;
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

/**
 * @brief Initialize partitionings needed to partition the SparseMatrix.
 *
 * This partitionings are needed to partition the matrix and to be able to exchange values via LAIK
 *
 * @param A SparseMatrix
 */
void init_SPM_partitionings(SparseMatrix &A)
{
    // Init spaces
    if (A.space)
        exit_hpcg_run("Something's wrong. Value should be NULL.");
    if (A.space2d)
        exit_hpcg_run("Something's wrong. Value should be NULL.");

    A.mapping_ = new uint64_t[A.localNumberOfRows];
    A.offset_ = -1;

    
    A.space = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows);
    A.space2d = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows * numberOfNonzerosPerRow);

    A.nonzerosInRow_d = laik_new_data(A.space, laik_UChar);
    A.matrixDiagonal_d = laik_new_data(A.space, laik_Double);
    A.mtxIndG_d = laik_new_data(A.space2d, laik_UInt64);
    A.matrixValues_d = laik_new_data(A.space2d, laik_Double);

    Laik_Partitioner *partitioner_1d = laik_new_partitioner("partitioner_1d_members_of_A", partitioner_1d_members_of_A, (void *)&A, LAIK_PF_None);
    Laik_Partitioner *partitioner_2d = laik_new_partitioner("partitioner_2d_members_of_A", partitioner_2d_members_of_A, (void *)&A, LAIK_PF_None);

    A.partitioning_1d = laik_new_partitioning(partitioner_1d, world, A.space, NULL);
    A.partitioning_2d = laik_new_partitioning(partitioner_2d, world, A.space2d, NULL);

    // Split data
    laik_switchto_partitioning(A.nonzerosInRow_d, A.partitioning_1d, LAIK_DF_None, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixDiagonal_d, A.partitioning_1d, LAIK_DF_None, LAIK_RO_None);
    laik_switchto_partitioning(A.mtxIndG_d, A.partitioning_2d, LAIK_DF_None, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixValues_d, A.partitioning_2d, LAIK_DF_None, LAIK_RO_None);
}

/**
 * @brief Map the local index (mtxIndL) to the corresponding global index. Then map this global index to the allocation index.
 *
 * This extra step has to be done since Laik makes use of the lex_layout under the hood.
 * v1.2 will hopefully implement the custom hpcg_layout, so that there is no need for this anymore
 *
 * @param[in] SparseMatrix contains offset_ used to calculate the allocation index
 * @param[in] localIndex which will be mapped to the corresponding allocation index
 *
 * @return index to the allocation buffer
 *
 * @see L2A_map
 */
allocation_int_t map_l2a_A(const SparseMatrix &A, local_int_t localIndex)
{
    if(localIndex < 0 || localIndex >= A.localNumberOfRows)
        exit_hpcg_run("LocalIndex not mapped to Allocation Index; Fatal Error.");

    return A.mapping_[localIndex];
}
 
 
/**
 * @brief Replace specified diagonal values of the matrix after calling "ReplaceMatrixDiagonal(A, exaggeratedDiagA)".
 * 
 * This is needed in the case of enabling repartition, because the Laik data container does not store the pointer to matrixValues.
 * So, we need to update the diagonal values in matrixValues as well.
 * 
 *  @param[inout] A The system matrix.
 */
void replaceMatrixValues(SparseMatrix &A)
{
    // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
    // below may result in global range values.
    global_int_t nx = A.geom->nx;
    global_int_t ny = A.geom->ny;
    global_int_t nz = A.geom->nz;
    global_int_t gnx = A.geom->gnx;
    global_int_t gny = A.geom->gny;
    global_int_t gnz = A.geom->gnz;
    global_int_t gix0 = A.geom->gix0;
    global_int_t giy0 = A.geom->giy0;
    global_int_t giz0 = A.geom->giz0;

    double *matrixDiagonal;
    laik_get_map_1d(A.matrixDiagonal_d, 0, (void **)&matrixDiagonal, 0);

    double *matrixValues;
    laik_get_map_1d(A.matrixValues_d, 0, (void **)&matrixValues, 0);

    for (local_int_t iz = 0; iz < nz; iz++)
    {
        global_int_t giz = giz0 + iz;
        for (local_int_t iy = 0; iy < ny; iy++)
        {
            global_int_t giy = giy0 + iy;
            for (local_int_t ix = 0; ix < nx; ix++)
            {
                global_int_t gix = gix0 + ix;
                local_int_t currentLocalRow = iz * nx * ny + iy * nx + ix;
                global_int_t currentGlobalRow = giz * gnx * gny + giy * gnx + gix;
                uint64_t currentValuePointer_index = -1;  // Index to current value in current row
                for (int sz = -1; sz <= 1; sz++)
                {
                    if (giz + sz > -1 && giz + sz < gnz)
                    {
                        for (int sy = -1; sy <= 1; sy++)
                        {
                            if (giy + sy > -1 && giy + sy < gny)
                            {
                                for (int sx = -1; sx <= 1; sx++)
                                {
                                    if (gix + sx > -1 && gix + sx < gnx)
                                    {
                                        global_int_t curcol = currentGlobalRow + sz * gnx * gny + sy * gnx + sx;
                                        currentValuePointer_index++;
                                        
                                        if (curcol == currentGlobalRow)
                                        {
                                            matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + currentValuePointer_index] = matrixDiagonal[map_l2a_A(A, currentLocalRow)]; 
                                        }
                                       
                                    } // end x bounds test
                                }     // end sx loop
                            }         // end y bounds test
                        }             // end sy loop
                    }                 // end z bounds test
                }                     // end sz loop
            } // end ix loop
        }     // end iy loop
    }         // end iz loop
    return;
}


/*
    Needed functions/variables for shrink/expand feature -END
*/

/*
    Helper functions for shrink/expand feature
        Forward declare them in the beginning of the file
*/

/**
 * @brief localToGlobalMap and globalToLocalMap needs to be updated after repartitioning
 *
 * This code is from GenerateProblem_ref.cpp.
 *
 * @param A
 */
void update_Maps(SparseMatrix &A)
{

    // Removed procs do not need this upate
    if(A.geom->rank < 0)
        return;
    
    A.localToGlobalMap.clear();
    A.globalToLocalMap.clear();

    // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
    // below may result in global range values.
    global_int_t nx = A.geom->nx;
    global_int_t ny = A.geom->ny;
    global_int_t nz = A.geom->nz;
    global_int_t gnx = A.geom->gnx;
    global_int_t gny = A.geom->gny;
    // global_int_t gnz = A.geom->gnz;
    global_int_t gix0 = A.geom->gix0;
    global_int_t giy0 = A.geom->giy0;
    global_int_t giz0 = A.geom->giz0;

    for (local_int_t iz = 0; iz < nz; iz++)
    {
        global_int_t giz = giz0 + iz;
        for (local_int_t iy = 0; iy < ny; iy++)
        {
            global_int_t giy = giy0 + iy;
            for (local_int_t ix = 0; ix < nx; ix++)
            {
                global_int_t gix = gix0 + ix;
                local_int_t currentLocalRow = iz * nx * ny + iy * nx + ix;
                global_int_t currentGlobalRow = giz * gnx * gny + giy * gnx + gix;
#ifndef HPCG_NO_OPENMP
                // C++ std::map is not threadsafe for writing
#pragma omp critical
#endif
                A.globalToLocalMap[currentGlobalRow] = currentLocalRow;

                A.localToGlobalMap[currentLocalRow] = currentGlobalRow;
            }
        }
    }

    // or just call ComputeRAnk... ?
    /* DEBUG */
    /* if all rows are mapped on proc 0 */
    for (int i = 0; i < A.totalNumberOfRows; i++)
    {
        int rank = ComputeRankOfMatrixRow(*(A.geom), i);
        // TEST CASE. -np 2 => shrinking to world size 1
        if (rank != 0)
        {
            std::string str{"Rank not equal to 0: " + std::to_string(rank)};
            std::cout << str;
            while (1)
                ;
        }
    }
    /* DEBUG */

    return;
}

/**
 * @brief HPCG_Params and Geometry need to be updated after repartitioning
 *
 * New joining processes also need to get the values of hpcg_params, if broadcast was done in HPCG_init.
 *
 * @param[in] A SparseMatrix
 */
void update_Geometry(SparseMatrix &A)
{
    /*
        Old local problem size is needed to calculate new local problem size.
        Because, we want to have the same total problem size.
    */
    global_int_t old_gnx = A.geom->gnx;
    global_int_t old_gny = A.geom->gny;
    global_int_t old_gnz = A.geom->gnz;

    // Delete old setup.
    DeleteGeometry(*(A.geom));

    // TODO. Send params_struct to new procs; old procs already have them
    // Broadcast updated values. New joining processes will need them, if broadcast was done in init
    // laik_broadcast(iparams, iparams, nparams, laik_Int32);
    
    // Update parameters
    hpcg_params.comm_rank = laik_myid(world);
    hpcg_params.comm_size = laik_size(world);

    // Construct the new geometry and linear system
    Geometry *new_geom = new Geometry;

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

    // Additionally, update localNumberOfRows
    A.localNumberOfRows = A.geom->nx * A.geom->ny * A.geom->nz;

    return;
}

/**
 * @brief Geometry of coarse need to be updated as well
 *
 * @param[in] Af Current Level Matrix Layer 
 */
void update_coarse_Geometry(SparseMatrix &Af)
{
    // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
    // below may result in global range values.
    global_int_t nxf = Af.geom->nx;
    global_int_t nyf = Af.geom->ny;
    global_int_t nzf = Af.geom->nz;

    local_int_t nxc, nyc, nzc; // Coarse nx, ny, nz
    assert(nxf % 2 == 0); assert(nyf % 2 == 0); assert(nzf % 2 == 0); // Need fine grid dimensions to be divisible by 2
    nxc = nxf / 2; nyc = nyf / 2; nzc = nzf / 2;

    // Construct the geometry and linear system
    Geometry *geomc = new Geometry;

    local_int_t zlc = 0; // Coarsen nz for the lower block in the z processor dimension
    local_int_t zuc = 0; // Coarsen nz for the upper block in the z processor dimension
    int pz = Af.geom->pz;
    if (pz > 0)
    {
        zlc = Af.geom->partz_nz[0] / 2; // Coarsen nz for the lower block in the z processor dimension
        zuc = Af.geom->partz_nz[1] / 2; // Coarsen nz for the upper block in the z processor dimension
    }

    // TEST AFTERWARDS
    global_int_t old_gnx = Af.Ac->geom->gnx;
    global_int_t old_gny = Af.Ac->geom->gny;
    global_int_t old_gnz = Af.Ac->geom->gnz;

    DeleteGeometry(*Af.Ac->geom);
    GenerateGeometry(Af.geom->size, Af.geom->rank, Af.geom->numThreads, Af.geom->pz, zlc, zuc, nxc, nyc, nzc, Af.geom->npx, Af.geom->npy, Af.geom->npz, geomc);
    Af.Ac->geom = geomc;
    Af.Ac->localNumberOfRows = nxc * nyc * nzc;

    // should be the same problem size
    assert(old_gnx == geomc->gnx);
    assert(old_gny == geomc->gny);
    assert(old_gnz == geomc->gnz);
    return;
}

/**
 * @brief Need to delete old mtxIndL as it will be re-init'ed in SetupHalo(A).
 *
 * We need the old size of localNumberOfRows,thus, we call this function before every other function in repartition_SparseMatrix(A)
 *
 * @param A SparseMatrix
 */
void delete_mtxIndL(SparseMatrix &A)
{
#ifndef HPCG_CONTIGUOUS_ARRAYS
    for (local_int_t i = 0; i < A.localNumberOfRows; ++i)
        delete[] A.mtxIndL[i];

    // if (A.mtxIndL) delete[] A.mtxIndL;

#else
    delete[] A.mtxIndL[0];
#endif
    return;
}

/**
 * @brief Helper function for  update_partitionings_x(SparseMatrix &A)
 *
 * @param A
 */
void re_init_mtxIndL(SparseMatrix &A)
{

    local_int_t **mtxIndL = new local_int_t *[A.localNumberOfRows];

    for (local_int_t i = 0; i < A.localNumberOfRows; ++i)
        mtxIndL[i] = 0;

#ifndef HPCG_CONTIGUOUS_ARRAYS
    // Now allocate the arrays pointed to
    for (local_int_t i = 0; i < A.localNumberOfRows; ++i)
      mtxIndL[i] = new local_int_t[numberOfNonzerosPerRow];
#else
    // Now allocate the arrays pointed to
    mtxIndL[0] = new local_int_t[localNumberOfRows * numberOfNonzerosPerRow];
    for (local_int_t i = 1; i < localNumberOfRows; ++i)
        mtxIndL[i] = mtxIndL[0] + i * numberOfNonzerosPerRow;
#endif

    A.mtxIndL= mtxIndL;
    
    return;
}

/**
 * @brief Partitionings for Laik Vectors need to be updated after repartitioning
 *
 * @param A SparseMatrix
 */
void update_partitionings_x(SparseMatrix &A)
{
    // removed procs only need "local" partitioning
    if(A.geom->rank < 0)
    {
        partition_d *pt_local = (partition_d *)malloc(sizeof(partition_d));

        pt_local->halo = false;
        pt_local->geom = A.geom;
        pt_local->size = A.totalNumberOfRows;
        pt_local->offset = 0; // don't need this
        /* pt_ext and the rest of the information in pt_data is not needed. */

        Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR_temporary", partitioner_alg_for_x_vector, (void *)pt_local, LAIK_PF_None);

        A.ext = NULL; // Partitioning is deleted later
        A.local = laik_new_partitioning(x_localPR, world, A.space, NULL);

      
        return;
    }

    // Existing procs will call setuphalo, as everything is calculcated there
    // As mtxIndL is recalculated:
    re_init_mtxIndL(A);

    // clear old variables
    free_L2A_map(A.mapping);
    // TODO other vars


    SetupHalo(A);
 
    return;
}

/**
 * @brief Partitionings for the SparseMatrix need to be updated after repartitioning.
 * 
 * Old ressources will be freed in the end
 *
 * @param A SparseMatrix
 */
void update_partitionings_A(SparseMatrix &A)
{

    // Removed procs do not need to recalculate this
    if (A.geom->rank >= 0)
    {
        if (A.mapping_) delete[] A.mapping_;
        A.mapping_ = new uint64_t[A.localNumberOfRows];
        A.offset_ = -1;
    }
    A.repartitioned = true;

    Laik_Partitioner * partitioner_1d = laik_new_partitioner("partitioner_1d_members_of_A", partitioner_1d_members_of_A, (void *)&A, LAIK_PF_None);
    Laik_Partitioner * partitioner_2d = laik_new_partitioner("partitioner_2d_members_of_A", partitioner_2d_members_of_A, (void *)&A, LAIK_PF_None);

    Laik_Partitioning * new_partitioning_1d = laik_new_partitioning(partitioner_1d, world, A.space, NULL);
    Laik_Partitioning * new_partitioning_2d = laik_new_partitioning(partitioner_2d, world, A.space2d, NULL);

    /* Debug */
    if (A.geom->rank >= 0)
    {
        for (size_t i = 0; i < A.localNumberOfRows; i++)
        {
            assert(map_l2a_A(A, i) <= A.totalNumberOfRows);
        }
        
    }

    // Split data
    laik_switchto_partitioning(A.nonzerosInRow_d, new_partitioning_1d, LAIK_DF_Preserve, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixDiagonal_d, new_partitioning_1d, LAIK_DF_Preserve, LAIK_RO_None);
    laik_switchto_partitioning(A.mtxIndG_d, new_partitioning_2d, LAIK_DF_Preserve, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixValues_d, new_partitioning_2d, LAIK_DF_Preserve, LAIK_RO_None);

    // Last coarse matrix does not have mgData
    if(A.mgData != NULL)
    {

        laik_switchto_partitioning(A.mgData->f2cOperator_d, new_partitioning_1d, LAIK_DF_Preserve, LAIK_RO_None);

        // /* Debug */
        // if (A.totalNumberOfRows == 8192 && A.geom->rank == 0)
        // {
        //     printf("DEBUGING START\n");
        //     std::string debug{"LAIK "};

        //     debug += std::to_string(A.geom->rank) + "\tPrinting f2c values\n";
        //     local_int_t *f2c;
        //     laik_get_map_1d(A.mgData->f2cOperator_d, 0, (void **)&f2c, 0);
        //     uint64_t count = 0;

        //     for (size_t currentCoarseRow = 0; currentCoarseRow < 1024; currentCoarseRow++)
        //     {
        //         count++;
        //         // debug += "f2c[" + std::to_string(map_l2a_A(A, currentCoarseRow)) + "] = " + std::to_string(f2c[map_l2a_A(A, currentCoarseRow)]) + "\n";
        //         // printf("f2c[%ld (AI:%lld) (GI: %lld)] = %d\n", currentCoarseRow, map_l2a_A(A, currentCoarseRow), map_l2a_A(A, currentCoarseRow) + A.offset_, f2c[map_l2a_A(A, currentCoarseRow)]);
        //         printf("f2c[%ld (AI:%lld) (GI: %lld)] = %d\n", currentCoarseRow, map_l2a_A(A, currentCoarseRow), map_l2a_A(A, currentCoarseRow) + A.offset_, f2c[map_l2a_A(A, currentCoarseRow)]);

        //         // printf("AI:%lld\n",map_l2a_A(A, currentCoarseRow));
        //     }

        //     debug += "\n I was " + std::to_string(count) + " times in the loop\n";

        //     debug += "DONE DEBUGGING #########\n";
        //     std::cout << debug;
        // }

        // /* Debug */
    }

    // while (1)
    //     ;

    laik_free_partitioning(A.partitioning_1d);
    laik_free_partitioning(A.partitioning_2d);

    A.partitioning_1d = new_partitioning_1d;
    A.partitioning_2d = new_partitioning_2d; 
}

/**
 * @brief Local values of the SparseMatrix need to be updated after repartitioning
 *
 * @param A SparseMatrix
 */
void update_Local_Values(SparseMatrix &A)
{
    // Removed procs do not need this upate
    if (A.geom->rank < 0)
        return;
        
    char *nonzerosInRow;
    laik_get_map_1d(A.nonzerosInRow_d, 0, (void **)&nonzerosInRow, 0);

    local_int_t localNumberOfNonZeros = 0;
    for (local_int_t i = 0; i < A.localNumberOfRows; i++)
        localNumberOfNonZeros += nonzerosInRow[map_l2a_A(A, i)];

    A.localNumberOfNonzeros = localNumberOfNonZeros;

    return;
}

/*
    Helper functions for shrink/expand feature -END
        Forward declare them in the beginning of the file
*/
#endif