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

    if (hpcg_params.comm_rank >= 0)
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
                    while (1)
                        ;
                }
            }
            /* DEBUG */

            pt_data *pt_local = (pt_data *)malloc(sizeof(pt_data));

            pt_local->halo = false;
            pt_local->geom = curLevelMatrix->geom;
            pt_local->size = curLevelMatrix->totalNumberOfRows;
            /* pt_ext and the rest of the information in pt_data is not needed. */

            Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR_temporary", partitioner_alg_for_x_vector, (void *)pt_local, LAIK_PF_None);

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
