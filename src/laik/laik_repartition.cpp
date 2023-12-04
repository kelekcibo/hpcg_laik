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
#include "ComputeOptimalShapeXYZ.hpp"

/*
    Includes -END
*/

/*
    Partitioner algorithm for SparseMatrix
*/
#ifdef REPARTITION

/**
 * @brief This partitioner is a little optimisation: 
 * 
 * if new processes join the computation, initial processes send several important local variables computed during the CG to old processes.
 * Therefore, it is enough, if only one initial process sends its values ONLY to newly joined processes.
 *
*/
void new_joining_procs(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    Laik_Space * space = p->space;
    int n = laik_space_size(space);

    int old_size = *((int *) laik_partitioner_data(p->partitioner));
    int new_size = laik_size(p->group);

    Laik_Range range;
    for (long long proc = old_size; proc < new_size; proc++)
    {
        laik_range_init_1d(&range, space, 0, n);
        laik_append_range(r, proc, &range, 1, 0);
    }

    // Partial broadcast is done by one initial process. So the initial process here is the process with ID 0.
    laik_range_init_1d(&range, space, 0, n);
    laik_append_range(r, 0, &range, 1, 0); // Processes with ID 0 will need to access it as well

    return;
}

/**
 * @brief Partitioner algorithm for 1d members of a SparseMatrix.
 *
 * The data which needs to be distributed/partitioned will be "splitted" here,
 *
 * For now, following data will be distributed:
 *  - nonzerosInRow (1d)
 *  - matrixDiagonal (1d)
 *  - matrixValues (2d)
 *  - mtxIndG (2d)
 *
 */
void partitioner_1d_members_of_A(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{
    SparseMatrix * A = (SparseMatrix *)laik_partitioner_data(p->partitioner);
    Laik_Space * A_space = p->space;

    assert(A->totalNumberOfRows == laik_space_size(A_space));

    Laik_Range range;
    for (long long i = 0; i < A->totalNumberOfRows; i++)
    {
        // assign every process its global part
        int proc = ComputeRankOfMatrixRow(*(A->geom), i);

        laik_range_init_1d(&range, A_space, i, i + 1);
        laik_append_range(r, proc, &range, 1, 0);
    }
    return;
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
    return;
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
void delete_mtxIndL(SparseMatrix &A);
void update_Geometry(SparseMatrix &A);
void update_coarse_Geometry(SparseMatrix &Af);
void update_partitionings_A(SparseMatrix &A);
void update_partitionings_x(SparseMatrix &A);
void broadcast_hpcg_params(void);
Geometry * calculate_old_geometry(SparseMatrix &A, int old_size);
Geometry * calculate_old_coarse_geometry(Geometry * old_geometry, int old_size);

/*
    Forw. Decl. of helper functions for repartition_SparseMatrix -END
*/

/**
 * @brief Redistribute data of the SparseMatrix according to the new world size.
 *
 * Following has to be updated:
 *  - Geometry
 *  - Data which is distributed
 *  - Local values of A 
 *
 * This has to be done for all "matrix-layers"
 * 
 * This function is only called by initial computing processes.
 * 
 * @param[inout] A the sparsematrix to be repartioned
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

    // Update Geometry of the SparseMatrices
    update_Geometry(A);
    curLevelMatrix = &A;
    for (int level = 1; level < numberOfMgLevels; ++level)
    {
        update_coarse_Geometry(*curLevelMatrix);
        curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
    }

    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level)
    {
        // Update partitionings for the SparseMatrices
        update_partitionings_A(*curLevelMatrix);
        // Update local/global variables
        update_Values(*curLevelMatrix);
        // Update partitionings for the LAIK vectors
        update_partitionings_x(*curLevelMatrix);
        // Make the nextcoarse grid the next level
        curLevelMatrix = curLevelMatrix->Ac;
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
 * @brief Redistribute all LAIK Vectors according to the new world.
 *
 * @param[inout] A SparseMatrix
 * @param[in] list of LAIK Vectors to be redistributed
 */
void re_switch_LaikVectors(SparseMatrix &A, std::vector<Laik_Blob *> list)
{
    int numberOfMgLevels = 4; // Number of levels including first

    // #################################################################################
    // ######### New joining procs have no initial partitiong activated on Laik Data containers.
    // ######### Need to calculate old partitioning
    // ######### Thus, we need the old geometry to create the correct local partitioning of the previous world
    Geometry * old_geometry[numberOfMgLevels];
    Laik_Partitioning *local_old[numberOfMgLevels];
    Laik_Group *parent = laik_group_parent(world);
    int old_size = laik_size(parent);

    if (laik_myid(world) >= old_size)
    {
        // Initialize pointers to old local partitionings
        for (int i = 0; i < numberOfMgLevels; i++) local_old[i] = 0;
        // calculate old geometry for layer 0 (main level matrix)
        old_geometry[0] = calculate_old_geometry(A, old_size);
        partition_d pt_local;
        pt_local.halo = false; pt_local.geom = old_geometry[0]; pt_local.size = A.totalNumberOfRows;
        pt_local.offset = 0; // don't need this, as well as pt_ext and the rest of the information in pt_data 
        // Create partitioner
        Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR", partitioner_alg_for_x_vector, (void *)&pt_local, LAIK_PF_None);
        // partitionings for the other layers are created later
        local_old[0] = laik_new_partitioning(x_localPR, parent, A.space, NULL);
    }
    // ######### Preparation of necessary information of new joining processes done
    // #################################################################################


    // Redistribute data of laik vectors in list
    for (uint32_t i = 0; i < list.size(); i++)
    {
        Laik_Blob *elem = list[i];
        elem->localLength = A.localNumberOfRows; // Need to update local length, since it probably changed

        // create new layout data for all LAIK Vectors
        Laik_vector_data *layout_data = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));
        layout_data->localLength = A.localNumberOfRows;
        layout_data->numberOfExternalValues = elem->exchangesValues ? A.numberOfExternalValues : 0;
        layout_data->id = laik_myid(world);
        laik_data_set_layout_data(elem->values, layout_data);

        // New joining procs have no initial partitiong activated on Laik Data containers.
        if (laik_myid(world) >= old_size)
            laik_set_initial_partitioning(elem->values, local_old[0]); // local_old[0] previously calculated

        // if(laik_myid(world) != -1) /// NOOOO WRONG! WILL HAVE TO DO ONE COPY, => IMPLEMENT REUSE PROPERLY.
        //     laik_switchto_partitioning(elem->values, A.ext, LAIK_DF_None, LAIK_RO_None); // due to optimisation reasons as mentioned in init_blob()

        laik_switchto_partitioning(elem->values, A.local, LAIK_DF_Preserve, LAIK_RO_None);
    }

    // Redistribute data of laik vectors stored in MGData of all coarse grid layers of A
    SparseMatrix *curLevelMatrix = &A;
    MGData *curMGData;
    for (int level = 1; level < numberOfMgLevels; ++level)
    {
        curMGData = curLevelMatrix->mgData;
        
        // create new layout_data
        Laik_vector_data *layout_data_Axf = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));
        Laik_vector_data *layout_data_rc = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));
        Laik_vector_data *layout_data_xc = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));

        // prepare data for Axf
        layout_data_Axf->id = laik_myid(world);
        layout_data_Axf->localLength = curLevelMatrix->localNumberOfRows;
        layout_data_Axf->numberOfExternalValues = curMGData->Axf_blob->exchangesValues ? curLevelMatrix->numberOfExternalValues : 0; // will return curLevelMatrix->numberOfExternalValues
        // prepare data for rc
        layout_data_rc->id = laik_myid(world);
        layout_data_rc->localLength = curLevelMatrix->Ac->localNumberOfRows;
        layout_data_rc->numberOfExternalValues = curMGData->rc_blob->exchangesValues ? curLevelMatrix->numberOfExternalValues : 0; // will return 0
        // prepare data for xc
        layout_data_xc->id = laik_myid(world);
        layout_data_xc->localLength = curLevelMatrix->Ac->localNumberOfRows;
        layout_data_xc->numberOfExternalValues = curMGData->xc_blob->exchangesValues ? curLevelMatrix->numberOfExternalValues : 0; // will return curLevelMatrix->numberOfExternalValues
        // debug
        // laik_data_set_name(curMGData->Axf_blob->values, (char *)curMGData->Axf_blob->name);
        // laik_data_set_name(curMGData->rc_blob->values, (char *)curMGData->rc_blob->name);
        // laik_data_set_name(curMGData->xc_blob->values, (char *)"MG_Data xc71");

        // set new layout_data
        laik_data_set_layout_data(curMGData->Axf_blob->values, layout_data_Axf);
        laik_data_set_layout_data(curMGData->rc_blob->values, layout_data_rc);
        laik_data_set_layout_data(curMGData->xc_blob->values, layout_data_xc);
     
        // we previously calculated local_old[0], now calculate the other three old local partitionings
        if (laik_myid(world) >= old_size)
        {
            // calculate old geometry for layer <level> (coarse grid <level> of matrix A)
            old_geometry[level] = calculate_old_coarse_geometry(old_geometry[level - 1], old_size);
            partition_d pt_local;
            pt_local.halo = false; pt_local.geom = old_geometry[level]; pt_local.size = curLevelMatrix->Ac->totalNumberOfRows;
            pt_local.offset = 0; // don't need this, as well as pt_ext and the rest of the information in pt_data
            // Create partitioner
            Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR", partitioner_alg_for_x_vector, (void *)&pt_local, LAIK_PF_None);
            // partitionings for the other layers are created here
            local_old[level] = laik_new_partitioning(x_localPR, parent, curLevelMatrix->Ac->space, NULL);
            // New joining procs have no initial partitiong activated on Laik Data containers.
            laik_set_initial_partitioning(curMGData->Axf_blob->values, local_old[level - 1]);
            laik_set_initial_partitioning(curMGData->rc_blob->values, local_old[level]);
            laik_set_initial_partitioning(curMGData->xc_blob->values, local_old[level]);
        }

        laik_switchto_partitioning(curMGData->Axf_blob->values, curLevelMatrix->local, LAIK_DF_Preserve, LAIK_RO_None);
        laik_switchto_partitioning(curMGData->rc_blob->values, curLevelMatrix->Ac->local, LAIK_DF_Preserve, LAIK_RO_None);
        laik_switchto_partitioning(curMGData->xc_blob->values, curLevelMatrix->Ac->local, LAIK_DF_Preserve, LAIK_RO_None);
    
        // Update local lengths as well
        curMGData->Axf_blob->localLength = curLevelMatrix->localNumberOfRows;
        curMGData->xc_blob->localLength = curLevelMatrix->Ac->localNumberOfRows;
        curMGData->rc_blob->localLength = curLevelMatrix->Ac->localNumberOfRows;
        
        curLevelMatrix = curLevelMatrix->Ac; // Next level
    }

    if (laik_myid(world) >= old_size)
    {
        // Ddelete old geometries and free old partitionings
        for (int level = 0; level < numberOfMgLevels; level++)
        {
            DeleteGeometry(*(old_geometry[level]));
            laik_free_partitioning(local_old[level]);
        }
    }

    return;    
}


/**
 * @brief Initialize partitionings needed to partition the SparseMatrix.
 *
 * This partitionings are needed to partition the matrix and to be able to exchange values via LAIK
 *
 * @param[inout] A SparseMatrix
 */
void init_SPM_partitionings(SparseMatrix &A)
{
    if (A.space) exit_hpcg_run("Something's wrong. Value should be NULL.", false);
    if (A.space2d) exit_hpcg_run("Something's wrong. Value should be NULL.", false);
    // Init spaces
    A.space = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows);
    A.space2d = laik_new_space_1d(hpcg_instance, A.totalNumberOfRows * numberOfNonzerosPerRow);
    // Init data containers
    A.nonzerosInRow_d = laik_new_data(A.space, laik_UChar);
    A.matrixDiagonal_d = laik_new_data(A.space, laik_Double);
    A.mtxIndG_d = laik_new_data(A.space2d, laik_UInt64);
    A.matrixValues_d = laik_new_data(A.space2d, laik_Double);
    // use custom layouts instead of the default lex layout
    laik_data_set_layout_factory(A.nonzerosInRow_d, laik_new_layout_vector);
    laik_data_set_layout_flag(A.nonzerosInRow_d, LAIK_Vector_Layout);
    laik_data_set_layout_factory(A.matrixDiagonal_d, laik_new_layout_vector);
    laik_data_set_layout_flag(A.matrixDiagonal_d, LAIK_Vector_Layout);

    Laik_vector_data *layout_data = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));
    layout_data->localLength = A.localNumberOfRows;
    layout_data->numberOfExternalValues = 0; // not used 
    layout_data->id = laik_myid(world);
    laik_data_set_layout_data(A.nonzerosInRow_d, layout_data);
    laik_data_set_layout_data(A.matrixDiagonal_d, layout_data); // same metadata

    // laik_data_set_layout_factory(A.mtxIndG_d, laik_new_layout_sparse);
    // laik_data_set_layout_factory(A.matrixValues_d, laik_new_layout_sparse);

    // New joining procs
    // Need geometry with old config. Store current config
    Geometry * geometry = A.geom;
    Laik_Partitioning *partitioning_1d_old = 0, *partitioning_2d_old = 0;
    Laik_Group *parent = laik_group_parent(world);
    if (laik_phase(hpcg_instance) > 0)
    {
        Geometry *old_geometry = calculate_old_geometry(A, laik_size(parent));

        A.geom = old_geometry;

        Laik_Partitioner *partitioner_1d = laik_new_partitioner("partitioner_1d_members_of_A", partitioner_1d_members_of_A, (void *)&A, LAIK_PF_None);
        Laik_Partitioner *partitioner_2d = laik_new_partitioner("partitioner_2d_members_of_A", partitioner_2d_members_of_A, (void *)&A, LAIK_PF_None);

        // partitionings before joining: empty for own process
        partitioning_1d_old = laik_new_partitioning(partitioner_1d, parent, A.space, NULL);
        partitioning_2d_old = laik_new_partitioning(partitioner_2d, parent, A.space2d, NULL);

        // Assign current geom again, we have partitioning of old world now
        DeleteGeometry(*old_geometry);
        A.geom = geometry;

        laik_set_initial_partitioning(A.nonzerosInRow_d, partitioning_1d_old);
        laik_set_initial_partitioning(A.matrixDiagonal_d, partitioning_1d_old);
        laik_set_initial_partitioning(A.mtxIndG_d, partitioning_2d_old);
        laik_set_initial_partitioning(A.matrixValues_d, partitioning_2d_old);
    }

    Laik_Partitioner *partitioner_1d = laik_new_partitioner("partitioner_1d_members_of_A", partitioner_1d_members_of_A, (void *)&A, LAIK_PF_None);
    Laik_Partitioner *partitioner_2d = laik_new_partitioner("partitioner_2d_members_of_A", partitioner_2d_members_of_A, (void *)&A, LAIK_PF_None);

    calculate_Mapping(A);

    A.partitioning_1d = laik_new_partitioning(partitioner_1d, world, A.space, NULL);
    A.partitioning_2d = laik_new_partitioning(partitioner_2d, world, A.space2d, NULL);

    Laik_DataFlow flow = LAIK_DF_None;

    // New joining procs need to preserve values
    if (laik_phase(hpcg_instance) > 0)  
        flow = LAIK_DF_Preserve;

    // Split data
    laik_switchto_partitioning(A.nonzerosInRow_d, A.partitioning_1d, flow, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixDiagonal_d, A.partitioning_1d, flow, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixValues_d, A.partitioning_2d, flow, LAIK_RO_None);
    laik_switchto_partitioning(A.mtxIndG_d, A.partitioning_2d, flow, LAIK_RO_None);

    if (laik_phase(hpcg_instance) > 0)
    {
        if(partitioning_1d_old) laik_free_partitioning(partitioning_1d_old);
        if(partitioning_2d_old) laik_free_partitioning(partitioning_2d_old);
    }

    return;
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
 */
allocation_int_t map_l2a_A(const SparseMatrix &A, local_int_t localIndex)
{
    if(localIndex < 0 || localIndex >= A.localNumberOfRows)
        exit_hpcg_run("LocalIndex not mapped to Allocation Index; Fatal Error.", false);

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

    double *matrixDiagonal; laik_get_map_1d(A.matrixDiagonal_d, 0, (void **)&matrixDiagonal, 0);
    double *matrixValues; laik_get_map_1d(A.matrixValues_d, 0, (void **)&matrixValues, 0);
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
                                            matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + currentValuePointer_index] = matrixDiagonal[currentLocalRow]; 
                                       
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
 * @param[in] A
 */
void update_Maps(SparseMatrix &A)
{

    // Removed procs do not need this upate
    if(A.geom->rank < 0) return;
    
    A.localToGlobalMap.clear();
    A.globalToLocalMap.clear();
    A.localToGlobalMap.resize(A.localNumberOfRows);

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

    // Update parameters
    hpcg_params.comm_rank = laik_myid(world);
    hpcg_params.comm_size = laik_size(world);

    // Construct the new geometry and linear system
    Geometry *new_geom = new Geometry;

    // need old gnx to calculate new nx, ny, nz. As npx, npy, npz may change due to new world size
    new_geom->gnx = old_gnx; new_geom->gny = old_gny; new_geom->gnz = old_gnz;

    int dynamicCalculation = hpcg_params.numThreads;
    // use numThreads param to jump into the if to calculate new nx,ny,nz for now
    GenerateGeometry(hpcg_params.comm_size, hpcg_params.comm_rank, -1, hpcg_params.pz, hpcg_params.zl, hpcg_params.zu, hpcg_params.nx, hpcg_params.ny, hpcg_params.nz, hpcg_params.npx, hpcg_params.npy, hpcg_params.npz, new_geom);
    new_geom->numThreads = dynamicCalculation;

    // should be the same problem size
    assert(old_gnx == new_geom->gnx); assert(old_gny == new_geom->gny); assert(old_gnz == new_geom->gnz);

    // Only new geom will be set, other variables are assigned in generateProblem and SetupHalo
    A.geom = new_geom;

    // Additionally, update localNumberOfRows
    A.localNumberOfRows = A.geom->nx * A.geom->ny * A.geom->nz;

    // New joining processes need this data, if broadcast was done in HPCG_init()
    broadcast_hpcg_params();

    return;
}

/**
 * @brief Geometry of coarse grid needs to be updated as well
 *
 * Additionaly, update f2cOperator in MGData 
 * 
 * @param[inout] Af Current coarse grid matrix Layer 
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

    // Update f2cOperator
    // Removed procs do not need this update
    if (Af.geom->rank >= 0)
    {
        local_int_t *f2cOperator = new local_int_t[Af.localNumberOfRows];
        local_int_t localNumberOfRows = nxc * nyc * nzc; // This is the size of our subblock
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
        for (local_int_t i = 0; i < localNumberOfRows; ++i)
            f2cOperator[i] = 0;
#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
        for (local_int_t izc = 0; izc < nzc; ++izc)
        {
            local_int_t izf = 2 * izc;
            for (local_int_t iyc = 0; iyc < nyc; ++iyc)
            {
                local_int_t iyf = 2 * iyc;
                for (local_int_t ixc = 0; ixc < nxc; ++ixc)
                {
                    local_int_t ixf = 2 * ixc;
                    local_int_t currentCoarseRow = izc * nxc * nyc + iyc * nxc + ixc;
                    local_int_t currentFineRow = izf * nxf * nyf + iyf * nxf + ixf;
                    f2cOperator[currentCoarseRow] = currentFineRow;
                } // end iy loop
            }     // end even iz if statement
        }         // end iz loop
        // Delete old f2cOperator
        delete[] Af.mgData->f2cOperator;
        // Assign new f2cOperator
        Af.mgData->f2cOperator= f2cOperator;
    }

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
    assert(old_gnx == geomc->gnx); assert(old_gny == geomc->gny); assert(old_gnz == geomc->gnz);
    return;
}

/**
 * @brief Need to delete old mtxIndL as it will be re-init'ed in SetupHalo(A).
 *
 * We need the old size of localNumberOfRows,thus, we call this function before every other function in repartition_SparseMatrix(A)
 *
 * @param[inout] A SparseMatrix
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
 * @brief Helper function for update_partitionings_x(SparseMatrix &A)
 *
 * @param[in] A SparseMatrix
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
 * @param[inout] A SparseMatrix
 */
void update_partitionings_x(SparseMatrix &A)
{
    // removed procs only need "local" partitioning
    if(A.geom->rank < 0)
    {
        partition_d *pt_local = (partition_d *)malloc(sizeof(partition_d));

        pt_local->halo = false; pt_local->geom = A.geom; pt_local->size = A.totalNumberOfRows;
        pt_local->offset = 0; // don't need this, as well as pt_ext and the rest of the information in pt_data
        // create partitioner
        Laik_Partitioner *x_localPR = laik_new_partitioner("x_localPR_temporary", partitioner_alg_for_x_vector, (void *)pt_local, LAIK_PF_None);
        A.ext = NULL; // Partitioning is deleted later
        // create local partitioning
        A.local = laik_new_partitioning(x_localPR, world, A.space, NULL);
        return;
    }

    // Existing procs will call setuphalo, as everything is calculcated there
    // As mtxIndL is recalculated:
    re_init_mtxIndL(A);
    // TODO clear old variables

    // Update maps of the SparseMatrix
    update_Maps(A);
    SetupHalo(A);
    return;
}

/**
 * @brief Calculate new mapping to due the lex layout
 * 
 * @param[inout] A SparseMatrix
 */
void calculate_Mapping(SparseMatrix &A)
{
    if (A.mapping_) delete[] A.mapping_;
    A.mapping_ = new uint64_t[A.localNumberOfRows]; A.offset_ = -1;

    int rank = A.geom->rank;
    uint64_t localIndex = 0;
    for (long long i = 0; i < A.totalNumberOfRows; i++)
    {
        int proc = ComputeRankOfMatrixRow(*(A.geom), i);

        /* Calculating the offset into the allocation buffer due to lex_layout needs to be done once for all LAIK containers in A */
        if (A.offset_ == -1 && proc == rank)
            A.offset_ = i;

        if (A.offset_ != -1 && proc == rank)
            A.mapping_[localIndex++] = i - A.offset_; // GlobalIndex - Offset = Allocation Index;
    }
    return;
}

/**
 * @brief Partitionings for the data containers in the SparseMatrix need to be updated after repartitioning.
 * 
 * Old ressources will be freed in the end
 *
 * @param[inout] A SparseMatrix
 */
void update_partitionings_A(SparseMatrix &A)
{
    // Removed procs do not need to recalculate this
    if (A.geom->rank >= 0)
        calculate_Mapping(A);

    Laik_Partitioner * partitioner_1d = laik_new_partitioner("partitioner_1d_members_of_A", partitioner_1d_members_of_A, (void *)&A, LAIK_PF_None);
    Laik_Partitioner * partitioner_2d = laik_new_partitioner("partitioner_2d_members_of_A", partitioner_2d_members_of_A, (void *)&A, LAIK_PF_None);

    Laik_Partitioning * new_partitioning_1d = laik_new_partitioning(partitioner_1d, world, A.space, NULL);
    Laik_Partitioning * new_partitioning_2d = laik_new_partitioning(partitioner_2d, world, A.space2d, NULL);

    // redistribute data
    // before, update localLength of layout_data stored in nonzerosInRow_d and matrixDiagonal_d
    Laik_vector_data *layout_data = (Laik_vector_data *)malloc(sizeof(Laik_vector_data));
    layout_data->localLength = A.localNumberOfRows;
    layout_data->numberOfExternalValues = 0; // not used
    layout_data->id = laik_myid(world);
    laik_data_set_layout_data(A.nonzerosInRow_d, layout_data);
    laik_data_set_layout_data(A.matrixDiagonal_d, layout_data); // same layout data
    // laik_data_set_name(A.nonzerosInRow_d, "DEBUG_DATA");
    // printf("LAIK %d \t GO TO 2 \n", laik_myid(world));
    // if (strncmp("DEBUG_DATA", d->name, sizeof("DEBUG_DATA")) == 0)
    // {
    //     printf("LAIK %lu \t Layouts after\n", laik_get_id_vector(fromList->layout));
    //     printf("Old layout %s\n", fromList->layout->describe(fromList->layout));
    //     if (laik_get_id_vector(fromList->layout) == 0)
    //         printf("New layout %s\n", toList->layout->describe(toList->layout));
    //     exit(1);
    // }
    laik_switchto_partitioning(A.nonzerosInRow_d, new_partitioning_1d, LAIK_DF_Preserve, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixDiagonal_d, new_partitioning_1d, LAIK_DF_Preserve, LAIK_RO_None);
    laik_switchto_partitioning(A.matrixValues_d, new_partitioning_2d, LAIK_DF_Preserve, LAIK_RO_None);
    laik_switchto_partitioning(A.mtxIndG_d, new_partitioning_2d, LAIK_DF_Preserve, LAIK_RO_None);

    laik_free_partitioning(A.partitioning_1d);
    laik_free_partitioning(A.partitioning_2d);

    A.partitioning_1d = new_partitioning_1d;
    A.partitioning_2d = new_partitioning_2d; 

    return;
}

/**
 * @brief Local/Global values of the SparseMatrix need to be updated after repartitioning
 *
 * @param[inout] A SparseMatrix
 */
void update_Values(SparseMatrix &A)
{
    // Removed procs do not need this upate
    if (A.geom->rank < 0) return;
        
    char *nonzerosInRow; laik_get_map_1d(A.nonzerosInRow_d, 0, (void **)&nonzerosInRow, 0);

    local_int_t localNumberOfNonZeros = 0;
    for (local_int_t i = 0; i < A.localNumberOfRows; i++)
        localNumberOfNonZeros += nonzerosInRow[i];

    A.localNumberOfNonzeros = localNumberOfNonZeros;

    laik_broadcast((void *)&A.totalNumberOfNonzeros, (void *)&A.totalNumberOfNonzeros, 1, laik_UInt64); // Send this to new procs
    return;
}

/**
 * @brief New joining processes need this data, if broadcast was done in HPCG_init()
 *
 */
void broadcast_hpcg_params(void)
{
    char cparams[][7] = {"--nx=", "--ny=", "--nz=", "--rt=", "--pz=", "--zl=", "--zu=", "--npx=", "--npy=", "--npz="};
    const int nparams = (sizeof cparams) / (sizeof cparams[0]);
    int * iparams = (int *)malloc(sizeof(int) * nparams);
    iparams[0] = hpcg_params.nx;
    iparams[1] = hpcg_params.ny;
    iparams[2] = hpcg_params.nz;
    iparams[3] = hpcg_params.runningTime;
    iparams[4] = hpcg_params.pz;
    iparams[5] = hpcg_params.zl;
    iparams[6] = hpcg_params.zu;
    // npx,npy,npz is 0 in the beginning
    iparams[7] = 0;
    iparams[8] = 0;
    iparams[9] = 0;
    laik_broadcast(iparams, iparams, nparams, laik_Int32);
    return;
}

/**
 * @brief Calculate old geometry. We need this to set initial partitioning for new joining processes
 * 
 * @param[in] A SparseMatrix
 * @param[in] old_size of parent world 
 * @return Pointer to geom with old config 
 */
Geometry * calculate_old_geometry(SparseMatrix &A, int old_size)
{
    Geometry *old_geometry = new Geometry;
    // // Calculate nx, ny, nz of old world
    int npx = 0; int npy = 0; int npz = 0;
    ComputeOptimalShapeXYZ(old_size, npx, npy, npz);
    local_int_t old_nx, old_ny, old_nz;
    old_nx = A.geom->gnx / npx; old_ny = A.geom->gny / npy; old_nz = A.geom->gnz / npz;
    // Test
    bool config_1 = A.geom->gnx % npx == 0; bool config_2 = A.geom->gny % npy == 0; bool config_3 = A.geom->gnz % npz == 0;
    // this tests are not needed, because we are asserting the previous config here...
    if (!config_1 || !config_2 || !config_3)
    {
        // This means, that expanding/shrinking will not work to the demanded size,
        assert(config_1 == true); // will fail
        exit_hpcg_run("It is not possible to expand/shrink the world as requested. Try other new sizes!", false);
    }
    GenerateGeometry(old_size, 0, hpcg_params.numThreads, hpcg_params.pz, hpcg_params.zl, hpcg_params.zu, old_nx, old_ny, old_nz, hpcg_params.npx, hpcg_params.npy, hpcg_params.npz, old_geometry);
    return old_geometry;
}

/**
 * @brief Calculate old coarse geometry. We need this to set initial partitioning for new joining processes
 *
 * @param[in] Af SparseMatrix
 * @param[in] old_size of parent world
 * @return Pointer to geom with old config
 */
Geometry * calculate_old_coarse_geometry(Geometry * old_geometry, int old_size)
{
    // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
    // below may result in global range values.
    global_int_t nxf = old_geometry->nx;
    global_int_t nyf = old_geometry->ny;
    global_int_t nzf = old_geometry->nz;

    local_int_t nxc, nyc, nzc; //Coarse nx, ny, nz
    assert(nxf%2==0); assert(nyf%2==0); assert(nzf%2==0); // Need fine grid dimensions to be divisible by 2
    nxc = nxf/2; nyc = nyf/2; nzc = nzf/2;

    // Construct the geometry and linear system
    Geometry * old_geomc = new Geometry;
    local_int_t zlc = 0; // Coarsen nz for the lower block in the z processor dimension
    local_int_t zuc = 0; // Coarsen nz for the upper block in the z processor dimension
    int pz = old_geometry->pz;
    if (pz > 0)
    {
        zlc = old_geometry->partz_nz[0] / 2; // Coarsen nz for the lower block in the z processor dimension
        zuc = old_geometry->partz_nz[1] / 2; // Coarsen nz for the upper block in the z processor dimension
    }
    GenerateGeometry(old_size, 0, old_geometry->numThreads, old_geometry->pz, zlc, zuc, nxc, nyc, nzc, old_geometry->npx, old_geometry->npy, old_geometry->npz, old_geomc);
    return old_geomc;
}

/*
    Helper functions for shrink/expand feature -END
        Forward declare them in the beginning of the file
*/
#endif