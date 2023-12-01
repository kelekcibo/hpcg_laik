/**
 * @file laik_debug.cpp
 * @brief Implementation of debug functions
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

#include "laik_debug.hpp"

/*
    Includes
*/

/*
    Implementation of debug functions
*/
/**
 * @brief Compare two double values x and y
 *
 * @param[in] x value
 * @param[in] y value
 * @param[in] doIO print
 * @param[in] curIndex index to the allocation buffer
 */
void compare2(double x, double y, bool doIO, allocation_int_t curIndex)
{
    double delta = std::abs(x - y);

    if (doIO) printf("map_l2a index: %lld\t| %.10f - %.10f | = %.10f\n", curIndex, x, y, delta);

    if (delta != 0.0)
    {
        if (doIO) printf("Difference is not tolerated: %.20f\nindex of map_l2a: %lld\n", delta, curIndex);
        exit_hpcg_run("delta does not equal 0!", false);
    }
}

/**
 * @brief Compare the two vectors x and y.
 *
 * @param[in] x vector
 * @param[in] y vector
 * @param[in] mapping due to the lex layout
 * @param[in] doIO print
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
        double delta = std::abs(xv[i] - yv[map_l2a_x(mapping, i, false)]);
        // if (doIO) printf("Index %lld: Delta %.10f\n", i, delta);
        if (doIO) printf("xv[%ld]=%.10f\tyv_blob[%lld]=%.10f\n", i, xv[i], map_l2a_x(mapping, i, false), yv[map_l2a_x(mapping, i, false)]);
        if (delta != 0)
        {
            if (doIO) printf("Difference is not tolerated: %.20f\n", delta);
            exit_hpcg_run("delta does not equal 0!", false);
        }
    }

    if (doIO) printf("Compare done\n");
}

/**
 * @brief Print the vector
 *
 * @param[in] x vector to be printed
 */
void printResultVector(Vector &x)
{
    if (laik_myid(world) == 0)
    {
        // printf("\n\nPrint result of vector\n");
        HPCG_fout << "\n\nPrint result of vector\n";

        double *xv = x.values;
        size_t length = x.localLength;

        // printf("Length = %ld\n", length);
        HPCG_fout << "Length = " << to_string(length) << "\n";

        for (size_t i = 0; i < length; i++)
            // printf("xv[%ld]=%.10f\n", i, xv[i]);
            HPCG_fout << "xv[" << to_string(i) << "]=" << to_string(xv[i]) << "\n";

        // printf("\nEnd of printing result of vector\n\n");
        HPCG_fout << "\nEnd of printing result of vector\n\n";
    }
}

/**
 * @brief Print the Laik Vector
 *
 * @param[in] x laik vector to be printed
 * @param[in] mapping due to the lex layout
 */
void printResultLaikVector(Laik_Blob *x)
{
    if (laik_myid(world) == 0)
        // HPCG_fout << "\n\nPrint result of vector\n";
    printf("\n\nPrint result of Laik-vector\n");

    double *xv;
    laik_get_map_1d(x->values, 0, (void **)&xv, 0);

    size_t localLength = x->localLength;

    if (laik_myid(world) == 0)
    {
        printf("localLength = %ld\n", localLength);
        for (size_t i = 0; i < localLength; i++)
            printf("xv[%ld]=%.10f\n", i, xv[i]);
    printf("\nEnd of printing result of vector\n\n");
    }

    // if (laik_myid(world) == 0)
    // {
    //     HPCG_fout << "localLength = " << to_string(localLength) << "\n";
    //     for (size_t i = 0; i < localLength; i++)
    //         HPCG_fout << "xv[" << to_string(i) << "]=" << to_string(xv[map_l2a_x(mapping, i, false)]) << "\n";
    //     HPCG_fout << "\nEnd of printing result of vector\n\n";
    // }
}

/**
 * @brief Print some information about the SparseMatrix
 *
 * @param[in] spm sparsematrix
 * @param[in] coarseLevel current matrix layer/level
 */
void printSPM(SparseMatrix *spm, int coarseLevel)
{

    // Global data
    HPCG_fout << "\n##################### Global stats #####################\n\n";

    HPCG_fout << "\nTotal # of rows " << spm->totalNumberOfRows
              << std::endl
              << "\nTotal # of nonzeros " << spm->totalNumberOfNonzeros
              << std::endl;

    // Local
    HPCG_fout << "\n##################### Local stats #####################\n\n";
    HPCG_fout << "\nLocal # of rows " << spm->localNumberOfRows
              << std::endl
              << "\nLocal # of nonzeros " << spm->localNumberOfNonzeros
              << std::endl;
    // mapping of rows
    HPCG_fout << "\n##################### Mapping of rows #####################\n\n";
    HPCG_fout << "\nLocal-to-global Map\n";
    HPCG_fout << "Local\tGlobal\n";
    for (int c = 0; c < spm->localToGlobalMap.size(); c++)
        HPCG_fout << c << "\t\t" << spm->localToGlobalMap[c] << std::endl;

    // other procs print on stdin
    if (spm->geom->rank != 0)
    {
        std::cout << "\n##################### My RANK (" << spm->geom->rank << ") #####################\n\n";
        // Global data
        std::cout << "\n##################### Global stats #####################\n\n";

        std::cout << "\nTotal # of rows " << spm->totalNumberOfRows
                  << std::endl
                  << "\nTotal # of nonzeros " << spm->totalNumberOfNonzeros
                  << std::endl;

        // Local
        std::cout << "\n##################### Local stats #####################\n\n";
        std::cout << "\nLocal # of rows " << spm->localNumberOfRows
                  << std::endl
                  << "\nLocal # of nonzeros " << spm->localNumberOfNonzeros
                  << std::endl;

        // mapping of rows
        std::cout << "\n##################### Mapping of rows #####################\n\n";
        std::cout << "\nLocal-to-global Map\n";
        std::cout << "Local\tGlobal\n";
        for (int c = 0; c < spm->localToGlobalMap.size(); c++)
            std::cout << c << "\t\t" << spm->localToGlobalMap[c] << std::endl;
        std::cout << "\n##################### My RANK (" << spm->geom->rank << ") END #####################\n\n";
    }

    return;
}

/**
 * @brief Print values of matrixValues and matrixDiagonal members of A.
 * 
 * This is to print the values, if repartitioning is enabled. Then we use Laik_Data members.
 *
 * @param[in] A SparseMatrix
 */
void printSPM_val(SparseMatrix &A)
{
#ifdef REPARTITION
    global_int_t nx = A.geom->nx;
    global_int_t ny = A.geom->ny;
    global_int_t nz = A.geom->nz;
    global_int_t gnx = A.geom->gnx;
    global_int_t gny = A.geom->gny;
    global_int_t gnz = A.geom->gnz;
    global_int_t gix0 = A.geom->gix0;
    global_int_t giy0 = A.geom->giy0;
    global_int_t giz0 = A.geom->giz0;

    const char *nonzerosInRow; laik_get_map_1d(A.nonzerosInRow_d, 0, (void **)&nonzerosInRow, 0);
    const double *matrixValues; laik_get_map_1d(A.matrixValues_d, 0, (void **)&matrixValues, 0);
    const double *matrixDiagonal; laik_get_map_1d(A.matrixDiagonal_d, 0, (void **)&matrixDiagonal, 0);
    double entry_val = 0.0; double entry_dia = 0.0;

    std::string debug{""};
    for (local_int_t iz = 0; iz < nz; iz++)
    {
        global_int_t giz = giz0 + iz;
        for (local_int_t iy = 0; iy < ny; iy++)
        {
            global_int_t giy = giy0 + iy;
            for (local_int_t ix = 0; ix < nx; ix++)
            {
                global_int_t gix = gix0 + ix; local_int_t currentLocalRow = iz * nx * ny + iy * nx + ix;
                global_int_t currentGlobalRow = giz * gnx * gny + giy * gnx + gix;

                debug += "Current Local Row (" + std::to_string(currentLocalRow) + ") " 
                      + "Current Global Row (" + std::to_string(currentGlobalRow) + ") " 
                      + "cur_nnz (" + std::to_string(nonzerosInRow[currentLocalRow]) + ") ";
                debug += "\nUsed Matrix values: ";
                uint64_t currentValuePointer_index = -1;      // Index to current value in current row
                global_int_t currentIndexPointerG_index = -1; // Index to current index in current row
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
                                        if (curcol == currentGlobalRow)
                                        {
                                            debug += std::to_string(matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + ++currentValuePointer_index]) + " [dia ";
                                            debug += std::to_string(matrixDiagonal[currentLocalRow]) + "], ";
                                        }
                                        else
                                        {
                                            debug += std::to_string(matrixValues[map_l2a_A(A, currentLocalRow) * numberOfNonzerosPerRow + ++currentValuePointer_index]) + ", ";
                                        }

                                    } // end x bounds test
                                }     // end sx loop
                            }         // end y bounds test
                        }             // end sy loop
                    }                 // end z bounds test
                }                     // end sz loop

                debug += "\n";
            } // end ix loop
        }     // end iy loop
    }         // end iz loop
    HPCG_fout << debug;
#endif
    return;
}

/**
 * @brief Print members of HPCG_Params Struct
 *
 * @param[in] params HPCG_Params
 * @param[in] doIO print
 */
void print_HPCG_PARAMS(HPCG_Params params, bool doIO)
{
    if (doIO)
    {
        std::string a{"####### PARAM values\n"};
        a += "rank: " + to_string(params.comm_rank) + "\n";
        a += "npx: " + to_string(params.npx) + "\n";
        a += "npy: " + to_string(params.npy) + "\n";
        a += "npz: " + to_string(params.npz) + "\n";
        a += "numThreads: " + to_string(params.numThreads) + "\n";
        a += "nx: " + to_string(params.nx) + "\n";
        a += "ny: " + to_string(params.ny) + "\n";
        a += "nz: " + to_string(params.nz) + "\n";
        a += "pz: " + to_string(params.pz) + "\n";
        a += "runningTime: " + to_string(params.runningTime) + "\n";
        a += "zl: " + to_string(params.zl) + "\n";
        a += "zu: " + to_string(params.zu) + "\n";
        std::cout << a;
    }

    return;
}

/**
 * @brief Print members of geometry struct
 *
 * @param[in] geom
 * @param[in] doIO
 */
void print_GEOMETRY(Geometry *geom, bool doIO)
{
    if (doIO)
    {
        std::string a{"####### Geometry values\t"};
        a += "printed by rank " + to_string(geom->rank) + "\n";
        a += "gix0: " + to_string(geom->gix0) + "\n";
        a += "giy0: " + to_string(geom->giy0) + "\n";
        a += "giz0: " + to_string(geom->giz0) + "\n";
        a += "gnx: " + to_string(geom->gnx) + "\n";
        a += "gny: " + to_string(geom->gny) + "\n";
        a += "gnz: " + to_string(geom->gnz) + "\n";
        a += "ipx: " + to_string(geom->ipx) + "\n";
        a += "ipy: " + to_string(geom->ipy) + "\n";
        a += "ipz: " + to_string(geom->ipz) + "\n";
        a += "npartz: " + to_string(geom->npartz) + "\n";
        a += "npx: " + to_string(geom->npx) + "\n";
        a += "npy: " + to_string(geom->npy) + "\n";
        a += "npz: " + to_string(geom->npz) + "\n";
        a += "numThreads: " + to_string(geom->numThreads) + "\n";
        a += "nx: " + to_string(geom->nx) + "\n";
        a += "ny: " + to_string(geom->ny) + "\n";
        a += "nz: " + to_string(geom->nz) + "\n";
        a += "pz: " + to_string(geom->pz) + "\n";
        a += "partz_ids: ";
        for (int i = 0; i < geom->npartz; ++i)
        {
            a += to_string(geom->partz_ids[i]);
            if (i != geom->npartz - 1) a += ", ";
        }
        a += "\n";
        a += "partz_nz: ";
        for (int i = 0; i < geom->npartz; ++i)
        {
            a += to_string(geom->partz_nz[i]);
            if (i != geom->npartz - 1) a += ", ";
        }
        a += "\n";
        a += "size: " + to_string(geom->size) + "\n";
        std::cout << a;
    }
    return;
}

/**
 * @brief Print information stored in a Laik_Blob
 *
 * @param[in] x Laik_Blob to be printed
 */
void print_LaikBlob(Laik_Blob *x)
{
    if (x != NULL)
    {
        std::string str{"####### Laik_Blob\n"};
        std::string str2{x->name};
        str += "LAIK " + to_string(laik_myid(world)) + "\n";
        str += "Vector name: " + str2 + "\n";
        str += "localLength: " + to_string(x->localLength) + "\n#######\n";
        std::cout << str;
    }
    return;
}

/**
 * @brief Debug function.
 *
 * Exit the program on all processes
 *
 * @param[in] msg will be printed by proc with ID 0 (No need to add "\\n")
 * @param[in] wait tells, if proc should wait instead of exiting
 */
void exit_hpcg_run(const char *msg, bool wait)
{
    if (laik_myid(world) == 0)
    {
        if (strcmp(msg, "") != 0) printf("\n\n####### %s\n####### Debug DONE -> Exiting #######\n", msg);
        else printf("\n\n####### Debug DONE -> Exiting #######\n");
    }
    if(wait) while(1) ;
    laik_finalize(hpcg_instance); exit(0);
    return;
}

/*
    Implementation of debug functions -END
*/