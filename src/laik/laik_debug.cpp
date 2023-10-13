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
#include <iostream>

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
 * @param x
 * @param y
 * @param doIO
 * @param curIndex
 */
void compare2(double x, double y, bool doIO, allocation_int_t curIndex)
{
    double delta = std::abs(x - y);

    if (doIO)
        printf("map_l2a index: %lld\t| %.10f - %.10f | = %.10f\n", curIndex, x, y, delta);

    if (delta != 0.0)
    {
        if (doIO)
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

    if (spm->geom->rank != 0)
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
    if (doIO)
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

    return;
}

/**
 * @brief Print members of geometry struct
 *
 * @param geom
 * @param doIO
 */
void print_GEOMETRY(Geometry *geom, bool doIO)
{
    if (doIO)
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

            if (i != geom->npartz - 1)
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

void print_LaikBlob(Laik_Blob *x)
{
    if (x != NULL)
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
void exit_hpcg_run(const char *msg)
{
    if (laik_myid(world) == 0)
        printf("\n####### %s\n####### Debug DONE -> Exiting #######\n", msg);

    exit(0);

    return;
}

/*
    Implementation of debug functions -END
*/