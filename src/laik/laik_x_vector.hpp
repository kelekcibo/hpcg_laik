/**
 * @file laik_x_vector.hpp
 * @brief Everything needed for using Laik Vectors
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef LAIK_X_VECTOR_HPP
#define LAIK_X_VECTOR_HPP

/*
    Includes
*/
#include <vector>
#include <set>
#include <map>

// forw. decl
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;
// forw. decl

#include "laik/hpcg_laik.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "Geometry.hpp"
#include "hpcg.hpp"
/*
    Includes -END
*/

/*
    Structs
*/
/**
 * @brief Struct (data) needed for the partitioner algorithm.
 *
 * if halo == true:
 *      This data is needed if a process needs to access external values.
 * else:
 *      This data is needed when a process should access only local values
 *
 * @param size                  size of vector x
 * @param elementsToSend        local index to vector x of elements to be sent (converting to global later)
 * @param geom                  Geometry of the problem
 * @param numberOfNeighbours    count neighbors we need to exchange data with
 * @param neighbors             Process ID's of neighbours
 * @param receiveLength         number of elements to be received/sent by/to neighbours
 * @param receiveList           Global Indices proc i needs to receive
 * @param localToGlobalMap      local-to-global mapping
 * @param halo                  x will be partioned such that indices to external values are owned by other procs as well
 * @param offset                offset to allocation buffer of partitioning with local values
 */
struct partition_data
{
    local_int_t size;                                  /* size of vector x */
    local_int_t *elementsToSend;                       /* local index to vector x of elements to be sent (converting to global later) */
    Geometry *geom;                                    /* Geometry of the problem */
    int numberOfNeighbours;                            /* count neighbors we need to exchange data with*/
    int *neighbors;                                    /* Process ID's of neighbours */
    local_int_t *receiveLength;                        /* number of elements to be received/sent by/to neighbours */
    std::map<int, std::set<global_int_t>> receiveList; /* Global Indices proc i needs to receive */
    std::vector<global_int_t> *localToGlobalMap;       /* local-to-global mapping */
    bool halo;                                         /* x will be partioned such that indices to external values are owned by other procs as well */
    int offset;                                        /* offset to allocation buffer of partitioning */
};
typedef partition_data partition_d;

/**
 * @brief Laik Vector
 * 
 * @param name for debugging
 * @param values LAIK container holding the data
 * @param localLength of the data
 */
struct Laik_Blob
{
    char *name; // name of the vector /* Debug

    Laik_Data *values;
    bool exchangesValues;
    mutable local_int_t localLength;
};

/*
    Structs -END
*/

/*
    Functions needed to exchange values via LAIK 
*/
extern void partitioner_alg_for_x_vector(Laik_RangeReceiver *r, Laik_PartitionerParams *p);
extern void init_partition_data(SparseMatrix &A, partition_d *local, partition_d *ext);
extern void init_partitionings(SparseMatrix &A, partition_d *local, partition_d *ext);
extern Laik_Blob *init_blob(const SparseMatrix &A, bool exchangesValues, char * name);
/*
    Functions needed to exchange values via LAIK -END
*/

/*
    Operations on laik vectors
*/
extern void fillRandomLaikVector(Laik_Blob *x);
extern void ZeroLaikVector(Laik_Blob *x);
extern void CopyLaikVectorToLaikVector(Laik_Blob *x, Laik_Blob *y);
extern void CopyVectorToLaikVector(Vector &v, Laik_Blob *x);
extern void CopyLaikVectorToVector(const Laik_Blob *x, Vector &v);
extern void CopyLaikVectorToVector(Laik_Blob *x, Vector &v);
extern void ScaleLaikVectorValue(Laik_Blob *v, local_int_t index, double value);
/*
    Operations on laik vectors -END
*/

/*
    clean up functions
*/
extern void DeleteLaikVector(Laik_Blob *x);
/*
    clean up functions -END
*/

#endif // LAIK_X_VECTOR_HPP