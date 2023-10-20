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
 * @brief Making use of the lex_layout, we need following mapping of indices:
 *
 * Local indices -> Global indices -> Allocation indices
 *
 * As LAIK allocates a buffer under the hood according to the lex layout, we need to specify a mapping for the corresponding
 * indices within that buffer.
 * 
 * This mapping is for Laik Vectors
 *
 */
struct Local2Allocation_map_x
{
    allocation_int_t offset_ext; /* Offset into allocation buffer with external values */
    allocation_int_t offset;     /* Offset into allocation buffer without external values */

    /* Mapping from Local to Global Indices */
    std::map<local_int_t, global_int_t> localToExternalMap; /* External global indices */
    std::vector<global_int_t> localToGlobalMap;             /* Owned global indices */
    local_int_t localNumberOfRows;                          /* Border between owned and external indices */
};
typedef Local2Allocation_map_x L2A_map;

/**
 * @brief Laik Vector
 * 
 * @param name for debugging
 * @param values LAIK container holding the data
 * @param localLength of the data
 */
struct Laik_Blob
{
    const char *name; // name of the vector /* Debug

    Laik_Data *values;
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
extern Laik_Blob *init_blob(const SparseMatrix &A);
extern allocation_int_t map_l2a_x(L2A_map *mapping, local_int_t local_index, bool halo);
/*
    Functions needed to exchange values via LAIK -END
*/

/*
    Operations on laik vectors
*/
extern void fillRandomLaikVector(Laik_Blob *x, L2A_map *mapping);
extern void ZeroLaikVector(Laik_Blob *x, L2A_map *mapping);
extern void CopyLaikVectorToLaikVector(Laik_Blob *x, Laik_Blob *y, L2A_map *mapping);
extern void CopyVectorToLaikVector(Vector &v, Laik_Blob *x_blob, L2A_map *mapping);
extern void CopyLaikVectorToVector(const Laik_Blob *x_blob, Vector &v, L2A_map *mapping);
extern void CopyLaikVectorToVector(Laik_Blob *x, Vector &v, L2A_map *mapping);
extern void ScaleLaikVectorValue(Laik_Blob *v, local_int_t index, double value, L2A_map *mapping);
/*
    Operations on laik vectors -END
*/

/*
    clean up functions
*/
extern void free_L2A_map(L2A_map *mapping);
extern void DeleteLaikVector(Laik_Blob *x);
/*
    clean up functions -END
*/

#endif // LAIK_X_VECTOR_HPP