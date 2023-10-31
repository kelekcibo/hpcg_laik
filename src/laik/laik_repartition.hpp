/**
 * @file laik_repartition.hpp
 * @brief Everything needed for enabling repartitioning
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef LAIK_REPARTITION_HPP
#define LAIK_REPARTITION_HPP

/*
    Includes
*/
#include <vector>

// forw. decl
struct Laik_Blob;
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
    Needed functions/variables for shrink/expand feature
*/
#ifdef REPARTITION
extern HPCG_Params hpcg_params;
const local_int_t numberOfNonzerosPerRow = 27; // We are approximating a 27-point finite element/volume/difference 3D stencil

extern allocation_int_t map_l2a_A(const SparseMatrix &A, local_int_t localIndex);
extern void init_SPM_partitionings(SparseMatrix &A);
extern void repartition_SparseMatrix(SparseMatrix &A);
extern void re_switch_LaikVectors(SparseMatrix &A, std::vector<Laik_Blob *> list);
extern void replaceMatrixValues(SparseMatrix &A);
extern void new_joining_procs(Laik_RangeReceiver *r, Laik_PartitionerParams *p);
extern void update_Maps(SparseMatrix &A);
extern void update_Values(SparseMatrix &A);
extern void calculate_Mapping(SparseMatrix &A);
extern void re_init_mtxIndL(SparseMatrix &A);

#endif
/*
    Needed functions/variables for shrink/expand feature -END
*/
#endif // LAIK_REPARTITION_HPP
