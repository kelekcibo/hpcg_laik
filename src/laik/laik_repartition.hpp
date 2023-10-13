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
#include "laik/hpcg_laik.hpp"
#include <vector>

// forw. decl
struct Laik_Blob;
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;
// forw. decl

#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "Geometry.hpp"
#include "hpcg.hpp"
/*
    Includes -END
*/

extern HPCG_Params hpcg_params;

extern void re_setup_problem(SparseMatrix &A);
extern void re_switch_LaikVectors(SparseMatrix &A, std::vector<Laik_Blob *> list);

#endif // LAIK_REPARTITION_HPP
