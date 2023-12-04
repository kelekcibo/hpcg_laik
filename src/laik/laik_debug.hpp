/**
 * @file laik_debug.hpp
 * @brief Debug functions
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef LAIK_DEBUG_HPP
#define LAIK_DEBUG_HPP

/*
    Includes 
*/
#include <iostream>

using std::to_string;

#include "laik_x_vector.hpp"

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

/*
    Debug functions
*/
extern void compare2(double x, double y, bool doIO, allocation_int_t curIndex);
extern void compareResult(Vector &x, Laik_Blob *y, bool doIO);
extern void printResultLaikVector(Laik_Blob *x);
extern void printResultVector(Vector &x);
extern void printSPM(SparseMatrix *A, int coarseLevel);
extern void printSPM_val(SparseMatrix &A);
extern void print_HPCG_PARAMS(HPCG_Params params, bool doIO);
extern void print_GEOMETRY(Geometry *geom, bool doIO);
extern void print_LaikBlob(Laik_Blob *x);
extern void exit_hpcg_run(const char *msg, bool wait);
/*
    Debug functions -END
*/

#endif // LAIK_DEBUG_HPP
