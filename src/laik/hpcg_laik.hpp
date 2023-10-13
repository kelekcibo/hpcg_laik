/**
 * @file hpcg_laik.hpp
 * 
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef HPCG_LAIK_HPP
#define HPCG_LAIK_HPP

/*
    Defs
*/
#ifndef REPARTITION
#define REPARTITION
#endif

typedef long long allocation_int_t;
/*
    Defs -END
*/

/*
    Includes
*/
extern "C"
{
#include <laik.h> // using laik library
}

// All other custom laik headers within src/laik
#include "laik_x_vector.hpp"
#include "laik_debug.hpp"
#include "laik_reductions.hpp"
#ifdef REPARTITION
#include "laik_repartition.hpp"
#endif
/*
    Includes -END
*/

/*
    Important global variables
*/
// should be initialized at the very beginning of the program. No use without init.
extern Laik_Instance *hpcg_instance;
extern Laik_Group *world;
/*
    Important global variables -END
*/

#endif // HPCG_LAIK_HPP