/**
 * @file hpcg_laik.cpp
 *
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */

/*
    Includes
*/
    #include "hpcg_laik.hpp"
/*
    Includes -END
*/

/*
    Important global variables
*/
// should be initialized at the very beginning of the program. No use without init.
Laik_Instance *hpcg_instance; // The instance the world is associated to
Laik_Group *world; // group of processes
/*
    Important global variables -END
*/