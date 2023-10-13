/**
 * @file laik_reductions.hpp
 * @brief Using LAIK for Allreduction-/Broadcastoperations
 * @version 1.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */
#ifndef LAIK_REDUCTIONS_HPP
#define LAIK_REDUCTIONS_HPP


/*
    Includes
*/
#include "laik/hpcg_laik.hpp"
/*
    Includes -END
*/

/*
    Functions
*/
extern void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type);
extern void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type);
extern void laik_barrier(void);
/*
    Functions -END
*/

#endif // LAIK_REDUCTIONS_HPP