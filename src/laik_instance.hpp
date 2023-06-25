/**
 * @brief Helper Functions + global variable instance + world (reason: using laik would result in more ifdef's)
 */
#ifndef LAIK_HELPER_HPP
#define LAIK_HELPER_HPP

extern "C"{
    #include <laik.h>
}

extern Laik_Instance *hpcg_instance;
extern Laik_Group * world;

extern void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type);
extern void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type);

#endif