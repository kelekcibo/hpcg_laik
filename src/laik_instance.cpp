/**
 * @brief Definiton for helper functions
 * 
 */

#include "laik_instance.hpp"

#include <iostream>
#include <cstring>
#include <cassert>


// should be initialized at the very beginning of the program. No use without init.
Laik_Instance *hpcg_instance; /* Laik instance during HPCG run */
Laik_Group *world;

/* Space and containers */
Laik_Space *x_space; /* space of vector x */
Laik_Data *x_vector; /* Proc i local portion */
Laik_Data *x_vector_halo; /* Proc i local portion + halo values */

/* Partitioning to switch between */
Laik_Partitioning *x_pt;  
Laik_Partitioning *x_halo_pt;


void partitioner_x(Laik_RangeReceiver *r, Laik_PartitionerParams *p)
{

    pt_data * data = (pt_data *) laik_partitioner_data(p->partitioner);
    int id1 = laik_myid(world);

    Laik_Range range;

    for(uint32_t id = 0; id < laik_size(world); id++)
    {
        // assign every process its local portion of the x vector
        laik_range_init_1d(&range, x_space, data->local_portion * id, data->local_portion * id + data->local_portion);
        laik_append_range(r, id, &range, 0, 0);
    }

    if (data->halo)
    {
        local_int_t index_elementsToSend = 0;

        // Process i needs to send local values
        for (int nb = 0; nb < data->numberOfNeighbours; nb++)
        {
            local_int_t numbOfElementsToSend = data->receiveLength[nb]; // due to symmetry, this is equal to sendLength[nb]

            for(local_int_t i = 0; i < numbOfElementsToSend; i++)
            {

                local_int_t j = data->elementsToSend[index_elementsToSend++]; // neighbour n needs this value from proc i
                j += data->local_portion * laik_myid(world);                  // offset to local portion of proc i

                printf("I (%d) need to give access to proc %d at local offset %d\n", id1, data->neighbors[nb], j - (data->local_portion * laik_myid(world)));

                laik_range_init_1d(&range, x_space, j, j+1);
                laik_append_range(r, data->neighbors[nb], &range, 0, 0);
            }
        }


        typedef std::set<global_int_t>::iterator set_iter;

        // process i needs to get external values
        for (int nb = 0; nb < data->numberOfNeighbours; nb++)
        {
            for (set_iter i = data->receiveList[data->neighbors[nb]].begin(); i != data->receiveList[data->neighbors[nb]].end(); ++i)
            {
                // TODO: How to get mapping from global to local of external values?
                int globalIndex = *i; // receiveList mit local indices speichern

                laik_range_init_1d(&range, x_space, globalIndex, globalIndex + 1);
                laik_append_range(r, id1, &range, 0, 0);

                printf("I (%d) need to have access to local offset %d of proc %d \n", id1, globalIndex, data->neighbors[nb]);
            }       
        }
    }

    // printf("Partitioning with halo (%s) intialized!\n", (data->halo == 1 ? "true":"false"));
}


/**
 * @brief Initializing the two partitiongs to switch between them (exchanging data between processes)
 * 
 * @param data - necessary information needed for the partitioner algorithm
 */
void init_partitionings(pt_data * data)
{
    x_space = laik_new_space_1d(hpcg_instance, data->size);
    x_vector = laik_new_data(x_space, laik_Double);

    // Initialize partitioners to partition data accordingly
    data->halo = false;
    Laik_Partitioner *x_pr = laik_new_partitioner("x_pt", partitioner_x, (void *)data, LAIK_PF_None);

    pt_data data2;
    std::memcpy((void *) &data2, (void *)data, sizeof(pt_data));

    data2.halo = true;
    Laik_Partitioner *x_halo_pr = laik_new_partitioner("x_pt_halo", partitioner_x, (void *) &data2, LAIK_PF_None);

    x_halo_pt = laik_new_partitioning(x_halo_pr, world, x_space, NULL);
    x_pt = laik_new_partitioning(x_pr, world, x_space, NULL); 
}


/**
 * @brief Switch between partitionings to exchange data
 * 
 * @param halo allow access to external values, if true
 */
void exchangeValues(bool halo){

    if(halo)
    {
        laik_switchto_partitioning(x_vector, x_halo_pt, LAIK_DF_Preserve, LAIK_RO_None);
        return;
    }

    laik_switchto_partitioning(x_vector, x_pt, LAIK_DF_Preserve, LAIK_RO_None);
}

/**
 * @brief Helper function implements Allreduce / Broadcast
 * 
 * @param partitioner needed to specify Allreduce or Broadcast
 */
void laik_helper(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type, Laik_ReductionOperation ro_type, Laik_Partitioner * partitioner)
{

    // reducing some buffer occurs very often in hpcg code base
    // implement it here once and call it from the hpcg code base

    // Definition and initialization of laik-specific data
    Laik_Space * space = laik_new_space_1d(hpcg_instance, n);
    Laik_Data * data = laik_new_data(space, data_type);

    // If broadcast: Root process has only access to container
    // else: everyone
    laik_switchto_new_partitioning(data, world, partitioner, LAIK_DF_None, LAIK_RO_None);

    // fill sendBuf into container
    uint64_t count;
    if(data_type == laik_Int32)
    {
        int * sendBuf_tmp = (int *) sendBuf;
        int * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            base[i] = sendBuf_tmp[i];

    }else if(data_type == laik_Double)
    { 
        double * sendBuf_tmp = (double *)sendBuf;
        double * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            base[i] = sendBuf_tmp[i];
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t * sendBuf_tmp = (uint64_t *)sendBuf;
        uint64_t * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            base[i] = sendBuf_tmp[i];

    }

    // Broadcast/Allreduce data with ro_type as reduction operation
    laik_switchto_new_partitioning(data, world, laik_All, LAIK_DF_Preserve, ro_type);

    // Store received result in recvBuf
    if (data_type == laik_Int32)
    {
        int * recvBuf_tmp = (int *)sendBuf;
        int * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            recvBuf_tmp[i] = base[i];

        assert(count == n); // elements in buffer; should be equal
        std::memcpy(recvBuf, recvBuf_tmp, count * sizeof(int32_t));
    }
    else if (data_type == laik_Double)
    {
        double * recvBuf_tmp = (double *)sendBuf;
        double * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            recvBuf_tmp[i] = base[i];

        assert(count == n); // elements in buffer; should be equal
        std::memcpy(recvBuf, recvBuf_tmp, count * sizeof(double));
    }
    else if (data_type == laik_UInt64)
    {
        uint64_t * recvBuf_tmp = (uint64_t *)sendBuf;
        uint64_t * base;
        laik_get_map_1d(data, 0, (void **)&base, &count);

        for (uint64_t i = 0; i < count; ++i)
            recvBuf_tmp[i] = base[i];

        assert(count == n); // elements in buffer; should be equal
        std::memcpy(recvBuf, (void *)recvBuf_tmp, count * sizeof(uint64_t));
    }

    // Free ressources
    laik_free(data);
    laik_free_space(space);
    
    return;
}

/**
 * @brief "Allreduce" a buffer based on the ro_type.
 *
 * @param sendBuf send buffer
 * @param recvBuf receive buffer
 * @param n size of send/recv buffer
 * @param data_type of the send buffer to be sent
 * @param ro_type type of reduction to be applied
 * @param allreduce if true, use allreduce
 */
void laik_allreduce(const void * sendBuf, void * recvBuf, uint64_t n, Laik_Type * data_type, Laik_ReductionOperation ro_type){
    laik_helper(sendBuf, recvBuf, n, data_type, ro_type, laik_All);
    return;
}

/**
 * @brief Broadcast a buffer to all processes in the world (Broadcast done by root process (rank==0)
 *
 * @param sendBuf send buffer
 * @param recvBuf receive buffer
 * @param n size of buffer
 * @param data_type of the buffer to be broadcasted
 */
void laik_broadcast(const void *sendBuf, void *recvBuf, uint64_t n, Laik_Type *data_type)
{
    laik_helper(sendBuf, recvBuf, n, data_type, LAIK_RO_None, laik_Master);
    return;
}

/**
 * @brief Synchronize all processes
 * 
 */
void laik_barrier()
{
    // arbitrary data
    int32_t data = 71; 

    // Synchronize all processes by making use of an All-to-All-Reduction
    laik_helper(&data, &data, 1, laik_Int32, LAIK_RO_None, laik_All);
    return;
}