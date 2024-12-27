/*! @file PusherBoris_2D.h

 @brief PusherBoris_2D.h  generic class for the particle pusher of Boris in 2D, no z momentum.

 @date 2022-02-04
 */

#ifndef PUSHERBORIS_2D_H
#define PUSHERBORIS_2D_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris_D
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_2D : public Pusher
{
public:
    //! Creator for Pusher
    PusherBoris_2D( Params &params, Species *species );
    ~PusherBoris_2D();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
    
};

#endif

