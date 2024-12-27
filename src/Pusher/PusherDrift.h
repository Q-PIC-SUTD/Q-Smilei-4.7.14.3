/*! @file PusherDrift.h

 @brief PusherDrift.h  generic class for the particle pusher of Boris.

 @date 2013-02-15
 */

#ifndef PUSHERDRIFT_H
#define PUSHERDRIFT_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherDrift
//  --------------------------------------------------------------------------------------------------------------------
class PusherDrift : public Pusher
{
public:
    //! Creator for Pusher
    PusherDrift( Params &params, Species *species );
    ~PusherDrift();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
    
};

#endif

