/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERBORISWNSURF_H
#define PUSHERBORISWNSURF_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBorisWnSurfV
//  --------------------------------------------------------------------------------------------------------------------
class PusherBorisWnSurf : public Pusher
{
public:
    //! Creator for Pusher
    PusherBorisWnSurf( Params &params, Species *species );
    ~PusherBorisWnSurf();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
    
};

#endif
