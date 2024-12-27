/*! @file PusherBunch.h

 @brief PusherBunch.h  generic class for the particle pusher of Boris.

 @date 2013-02-15
 */

#ifndef PUSHERBUNCH_H
#define PUSHERBUNCH_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBunch
//  --------------------------------------------------------------------------------------------------------------------
class PusherBunch : public Pusher
{
public:
    //! Creator for Pusher
    PusherBunch( Params &params, Species *species );
    ~PusherBunch();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref = 0 );
    
};

#endif

