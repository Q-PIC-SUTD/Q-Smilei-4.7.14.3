/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERMASSLESS_H
#define PUSHERMASSLESS_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherMasslessV
//  --------------------------------------------------------------------------------------------------------------------
class PusherMassless : public Pusher
{
public:
    //! Creator for Pusher
    PusherMassless( Params &params, Species *species );
    ~PusherMassless();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
    
};

#endif
