/*! @file PusherPrescribed.h

 @brief PusherPrescribed.h  generic class for the particle pusher of Boris.

 @date 2013-02-15
 */

#ifndef PUSHERPRESCRIBED_H
#define PUSHERPRESCRIBED_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherPrescribed
//  --------------------------------------------------------------------------------------------------------------------
class PusherPrescribed : public Pusher
{
public:
    //! Creator for Pusher
    PusherPrescribed( Params &params, Species *species );
    ~PusherPrescribed();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time =0);
    
};

#endif

