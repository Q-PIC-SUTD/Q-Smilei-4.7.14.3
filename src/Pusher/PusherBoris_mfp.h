/*! @file PusherBoris_mfp.h

 @brief PusherBoris_mfp.h  generic class for the particle pusher of Boris.

 @date 2013-02-15
 */

#ifndef PUSHERBORIS_MFP_H
#define PUSHERBORIS_MFP_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris_mfp
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_mfp : public Pusher
{
public:
    //! Creator for Pusher
    PusherBoris_mfp( Params &params, Species *species );
    ~PusherBoris_mfp();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
   
};

#endif

