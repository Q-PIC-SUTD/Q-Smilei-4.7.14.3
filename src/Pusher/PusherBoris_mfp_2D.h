/*! @file PusherBoris_mfp_2D.h

 @brief PusherBoris_mfp_2D.h  generic class for the particle pusher of Boris with mfp in 2D.

 @date 2022-01-26
 */

#ifndef PUSHERBORIS_MFP_2D_H
#define PUSHERBORIS_MFP_2D_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris_mfp_2D
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_mfp_2D : public Pusher
{
public:
    //! Creator for Pusher
    PusherBoris_mfp_2D( Params &params, Species *species );
    ~PusherBoris_mfp_2D();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );;
    
};

#endif

