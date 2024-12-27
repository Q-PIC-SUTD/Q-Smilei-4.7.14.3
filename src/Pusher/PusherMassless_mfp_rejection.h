/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERMASSLESS_mfp_rejection_H
#define PUSHERMASSLESS_mfp_rejection_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherMassless_mfp_rejectionV
//  --------------------------------------------------------------------------------------------------------------------
class PusherMassless_mfp_rejection : public Pusher
{
public:
    //! Creator for Pusher
    PusherMassless_mfp_rejection( Params &params, Species *species );
    ~PusherMassless_mfp_rejection();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
    
};

#endif
