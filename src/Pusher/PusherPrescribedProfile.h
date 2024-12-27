/*! @file PusherPrescribedProfile.h

 @brief PusherPrescribedProfile.h  generic class for the particle pusher of Boris.

 @date 2013-02-15
 */

#ifndef PUSHERPRESCRIBEDPROFILE_H
#define PUSHERPRESCRIBEDPROFILE_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherPrescribedProfile
//  --------------------------------------------------------------------------------------------------------------------
class PusherPrescribedProfile : public Pusher
{
public:
    //! Creator for Pusher
    PusherPrescribedProfile( Params &params, Species *species );
    ~PusherPrescribedProfile();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time =0);
    
};

#endif

