/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERMASSLESSQUASISTATIC_H
#define PUSHERMASSLESSQUASISTATIC_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherMasslessQuasistaticV
//  --------------------------------------------------------------------------------------------------------------------
class PusherMasslessQuasistatic : public Pusher
{
public:
    //! Creator for Pusher
    PusherMasslessQuasistatic( Params &params, Species *species );
    ~PusherMasslessQuasistatic();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
    
};

#endif
