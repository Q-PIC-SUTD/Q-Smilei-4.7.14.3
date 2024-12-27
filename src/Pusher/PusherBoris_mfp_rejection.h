/*! @file PusherBoris_mfp_rejection.h

 @brief PusherBoris_mfp_rejection.h  generic class for the particle pusher of Boris.
https://mathworld.wolfram.com/SpherePointPicking.html
Muller, M. E. "A Note on a Method for Generating Points Uniformly on N-Dimensional Spheres." Comm. Assoc. Comput. Mach. 2, 19-20, Apr. 1959.
Marsaglia, G. "Choosing a Point from the Surface of a Sphere." Ann. Math. Stat. 43, 645-646, 1972.
 @date 2013-02-15
 */

#ifndef PUSHERBORIS_MFP_REJECTION_H
#define PUSHERBORIS_MFP_REJECTION_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris_mfp
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_mfp_rejection : public Pusher
{
public:
    //! Creator for Pusher
    PusherBoris_mfp_rejection( Params &params, Species *species );
    ~PusherBoris_mfp_rejection();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 );
   
};

#endif

