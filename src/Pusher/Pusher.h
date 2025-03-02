#ifndef PUSHER_H
#define PUSHER_H

#include "Params.h"
#include "Field.h"
#include "Random.h"
class Particles;


//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class Pusher
{

public:
    //! Creator for Pusher
    Pusher( Params &params, Species *species );
    virtual ~Pusher();
    
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0, double time=0 ) = 0;
    
protected:
    double dt, dts2, dts4;
    //! \todo Move mass_ in Particles_
    // mass_ relative to Species but used in the particle pusher
    double mass_;
    double mfp_;
    double prescribed_freq;
    double prescribed_amp;
    double V_const;
    double one_over_mass_;
    double time_frozen_;
    Profile *prescribed_profile;
    int nDim_;
    //! min_loc_vec (copy from picparams)
    std::vector<double> min_loc_vec;
    double dx_inv_[3];
    unsigned int nspace[3];
    bool vecto;
    Random *rand_pusher;
};//END class

#endif
