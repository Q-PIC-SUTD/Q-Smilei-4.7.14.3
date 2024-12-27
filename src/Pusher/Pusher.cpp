#include "Pusher.h"
#include "Params.h"
#include "Species.h"
#include "Random.h"
Pusher::Pusher( Params &params, Species *species ) :
    min_loc_vec( species->min_loc_vec ),
    vecto( params.vectorization_mode=="on" || params.vectorization_mode=="adaptive_mixed_sort" || params.vectorization_mode=="adaptive" || params.cell_sorting )
{
    for( unsigned int ipos=0; ipos < params.nDim_particle ; ipos++ ) {
        dx_inv_[ipos] = species->dx_inv_[ipos];
    }
    
    rand_pusher = new Random(params.random_seed);

    nspace[0] = 0;
    nspace[1] = params.n_space[1]+1;
    nspace[2] = params.n_space[2]+1;
    mfp_ = species->mfp_;
    if(species->pusher_name_=="prescribed"){
        prescribed_freq = species->prescribed_freq;
        prescribed_amp = species->prescribed_amp;

    }
    if(species->pusher_name_=="prescribed_profile"){
        prescribed_freq = species->prescribed_freq;
        prescribed_profile = species->prescribed_profile;
        prescribed_amp = species->prescribed_amp;
    }
    //if(species->pusher_name_=="massless"){
    V_const = species->V_const;

    //}    
    mass_          = species->mass_;
    if( mass_ > 0. ) {
        one_over_mass_ = 1.0/mass_;
    } else {
        one_over_mass_ = 0.;
    }
    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
    
    nDim_          = params.nDim_particle;
    time_frozen_ = species->time_frozen_;
    
    
}

Pusher::~Pusher()
{
}
