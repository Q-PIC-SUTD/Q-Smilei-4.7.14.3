// --------------------------------------------------------------------------------------------------------------------
//
//! \file PusherFactory.h
//
//! \brief Class PusherFactory that manage the pusher association to the species.
//
// --------------------------------------------------------------------------------------------------------------------

#ifndef PUSHERFACTORY_H
#define PUSHERFACTORY_H

#include "Pusher.h"
#include "PusherBoris.h"
#include "PusherPonderomotiveBoris.h"
#include "PusherPonderomotivePositionBoris.h"
#include "PusherVay.h"
#include "PusherBorisNR.h"
#include "PusherHigueraCary.h"
#include "PusherPhoton.h"
#include "PusherBoris_mfp.h"
#include "PusherBorisWnSurf.h"
#include "PusherBoris_mfp_normal.h"
#include "PusherBoris_mfp_rejection.h"
#include "PusherBoris_mfp_2D.h"
#include "PusherBoris_2D.h"
#include "PusherDrift.h"
#include "PusherPrescribed.h"
#include "PusherPrescribedProfile.h"
#include "PusherMassless.h"
#include "PusherMassless_mfp_rejection.h"
#include "PusherMasslessQuasistatic.h"
#include "PusherBunch.h"
// #ifdef _VECTO
// #include "PusherPonderomotiveBorisV.h"
// #include "PusherPonderomotivePositionBorisV.h"
// #endif

#include "Params.h"
#include "Species.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherFactory
//
//! \brief Class PusherFactory that manage the pusher association to the species.
//  --------------------------------------------------------------------------------------------------------------------
class PusherFactory
{
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate pusher for the species ispec
    //! \param ispec SpeciesId
    //! \param params Parameters
    //  --------------------------------------------------------------------------------------------------------------------
    static Pusher *create( Params &params, Species *species )
    {
        Pusher *Push = NULL;
        
        // Particle of matter
        if( species->mass_ > 0 ) {
            // assign the correct Pusher to Push
            // Pusher of Boris
            if( species->pusher_name_ == "boris" ) {
                /* if(params.with_nSurf) {
                    Push = new PusherBorisWnSurf( params, species );
                    // MESSAGE( "Apply nSurf for species" << species->name_ );
                } else {
                    Push = new PusherBoris( params, species );
                }*/
                    Push = new PusherBoris( params, species );
            } else if( species->pusher_name_ == "ponderomotive_boris" ) {
            
                int n_envlaser = params.Laser_Envelope_model;
                if( n_envlaser <1 ) {
                    ERROR_NAMELIST( "No Laser Envelope present." 
                                <<  " The pusher ponderomotive_boris can be used only in presence of a Laser Envelope.",
                    LINK_NAMELIST + std::string("#pusher"));
                }
                
                Push = new PusherPonderomotiveBoris( params, species );
            // Non-relativistic Boris pusher
            } else if( species->pusher_name_ == "borisnr" ) {
                Push = new PusherBorisNR( params, species );
            }
            // Pusher of J.L. Vay
            else if( species->pusher_name_ == "vay" ) {
                Push = new PusherVay( params, species );
            // Pusher of Higuera Cary
            } else if( species->pusher_name_ == "higueracary" ) {
                Push = new PusherHigueraCary( params, species );
            } else if( species->pusher_name_ == "boris_mfp" ) {               
                if( species->mfp_ > 0.) {
                    Push = new PusherBoris_mfp( params, species );                   
                }               
                else {
                    ERROR( "Invalid mean free path" );
                    MESSAGE( "species:" << species->name_ << " - mfp" << species->mfp_ );
                }
            } else if( species->pusher_name_ == "boris_mfp_normal" ) {               
                if( species->mfp_ > 0.) {
                    Push = new PusherBoris_mfp_normal( params, species );                   
                }               
                else {
                    ERROR( "Invalid mean free path" );
                    MESSAGE( "species:" << species->name_ << " - mfp" << species->mfp_ );
                }
            } else if( species->pusher_name_ == "boris_mfp_rejection" ) {               
                if( species->mfp_ > 0.) {
                    Push = new PusherBoris_mfp_rejection( params, species );                   
                }               
                else {
                    ERROR( "Invalid mean free path" );
                    MESSAGE( "species:" << species->name_ << " - mfp" << species->mfp_ );
                }
            } else if(species->pusher_name_ == "boris_2D") {
                Push = new PusherBoris_2D( params, species );
            } else if( species->pusher_name_ == "boris_mfp_2D" ) {
                if( species->mfp_ > 0.) {
                    Push = new PusherBoris_mfp_2D( params, species );
                }
                else {
                    ERROR( "Invalid mean free path" );
                    MESSAGE( "species:" << species->name_ << " - mfp" << species->mfp_ );
                }
            } else if( species->pusher_name_ == "drift" ) {
                Push = new PusherDrift( params, species );
            } else if(species->pusher_name_ == "prescribed") {
                Push = new PusherPrescribed( params, species );  
            } else if(species->pusher_name_ == "prescribed_profile") {
                Push = new PusherPrescribedProfile( params, species );  
            } else if(species->pusher_name_ == "massless") {
                if(params.quasistatic){
                    Push = new PusherMasslessQuasistatic( params, species );   
                } else {
                    Push = new PusherMassless( params, species );   
                }
                       
            } else if(species->pusher_name_ == "massless_mfp_rejection") {
                
                    Push = new PusherMassless_mfp_rejection( params, species );   
                
                       
            } else {
                ERROR_NAMELIST( "For species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`",
                    LINK_NAMELIST + std::string("#pusher") );
            }
        }
        // Photon
        else if( species->mass_ == 0 ) {
            if( species->pusher_name_ == "norm" ) {
                Push = new PusherPhoton( params, species );
 
            } else {
                ERROR_NAMELIST( "For photon species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`",
                   LINK_NAMELIST + std::string("#pusher") );
            }
        }
        
        if( params.Laser_Envelope_model ) {
            if( species->pusher_name_ != "ponderomotive_boris" ) {
                ERROR_NAMELIST( "For species " << species->name_ 
                << " the only pusher available to interact with the envelope is ponderomotive_boris",
                LINK_NAMELIST + std::string("#pusher") );
            }
        }
        return Push;
    }
    
    static Pusher *create_ponderomotive_position_updater( Params &params, Species *species )
    {
        Pusher *Push_ponderomotive_position = NULL;
        
        // Particle of matter
        if( species->mass_ > 0 ) {
            // assign the correct Pusher to Push_ponderomotive_position
            if( species->pusher_name_ == "ponderomotive_boris" ) {
                    Push_ponderomotive_position = new PusherPonderomotivePositionBoris( params, species );
            }
            
            else {
                ERROR_NAMELIST( "For species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`",
                   LINK_NAMELIST + std::string("#pusher") );
            }
        } else {
            ERROR_NAMELIST( "Ponderomotive pusher is not a valid choice for photons",
            LINK_NAMELIST + std::string("#pusher"));
        }
        return Push_ponderomotive_position;
    }
    
};

#endif
