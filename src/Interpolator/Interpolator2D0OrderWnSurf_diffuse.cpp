#include "Interpolator2D0OrderWnSurf_diffuse.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "LaserEnvelope.h"
#include "userFunctions.h"
using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D0OrderWnSurf_diffuse
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D0OrderWnSurf_diffuse::Interpolator2D0OrderWnSurf_diffuse( Params &params, Patch *patch ) : Interpolator2D( params, patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D0OrderWnSurf_diffuse::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc, double *nSurfLoc )
{

    // Static cast of the electromagnetic fields
    Field2D *nSurfx2D = static_cast<Field2D *>( EMfields->nSurfx_ );
    Field2D *nSurfy2D = static_cast<Field2D *>( EMfields->nSurfy_ );
    Field2D *nSurfz2D = static_cast<Field2D *>( EMfields->nSurfz_ );

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_x_old = particles.getPtrPosition_old(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ position_y_old = particles.getPtrPosition_old(1);

    double * __restrict__ momentum_x = particles.getPtrMomentum(0);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Normalized particle position
    double xpn = position_x[ipart]*d_inv_[0];
    double ypn = position_y[ipart]*d_inv_[1];
    // Calculate coeffs
    coeffs( xpn, ypn );

    // Interpolation of Ex^(d,p)
    double nSurfLocx = compute( &coeffxd_[1], &coeffyp_[1], nSurfx2D, ip_, jp_ );
    double nSurfLocy = compute( &coeffxp_[1], &coeffyd_[1], nSurfy2D, ip_, jp_ );
    double nSurfLocz = compute( &coeffxp_[1], &coeffyp_[1], nSurfz2D, ip_, jp_ );
    // Decide whether the particle is reflected
    double pn = momentum_x[ipart]*nSurfLocx+momentum_y[ipart]*nSurfLocy+momentum_z[ipart]*nSurfLocz;
    double nSurf2 = (nSurfLocx)*(nSurfLocx)+(nSurfLocy)*(nSurfLocy)+(nSurfLocz)*(nSurfLocz);
    bool reflected = false;
    if(nSurf2>0 && pn>0) {
        double prob = rand_->uniform();
        if(prob > specularity) {
            reflected = true;
            double u = rand_->uniform(); // CDF(x) = 1/2(sinx+1) = u uniformly [0,1]
            double angle = asin(sqrt(u)); // asin [-pi/2, pi/2]
            double mometum_magnitude = sqrt(momentum_x[ipart]*momentum_x[ipart]+momentum_y[ipart]*momentum_y[ipart]);
            momentum_x[ipart] = (nSurfLocx * cos(angle) - nSurfLocy * sin(angle))*(-1); // rotate around (-1)*nsurf
            momentum_y[ipart] = (nSurfLocx * sin(angle) + nSurfLocy * cos(angle))*(-1);
            double mometum_magnitude2 = sqrt(momentum_x[ipart]*momentum_x[ipart]+momentum_y[ipart]*momentum_y[ipart]);
            momentum_x[ipart] = momentum_x[ipart]/mometum_magnitude2*mometum_magnitude; // re-normalize momentum magnitude
            momentum_y[ipart] = momentum_y[ipart]/mometum_magnitude2*mometum_magnitude;
            
            // Assuming point of reflection is at the middle of positon - position_old (valid for distance travel within 1 timestep << cell size)
            double xn = (position_x[ipart]-position_x_old[ipart])*nSurfLocx+(position_y[ipart]-position_y_old[ipart])*nSurfLocy;
            position_x[ipart] -= nSurfLocx*xn/nSurf2;
            position_y[ipart] -= nSurfLocy*xn/nSurf2 ;
        } else{
            reflected = true;
            momentum_x[ipart] -= 2*nSurfLocx*pn/nSurf2;
            momentum_y[ipart] -= 2*nSurfLocy*pn/nSurf2;
            momentum_z[ipart] -= 2*nSurfLocz*pn/nSurf2;
            // Assuming point of reflection is at the middle of positon - position_old (valid for distance travel within 1 timestep << cell size)
            double xn = (position_x[ipart]-position_x_old[ipart])*nSurfLocx+(position_y[ipart]-position_y_old[ipart])*nSurfLocy;
                position_x[ipart] -= nSurfLocx*xn/nSurf2;
                position_y[ipart] -= nSurfLocy*xn/nSurf2 ;
        }
        
        // Interpolation of nSurf at new location if is used in subsequent moves, likely not
        /*
        
            double xpn = position_x[ipart]*d_inv_[0];
            double ypn = position_y[ipart]*d_inv_[1];
            // Calculate coeffs
            coeffs( xpn, ypn );
            double nSurfLocx = compute( &coeffxd_[1], &coeffyp_[1], nSurfx2D, ip_, jp_ );
            double nSurfLocy = compute( &coeffxp_[1], &coeffyd_[1], nSurfy2D, ip_, jp_ );
            double nSurfLocz = compute( &coeffxp_[1], &coeffyp_[1], nSurfz2D, ip_, jp_ );
        */
        // here take the assumption that the distance move in  1 timestep is small, change in interpolated field value at surface reflection is negligible
        // MESSAGE( "reflected");
    }
    *( nSurfLoc+0*nparts ) = nSurfLocx;
    // Interpolation of Ey^(p,d)
    *( nSurfLoc+1*nparts ) = nSurfLocy;
    // Interpolation of Ez^(p,p)
    *( nSurfLoc+2*nparts ) = nSurfLocz;
 
} // END Interpolator2D0OrderWnSurf_diffuse

// -----------------------------------------------------------------------------
//
//! Interpolation of all fields and currents for a single particles
//! located at istart.
//! This version is not vectorized.
//! The input parameter iend not used for now, probes are interpolated one by one for now.
//
// -----------------------------------------------------------------------------
void Interpolator2D0OrderWnSurf_diffuse::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    /*
    
    
    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    // Interpolate E, B
    // Compute coefficient for ipart position
    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );
    Field2D *Jx2D = static_cast<Field2D *>( EMfields->Jx_ );
    Field2D *Jy2D = static_cast<Field2D *>( EMfields->Jy_ );
    Field2D *Jz2D = static_cast<Field2D *>( EMfields->Jz_ );
    Field2D *Rho2D= static_cast<Field2D *>( EMfields->rho_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    // Calculate coeffs
    coeffs( xpn, ypn );

    int nparts( particles.size() );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) = compute( &coeffxd_[1], &coeffyp_[1], Ex2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp_[1], &coeffyd_[1], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp_[1], &coeffyp_[1], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[1], &coeffyd_[1], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd_[1], &coeffyp_[1], By2D, id_, jp_ );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd_[1], &coeffyd_[1], Bz2D, id_, jd_ );
    // Interpolation of Jx^(d,p)
    JLoc->x = compute( &coeffxd_[1], &coeffyp_[1], Jx2D, id_, jp_ );
    // Interpolation of Jy^(p,d)
    JLoc->y = compute( &coeffxp_[1], &coeffyd_[1], Jy2D, ip_, jd_ );
    // Interpolation of Jz^(p,p)
    JLoc->z = compute( &coeffxp_[1], &coeffyp_[1], Jz2D, ip_, jp_ );
    // Interpolation of Rho^(p,p)
    ( *RhoLoc ) = compute( &coeffxp_[1], &coeffyp_[1], Rho2D, ip_, jp_ );
    */
}

//! Interpolator on another field than the basic ones
void Interpolator2D0OrderWnSurf_diffuse::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    /*
    
    Field2D *F = static_cast<Field2D *>( *field );
    double *coeffx = F->isDual( 0 ) ? &coeffxd_[1] : &coeffxp_[1];
    double *coeffy = F->isDual( 1 ) ? &coeffyd_[1] : &coeffyp_[1];
    int *i = F->isDual( 0 ) ? &id_ : &ip_;
    int *j = F->isDual( 1 ) ? &jd_ : &jp_;

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        coeffs( xpn, ypn );
        FieldLoc[ipart] = compute( coeffx, coeffy, F, *i, *j );
    }
    */
}

// -----------------------------------------------------------------------------
//! Wrapper called by the particle dynamics section
// -----------------------------------------------------------------------------
void Interpolator2D0OrderWnSurf_diffuse::fieldsWrapper(   ElectroMagn *EMfields,
                                            Particles &particles,
                                            SmileiMPI *smpi,
                                            int *istart,
                                            int *iend,
                                            int ithread,
                                            int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *nSurfpart = &( smpi->dynamics_nSurfpart[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] , &( *nSurfpart )[ipart] );
        //Buffering of iol and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;
    }

}

// -----------------------------------------------------------------------------
//! Interpolator specific to tracked particles.
//! A selection of particles may be provided
// -----------------------------------------------------------------------------
void Interpolator2D0OrderWnSurf_diffuse::fieldsSelection( ElectroMagn *EMfields,
                                            Particles &particles,
                                            double *buffer,
                                            int offset,
                                            vector<unsigned int> *selection )
{
/*
    if( selection ) {

        int nsel_tot = selection->size();
        for( int isel=0 ; isel<nsel_tot; isel++ ) {
            fields( EMfields, particles, ( *selection )[isel], offset, buffer+isel, buffer+isel+3*offset );
        }

    } else {

        int npart_tot = particles.size();
        for( int ipart=0 ; ipart<npart_tot; ipart++ ) {
            fields( EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset );
        }

    }*/
}

// -----------------------------------------------------------------------------
//! Interpolator specific to the envelope model
// -----------------------------------------------------------------------------
void Interpolator2D0OrderWnSurf_diffuse::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
   /* // Static cast of the envelope fields
    Field2D *Phi2D = static_cast<Field2D *>( EMfields->envelope->Phi_ );
    Field2D *GradPhix2D = static_cast<Field2D *>( EMfields->envelope->GradPhix_ );
    Field2D *GradPhiy2D = static_cast<Field2D *>( EMfields->envelope->GradPhiy_ );
    Field2D *GradPhiz2D = static_cast<Field2D *>( EMfields->envelope->GradPhiz_ );

    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );


        // -------------------------
        // Interpolation of Phi^(p,p)
        // -------------------------
        ( *PHIpart )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], Phi2D, ip_, jp_ );

        // -------------------------
        // Interpolation of GradPhix^(p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhix2D, ip_, jp_ );

        // -------------------------
        // Interpolation of GradPhiy^(p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiy2D, ip_, jp_ );

        // -------------------------
        // Interpolation of GradPhiz^(p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiz2D, ip_, jp_ );


        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;

    }
*/

} // END Interpolator2D0OrderWnSurf_diffuse


void Interpolator2D0OrderWnSurf_diffuse::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    /*
    // Static cast of the envelope fields
    Field2D *Phi_m2D = static_cast<Field2D *>( EMfields->envelope->Phi_m );

    Field2D *GradPhix_m2D = static_cast<Field2D *>( EMfields->envelope->GradPhix_m );
    Field2D *GradPhiy_m2D = static_cast<Field2D *>( EMfields->envelope->GradPhiy_m );
    Field2D *GradPhiz_m2D = static_cast<Field2D *>( EMfields->envelope->GradPhiz_m );

    std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];


        // Indexes of the central nodes
        ip_ = round( xpn );
        jp_ = round( ypn );

        // Declaration and calculation of the coefficient for interpolation
        double delta2;


        deltax   = xpn - ( double )ip_;
        delta2  = deltax*deltax;
        coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );


        deltay   = ypn - ( double )jp_;
        delta2  = deltay*deltay;
        coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );



        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;




        // -------------------------
        // Interpolation of Phiold^(p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], Phi_m2D, ip_, jp_ );

        // -------------------------
        // Interpolation of GradPhixold^(p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhix_m2D, ip_, jp_ );

        // -------------------------
        // Interpolation of GradPhiyold^(p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiy_m2D, ip_, jp_ );

        // -------------------------
        // Interpolation of GradPhizold^(p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiz_m2D, ip_, jp_ );

        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;


    }*/

} // END Interpolator2D0OrderWnSurf_diffuse


void Interpolator2D0OrderWnSurf_diffuse::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    /*
    // Static cast of the electromagnetic fields
    Field2D *Env_A_abs_2D = static_cast<Field2D *>( EMfields->Env_A_abs_ );
    Field2D *Env_Chi_2D = static_cast<Field2D *>( EMfields->Env_Chi_ );
    Field2D *Env_E_abs_2D = static_cast<Field2D *>( EMfields->Env_E_abs_ );
    Field2D *Env_Ex_abs_2D = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];



    // Indexes of the central nodes
    ip_ = round( xpn );
    jp_ = round( ypn );



    // Declaration and calculation of the coefficient for interpolation
    double delta2;


    deltax   = xpn - ( double )ip_;
    delta2  = deltax*deltax;
    coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );

    deltay   = ypn - ( double )jp_;
    delta2  = deltay*deltay;
    coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );



    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;


    // -------------------------
    // Interpolation of Env_A_abs_^(p,p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_A_abs_2D, ip_, jp_ );

    // -------------------------
    // Interpolation of Env_Chi_^(p,p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_Chi_2D, ip_, jp_ );

    // -------------------------
    // Interpolation of Env_E_abs_^(p,p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_E_abs_2D, ip_, jp_ );

    // -------------------------
    // Interpolation of Env_Ex_abs_^(p,p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_Ex_abs_2D, ip_, jp_ );
*/
} // END Interpolator2D0OrderWnSurf_diffuse

void Interpolator2D0OrderWnSurf_diffuse::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    /*
    // Static cast of the envelope fields
    Field2D *EnvEabs  = static_cast<Field2D *>( EMfields->Env_E_abs_ );
    Field2D *EnvExabs = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );

    std::vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
    std::vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );

    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];

        // Indexes of the central nodes
        ip_ = round( xpn );
        jp_ = round( ypn );

        // Declaration and calculation of the coefficient for interpolation
        double delta2;

        deltax   = xpn - ( double )ip_;
        delta2  = deltax*deltax;
        coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );

        deltay   = ypn - ( double )jp_;
        delta2  = deltay*deltay;
        coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );

        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;

        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], EnvEabs, ip_, jp_ );

        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], EnvExabs, ip_, jp_ );


    }
*/

} // END Interpolator2D0OrderWnSurf_diffuse
