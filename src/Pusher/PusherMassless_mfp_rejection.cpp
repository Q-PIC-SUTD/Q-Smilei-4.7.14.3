#include "PusherMassless_mfp_rejection.h"

#include <iostream>
#include <cmath>
#ifdef _GPU
#include <accelmath.h>
#endif

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherMassless_mfp_rejection::PusherMassless_mfp_rejection( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherMassless_mfp_rejection::~PusherMassless_mfp_rejection()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Massless) scheme
***********************************************************************/

void PusherMassless_mfp_rejection::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset , double time)

{
    std::vector<double> * Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> * Bpart = &( smpi->dynamics_Bpart[ithread] );
    double * __restrict__ invgf = &( smpi->dynamics_invgf[ithread][0] );

    double pxsm, pysm, pzsm;
    double local_invgf;
    double umx, umy, umz; //upx, upy, upz;
    double charge_over_mass_dts2;
    double inv_det_T, Tx, Ty, Tz;
    //double Tx2, Ty2, Tz2;
    //double TxTy, TyTz, TzTx;
    //double alpha;

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = NULL;
    double * __restrict__ position_z = NULL;
    if (nDim_>1) {
        position_y = particles.getPtrPosition(1);
        if (nDim_>2) {
            position_z = particles.getPtrPosition(2);
        }
    }
    double * __restrict__ momentum_x = particles.getPtrMomentum(0);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);

    short * __restrict__ charge = particles.getPtrCharge();
    double * __restrict__ acc_path = particles.getPtrAccPath();
    int nparts;
    if (vecto) {
        nparts = Epart->size()/3;
    } else {
        nparts = particles.size();
    }
    double * __restrict__ Ex = &( ( *Epart )[0*nparts] );
    double * __restrict__ Ey = &( ( *Epart )[1*nparts] );
    double * __restrict__ Ez = &( ( *Epart )[2*nparts] );
    double * __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    double * __restrict__ By = &( ( *Bpart )[1*nparts] );
    double * __restrict__ Bz = &( ( *Bpart )[2*nparts] );
    //MESSAGE( 1, "Pusher Boris" );
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        int ipart2 = ipart - ipart_buffer_offset;


        // Decide whether to scatter the particle
        double mometum_magnitude; 
        mometum_magnitude = sqrt(momentum_x[ipart]*momentum_x[ipart]+momentum_y[ipart]*momentum_y[ipart]+momentum_z[ipart]*momentum_z[ipart]);    
        acc_path[ipart]  += dt;
        double collision_prob = dt/mfp_;
        double coin = rand_pusher->uniform();
        
        if(coin < collision_prob){
            double tempx = rand_pusher->uniform2();
            double tempy = rand_pusher->uniform2();
            // double tempz = rand_pusher->uniform2();
            double temp_mag = sqrt(tempx*tempx+tempy*tempy);
            while(temp_mag>1.-1.e-3){
                tempx = rand_pusher->uniform2();
                tempy = rand_pusher->uniform2();
                // tempz = rand_pusher->uniform2();
                temp_mag = sqrt(tempx*tempx+tempy*tempy);
            }

            momentum_x[ipart] = mometum_magnitude* tempx/temp_mag;
            momentum_y[ipart] = mometum_magnitude* tempy/temp_mag;
            // momentum_z[ipart] = mometum_magnitude* tempz/temp_mag;

            acc_path[ipart] = 0; 

        }
        /*
        if(coin < collision_prob){
            // MESSAGE( "deflection event after acc_path" << acc_path[ipart] << "collision_prob" << collision_prob << "coin" << coin);
            // double theta;
            double phi;
            //theta = ((float) rand()/RAND_MAX)*M_PI;
            phi = rand_pusher->uniform()*2*M_PI;
            momentum_x[ipart] = mometum_magnitude*sin(phi);
            momentum_y[ipart] = mometum_magnitude*cos(phi);
            //momentum[2][ipart] = mometum_magnitude*cos(theta);
            acc_path[ipart] = 0; 
        }*/
        // Change velocity in E and B field

        charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;

        // init Half-acceleration in the electric field
        pxsm = charge_over_mass_dts2*( Ex[ipart2] );
        pysm = charge_over_mass_dts2*( Ey[ipart2] );
        pzsm = 0; //charge_over_mass_dts2*( Ez[ipart2] );

        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
        umx = momentum_x[ipart] + pxsm;
        umy = momentum_y[ipart] + pysm;
        umz = 0; //momentum_z[ipart] + pzsm;

        // Rotation in the magnetic field
        local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
        Tx    = local_invgf * ( Bx[ipart2] );
        Ty    = local_invgf * ( By[ipart2] );
        Tz    = local_invgf * ( Bz[ipart2] );
        inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );

        pxsm += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* umx  +      2.0*( Tx*Ty+Tz )* umy  +      2.0*( Tz*Tx-Ty )* umz )*inv_det_T;
        pysm += ( 2.0*( Tx*Ty-Tz )* umx  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* umy  +      2.0*( Ty*Tz+Tx )* umz )*inv_det_T;
        pzsm += 0;//( 2.0*( Tz*Tx+Ty )* umx  +      2.0*( Ty*Tz-Tx )* umy  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* umz )*inv_det_T;

        // finalize Half-acceleration in the electric field
        local_invgf = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        invgf[ipart2] = local_invgf; //1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );

        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = 0;//pzsm;

        // Renormalize momentum to V_const
        
        mometum_magnitude = sqrt(momentum_x[ipart]*momentum_x[ipart]+momentum_y[ipart]*momentum_y[ipart]+momentum_z[ipart]*momentum_z[ipart]);
        momentum_x[ipart] = momentum_x[ipart]/mometum_magnitude*V_const;
        momentum_y[ipart] = momentum_y[ipart]/mometum_magnitude*V_const;
        // momentum_z[ipart] = 0; //momentum_z[ipart]/mometum_magnitude*V_const; // for 2D electrons in graphene
        // Move the particle
        local_invgf *= dt;
        //position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart2];
        position_x[ipart] += momentum_x[ipart]*local_invgf;
        if (nDim_>1) {
            position_y[ipart] += momentum_y[ipart]*local_invgf;
            /*if (nDim_>2) {
                position_z[ipart] += pzsm*local_invgf;
            }*/
        }
    }

    // if (nDim_>1) {
    //     #pragma omp simd
    //     for( int ipart=istart ; ipart<iend; ipart++ ) {
    //         position_y[ipart] += momentum_y[ipart]*invgf[ipart-ipart_buffer_offset]*dt;
    //     }
    // }
    // 
    // if (nDim_>2) {
    //     #pragma omp simd
    //     for( int ipart=istart ; ipart<iend; ipart++ ) {
    //         position_z[ipart] += momentum_z[ipart]*invgf[ipart-ipart_buffer_offset]*dt;
    //     }
    // }
}