#include "PusherBorisWnSurf.h"

#include <iostream>
#include <cmath>
#ifdef _GPU
#include <accelmath.h>
#endif

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherBorisWnSurf::PusherBorisWnSurf( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherBorisWnSurf::~PusherBorisWnSurf()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBorisWnSurf::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset , double time)
{
    
    std::vector<double> * Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> * Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> * nSurfpart = &( smpi->dynamics_nSurfpart[ithread] );    
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
    double * __restrict__ position_x_old = particles.getPtrPosition_old(0);
    double * __restrict__ position_y = NULL;
    double * __restrict__ position_y_old = NULL;
    double * __restrict__ position_z = NULL;
    double * __restrict__ position_z_old = NULL;
    if (nDim_>1) {
        position_y = particles.getPtrPosition(1);
        position_y_old = particles.getPtrPosition_old(1);
        if (nDim_>2) {
            position_z = particles.getPtrPosition(2);
            position_z_old = particles.getPtrPosition_old(2);
        }
    }
    double * __restrict__ momentum_x = particles.getPtrMomentum(0);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);

    short * __restrict__ charge = particles.getPtrCharge();

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
    double * __restrict__ nSurfx = &( ( *nSurfpart )[0*nparts] );
    double * __restrict__ nSurfy = &( ( *nSurfpart )[1*nparts] );
    double * __restrict__ nSurfz = &( ( *nSurfpart )[2*nparts] );
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        
        int ipart2 = ipart - ipart_buffer_offset;
        bool reflected = false;
        charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;
        double pn = momentum_x[ipart]*nSurfx[ipart2]+momentum_y[ipart]*nSurfy[ipart2]+momentum_z[ipart]*nSurfz[ipart2];
        double nSurf2 = (nSurfx[ipart2])*(nSurfx[ipart2])+(nSurfy[ipart2])*(nSurfy[ipart2])+(nSurfz[ipart2])*(nSurfz[ipart2]);
        // reflect off surface if particle is outside and dot(n,p)>0. Momentum reflected specularly, position kept
        if(nSurf2>0 && pn>0) {
            reflected = true;
            momentum_x[ipart] -= 2*nSurfx[ipart2]*pn/nSurf2;
            momentum_y[ipart] -= 2*nSurfy[ipart2]*pn/nSurf2;
            momentum_y[ipart] -= 2*nSurfz[ipart2]*pn/nSurf2;
            if(nDim_>2){
                // Assuming point of reflection is at the middle of positon - position_old (valid for distance travel within 1 timestep << cell size)
                // May move this part to interpolation
                double xn = (position_x[ipart]-position_x_old[ipart])*nSurfx[ipart2]+(position_y[ipart]-position_y_old[ipart])*nSurfy[ipart2]+(position_z[ipart]-position_z_old[ipart])*nSurfz[ipart2];
                position_x[ipart] -= nSurfx[ipart2]*xn/nSurf2;
                position_y[ipart] -= nSurfy[ipart2]*xn/nSurf2 ;
                position_z[ipart] -= nSurfz[ipart2]*xn/nSurf2 ;
            } else if (nDim_>1){
                double xn = (position_x[ipart]-position_x_old[ipart])*nSurfx[ipart2]+(position_y[ipart]-position_y_old[ipart])*nSurfy[ipart2];
                position_x[ipart] -= nSurfx[ipart2]*xn/nSurf2;
                position_y[ipart] -= nSurfy[ipart2]*xn/nSurf2 ;
            } else {
                double xn = (position_x[ipart]-position_x_old[ipart])*nSurfx[ipart2];
                position_x[ipart] -= nSurfx[ipart2]*xn/nSurf2;

            }            
            // here take the assumption that the distance move in  1 timestep is small, change in interpolated field value at surface reflection is negligible
            // MESSAGE( "reflected");
        }

        // init Half-acceleration in the electric field
        pxsm = charge_over_mass_dts2*( Ex[ipart2] );
        pysm = charge_over_mass_dts2*( Ey[ipart2] );
        pzsm = charge_over_mass_dts2*( Ez[ipart2] );

        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
        umx = momentum_x[ipart] + pxsm;
        umy = momentum_y[ipart] + pysm;
        umz = momentum_z[ipart] + pzsm;

        // Rotation in the magnetic field
        local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
        Tx    = local_invgf * ( Bx[ipart2] );
        Ty    = local_invgf * ( By[ipart2] );
        Tz    = local_invgf * ( Bz[ipart2] );
        inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );

        pxsm += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* umx  +      2.0*( Tx*Ty+Tz )* umy  +      2.0*( Tz*Tx-Ty )* umz )*inv_det_T;
        pysm += ( 2.0*( Tx*Ty-Tz )* umx  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* umy  +      2.0*( Ty*Tz+Tx )* umz )*inv_det_T;
        pzsm += ( 2.0*( Tz*Tx+Ty )* umx  +      2.0*( Ty*Tz-Tx )* umy  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* umz )*inv_det_T;

        // finalize Half-acceleration in the electric field
        local_invgf = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        invgf[ipart2] = local_invgf; //1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );

        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = pzsm;

        // Move the particle
        local_invgf *= dt;
        //position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart2];
        position_x[ipart] += pxsm*local_invgf;
        if (nDim_>1) {
            position_y[ipart] += pysm*local_invgf;
            if (nDim_>2) {
                position_z[ipart] += pzsm*local_invgf;
            }
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
