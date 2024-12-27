#include "PusherMasslessQuasistatic.h"

#include <iostream>
#include <cmath>
#ifdef _GPU
#include <accelmath.h>
#endif

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherMasslessQuasistatic::PusherMasslessQuasistatic( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherMasslessQuasistatic::~PusherMasslessQuasistatic()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (MasslessQuasistatic) scheme
***********************************************************************/

void PusherMasslessQuasistatic::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset , double time)
{

    std::vector<double> * Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> * Bpart = &( smpi->dynamics_Bpart[ithread] );
    double * __restrict__ invgf = &( smpi->dynamics_invgf[ithread][0] );

    double pxsm, pysm, pzsm;
    double local_invgf;
    double umx, umy, umz; //upx, upy, upz;
    double charge_over_mass_dts;
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

    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        int ipart2 = ipart - ipart_buffer_offset;

        charge_over_mass_dts = ( double )( charge[ipart] )*dt;

        // acceleration in the electric field
        pxsm = charge_over_mass_dts*( Ex[ipart2] );
        pysm = charge_over_mass_dts*( Ey[ipart2] );
        pzsm = charge_over_mass_dts*( Ez[ipart2] );

        momentum_x[ipart] += pxsm;
        momentum_y[ipart] += pysm;
        momentum_z[ipart] += 0. ;//pzsm;

        double momentum_mag=sqrt(momentum_x[ipart]*momentum_x[ipart]+momentum_y[ipart]*momentum_y[ipart]+momentum_z[ipart]*momentum_z[ipart]);

        //position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart2];
        position_x[ipart] += V_const*momentum_x[ipart]/momentum_mag*dt;
        if (nDim_>1) {
            position_y[ipart] += V_const*momentum_y[ipart]/momentum_mag*dt;
            if (nDim_>2) {
                position_z[ipart] += V_const*momentum_z[ipart]/momentum_mag*dt;
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
