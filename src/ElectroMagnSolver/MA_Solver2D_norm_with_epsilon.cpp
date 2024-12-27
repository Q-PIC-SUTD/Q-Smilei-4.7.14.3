
#include "MA_Solver2D_norm_with_epsilon.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MA_Solver2D_norm_with_epsilon::MA_Solver2D_norm_with_epsilon( Params &params )
    : Solver2D( params )
{
}

MA_Solver2D_norm_with_epsilon::~MA_Solver2D_norm_with_epsilon()
{
}

void MA_Solver2D_norm_with_epsilon::operator()( ElectroMagn *fields )
{

    // Static-cast of the fields
    Field2D *Ex2D = static_cast<Field2D *>( fields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( fields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    Field2D *Jx2D = static_cast<Field2D *>( fields->Jx_ );
    Field2D *Jy2D = static_cast<Field2D *>( fields->Jy_ );
    Field2D *Jz2D = static_cast<Field2D *>( fields->Jz_ );
    Field2D *epsilon2D = static_cast<Field2D *>( fields->epsilon_ );
    Field2D *epsilonx2D = static_cast<Field2D *>( fields->epsilonx_ );
    Field2D *epsilony2D = static_cast<Field2D *>( fields->epsilony_ );
    Field2D *epsilonz2D = static_cast<Field2D *>( fields->epsilonz_ );
    // Electric field Ex^(d,p)
    
    DEBUG( "nx_d:" << nx_d );
    DEBUG( "ny_p:" << ny_p );
    DEBUG( "nx_p:" << nx_p );
    DEBUG( "ny_d:" << ny_d );
    DEBUG( "Ex size 0:" << Ex2D->dims_[0]);
    DEBUG( "Ex size 1:" << Ex2D->dims_[1]);
    DEBUG( "epsilon size :" << epsilon2D->dims_[0] << " " << epsilon2D->dims_[1] );
    DEBUG( "epsilonx size :" << epsilonx2D->dims_[0] << " " << epsilonx2D->dims_[1] );
    DEBUG( "epsilony size :" << epsilony2D->dims_[0] << " " << epsilony2D->dims_[1] );
    DEBUG( "epsilonz size :" << epsilonz2D->dims_[0] << " " << epsilonz2D->dims_[1] );

    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            ( *Ex2D )( i, j ) += 1./(*epsilonx2D)(i,j)*(-dt*( *Jx2D )( i, j ) + dt_ov_dy * ( ( *Bz2D )( i, j+1 ) - ( *Bz2D )( i, j ) ));
            //( *Ex2D )( i, j ) += 1.*(-dt*( *Jx2D )( i, j ) + dt_ov_dy * ( ( *Bz2D )( i, j+1 ) - ( *Bz2D )( i, j ) ));
            // MESSAGE( "Ex j:" << j );
        }
        //MESSAGE( "Ex i:" << i );
    }
    
    // Electric field Ey^(p,d)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            ( *Ey2D )( i, j ) += 1./(*epsilony2D)(i,j)*(-dt*( *Jy2D )( i, j ) - dt_ov_dx * ( ( *Bz2D )( i+1, j ) - ( *Bz2D )( i, j ) ));
            //( *Ey2D )( i, j ) += 1.*(-dt*( *Jy2D )( i, j ) - dt_ov_dx * ( ( *Bz2D )( i+1, j ) - ( *Bz2D )( i, j ) ));
            
        }
        //MESSAGE( "Ey i:" << i );
    }
    
    // Electric field Ez^(p,p)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            
             (*Ez2D )( i, j ) += 1./(*epsilonz2D)(i,j)*(-dt*( *Jz2D )( i, j )
                                 +               dt_ov_dx * ( ( *By2D )( i+1, j ) - ( *By2D )( i, j ) )
                                 -               dt_ov_dy * ( ( *Bx2D )( i, j+1 ) - ( *Bx2D )( i, j ) ));
                 
            /*( *Ez2D )( i, j ) += 1.*(-dt*( *Jz2D )( i, j )
                                 +               dt_ov_dx * ( ( *By2D )( i+1, j ) - ( *By2D )( i, j ) )
                                 -               dt_ov_dy * ( ( *Bx2D )( i, j+1 ) - ( *Bx2D )( i, j ) ));*/

        }
        //MESSAGE( "Ez i:" << i );
    }
    
}

