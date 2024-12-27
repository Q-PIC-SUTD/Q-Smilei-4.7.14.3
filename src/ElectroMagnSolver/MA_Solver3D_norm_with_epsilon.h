#ifndef MA_SOLVER3D_NORM_with_epsilon_H
#define MA_SOLVER3D_NORM_with_epsilon_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver3D_norm_with_epsilon : public Solver3D
{

public:
    MA_Solver3D_norm_with_epsilon( Params &params );
    virtual ~MA_Solver3D_norm_with_epsilon();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:

};//END class

#endif

