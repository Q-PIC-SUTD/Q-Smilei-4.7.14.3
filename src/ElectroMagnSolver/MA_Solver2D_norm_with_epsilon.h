#ifndef MA_SOLVER2D_NORM_with_epsilon_H
#define MA_SOLVER2D_NORM_with_epsilon_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver2D_norm_with_epsilon : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Yee
    MA_Solver2D_norm_with_epsilon( Params &params );
    virtual ~MA_Solver2D_norm_with_epsilon();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:

};//END class

#endif

