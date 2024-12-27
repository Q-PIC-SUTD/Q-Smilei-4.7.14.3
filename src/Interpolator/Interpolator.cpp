#include "Interpolator.h"

#include <cmath>
#include <iostream>

#include "Params.h"
#include "Patch.h"
#include "Random.h"
using namespace std;

Interpolator::Interpolator( Params &params, Patch *patch )
{
    rand_ = new Random(params.random_seed);
}

