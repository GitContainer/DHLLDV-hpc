#include "homogeneous.h"
#include "fluidutils.h"
#include "particleutils.h"

#include <cmath>
#include <boost/units/cmath.hpp>

//#include <iostream>


dpdx Homogeneous::pressureLoss(quantity<velocity> v,
        quantity<length> D,
        quantity<length> d,
        quantity<length> eps,
        quantity<kinematic_viscosity> nu,
        quantity<mass_density> rhow,
        quantity<mass_density> rhos,
        quantity<dimensionless> Cvs)
{
    return relativeExcessGradient(v, D, d, eps, nu, rhow, rhos, Cvs, true) * relativeDensity(rhos, rhow) * Cvs
        + fluidPressureLoss(v, D, eps, nu*rhow, rhow);
}

<<<<<<< HEAD
dpdx Homogeneous::relativeExcessGradient(
=======
quantity<pressure_gradient> Homogeneous::relativeExcessGradient(
>>>>>>> e17fb798a78abf6777cf3b522b808c142933f340
        quantity<velocity> v,
        quantity<length> D,
        quantity<length> d,
        quantity<length> eps,
        quantity<kinematic_viscosity> nu,
        quantity<mass_density> rhow,
        quantity<mass_density> rhos,
        quantity<dimensionless> Cvs,
        bool use_sf = true )
{
    quantity<dimensionless> Re = reynolds(v, D, nu);
    quantity<dimensionless> lambda = frictionfactor(Re, D, eps);

    quantity<mass_density> rhom = rhow + Cvs*(rhos - rhow);

    quantity<dimensionless> ratio_dv_D = 11.6 * nu / ( sqrt(lambda/8.0) * v * d );
    if ( ratio_dv_D > 1.0 )
    {
        ratio_dv_D = 1.0;
    }

<<<<<<< HEAD
    quantity<dimensionless> sb = pow<2>( (ACv/kvK) * log( rhom/rhow ) * sqrt(lambda/8.0) + 1);
    quantity<dimensionless> top = 1.0 + relativeDensity(rhos, rhow)*Cvs - sb;
    quantity<dimensionless> bottom = relativeDensity(rhos, rhow)*Cvs*sb;
    
    dpdx il = fluidPressureLoss(v, D, eps, nu*rhow, rhow);
    quantity<dimensionless> f = d / ( particle_ratio * D );
    
    if( !use_sf || (f < 1.0) )
    {
        return il*(1.0 - (1.0-top/bottom)*(1.0-ratio_dv_D));
    }
    else
    {
        // Sliding flow as per equation 8.8-5
        return (il * (1.0 - (1.0 - top/bottom)*(1.0 - ratio_dv_D)) + (f-1)*musf)/f;
    }
=======
    quantity<dimensionless> sb = pow<2>( (ACV/KVK) * log( rhom/rhow ) * sqrt(lambda/8.0) + 1);

    quantity<pressure_gradient> il = fluidPressureLoss(v, D, eps, dynamicViscosity(nu, rhow), rhow);
    quantity<dimensionless> f = d / ( PARTICLE_RATIO * D );

    quantity<pressure_gradient> result;
    quantity<dimensionless> frac = ( 1 + relativeDensity(rhos,rhow)*Cvs - sb ) / ( relativeDensity(rhos,rhow)*Cvs*sb );

    if ( (!use_sf) || ( f < 1.0 ) )
    {
        result = il * (1 - (1 - frac)*(1 - ratio_dv_D));
    } else {
        result = ( il * ( 1 - (1 - frac)*(1 - ratio_dv_D) ) + (f-1)*MUSF ) / f;
    }

    return result;
>>>>>>> e17fb798a78abf6777cf3b522b808c142933f340
}
