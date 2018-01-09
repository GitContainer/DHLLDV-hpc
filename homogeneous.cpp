#include "homogeneous.h"
#include "fluidutils.h"

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

dpdx Homogeneous::relativeExcessGradient(
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
}
