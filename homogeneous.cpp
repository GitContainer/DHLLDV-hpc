#include "homogeneous.h"
#include "fluidutils.h"

#include <cmath>
#include <boost/units/cmath.hpp>

quantity<pressure_gradient> Homogeneous::pressureLoss(quantity<velocity> v,
                                                      quantity<length> D,
                                                      quantity<length> d,
                                                      quantity<length> eps,
                                                      quantity<kinematic_viscosity> nu,
                                                      quantity<mass_density> rhow,
                                                      quantity<mass_density> rhos,
                                                      quantity<dimensionless> Cvs)
{
    return 1.0*pascals_per_meter;
}

quantity<dimensionless> Homogeneous::relativeExcessGradient(
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



    return 1.0;
}
