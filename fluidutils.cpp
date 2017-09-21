#include "fluidutils.h"

#include <boost/units/cmath.hpp>


// Reynolds number with dynamic viscosity
quantity<dimensionless> reynolds(quantity<velocity> v, quantity<length> L, quantity<dynamic_viscosity> mu , quantity<mass_density> rho)
{
    return ( rho * v * L / mu );
}

// Reynolds number with kinematic viscosity
quantity<dimensionless> reynolds( quantity<velocity> v, quantity<length> L, quantity<kinematic_viscosity> nu)
{
    return ( v * L / nu );
}

// Approximation of the Darcy Weisbach friction factor using the formula of S.W. Churchill.
quantity<dimensionless> frictionfactor( quantity<dimensionless> Re, quantity<length> D, quantity<length> eps )
{
    /* Note that the implementation does not use boost quantities internally. This is due to the limited
     * availability of power functions for boost quantities. Input and output are fully in boost quantities.
     */

    quantity<dimensionless> A = pow<16>( 2.457 * log( pow( pow(7.0/Re, 0.9) + 0.27 * eps/D, -1) ) );
    quantity<dimensionless> B = pow<16>( 37530/Re );

    quantity<dimensionless> result = 8.0*pow( pow<12>( 8.0/Re ) + pow(A+B, -1.5), 1.0/12.0 );

    return result;
}

// Pressure loss per meter for single phase fluids
quantity<pressure_gradient> fluidPressureLoss( quantity<velocity> v, quantity<length> D, quantity<length> eps, quantity<dynamic_viscosity> mu, quantity<mass_density> rho )
{
    quantity<dimensionless> Re = reynolds(v, D, mu, rho);
    quantity<dimensionless> lambda = frictionfactor( Re, D, eps);

    return ( 0.5 * lambda / D * rho * pow<2>(v) );
}

