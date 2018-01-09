#include "fluidutils.h"

#include <boost/units/cmath.hpp>
#include <iostream>

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
    quantity<dimensionless> A = pow<16>( 2.457 * log( pow<-1>( pow<static_rational<9, 10> >(7.0/Re) + 0.27 * eps/D) ) );
    quantity<dimensionless> B = pow<16>( 37530/Re );

    quantity<dimensionless> result = 8.0*pow<static_rational<1, 12> >( pow<12>( 8.0/Re ) + pow<static_rational<-3, 2> >(A+B));

    return result;
}

// Pressure loss per meter for single phase fluids
quantity<pressure_gradient> fluidPressureLoss( quantity<velocity> v, quantity<length> D, quantity<length> eps, quantity<dynamic_viscosity> mu, quantity<mass_density> rho )
{
    quantity<dimensionless> Re = reynolds(v, D, mu, rho);
    quantity<dimensionless> lambda = frictionfactor( Re, D, eps);

    return ( 0.5 * lambda / D * rho * pow<2>(v) );
}
<<<<<<< HEAD
=======

// Kinematic viscosity from dynamic viscosity and density
quantity<kinematic_viscosity> kinematicViscosity( quantity<dynamic_viscosity> mu, quantity<mass_density> rho )
{
    return (mu / rho);
}

// Dynamic viscosity from kinemati viscosity and density
quantity<dynamic_viscosity> dynamicViscosity( quantity<kinematic_viscosity> nu, quantity<mass_density> rho )
{
    return ( nu * rho );
}
>>>>>>> e17fb798a78abf6777cf3b522b808c142933f340
