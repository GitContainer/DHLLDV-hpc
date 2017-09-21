#ifndef FLUIDUTILS_H
#define FLUIDUTILS_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dynamic_viscosity.hpp>
#include <boost/units/systems/si/kinematic_viscosity.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

using namespace boost::units;
using namespace boost::units::si;

// Reynolds number with dynamic viscosity
quantity<dimensionless> reynolds( quantity<velocity> v, quantity<length> L, quantity<dynamic_viscosity> mu, quantity<mass_density> rho );

// Reynolds number with kinematic viscosity
quantity<dimensionless> reynolds( quantity<velocity> v, quantity<length> L, quantity<kinematic_viscosity> nu);

// Approximation of the Darcy Weisbach friction factor using the formula of S.W. Churchill.
quantity<dimensionless> frictionfactor( quantity<dimensionless> Re, quantity<length> D, quantity<length> eps );

// Pressure loss per meter for single phase fluids
quantity<pressure_gradient> fluidPressureLoss( quantity<velocity> v, quantity<length> D, quantity<length> eps, quantity<dynamic_viscosity> mu, quantity<mass_density> rho );

#endif
