#ifndef PARTICLEUTILS_H
#define PARTICLEUTILS_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dynamic_viscosity.hpp>
#include <boost/units/systems/si/kinematic_viscosity.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

using namespace boost::units;
using namespace boost::units::si;

// Relative density of solids in fluid
quantity<dimensionless> relativeDensity( quantity<mass_density> rhos, quantity<mass_density> rhol );

#endif
