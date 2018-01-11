#ifndef WATERDATA_H
#define WATERDATA_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dynamic_viscosity.hpp>
#include <boost/units/systems/si/kinematic_viscosity.hpp>
#include <boost/units/systems/si/io.hpp>

#include <boost/units/systems/temperature/celsius.hpp>

using namespace boost::units;
using namespace boost::units::si;
using namespace boost::units::celsius;

quantity<mass_density> waterDensity( quantity<celsius::temperature> T );

quantity<kinematic_viscosity> waterKinematicViscosity( quantity<celsius::temperature> T );

#endif
