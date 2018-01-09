#ifndef HOMOGENEOUS_H
#define HOMOGENEOUS_H

#include "regime.h"

//#include <boost/units/systems/si/dimensionless.hpp>
//#include <boost/units/systems/si/length.hpp>
//#include <boost/units/systems/si/mass_density.hpp>
//#include <boost/units/systems/si/dynamic_viscosity.hpp>
//#include <boost/units/systems/si/kinematic_viscosity.hpp>
//#include <boost/units/systems/si/io.hpp>

//using namespace boost::units;
//using namespace boost::units::si;

#define ACv (3.0)
#define kvK (0.5)
#define particle_ratio (0.015)
#define musf (0.415 * pascals / meter )

class Homogeneous : public Regime {

    dpdx pressureLoss(
            quantity<velocity> v,
            quantity<length> D,
            quantity<length> d,
            quantity<length> eps,
            quantity<kinematic_viscosity> nu,
            quantity<mass_density> rhow,
            quantity<mass_density> rhos,
            quantity<dimensionless> Cvs);

    dpdx relativeExcessGradient(
            quantity<velocity> v,
            quantity<length> D,
            quantity<length> d,
            quantity<length> eps,
            quantity<kinematic_viscosity> nu,
            quantity<mass_density> rhow,
            quantity<mass_density> rhos,
            quantity<dimensionless> Cvs,
            bool use_sf );
};

#endif
