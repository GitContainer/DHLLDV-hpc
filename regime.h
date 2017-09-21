#ifndef REGIME_H
#define REGIME_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/pressure.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

using namespace boost::units;
using namespace boost::units::si;

class Regime
{
public:

    // Pressure loss per meter
    virtual quantity<pressure_gradient> pressureLoss(
            quantity<velocity> v,
            quantity<length> D,
            quantity<length> d,
            quantity<length> eps,
            quantity<kinematic_viscosity> nu,
            quantity<mass_density> rhow,
            quantity<mass_density> rhos,
            quantity<dimensionless> Cvs) = 0;

    // relative excess hydraulic gradient
    virtual quantity<dimensionless> relativeExcessGradient(
            quantity<velocity> v,
            quantity<length> D,
            quantity<length> d,
            quantity<length> eps,
            quantity<kinematic_viscosity> nu,
            quantity<mass_density> rhow,
            quantity<mass_density> rhos,
            quantity<dimensionless> Cvs,
            bool use_sf = true
            ) = 0;
};

#endif
