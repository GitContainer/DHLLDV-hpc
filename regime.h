#ifndef REGIME_H
#define REGIME_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/pressure.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

using namespace boost::units;
using namespace boost::units::si;

typedef quantity<pressure_gradient> dpdx;

class Regime
{
public:

    // Pressure loss per meter
    virtual dpdx pressureLoss(
            quantity<velocity> v,
            quantity<length> D,
            quantity<length> d,
            quantity<length> eps,
            quantity<kinematic_viscosity> nu,
            quantity<mass_density> rhow,
            quantity<mass_density> rhos,
            quantity<dimensionless> Cvs) = 0;

    // relative excess hydraulic gradient
<<<<<<< HEAD
    virtual dpdx relativeExcessGradient(
=======
    virtual quantity<pressure_gradient> relativeExcessGradient(
>>>>>>> e17fb798a78abf6777cf3b522b808c142933f340
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
