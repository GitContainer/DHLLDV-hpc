#ifndef LIMITDEPOSITVELOCITY_H
#define LIMITDEPOSITVELOCITY_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dynamic_viscosity.hpp>
#include <boost/units/systems/si/kinematic_viscosity.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

// Tolerance value for the iterative calculation
#define TOL (0.01e-3 * meter_per_second)
#define MAX_ITER 20

#define PI (4.0*atan(1.0))

using namespace boost::units;
using namespace boost::units::si;

namespace DHLLDV {
    namespace LDV {
        quantity<velocity> limitDepositVelocity(quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> gravity, quantity<dimensionless> musf, quantity<dimensionless> Cvb);

        quantity<velocity> verySmallParticles(quantity<kinematic_viscosity> nu, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<acceleration> gravity);

        quantity<velocity> smallParticles(quantity<kinematic_viscosity> nu,
            quantity<dimensionless> Cvs,
            quantity<dimensionless> Rsd, quantity<dimensionless> ap, quantity<velocity> vt, quantity<dimensionless> beta, quantity<dimensionless> KC, quantity<length> D,
            quantity<length> eps, quantity<length> d,
            quantity<acceleration> gravity);

        quantity<velocity> largeParticles(quantity<kinematic_viscosity> nu,
            quantity<dimensionless> Cvs, quantity<dimensionless> Rsd, quantity<dimensionless> ap, quantity<velocity> vt, quantity<dimensionless> beta, quantity<dimensionless> KC, quantity<length> D,
            quantity<length> eps, quantity<length> d,
            quantity<acceleration> gravity, quantity<dimensionless> musf, quantity<dimensionless> Cvb);

        quantity<velocity> lowerLimit(quantity<kinematic_viscosity> nu, quantity<length> D, quantity<length> d, quantity<length> eps, quantity<velocity> vt, quantity<dimensionless> Cvs, quantity<dimensionless> beta, quantity<dimensionless> KC, quantity<acceleration> g, quantity<dimensionless> musf, quantity<velocity> v_s, quantity<velocity> v_r);

        quantity<velocity> upperLimit(quantity<length> d, quantity<dimensionless> Rsd, quantity<velocity> v_s, quantity<velocity> v_r);
    }
}

#endif
