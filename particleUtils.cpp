#include "particleUtils.h"
#include "fluidutils.h"

#include <boost/units/cmath.hpp>

#include <iostream>

quantity<velocity> terminalSettlingRuby(quantity<kinematic_viscosity> nu, quantity<length> d, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<acceleration> g)
{
    return 10.0 * nu / d * (sqrt(1.0 + 0.01 * relativeDensity(rhos,rhol) * g * pow<3>(d) / pow<2>(nu) ) - 1.0);
}

// Relative density of solids in fluid
quantity<dimensionless> relativeDensity( quantity<mass_density> rhos, quantity<mass_density> rhol )
{
    return ((rhos-rhol)/rhol);
}
