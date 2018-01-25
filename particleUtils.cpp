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

quantity<dimensionless> sqrtCx(quantity<length> d, quantity<kinematic_viscosity> nu, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<acceleration> g)
{
    quantity<dimensionless> factor = 0.6, smallFactor = 1.8;

    quantity<dimensionless> Fr = terminalSettlingRuby(nu,d,rhos,rhol,g) / sqrt(g*d);
    quantity<dimensionless> wilson = 0.226 * pow<static_rational<1,6> >(g / d * second * second);
    quantity<dimensionless> gibert = pow<static_rational<-10, 9> >(Fr);

    if (gibert > smallFactor)
    {
        gibert = pow<static_rational<1,4> >(smallFactor) * pow<static_rational<3,4> >(gibert);
    }

    if (gibert < wilson)
    {
        gibert = gibert * factor + wilson * (1.0 - factor);
    }

    return gibert;
}
