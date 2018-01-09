#include "particleutils.h"

#include <iostream>

#include <cmath>
#include <boost/units/cmath.hpp>

// Relative density of solids in fluid
quantity<dimensionless> relativeDensity( quantity<mass_density> rhos, quantity<mass_density> rhol )
{
    return ((rhos-rhol)/rhol);
}

quantity<velocity> terminalSettlingVelocity( quantity<length> diameter, quantity<mass_density>
                                            rhos, quantity<mass_density> rhol)
{
    quantity<velocity> vt;
    quantity<dimensionless> Rsd = relativeDensity(rhos, rhol);
    
    // Convert boost units to doubles for internal use
    double d = 1000.0 * diameter/meter;
    
    if (diameter < 0.1e-3*meter)
    {
        // Stokes region
        vt = 0.424 * Rsd * pow(d, 2) * meter_per_second;
    }
    else if (diameter > 1e-3*meter)
    {
        // Rittinger equation
        vt = 0.087 * sqrt( Rsd * d ) * meter_per_second;
    }
    else
    {
        // Transition zone, Budryck equation
        vt = 0.008925/d * ( sqrt( 1.0 + 95.0 * Rsd * pow(d, 3) ) - 1.0 ) * meter_per_second;
    }
    
    return vt;
}
