#include "particleutils.h"

// Relative density of solids in fluid
quantity<dimensionless> relativeDensity( quantity<mass_density> rhos, quantity<mass_density> rhol )
{
    return ((rhos-rhol)/rhol);
}
