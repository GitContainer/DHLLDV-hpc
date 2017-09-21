#include "fluidutils.h"
#include "homogeneous.h"

//#include "core/units/systems/si/pressure_gradient.hpp"

#include <iostream>

using namespace boost::units;
using namespace boost::units::si;

int main()
{
    quantity<dynamic_viscosity> mu = 1.002e-3 *pascals*seconds;
    quantity<mass_density> rho = 1000 * kilogrammes_per_cubic_metre;
    quantity<velocity> v = 5.0 * meters_per_second;
    quantity<length> D = 0.5 * meter;
    quantity<length> L = 1200.0 * meter;

    quantity<dimensionless> Re = reynolds(v, D, mu, rho);
    //quantity<kinematic_viscosity> nu = mu/rho;

    std::cout << "Reynolds number: " << Re << std::endl;
    //std::cout << reynolds(v, D, mu/rho) << std::endl;
    //std::cout << reynolds(v, D, nu) << std::endl;

    std::cout << "Friction factor: " << frictionfactor( Re, D, 0.0*meter ) << std::endl;

    Regime *ho = new Homogeneous();
    std::cout << ho->pressureLoss(1.0*meter_per_second, 1.0*meter, 1.0*meter, 0.0*meter, 1.0*meter*meter/second, 1000.0*kilogrammes_per_cubic_metre, 2650.0*kilogrammes_per_cubic_metre, 0.2) << std::endl;

    std::cout << fluidPressureLoss(v, D, 0.0*meter, mu, rho) * L << std::endl;

    return 0;
}
