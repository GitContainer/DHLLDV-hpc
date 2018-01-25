#ifndef SLIPRATIO_H
#define SLIPRATIO_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dynamic_viscosity.hpp>
#include <boost/units/systems/si/kinematic_viscosity.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

// Tolerance value for the iterative calculation
//#define TOL (0.01e-3 * meter_per_second)
//#define MAX_ITER 20

//#define PI (4.0*atan(1.0))

using namespace boost::units;
using namespace boost::units::si;

namespace DHLLDV {
    quantity<dimensionless> slipRatio(quantity<velocity> v, quantity<length> D, quantity<length> d, quantity<length> eps, quantity<kinematic_viscosity> nu, quantity<mass_density> rhol, quantity<mass_density> rhos, quantity<dimensionless> Cvs, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb);

    quantity<dimensionless> slipRatio_FB(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb);

    quantity<dimensionless> slipRatio_LDV(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb);

    quantity<dimensionless> slipRatio_HeHo(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<length> D, quantity<length> d, quantity<length> eps, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<acceleration> g);

    quantity<dimensionless> slipRatio_tan(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb);
}


#endif
