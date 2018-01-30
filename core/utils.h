#ifndef UTIL_H
#define UTIL_H

#define PI (4.0*atan(1.0))

#include <boost/units/systems/si/plane_angle.hpp>
#include <boost/units/cmath.hpp>

using namespace boost::units;
using namespace boost::units::si;

//#include <iostream>

quantity<plane_angle> a2beta( quantity<dimensionless> Arel );
quantity<dimensionless> beta( quantity<dimensionless> Re );



#endif
