#ifndef UTIL_H
#define UTIL_H

#define PI (4.0*atan(1.0))

#include <boost/units/systems/si/plane_angle.hpp>
#include <boost/units/cmath.hpp>

using namespace boost::units;
using namespace boost::units::si;

#include <iostream>

quantity<plane_angle> a2beta( quantity<dimensionless> Arel )
{
    int it = 0, MAX_ITER = 20;
    quantity<plane_angle> beta = 0.25 * PI * radians, TOL = 1e-6 * radians, db = 10.0*TOL;

    while( (abs(db) > TOL) && (it++ < MAX_ITER) )
    {
        db = (Arel * PI * radians - beta + (0.5 * sin(2.0*beta) ) * radians) / (-1.0 + cos(2.0*beta) );
        beta = beta - db;
    }

    return beta;
}

quantity<dimensionless> beta( quantity<dimensionless> Re )
{
    return ((4.7 + 0.41 * pow<static_rational<3,4> >(Re)) / (1.0 + 0.175 * pow<static_rational<3,4> >(Re)));
}



#endif
