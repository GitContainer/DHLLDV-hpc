#include "utils.h"

quantity<plane_angle> a2beta( quantity<dimensionless> Arel )
{
    int it = 0, _MAX_ITER = 20;
    quantity<plane_angle> beta = 0.25 * PI * radians, _TOL = 1e-6 * radians, db = 10.0*_TOL;

    while( (abs(db) > _TOL) && (it++ < _MAX_ITER) )
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
