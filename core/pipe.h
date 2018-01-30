#ifndef PIPE_H
#define PIPE_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/plane_angle.hpp>
#include <boost/units/systems/si/io.hpp>

using namespace boost::units;
using namespace boost::units::si;

class Pipe {
public:
    quantity<length> D() { return diameter; };
    void setD( quantity<length> D ) { diameter = D; };

    quantity<length> eps() { return roughness; };
    void setEps( quantity<length> e ) { roughness = e; };

    quantity<plane_angle> a() { return inclination; };
    void setA( quantity<plane_angle> a ) { inclination = a; };

private:
    quantity<length> diameter = 0.762 * meter;
    quantity<length> roughness = 0.000001 * meter;
    quantity<plane_angle> inclination = 0.0 * radians;
};

#endif

