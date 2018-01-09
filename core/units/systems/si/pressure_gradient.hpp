#ifndef BOOST_UNITS_SI_PRESSURE_GRADIENT_HPP
#define BOOST_UNITS_SI_PRESSURE_GRADIENT_HPP

#include <boost/units/systems/si/base.hpp>
#include "../../physical_dimensions/pressure_gradient.hpp"

namespace boost {

namespace units {

namespace si {

typedef unit<pressure_gradient_dimension,si::system>      pressure_gradient;

BOOST_UNITS_STATIC_CONSTANT(pascal_per_meter,pressure_gradient);
BOOST_UNITS_STATIC_CONSTANT(pascal_per_meters,pressure_gradient);
BOOST_UNITS_STATIC_CONSTANT(pascals_per_meter,pressure_gradient);
BOOST_UNITS_STATIC_CONSTANT(pascals_per_meters,pressure_gradient);

} // namespace si

} // namespace units

} // namespace boost

#endif // BOOST_UNITS_SI_VOLUME_FLOW_HPP
