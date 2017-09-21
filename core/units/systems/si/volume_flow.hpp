#ifndef BOOST_UNITS_SI_VOLUME_FLOW_HPP
#define BOOST_UNITS_SI_VOLUME_FLOW_HPP

#include <boost/units/systems/si/base.hpp>
#include "../../physical_dimensions/volume_flow.hpp"

namespace boost {

namespace units {

namespace si {

typedef unit<volume_flow_dimension,si::system>      volume_flow;

BOOST_UNITS_STATIC_CONSTANT(cubic_meter_per_second,volume_flow);
BOOST_UNITS_STATIC_CONSTANT(cubic_meters_per_second,volume_flow);
BOOST_UNITS_STATIC_CONSTANT(cubic_metre_per_second,volume_flow);
BOOST_UNITS_STATIC_CONSTANT(cubic_metres_per_second,volume_flow);

} // namespace si

} // namespace units

} // namespace boost

#endif // BOOST_UNITS_SI_VOLUME_FLOW_HPP
