#ifndef BOOST_UNITS_VOLUME_FLOW_DERIVED_DIMENSION_HPP
#define BOOST_UNITS_VOLUME_FLOW_DERIVED_DIMENSION_HPP

#include <boost/units/derived_dimension.hpp>
#include <boost/units/physical_dimensions/length.hpp>
#include <boost/units/physical_dimensions/time.hpp>

namespace boost {

namespace units {

/// derived dimension for volume_flow : L^3 T^-1
typedef derived_dimension<length_base_dimension,3,
                          time_base_dimension,-1>::type volume_flow_dimension;

} // namespace units

} // namespace boost

#endif // BOOST_UNITS_VOLUME_FLOW_DERIVED_DIMENSION_HPP
