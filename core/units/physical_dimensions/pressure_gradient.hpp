#ifndef BOOST_UNITS_PRESSURE_GRADIENT_DERIVED_DIMENSION_HPP
#define BOOST_UNITS_PRESSURE_GRADIENT_DERIVED_DIMENSION_HPP

#include <boost/units/derived_dimension.hpp>
#include <boost/units/physical_dimensions/mass.hpp>
#include <boost/units/physical_dimensions/length.hpp>
#include <boost/units/physical_dimensions/time.hpp>

namespace boost {

namespace units {

/// derived dimension for pressure_gradient : M L^-2 S^-2
typedef derived_dimension<length_base_dimension,-2,
                          mass_base_dimension, 1,
                          time_base_dimension,-2>::type pressure_gradient_dimension;

} // namespace units

} // namespace boost

#endif // BOOST_UNITS_VOLUME_FLOW_DERIVED_DIMENSION_HPP
