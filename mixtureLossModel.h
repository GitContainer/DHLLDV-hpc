#ifndef MIXTURELOSSMODEL_H
#define MIXTURELOSSMODEL_H

#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dynamic_viscosity.hpp>
#include <boost/units/systems/si/kinematic_viscosity.hpp>
#include <boost/units/systems/si/io.hpp>

#include "core/units/systems/si/pressure_gradient.hpp"

using namespace boost::units;
using namespace boost::units::si;

typedef quantity<pressure_gradient> dpdx;

class mixtureLossModel
{
public:
    virtual dpdx pressureLoss( quantity<velocity> v ) = 0;

    void setRhoLiquid( quantity<mass_density> rho ) { rhol = rho; };
    quantity<mass_density> getRhoLiquid() { return rhol; };

    void setKinematicViscosity( quantity<kinematic_viscosity> n ) { nu = n; };
    quantity<kinematic_viscosity> getKinematicViscosity() { return nu; };

    void setPipeDiameter( quantity<length> Dp ) { D = Dp; };
    quantity<length> getPipeDiameter() { return D; };

    void setPipeRoughness( quantity<length> e ) { eps = e; };
    quantity<length> getPipeRoughness() { return eps; };

    void setRhoSolids( quantity<mass_density> rho ) { rhos = rho; };
    quantity<mass_density> getRhoSolids() { return rhos; };

    void setSpatialVolumetricConcentration( quantity<dimensionless> cv ) { Cvs = cv; };
    quantity<dimensionless> getSpatialVolumentricConcentration() { return Cvs; };

    quantity<mass_density> getRhoMixture() { return ( Cvs * (rhos - rhol) + rhol ); };

protected:
    quantity<mass_density> rhol = 1025.0 * kilogrammes_per_cubic_metre;
    quantity<kinematic_viscosity> nu = 1.3e-6 * meter * meter_per_second;
    quantity<length> D = 0.762 * meter;
    quantity<length> eps = 0.000001 * meter;

    // Something with a soil class for particle diameter ...
    quantity<length> d = 0.001 * meter;

    quantity<mass_density> rhos = 2585.0 * kilogrammes_per_cubic_metre;
    quantity<dimensionless> Cvs = 0.175;

};

#endif
