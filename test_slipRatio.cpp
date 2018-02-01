#include "limitDepositVelocity.h"
#include "slipRatio.h"
//#include "fluidutils.h"
//#include "particleUtils.h"
#include "waterData.h"

#include <iostream>
#include <fstream>

#include <boost/units/cmath.hpp>

using namespace dhlldv::LDV;

int main()
{
    // Input -------------------------------------------------------------------

    // Temperature
    quantity<celsius::temperature> T = 20.0 * celsius::degree;

    // Liquid density
    quantity<mass_density> rhow = waterDensity(T);
    rhow = 1025.0 * kilogrammes_per_cubic_metre;

    // Liquid viscosity
    quantity<kinematic_viscosity> nu = waterKinematicViscosity(T);
    nu = 0.0000013 * meter * meter_per_second;

    // Pipe diameter
    quantity<length> D = 0.762 * meter;

    // Pipe roughness
    quantity<length> eps = 4.5e-05 * meter;
    eps = 0.000001 * meter;

    // Particle diameter
    quantity<length> d = 0.001 * meter;
//    quantity<length> d = 0.000355 * meter;

    // Particle density
    quantity<mass_density> rhos = 2650.0 * kilogrammes_per_cubic_metre;

    // Spatial volumetric concentration
    quantity<dimensionless> Cvs = 0.175;
    quantity<dimensionless> Cvb = 0.6;
    Cvb = 0.55;

    // Sliding friction coefficient
    quantity<dimensionless> musf = 0.416;

    // Gravitational acceleration
    quantity<acceleration> g = 9.80665 * meters_per_second_squared;


    std::ofstream out;
    out.open("../output/slipRatio.txt");

    quantity<velocity> v;

    for (int i = 0; i <= 1000; i++ )
    {
        v = 0.01 * i * meter_per_second;

        out << v.value() << ", "
            << (dhlldv::slipRatio_HeHo(v, nu, D, d, eps, rhos, rhow, g)).value() << ", "
            << (dhlldv::slipRatio_LDV(v, nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << ", "
            << (dhlldv::slipRatio_FB(v, nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << ", "
            << (dhlldv::slipRatio_tan(v, nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << ", "
            << (dhlldv::slipRatio(v,D,d,eps,nu,rhow,rhos,Cvs,g,musf,Cvb)).value()
            << std::endl;
    }

/*

    // quantities which are shared between methods are calculated first
    quantity<dimensionless> Rsd = relativeDensity(rhos,rhow);
    // TODO: change factor to 3.2 for production
    quantity<dimensionless> ap = 3.4 * pow<static_rational<2,9> >(1.65 / Rsd);
    quantity<velocity> vt;
    quantity<dimensionless> Re_p;
    quantity<dimensionless> beta;
    quantity<dimensionless> KC;

    for (int i = -500; i <= -100; i = i + 7)
    {
        d = pow(10.0, 1.0*i/100.0) * meter;
        vt = terminalSettlingRuby(nu, d, rhos, rhow, g);
        Re_p = vt * d / nu;
        beta = ( 4.7 + 0.41 * pow<static_rational<3,4> >(Re_p) ) / ( 1.0 + 0.175 * pow<static_rational<3,4> >(Re_p) );
        KC = 0.175*(1.0 + beta);

        out << d.value() << ", "
            << (verySmallParticles(nu,rhos,rhow,D,eps,g)).value() << ", "
            << (smallParticles(nu,Cvs,Rsd,ap,vt,beta,KC,D,eps,d,g)).value() << ", "
            << (largeParticles(nu,Cvs,Rsd,ap,vt,beta,KC,D,eps,d,g,musf,Cvb)).value() << ", "
            << (limitDepositVelocity(nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << ", "
            << (lowerLimit(nu,D,d,eps,vt,Cvs,beta,KC,g,musf)).value() << std::endl;
    }*/

    out.close();

    return 0;
}
