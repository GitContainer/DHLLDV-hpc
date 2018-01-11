#include "limitDepositVelocity.h"
#include "fluidutils.h"
#include "particleUtils.h"

#include <boost/units/cmath.hpp>
#include <boost/units/systems/si/dose_equivalent.hpp>

#include <iostream>

quantity<velocity> DHLLDV::LDV::limitDepositVelocity(quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb)
{
    quantity<velocity> v_r, v_s;

    // quantities which are shared between methods are calculated first
    quantity<dimensionless> Rsd = relativeDensity(rhos,rhol);
    // TODO: change factor to 3.2 for production
    quantity<dimensionless> ap = 3.4 * pow<static_rational<2,9> >(1.65 / Rsd);
    quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g);
    quantity<dimensionless> Re_p = vt * d / nu;
    quantity<dimensionless> beta = ( 4.7 + 0.41 * pow<static_rational<3,4> >(Re_p) ) / ( 1.0 + 0.175 * pow<static_rational<3,4> >(Re_p) );
    quantity<dimensionless> KC = 0.175*(1.0 + beta);

    // Comparison between "very small particles" and "small particles", result saved in v_s
    v_r = verySmallParticles(nu, rhos, rhol, D, eps, g);
    v_s = smallParticles(nu, Cvs, Rsd, ap, vt, beta, KC, D, eps, d, g);

    if ( v_r > v_s )
    {
        v_s = v_r;
    }

    // The "large particles" and "very large particles", saved in v_r
    v_r = largeParticles(nu, Cvs, Rsd, ap, vt, beta, KC, D, eps, d, g, musf, Cvb);

    // The upper limit
    v_r = upperLimit(d, Rsd, v_s, v_r);

    // The lower limit
    v_s = lowerLimit(nu, D, d, eps, vt, Cvs, beta, KC, g, musf, v_s, v_r);

    return (v_s > v_r) ? v_s : v_r;
}

quantity<velocity> DHLLDV::LDV::verySmallParticles(quantity<kinematic_viscosity> nu, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<acceleration> gravity)
{
    quantity<velocity> result = 10.0 * TOL;
    quantity<velocity> v0 = -result;

    int it = 0;

    while ( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;
        result = 1.4 * pow<static_rational<1,3> >(nu * relativeDensity(rhos, rhol) * gravity)
            * sqrt( 8.0 / frictionfactor(reynolds(v0, D, nu), D, eps) );
    }

    return result;
}

quantity<velocity> DHLLDV::LDV::smallParticles(quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<dimensionless> Rsd, quantity<dimensionless> ap, quantity<velocity> vt, quantity<dimensionless> beta,
                                               quantity<dimensionless> KC, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> gravity)
{
    quantity<velocity> result = 10.0 * TOL;
    quantity<velocity> v0 = -result;

    int it = 0;

    while ( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;
        result = ap * pow<static_rational<1,3> >(vt * pow(1-Cvs/KC, beta) * Cvs * (2.0 * gravity * Rsd * D) / frictionfactor(reynolds(v0, D, nu), D, eps));
    }

    return result;
}

quantity<velocity> DHLLDV::LDV::largeParticles(quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvs, quantity<dimensionless> Rsd, quantity<dimensionless> ap, quantity<velocity> vt, quantity<dimensionless> beta,
                                               quantity<dimensionless> KC, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> gravity, quantity<dimensionless> musf, quantity<dimensionless> Cvb)
{
    // Distuingish between large particles and very large particles
    quantity<dimensionless> Cvr_ldv;
    if ( d <= 0.015 * D )
    {
        Cvr_ldv = 0.0065 / ( 2.0 * gravity * Rsd * D ) * pow<2>(meters_per_second);
    }
    else
    {
        Cvr_ldv = 0.053 / (2.0 * gravity * Rsd * D) * sqrt(d/D) * pow<2>(meters_per_second);
    }

    quantity<velocity> result = 10.0 * TOL;
    quantity<velocity> v0 = -result;

    int it = 0;

    while ( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;
        result = ap * pow<static_rational<1,3> >( pow(1.0-Cvs/KC, beta) * Cvs * sqrt(musf * Cvb * PI/8.0) * sqrt(Cvr_ldv) / frictionfactor(reynolds(v0, D, nu), D, eps) ) * sqrt( 2.0 * gravity * Rsd * D);
    }

    return result;
}

quantity<velocity> DHLLDV::LDV::lowerLimit(quantity<kinematic_viscosity> nu, quantity<length> D, quantity<length> d, quantity<length> eps, quantity<velocity> vt, quantity<dimensionless> Cvs, quantity<dimensionless> beta,
                                           quantity<dimensionless> KC, quantity<acceleration> g, quantity<dimensionless> musf, quantity<velocity> v_s, quantity<velocity> v_r)
{
    quantity<velocity> result = 10.0 * TOL;
    quantity<velocity> v0 = - result;

    quantity<velocity> B;
    quantity<dose_equivalent> C;

    int it = 0;

    while ( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;

        B = vt * pow(1.0-Cvs/KC, beta);
        C = 72.25 / frictionfactor( reynolds(v0,D,nu), D, eps ) * pow<static_rational<10,3> >(vt/sqrt(g*d)) * pow<static_rational<2,3> >(nu*g) / musf;

        result = 0.5*B + 0.5*sqrt( pow<2>(B) + 4.0 * C );
    }

    return result;
}

quantity<velocity> DHLLDV::LDV::upperLimit(quantity<length> d, quantity<dimensionless> Rsd, quantity<velocity> v_s, quantity<velocity> v_r)
{
    // Settings for the method
    quantity<length> drough = 0.002 * meter; // valid for sand with Rsd=1.65

    if (!(d > drough))
    {
        if( v_s < v_r )
        {
            v_r = v_s;
        }
        else
        {
            quantity<length> d0 = 0.0005 * 1.65 / Rsd * meter;
            v_r = v_s * exp(-d/d0) + v_r * (1.0 - exp(-d/d0) );
        }
    }

    return v_r;
}
