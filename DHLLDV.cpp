#include "DHLLDV.h"

#include "core/utils.h"

#include <boost/units/cmath.hpp>

#include <iostream>

DHLLDV::DHLLDV()
{

}

DHLLDV::~DHLLDV()
{

}

dpdx DHLLDV::pressureLoss( quantity<velocity> v )
{
    dpdx result = 1.0 * pascals_per_meter;

    return result;
}

quantity<dimensionless> DHLLDV::homogeneousHeadLoss( quantity<velocity> v )
{
    return frictionfactor(reynolds(v, D, nu), D, eps) * pow<2>(v) / (2.0 * g * D);
}

quantity<dimensionless> DHLLDV::fixedBedErhg( quantity<velocity> v )
{
    return (fixedBedHeadLoss(v) - homogeneousHeadLoss(v)) / ( relativeDensity(rhos,rhol) * Cvs );
}

quantity<dimensionless> DHLLDV::fixedBedHeadLoss( quantity<velocity> v )
{
    return (fixedBedPressureLoss(v) / (rhol * g) );
}

dpdx DHLLDV::fixedBedPressureLoss(quantity<velocity> v)
{
    quantity<area> Ap, A2, A1;
    Ap = 0.25 * PI * pow<2>(D);
    A2 = Ap * Cvs / Cvb;
    A1 = Ap - A2;

    quantity<plane_angle> beta = a2beta( Cvs/Cvb );

    quantity<length> O1, O12;
    O1 = (PI*radians - beta)/radians * D;
    O12 = D * sin(beta);

    quantity<length> DH1 = 4.0 * A1 / (O1 + O12);
    quantity<velocity> v1 = v * Ap / A1, v2 = 0.0 * meter_per_second;

    quantity<dimensionless> lambda1 = frictionfactor( reynolds(v,DH1,nu), DH1, eps );
    quantity<pressure_gradient> F1_l = (lambda1 * rhol * pow<2>(v1) * 0.125) * O1 / A1, F12_l;

    quantity<dimensionless> lambda12 = frictionfactor( reynolds(v1-v2,DH1,nu), DH1, eps ), lambda12_sf;

    lambda12_sf = 0.83 * lambda1 + 0.37 * pow<static_rational<273, 100> >( (v1-v2)/sqrt(2.0*g*DH1*relativeDensity(rhos,rhol)) ) * pow<static_rational<94, 1000> >( rhos*PI/6.0 * pow<3>(d/meter)/rhol );
    lambda12 = (lambda12 > lambda12_sf) ? lambda12 : lambda12_sf;

    F12_l = (lambda12 * rhol * pow<2>(v) * 0.125) * O12 / A1;

    return (F1_l + F12_l);
}

quantity<dimensionless> DHLLDV::slidingBedErhg( quantity<velocity> v )
{
    return musf;
}

quantity<dimensionless> DHLLDV::slidingBedHeadLoss( quantity<velocity> v )
{
    return slidingBedErhg(v) * relativeDensity(rhos,rhol) * Cvs + homogeneousHeadLoss(v);
}

quantity<dimensionless> DHLLDV::heterogeneousErhg( quantity<velocity> v )
{
    quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g);
    quantity<dimensionless> b = beta( reynolds(vt,d,nu) );
    quantity<dimensionless> Kc = 0.175 * (1.0 + b);
    quantity<dimensionless> lambda = frictionfactor( reynolds(v,D,nu), D, eps );

    quantity<dimensionless> result = vt * pow(1.0 - Cvs/Kc, b) / v + 72.25 * (1.0/lambda) * pow<static_rational<10,3> >(vt / sqrt(g*d)) * pow<2>( pow<static_rational<1,3> >(nu * g) / v ) ;

    return result;
}

quantity<dimensionless> DHLLDV::heterogeneousHeadLoss( quantity<velocity> v )
{
    return (heterogeneousErhg(v) * relativeDensity(rhos,rhol) * Cvs + homogeneousHeadLoss(v) );
}

dpdx DHLLDV::heterogeneousPressureLoss( quantity<velocity> v )
{
    return heterogeneousHeadLoss(v) * g * rhol;
}
