#include "DHLLDV.h"
#include "core/utils.h"

#include <boost/units/cmath.hpp>

#include <iostream>
#include <algorithm>

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

quantity<dimensionless> DHLLDV::carrierHeadLoss( quantity<velocity> v )
{
    quantity<velocity> vt = terminalSettlingRuby(nu,d,rhos,rhol,g,shapeFactor);
    quantity<velocity> hvt = vt * pow( 1.0 - Cvs, beta(Re(vt,d,nu)) );

    return ff(Re(v, p.D(), nu), p.D(), p.eps()) * pow<2>(v + hvt * sin(p.a()) * Cvs ) / (2.0 * g * p.D()) + sin(p.a());
}

quantity<dimensionless> DHLLDV::homogeneousHeadLoss( quantity<velocity> v )
{
    quantity<dimensionless> result;

    result = ((rhos/rhol - 1.0) * carrierHeadLoss(v) / (Cvs * relativeDensity(rhos,rhol) ) - sin(p.a()) ) * ErhgELMFactor(v);
    result *= (HeHoTransition) ? HoMobilisationFactor(v) : (quantity<dimensionless>)1.0;
    result += sin(p.a());

    result = result * Cvs * relativeDensity(rhos,rhol) + carrierHeadLoss(v);

    return result;
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
    Ap = 0.25 * PI * pow<2>(p.D());
    A2 = Ap * Cvs / Cvb;
    A1 = Ap - A2;

    quantity<plane_angle> beta = a2beta( Cvs/Cvb );

    quantity<length> O1, O12;
    O1 = (PI*radians - beta)/radians * p.D();
    O12 = p.D() * sin(beta);

    quantity<length> DH1 = 4.0 * A1 / (O1 + O12);
    quantity<velocity> v1 = v * Ap / A1, v2 = 0.0 * meter_per_second;

    quantity<dimensionless> lambda1 = ff( Re(v,DH1,nu), DH1, p.eps() );
    quantity<pressure_gradient> F1_l = (lambda1 * rhol * pow<2>(v1) * 0.125) * O1 / A1, F12_l;

    quantity<dimensionless> lambda12 = ff( Re(v1-v2,DH1,nu), DH1, p.eps() ), lambda12_sf;

    lambda12_sf = 0.83 * lambda1 + 0.37 * pow<static_rational<273, 100> >( (v1-v2)/sqrt(2.0*g*DH1*relativeDensity(rhos,rhol)) ) * pow<static_rational<94, 1000> >( rhos*PI/6.0 * pow<3>(d/meter)/rhol );
    lambda12 = (lambda12 > lambda12_sf) ? lambda12 : lambda12_sf;

    F12_l = (lambda12 * rhol * pow<2>(v1) * 0.125) * O12 / A1;

    return (F1_l + F12_l + sin(p.a())*rhol*g);
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
    quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g, shapeFactor);
    quantity<dimensionless> b = beta( Re(vt,d,nu) );
    quantity<dimensionless> Kc = 0.175 * (1.0 + b);
    quantity<dimensionless> lambda = ff( Re(v,p.D(),nu), p.D(), p.eps() );

    quantity<velocity> ustar = sqrt(0.125*lambda) * v;

    if ( (5.8 * ustar) < (vt * sin(p.a())) )
    {
        ustar = vt * sin(p.a()) / 5.8;
    }

    quantity<dimensionless> cFac = pow<static_rational<4,3> >((vt*cos(p.a()) / (11.6*ustar - vt*sin(p.a()))) / (vt / (11.6*ustar)) );
    quantity<dimensionless> result = vt * cos(p.a()) * pow(1.0 - Cvs/Kc, b) / v;

    quantity<dimensionless> mobiFac = HeHoTransition ? HeMobilisationFactor(v) : (quantity<dimensionless>)1.0;

    result += cFac * pow<2>(1.25*cHe) / lambda / pow<3>(sqrtCx(d, nu, rhos, rhol, g)) * pow<2>(pow<static_rational<1,3> >(nu * g) / v) * mobiFac;

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

quantity<dimensionless> DHLLDV::ErhgELMFactor(quantity<velocity> v)
{
    quantity<dimensionless> lambda = ff(Re(v,p.D(),nu), p.D(), p.eps());
    quantity<velocity> ustar = sqrt( 0.125 * lambda ) * v;
    quantity<length> deltaV = 11.6 * nu / ustar;
    quantity<dimensionless> ratio = deltaV / d;

    ratio = std::min(ratio, (quantity<dimensionless>)1.0);

    quantity<dimensionless> Rsd = relativeDensity(rhos,rhol);

    quantity<dimensionless> factor = pow<2>(Acv / 0.4 * log(1.0 + Rsd*Cvs) * sqrt(lambda * 0.125) + 1.0 );
    quantity<dimensionless> reduction = (1.0 + Rsd*Cvs - factor) / (Rsd*Cvs*factor);

    return (1.0 - (1.0 - reduction) * (1 - ratio));
}

quantity<dimensionless> DHLLDV::particleLiftRatio(quantity<velocity> v)
{
    quantity<dimensionless> Cvr = Cvs / Cvb;
    quantity<dimensionless> lambda = ff( Re(v,p.D(),nu), p.D(), p.eps() );
    quantity<velocity> ustar = sqrt(0.125*lambda) * v, vth, vt = terminalSettlingRuby(nu,d,rhos,rhol,g,shapeFactor);
    quantity<dimensionless> b = beta(Re(vt,d,nu));

    vth = vt * pow(1.0 - Cvs / (0.175 * (1.0 + b)), b);

    quantity<force> FL, FK, FG;

    FL = CL * rhol * pow<2>(ustar) * PI * pow<2>(d) * (1.0 - Cvr) * 0.125;
    FK = shapeFactor * 0.5 * rhos * PI / 6.0 * pow<3>(d) * pow<2>(vth) * ustar / (nTimesThickness * 11.6 * nu);
    FG = shapeFactor * (rhos - rhol) * PI * pow<3>(d) * g;

    return ( FL / (FK + FG) );
}

quantity<dimensionless> DHLLDV::HeMobilisationFactor(quantity<velocity> v)
{
    quantity<dimensionless> result, smoothing = 0.5, LR;
    LR = particleLiftRatio(v);

    if (pow<2>(LR) > smoothing)
    {
        result = (1.0 - smoothing) * (smoothing / pow<2>(LR) );
    }
    else
    {
        result = (1.0 - pow<2>(LR));
    }

    return std::max(result, (quantity<dimensionless>)1e-6);
}

quantity<dimensionless> DHLLDV::HoMobilisationFactor(quantity<velocity> v) // TODO
{
    quantity<dimensionless> Cvr, alphaSM, CvsLDV, b, RZLDV, lambdaLDV, lambdaVls, CvBottom, CvsVls, RZVls;
    Cvr = Cvs / Cvb;
    alphaSM = 0.9847 + 0.304 * Cvr - 1.196 * pow<2>(Cvr) - 0.5564 * pow<3>(Cvr) + 0.47 * pow<4>(Cvr);

    CvsLDV = Cvb * exp( -0.5 * alphaSM / Cvr );
    quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g, shapeFactor), vstar, vstar_v;
    b = beta( Re(vt, d, nu) );

    RZLDV = pow(1.0 - CvsLDV, b);
    lambdaLDV = ff( Re(LDV(), p.D(), nu), p.D(), p.eps() );
    vstar = LDV() * sqrt( 0.125 * lambdaLDV );

    lambdaVls = ff( Re(v, p.D(), nu), p.D(), p.eps() );
    vstar_v = v * sqrt(0.125 * lambdaVls);

    CvBottom = Cvb * alphaSM * (vstar / vstar_v) / ( 1.0 - exp(-alphaSM / Cvr * vstar / vstar_v) );
    CvsVls = CvBottom * exp(-alphaSM / Cvr * vstar / vstar_v * 0.5);
    RZVls = pow( 1.0 - CvsVls, b );

    return exp( -dFraction * alphaSM / Cvr * vstar / vstar_v * RZVls / RZLDV );
}

quantity<velocity> DHLLDV::LDVsmooth()
{
    quantity<dimensionless> A1, CvMax, kappa, lambda, Rsd, b;

    Rsd = relativeDensity(rhos,rhol);

    A1 = alphap * pow( 1.65 / Rsd, 2.0/9.0 );
    CvMax = 1.0 / (1.0 + 4.7);

    quantity<velocity> TOL = 1e-6 * meter_per_second, result = 10.0 * TOL, v0 = -result, vt, lowlim;
    vt = terminalSettlingRuby(nu,d,rhos,rhol,g,shapeFactor);

    b = beta(Re(vt,d,nu));
    kappa = CvMax * (1.0 + b);

    int it = 0, MAX_ITER = 20;

    while( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;
        lambda = ff(Re(result, p.D(), nu), p.D(), p.eps());
        result = A1 * pow<static_rational<1,3> >(vt * pow(1.0 - Cvs/kappa, b) * Cvs / (lambda * sqrt(2.0 * p.D() * g * Rsd)) ) * sqrt(2.0 * g * Rsd * p.D());
    }

    lowlim = 0.85 * 2.2637 * pow<static_rational<1,3> >( musf * g * Rsd * nu * Cvb ) * sqrt(lambda * 0.125);

    return std::max(lowlim, result);
}

quantity<velocity> DHLLDV::LDVrough()
{
    quantity<dimensionless> A1, CvMax, kappa, lambda, Rsd, b, Cvrldv;

    Rsd = relativeDensity(rhos,rhol);

    A1 = alphap * pow(1.65 / Rsd, 2.0/9.0);
    CvMax = 1.0 / (1.0 + 4.7);

    quantity<velocity> TOL = 1e-6 * meter_per_second, result = 1e6 * TOL, vt, v0 = -result;
    vt = terminalSettlingRuby(nu,d,rhos,rhol,g,shapeFactor);
    b = beta(Re(vt,d,nu));
    kappa = CvMax * (1.0 + b);

    Cvrldv = 0.00017 / ( p.D()/meter * Rsd/1.65 );
    if ( d > ratioDd*p.D() )
    {
        Cvrldv *= sqrt( d / (p.D() * ratioDd) );
    }

    int it = 0, MAX_ITER = 20;

    while( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;
        lambda = ff(Re(result, p.D(), nu), p.D(), p.eps());
        result = A1 * pow<static_rational<1,3> >(pow(1.0 - Cvs/kappa, b) * Cvs * sqrt(musf * Cvb * PI * 0.125) * sqrt(Cvrldv) / lambda ) * sqrt(2.0 * g * Rsd * p.D());
    }

    return result;
}

quantity<velocity> DHLLDV::LDVsliding()
{
    quantity<dimensionless> b, kappa, lambda;

    quantity<velocity> TOL = 1e-6 * meter_per_second, result = 1e6 * TOL, vt, v0 = -result, B;
    quantity<dose_equivalent> C;
    vt = terminalSettlingRuby(nu,d,rhos,rhol,g,shapeFactor);
    b = beta(Re(vt,d,nu));
    kappa = 0.175 * (1.0 + b);

    int it = 0, MAX_ITER = 20;

    while( (abs(v0 - result) > TOL) && (it++ < MAX_ITER) )
    {
        v0 = result;
        lambda = ff(Re(result, p.D(), nu), p.D(), p.eps());

        B = vt * pow(1.0-Cvs/kappa, b) / musf;
        C = pow<2>(1.25*cHe) / musf / lambda * pow<3>(1.0 / sqrtCx(d,nu,rhos,rhol,g)) * pow<static_rational<2,3> >(nu * g);

        result = facTransSbHe * 0.5 * ( B + sqrt(pow<2>(B) + 4.0 * C) );
    }

    return result;
}

quantity<velocity> DHLLDV::LDV()
{
    quantity<dimensionless> Rsd = relativeDensity(rhos,rhol);
    quantity<length> d0 = dTransSmoothRough * (1.65 / Rsd);
    quantity<velocity> result, smooth, rough, trans;

    smooth = LDVsmooth();
    rough = LDVrough();
    trans = LDVsliding();

    if ( smooth < rough )
    {
        result = smooth;
    } else {
        result = smooth * exp(-d / d0) + rough * (1.0 - exp(-d / d0));
    }

    result = (trans > result) ? trans : result;
    result *= pow<static_rational<1, 3> >(cos(p.a()));

    return result;
}








