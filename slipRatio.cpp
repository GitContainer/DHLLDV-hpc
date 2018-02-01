#include "slipRatio.h"
#include "particleUtils.h"
#include "fluidutils.h"
#include "limitDepositVelocity.h"

#include "boost/units/cmath.hpp"
#include <algorithm>

#include <iostream>

quantity<dimensionless> dhlldv::slipRatio( quantity<velocity> v, quantity<length> D, quantity<length> d, quantity<length> eps, quantity<kinematic_viscosity> nu, quantity<mass_density> rhol, quantity<mass_density> rhos,
                                           quantity<dimensionless> Cvt, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb )
{
    quantity<dimensionless> Rsd = relativeDensity(rhos, rhol);
    quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g, 0.77);
    quantity<dimensionless> Re_p = Re(vt, d, nu);

    quantity<dimensionless> beta = (4.7 + 0.41 * pow<static_rational<3,4> >(Re_p)) / ( 1.0 + 0.175 * pow<static_rational<3,4> >(Re_p) );
    quantity<dimensionless> KC = 0.175 * (1.0 + beta);

    quantity<velocity> vls_ldv = LDV::limitDepositVelocity(nu, Cvt, rhos, rhol, D, eps, d, g, musf, Cvb);
    quantity<dimensionless> CD = 4.0/3.0 * g*Rsd*d/pow<2>(vt);

    quantity<dimensionless> Xi_ldv = 1.0/(2.0*CD) * pow(1.0-Cvt/KC, beta) * vls_ldv / v;
    quantity<velocity> vs_ldv = vls_ldv * Xi_ldv;

    quantity<dimensionless> kappa = 1.0 / (1.0 - Xi_ldv);
    quantity<dimensionless> Xi_fb = 1.0 - Cvt * vls_ldv / (((Cvb - kappa * Cvt)*(vls_ldv - v)) + kappa * Cvt * vls_ldv);
    quantity<dimensionless> Xi_th = (Xi_fb < Xi_ldv) ? Xi_fb : Xi_ldv;

    quantity<velocity> vls_t = vls_ldv * pow( 2.5 / CD * pow(1-Cvt/KC,beta) * pow(1-Cvt/Cvb,-1), 0.25);
    quantity<dimensionless> Xi_t = (1.0 - Cvt/Cvb) * (1 - 0.8 * v/vls_t);

    double alpha_xi = 0.0;

    quantity<dimensionless> Xi = Xi_th * pow(1.0 - v/vls_t, alpha_xi) + Xi_t * pow(v/vls_t, alpha_xi);
    quantity<dimensionless> compare = 0.0;
    Xi = std::max(Xi, compare);
    compare = 1.0;
    Xi = std::min(Xi, compare);

    return Xi;
}

quantity<dimensionless> dhlldv::slipRatio_FB(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvt, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb)
{
    quantity<velocity> v_ldv = LDV::limitDepositVelocity(nu, Cvt, rhos, rhol, D, eps, d, g, musf, Cvb);
    quantity<dimensionless> kappa = 1.0 / ( 1.0 - slipRatio_LDV(v, nu, Cvt, rhos, rhol, D, eps, d, g, musf, Cvb) );

    quantity<dimensionless> result = 1.0 - Cvt * v_ldv / (((Cvb - kappa * Cvt)*(v_ldv - v)) + kappa * Cvt * v_ldv);
    if (result < 0.0)
    {
        result = 0.0;
    }

    return result;

}

/*
def slip_ratio(vls, Dp,  d, epsilon, nu, rhol, rhos, Cvt):
    """
    Return the slip ratio (Xi) for the given slurry.
    Dp = Pipe diameter (m)
    d = Particle diameter (m)
    epsilon = absolute pipe roughness (m)
    nu = fluid kinematic viscosity in m2/sec
    rhol = density of the fluid (ton/m3)
    rhos = particle density (ton/m3)
    Cvt = transport volume concentration
    """
    Rsd = (rhos-rhol)/rhol
    vt = heterogeneous.vt_ruby(d, Rsd, nu)  # particle shape factor assumed for sand for now
    CD = (4/3.)*((gravity*Rsd*d)/vt**2)      # eqn 4.4-6 without the shape factor
    Rep = vt*d/nu  # eqn 4.2-6
    top = 4.7 + 0.41*Rep**0.75
    bottom = 1. + 0.175*Rep**0.75
    beta = top/bottom  # eqn 4.6-4
    KC = 0.175*(1+beta)
    vls_ldv = LDV(vls, Dp, d, epsilon, nu, rhol, rhos, Cvt)
    Xi_ldv = (1/(2*CD)) * (1-Cvt/KC)**beta * (vls_ldv/vls)  # eqn 8.12-1

    vs_ldv = vls_ldv * Xi_ldv
    Kldv = 1/(1 - Xi_ldv)   # eqn 7.9-14
    Xi_fb = 1-((Cvt*vs_ldv)/(stratified.Cvb-Kldv*Cvt)*(vs_ldv-vls)+Kldv*Cvt*vs_ldv)  # eqn 8.12-2
    Xi_th = min(Xi_fb, Xi_ldv)  # eqn 8.12-3

    vls_t = vls_ldv*(5 * (1/(2*CD)) * (1-Cvt/KC)**beta * (1-Cvt/stratified.Cvb)**(-1))**(1./4) # eqn8.12-4
    Xi_t = (1-Cvt/stratified.Cvb) * (1-(4./5)*(vls/vls_t))  # 8.12-5

    Xi = Xi_th*(1-(vls/vls_t)**alpha_xi) + Xi_t *(vls/vls_t)**alpha_xi    # eqn 8.12-7
    return min(max(Xi, 0.0),1.0)
*/

quantity<dimensionless> dhlldv::slipRatio_LDV(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvt, quantity<mass_density> rhos, quantity<mass_density> rhol,
                                              quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb)
{
    quantity<dimensionless> Cvr = Cvt / Cvb;
    quantity<dimensionless> alpha = 0.58 * pow<static_rational<-42, 100> >( Cvr );

    quantity<velocity> v_ldv = LDV::limitDepositVelocity(nu, Cvt, rhos, rhol, D, eps, d, g, musf, Cvb);

    quantity<dimensionless> power =
            -(0.83 + musf/4.0 + pow<2>(Cvr - 0.5 - 0.075*D/meter) + 0.025 * D/meter ) * pow<static_rational<25,1000> >(D/meter) * pow(v/v_ldv, alpha) * pow<static_rational<65,100> >(Cvr);

    return (1.0 - Cvr) * exp(power) * pow<4>( v_ldv / v );

    /*quantity<velocity> vt = terminalSettlingRuby(nu,d,rhos,rhol,g);
    quantity<dimensionless> Cd = (4.0/3.0)*( g * relativeDensity(rhos,rhol) * d / pow<2>(vt) );

    quantity<dimensionless> Re_p = reynolds(vt,d,nu);
    quantity<dimensionless> beta = (4.7 + 0.41*pow(Re_p,0.75)) / (1.0 + 0.175*pow(Re_p,0.75));
    quantity<dimensionless> Kc = 0.175*(1.0 + beta);

    quantity<velocity> v_ldv = LDV::limitDepositVelocity(nu, Cvt, rhos, rhol, D, eps, d, g, musf, Cvb);

    return (1.0/(2.0*Cd) * pow(1.0-Cvt/Kc, beta) * v_ldv / v );*/
}

quantity<dimensionless> dhlldv::slipRatio_HeHo(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<length> D, quantity<length> d, quantity<length> eps, quantity<mass_density> rhos, quantity<mass_density> rhol,
                                               quantity<acceleration> g)
{
    /*quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g);

    return (8.5 / sqrt(ff(Re(v, D, nu), D, eps)) * pow<static_rational<5,3> >(vt / sqrt(g * d)) * pow<static_rational<1,3> >(nu * g) * vt / pow<2>(v));*/

    quantity<plane_angle> ia = 0.0 * radians;

    quantity<velocity> vt = terminalSettlingRuby(nu, d, rhos, rhol, g, 0.77);
    quantity<velocity> vstar = sqrt( ff(Re(v,D,nu),D,eps)/8.0) * v;
    if ( (5.8*vstar) < vt*sin(ia) )
    {
        vstar = vt * sin(ia) / 5.8;
    }

    quantity<dimensionless> corrFactor = pow<static_rational<4,3> >(vt * cos(ia) / (11.6 * vstar - vt * sin(ia)) / vt * (11.6 * vstar));

    return corrFactor * pow<2>(1.25*6.80) / ff(Re(v,D,nu),D,eps) / pow<3>( sqrtCx(d,nu,rhos,rhol,g) ) * pow<2>(pow<static_rational<1,3> >(nu*g) / v);
}

quantity<dimensionless> dhlldv::slipRatio_tan(quantity<velocity> v, quantity<kinematic_viscosity> nu, quantity<dimensionless> Cvt, quantity<mass_density> rhos, quantity<mass_density> rhol, quantity<length> D, quantity<length> eps, quantity<length> d, quantity<acceleration> g, quantity<dimensionless> musf, quantity<dimensionless> Cvb)
{
    quantity<dimensionless> Cvr = Cvt / Cvb;
    quantity<dimensionless> alpha = 0.58 * pow<static_rational<-42, 100> >( Cvr );

    quantity<velocity> v_ldv = LDV::limitDepositVelocity(nu,Cvt,rhos,rhol,D,eps,d,g,musf,Cvb);

    quantity<dimensionless> power =
            -( (0.83 + musf/4.0 + pow<2>(Cvr - 0.5 - 0.075*D/meter) + 0.025 * D/meter ) ) * pow<static_rational<25,1000> >(D/meter) * pow(v/v_ldv, alpha) * pow<static_rational<65,100> >(Cvr);

    quantity<velocity> vlst = pow<static_rational<1,4> >( 5.0 * exp(power) ) * v_ldv;

    return ( (1.0- Cvt / Cvb)*( 1.0 - 0.8*(v / vlst) ) );
}

/*
def slip_ratio(vls, Dp,  d, epsilon, nu, rhol, rhos, Cvt):
    """
    Return the slip ratio (Xi) for the given slurry.
    Dp = Pipe diameter (m)
    d = Particle diameter (m)
    epsilon = absolute pipe roughness (m)
    nu = fluid kinematic viscosity in m2/sec
    rhol = density of the fluid (ton/m3)
    rhos = particle density (ton/m3)
    Cvt = transport volume concentration
    """
    Rsd = (rhos-rhol)/rhol
    vt = heterogeneous.vt_ruby(d, Rsd, nu)  # particle shape factor assumed for sand for now
    CD = (4/3.)*((gravity*Rsd*d)/vt**2)      # eqn 4.4-6 without the shape factor
    Rep = vt*d/nu  # eqn 4.2-6
    top = 4.7 + 0.41*Rep**0.75
    bottom = 1. + 0.175*Rep**0.75
    beta = top/bottom  # eqn 4.6-4
    KC = 0.175*(1+beta)
    vls_ldv = LDV(vls, Dp, d, epsilon, nu, rhol, rhos, Cvt)
    Xi_ldv = (1/(2*CD)) * (1-Cvt/KC)**beta * (vls_ldv/vls)  # eqn 8.12-1

    vs_ldv = vls_ldv * Xi_ldv
    Kldv = 1/(1 - Xi_ldv)   # eqn 7.9-14
    Xi_fb = 1-((Cvt*vs_ldv)/(stratified.Cvb-Kldv*Cvt)*(vs_ldv-vls)+Kldv*Cvt*vs_ldv)  # eqn 8.12-2
    Xi_th = min(Xi_fb, Xi_ldv)  # eqn 8.12-3

    vls_t = vls_ldv*(5 * (1/(2*CD)) * (1-Cvt/KC)**beta * (1-Cvt/stratified.Cvb)**(-1))**(1./4) # eqn8.12-4
    Xi_t = (1-Cvt/stratified.Cvb) * (1-(4./5)*(vls/vls_t))  # 8.12-5

    Xi = Xi_th*(1-(vls/vls_t)**alpha_xi) + Xi_t *(vls/vls_t)**alpha_xi    # eqn 8.12-7
return min(max(Xi, 0.0),1.0)*/
