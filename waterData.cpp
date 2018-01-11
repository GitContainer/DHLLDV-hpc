#include "waterData.h"

#include <boost/units/cmath.hpp>

#include <iostream>

quantity<mass_density> waterDensity( quantity<celsius::temperature> T )
{
    // Formula from view-source:https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html

/*    msg +=  + Math.round((-0.00000000009204453627*vA*vA*vA*vA
 * +0.00000003420742008672*vA*vA*vA
 * -0.00000708919807166417*vA*vA
 * +0.00004375294545181970*vA
 * +0.99988826440573500000
 * )*100000)/100000 + " g/cm3\n";
    */

    return (0.00000009204453627*pow<4>(T/(celsius::degree)) + 0.00003420742008672*pow<3>(T/(celsius::degree)) - 0.00708919807166417*pow<2>(T/(celsius::degree))
             + 0.04375294545181970 * T/(celsius::degree) + 999.88826440573500000) * kilogrammes_per_cubic_metre;
}

quantity<kinematic_viscosity> waterKinematicViscosity( quantity<celsius::temperature> T )
{
    // Formula from view-source:https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
/*     msg +=  + Math.round((0.00000000000282244333*vA*vA*vA*vA*vA*vA
 * -0.00000000126441088087*vA*vA*vA*vA*vA
 * +0.00000023336659710795*vA*vA*vA*vA
 * -0.0000234079044336466*vA*vA*vA
 * +0.00144686943485654*vA*vA
 * -0.0607310297913931*vA
 * +1.79194000343777)*0.000001*10000000000)/10000000000 + " m2/s\n";
*/
    return ( 0.00000000000282244333 * pow<6>(T/celsius::degree) - 0.00000000126441088087 * pow<5>(T/celsius::degree) + 0.00000023336659710795 * pow<4>(T/celsius::degree) - 0.0000234079044336466 * pow<3>(T/celsius::degree)
                   + 0.00144686943485654 * pow<2>(T/celsius::degree) - 0.0607310297913931 * T/celsius::degree + 1.79194000343777 ) * 1e-6 * meter * meter_per_second;
}

