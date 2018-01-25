#include "fluidutils.h"
#include "homogeneous.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace boost::units;
using namespace boost::units::si;

int main()
{
    quantity<length> Dp = 0.762 * meter;
    quantity<length> d = 2.0/1000.0 * meter;
    quantity<length> eps = 4.5e-05 * meter;
    quantity<mass_density> rhow = 1000 * kilogrammes_per_cubic_metre;
    quantity<acceleration> g = 9.80665 * meter/second/second;
    quantity<kinematic_viscosity> nu = 1.0065122620717291e-6 * meter * meter / second;

    quantity<temperature> T = (20.0 + 273.15) * kelvin;

    quantity<mass_density> rhos = 2650.0 * kilogrammes_per_cubic_metre;

    // Relative solids density
    quantity<dimensionless> Rsd = relativeDensity(rhos, rhow);
    quantity<dimensionless> Cvs = 0.1;

    quantity<mass_density> rhom = Cvs * (rhos - rhow) + rhow;


    Regime *ho = new Homogeneous();

    std::vector< quantity<velocity> > vel;
    for(int i=0; i<=200; i++)
    {
        vel.push_back( i/10.0 * meter_per_second );
    }

    std::vector< quantity<dimensionless> > il;

    for(int i=0; i<=200; i++)
    {
        il.push_back( fluidPressureLoss(vel[i], Dp, eps, nu*rhow, rhow) / (g * rhow) );
    }

    std::ofstream outFile;
    outFile.open("output/ELM.txt");
    for(quantity<dimensionless> dp : il )
    {
        outFile << dp << std::endl;
    }
    outFile.close();

    //== Fixed bed =============================================================


    //== Sliding bed ===========================================================
    std::vector< quantity<dimensionless> > SB_Erhg;
    for(int i = 0; i <= 200; i++)
    {
        SB_Erhg.push_back( musf * meter / pascals );
    }

    std::ofstream SB_file;
    SB_file.open("output/SB_Erhg.txt");
    for( quantity<dimensionless> dp : SB_Erhg )
    {
        SB_file << dp << std::endl;
    }
    SB_file.close();


    //== Heterogeneous =========================================================


    //== ELM ===================================================================


    //== Homogeneous ===========================================================





    return 0;
}
