#include <iostream>
#include <fstream>
#include <vector>

#include "DHLLDV.h"
#include "waterData.h"

int main()
{
    quantity<celsius::temperature> T = 20.0 * celsius::degree;

    DHLLDV *mlm = new DHLLDV();

    mlm->setRhoLiquid( waterDensity(T) );
    mlm->setKinematicViscosity( waterKinematicViscosity(T) );

    std::vector<quantity<velocity> > v;

    for(int i = 0; i <= 300; i++)
    {
        v.push_back( 1.0 * i / 10.0 * meter_per_second );
    }


    // Equivalent liquid model
    std::ofstream out;
    out.open("../output/fig8.16.txt");

    mlm->fixedBedPressureLoss(5.0*meter_per_second);

    for( auto vs : v )
    {
        out << vs.value() << ", "
            << (mlm->homogeneousHeadLoss(vs)).value() << ", "       // Equivalent liquid model
            << (mlm->fixedBedErhg(vs)).value() << ", "              // Fixed bed
            << (mlm->fixedBedHeadLoss(vs)).value() << ", "          // Fixed bed
            << (mlm->slidingBedErhg(vs)).value() << ", "            // Sliding bed
            << (mlm->slidingBedHeadLoss(vs)).value() << ", "                // Sliding bed
            << (mlm->heterogeneousErhg(vs)).value() << ", "
            << (mlm->heterogeneousHeadLoss(vs)).value()
            << std::endl;
    }
    out.close();



    return 0;
}

