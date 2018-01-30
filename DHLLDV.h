#ifndef DHLLDV_H
#define DHLLDV_H

#include "mixtureLossModel.h"
#include "fluidutils.h"
#include "particleUtils.h"

class DHLLDV : public mixtureLossModel
{
public:
    DHLLDV();
    ~DHLLDV();

    dpdx pressureLoss(quantity<velocity> v);

    quantity<dimensionless> homogeneousHeadLoss( quantity<velocity> v );

    quantity<dimensionless> fixedBedErhg( quantity<velocity> v );
    quantity<dimensionless> fixedBedHeadLoss( quantity<velocity> v );
    dpdx fixedBedPressureLoss(quantity<velocity> v);

    quantity<dimensionless> slidingBedErhg( quantity<velocity> v );
    quantity<dimensionless> slidingBedHeadLoss( quantity<velocity> v );

    quantity<dimensionless> heterogeneousErhg( quantity<velocity> v );
    quantity<dimensionless> heterogeneousHeadLoss( quantity<velocity> v );
    dpdx heterogeneousPressureLoss( quantity<velocity> v );

    quantity<acceleration> getGravity() { return g; };
    void setGravity( quantity<acceleration> gravity ) { g = gravity; };

private:
    quantity<acceleration> g = 9.80665 * meters_per_second_squared;
    quantity<dimensionless> musf = 0.416;
    quantity<dimensionless> Cvb = 0.6;


};

#endif
