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
    
    quantity<dimensionless> headLoss( quantity<velocity> v );
    quantity<dimensionless> headLossCvt( quantity<velocity> v );

    quantity<dimensionless> homogeneousHeadLoss( quantity<velocity> v );

    quantity<dimensionless> fixedBedErhg( quantity<velocity> v );
    quantity<dimensionless> fixedBedHeadLoss( quantity<velocity> v );
    dpdx fixedBedPressureLoss(quantity<velocity> v);

    quantity<dimensionless> carrierHeadLoss( quantity<velocity> v );

    quantity<dimensionless> slidingBedErhg( quantity<velocity> v );
    quantity<dimensionless> slidingBedHeadLoss( quantity<velocity> v );

    quantity<dimensionless> heterogeneousErhg( quantity<velocity> v );
    quantity<dimensionless> heterogeneousHeadLoss( quantity<velocity> v );
    dpdx heterogeneousPressureLoss( quantity<velocity> v );

    quantity<acceleration> getGravity() { return g; };
    void setGravity( quantity<acceleration> gravity ) { g = gravity; };

    quantity<dimensionless> getBedConcentration() { return Cvb; };
    void setBedConcentration( quantity<dimensionless> c ) { Cvb = c; };

    quantity<dimensionless> getCHe() { return cHe; };
    void setCHe( quantity<dimensionless> c ) { cHe = c; };

    quantity<dimensionless> getShapeFactor() { return shapeFactor; };
    void setShapeFactor( quantity<dimensionless> f ) { shapeFactor = f; };

    quantity<dimensionless> particleLiftRatio(quantity<velocity> v);
    quantity<dimensionless> HoMobilisationFactor(quantity<velocity> v);
    quantity<dimensionless> HeMobilisationFactor(quantity<velocity> v);
    quantity<dimensionless> ErhgELMFactor(quantity<velocity> v);

    quantity<velocity> LDVsmooth();
    quantity<velocity> LDVrough();
    quantity<velocity> LDVsliding();
    quantity<velocity> LDV();

    quantity<velocity> LSDV(quantity<velocity> v);

    quantity<dimensionless> heterogeneousSlip( quantity<velocity> v );
    quantity<dimensionless> aroundLDVSlip( quantity<velocity> v );
    quantity<dimensionless> threeLMSlip( quantity<velocity> v );
    quantity<dimensionless> fixedBedSlip( quantity<velocity> v );
    quantity<dimensionless> slipTangent(quantity<velocity> v);
    quantity<dimensionless> slipRatio(quantity<velocity> v);

private:
    quantity<acceleration> g = 9.80665 * meters_per_second_squared;
    quantity<dimensionless> musf = 0.416;
    quantity<dimensionless> Cvb = 0.55;
    quantity<dimensionless> cHe = 6.80;
    quantity<dimensionless> CL = 0.270;
    quantity<dimensionless> shapeFactor = 0.77;
    quantity<dimensionless> nTimesThickness = 10.0;
    quantity<dimensionless> Acv = 3.0;
    quantity<dimensionless> alphap = 3.4;
    quantity<dimensionless> ratioDd = 0.015;
    quantity<dimensionless> facTransSbHe = 0.825;
    quantity<length> dTransSmoothRough = 0.0005 * meter;
    quantity<dimensionless> dFraction = 0.1;
    quantity<dimensionless> slipRatioPower = 0.5;

    bool HeHoTransition = true;
    
    quantity<dimensionless> srs(quantity<velocity> v);

    // Friction factor with bed in sliding flow
    quantity<dimensionless> ff_bed(quantity<length> D, quantity<length> d, quantity<velocity> v1, quantity<velocity> v2);
};

#endif
