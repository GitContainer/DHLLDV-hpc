#include "limitDepositVelocity.h"
#include "slipRatio.h"
#include "fluidutils.h"
#include "particleUtils.h"
#include "waterData.h"

#include "xlsxwriter.h"

#include <iostream>
#include <string>

#include <boost/units/cmath.hpp>

#define TOL_V (1.0e-10 * meter_per_second)
//#define MAX_ITER 30;

template <typename T>
void write_parameter(lxw_worksheet* w, int row, int col, T q, std::string desc)
{
    worksheet_write_string(w, row, col  , desc.c_str(), NULL);
    worksheet_write_number(w, row, col+1, q.value(), NULL);
    worksheet_write_string(w, row, col+2, "TODO", NULL);    // TODO: get symbol from quantity
}

quantity<velocity> maximumLineSpeed(quantity<length> D, quantity<length> epsilon, quantity<kinematic_viscosity> nu)
{
    quantity<velocity> v0 = 10.0 * TOL_V, result = -v0;
    int it = 0;

    quantity<dimensionless> lambda;
    quantity<dimensionless> head;

    while ( (abs(v0 - result) > TOL_V) && (it++ < MAX_ITER) )
    {
//        lambda = ff( Re(result, D, nu) );

    }

    return result;
}

/*
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
*/

int main()
{
    lxw_workbook *workbook = workbook_new("../output/imitate.xlsx");
    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, "Head losses single particle");

    //--------------------------------------------------------------------------
    // Carrier liquid input parameters
    worksheet_write_string(worksheet, 2, 1, "Carrier liquid input parameters", NULL);

    // Density of carrier fluid
    quantity<mass_density> rhol = 1025.0 * kilogrammes_per_cubic_metre;
    write_parameter(worksheet, 3, 1, rhol, "Density of carrier liquid rho_l");

    // Kinematic viscosity of carrier fluid
    quantity<kinematic_viscosity> nu = 0.0000013 * meter * meter_per_second;
    write_parameter(worksheet, 4, 1, nu, "Kinematic viscosity nul");

    //--------------------------------------------------------------------------
    // Pipe input parameters
    worksheet_write_string(worksheet, 6, 1, "Pipe input parameters", NULL);

    quantity<length> D = 0.7620 * meter;
    write_parameter(worksheet, 7, 1, D, "Pipe diameter D");

    quantity<length> eps = 0.000001 * meter;
    write_parameter(worksheet, 8, 1, eps, "Wall roughness epsilon");

    //quantity<angle> incl_angle = 0.0 * radians; // TODO
    //write_parameter(worksheet, 9, 1, incl_angle, "Pipe inclination angle"); // TODO



    //--------------------------------------------------------------------------
    // Line speed
    worksheet_write_string(worksheet,1,5,"Line speed", NULL);
    worksheet_write_string(worksheet,2,5,"m/s", NULL);

    // Maximum linespeed where the hydraulic gradient is 1 for the carrier liquid
    quantity<velocity> vmax = maximumLineSpeed(D, eps, nu);

    std::cout << vmax << std::endl;

    return workbook_close(workbook);
}
