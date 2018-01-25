#include "xlsxwriter.h"

#include "limitDepositVelocity.h"
#include "slipRatio.h"
#include "fluidutils.h"
#include "particleUtils.h"
#include "waterData.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/units/cmath.hpp>

using namespace DHLLDV::LDV;

int main()
{
    // Input -------------------------------------------------------------------

    // Temperature
    quantity<celsius::temperature> T = 20.0 * celsius::degree;

    // Liquid density
    quantity<mass_density> rhow = waterDensity(T);
    rhow = 1025.0 * kilogrammes_per_cubic_metre;

    // Liquid viscosity
    quantity<kinematic_viscosity> nu = waterKinematicViscosity(T);
    nu = 0.0000013 * meter * meter_per_second;

    // Pipe diameter
    quantity<length> D = 0.762 * meter;
    D = 0.0254 * meter;

    // Pipe roughness
    quantity<length> eps = 4.5e-05 * meter;
    eps = 0.000001 * meter;

    // Particle diameter
    quantity<length> d = 0.001 * meter;
//    quantity<length> d = 0.000355 * meter;

    // Particle density
    quantity<mass_density> rhos = 2650.0 * kilogrammes_per_cubic_metre;

    // Spatial volumetric concentration
    quantity<dimensionless> Cvs = 0.175;
    quantity<dimensionless> Cvb = 0.6;
    Cvb = 0.55;

    // Sliding friction coefficient
    quantity<dimensionless> musf = 0.416;

    // Gravitational acceleration
    quantity<acceleration> g = 9.80665 * meters_per_second_squared;

    lxw_workbook *wb = workbook_new("2LM-3LM-Slager.xlsx");
    lxw_worksheet *ws1 = workbook_add_worksheet(wb, "Manual");
    lxw_worksheet *ws2 = workbook_add_worksheet(wb, "Calculations");
    lxw_worksheet *ws3 = workbook_add_worksheet(wb, "Graphs");
    lxw_worksheet *ws4 = workbook_add_worksheet(wb, "Comparison");

    std::stringstream ss;
    ss << "Dp = " << D;

    worksheet_write_string(ws2, 0, 0, (ss.str()).c_str(), NULL);

    workbook_close(wb);

/*    // Create a workbook and add a worksheet.
    lxw_workbook  *workbook  = workbook_new("tutorial01.xlsx");
    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, NULL);
    // Start from the first cell. Rows and columns are zero indexed.
    int row = 0;
    int col = 0;
    // Iterate over the data and write it out element by element.
    for (row = 0; row < 4; row++) {
        worksheet_write_string(worksheet, row, col,     expenses[row].item, NULL);
        worksheet_write_number(worksheet, row, col + 1, expenses[row].cost, NULL);
    }
    // Write a total using a formula.
    worksheet_write_string (worksheet, row, col,     "Total",       NULL);
    worksheet_write_formula(worksheet, row, col + 1, "=SUM(B1:B4)", NULL);
    // Save the workbook and free any allocated memory.
    return workbook_close(workbook);*/


/*
    std::ofstream out;
    out.open("output/LDV.txt");

    // quantities which are shared between methods are calculated first
    quantity<dimensionless> Rsd = relativeDensity(rhos,rhow);
    // TODO: change factor to 3.2 for production
    quantity<dimensionless> ap = 3.4 * pow<static_rational<2,9> >(1.65 / Rsd);
    quantity<velocity> vt;
    quantity<dimensionless> Re_p;
    quantity<dimensionless> beta;
    quantity<dimensionless> KC;

    for (int i = -500; i <= -100; i = i + 7)
    {
        d = pow(10.0, 1.0*i/100.0) * meter;
        vt = terminalSettlingRuby(nu, d, rhos, rhow, g);
        Re_p = vt * d / nu;
        beta = ( 4.7 + 0.41 * pow<static_rational<3,4> >(Re_p) ) / ( 1.0 + 0.175 * pow<static_rational<3,4> >(Re_p) );
        KC = 0.175*(1.0 + beta);

        out << d.value() << ", "
            << (verySmallParticles(nu,rhos,rhow,D,eps,g)).value() << ", "
            << (smallParticles(nu,Cvs,Rsd,ap,vt,beta,KC,D,eps,d,g)).value() << ", "
            << (largeParticles(nu,Cvs,Rsd,ap,vt,beta,KC,D,eps,d,g,musf,Cvb)).value() << ", "
            << (limitDepositVelocity(nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << ", "
            << (lowerLimit(nu,D,d,eps,vt,Cvs,beta,KC,g,musf)).value() << std::endl;
    }

    out.close();
*/

    /*std::ofstream out;
    out.open("output/LDV_0.1524m.txt");

    for (int i = -500; i <= -100; i++)
    {
        d = pow(10.0, 1.0*i/100.0) * meter;
        out << d.value() << ", " << (limitDepositVelocity(nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << std::endl;
    }

    out.close();


    D = 0.762 * meter;
    out.open("output/LDV_0.762m.txt");

    for (int i = -500; i <= -100; i++)
    {
        d = pow(10.0, 1.0*i/100.0) * meter;
        out << d.value() << ", " << (limitDepositVelocity(nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << std::endl;
    }

    out.close();

    D = 0.3 * meter;
    out.open("output/LDV_0.3m.txt");

    for (int i = -500; i <= -100; i++)
    {
        d = pow(10.0, 1.0*i/100.0) * meter;
        out << d.value() << ", " << (limitDepositVelocity(nu, Cvs, rhos, rhow, D, eps, d, g, musf, Cvb)).value() << std::endl;
    }

    out.close();


    std::cout << DHLLDV::slipRatio(5.0*meter_per_second, D, d, eps, nu, rhow, rhos, Cvs, g, musf, Cvb) << std::endl;*/

    return 0;
}


/*#include "xlsxwriter.h"
// Some data we want to write to the worksheet.
struct expense {
    char item[32];
    int  cost;
};
struct expense expenses[] = {
    {"Rent", 1000},
    {"Gas",   100},
    {"Food",  300},
    {"Gym",    50},
};
int main() {
    // Create a workbook and add a worksheet.
    lxw_workbook  *workbook  = workbook_new("tutorial01.xlsx");
    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, NULL);
    // Start from the first cell. Rows and columns are zero indexed.
    int row = 0;
    int col = 0;
    // Iterate over the data and write it out element by element.
    for (row = 0; row < 4; row++) {
        worksheet_write_string(worksheet, row, col,     expenses[row].item, NULL);
        worksheet_write_number(worksheet, row, col + 1, expenses[row].cost, NULL);
    }
    // Write a total using a formula.
    worksheet_write_string (worksheet, row, col,     "Total",       NULL);
    worksheet_write_formula(worksheet, row, col + 1, "=SUM(B1:B4)", NULL);
    // Save the workbook and free any allocated memory.
    return workbook_close(workbook);
}

*/
