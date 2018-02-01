#include "limitDepositVelocity.h"
#include "slipRatio.h"
#include "fluidutils.h"
#include "particleUtils.h"
#include "waterData.h"
#include "core/utils.h"
#include "DHLLDV.h"

#include "xlsxwriter.h"

#include <iostream>
#include <string>
#include <algorithm>

#include <boost/units/cmath.hpp>

#define TOL_V (1.0e-14 * meter_per_second)
//#define MAX_ITER 30;

using namespace boost::units;
using namespace boost::units::si;

template <typename Unit, typename T>
std::string unitSymbol(const quantity<Unit, T> &)
{
    return symbol_string(Unit());
}

template <typename T>
void write_parameter(lxw_worksheet* w, int row, int col, T q, std::string desc, lxw_format *format, lxw_format *descfor)
{
    worksheet_write_string(w, row, col  , desc.c_str(), descfor);
    worksheet_write_number(w, row, col+1, q.value(), format);
    worksheet_write_string(w, row, col+2, unitSymbol(q).c_str(), descfor);
}

quantity<velocity> maximumLineSpeed(quantity<length> D, quantity<length> epsilon, quantity<kinematic_viscosity> nu, quantity<acceleration> g)
{
    quantity<velocity> v0 = 10.0 * TOL_V;
    int it = 0;

    quantity<dimensionless> lambda;

    quantity<velocity> vl = 0.0 * meter_per_second, vr = 50.0 * meter_per_second, vc;

    while ( ((vr - vl) > TOL_V) && (it++ < MAX_ITER) )
    {
        vc = 0.5 * (vl + vr);
        lambda = ff( Re(vc, D, nu), D, epsilon );

        if( (0.5*lambda * pow<2>(vc) / g / D) > 1.0 )
        {
            vr = vc;
        }
        else
        {
            vl = vc;
        }
    }

    return 0.5*(vl+vr);
}

int main()
{
    int n = 500;

    DHLLDV *mlm = new DHLLDV();

    lxw_workbook *workbook = workbook_new("../output/imitate.xlsx");
    lxw_worksheet *ws = workbook_add_worksheet(workbook, "Head losses single particle");

    lxw_format *header = workbook_add_format(workbook);
    format_set_bold(header);
    format_set_font_size(header, 10.0);
//    format_set_font_name(header, "Times New Roman");

    lxw_format *pm = workbook_add_format(workbook);
    format_set_font_color(pm, LXW_COLOR_GREEN);
    format_set_bold(pm);
    format_set_align(pm, LXW_ALIGN_CENTER);
    format_set_font_size(pm, 10.0);

    lxw_format *pm_desc = workbook_add_format(workbook);
    format_set_font_size(pm_desc, 10.0);

    lxw_format *clhead = workbook_add_format(workbook);
    format_set_font_color(clhead, LXW_COLOR_GREEN);
    format_set_font_size(clhead, 10.0);
    format_set_bold(clhead);
    format_set_align(clhead, LXW_ALIGN_CENTER);
    format_set_align(clhead, LXW_ALIGN_VERTICAL_CENTER);
    format_set_text_wrap(clhead);

    lxw_format *clnum = workbook_add_format(workbook);
    format_set_font_color(clnum, LXW_COLOR_GREEN);
    format_set_font_size(clnum, 10.0);
    format_set_align(clnum, LXW_ALIGN_CENTER);
    format_set_num_format(clnum, "0.0000");

    lxw_format *dpm = workbook_add_format(workbook);
    format_set_font_color(dpm, LXW_COLOR_RED);
    format_set_font_size(dpm, 10.0);
    format_set_bold(dpm);
    format_set_align(dpm, LXW_ALIGN_CENTER);

    worksheet_set_column_opt(ws, 1, 1, 50.0, NULL, NULL);
    worksheet_set_column_opt(ws, 2, 2, 20.0, NULL, NULL);
    worksheet_set_column_opt(ws, 3, 3, 10.0, NULL, NULL);
    worksheet_set_column_opt(ws, 4, 4, 136.86, NULL, NULL);
    worksheet_set_column_opt(ws, 5, 60, 16.14, NULL, NULL);

//    worksheet_set_row_opt(ws, 2, 17.25, NULL, NULL);
    worksheet_set_row_opt(ws, 0, 49.50, NULL, NULL);
    worksheet_set_row_opt(ws, 1, 99.75, NULL, NULL);

    //--------------------------------------------------------------------------
    // Carrier liquid input parameters
    worksheet_write_string(ws, 2, 1, "Carrier liquid input parameters", header);
    write_parameter(ws, 3, 1, mlm->getRhoLiquid(), "Density of carrier liquid rho_l", pm, pm_desc);
    write_parameter(ws, 4, 1, mlm->getKinematicViscosity(), "Kinematic viscosity nul", pm, pm_desc);

    //--------------------------------------------------------------------------
    // Pipe input parameters
    worksheet_write_string(ws, 6, 1, "Pipe input parameters", header);
    write_parameter(ws, 7, 1, mlm->getPipe().D(), "Pipe diameter D", pm, pm_desc);
    write_parameter(ws, 8, 1, mlm->getPipe().eps(), "Wall roughness epsilon", pm, pm_desc);
    write_parameter(ws, 9, 1, mlm->getPipe().a(), "Pipe inclination angle", pm, pm_desc);

    //--------------------------------------------------------------------------
    // Solids input parameters
    worksheet_write_string(ws, 11, 1, "Solids input parameters", header);
    write_parameter(ws, 12, 1, mlm->getRhoSolids(), "Density of solids", pm, pm_desc);
    // TODO: d50

    write_parameter(ws, 15, 1, -mlm->getBedConcentration() + (quantity<dimensionless>)1.0, "Bed porosity", pm, pm_desc);

    write_parameter(ws, 17, 1, mlm->getSpatialVolumentricConcentration(), "Volumetric concentration Cvs or Cvt", pm, pm_desc);

    //--------------------------------------------------------------------------
    // DHLLDV Framework model parameters
    worksheet_write_string(ws, 22, 1, "DHLLDV Framework model parameters", header);
    write_parameter(ws, 24, 1, mlm->getCHe(), "Homogeneous factor c", pm, pm_desc);


    //--------------------------------------------------------------------------
    // Derived quantities carrier liquid
    worksheet_write_string(ws, 38, 1, "Derived quantities carrier liquid", header);
    write_parameter(ws, 39, 1, relativeDensity(mlm->getRhoSolids(), mlm->getRhoLiquid()), "Relative submerged density", dpm, pm_desc);

    quantity<length> d = 0.001 * meter;
    quantity<velocity> vt = terminalSettlingRuby( mlm->getKinematicViscosity(), d, mlm->getRhoSolids(), mlm->getRhoLiquid(), mlm->getGravity(), 0.77 );
    write_parameter(ws, 40, 1, mlm->getRhoMixture(), "Mixture density", dpm, pm_desc);
    write_parameter(ws, 41, 1, vt, "Terminal settling velocity", dpm, pm_desc);
    write_parameter(ws, 42, 1, sqrtCx(d, mlm->getKinematicViscosity(), mlm->getRhoSolids(), mlm->getRhoLiquid(), mlm->getGravity()), "Durand & Condolios drag coefficient sqrtCx", dpm, pm_desc);
    write_parameter(ws, 43, 1, beta( Re(vt, d, mlm->getKinematicViscosity()) ), "Richardson & Zaki power beta", dpm, pm_desc );
    write_parameter(ws, 44, 1, mlm->LDVsmooth(), "Limit deposit velocity smooth bed", dpm, pm_desc);
    write_parameter(ws, 45, 1, mlm->LDVrough(), "Limit deposit velocity rough bed", dpm, pm_desc);
    write_parameter(ws, 46, 1, mlm->LDVsliding(), "Limit deposit velocity sliding bed minimum", dpm, pm_desc);
    write_parameter(ws, 47, 1, mlm->LDV(), "Resulting limit deposit velocity", dpm, pm_desc);


    //--------------------------------------------------------------------------
    // Derived quantities pseudo liquid
    worksheet_write_string(ws, 55, 1, "Derived quantities pseudo liquid", header);



    //--------------------------------------------------------------------------
    // Additional data required for graphs
    worksheet_write_string(ws, 79, 1, "Additional data required for graphs", header);


    //--------------------------------------------------------------------------
    // Times for the different calculations
    worksheet_write_string(ws, 85, 1, "Times for the different calculations", header);




    //--------------------------------------------------------------------------
    // Line speed
    worksheet_write_string(ws,1,5,"Line speed", clhead);
    worksheet_write_string(ws,0,6,"Carrier liquid", clhead);
    worksheet_write_string(ws,1,6,"Lambda carrier liquid", clhead);
    worksheet_merge_range(ws,0,7,0,9,"Carrier liquid without fines", clhead);
    worksheet_write_string(ws,1,7,"Hydraulic gradient carrier liquid", clhead);
    worksheet_write_string(ws,1,8,"Hydraulic gradient carrier ELM liquid", clhead);
    worksheet_write_string(ws,1,9,"Solids effect ELM carrier liquid", clhead);

    worksheet_write_string(ws,1,19,"Hydraulic gradient FB carrier liquid", clhead);

    worksheet_write_string(ws,1,26,"Hydraulic gradient SB carrier liquid", clhead);

    worksheet_write_string(ws,1,33,"Hydraulic gradient He carrier liquid", clhead);

    worksheet_write_string(ws,1,40,"Hydraulic gradient Ho carrier liquid", clhead);

    worksheet_write_string(ws,1,54,"Heterogeneous slip curve", clhead);

    // Maximum linespeed where the hydraulic gradient is 1 for the carrier liquid
    quantity<velocity> vmax = maximumLineSpeed(mlm->getPipe().D(), mlm->getPipe().eps(),
                                               mlm->getKinematicViscosity(), mlm->getGravity() );
    quantity<velocity> dv = vmax/((double)n), vtemp;
    quantity<dimensionless> lambda;

    for(int i = 0; i <= n; i++)
    {
        vtemp = ((double)i)*dv;
        vtemp = std::max( vtemp, 0.001*meter_per_second );
        worksheet_write_number(ws, 3+i, 5, vtemp.value(), clnum);

        lambda = ff( Re(vtemp, mlm->getPipe().D(), mlm->getKinematicViscosity() ), mlm->getPipe().D(), mlm->getPipe().eps() );
        worksheet_write_number(ws, 3+i, 6, lambda.value(), clnum);
        worksheet_write_number(ws, 3+i, 7, mlm->carrierHeadLoss(vtemp), clnum);
        worksheet_write_number(ws, 3+i, 8, mlm->carrierHeadLoss(vtemp)*mlm->getRhoMixture()/mlm->getRhoLiquid(), clnum );

        worksheet_write_number(ws, 3+i, 19, mlm->fixedBedHeadLoss(vtemp), clnum);

        worksheet_write_number(ws, 3+i, 26, mlm->slidingBedHeadLoss(vtemp), clnum);

        worksheet_write_number(ws, 3+i, 33, mlm->heterogeneousHeadLoss(vtemp), clnum);

        worksheet_write_number(ws, 3+i, 40, mlm->homogeneousHeadLoss(vtemp), clnum);


        worksheet_write_number(ws, 3+i, 54, 0.0, clnum);
    }

    // Dimension headers of the columns
    worksheet_write_string(ws,2,5, unitSymbol(vmax).c_str(), clhead);
    worksheet_write_string(ws,2,6, unitSymbol(lambda).c_str(), clhead);
    worksheet_write_string(ws,2,7, unitSymbol(mlm->homogeneousHeadLoss(vmax)).c_str(), clhead);
    worksheet_write_string(ws,2,8, unitSymbol(mlm->homogeneousHeadLoss(vmax)*mlm->getRhoMixture()/mlm->getRhoLiquid()).c_str(), clhead);

    worksheet_write_string(ws,2,19, unitSymbol(mlm->fixedBedHeadLoss(vmax)).c_str(), clhead);

    worksheet_write_string(ws,2,26, unitSymbol(mlm->slidingBedHeadLoss(vmax)).c_str(), clhead);

    worksheet_write_string(ws,2,33, unitSymbol(mlm->heterogeneousHeadLoss(vmax)).c_str(), clhead);

    worksheet_write_string(ws,2,40, unitSymbol(mlm->homogeneousHeadLoss(vmax)).c_str(), clhead);

    worksheet_write_string(ws,2,54, "", clhead);

    //--------------------------------------------------------------------------
    // Charts
    lxw_chart *cl_chart = workbook_add_chart(workbook, LXW_CHART_SCATTER_STRAIGHT);
    lxw_chart_series *series;
    lxw_image_options ch_options = { .x_offset = 0, .y_offset = 0, .x_scale = 2.0160018, .y_scale = 1.8 };

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$H$4:$H$503");
    chart_series_set_name(series, "Carrier liquid");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$I$4:$I$503");
    chart_series_set_name(series, "ELM carrier liquid");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$T$4:$T$503");
    chart_series_set_name(series, "Fixed bed regime");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$AA$4:$AA$503");
    chart_series_set_name(series, "Sliding bed regime");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$AH$4:$AH$503");
    chart_series_set_name(series, "Heterogeneous regime");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$AO$4:$AO$503");
    chart_series_set_name(series, "Homogeneous regime");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$AV$4:$AV$503");
    chart_series_set_name(series, "Resulting curve Cvs = c");

    series = chart_add_series(cl_chart, "='Head losses single particle'!$F$4:$F$503", "='Head losses single particle'!$BI$4:$BI$503");
    chart_series_set_name(series, "Resulting curve Cvt = c");

    chart_axis_set_min(cl_chart->x_axis, 0.0);
    chart_axis_set_max(cl_chart->x_axis, 10.0);

    chart_axis_set_min(cl_chart->y_axis, 0.0);
    chart_axis_set_max(cl_chart->y_axis, 0.4);

    chart_axis_major_gridlines_set_visible(cl_chart->x_axis, LXW_TRUE);
    chart_axis_set_minor_unit(cl_chart->x_axis, 0.1);
    chart_axis_minor_gridlines_set_visible(cl_chart->x_axis, LXW_TRUE);
    chart_axis_set_minor_unit(cl_chart->y_axis, 0.01);
    chart_axis_minor_gridlines_set_visible(cl_chart->y_axis, LXW_TRUE);

    chart_axis_set_name(cl_chart->x_axis, "Line speed [m/s]");
    chart_axis_set_name(cl_chart->y_axis, "Carrier liquid based hydraulic gradient [m/m]");

    chart_title_set_name(cl_chart, "Carrier liquid based hydraulic gradients without fines");
    lxw_chart_font font = { name: NULL, size: 10, bold: LXW_FALSE, italic: LXW_FALSE, underline: LXW_FALSE, rotation: 0, color: LXW_COLOR_BLACK };
    chart_title_set_name_font(cl_chart, &font);

    worksheet_insert_chart_opt(ws,2,4,cl_chart, &ch_options);

    return workbook_close(workbook);
}
