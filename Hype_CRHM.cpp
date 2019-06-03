// 05/31/19
//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "Hype_CRHM.h"
#include "NewModules.h"
#include "DefCRHMGlobal.h"
#include <sstream>

#include <math.h>
#include <stdlib.h>
//---------------------------------------------------------------------------
#pragma package(smart_init)

//---------------------------------------------------------------------------

using namespace std;
using namespace CRHM;

extern float fLimit;

//---------------------------------------------------------------------------

// const long i_in = 0; // inorganic nitrogen
// const long i_on = 1; // organic nitrogen
// const long i_sp = 2; // soluble (reactive) phosphorus, i.e. phosphate
// const long i_pp = 3; // particulate phosphorus
// const long i_oc = 4; // (dissolved) organic carbon

// kg/km2 = (mg/l)*mm = mg/m2
// Kg == mm
// mg/l == Kg/Km2 = 1/mm
// l/m2 == mm = Kg/m2
//==============================================================================

ClassWQ_Soil* ClassWQ_Soil::klone(string name) const{
  return new ClassWQ_Soil(name);
}

void ClassWQ_Soil::decl(void) {

  Description = "'Handles soil moisture throughout the year.' \
                 'Standard version,' \
                 'Version with Culvert limited runoff.' \
                 'Version with Tile drainage calculated ssr.'";

  variation_set = VARIATION_1;

  declvar("culvert_Q", NHRU, "flow in culvert.", "(m^3/s)", &culvert_Q);

  declvar("culvert_water_H", NHRU, "depth of pond at culvert inlet.", "(m)", &culvert_water_H);

  declvar("culvert_water_A", NHRU, "surface area of culvert pond.", "(m^2)", &culvert_water_A);

  declvar("culvert_water_V", NHRU, "volume of water in culvert pond.", "(m^3)", &culvert_water_V);

  declvar("culvert_over_Q", NHRU, "flow over the road.", "(m^3/s)", &culvert_over_Q);

  declvar("culvert_evap", NHRU, "Depth of water evaporating from culvert pond.", "(mm/int)", &culvert_evap);

  declstatdiag("cum_culvert", NHRU, "Cumulative culvert HRU flow.", "(m^3)", &cum_culvert);

  declstatdiag("cum_culvert_over", NHRU, "Cumulative culvert HRU overflow.", "(m^3)", &cum_culvert_over);

  decldiag("HD", NHRU, "ratio of depth of water at culvert/culvert diameter.", "()", &HD);

  declparam("channel_slope", NHRU, "[0.002]", "0.0001", "0.01", "soil slope to culvert.", "()", &channel_slope);

  declparam("side_slope", NHRU, "[0.02]", "0.0001", "0.01", "side soil slope mormal to culvert slope.", "()", &side_slope);

  declparam("culvert_diam", NHRU, "[0.5]", "0.1", "5.0", "culvert diameter.", "(m)", &culvert_diam);

  declparam("culvert_water_Dmax", NHRU, "[2.0]", "0.1", "10.0", "maximum depth of pond at culvert inlet.", "(m)", &culvert_water_Dmax);

  declparam("number_culverts", NHRU, "[1.0]", "0.0", "10.0", "number of culverts and efficiency factor. Zero = no culvert.", "()", &number_culverts);

  declparam("culvert_type", NHRU, "[0]", "0", "4", "0- thin walled projection, 1- square edged (flush) inlet, 2- socket and concrete pipe, 3- 45 degree beveled inlet, 4- well-rounded (streamlined) inlet.", "()", &culvert_type);


  variation_set = VARIATION_2;

  declvar("tile_flow", NHRU, "flow in tile drainage.", "(m^3/s)", &tile_flow);

  declvar("Dequiv_depth", NHRU, "closed-form expression for equivalent depth.", "(m)", &Dequiv_depth);

  declvar("x", NHRU, ".", "(m)", &x);

  declvar("Dw", NHRU, "steady state depth of the watertable midway between the drain.", "(m)", &Dw);

  declparam("Ka", NHRU, "[1]", "0.1", "4", "hydraulic conductivity of the soil above drain level.", "(m/d)", &Ka);

  declparam("Kb", NHRU, "[1]", "0.1", "4", "hydraulic conductivity of the soil below drain level.", "(m/d)", &Kb);

  declparam("Di", NHRU, "[2]", "0", "10", "depth of the impermeable layer below drain level.", "(m)", &Di);

  declparam("Dd", NHRU, "[0.5]", "0", "4", "depth of drains.", "(m)", &Dd);

  declparam("soil_poro_moist", NHRU, "[0.5]", "0", "1", "soil porosity of moist layer.", "()", &soil_poro_moist);

  declparam("L", NHRU, "[2]", "0", "20", "spacing between drains.", "(m)", &L);

  declparam("r", NHRU, "[0.1]", "0.01", "1", "drain radius.", "(m)", &r);


  variation_set = VARIATION_ORG;

  if(Global::nlay < 2){
    Global::nlay = 2;
    Global::maxlay = 2;
  }

  declvar("current_getstep", BASIN, "current getstep", "()", &current_getstep);

  declvar("redirected_residual", NHRU, "redirected residual after topping up Sd and soil_rechar in Netroute/_D/_M/_M_D.", "(mm*km^2/int)", &redirected_residual);

  declvar("redirected_residual_conc", NDEFN, "Concentration:: redirected residual after topping up Sd and soil_rechar in Netroute/_D/_M/_M_D.", "(mg/l)", &redirected_residual_conc, &redirected_residual_conc_lay, numsubstances);

  declstatdiag("cum_redirected_residual", NHRU, "cumulative HRU redirected_residual to another HRU.", "(mm*km^2)", &cum_redirected_residual);

  declstatdiag("cum_redirected_residual_mWQ", NDEFN, "mass of solute redirected_residual to another HRU.", "(mg*km^2)", &cum_redirected_residual_mWQ, &cum_redirected_residual_mWQ_lay, numsubstances);

  declstatvar("Sd", NHRU, "Depression storage.", "(mm)", &Sd);

  declstatvar("Sd_conc", NDEFN, "Concentration: Depression storage.", "(mg/l)", &Sd_conc, &Sd_conc_lay, numsubstances);

  declstatvar("gw", NHRU, "ground water storage.", "(mm)", &gw);

  declstatvar("gw_conc", NDEFN, "Concentration: ground water storage.", "(mg/l)", &gw_conc, &gw_conc_lay, numsubstances);

  declvar("solute_deposit", NHRU, "solute deposit left by evaporation.", "(mg)", &solute_deposit);

  declstatvar("cum_solute_deposit", NHRU, "cumulative solute deposit left by evaporation.", "(mg)", &cum_solute_deposit);

  declstatvar("soil_rechr", NHRU, "moisture content of soil recharge zone, ie, the"//
    "portion of the soil profile from which evaporation can take place.", "(mm)", &soil_rechr);

  declstatvar("soil_moist", NHRU, "moisture content of soil profile to the depth"//
    "of the rooting zone of the major vegetation type on the HRU.", "(mm)", &soil_moist);

  declstatvar("soil_moist_conc", NDEFN, "Concentration: moisture content of soil profile to the depth.", "(mg/l)", &soil_moist_conc, &soil_moist_conc_lay, numsubstances);

  declstatvar("potential", NHRU, ".", "(mm)", &potential);

  declstatvar("direct_excs", NHRU, ".", "(mm)", &direct_excs);

  declstatvar("potential_mWQ", NDEFN, ".", "(mg/l)", &potential_mWQ, &potential_mWQ_lay, numsubstances);

  declstatvar("direct_excs_mWQ", NDEFN, ".", "(mg/l)", &direct_excs_mWQ, &direct_excs_mWQ_lay, numsubstances);

  declstatvar("soil_moist_conc", NDEFN, "Concentration: moisture content of soil profile to the depth.", "(mg/l)", &soil_moist_conc, &soil_moist_conc_lay, numsubstances);

  declstatvar("soil_lower_max", NHRU, "maximum moisture content of lower soil profile to the depth"//
    "of the rooting zone of the major vegetation type on the HRU. (N.B. not Hype lower layer)", "(mm)", &soil_lower_max);

  declstatvar("soil_lower", NHRU, "moisture content of lower soil profile to the depth"//
    "of the rooting zone of the major vegetation type on the HRU. (N.B. not Hype lower layer)", "(mm)", &soil_lower);

  decllocal("cum_hru_condense", NHRU, "cumulative condensation over HRU.", "(mm)", &cum_hru_condense);

  declvar("cum_Sd_evap", NHRU, "cumulative Sd evaporation(included in hru_actet).", "(mm)", &cum_Sd_evap);

  declstatvar("cum_Sd_ssr", NHRU, "Accumulation of Sd excess from a HRU to ssr.", "(mm)", &cum_Sd_ssr);

  declstatvar("cum_Sd_gw", NHRU, "Accumulation of Sd excess from a HRU to gw.", "(mm)", &cum_Sd_gw);

  declstatvar("cum_lower_ssr", NHRU, "Accumulation of Sd excess from a HRU to ssr.", "(mm)", &cum_lower_ssr);

  declvar("soil_gw", NHRU, "Portion of excess soil water from a HRU that enters groundwater reservoirs.", "(mm/int)", &soil_gw);

  declvar("soil_gw_conc", NDEFN, "Concentration: Portion of excess soil water from a HRU that enters groundwater reservoirs.", "(mg/l)", &soil_gw_conc, &soil_gw_conc_lay, numsubstances);

  declvar("soil_gw_D", NHRU, "Portion of excess soil water from a HRU that enters groundwater reservoirs.", "(mm/d)", &soil_gw_D);

  declvar("gw_flow", NHRU, "Drainage from HRU ground water reservoir.", "(mm/int)", &gw_flow);

  declvar("gw_flow_conc", NDEFN, "Concentration: Drainage from HRU ground water reservoir.", "(mg/l)", &gw_flow_conc, &gw_flow_conc_lay, numsubstances);

  declvar("gw_flow_D", NHRU, "Daily drainage from HRU ground water reservoir.", "(mm/d)", &gw_flow_D);

  declvar("infil_act", NHRU, "Actual amount of water infiltrating the soil on each HRU.", "(mm/int)", &infil_act);

  declvar("infil_act_conc", NDEFN, "Concentration: Actual amount of water infiltrating the soil on each HRU.", "(mm/int)", &infil_act_conc, &infil_act_conc_lay, numsubstances);

  declvar("cum_infil_act", NHRU, "Accumulation of the actual amount of water infiltrating the soil on each HRU.", "(mm)", &cum_infil_act);

  declvar("cum_infil_act_mWQ", NDEFN, "mass of solute infiltrating the soil on each HRU.", "(mg)", &cum_infil_act_mWQ, &cum_infil_act_mWQ_lay, numsubstances);

  declvar("infil_act_D", NHRU, "Daily actual amount of water infiltrating the soil on each HRU.", "(mm/d)", &infil_act_D);

  declstatvar("cum_gw_flow", NHRU, "Accumulation of excess soil water from a HRU that enters groundwater reservoirs.", "(mm)", &cum_gw_flow);

  declstatvar("cum_gw_flow_mWQ", NDEFN, "mass of solute from a HRU that enters groundwater reservoirs.", "(mg)", &cum_gw_flow_mWQ, &cum_gw_flow_mWQ_lay, numsubstances);

  declvar("soil_ssr", NHRU, "Portion of soil moisture and recharge excess from a HRU that enters subsurface reservoirs.", "(mm/int)", &soil_ssr);

  declvar("soil_ssr_conc", NDEFN, "Concentration: Portion of soil moisture and recharge excess from a HRU that enters subsurface reservoirs.", "(mg/l)", &soil_ssr_conc, &soil_ssr_conc_lay, numsubstances);

  declvar("rechr_ssr", NHRU, "Portion of excess soil water from a HRU that enters subsurface reservoirs.", "(mm/int)", &rechr_ssr);

  declvar("rechr_ssr_conc", NDEFN, "Concentration: Portion of excess soil water from a HRU that enters subsurface reservoirs.", "(mg/l)", &rechr_ssr_conc, &rechr_ssr_conc_lay, numsubstances);

  declvar("rechr_ssr_conc", NDEFN, "Concentration: Portion of excess soil water from a HRU that enters subsurface reservoirs.", "(mg/l)", &rechr_ssr_conc, &rechr_ssr_conc_lay, numsubstances);

  declvar("funct_unmixed_layer", NHRU, "unmixed_layer function values.", "()", &funct_unmixed_layer);


  declstatvar("cum_soil_ssr", NHRU, "Accumulation of soil moisture from a HRU to ssr.", "(mm)", &cum_soil_ssr);

  declstatvar("cum_soil_ssr_mWQ", NDEFN, "mass of solute from a HRU to ssr.", "(mg)", &cum_soil_ssr_mWQ, &cum_soil_ssr_mWQ_lay, numsubstances);

  declstatvar("cum_rechr_ssr", NHRU, "Accumulation of Portion of excess from a HRU to ssr.", "(mm)", &cum_rechr_ssr);

  declstatvar("cum_rechr_ssr_mWQ", NDEFN, "mass of solute portion of excess from a HRU to ssr.", "(mg)", &cum_rechr_ssr_mWQ, &cum_rechr_ssr_mWQ_lay, numsubstances);

  declvar("soil_ssr_D", NHRU, "Portion of excess soil water from a HRU that enters subsurface reservoirs.", "(mm/d)", &soil_ssr_D);

  declvar("soil_runoff", NHRU, "Portion of excess soil water from a HRU to runoff.", "(mm/int)", &soil_runoff);

  declvar("soil_runoff_conc", NDEFN, "Concentration: Portion of excess soil water from a HRU to runoff.", "()", &soil_runoff_conc, &soil_runoff_conc_lay, numsubstances);

  declstatvar("cum_soil_runoff", NHRU, "Accumulation of Portion of excess soil water from a HRU to runoff.", "(mm)", &cum_soil_runoff);

  declstatvar("cum_soil_runoff_mWQ", NDEFN, "mass of solute of portion of excess soil water from a HRU to runoff.", "(mg)", &cum_soil_runoff_mWQ, &cum_soil_runoff_mWQ_lay, numsubstances);

  declvar("soil_runoff_D", NHRU, "Portion of excess soil water from a HRU that enters groundwater reservoirs.", "(mm/d)", &soil_runoff_D);

  declstatvar("cum_runoff_to_Sd", NHRU, "Cumulative portion of runoff to depression storage.", "(mm/int)", &cum_runoff_to_Sd);

  declstatvar("cum_runoff_to_Sd_mWQ", NDEFN, "mass of solute portion of runoff to depression storage.", "(mg)", &cum_runoff_to_Sd_mWQ, &cum_runoff_to_Sd_mWQ_lay, numsubstances);

  declstatvar("cum_runoff_to_ssr", NHRU, "Cumulative portion of runoff to interflow(ssr).", "(mm/int)", &cum_runoff_to_ssr);

  declstatvar("cum_soil_gw", NHRU, "Accumulation of excess soil water from a HRU that enters groundwater reservoirs.", "(mm)", &cum_soil_gw);

  declstatvar("cum_soil_gw_mWQ", NDEFN, "mass of solute of excess soil water from a HRU that enters groundwater reservoirs.", "(mg)", &cum_soil_gw_mWQ, &cum_soil_gw_mWQ_lay, numsubstances);


  declstatvar("soil_rechr_change_mWQ", NDEFN, "mass of solute soil_rechr change.", "(mg)", &soil_rechr_change_mWQ, &soil_rechr_change_mWQ_lay, numsubstances);

  declstatvar("soil_moist_change_mWQ", NDEFN, "mass of solute soil_moist change.", "(mg)", &soil_moist_change_mWQ, &soil_moist_change_mWQ_lay, numsubstances);

  declstatvar("soil_bottom_change_mWQ", NDEFN, "mass of solute soil_bottom change.", "(mg)", &soil_bottom_change_mWQ, &soil_bottom_change_mWQ_lay, numsubstances);

  declstatvar("Sd_change_mWQ", NDEFN, "mass of solute Sd change.", "(mg)", &Sd_change_mWQ, &Sd_change_mWQ_lay, numsubstances);

  declstatvar("gw_change_mWQ", NDEFN, "mass of solute gw change.", "(mg)", &gw_change_mWQ, &gw_change_mWQ_lay, numsubstances);


  decllocal("snowinfil_buf", NHRU, "buffer snow infiltration.", "(mm/d)", &snowinfil_buf);

  decllocal("runoff_buf", NHRU, "buffer runoff.", "(mm/d)", &runoff_buf);

  decllocal("meltrunoff_buf", NHRU, "buffer melt runoff.", "(mm/d)", &meltrunoff_buf);

  decllocal("hru_evap_buf", NHRU, "buffer evaporation.", "(mm/d)", &hru_evap_buf);

  decllocal("soil_rechr_Init", NHRU, "initial soil recharge.", "(mm)", &soil_rechr_Init);

  decllocal("soil_moist_Init", NHRU, "initial soil moisture.", "(mm)", &soil_moist_Init);

  decllocal("soil_bottom_Init", NHRU, "initial bottom soil moisture.", "(mm)", &soil_bottom_Init);

  decllocal("Sd_Init", NHRU, "initial Depression storage.", "(mm)", &Sd_Init);

  decllocal("gw_Init", NHRU, "initial ground water storage.", "(mm)", &gw_Init);

  declvar("soil_moist_conc_Init", NDEFN, "initial soil moisture conc.", "(mg/l)", &soil_moist_conc_Init, &soil_moist_conc_Init_lay, numsubstances);

  declvar("soil_bottom_conc_Init", NDEFN, "initial bottom soil moisture conc.", "(mg/l)", &soil_bottom_conc_Init, &soil_bottom_conc_Init_lay, numsubstances);

  declvar("soil_top_conc_Init", NDEFN, "initial top soil moisture conc.", "(mg/l)", &soil_top_conc_Init, &soil_top_conc_Init_lay, numsubstances);

  declvar("Sd_conc_Init", NDEFN, "initial concentration of nutrient species �layer� in the initial depression storage.", "(mg/l)", &Sd_conc_Init, &Sd_conc_Init_lay, numsubstances);

  declvar("gw_conc_Init", NDEFN, "initial concentration of nutrient species �layer�  in the groundwater reservoir.", "(mg/l)", &gw_conc_Init, &gw_conc_Init_lay, numsubstances);


  declparam("calcN", NHRU, "[0]", "0", "1", "flag for nitrogen simulation", "()", &calcN);

  declparam("calcP", NHRU, "[0]", "0", "1", "flag for phosphorus simulation", "()", &calcP);

  declparam("calcC", NHRU, "[0]", "0", "1", "flag for carbon simulation", "()", &calcC);

  declparam("unmixed_surface_layer", NHRU, "[10]", "0", "1000", "The depth of this mixed layer depends on tillage practices, but the top 5cm is mixed with most tillage events, and most of the top 15cm is at least partially mixed with most events. The 60cm is primarily for nitrates given that they are more likely to leach downward. In agricultural systems depths of 5cm, 15cm, and 60cm correspond with typical depths for soil testing.", "(mm)", &unmixed_surface_layer);

  declparam("CV_SWE", NHRU, "[1]", "0", "1", "Coefficient of variation: values can be taken from Gray, D. M., Toth, B., Zhao, L., Pomeroy, J. W., & Granger, R. J. (2001). Estimating areal snowmelt infiltration into frozen soils. Hydrological Processes, 15(16), 3095�3111. https://doi.org/10.1002/hyp.320", "()", &CV_SWE);

  declparam("basin_area", BASIN, "3", "1e-6", "1e+09", "total basin area.", "(km^2)", &basin_area);

  declparam("hru_area", NHRU, "[1]", "1e-6", "1e+09", "hru area.", "(km^2)", &hru_area);

  declparam("Sdmax", NHRU, "[0]", "0.0", "5000.0","Maximum depression storage.", "(mm)", &Sdmax);

  declparam("Sdinit", NHRU, "[0]", "0.0", "5000.0","Initial depression storage.", "(mm)", &Sdinit);

  declparam("Sd_conc_init", NDEFN, "[0]", "0.0", "10.0","Initial depression storage.", "(mg/l)", &Sd_conc_init, &Sd_conc_init_lay, numsubstances);

  declparam("soil_rechr_max", NHRU, "[60.0]", "0.0", "350.0",
    "Maximum value for soil recharge zone (upper portion of soil_moist where losses occur as both evaporation "//
    "and transpiration).  Must be less than or equal to soil_moist.","(mm)", &soil_rechr_max);

  declparam("soil_rechr_init", NHRU, "[30.0]", "0.0", "250.0", "Initial value for soil recharge zone (upper part of "//
    "soil_moist).  Must be less than or equal to soil_moist_init.", "(mm)", &soil_rechr_init);

  declparam("soil_moist_max", NHRU, "[375.0]", "0.0", "5000.0", "Maximum available water holding capacity of soil profile."//
    "Soil profile is surface to bottom of rooting zone.", "(mm)", &soil_moist_max);

  declparam("soil_moist_init", NHRU, "[187.0]", "0.0", "5000.0",
    "Initial value of available water in soil profile.", "(mm)", &soil_moist_init);

  declparam("soil_gw_K", NHRU, "[0.0]", "0.", "100.0", "The maximum amount of the soil water excess for an HRU "//
    "that is routed directly to the associated groundwater reservoir each day.", "(mm/d)", &soil_gw_K);

  declparam("gw_max", NHRU, "[375.0]", "0.0", "5000.0", "Maximum available water holding capacity of ground water reservoir.", "(mm)", &gw_max);

  declparam("gw_init", NHRU, "[187.0]", "0.0", "5000.0", "Initial value of available water in ground water reservoir.", "(mm)", &gw_init);

  declparam("gw_conc_init", NDEFN, "[1]", "0.0", "1.0", "Initial value of available water in ground water reservoir.", "(mg/l)", &gw_conc_init, &gw_conc_init_lay, numsubstances);

  declparam("gw_K", NHRU, "[0.0]", "0.", "100.0", "daily ground water drainage from gw reservoir.", "(mm/d)", &gw_K);

  declparam("rechr_ssr_K", NHRU, "[0.0]", "0.", "100.0", "daily ssr drainage from recharge.", "(mm/d)", &rechr_ssr_K);

  declparam("lower_ssr_K", NHRU, "[0.0]", "0.", "100.0", "daily ssr drainage from soil column.", "(mm/d)", &lower_ssr_K);

  declparam("Sd_ssr_K", NHRU, "[0.0]", "0.", "100.0", "daily depression storage ssr drainage factor.", "(mm/d)", &Sd_ssr_K);

  declparam("Sd_gw_K", NHRU, "[0.0]", "0.", "100.0", "daily depression storage gw drainage.", "(mm/d)", &Sd_gw_K);

  declparam("cov_type", NHRU,
    "[1]", "0", "2", "Vegetation evaporation type designation for HRU:  "//
    "0 = bare soil (no evaporation), 1 = crops (recharge layer), 2 = grasses & shrubs (all soil moisture).", "()", &cov_type);

  declparam("transp_limited", NHRU, "[0]", "0", "1", "limit transpiration to recharge layer only  on-1/off-0.", "()", &transp_limited);

  declparam("soil_ssr_runoff", NHRU, "[1]", "0", "1", "soil column excess to interflow(ssr)/runoff (and possibly Sd)  interflow-0/runoff-1.", "()", &soil_ssr_runoff);

  declparam("inhibit_evap", NHRU, "[0]", "0", "1", "inhibit evapatation, 1 -> inhibit.", "()", &inhibit_evap);

  declparam("rain_conc", NDEFN, "0", "0", "1000", "rain solute concentration", "(mg/l)", &rain_conc, &rain_conc_lay, numsubstances);

  declparam("Atmos_mWQ", NDEFN, "0", "0", "3", "atmospheric solute deposit", "(mg/int)", &Atmos_mWQ, &Atmos_mWQ_lay, numsubstances);

  declparam("soil_withdrawal", NDEFN, "[3]", "1", "4",
      "Select water withdrawal function for soil type: 1 = sand, 2 = loam, 3 = clay, 4 = organic. soil_withdrawal[1] - rechr layer, soil_withdrawal[2] - lower layer", "()",
      &soil_withdrawal, &soil_withdrawal_Tables, 2);


  declputvar("*", "hru_actet", "(mm/int)", &hru_actet);

  declputvar("*", "hru_cum_actet", "(mm)", &hru_cum_actet);


  evapDiv = declgetvar("*", "hru_evap", "(mm/int)", &hru_evap);


  declgetvar("*", "SWE", "(mm)", &SWE);

  declgetvar("*", "SWE_max", "(mm)", &SWE_max);

  declgetvar("*", "SWE_conc", "(mg/l)", &SWE_conc, &SWE_conc_lay);

  declputvar("*", "conc_top", "(mg/l)", &conc_top, &conc_top_lay);

  declputvar("*", "conc_bottom", "(mg/l)", &conc_bottom, &conc_bottom_lay);

  declputvar("*", "conc_below", "(mg/l)", &conc_below, &conc_below_lay);

  declgetvar("*", "infil", "(mm/int)", &infil);

  snowinfilDiv = declgetvar("*", "snowinfil", "(mm/int)", &snowinfil);

  runoffDiv = declgetvar("*", "runoff", "(mm/int)", &runoff);

  meltrunoffDiv = declgetvar("*", "meltrunoff", "(mm/int)", &meltrunoff);


  decllocal("redirected_residual_0", NHRU, "", "", &redirected_residual_0);

  decllocal("Sd_0", NHRU, "Depression storage.", "(mm)", &Sd_0);

  decllocal("gw_0", NHRU, "ground water storage.", "(mm)", &gw_0);

  decllocal("soil_rechr_0", NHRU, "moisture content of soil recharge zone.", "(mm)", &soil_rechr_0);

  decllocal("soil_moist_0", NHRU, "moisture content of soil profile to the depth.", "(mm)", &soil_moist_0);

  decllocal("soil_lower_0", NHRU, "moisture content of soil profile to the depth.", "(mm)", &soil_lower_0);

  decllocal("gw_flow_0", NHRU, "Drainage from HRU ground water reservoir.", "(mm/int)", &gw_flow_0);

  decllocal("hru_cum_actet_0", NHRU, "cumulative HRU redirected_residual to another HRU.", "(mm*km^2)", &hru_cum_actet_0);

  decllocal("cum_soil_runoff_0", NHRU, "cumulative HRU redirected_residual to another HRU.", "(mm*km^2)", &cum_soil_runoff_0);

  decllocal("cum_redirected_residual_0", NHRU, "cumulative HRU redirected_residual to another HRU.", "(mm*km^2)", &cum_redirected_residual_0);

  decllocal("cum_soil_ssr_0", NHRU, "Accumulation of soil moisture from a HRU to ssr.", "(mm)", &cum_soil_ssr_0);

  decllocal("cum_rechr_ssr_0", NHRU, "Accumulation of Portion of excess from a HRU to ssr.", "(mm)", &cum_rechr_ssr_0);

  decllocal("hru_actet_0", NDEFN, "","(mm/int)", &hru_actet_0);

  decllocal("cum_hru_condense_0", NHRU, "cumulative condensation over HRU.", "(mm)", &cum_hru_condense_0);

  decllocal("cum_Sd_evap_0", NHRU, "cumulative Sd evaporation(included in hru_actet).", "(mm)", &cum_Sd_evap_0);

  decllocal("cum_Sd_ssr_0", NHRU, "Accumulation of Sd excess from a HRU to ssr.", "(mm)", &cum_Sd_ssr_0);

  declstatvar("cum_Sd_gw_0", NHRU, "Accumulation of Sd excess from a HRU to gw.", "(mm)", &cum_Sd_gw_0);

  decllocal("cum_lower_ssr_0", NHRU, "Accumulation of Sd excess from a HRU to ssr.", "(mm)", &cum_lower_ssr_0);

  decllocal("cum_infil_act_0", NHRU, "Accumulation of the actual amount of water infiltrating the soil on each HRU.", "(mm)", &cum_infil_act_0);

  decllocal("cum_gw_flow_0", NHRU, "Accumulation of excess soil water from a HRU that enters groundwater reservoirs.", "(mm)", &cum_gw_flow_0);

  decllocal("cum_soil_runoff_0", NHRU, "Accumulation of Portion of excess soil water from a HRU to runoff.", "(mm)", &cum_soil_runoff_0);

  decllocal("cum_runoff_to_Sd_0", NHRU, "Cumulative portion of runoff to depression storage.", "(mm/int)", &cum_runoff_to_Sd_0);

  decllocal("cum_runoff_to_ssr_0", NHRU, "Cumulative portion of runoff to interflow(ssr).", "(mm/int)", &cum_runoff_to_ssr_0);

  decllocal("cum_soil_gw_0", NHRU, "Accumulation of excess soil water from a HRU that enters groundwater reservoirs.", "(mm)", &cum_soil_gw_0);

  decllocal("cum_solute_deposit_0", NHRU, "cumulative solute deposit left by evaporation.", "(mg)", &cum_solute_deposit_0);
}

void ClassWQ_Soil::init(void) {

  if(Global::nlay < numsubstances){
    Global::nlay = numsubstances;
    Global::maxlay = numsubstances;
  }

  FaultsAllowed = 0;
  nhru = getdim(NHRU);
  nlay = getdim(NLAY);

  if(snowinfilDiv > 1){
    String S = "Soil:  \"snowinfil\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(evapDiv > 1){
    String S = "Soil:  \"hru_evap\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(meltrunoffDiv > 1){
    String S = "Netroute:  \"meltrunoff\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(runoffDiv > 1){
    String S = "Netroute:  \"runoff\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  for(Sub = 0; Sub < numsubstances; Sub++){
    for(hh = 0; hh < nhru; ++hh){
      if(Sub == 0){
        soil_rechr[hh] = soil_rechr_init[hh];
        soil_moist[hh] = soil_moist_init[hh];

        soil_lower_max[hh] = soil_moist_max[hh] - soil_rechr_max[hh];
        soil_lower[hh] = soil_moist[hh] - soil_rechr[hh];

        soil_moist_conc_lay[Sub][hh] = 0.0;

        solute_deposit[hh] = 0.0;
        cum_solute_deposit[hh] = 0.0;

        soil_runoff_D[hh] = 0.0;
        soil_ssr_D[hh] = 0.0;
        soil_gw_D[hh] = 0.0;
        gw_flow_D[hh] = 0.0;
        infil_act_D[hh] = 0.0;

        hru_cum_actet[hh] = 0.0;
        cum_hru_condense[hh] = 0.0;
        cum_Sd_evap[hh] = 0.0;

        cum_Sd_ssr[hh] = 0.0;
        cum_Sd_gw[hh] = 0.0;
        cum_lower_ssr[hh] = 0.0;
        cum_runoff_to_ssr[hh] = 0.0;

        if(soil_rechr[hh] > soil_moist[hh]) {
          soil_rechr[hh] = soil_moist[hh];
          string S = string("'") + Name + " (Soil)' Soil_rechr greater than soil_moist, soil_rechr set to soil_moist, hru = " + String(hh).c_str();
          CRHMException TExcept(S.c_str(), WARNING);
          LogError(TExcept);
          throw TExcept;
        }

        if(soil_rechr_max[hh] > soil_moist_max[hh]) {
          string S = string("'") + Name + " (Soil)' Soil_rechr_max cannot be greater than soil_moist_max in hru = " + String(hh+1).c_str();
          CRHMException TExcept(S.c_str(), TERMINATE);
          LogError(TExcept);
          throw TExcept;
        }

        if(Sdinit[hh] > Sdmax[hh]) {
          string S = string("'") + Name + " (Soil)' Initial value of depression storage is greater than the maximum value in hru = " + String(hh+1).c_str();
          CRHMException Except(S.c_str() ,TERMINATE);
          LogError(Except);
          throw Except;
        }

        if(variation == VARIATION_1){
          if(culvert_water_Dmax[hh]/culvert_diam[hh] > 2.5){
            String S = "soil: " + String(Name.c_str()) +  " ratio of H/D > 2.5 in HRU " + String(hh+1);
            CRHMException TExcept(S.c_str(), WARNING);
            LogError(TExcept);
          }
          culvert_water_V[hh] = 0.0;
          culvert_water_H[hh] = 0.0;
          culvert_water_A[hh] = 0.0;
          culvert_over_Q[hh] = 0.0;
          culvert_Q[hh] = 0.0;
          culvert_evap[hh] = 0.0;
          cum_culvert[hh] = 0.0;
          cum_culvert_over[hh] = 0.0;
        }

        if(variation == VARIATION_2){
          tile_flow[hh] = 0.0;
          x[hh] = 2.0*M_PI*(Di[hh] - Dd[hh])/L[hh];
          Dequiv_depth[hh] =  M_PI*L[hh]/8.0*(M_PI*log(L[hh]/(r[hh])) + FunctX(x[hh]));
          Dw[hh] = soil_moist_init[hh]/soil_poro_moist[hh];
        }
      }

      Reset_WQ(hh, infil_act, infil_act_conc_lay[Sub]);

      Reset_WQ(hh, redirected_residual, redirected_residual_conc_lay[Sub]);
      Reset_WQ(hh, cum_redirected_residual, cum_redirected_residual_mWQ_lay[Sub]);
      Reset_WQ(hh, cum_infil_act, cum_infil_act_mWQ_lay[Sub]);


      Reset_WQ(hh, cum_soil_runoff, cum_soil_runoff_mWQ_lay[Sub]);
      Reset_WQ(hh, cum_soil_ssr, cum_soil_ssr_mWQ_lay[Sub]);
      Reset_WQ(hh, cum_rechr_ssr, cum_rechr_ssr_mWQ_lay[Sub]);
      Reset_WQ(hh, cum_soil_gw, cum_soil_gw_mWQ_lay[Sub]);
      Reset_WQ(hh, cum_gw_flow, cum_gw_flow_mWQ_lay[Sub]);
      Reset_WQ(hh, cum_runoff_to_Sd, cum_runoff_to_Sd_mWQ_lay[Sub]);

      Set_WQ(hh, Sd, Sd_conc_lay[Sub], Sdinit[hh], Sd_conc_init_lay[Sub][hh]);
      Set_WQ(hh, gw, gw_conc_lay[Sub], gw_init[hh], gw_conc_init_lay[Sub][hh]);
   } // hh
  } // Sub
}

void ClassWQ_Soil::run(void) {

try{
  float excs, excs_mWQ, condense, et;

  long step = getstep();
  current_getstep[0] = step;
  long nstep = step%Global::Freq;

  for(Sub = 0; Sub < numsubstances; Sub++){

    if(Sub == 0) // saves all HRUs, non WQ variables the first time
      Save();

    if(step == 1){ // First Time only
      String S = String("Initial Substance# ") + (Sub+1);
      LogDebug(S.c_str());
      LogMessage(" ");

      for(hh = 0; chkStruct(); ++hh) {
        soil_top_conc_Init_lay[Sub][hh] = conc_top_lay[Sub][hh];
        soil_bottom_conc_Init_lay[Sub][hh] = conc_bottom_lay[Sub][hh];
        if(soil_rechr[hh] + soil_lower[hh] > 0.0)
          soil_moist_conc_Init_lay[Sub][hh] = (soil_rechr[hh]*conc_top_lay[Sub][hh] + soil_lower[hh]*conc_bottom_lay[Sub][hh])/(soil_rechr[hh] + soil_lower[hh]); // mg/l
        else
          soil_moist_conc_Init_lay[Sub][hh] = 0.0; // mg/l

        Sd_Init[hh] = Sd[hh];
        Sd_conc_Init_lay[Sub][hh] = Sd_conc_lay[Sub][hh];

        gw_Init[hh] = gw[hh];
        gw_conc_Init_lay[Sub][hh] = gw_conc_lay[Sub][hh];

// soil_rechr_conc, soil_lower_conc and soil_moist_conc use the Hype initial values

        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr_init     (mm) (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr[hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr_init_conc     (mm) (mm*hru) (mm*hru/basin): ").c_str(), conc_top_lay[Sub][hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist_init     (mm) (mm*hru) (mm*hru/basin): ").c_str(), soil_moist[hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist_init_conc     (mm) (mm*hru) (mm*hru/basin): ").c_str(), soil_moist_conc_Init_lay[Sub][hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_lower_init     (mm) (mm*hru) (mm*hru/basin): ").c_str(), soil_lower[hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_lower_init_conc     (mm) (mm*hru) (mm*hru/basin): ").c_str(), conc_bottom_lay[Sub][hh], hru_area[hh], basin_area[0]);

        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sdinit              (mm) (mm*hru) (mm*hru/basin): ").c_str(), Sd[hh],         hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd_init_conc             (mm) (mm*hru) (mm*hru/basin): ").c_str(), Sd_conc_lay[Sub][hh],         hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw_init             (mm) (mm*hru) (mm*hru/basin): ").c_str(), gw[hh],         hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw_init_conc             (mm) (mm*hru) (mm*hru/basin): ").c_str(), gw_conc_lay[Sub][hh],         hru_area[hh], basin_area[0]);
        LogDebug(" ");
      } // for
    } // if step == 1

    for(hh = 0; chkStruct(); ++hh){

      if(Sub != 0)
        Restore(hh);

      if(snowinfilDiv == 1) // interval value
         snowinfil_buf[hh] = snowinfil[hh];

      if(runoffDiv == 1) // interval value
         runoff_buf[hh] = runoff[hh];

      if(meltrunoffDiv == 1) // interval value
         meltrunoff_buf[hh] = meltrunoff[hh];

      if(evapDiv == 1) // interval value
         hru_evap_buf[hh] = hru_evap[hh];

//      if(nstep == 1){ // beginning of every day

// update from Hype

        if(soil_rechr[hh] + soil_lower[hh] > 0.0)
          soil_moist_conc_lay[Sub][hh] = (soil_rechr[hh]*conc_top_lay[Sub][hh] + soil_lower[hh]*conc_bottom_lay[Sub][hh])/(soil_rechr[hh] + soil_lower[hh]); // mg/l
        else
          soil_moist_conc_lay[Sub][hh] = 0.0; // mg/l

        soil_runoff_D[hh] = 0.0;
        soil_ssr_D[hh] = 0.0;
        soil_gw_D[hh] = 0.0;
        gw_flow_D[hh] = 0.0;
        infil_act_D[hh] = 0.0;
        gw_conc_lay[Sub][hh] = 0.0;
//      } // if

      solute_deposit[hh] = 0.0;
      hru_actet[hh] = 0.0;

      Reset_WQ(hh, soil_gw, soil_gw_conc_lay[Sub]);
      Reset_WQ(hh, soil_ssr, soil_ssr_conc_lay[Sub]);
      Reset_WQ(hh, rechr_ssr, rechr_ssr_conc_lay[Sub]);
      Reset_WQ(hh, infil_act, infil_act_conc_lay[Sub]);
      Reset_WQ(hh, soil_runoff, soil_runoff_conc_lay[Sub]);

      if(hru_evap_buf[hh] < 0.0) {
        condense = -hru_evap_buf[hh];
        cum_hru_condense[hh] += condense;
        hru_evap_buf[hh] = 0.0;
      }
      else
        condense = 0.0;

  //******Add infiltration to soil and compute excess

      if(soil_moist_max[hh] > 0.0){

        potential[hh] = infil[hh] + snowinfil_buf[hh] + condense;

        if(potential[hh] > 0.0){
          potential_mWQ_lay[Sub][hh] = (infil[hh] + condense)*rain_conc_lay[Sub][hh] + snowinfil_buf[hh]*SWE_conc_lay[Sub][hh];
          if(!inhibit_evap[hh]) // only when no snowcover
            potential_mWQ_lay[Sub][hh] += Atmos_mWQ_lay[Sub][hh];
        }
        else
          potential_mWQ_lay[Sub][hh] = 0.0;

        if(soil_moist[hh] + potential[hh] > soil_moist_max[hh]){
          direct_excs[hh] = soil_moist[hh] + potential[hh] - soil_moist_max[hh];
          direct_excs_mWQ_lay[Sub][hh] = potential_mWQ_lay[Sub][hh]*direct_excs[hh]/potential[hh];
          potential[hh] -= direct_excs[hh];
          potential_mWQ_lay[Sub][hh] -=  direct_excs_mWQ_lay[Sub][hh];
        }
        else{
          direct_excs[hh] = 0.0;
          direct_excs_mWQ_lay[Sub][hh] = 0.0;
        }

        soil_lower[hh] = soil_moist[hh] - soil_rechr[hh];

        if(soil_lower[hh])
          conc_bottom_lay[Sub][hh] = (soil_moist_conc_lay[Sub][hh]*soil_moist[hh] - conc_top_lay[Sub][hh]*soil_rechr[hh])/soil_lower[hh];
        else
          conc_bottom_lay[Sub][hh] = 0.0;

        if(soil_rechr[hh] + potential[hh] > soil_rechr_max[hh]){
          excs = soil_rechr[hh] + potential[hh] - soil_rechr_max[hh];
          conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + potential_mWQ_lay[Sub][hh];
          conc_top_lay[Sub][hh] /= (soil_rechr[hh] + excs);
          excs_mWQ = conc_top_lay[Sub][hh]*excs;
          soil_rechr[hh] = soil_rechr_max[hh];
        }
        else{
          excs = 0.0;
          excs_mWQ = 0.0;
          soil_rechr[hh] = soil_rechr[hh] + potential[hh];

          conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + potential_mWQ_lay[Sub][hh];
          if(soil_rechr[hh] > 0.0)
            conc_top_lay[Sub][hh] /= soil_rechr[hh];
          else
            conc_top_lay[Sub][hh] = 0.0;

          soil_moist[hh] = soil_lower[hh] + soil_rechr[hh];
          soil_moist_conc_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + conc_top_lay[Sub][hh]*soil_rechr[hh])/soil_moist[hh]; // amount used
        }

        if(excs > 0.0){
          conc_bottom_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + excs_mWQ)/soil_lower[hh]; // put remaning excs_mWQ in lower (some already in rechr)
          soil_moist[hh] = soil_lower[hh] + soil_rechr[hh] + excs;
          soil_moist_conc_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + conc_top_lay[Sub][hh]*soil_rechr[hh] + excs_mWQ)/soil_moist[hh];

          if(soil_moist[hh] > soil_moist_max[hh]){
            float excs0 = soil_moist[hh] - soil_moist_max[hh];
            soil_moist_conc_lay[Sub][hh] = (soil_moist_conc_lay[Sub][hh]*soil_moist[hh] - excs_mWQ*excs0/excs); // amount used

            conc_bottom_lay[Sub][hh] = (soil_moist_conc_lay[Sub][hh] - conc_top_lay[Sub][hh]*soil_rechr[hh]); // amount used

            soil_moist[hh] = soil_moist_max[hh];
            soil_moist_conc_lay[Sub][hh] /= soil_moist_max[hh]; // amount used

            soil_lower[hh] = soil_moist[hh] - soil_rechr[hh];
            conc_bottom_lay[Sub][hh] /= soil_lower[hh]; // amount used

            excs_mWQ = excs_mWQ*excs0/excs;
            excs = excs0;
          }
          else{
            soil_lower[hh] = soil_moist[hh] - soil_rechr[hh];
            conc_bottom_lay[Sub][hh] = (soil_moist_conc_lay[Sub][hh]*soil_moist[hh] - conc_top_lay[Sub][hh]*soil_rechr[hh])/soil_lower[hh]; // amount used

            excs = 0.0;
            excs_mWQ = 0.0;
          }
        }

        infil_act_conc_lay[Sub][hh] = potential_mWQ_lay[Sub][hh];
        cum_infil_act_mWQ_lay[Sub][hh] += infil_act_conc_lay[Sub][hh];

        infil_act[hh] = potential[hh];
        cum_infil_act[hh] += infil_act[hh];
        infil_act_D[hh] += infil_act[hh];

  //  Handle subsurface runoff - does not affect WQ - concentrates

        if(!inhibit_evap[hh]){
          if(variation == VARIATION_2){ // handle tile drainage
            Dw[hh] = soil_moist_init[hh]/soil_poro_moist[hh]/1000;  // convert mm -> m
            tile_flow[hh] = 8.0*Kb[hh]*Dequiv_depth[hh]*(Dd[hh] - Dw[hh]) + 4.0*Ka[hh]*sqr(Dd[hh] - Dw[hh])/sqr(L[hh]);
            tile_flow[hh] /= Global::Freq;

            if(tile_flow[hh] > 0.0){
              float Actual = soil_rechr[hh]/soil_rechr_max[hh]*tile_flow[hh]; // take out of recharge layer

              if(Actual > soil_rechr[hh])
                rechr_ssr[hh] = soil_rechr[hh];
              else
                rechr_ssr[hh] = Actual;

              soil_rechr[hh] -= rechr_ssr[hh];
              soil_moist[hh] -= rechr_ssr[hh];
              rechr_ssr_conc_lay[Sub][hh] = conc_top_lay[Sub][hh];
              soil_ssr[hh] = rechr_ssr[hh];

             float remainder = Actual - rechr_ssr[hh];
             if(remainder > soil_moist[hh])
               remainder = soil_moist[hh];

             soil_ssr[hh] += remainder;
             soil_moist[hh] -= remainder;
             
             soil_lower[hh] = soil_moist[hh] - soil_rechr[hh];
             rechr_ssr_conc_lay[Sub][hh] = (rechr_ssr_conc_lay[Sub][hh]*rechr_ssr[hh] + conc_bottom_lay[Sub][hh]*remainder)/(rechr_ssr[hh] + remainder);
             rechr_ssr[hh] += remainder;
            }
          } // handle tile drainage
          else{
            rechr_ssr[hh] = soil_rechr[hh]/soil_rechr_max[hh]*rechr_ssr_K[hh]/Global::Freq;
            if(rechr_ssr[hh] > 0.0){
              soil_ssr[hh] = rechr_ssr[hh];
              soil_ssr_conc_lay[Sub][hh] = conc_top_lay[Sub][hh];

              soil_rechr[hh] -= rechr_ssr[hh];
              rechr_ssr_conc_lay[Sub][hh] = conc_top_lay[Sub][hh];

              if(soil_rechr[hh] < 0.0){
                soil_moist[hh] -= soil_rechr[hh];
                soil_rechr[hh] = 0.0;
                conc_top_lay[Sub][hh] = 0.0;
              }
              else
                soil_moist[hh] -= rechr_ssr[hh];

              soil_lower[hh] = soil_moist[hh] - soil_rechr[hh];
              soil_moist_conc_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + soil_rechr[hh]*conc_top_lay[Sub][hh])/soil_moist[hh];
            } // handle non variation ssr
          }

          cum_rechr_ssr[hh] += rechr_ssr[hh];
          cum_rechr_ssr_mWQ_lay[Sub][hh] += rechr_ssr_conc_lay[Sub][hh]*rechr_ssr[hh];
        } // !inhibit_evap

  //  Handle excess to gw

        float s2gw_k = soil_gw_K[hh]/Global::Freq;
        if(s2gw_k > 0)
          if(direct_excs[hh] >= s2gw_k) { // to gw 03/04/10 changed from >
            soil_gw[hh] = s2gw_k; // soil_gw[hh] = s2gw_k;
            soil_gw_conc_lay[Sub][hh] = direct_excs_mWQ_lay[Sub][hh];
            direct_excs_mWQ_lay[Sub][hh] -= soil_gw_conc_lay[Sub][hh]*s2gw_k;
            direct_excs[hh] -= s2gw_k;
          }
          else { // to gw
            soil_gw[hh] = direct_excs[hh];
            soil_gw_conc_lay[Sub][hh] = direct_excs_mWQ_lay[Sub][hh];
            direct_excs[hh] = 0.0;
            direct_excs_mWQ_lay[Sub][hh] = 0.0;
          }

  //  Handle excess to interflow or runoff

        if(!soil_ssr_runoff[hh] && direct_excs[hh] > 0.0){ // to interflow
          soil_ssr_conc_lay[Sub][hh] = soil_ssr_conc_lay[Sub][hh]*soil_ssr[hh] + direct_excs_mWQ_lay[Sub][hh];
          soil_ssr[hh] += direct_excs[hh];
          if(soil_ssr[hh] > 0.0)
            soil_ssr_conc_lay[Sub][hh] /= soil_ssr[hh];
          else
            soil_ssr_conc_lay[Sub][hh] = 0.0;
          direct_excs[hh] = 0.0;
          direct_excs_mWQ_lay[Sub][hh] = 0.0;
        }
      }
      else{ // soil_moist_max <= 0.0, i.e. Pond
        excs = infil[hh] + snowinfil_buf[hh] + condense;
          excs_mWQ = (infil[hh] + condense)*rain_conc_lay[Sub][hh] + snowinfil_buf[hh]*SWE_conc_lay[Sub][hh];
          if(!inhibit_evap[hh]) // only when no snowcover
            excs_mWQ += Atmos_mWQ_lay[Sub][hh];
      }

      float runoff_to_Sd = 0.0;

      soil_runoff[hh] = direct_excs[hh] + meltrunoff_buf[hh] + runoff_buf[hh] + redirected_residual[hh]/hru_area[hh]; // last term (mm*km^2/int)

      if(soil_runoff[hh] > 0.0){
        if(SWE_max[hh] > 0.0)
          funct_unmixed_layer[hh] = tanh(1.26*SWE[hh]/(CV_SWE[hh]*SWE_max[hh]));
        else
          funct_unmixed_layer[hh] = 1.0;

        float amount = conc_top_lay[Sub][hh]*unmixed_surface_layer[hh]*(1.0 - funct_unmixed_layer[hh]);

        if(amount > conc_top_lay[Sub][hh]*soil_rechr[hh])
          conc_top_lay[Sub][hh] = 0.0;
        else
          conc_top_lay[Sub][hh] = (conc_top_lay[Sub][hh]*soil_rechr[hh] - amount)/soil_rechr[hh];

        soil_runoff_conc_lay[Sub][hh] = (excs_mWQ + direct_excs_mWQ_lay[Sub][hh] + meltrunoff_buf[hh]*SWE_conc_lay[Sub][hh] +
          runoff_buf[hh]*rain_conc_lay[Sub][hh] + redirected_residual_conc_lay[Sub][hh]*redirected_residual[hh]/hru_area[hh]
          + amount)/soil_runoff[hh]; // last term (mm*km^2/int)
      }
      else
        soil_runoff_conc_lay[Sub][hh] = 0.0;

      cum_redirected_residual[hh] += redirected_residual[hh];

      redirected_residual[hh] = 0;

      if(soil_runoff[hh] > 0.0 && Sdmax[hh] > 0.0){
        float Fix = -12.0;
        if(soil_runoff[hh]/Sdmax[hh] < 12.0)
          Fix = -soil_runoff[hh]/Sdmax[hh];
        float Ds = (Sdmax[hh] - Sd[hh])*(1 - exp(Fix));

        if(soil_moist_max[hh] <= 0.0) // handle pond
          Ds = Sdmax[hh] - Sd[hh];

        if(Ds > 0.0){
          if(soil_runoff[hh] > Ds){
            soil_runoff[hh] -= Ds;
            if(soil_runoff[hh] < 0.0){
              soil_runoff[hh] = 0.0;
              soil_runoff_conc_lay[Sub][hh] = 0.0;
            }
            Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + soil_runoff_conc_lay[Sub][hh]*Ds;
            Sd[hh] += Ds;
            Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]/Sd[hh];
            runoff_to_Sd += Ds;
          }
          else{
            Sd[hh] += soil_runoff[hh];
            Sd_conc_lay[Sub][hh] += soil_runoff_conc_lay[Sub][hh];
            runoff_to_Sd += soil_runoff[hh];
            soil_runoff[hh] = 0.0;
            soil_runoff_conc_lay[Sub][hh] = 0.0;
          }
        }
      }

      if(variation == VARIATION_1){

        float culvert_C[5] = {0.5, 0.6, 0.7, 0.75, 0.97};

        culvert_water_H[hh] = 0.0;
        culvert_water_A[hh] = 0.0;
        culvert_over_Q[hh] = 0.0;
        culvert_Q[hh] = 0.0;
        culvert_evap[hh] = 0.0;

        if((soil_runoff[hh] > 0.0 || culvert_water_V[hh] > 0.0) && number_culverts[hh] > 0.0){ // culvert addition. Inputs are in (mm)
          culvert_water_V[hh] += soil_runoff[hh]*(hru_area[hh]*1000.0); // convert mm to m3
          soil_runoff[hh] = 0.0;

          culvert_water_H[hh] = pow(3.0*culvert_water_V[hh]*channel_slope[hh]*side_slope[hh], 1.0/3.0); // (m)

          if(culvert_water_H[hh] > 0.0){

  // calculate overflow if culvert level exceeds Dmax

            culvert_water_H[hh] = pow(3.0*culvert_water_V[hh]*channel_slope[hh]*side_slope[hh], 1.0/3.0); // (m)

            if(culvert_water_H[hh] > culvert_water_Dmax[hh]){ // (m) overflow over road
              culvert_water_H[hh] = culvert_water_Dmax[hh]; // (m)
              float maxVol = pow(culvert_water_Dmax[hh], 3.0)/(3.0*channel_slope[hh]*side_slope[hh]); // (m3)

              culvert_over_Q[hh] = (culvert_water_V[hh] - maxVol)/86400.0*Global::Freq; // (m3) to (m3/int)
              culvert_water_V[hh] = maxVol; // (m3)

              cum_culvert_over[hh] += culvert_over_Q[hh]*86400.0/Global::Freq; // (m3/s) to (m3)
              soil_runoff[hh] += culvert_over_Q[hh]*86400.0/Global::Freq/(hru_area[hh]*1000.0); // (m3/s) to (mm/int)
            }
            HD[hh] = culvert_water_H[hh]/culvert_diam[hh];

  // calculate flow through culvert

            if(HD[hh] <= 0.0)
              culvert_Q[hh] = 0.0;
            else if(HD[hh] < 1.5)
              culvert_Q[hh] = max <float>((-0.544443*pow(HD[hh], 4.0) + 0.221892*pow(HD[hh], 3.0) + 2.29756*pow(HD[hh], 2.0)
                   + 0.159413*HD[hh] + 0.00772254)*culvert_C[culvert_type[hh]]*number_culverts[hh]*pow(culvert_diam[hh], 2.5), 0.0); // (m3/s)
            else
              culvert_Q[hh] = culvert_C[culvert_type[hh]]*number_culverts[hh]*0.785*pow(culvert_diam[hh], 2.5)*sqrt(2.0*9.81*(HD[hh] - 0.5));

            if(culvert_water_V[hh] > culvert_Q[hh]*86400.0/Global::Freq) // (m3/s) to (m3))
              culvert_water_V[hh] -= culvert_Q[hh]*86400.0/Global::Freq; // (m3/s) to (m3)
            else{
              culvert_Q[hh] = culvert_water_V[hh]*Global::Freq/86400.0;  // (m3) to (m3/s)
              culvert_water_V[hh] = 0.0;
            }

            cum_culvert[hh] += culvert_Q[hh]*86400.0/Global::Freq; // (m3/s) to (m3/int)
            soil_runoff[hh] += culvert_Q[hh]*86400.0/Global::Freq/(hru_area[hh]*1000.0); // (m3/s) to (mm/int)
          }
          culvert_water_A[hh] = sqr(culvert_water_H[hh])/(channel_slope[hh]*side_slope[hh]); // used for evaporation
        } // culvert addition
      }

      soil_runoff_D[hh] += soil_runoff[hh];
      cum_soil_runoff[hh] += soil_runoff[hh];
      cum_soil_runoff_mWQ_lay[Sub][hh] += soil_runoff[hh]*soil_runoff_conc_lay[Sub][hh];
      cum_runoff_to_Sd[hh] += runoff_to_Sd;
      cum_runoff_to_Sd_mWQ_lay[Sub][hh] += runoff_to_Sd*soil_runoff_conc_lay[Sub][hh];

      if(Sd[hh] > 0.0 && Sd_gw_K[hh] > 0.0){
        float Sd2gw_k = Sd_gw_K[hh]/Global::Freq;
        if(Sd2gw_k > Sd[hh])
          Sd2gw_k = Sd[hh];

        soil_gw_conc_lay[Sub][hh] = soil_gw_conc_lay[Sub][hh]*soil_gw[hh] + Sd_conc_lay[Sub][hh]*Sd2gw_k;
        soil_gw[hh] += Sd2gw_k;

        if(soil_gw[hh] > 0.0)
          soil_gw_conc_lay[Sub][hh] /= soil_gw[hh];
        else
          soil_gw_conc_lay[Sub][hh] = 0.0;

        Sd[hh] -= Sd2gw_k;

        if(Sd[hh] < 0.0){
          Sd[hh] = 0.0;
          Sd_conc_lay[Sub][hh] = 0.0;
        }
        cum_Sd_gw[hh] += Sd2gw_k;
      }

      soil_gw_D[hh] += soil_gw[hh];
      cum_soil_gw[hh] += soil_gw[hh];
      cum_soil_gw_mWQ_lay[Sub][hh] += soil_gw[hh]*Sd_conc_lay[Sub][hh];

      gw[hh] += soil_gw[hh];
      gw_flow[hh] = 0.0;
      gw_flow_conc_lay[Sub][hh] = 0.0;

      if(gw[hh] > gw_max[hh]){
        gw_flow_conc_lay[Sub][hh] = gw_conc_lay[Sub][hh];
        gw_flow[hh] += gw[hh] - gw_max[hh];
        gw[hh] = gw_max[hh];
      }

      if(gw_max[hh] > 0.0){ // prevents divide by zero error
        float spill = gw[hh]/gw_max[hh]*gw_K[hh]/Global::Freq;
        gw[hh] -= spill;
        gw_flow_conc_lay[Sub][hh] = gw_flow_conc_lay[Sub][hh]*gw_flow[hh] - gw_conc_lay[Sub][hh]*spill;
        gw_flow[hh] += spill;
        if(gw_flow[hh] > 0.0)
          gw_flow_conc_lay[Sub][hh] /= gw_flow[hh];
        else
          gw_flow_conc_lay[Sub][hh] = 0.0;
      }

      gw_flow_D[hh] += gw_flow[hh];
      cum_gw_flow[hh] += gw_flow[hh];
      cum_gw_flow_mWQ_lay[Sub][hh] += gw_flow[hh]*gw_conc_lay[Sub][hh];

      if(Sd[hh] > 0.0 && Sd_ssr_K[hh] > 0.0){
        float Sd2ssr_k = Sd_ssr_K[hh]/Global::Freq; // WHY not proportional?
        if(Sd2ssr_k >= Sd[hh])
          Sd2ssr_k = Sd[hh];

        soil_ssr_conc_lay[Sub][hh] = soil_ssr_conc_lay[Sub][hh]*soil_ssr[hh] + Sd_conc_lay[Sub][hh]*Sd2ssr_k;
        soil_ssr[hh] += Sd2ssr_k;
        if(soil_ssr[hh] > 0.0)
          soil_ssr_conc_lay[Sub][hh] /= soil_ssr[hh];
        else
          soil_ssr_conc_lay[Sub][hh] = 0.0;

        if(Sd[hh] - Sd2ssr_k < 0.0){
          Sd[hh] = 0.0;
          Sd_conc_lay[Sub][hh] = 0.0;
        }
        else{
          Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] -  Sd_conc_lay[Sub][hh]*Sd2ssr_k;
          Sd[hh] -= Sd2ssr_k;
          if(Sd[hh] > 0.0)
            Sd_conc_lay[Sub][hh] /= Sd[hh];
          else
            Sd_conc_lay[Sub][hh] = 0.0;
        }
      }

      float s2ssr_k = lower_ssr_K[hh]/Global::Freq;
      if(s2ssr_k > 0.00001){
        float avail = soil_lower[hh];
        if(s2ssr_k >= avail)
          s2ssr_k = avail;

        soil_lower[hh] -= s2ssr_k; // soil lower concentration stays the same

        soil_ssr_conc_lay[Sub][hh] = soil_ssr_conc_lay[Sub][hh]*soil_ssr[hh] + conc_bottom_lay[Sub][hh]*s2ssr_k;
        soil_ssr[hh] += s2ssr_k;

        if(soil_ssr[hh] > 0.0)
          soil_ssr_conc_lay[Sub][hh] /= soil_ssr[hh];
        else
          soil_ssr_conc_lay[Sub][hh] = 0.0;

        soil_moist[hh] = soil_lower[hh] + soil_rechr[hh];
        soil_moist_conc_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + soil_rechr[hh]*conc_top_lay[Sub][hh])/soil_moist[hh];

        cum_lower_ssr[hh] += s2ssr_k;
      }

      cum_soil_ssr[hh] += soil_ssr[hh];
      cum_soil_ssr_mWQ_lay[Sub][hh] += soil_ssr[hh]*soil_moist_conc_lay[Sub][hh];
      soil_ssr_D[hh] += soil_ssr[hh];

  //******Compute actual evapotranspiration

      float culvert_pond = 0.0; // convert m3 to mm

      float culvert_evapL = 0;

      if(variation == VARIATION_1 && culvert_water_V[hh] > 0.0 && hru_evap_buf[hh] > 0.0){ // conditions for culvert evaporation
        culvert_pond = culvert_water_V[hh]/(hru_area[hh]*1000.0); // convert m3 to mm over HRU area
        culvert_evapL = hru_evap_buf[hh]*culvert_water_A[hh]/(hru_area[hh]*1e6); // calculate potential evaporation from pond

        if(culvert_evapL > culvert_pond) // limit to amount available
          culvert_evapL = culvert_pond;

        culvert_evap[hh] = culvert_evapL;
        hru_actet[hh] += culvert_evapL;
        culvert_water_V[hh] = (culvert_pond - culvert_evapL)*(hru_area[hh]*1000.0); // remove evaporation from culvert pond and convert to volume
      }

      float avail_evap = hru_evap_buf[hh] - culvert_evapL;
      if(Sd[hh] + soil_moist[hh] + culvert_pond > 0.0)
        avail_evap *= (Sd[hh]/(Sd[hh] + soil_moist[hh]));
      else
        avail_evap = 0.0;

      if(Sd[hh] > 0.0 && avail_evap > 0.0){
        if(Sd[hh] >= avail_evap){
          if(Sd[hh] - avail_evap < 0.0){
            Sd[hh] -= avail_evap;
            Sd[hh] = 0.0;
            Sd_conc_lay[Sub][hh] = 0.0;
          }
          else{
            Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] -  Sd_conc_lay[Sub][hh]*avail_evap;
            Sd[hh] -= avail_evap;

            if(Sd[hh] > 0.0)
              Sd_conc_lay[Sub][hh] /= Sd[hh];
            else
              Sd_conc_lay[Sub][hh] = 0.0;
          }
        }
        else {
          avail_evap = Sd[hh];
          Sd[hh] = 0.0;
          Sd_conc_lay[Sub][hh] = 0.0;
        }
        cum_Sd_evap[hh] += avail_evap;
        hru_actet[hh] += avail_evap;
      }
      else
        avail_evap = 0.0;

      avail_evap = hru_evap_buf[hh] - avail_evap - culvert_evapL;

      if(avail_evap > 0.0 && soil_moist[hh] > 0.0 && cov_type[hh] > 0){

        float pctl, pctr, etl, ets, etr;

        if((soil_moist_max[hh] - soil_rechr_max[hh]) > 0.0)
          pctl = (soil_moist[hh] - soil_rechr[hh])/(soil_moist_max[hh] - soil_rechr_max[hh]);
        else
          pctl = 0.0;

        pctr = soil_rechr[hh]/soil_rechr_max[hh];

        etr = avail_evap; // default

        switch (soil_withdrawal_Tables[0][hh]){  // handle recharge layer
          case 1: //******sandy soil
            if(pctr < 0.25)
              etr = 0.5*pctr*avail_evap;
          break;
          case 2: //******loam soil
            if(pctr < 0.5)
              etr = pctr*avail_evap;
          break;
          case 3: //******clay soil
            if(pctr <= 0.33)
              etr = 0.5*pctr*avail_evap;
            else if(pctr < 0.67)
              etr = pctr*avail_evap;
          break;
          case 4: //******organic soil
  //         use default above etr = avail_evap;
          break;
        } // recharge switch

        if(etr > avail_evap){
          etl = 0.0; // default value
          etr = avail_evap;
        }
        else
          etl = avail_evap - etr; // default value

        switch (soil_withdrawal_Tables[1][hh]){  // handle lower layer
          case 1: //******sandy soil
            if(pctl < 0.25)
              etl = 0.5*pctl*etl;
          break;
          case 2: //******loam soil
            if(pctl < 0.5)
              etl = pctl*etl;
          break;
          case 3: //******clay soil
            if(pctl <= 0.33)
              etl = 0.5*pctl*etl;
            else if(pctr < 0.67)
              etl = pctl*etl;
          break;
          case 4: //******organic soil
  //         use default above etl = avail_evap - etr;
          break;
        } // lower switch

    //******Soil moisture accounting  Remember that soil_rechr[hh] includes soil_moist[hh]

        long et_type = cov_type[hh];

  //****** et_type = 0 - no evaporation, bare soil - cov_type = 0
  //****** et_type = 1 - recharge layer only, crops - cov_type = 1
  //****** et_type = 2 - all soil moisture, grasses & shrubs - cov_type = 2

        if(transp_limited[hh] == 1 && et_type == 2)
          et_type = 1;

        et = 0.0;

        switch (et_type){  // handle evaporation
          case 0, -1:  // avail_evap <= 0
            break;
          case 1:
            if(etr > soil_rechr[hh]) {
              soil_rechr[hh] = 0.0;
              et = soil_rechr[hh];
            }
            else{
              soil_rechr[hh] = soil_rechr[hh] - etr;
              et = etr;
            }
            soil_moist[hh] = soil_moist[hh] - et;
            break;
          case 2:
            if(etr + etl >= soil_moist[hh]) {
              et = soil_moist[hh];
              soil_moist[hh] = 0.0;
              soil_rechr[hh] = 0.0;
              soil_lower[hh] = 0.0;
            }
            else {
              et = etr + etl;
              soil_moist[hh] = soil_moist[hh] - et;

              if(etr > soil_rechr[hh]) {
                soil_lower[hh] = soil_lower[hh] - (et - soil_rechr[hh]);
                soil_rechr[hh] = 0.0;
              }
              else{
                soil_rechr[hh] = soil_rechr[hh] - etr;
                soil_lower[hh] = soil_lower[hh] - etl;
              }
            }
            break;
        } // switch

        hru_actet[hh] += et;

      } // soil_moist[hh] > 0.0

      if(soil_moist_max[hh] <= 0.0 && Sdmax[hh] <= 0.0) // assume lake evaporation
        hru_actet[hh] = hru_evap_buf[hh]; // evaporate all

      hru_cum_actet[hh] += hru_actet[hh];
      cum_solute_deposit[hh] += solute_deposit[hh];
    } // for hh
  } // for Sub

  for(hh = 0; chkStruct(); ++hh) {

    if(nstep == 0){ // end of every day
      if(snowinfilDiv > 1) // daily value - ready for next day
         snowinfil_buf[hh] = snowinfil[hh]/snowinfilDiv;

      if(runoffDiv > 1) // daily value - ready for next day
         runoff_buf[hh] = runoff[hh]/runoffDiv;

      if(meltrunoffDiv > 1) // daily value - ready for next day
         meltrunoff_buf[hh] = meltrunoff[hh]/meltrunoffDiv;

      if(evapDiv > 1) // daily value - ready for next day
         hru_evap_buf[hh] = hru_evap[hh]/evapDiv;
    }
  } // for hh
      }  // try
      catch(Sysutils::Exception &E)
      {
        String S = E.Message + " at " +
          Global::DTnow.FormatString("yyyy'/'m'/'d hh':'nn") + " (" + getstep() + ") in '" + Global::OurModulesList->Strings[Global::CurrentModuleRun] +
            "'" + " (" + FloatToStrF(Global::DTnow, ffGeneral, 10, 6) + ") hh = " + FloatToStrF(hh, ffGeneral, 6, 0) + " Sub = " + FloatToStrF(Sub, ffGeneral, 6, 0);
        LogError(S, WARNING);
        if(++FaultsAllowed == 1)
         throw;
      }
}

void ClassWQ_Soil::finish(bool good) {

    float Allcum_soil_runoff = 0.0;
    float Allcum_soil_ssr = 0.0;
    float Allcum_rechr_ssr = 0.0;
    float Allcum_soil_gw = 0.0;
    float Allcum_gw_flow = 0.0;
    float Allcum_infil_act = 0.0;
    float Allcum_soil_moist_change = 0.0;
    float Allcum_Sd_change = 0.0;
    float Allcum_gw_change = 0.0;
    float Allcum_evap = 0.0;
    float Allcum_solute_deposit = 0.0;

    float Allcum_soil_runoff_mWQ = 0.0;
    float Allcum_soil_ssr_mWQ = 0.0;
    float Allcum_rechr_ssr_mWQ = 0.0;
    float Allcum_soil_gw_mWQ = 0.0;
    float Allcum_gw_flow_mWQ = 0.0;
    float Allcum_infil_act_mWQ = 0.0;
    float Allcum_soil_moist_change_mWQ = 0.0;
    float Allcum_Sd_change_mWQ = 0.0;
    float Allcum_gw_change_mWQ = 0.0;

    float AllTotal = 0.0;

    String S = String("H2O");
    LogDebug(S.c_str());
    LogMessage(" ");

    for(hh = 0; chkStruct(); ++hh) {
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr[hh], hru_area[hh], basin_area[0], " *** information only - already included in 'soil_moist'.");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr_change           (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr[hh] - soil_rechr_Init[hh], hru_area[hh], basin_area[0], " *** information only - already included in 'soil_moist'.");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_rechr_ssr               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_rechr_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_moist[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist_change           (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_moist[hh] - soil_moist_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_lower                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_lower[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_lower_change           (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_lower[hh] - soil_bottom_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' hru_cum_actet               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), hru_cum_actet[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_Sd_evap                 (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_Sd_evap[hh], hru_area[hh], basin_area[0], " *** information only - already included in 'hru_cum_actet'.");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_hru_condense            (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_hru_condense[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_infil_act(all)          (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_infil_act[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_gw                 (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_soil_gw[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_runoff             (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_soil_runoff[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_rechr_ssr               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_rechr_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_ssr                (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_soil_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd                          (mm)   (mm*hru) (mm*hru/basin): ").c_str(), Sd[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd_change                   (mm)   (mm*hru) (mm*hru/basin): ").c_str(), Sd[hh] - Sd_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_runoff_to_Sd            (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_runoff_to_Sd[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_Sd_ssr                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_Sd_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_lower_ssr               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_lower_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw                          (mm)   (mm*hru) (mm*hru/basin): ").c_str(), gw[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw_change                   (mm)   (mm*hru) (mm*hru/basin): ").c_str(), gw[hh] - gw_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_gw_flow                 (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_gw_flow[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_redirected_residual     (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_redirected_residual[hh]/hru_area[hh], hru_area[hh], basin_area[0], " *** Added to this HRU surface runoff");
      LogDebug(" ");

      float Total = cum_soil_runoff[hh] + cum_soil_ssr[hh] + cum_soil_gw[hh]
                     + cum_runoff_to_Sd[hh] + cum_gw_flow[hh]
                     + (soil_moist[hh] - soil_moist_Init[hh]) + (Sd[hh] - Sd_Init[hh]) + (gw[hh] - gw_Init[hh]) + hru_cum_actet[hh]
                     - cum_redirected_residual[hh]/hru_area[hh];
      AllTotal += Total;
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Total                       (mm) (mm*hru) (mm*hru/basin): ").c_str(), Total/hru_area[hh], hru_area[hh], basin_area[0], " *** HRU mass balance");
      LogDebug(" ");

      if(variation == VARIATION_1){
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_culvert      (m3) (m3*hru) (m3*hru/basin): ").c_str(), cum_culvert[hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_culvert_over (m3) (m3*hru) (m3*hru/basin): ").c_str(), cum_culvert_over[hh], hru_area[hh], basin_area[0]);
        LogDebug(" ");
      }


      LogDebug(" ");

      Allcum_soil_runoff += cum_soil_runoff[hh]*hru_area[hh];
      Allcum_soil_ssr += cum_soil_ssr[hh]*hru_area[hh];
      Allcum_rechr_ssr += cum_rechr_ssr[hh]*hru_area[hh];
      Allcum_soil_gw += cum_soil_gw[hh]*hru_area[hh];
      Allcum_gw_flow += cum_gw_flow[hh]*hru_area[hh];
      Allcum_infil_act += cum_infil_act[hh]*hru_area[hh];
      Allcum_soil_moist_change += (soil_moist[hh] - soil_moist_Init[hh])*hru_area[hh];
      Allcum_Sd_change += (Sd[hh] - Sd_Init[hh])*hru_area[hh];
      Allcum_gw_change += (gw[hh] - gw_Init[hh])*hru_area[hh];

      Allcum_evap += hru_cum_actet[hh]*hru_area[hh];
    } // hh

    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_runoff (mm*basin):           ").c_str(), Allcum_soil_runoff);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_ssr (mm*basin):              ").c_str(), Allcum_soil_ssr);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_rechr_ssr (mm*basin):             ").c_str(), Allcum_rechr_ssr);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_gw (mm*basin):               ").c_str(), Allcum_soil_gw);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_gw_flow (mm*basin):               ").c_str(), Allcum_gw_flow);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_infil_act (mm*basin):             ").c_str(), Allcum_infil_act);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_moist_change (mm*basin):     ").c_str(), Allcum_soil_moist_change);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_Sd_change (mm*basin):             ").c_str(), Allcum_Sd_change);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_gw_change (mm*basin):             ").c_str(), Allcum_gw_change);
    LogDebug(" ");
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_actet (mm*basin):                 ").c_str(), Allcum_evap);
    LogDebug(" ");

    LogMessage(string("'" + Name + " (Soil_WQ)' AllTotal              (mm*basin):        ").c_str(), AllTotal);
    LogDebug(" ");

  for(Sub = 0; Sub < numsubstances; Sub++){

    String S = String("Substance# ") + (Sub+1);
    LogDebug(S.c_str());
    LogMessage(" ");

    Allcum_soil_runoff = 0.0;
    Allcum_soil_ssr = 0.0;
    Allcum_rechr_ssr = 0.0;
    Allcum_soil_gw = 0.0;
    Allcum_gw_flow = 0.0;
    Allcum_infil_act = 0.0;
    Allcum_soil_moist_change = 0.0;
    Allcum_Sd_change = 0.0;
    Allcum_gw_change = 0.0;
    Allcum_evap = 0.0;
    Allcum_solute_deposit = 0.0;

    Allcum_soil_runoff_mWQ = 0.0;
    Allcum_soil_ssr_mWQ = 0.0;
    Allcum_rechr_ssr_mWQ = 0.0;
    Allcum_soil_gw_mWQ = 0.0;
    Allcum_gw_flow_mWQ = 0.0;
    Allcum_infil_act_mWQ = 0.0;
    Allcum_soil_moist_change_mWQ = 0.0;
    Allcum_Sd_change_mWQ = 0.0;
    Allcum_gw_change_mWQ = 0.0;

    float AllTotal = 0.0;

    for(hh = 0; chkStruct(); ++hh) {
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr[hh], hru_area[hh], basin_area[0], " *** information only - already included in 'soil_moist'.");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr_change           (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr[hh] - soil_rechr_Init[hh], hru_area[hh], basin_area[0], " *** information only - already included in 'soil_moist'.");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr_conc             (mg/l) (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr[hh]*conc_top_lay[Sub][hh], hru_area[hh], basin_area[0]);

      soil_rechr_change_mWQ_lay[Sub][hh] = soil_rechr[hh]*conc_top_lay[Sub][hh] - soil_rechr_Init[hh]*soil_top_conc_Init_lay[Sub][hh];
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_rechr_mWQ_change       (mg)   (mm*hru) (mm*hru/basin): ").c_str(), soil_rechr_change_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_rechr_ssr               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_rechr_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_rechr_ssr_mWQ           (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_rechr_ssr_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_moist[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist_change           (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_moist[hh] - soil_moist_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist_conc             (mg/l) (mm*hru) (mm*hru/basin): ").c_str(), soil_moist[hh]*soil_moist_conc_lay[Sub][hh], hru_area[hh], basin_area[0]);

      soil_moist_change_mWQ_lay[Sub][hh] = soil_moist[hh]*soil_moist_conc_lay[Sub][hh] - soil_moist_Init[hh]*soil_moist_conc_Init_lay[Sub][hh];
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_moist_mWQ_change       (mg)   (mm*hru) (mm*hru/basin): ").c_str(), soil_moist_change_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_lower                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_lower[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_lower_change           (mm)   (mm*hru) (mm*hru/basin): ").c_str(), soil_lower[hh] - soil_bottom_Init[hh], hru_area[hh], basin_area[0]);

      soil_bottom_change_mWQ_lay[Sub][hh] = soil_lower[hh]*conc_bottom_lay[Sub][hh] - soil_bottom_Init[hh]*soil_bottom_conc_Init_lay[Sub][hh];
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' soil_bottom_mWQ_change      (mg)   (mm*hru) (mm*hru/basin): ").c_str(), soil_bottom_change_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' hru_cum_actet               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), hru_cum_actet[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_Sd_evap                 (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_Sd_evap[hh], hru_area[hh], basin_area[0], " *** information only - already included in 'hru_cum_actet'.");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_hru_condense            (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_hru_condense[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_infil_act(all)          (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_infil_act[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_infil_act_mWQ(all)      (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_infil_act_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_gw                 (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_soil_gw[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_gw_mWQ             (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_soil_gw_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_runoff             (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_soil_runoff[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_runoff_mWQ         (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_soil_runoff_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_rechr_ssr               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_rechr_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_rechr_ssr_mWQ           (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_rechr_ssr_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_ssr                (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_soil_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_soil_ssr_mWQ            (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_soil_ssr_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd                          (mm)   (mm*hru) (mm*hru/basin): ").c_str(), Sd[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd_change                   (mm)   (mm*hru) (mm*hru/basin): ").c_str(), Sd[hh] - Sd_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd_conc                     (mm)   (mm*hru) (mm*hru/basin): ").c_str(), Sd[hh]*Sd_conc_lay[Sub][hh], hru_area[hh], basin_area[0]);

      Sd_change_mWQ_lay[Sub][hh] = Sd[hh]*Sd_conc_lay[Sub][hh] - Sd_Init[hh]*Sd_conc_Init[hh];
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Sd_mWQ_change               (mg)   (mm*hru) (mm*hru/basin): ").c_str(), Sd_change_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_runoff_to_Sd            (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_runoff_to_Sd[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_runoff_to_Sd_mWQ        (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_runoff_to_Sd_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_Sd_ssr                  (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_Sd_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_lower_ssr               (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_lower_ssr[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw                          (mm)   (mm*hru) (mm*hru/basin): ").c_str(), gw[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw_change                   (mm)   (mm*hru) (mm*hru/basin): ").c_str(), gw[hh] - gw_Init[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw_conc                     (mg/l) (mm*hru) (mm*hru/basin): ").c_str(), gw[hh]*gw_conc_lay[Sub][hh], hru_area[hh], basin_area[0]);

      gw_change_mWQ_lay[Sub][hh] = gw[hh]*gw_conc_lay[Sub][hh] - gw_Init[hh]*gw_conc_Init[hh];
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' gw_mWQ_change               (mg)   (mm*hru) (mm*hru/basin): ").c_str(), gw_conc_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_gw_flow                 (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_gw_flow[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_gw_flow_mWQ             (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_gw_flow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_redirected_residual     (mm)   (mm*hru) (mm*hru/basin): ").c_str(), cum_redirected_residual[hh]/hru_area[hh], hru_area[hh], basin_area[0], " *** Added to this HRU surface runoff");
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_redirected_residual_mWQ (mg)   (mg*hru) (mg*hru/basin): ").c_str(), cum_redirected_residual_mWQ_lay[Sub][hh]/hru_area[hh], hru_area[hh], basin_area[0], " *** Added to this HRU surface runoff");
      LogDebug(" ");

      float Total = cum_soil_runoff[hh] + cum_soil_ssr[hh] + cum_soil_gw[hh]
                     + cum_runoff_to_Sd[hh] + cum_gw_flow[hh]
                     + (soil_moist[hh] - soil_moist_Init[hh]) + (Sd[hh] - Sd_Init[hh]) + (gw[hh] - gw_Init[hh]) + hru_cum_actet[hh]
                     - cum_redirected_residual[hh]/hru_area[hh];
      AllTotal += Total;
      LogMessageA(hh, string("'" + Name + " (Soil_WQ)' Total                       (mm) (mm*hru) (mm*hru/basin): ").c_str(), Total/hru_area[hh], hru_area[hh], basin_area[0], " *** HRU mass balance");
      LogDebug(" ");

      if(variation == VARIATION_1){
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_culvert      (m3) (m3*hru) (m3*hru/basin): ").c_str(), cum_culvert[hh], hru_area[hh], basin_area[0]);
        LogMessageA(hh, string("'" + Name + " (Soil_WQ)' cum_culvert_over (m3) (m3*hru) (m3*hru/basin): ").c_str(), cum_culvert_over[hh], hru_area[hh], basin_area[0]);
        LogDebug(" ");
      }


      LogDebug(" ");

      Allcum_soil_runoff += cum_soil_runoff[hh]*hru_area[hh];
      Allcum_soil_ssr += cum_soil_ssr[hh]*hru_area[hh];
      Allcum_rechr_ssr += cum_rechr_ssr[hh]*hru_area[hh];
      Allcum_soil_gw += cum_soil_gw[hh]*hru_area[hh];
      Allcum_gw_flow += cum_gw_flow[hh]*hru_area[hh];
      Allcum_infil_act += cum_infil_act[hh]*hru_area[hh];
      Allcum_soil_moist_change += (soil_moist[hh] - soil_moist_Init[hh])*hru_area[hh];
      Allcum_Sd_change += (Sd[hh] - Sd_Init[hh])*hru_area[hh];
      Allcum_gw_change += (gw[hh] - gw_Init[hh])*hru_area[hh];

      Allcum_soil_runoff_mWQ += cum_soil_runoff_mWQ_lay[Sub][hh]*hru_area[hh];
      Allcum_soil_ssr_mWQ += cum_soil_ssr_mWQ_lay[Sub][hh]*hru_area[hh];
      Allcum_rechr_ssr_mWQ += cum_rechr_ssr_mWQ_lay[Sub][hh]*hru_area[hh];
      Allcum_soil_gw_mWQ += cum_soil_gw_mWQ_lay[Sub][hh]*hru_area[hh];
      Allcum_gw_flow_mWQ += cum_gw_flow_mWQ_lay[Sub][hh]*hru_area[hh];
      Allcum_infil_act_mWQ += cum_infil_act_mWQ_lay[Sub][hh]*hru_area[hh];
      Allcum_soil_moist_change_mWQ += (soil_moist[hh]*soil_moist_conc_lay[Sub][hh] - soil_moist_Init[hh]*soil_moist_conc_Init[hh])*hru_area[hh];
      Allcum_Sd_change_mWQ += (Sd[hh]*Sd_conc_lay[Sub][hh] - Sd_Init[hh]*Sd_conc_Init_lay[Sub][hh])*hru_area[hh];
      Allcum_gw_change_mWQ += (gw[hh]*gw_conc_lay[Sub][hh] - gw_Init[hh]*gw_conc_Init_lay[Sub][hh])*hru_area[hh];
      Allcum_evap += hru_cum_actet[hh]*hru_area[hh];
      Allcum_solute_deposit += cum_solute_deposit[hh]*hru_area[hh];
    }

    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_runoff (mm*basin):           ").c_str(), Allcum_soil_runoff);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_runoff_mWQ (mg*basin):       ").c_str(), Allcum_soil_runoff_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_ssr (mm*basin):              ").c_str(), Allcum_soil_ssr);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_ssr_mWQ (mg*basin):          ").c_str(), Allcum_soil_ssr_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_rechr_ssr (mm*basin):             ").c_str(), Allcum_rechr_ssr);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_rechr_ssr_mWQ (mg*basin):         ").c_str(), Allcum_rechr_ssr_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_gw (mm*basin):               ").c_str(), Allcum_soil_gw);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_gw_mWQ (mg*basin):           ").c_str(), Allcum_soil_gw_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_gw_flow (mm*basin):               ").c_str(), Allcum_gw_flow);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_gw_flow_mWQ (mg*basin):           ").c_str(), Allcum_gw_flow_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_infil_act (mm*basin):             ").c_str(), Allcum_infil_act);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_infil_act_mWQ (mg*basin):         ").c_str(), Allcum_infil_act_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_moist_change (mm*basin):     ").c_str(), Allcum_soil_moist_change);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_soil_moist_change_mWQ (mg*basin): ").c_str(), Allcum_soil_moist_change_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_Sd_change (mm*basin):             ").c_str(), Allcum_Sd_change);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_Sd_change_mWQ (mg*basin):         ").c_str(), Allcum_Sd_change_mWQ);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_gw_change (mm*basin):             ").c_str(), Allcum_gw_change);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_gw_change_mWQ (mg*basin):         ").c_str(), Allcum_gw_change_mWQ);
    LogDebug(" ");
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_actet (mm*basin):                 ").c_str(), Allcum_evap);
    LogMessage(string("'" + Name + " (Soil_WQ)' Allcum_solute_deposit (mg*basin):        ").c_str(), Allcum_solute_deposit);
    LogDebug(" ");

    LogMessage(string("'" + Name + " (Soil_WQ)' AllTotal              (mm*basin):        ").c_str(), AllTotal);
    LogDebug(" ");
  } // sub
}

float ClassWQ_Soil::FunctX(const float x){
  float X = 0.0, F;
  for( long n = 1; n < 100; ++n){
    float y = -2.0*n*x;
    F = 4.0*exp(y)/(n*(1.0-exp(y)));
    X += F;
    if(fabs(F) < 0.001)
      return X;
  }
  return X;
}

void ClassWQ_Soil::Set_WQ(const long hru, float *var, float *var_cWQ, float amount, float amount_cWQ){

  var[hru] = amount;
  if(amount > 0.0)
    var_cWQ[hru] = amount_cWQ;
  else
    var_cWQ[hru] = 0.0;
}

void ClassWQ_Soil::copy_array(float *from, float *to){
  for(hh = 0; chkStruct(); ++hh)
    to[hh] = from[hh];
}
/*
void ClassWQ_Soil::copy_basin(float *from, float *to){
  to[0] = from[0];
}*/

void ClassWQ_Soil::restore_hru(float *from, float *to, long hh){
  to[hh] = from[hh];
}

void ClassWQ_Soil::Restore(long hh){

  restore_hru(redirected_residual_0, redirected_residual, hh);
  restore_hru(Sd_0, Sd, hh);
  restore_hru(gw_0, gw, hh);
  restore_hru(soil_rechr_0, soil_rechr, hh);
  restore_hru(soil_moist_0, soil_moist, hh);
  restore_hru(soil_lower_0, soil_lower, hh);
  restore_hru(gw_flow_0, gw_flow, hh);
  restore_hru(hru_actet_0, hru_actet, hh);
  restore_hru(hru_cum_actet_0, hru_cum_actet, hh);
  restore_hru(cum_soil_runoff_0, cum_soil_runoff, hh);
  restore_hru(cum_redirected_residual_0, cum_redirected_residual, hh);
  restore_hru(cum_hru_condense_0, cum_hru_condense, hh);
  restore_hru(cum_Sd_evap_0, cum_Sd_evap, hh);
  restore_hru(cum_Sd_ssr_0, cum_Sd_ssr, hh);
  restore_hru(cum_Sd_gw_0, cum_Sd_gw, hh);
  restore_hru(cum_lower_ssr_0, cum_lower_ssr, hh);
  restore_hru(cum_infil_act_0, cum_infil_act, hh);
  restore_hru(cum_gw_flow_0, cum_gw_flow, hh);
  restore_hru(cum_soil_ssr_0, cum_soil_ssr, hh);
  restore_hru(cum_rechr_ssr_0, cum_rechr_ssr, hh);
  restore_hru(cum_runoff_to_Sd_0, cum_runoff_to_Sd, hh);
  restore_hru(cum_runoff_to_ssr_0, cum_runoff_to_ssr, hh);
  restore_hru(cum_soil_gw_0, cum_soil_gw, hh);
  restore_hru(cum_infil_act_0, cum_infil_act, hh);
  restore_hru(cum_solute_deposit_0, cum_solute_deposit, hh);
}

void ClassWQ_Soil::Save(){

  copy_array(redirected_residual, redirected_residual_0);
  copy_array(Sd, Sd_0);
  copy_array(gw, gw_0);
  copy_array(soil_rechr, soil_rechr_0);
  copy_array(soil_moist, soil_moist_0);
  copy_array(soil_lower, soil_lower_0);
  copy_array(gw_flow, gw_flow_0);
  copy_array(hru_actet, hru_actet_0);
  copy_array(hru_cum_actet, hru_cum_actet_0);
  copy_array(cum_soil_runoff, cum_soil_runoff_0);
  copy_array(cum_redirected_residual, cum_redirected_residual_0);
  copy_array(cum_hru_condense, cum_hru_condense_0);
  copy_array(cum_Sd_evap, cum_Sd_evap_0);
  copy_array(cum_Sd_ssr, cum_Sd_ssr_0);
  copy_array(cum_Sd_gw, cum_Sd_gw_0);
  copy_array(cum_lower_ssr, cum_lower_ssr_0);
  copy_array(cum_gw_flow, cum_gw_flow_0);
  copy_array(cum_soil_ssr, cum_soil_ssr_0);
  copy_array(cum_rechr_ssr, cum_rechr_ssr_0);
  copy_array(cum_runoff_to_Sd, cum_runoff_to_Sd_0);
  copy_array(cum_runoff_to_ssr, cum_runoff_to_ssr_0);
  copy_array(cum_soil_gw, cum_soil_gw_0);
  copy_array(cum_infil_act, cum_infil_act_0);
  copy_array(cum_solute_deposit, cum_solute_deposit_0);
}

void ClassWQ_Soil::Reset_Basin_WQ(long hru, float *var, float *var_cWQ){
  var[hru] = 0.0;
  var_cWQ[hru] = 0.0;
}

void ClassWQ_Soil::Reset_WQ(long hru, float *var, float *var_cWQ){
  var[hru] = 0.0;
  var_cWQ[hru] = 0.0;
}

void ClassWQ_Soil::Reset_WQ(long hru, float *var, float **var_cWQ_lay){
  var[hru] = 0.0;
  for(long Sub = 0; Sub < numsubstances; Sub++){
    var_cWQ_lay[Sub][hru] = 0.0;
  }
}

ClassWQ_Netroute* ClassWQ_Netroute::klone(string name) const{
  return new ClassWQ_Netroute(name);
}

void ClassWQ_Netroute::decl(void) {

// kg/km2 = (mg/l)*mm = mg/m2

  Description = "'Handles the routing of surface runoff, subsurface runoff and HRU routing using the lag and route method.'\
                    'uses Muskingum,' \
                    'uses Clark.'";

  declvar("inflow", NHRU, "inflow from other HRUs", "(mm*km^2/int)", &inflow);

  declvar("inflow_mWQ", NDEFN, "Mass: inflow from other HRUs", "(mg)", &inflow_mWQ, &inflow_mWQ_lay, numsubstances);

  declstatvar("cuminflow", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &cuminflow);

  declstatvar("cuminflow_mWQ", NDEFN, "cumulative mass of solute inflow from other HRUs", "(mg)", &cuminflow_mWQ, &cuminflow_mWQ_lay, numsubstances);

  declvar("outflow", NHRU, "HRU outflow", "(mm*km^2/int)", &outflow);

  declvar("outflow_mWQ", NDEFN, "Concentration: HRU outflow", "(mg/int)", &outflow_mWQ, &outflow_mWQ_lay, numsubstances);

  declstatvar("cumoutflow", NHRU, "cumulative HRU outflow", "(mm*km^2)", &cumoutflow);

  declstatvar("cumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU outflow", "(mg)", &cumoutflow_mWQ, &cumoutflow_mWQ_lay, numsubstances);

  decldiag("outflow_diverted", NHRU, "HRU outflow diverted to another HRU", "(mm*km^2/int)", &outflow_diverted);

  decldiag("outflow_diverted_conc", NDEFN, "Concentration: HRU outflow diverted to another HRU", "(mg/l)", &outflow_diverted_conc, &outflow_diverted_conc_lay, numsubstances);

  declstatvar("cumoutflow_diverted", NHRU, "cumulative HRU outflow diverted to another HRU", "(mm*km^2/int)", &cumoutflow_diverted);

  declstatvar("cumoutflow_diverted_mWQ", NDEFN, "cumulative mass of solute HRU outflow diverted to another HRU", "(mg)", &cumoutflow_diverted_mWQ, &cumoutflow_diverted_mWQ_lay, numsubstances);

  declstatvar("cum_to_Sd", NHRU, "cumulative other HRU to depressional storage (Sd) of this HRU", "(mm)", &cum_to_Sd);

  declstatvar("cum_to_Sd_mWQ", NDEFN, "cumulative mass of solute from other HRU to depressional storage (Sd) of this HRU", "(mg)", &cum_to_Sd_mWQ, &cum_to_Sd_mWQ_lay, numsubstances);

  declstatvar("cum_to_soil_rechr", NHRU, "cumulative other HRU to soil_rechr of this HRU", "(mm)", &cum_to_soil_rechr);

  declstatvar("cum_to_soil_rechr_mWQ", NDEFN, "cumulative mass of solute from other HRU to soil_rechr of this HRU", "(mg)", &cum_to_soil_rechr_mWQ, &cum_to_soil_rechr_mWQ_lay, numsubstances);

  declvar("gwinflow", NHRU, "ground water inflow", "(mm*km^2/int)", &gwinflow);

  declvar("gwinflow_mWQ", NDEFN, "Concentration: ground water inflow", "(mg/int)", &gwinflow_mWQ, &gwinflow_mWQ_lay, numsubstances);

  declstatvar("gwcuminflow", NHRU, "cumulative gw inflow", "(mm*km^2)", &gwcuminflow);

  declstatvar("gwcuminflow_mWQ", NDEFN, "cumulative mass of solute gw inflow", "(mg)", &gwcuminflow_mWQ, &gwcuminflow_mWQ_lay, numsubstances);

  declvar("gwoutflow", NHRU, "HRU gw outflow", "(mm*km^2/int)", &gwoutflow);

  declvar("gwoutflow_mWQ", NDEFN, "Concentration: HRU gw outflow", "(mg/int)", &gwoutflow_mWQ, &gwoutflow_mWQ_lay, numsubstances);

  declstatvar("gwcumoutflow", NHRU, "cumulative HRU gw outflow", "(mm*km^2)", &gwcumoutflow);

  declstatvar("gwcumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU gw outflow", "(mg)", &gwcumoutflow_mWQ, &gwcumoutflow_mWQ_lay, numsubstances);

  decldiag("gwoutflow_diverted", NHRU, "HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwoutflow_diverted);

  decldiag("gwoutflow_diverted_conc", NDEFN, "HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwoutflow_diverted_conc, &gwoutflow_diverted_conc_lay, numsubstances);

  declstatvar("gwcumoutflow_diverted", NHRU, "cumulative HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwcumoutflow_diverted);

  declstatvar("gwcumoutflow_diverted_mWQ", NDEFN, "cumulative mass of solute HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwcumoutflow_diverted_mWQ, &gwcumoutflow_diverted_mWQ_lay, numsubstances);

  declvar("ssrinflow", NHRU, "inflow from other HRUs", "(mm*km^2/int)", &ssrinflow);

  declvar("ssrinflow_mWQ", NDEFN, "Concentration: inflow from other HRUs", "(mg/l)", &ssrinflow_mWQ, &ssrinflow_mWQ_lay, numsubstances);

  declstatvar("ssrcuminflow", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &ssrcuminflow);

  declstatvar("ssrcuminflow_mWQ", NDEFN, "cumulative mass of solute of inflow from other HRUs", "(mg)", &ssrcuminflow_mWQ, &ssrcuminflow_mWQ_lay, numsubstances);

  declvar("ssroutflow", NHRU, "HRU outflow", "(mm*km^2/int)", &ssroutflow);

  declvar("ssroutflow_mWQ", NDEFN, "Concentration: HRU outflow", "(mg/int)", &ssroutflow_mWQ, &ssroutflow_mWQ_lay, numsubstances);

  declstatvar("ssrcumoutflow", NHRU, "cumulative HRU outflow", "(mm*km^2)", &ssrcumoutflow);

  declstatvar("ssrcumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU outflow", "(mg)", &ssrcumoutflow_mWQ, &ssrcumoutflow_mWQ_lay, numsubstances);

  declstatvar("HRU_cumbasinflow", NHRU, "cumulative HRU to basinflow", "(mm*km^2)", &HRU_cumbasinflow);

  declstatvar("HRU_cumbasinflow_mWQ", NDEFN, "cumulative HRU to basinflow", "(mg*km^2)", &HRU_cumbasinflow_mWQ, &HRU_cumbasinflow_mWQ_lay, numsubstances);

  declvar("runinflow", NHRU, "inflow from other HRUs", "(mm*km^2/int)", &runinflow);

  declvar("runinflow_mWQ", NDEFN, "Concentration: inflow from other HRUs", "(mg/int)", &runinflow_mWQ, &runinflow_mWQ_lay, numsubstances);

  declstatvar("runcuminflow", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &runcuminflow);

  declstatvar("runcuminflow_mWQ", NDEFN, "cumulative mass of solute inflow from other HRUs", "(mg)", &runcuminflow_mWQ, &runcuminflow_mWQ_lay, numsubstances);

  declvar("runoutflow", NHRU, "HRU outflow", "(mm*km^2/int)", &runoutflow);

  declvar("runoutflow_mWQ", NDEFN, "Concentration: HRU outflow", "(mg/l)", &runoutflow_mWQ, &runoutflow_mWQ_lay, numsubstances);

  declstatvar("runcumoutflow", NHRU, "cumulative HRU outflow", "(mm*km^2)", &runcumoutflow);

  declstatvar("runcumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU outflow", "(mg)", &runcumoutflow_mWQ, &runcumoutflow_mWQ_lay, numsubstances);

  declstatvar("cum_preferential_flow_to_gw", NHRU, "cumulative other HRU's runoff to gw of this HRU via preferential flow path", "(mm)", &cum_preferential_flow_to_gw);


  declvar("basinflow", BASIN, "basin surface and sub-surface outflow", "(m^3/int)", &basinflow);

  declvar("basinflow_conc", NDEF, "basin surface and sub-surface outflow", "(mg/l)", &basinflow_conc, &basinflow_conc_lay, numsubstances);

  decldiag("basinflow_s", BASIN, "basin surface and sub-surface outflow", "(m^3/s)", &basinflow_s);

  declvar("cumbasinflow", BASIN, "cumulative basin surface and sub-surface outflow", "(m^3)", &cumbasinflow);

  declvar("cumbasinflow_mWQ", NDEF, "cumulative mass of solute basin surface and sub-surface outflow", "(m^3)", &cumbasinflow_mWQ, &cumbasinflow_mWQ_lay, numsubstances);

  declvar("basingw", BASIN, "cumulative basin groundwater outflow", "(m^3/int)", &basingw);

  declvar("basingw_conc", NDEF, "cumulative basin groundwater outflow", "(m^3/int)", &basingw_conc, &basingw_conc_lay, numsubstances);

//  declvar("Used_mWQ", NDEF, "directed to basinbasin surface and sub-surface outflow", "(mg/int)", &Used_mWQ, &Used_mWQ_lay, numsubstances);

  decldiag("basingw_s", BASIN, "cumulative basin groundwater outflow", "(m^3/s)", &basingw_s);

  declstatvar("cumbasingw", BASIN, "cumulative basin groundwater outflow", "(m^3)", &cumbasingw);

  declvar("cumbasingw_mWQ", NDEF, "cumulative mass of solute basin groundwater outflow", "(m^3)", &cumbasingw_mWQ, &cumbasingw_mWQ_lay, numsubstances);


  decllocal("soil_ssr_Buf", NHRU, "buffer subsurface runoff", "(mm/d)", &soil_ssr_Buf);

  declvar("soil_ssr_Buf_conc", NDEFN, "buffer subsurface runoff", "(mm/d)", &soil_ssr_Buf_conc, &soil_ssr_Buf_conc_lay, numsubstances);

  decllocal("soil_runoff_Buf", NHRU, "buffer rain runoff", "(mm/d)", &soil_runoff_Buf);

  declvar("soil_runoff_Buf_conc", NDEFN, "buffer rain runoff", "(mm/d)", &soil_runoff_Buf_conc, &soil_runoff_Buf_conc_lay, numsubstances);

  decllocal("soil_gw_Buf", NHRU, "buffer soil_gw(gw_flow) runoff", "(mm/d)", &soil_gw_Buf);

  declvar("soil_gw_Buf_conc", NDEFN, "buffer soil_gw(gw_flow) runoff", "(mm/d)", &soil_gw_Buf_conc, &soil_gw_Buf_conc_lay, numsubstances);


  declparam("basin_area", BASIN, "3", "1e-6", "1e09", "Total basin area", "(km^2)", &basin_area);

  declparam("hru_area", NHRU, "[1]", "1e-6", "1e09", "HRU area", "(km^2)", &hru_area);

  declparam("Kstorage", NHRU, "[0.0]", "0.0","200.0", "aggregated storage constant", "(d)", &Kstorage);

  declparam("Lag", NHRU, "[0.0]", "0.0","1.0E4.0", "aggregated lag delay", "(h)", &Lag);

  declparam("ssrKstorage", NHRU, "[0.0]", "0.0","200.0", "subsurface runoff storage constant", "(d)", &ssrKstorage);

  declparam("ssrLag", NHRU, "[0.0]", "0.0","1.0E4.0", "subsurface runoff lag delay", "(h)", &ssrLag);

  declparam("runKstorage", NHRU, "[0.0]", "0.0","200.0", "runoff storage constant", "(d)", &runKstorage);

  declparam("runLag", NHRU, "[0.0]", "0.0","1.0E4", "runoff lag delay", "(h)", &runLag);

  declparam("gwKstorage", NHRU, "[0.0]", "0.0","200.0", "gw storage constant", "(d)", &gwKstorage);

  declparam("gwLag", NHRU, "[0.0]", "0.0","1.0E4", "gw lag delay", "(h)", &gwLag);

  declparam("order", NHRU, "[1,2,3,4,5!]", "1","1000", "HRU routing process order", "()", &order);

  declparam("whereto", NHRU, "[0]", "0", "1000", "send to; 0 - basin outflow, or HRU input", "()", &whereto);

  declparam("gwwhereto", NHRU, "[0]", "-1000", "1000", "send to: 0 - basingw, >0 - other HRU surface input <0 - other abs(-HRU) gw input, or (< -HRUmax or > +HRUmax) - surface basinflow", "()", &gwwhereto);

  declparam("Sdmax", NHRU, "[0]", "0.0", "1000.0","Maximum depression storage", "(mm)", &Sdmax);

  declparam("soil_rechr_max", NHRU, "[60.0]", "0.0", "350.0", "soil recharge maximum (<= soil_moist_max).", "(mm)", &soil_rechr_max);

  declparam("Sd_ByPass", NHRU, "[0]", "0", "1","0 - normal, 1 - Bypass Pond/Depressional storage (i.e. Sd).", "()", &Sd_ByPass);

  declparam("soil_rechr_ByPass", NHRU, "[1]", "0", "1","0 - normal, 1 - Bypass recharge layer (i.e. soil_rechr).", "()", &soil_rechr_ByPass);

  declparam("preferential_flow", NHRU, "[0]", "0", "1","0 - no preferential and remain as runoff routing to other HRU, 1 - preferential flow and route runoff to other HRU's gw.", "()", &preferential_flow);


  soil_gwDiv = declgetvar("*", "gw_flow", "(mm/int)", &soil_gw);

  soil_ssrDiv = declgetvar("*", "soil_ssr", "(mm/int)", &soil_ssr);

  declgetvar("*", "soil_ssr_conc", "(mg)", &soil_ssr_conc, &soil_ssr_conc_lay);

  soil_runoffDiv = declgetvar("*", "soil_runoff", "(mm/int)", &soil_runoff);

  declgetvar("*", "soil_gw_conc", "(mg)", &soil_gw_conc, &soil_gw_conc_lay);

  declgetvar("*", "soil_ssr_conc", "(mg)", &soil_ssr_conc, &soil_ssr_conc_lay);

  declgetvar("*", "soil_runoff_conc", "(mg)", &soil_runoff_conc, &soil_runoff_conc_lay);


  declputvar("*", "Sd", "(mm)", &Sd);

  declputvar("*", "Sd_conc", "(mg)", &Sd_conc, &Sd_conc_lay);

  declputvar("*", "soil_moist", "(mm)", &soil_moist);

  declputvar("*", "soil_moist_conc", "(mg)", &soil_moist_conc, &soil_moist_conc_lay);

  declputvar("*", "soil_lower", "(mm)", &soil_lower);

  declputvar("*", "soil_rechr", "(mm)", &soil_rechr);

  declputvar("*", "redirected_residual", "(mg)", &redirected_residual);

  declputvar("*", "redirected_residual_conc", "(mm*km^2/int)", &redirected_residual_conc, &redirected_residual_conc_lay);

  declputvar("*", "cum_redirected_residual", "(mg)", &cum_redirected_residual);

  declputvar("*", "cum_redirected_residual_mWQ", "(mg)", &cum_redirected_residual_mWQ, &cum_redirected_residual_mWQ_lay);

  declputvar("*", "gw", "(mm)", &gw);

  declputvar("*", "gw_conc", "(mg)", &gw_conc, &gw_conc_lay);

  declputvar("*", "conc_top", "(mg/l)", &conc_top, &conc_top_lay);

  declputvar("*", "conc_bottom", "(mg/l)", &conc_bottom, &conc_bottom_lay);


  variation_set = VARIATION_0;

  decllocal("Ktravel", NHRU, "travel time", "(d)", &Ktravel);

  declparam("route_n", NHRU, "[0.025]", "0.016","0.2", "Manning roughness coefficient", "()", &route_n);

  declparam("route_R", NHRU, "[0.5]", "0.01","1.0E4", "hydraulic radius", "(m)", &route_R);

  declparam("route_S0", NHRU, "[1e-3]", "1e-6","1.0", "longitudinal channel slope", "()", &route_S0);

  declparam("route_L", NHRU, "[200.0]", "0.01","1.0E10", "routing length", "(m)", &route_L);

  declparam("route_X_M", NHRU, "[0.25]", "0.0","0.5", "dimensionless weighting factor", "()", &route_X_M);

  declparam("Channel_shp", NHRU, "[0]", "0", "2", "rectangular - 0/parabolic - 1/triangular - 2", "()", &route_Cshp);


  variation_set = VARIATION_1;

  declparam("Kstorage", NHRU, "[0.0]", "0.0","200.0", "aggregated storage constant", "(d)", &Kstorage);


  variation_set = VARIATION_ORG;

  decllocal("outflow_0", NHRU, "HRU outflow", "(mm*km^2/int)", &outflow_0);

  declvar("gwoutflow_0", NHRU, "HRU gw outflow", "(mm*km^2/int)", &gwoutflow_0);

  decllocal("redirected_residual_0", NHRU, "", "", &redirected_residual_0);

  decllocal("Sd_0", NHRU, "Depression storage.", "(mm)", &Sd_0);

  decllocal("gw_0", NHRU, "ground water storage.", "(mm)", &gw_0);

  decllocal("soil_rechr_0", NHRU, "moisture content of soil recharge zone.", "(mm)", &soil_rechr_0);

  decllocal("soil_moist_0", NHRU, "moisture content of soil profile to the depth.", "(mm)", &soil_moist_0);

  decllocal("soil_lower_0", NHRU, "moisture content of soil profile to the depth.", "(mm)", &soil_lower_0);

  decllocal("cum_to_Sd_0", NHRU, "cumulative other HRU to depressional storage (Sd) of this HRU", "(mm)", &cum_to_Sd_0);

  decllocal("basinflow_0", BASIN, "basin surface and sub-surface outflow", "(m^3/int)", &basinflow_0);

  decllocal("basingw_0", BASIN, "cumulative basin groundwater outflow", "(m^3/int)", &basingw_0);

  decllocal("cumbasinflow_0", BASIN, "cumulative basin surface and sub-surface outflow", "(m^3)", &cumbasinflow_0);

  decllocal("cumbasingw_0", BASIN, "cumulative basin groundwater outflow", "(m^3)", &cumbasingw_0);

   decllocal("cuminflow_0", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &cuminflow_0);

  decllocal("cumoutflow_0", NHRU, "cumulative HRU outflow", "(mm*km^2)", &cumoutflow_0);

  decllocal("cum_to_Sd_0", NHRU, "cumulative other HRU to depressional storage (Sd) of this HRU", "(mm)", &cum_to_Sd_0);

  decllocal("cum_to_soil_rechr_0", NHRU, "cumulative other HRU to soil_rechr of this HRU", "(mm)", &cum_to_soil_rechr_0);

  decllocal("gwcuminflow_0", NHRU, "cumulative gw inflow", "(mm*km^2)", &gwcuminflow_0);

  decllocal("gwcumoutflow_0", NHRU, "cumulative HRU gw outflow", "(mm*km^2)", &gwcumoutflow_0);

  decllocal("gwcumoutflow_diverted_0", NHRU, "cumulative HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwcumoutflow_diverted_0);

  decllocal("ssrcuminflow_0", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &ssrcuminflow_0);

  decllocal("cumoutflow_diverted_0", NHRU, "cumulative HRU outflow diverted to another HRU", "(mm*km^2/int)", &cumoutflow_diverted_0);

  decllocal("ssrcumoutflow_0", NHRU, "cumulative HRU outflow", "(mm*km^2)", &ssrcumoutflow_0);

  decllocal("runcuminflow_0", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &runcuminflow_0);

  decllocal("runcumoutflow_0", NHRU, "cumulative HRU outflow", "(mm*km^2)", &runcumoutflow_0);

  decllocal("cum_preferential_flow_to_gw_0", NHRU, "cumulative other HRU's runoff to gw of this HRU via preferential flow path", "(mm)", &cum_preferential_flow_to_gw_0);

  decllocal("cum_redirected_residual_0", NHRU, "cumulative HRU redirected_residual to another HRU.", "(mm*km^2)", &cum_redirected_residual_0);

  decllocal("HRU_cumbasinflow_0", NHRU, "cumulative HRU to basinflow.", "(mm*km^2)", &HRU_cumbasinflow_0);

}

void ClassWQ_Netroute::init(void) {

  nhru = getdim(NHRU);

  try {
    hruDelay_cWQ = new ClassMuskingum*[numsubstances]; // [numsubstances][nhru] handled locally
    Clark_hruDelay_cWQ = new ClassClark*[numsubstances];
    ssrDelay_cWQ = new ClassClark*[numsubstances];
    runDelay_cWQ = new ClassClark*[numsubstances];
    gwDelay_cWQ = new ClassClark*[numsubstances];
  }
  catch (std::bad_alloc) {
    CRHMException Except("Could not allocate in module CRACK." ,TERMINATE);
    LogError(Except);
    throw Except;
  }

  if(soil_ssrDiv > 1){
    String S = "Netroute:  \"soil_ssr\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(soil_runoffDiv > 1){
    String S = "Netroute:  \"soil_runoff\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(soil_gwDiv > 1){
    String S = "Netroute:  \"gw_flow\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(variation == VARIATION_ORG){
    const float Vw[3] = {1.67, 1.44, 1.33}; // rectangular - 0/parabolic - 1/triangular - 2

    for(hh = 0; hh < nhru; ++hh){
      float Vavg = (1.0/route_n[hh])*pow(route_R[hh], 2.0/3.0)*pow(route_S0[hh], 0.5f); // (m/s)
      Ktravel[hh] = route_L[hh]/(Vw[route_Cshp[hh]]*Vavg)/86400.0; // (d)
    }

    hruDelay = new ClassMuskingum(inflow, outflow, Ktravel, route_X_M, Lag, nhru);
    for(hh = 0; hh < nhru; ++hh){
      if(Ktravel[hh] >= (Global::Interval/(2.0*route_X_M[hh]))){
        String S = string("'" + Name + " (Netroute_M_D) Muskingum coefficient negative in HRU ").c_str() + IntToStr(hh+1);
        CRHMException TExcept(S.c_str(), WARNING);
        LogError(TExcept);
      }

      if(Ktravel[hh] < (Global::Interval/(2.0*(1.0-route_X_M[hh])))){ //    if(hruDelay->c0[hh] < 0.0)
        hruDelay->c0[hh] = 0.0;
        hruDelay->c1[hh] = 1.0;
        hruDelay->c2[hh] = 0.0;
      }
    }
  }
  else if(variation == VARIATION_1)
    Clark_hruDelay = new ClassClark(inflow, outflow, Kstorage, Lag, nhru);

  ssrDelay = new ClassClark(ssrinflow, ssroutflow, ssrKstorage, ssrLag, nhru);
  runDelay = new ClassClark(runinflow, runoutflow, runKstorage, runLag, nhru);
  gwDelay = new ClassClark(gwinflow, gwoutflow, gwKstorage, gwLag, nhru);

  Reset_Basin_WQ(0, basinflow, basinflow_conc);
  Reset_Basin_WQ(0, cumbasinflow, cumbasinflow_mWQ);
  Reset_Basin_WQ(0, basingw, basingw_conc);
  Reset_Basin_WQ(0, cumbasingw, cumbasingw_mWQ);

  basinflow_s[0] = 0.0;
  basingw_s[0] = 0.0;

  for(Sub = 0; Sub < numsubstances; ++Sub){

    if(variation == VARIATION_1)
      Clark_hruDelay_cWQ[Sub] = new ClassClark(inflow_mWQ_lay[Sub], outflow_mWQ_lay[Sub], Kstorage, Lag, nhru);
    else // ClassMuskingum
      hruDelay_cWQ[Sub] = new ClassMuskingum(inflow_mWQ_lay[Sub], outflow_mWQ_lay[Sub], Ktravel, route_X_M, Lag, nhru);

    ssrDelay_cWQ[Sub] = new ClassClark(ssrinflow_mWQ_lay[Sub], ssroutflow_mWQ_lay[Sub], ssrKstorage, ssrLag, nhru, -1);
    runDelay_cWQ[Sub] = new ClassClark(runinflow_mWQ_lay[Sub], runoutflow_mWQ_lay[Sub], runKstorage, runLag, nhru, -1);
    gwDelay_cWQ[Sub] = new ClassClark(gwinflow_mWQ_lay[Sub], gwoutflow_mWQ_lay[Sub], gwKstorage, gwLag, nhru, -1);
  } // for Sub

  for(hh = 0; hh < nhru; ++hh) {
    Reset_WQ(hh, inflow, inflow_mWQ_lay);
    Reset_WQ(hh, cuminflow, cuminflow_mWQ_lay);

    Reset_WQ(hh, outflow, outflow_mWQ_lay);
    Reset_WQ(hh, cumoutflow, cumoutflow_mWQ_lay);

    Reset_WQ(hh, gwinflow, gwinflow_mWQ_lay); ;
    Reset_WQ(hh, gwcuminflow, gwcuminflow_mWQ_lay); ;

    Reset_WQ(hh, gwoutflow, gwoutflow_mWQ_lay);
    Reset_WQ(hh, gwcumoutflow, gwcumoutflow_mWQ_lay);

    Reset_WQ(hh, ssrinflow, ssrinflow_mWQ_lay);
    Reset_WQ(hh, ssrcuminflow, ssrcuminflow_mWQ_lay);

    Reset_WQ(hh, ssroutflow, ssroutflow_mWQ_lay);
    Reset_WQ(hh, ssrcumoutflow, ssrcumoutflow_mWQ_lay);

    Reset_WQ(hh, runinflow, runinflow_mWQ_lay);
    Reset_WQ(hh, runcuminflow, runcuminflow_mWQ_lay);

    Reset_WQ(hh, runoutflow, runoutflow_mWQ_lay);
    Reset_WQ(hh, runcumoutflow, runcumoutflow_mWQ_lay);

    Reset_WQ(hh, outflow_diverted, outflow_diverted_conc_lay);
    Reset_WQ(hh, cumoutflow_diverted, cumoutflow_diverted_mWQ_lay);

    Reset_WQ(hh, gwoutflow_diverted, gwoutflow_diverted_conc_lay);
    Reset_WQ(hh, gwcumoutflow_diverted, gwcumoutflow_diverted_mWQ_lay);

    Reset_WQ(hh, cum_to_Sd, cum_to_Sd_mWQ_lay);

    Reset_WQ(hh, cum_to_soil_rechr, cum_to_soil_rechr_mWQ_lay);

    Reset_WQ(hh, HRU_cumbasinflow, HRU_cumbasinflow_mWQ_lay);

    cum_preferential_flow_to_gw[hh] = 0.0;
    soil_ssr_Buf[hh] = 0.0;
    soil_runoff_Buf[hh] = 0.0;
    soil_gw_Buf[hh] = 0.0;


    boolean OK = false;
    for(long jj = 0; chkStruct(jj); ++jj)
      if(order[jj] - 1 == hh){
        OK = true;
        break;
      }

    if(!OK){
        string SS = string("'" + Name + " (Netroute)' the 'order' parameter does not have a unique value for each HRU");
        CRHMException Except(SS.c_str() ,ERR);
        LogError(Except);
        throw Except;
    }
  } // for hh
}

void ClassWQ_Netroute::run(void) {

  long step = getstep();
  long nstep = step% Global::Freq;

  float gw_Amount = 0.0;
  float gw_Amount_mWQ = 0.0;

  float Amount = 0.0;
  float Amount_mWQ = 0.0;

  Reset_WQ(0, basinflow, basinflow_conc_lay);
  Reset_WQ(0, basingw, basingw_conc_lay);

  for(long hh = 0; chkStruct(hh); ++hh){ // HRUs not in sequence
    for(Sub = 0; Sub < numsubstances; Sub++){

      if(soil_gwDiv == 1){ // interval value
         soil_gw_Buf[hh] = soil_gw[hh];
         soil_gw_Buf_conc_lay[Sub][hh] = soil_gw_conc_lay[Sub][hh];
      }
      if(soil_ssrDiv == 1){ // interval value
         soil_ssr_Buf[hh] = soil_ssr[hh];
         soil_ssr_Buf_conc_lay[Sub][hh] = soil_ssr_conc_lay[Sub][hh];
      }

      if(soil_runoffDiv == 1){ // interval value
         soil_runoff_Buf[hh] = soil_runoff[hh];
         soil_runoff_Buf_conc_lay[Sub][hh] = soil_runoff_conc_lay[Sub][hh];
      }
    } // Sub
  } // hh

  for(Sub = 0; Sub < numsubstances; ++Sub){
    for(long jj = 0; chkStruct(jj); ++jj){ // HRUs not in sequence
      for(hh = 0; chkStruct(hh); ++hh)
        if(order[hh] - 1 == jj)
          break;

      gw_Amount = 0.0;
      gw_Amount_mWQ = 0.0;

      Amount = 0.0;
      Amount_mWQ = 0.0;

      gwinflow[hh] = soil_gw_Buf[hh]*hru_area[hh];

      gwoutflow_diverted[hh] = 0.0;

      gwinflow_mWQ[hh] = soil_gw_conc_lay[Sub][hh]*soil_gw_Buf[hh]*hru_area[hh];
      gwoutflow_diverted_conc_lay[Sub][hh] = 0.0;

      for(long hhh = 0; chkStruct(hhh); ++hhh) {
        if(gwoutflow[hhh] > 0.0 && gwwhereto[hhh] && (abs(gwwhereto[hhh])-1 == hh || abs(gwwhereto[hhh]) > nhru)){ // handles "gwwhereto" <> 0
            gwoutflow_diverted[hhh] = gwoutflow[hhh]; // gwoutflow_diverted[hh] = gwoutflow[hhh];

            gw_Amount = gwoutflow[hhh]/hru_area[hh]; // units (mm*km^2/int)
            gw_Amount_mWQ = gwoutflow_mWQ_lay[Sub][hhh]/hru_area[hh]; // units (mm*km^2/int)

            gwcumoutflow_diverted[hhh] += gw_Amount;
            gwcumoutflow_diverted_mWQ_lay[Sub][hhh] += gw_Amount_mWQ;

            if(abs(gwwhereto[hhh]) <= nhru){
              if(gwwhereto[hhh] > 0){ // direct to HRU surface
                float free = soil_rechr_max[hh] - soil_rechr[hh];
                float free_mWQ = Amount_mWQ*free/gw_Amount;

                if(free > 0.0 && !soil_rechr_ByPass[hh]){
                  if(free > gw_Amount){ // units (mm*km^2/int)

                    conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + Amount_mWQ;
                    soil_moist_conc_lay[Sub][hh] = soil_moist_conc_lay[Sub][hh]*soil_moist[hh] + Amount_mWQ;

                    conc_top_lay[Sub][hh] /= (soil_rechr[hh] + gw_Amount);
                    soil_moist_conc_lay[Sub][hh] /= (soil_moist[hh] + gw_Amount);

                    cum_to_soil_rechr_mWQ_lay[Sub][hh] += gw_Amount_mWQ;

                    if(Sub == numsubstances-1){
                      soil_rechr[hh] += gw_Amount;
                      soil_moist[hh] += gw_Amount;
                      cum_to_soil_rechr[hh] += gw_Amount;
                    }

                    gw_Amount = 0.0;
                    gw_Amount_mWQ = 0.0;
                  }
                  else {
                    gw_Amount_mWQ = gw_Amount_mWQ*(gw_Amount - free)/Amount;
                    gw_Amount = (gw_Amount - free);

                    conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + free_mWQ;
                    soil_moist_conc_lay[Sub][hh] = soil_moist_conc_lay[Sub][hh]*soil_moist[hh] + free_mWQ;

                    conc_top_lay[Sub][hh] /= soil_rechr_max[hh];
                    soil_moist_conc_lay[Sub][hh] /= (soil_moist[hh] + gw_Amount);

                    cum_to_soil_rechr_mWQ_lay[Sub][hh] += free_mWQ;

                    if(Sub == numsubstances-1){
                      soil_rechr[hh] = soil_rechr_max[hh];
                      soil_moist[hh] = soil_moist[hh] + free_mWQ;
                      cum_to_soil_rechr[hh] += free;
                    }
                  }
                }

                free = Sdmax[hh] - Sd[hh];
                free_mWQ = Amount_mWQ*free/gw_Amount;

                if(free > 0.0 && !Sd_ByPass[hh] && gw_Amount > 0.0){
                  if(free > gw_Amount){ // units (mm*km^2/int)
                    Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + Amount_mWQ;

                    Sd_conc_lay[Sub][hh] /= (Sd[hh] + gw_Amount);
                    cum_to_Sd_mWQ_lay[Sub][hh] += Amount_mWQ;

                    if(Sub == numsubstances-1){
                      Sd[hh] += gw_Amount;
                      cum_to_Sd[hh] += gw_Amount;
                    }

                    gw_Amount = 0.0;
                    gw_Amount_mWQ = 0.0;
                  }
                  else {
                    gw_Amount_mWQ = gw_Amount_mWQ*(gw_Amount - free)/Amount;

                    Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + free_mWQ;
                    cum_to_Sd_mWQ_lay[Sub][hh] += free_mWQ;
                    Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]/Sdmax[hh];

                    if(Sub == numsubstances-1){
                      gw_Amount = (gw_Amount - free);
                      Sd[hh] = Sdmax[hh];
                      cum_to_Sd[hh] += free;
                    }
                  }
                }
              } // hh > 0
              else{ // hh < 0
                gw_conc_lay[Sub][hh] = gw_conc_lay[Sub][hh]*gw[hh] + gw_Amount_mWQ*gw_Amount;
                gw[hh] += gw_Amount;
                gw_conc_lay[Sub][hh] /= gw[hh];

                gw_Amount = 0.0;
                gw_Amount_mWQ = 0.0;
              }
            }
            else if(gwwhereto[hh] == 0){ // move to basin gw
              basingw_conc_lay[Sub][0] = basingw_conc_lay[Sub][0]*basingw[0] + gwoutflow_mWQ_lay[Sub][hh]*1000;

              if(basingw[0] > 0.0)
                basingw_conc_lay[Sub][0] /= (basingw[0] + gw_Amount*1000);
              else
                basingw_conc_lay[Sub][0] = 0.0;

              gwcumoutflow_mWQ_lay[Sub][hh] += gw_Amount_mWQ;

              if(Sub == numsubstances-1){
                basingw[0] += gw_Amount*1000; // (m3) end of every day
                gwcumoutflow[hh] += gw_Amount;
              }

              gw_Amount = 0.0;
              gw_Amount_mWQ = 0.0;
            }
            else{
              if(Sub == numsubstances-1){
               HRU_cumbasinflow[hh] +=  gw_Amount;
                basinflow[0] += gw_Amount*hru_area[hh]*1000; // (m3)
                cumoutflow[hh] += gw_Amount*hru_area[hh];
              }
              gw_Amount = 0.0;
              gw_Amount_mWQ = 0.0;
            } // is HRU in range
          } // handles "gwwhereto" <> 0
      } // for hhh

      gwcuminflow[hh] += gwinflow[hh];

      runinflow[hh] = soil_runoff_Buf[hh]*hru_area[hh];
      ssrinflow[hh] = soil_ssr_Buf[hh]*hru_area[hh];
      inflow[hh] = runoutflow[hh] + ssroutflow[hh] + gw_Amount; // add this HRU runoff and subsurface flow

      runcuminflow[hh] += runinflow[hh];
      runcumoutflow[hh] += runoutflow[hh];
      ssrcuminflow[hh] += ssrinflow[hh];
      ssrcumoutflow[hh] += ssroutflow[hh];
      cuminflow[hh] += inflow[hh];

      runinflow_mWQ_lay[Sub][hh] = soil_runoff_Buf_conc_lay[Sub][hh]*soil_runoff_Buf[hh]*hru_area[hh];
      runcuminflow_mWQ_lay[Sub][hh] += runinflow_mWQ_lay[Sub][hh];
      runcumoutflow_mWQ_lay[Sub][hh] += runoutflow_mWQ_lay[Sub][hh];

      ssrinflow_mWQ_lay[Sub][hh] = soil_ssr_Buf_conc_lay[Sub][hh]*soil_ssr_Buf[hh]*hru_area[hh];
      ssrcuminflow_mWQ_lay[Sub][hh] += ssrinflow_mWQ_lay[Sub][hh];
      ssrcumoutflow_mWQ_lay[Sub][hh] += ssroutflow_mWQ_lay[Sub][hh];

      inflow_mWQ_lay[Sub][hh] = runoutflow_mWQ_lay[Sub][hh] + ssroutflow_mWQ_lay[Sub][hh] + Amount_mWQ; // add this HRU runoff and subsurface flow

      for(long hhh = 0; chkStruct(hhh); ++hhh) {
          if(outflow[hhh] > 0.0){
            Amount = outflow[hhh]/hru_area[hh]; // outflow (mm*km^2/int)
            if(whereto[hhh] && whereto[hhh]-1 == hh){
              Amount = outflow[hhh]/hru_area[hh]; // outflow (mm*km^2/int)
              Amount_mWQ = outflow_mWQ_lay[Sub][hhh]/hru_area[hh];

            if(preferential_flow[hh]){
              gw_conc_lay[Sub][hh] = gw_conc_lay[Sub][hh]*gw[hh] + Amount_mWQ;
              gw[hh] += Amount;
              gw_conc_lay[Sub][hh] /= gw[hh];
              cum_preferential_flow_to_gw[hh] += Amount;
              Reset_WQ(hh, outflow, outflow_mWQ_lay);

              if(Sub == numsubstances-1){
                gw[hh] += Amount;
                cum_preferential_flow_to_gw[hh] += Amount;
                Amount = 0.0;
                Amount_mWQ = 0.0;
              }
            }
            else if(!soil_rechr_ByPass[hh] && Amount > 0.0){
                if(soil_rechr[hh] + Amount >= soil_rechr_max[hh]){ // units (mm*km^2/int)
                  float Excess = soil_rechr[hh] + Amount - soil_rechr_max[hh];
                  float Free = Amount - Excess;

                  conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + Amount_mWQ;
                  conc_top_lay[Sub][hh] /= (soil_rechr[hh] + Amount);

                  soil_moist_conc_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + conc_top_lay[Sub][hh]*soil_rechr[hh] + Amount_mWQ*Free)/(soil_lower[hh] + soil_rechr[hh] + Free); // present mQW

                  if(Sub == numsubstances-1){
                    soil_rechr[hh] += Free;
                    soil_moist[hh] = soil_lower[hh] + soil_rechr[hh];
                  }
                  Amount = Excess;
                  Amount_mWQ = conc_top_lay[Sub][hh]*Excess/Amount;
                }
                else{
                  conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + Amount_mWQ;
                  if(soil_rechr[hh] + Amount > 0.0)
                    conc_top_lay[Sub][hh] /= (soil_rechr[hh] + Amount);
                  else
                    conc_top_lay[Sub][hh] = 0.0;

                  soil_moist_conc_lay[Sub][hh] = (conc_bottom_lay[Sub][hh]*soil_lower[hh] + conc_top_lay[Sub][hh]*(soil_rechr[hh] + Amount_mWQ))/(soil_lower[hh] + soil_rechr[hh] + Amount); // amount used

                  if(Sub == numsubstances-1){
                    soil_rechr[hh] = soil_rechr[hh] + Amount;
                    soil_moist[hh] = soil_lower[hh] + soil_rechr[hh];
                  }

                  Amount = 0.0;
                  Amount_mWQ = 0.0;
                } // else
              } // if !soil_rechr_ByPass
              else if(!Sd_ByPass[hh] && Amount > 0.0){
                if(Sd[hh] + Amount >= Sdmax[hh]){ // units (mm*km^2/int)
                  float Excess = Sd[hh] + Amount - Sdmax[hh];
                  float Free = Amount - Excess;

                  Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + Amount_mWQ;
                  Sd_conc_lay[Sub][hh] /= (Sd[hh] + Amount);

                  cum_to_Sd_mWQ_lay[Sub][hh] += Amount_mWQ;

                  if(Sub == numsubstances-1){
                    Sd[hh] += Free;
                    cum_to_Sd[hh] += Amount;
                  }

                  Amount = Excess;
                  Amount_mWQ = Sd_conc_lay[Sub][hh]*Excess;
                }
                else{
                  Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + Amount_mWQ;

                  if(Sd[hh] + Amount > 0.0)
                    Sd_conc_lay[Sub][hh] /= (Sd[hh] + Amount);
                  else
                    Sd_conc_lay[Sub][hh] = 0.0;

                  if(Sub == numsubstances-1){
                    Sd[hh] = Sd[hh] + Amount;
                    cum_to_Sd[hh] += Amount;
                  }

                  Amount = 0.0;
                  Amount_mWQ = 0.0;
                } // else
              } // if !Sd_ByPass
              else if(Amount > 0.0){ // handle excess
                redirected_residual_conc_lay[Sub][hh] = Amount_mWQ*hru_area[hh]/(redirected_residual[hh] + Amount*hru_area[hh]);

                if(Sub == numsubstances-1)
                  redirected_residual[hh] += Amount*hru_area[hh];
              }
            } // whereto defined
          } // outflow[hhh] > 0.0
      } // for hhh

      if(gw_Amount > 0.0 && gwwhereto[hh] == 0){ // move to basin gw
        float basingw0 = basingw[0] + gw_Amount;
        basingw_conc_lay[Sub][0] = basingw_conc_lay[Sub][0]*basingw[0] + gw_Amount_mWQ*1000;

          cumbasingw_mWQ_lay[Sub][0] += basingw[0]*gw_Amount_mWQ;
          gwcumoutflow_mWQ_lay[Sub][hh] += gw_Amount_mWQ;

          if(basingw[0] > 0.0)
            basingw_conc_lay[Sub][0] /= basingw0;
          else
            basingw_conc_lay[Sub][0] = 0.0;

          if(Sub == numsubstances-1){
            basingw[0] = basingw0;
            cumbasingw[0] += gw_Amount;
            gwcumoutflow[hh] += gw_Amount;
          }

          gw_Amount = 0.0;
          gw_Amount_mWQ = 0.0;
      } // if gw_Amount > 0.0

      if(nstep == 0){ // end of every day
        if(soil_ssrDiv > 1) // daily value - ready for next day
           soil_ssr_Buf[hh] = soil_ssr[hh]/soil_ssrDiv;

        if(soil_runoffDiv > 1) // daily value - ready for next day
           soil_runoff_Buf[hh] = soil_runoff[hh]/soil_runoffDiv;

        if(soil_gwDiv > 1) // daily value - ready for next day
           soil_gw_Buf[hh] = soil_gw[hh]/soil_gwDiv;
      } // end if
      
/*      if(variation == VARIATION_ORG)
        hruDelay->DoMuskingum(hh); // need to update for later HRUs
      else if(variation == VARIATION_1)
        Clark_hruDelay->DoClark(hh); // need to update for later HRUs*/
    } // for jj accessing hh
  } // Sub

  if(variation == VARIATION_ORG){
    hruDelay->DoMuskingum(); // need to update for later HRUs
    for(long Sub = 0; Sub < numsubstances; Sub++)
      hruDelay_cWQ[Sub]->DoMuskingum(); // need to update for later HRUs
  }
  else if(variation == VARIATION_1){
    Clark_hruDelay->DoClark(); // need to update for later HRUs
    for(long Sub = 0; Sub < numsubstances; Sub++)
      Clark_hruDelay_cWQ[Sub]->DoClark(); // need to update for later HRUs
  }

  runDelay->DoClark();
  ssrDelay->DoClark();
  gwDelay->DoClark();

  for(long Sub = 0; Sub < numsubstances; Sub++){
    runDelay_cWQ[Sub]->DoClark();
    ssrDelay_cWQ[Sub]->DoClark();
    gwDelay_cWQ[Sub]->DoClark();
  }

  basinflow_s[0] = basinflow[0]*Global::Freq/86400.0;
  basingw_s[0] = basingw[0]*Global::Freq/86400.0;

  cumbasinflow[0] += basinflow[0];
  cumbasingw[0] += basingw[0];
}

void ClassWQ_Netroute::finish(bool good) {

  float Allcuminflow = 0.0;
  float Allcumoutflow = 0.0;
  float Allcumoutflowdiverted = 0.0;

  float Allcuminflow_mWQ = 0.0;
  float Allcumoutflow_mWQ = 0.0;
  float Allcumoutflowdiverted_mWQ = 0.0;

  float Allgwcuminflow = 0.0;
  float Allgwcumoutflow = 0.0;
  float Allgwcumoutflowdiverted = 0.0;

  float Allgwcuminflow_mWQ = 0.0;
  float Allgwcumoutflow_mWQ = 0.0;
  float Allgwcumoutflowdiverted_mWQ = 0.0;

  float Allssrcuminflow = 0.0;
  float Allssrcumoutflow = 0.0;
  float Allruncuminflow = 0.0;
  float Allruncumoutflow = 0.0;

  float Allssrcuminflow_mWQ = 0.0;
  float Allssrcumoutflow_mWQ = 0.0;
  float Allruncuminflow_mWQ = 0.0;
  float Allruncumoutflow_mWQ = 0.0;

  float AllSdcuminflow = 0.0;
  float Allrechrcuminflow = 0.0;

  float AllSdcuminflow_mWQ = 0.0;
  float Allrechrcuminflow_mWQ = 0.0;
  float AllTotal = 0.0;
  float Total;

  String S = String("H2o");
  LogDebug(S.c_str());
  LogMessage(" ");

  for(hh = 0; chkStruct(); ++hh) {
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cuminflow                   (mm) (mm*km^2) (mm*basin): ").c_str(), cuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cumoutflow                  (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cumoutflow_diverted         (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' hruDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), hruDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrcuminflow                (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrcumoutflow               (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), ssrDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runoffcuminflow             (mm) (mm*km^2) (mm*basin): ").c_str(), runcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runoffcumoutflow            (mm) (mm*km^2) (mm*basin): ").c_str(), runcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), runDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcuminflow                 (mm) (mm*km^2) (mm*basin): ").c_str(), gwcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcumoutflow                (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcumoutflow_diverted       (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwDelay_in_storage          (mm) (mm*km^2) (mm*basin): ").c_str(), gwDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_to_Sd                   (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_Sd[hh], hru_area[hh], basin_area[0], " *** Added to this HRU Sd");
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_to_soil_rechr           (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_soil_rechr[hh], hru_area[hh], basin_area[0], " *** Added to this HRU recharge");
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_redirected_residual     (mm) (mm*km^2) (mm*basin): ").c_str(), cum_redirected_residual[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' HRU_cumbasinflow            (mm) (mm*km^2) (mm*basin): ").c_str(), HRU_cumbasinflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");

//    Total = cumoutflow[hh] + gwcumoutflow[hh] - cumbasinflow[hh] - cum_to_Sd[hh] - cum_to_soil_rechr[hh] - gwcumoutflow[hh]
//            + cumoutflow_diverted[hh] + gwcumoutflow_diverted[hh] - cum_redirected_residual[hh];
//    AllTotal += Total;

//    LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' Total                       (mm) (mm*hru) (mm*hru/basin): ").c_str(), Total, hru_area[hh], basin_area[0], " *** HRU mass balance");
//    LogDebug(" ");
  } // hh

  for(Sub = 0; Sub < numsubstances; ++Sub){

    Allcuminflow = 0.0;
    Allcumoutflow = 0.0;
    Allcumoutflowdiverted = 0.0;

    Allcuminflow_mWQ = 0.0;
    Allcumoutflow_mWQ = 0.0;
    Allcumoutflowdiverted_mWQ = 0.0;

    Allgwcuminflow = 0.0;
    Allgwcumoutflow = 0.0;
    Allgwcumoutflowdiverted = 0.0;

    Allgwcuminflow_mWQ = 0.0;
    Allgwcumoutflow_mWQ = 0.0;
    Allgwcumoutflowdiverted_mWQ = 0.0;

    Allssrcuminflow = 0.0;
    Allssrcumoutflow = 0.0;
    Allruncuminflow = 0.0;
    Allruncumoutflow = 0.0;

    Allssrcuminflow_mWQ = 0.0;
    Allssrcumoutflow_mWQ = 0.0;
    Allruncuminflow_mWQ = 0.0;
    Allruncumoutflow_mWQ = 0.0;

    AllSdcuminflow = 0.0;
    Allrechrcuminflow = 0.0;

    AllSdcuminflow_mWQ = 0.0;
    Allrechrcuminflow_mWQ = 0.0;

    for(hh = 0; chkStruct(); ++hh) {
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cuminflow                   (mm) (mm*km^2) (mm*basin): ").c_str(), cuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cuminflow_mWQ               (mm) (mm*km^2) (mm*basin): ").c_str(), cuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cumoutflow                  (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cumoutflow_mWQ              (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cumoutflow_diverted         (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cumoutflow_diverted_mWQ     (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_diverted_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' hruDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), hruDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrcuminflow                (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrcuminflow_mWQ            (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrcumoutflow               (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrcumoutflow_mWQ           (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' ssrDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), ssrDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runoffcuminflow             (mm) (mm*km^2) (mm*basin): ").c_str(), runcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runoffcuminflow_mWQ         (mm) (mm*km^2) (mm*basin): ").c_str(), runcuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runoffcumoutflow            (mm) (mm*km^2) (mm*basin): ").c_str(), runcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runoffcumoutflow_mWQ        (mm) (mm*km^2) (mm*basin): ").c_str(), runcumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' runDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), runDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcuminflow                 (mm) (mm*km^2) (mm*basin): ").c_str(), gwcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcuminflow_mWQ             (mm) (mm*km^2) (mm*basin): ").c_str(), gwcuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcumoutflow                (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcumoutflow_mWQ            (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcumoutflow_diverted       (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwcumoutflow_diverted_mWQ   (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_diverted_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' gwDelay_in_storage          (mm) (mm*km^2) (mm*basin): ").c_str(), gwDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_to_Sd                   (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_Sd[hh], hru_area[hh], basin_area[0], " *** Added to this HRU Sd");
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_to_Sd_mWQ               (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_Sd_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0], " *** Added to this HRU Sd");
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_to_soil_rechr           (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_soil_rechr[hh], hru_area[hh], basin_area[0], " *** Added to this HRU recharge");
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_to_soil_rechr_mWQ       (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_soil_rechr_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0], " *** Added to this HRU recharge");
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_redirected_residual     (mm) (mm*km^2) (mm*basin): ").c_str(), cum_redirected_residual[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' cum_redirected_residual_mWQ (mm) (mm*km^2) (mm*basin): ").c_str(), cum_redirected_residual_mWQ[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' HRU_cumbasinflow            (mm) (mm*km^2) (mm*basin): ").c_str(), HRU_cumbasinflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' HRU_cumbasinflow_mwq        (mm) (mm*km^2) (mm*basin): ").c_str(), HRU_cumbasinflow_mWQ_lay[Sub][hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogDebug(" ");

//      Total = cumoutflow[hh] + gwcumoutflow[hh] - cumbasinflow[hh] - cum_to_Sd_mWQ_lay[Sub][hh] - cum_to_soil_rechr[hh] - gwcumoutflow[hh]
//              + cumoutflow_diverted[hh] + gwcumoutflow_diverted[hh] - cum_redirected_residual[hh];
//      AllTotal += Total;

//      LogMessageA(hh, string("'" + Name + " (Netroute_WQ)' Total                       (mm) (mm*hru) (mm*hru/basin): ").c_str(), Total, hru_area[hh], basin_area[0], " *** HRU mass balance");
//      LogDebug(" ");

      if(Sub == 0){ // sum only once
        Allcuminflow += cuminflow[hh];
        Allcumoutflow += cumoutflow[hh];

        Allgwcuminflow += gwcuminflow[hh];
        Allgwcumoutflow += gwcumoutflow[hh];

        Allssrcumoutflow += ssrcumoutflow[hh];
        Allssrcuminflow += ssrcuminflow[hh];
        Allruncuminflow += runcuminflow[hh];
        Allruncumoutflow += runcumoutflow[hh];

        AllSdcuminflow += cum_to_Sd[hh];
        Allrechrcuminflow += cum_to_soil_rechr[hh];

        Allcumoutflowdiverted += cumoutflow_diverted[hh];
      }

      Allcuminflow_mWQ += cuminflow_mWQ_lay[Sub][hh];
      Allcumoutflow_mWQ += cumoutflow_mWQ_lay[Sub][hh];
      Allcumoutflowdiverted_mWQ += cumoutflow_diverted_mWQ_lay[Sub][hh];

      Allgwcuminflow_mWQ += gwcuminflow_mWQ_lay[Sub][hh];
      Allgwcumoutflow_mWQ += gwcumoutflow_mWQ_lay[Sub][hh];
      Allgwcumoutflowdiverted_mWQ += gwcumoutflow_diverted_mWQ_lay[Sub][hh];

      Allssrcumoutflow_mWQ += ssrcumoutflow_mWQ_lay[Sub][hh];
      Allssrcuminflow_mWQ += ssrcuminflow_mWQ_lay[Sub][hh];
      Allruncuminflow_mWQ += runcuminflow_mWQ_lay[Sub][hh];
      Allruncumoutflow_mWQ += runcumoutflow_mWQ_lay[Sub][hh];

      AllSdcuminflow_mWQ += cum_to_Sd_mWQ_lay[Sub][hh];
      Allrechrcuminflow_mWQ += cum_to_soil_rechr_mWQ_lay[Sub][hh];
      LogMessage(" ");
    } // for hh

    LogMessage(string("'" + Name + " (Netroute_WQ)' cumbasinflow(m^3):           ").c_str(), cumbasinflow[0]);
    LogMessage(string("'" + Name + " (Netroute_WQ)' cumbasinflow_mWQ_lay(mg):    ").c_str(), cumbasinflow_mWQ_lay[Sub][0]);
    LogMessage(string("'" + Name + " (Netroute_WQ)' cumbasingw(m^3):             ").c_str(), cumbasingw[0]);
    LogMessage(string("'" + Name + " (Netroute_WQ)' cumbasingw_mWQ_lay(mg):      ").c_str(), cumbasingw_mWQ_lay[Sub][0]);
    LogMessage(" ");

    LogMessage(string("'" + Name + " (Netroute_WQ)' Allgwcuminflow (mm*basin):              ").c_str(), Allgwcuminflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allgwcuminflow_mWQ (mm*basin):          ").c_str(), Allgwcuminflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allgwcumoutflow (mm*basin):             ").c_str(), Allgwcumoutflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allgwcumoutflow_mWQ (mm*basin):         ").c_str(), Allgwcumoutflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allgwcumoutflowdiverted (mm*basin):     ").c_str(), Allgwcumoutflowdiverted);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allgwcumoutflowdiverted_mWQ (mm*basin): ").c_str(), Allgwcumoutflowdiverted_mWQ);
    LogDebug(" ");

    LogMessage(string("'" + Name + " (Netroute_WQ)' Allcuminflow (mm*basin):                ").c_str(), Allcuminflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allcuminflow_mWQ (mm*basin):            ").c_str(), Allcuminflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allcumoutflow (mm*basin):               ").c_str(), Allcumoutflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allcumoutflow_mWQ (mm*basin):           ").c_str(), Allcumoutflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allcumoutflowdiverted (mm*basin):       ").c_str(), Allcumoutflowdiverted);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allcumoutflowdiverted_mWQ (mm*basin):   ").c_str(), Allcumoutflowdiverted_mWQ);
    LogDebug(" ");

    LogMessage(string("'" + Name + " (Netroute_WQ)' Allssrcuminflow (mm*basin):             ").c_str(), Allssrcuminflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allssrcuminflow_mWQ (mm*basin):         ").c_str(), Allssrcuminflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allssrcumoutflow (mm*basin):            ").c_str(), Allssrcumoutflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allssrcumoutflow_mWQ (mm*basin):        ").c_str(), Allssrcumoutflow_mWQ);
    LogDebug(" ");

    LogMessage(string("'" + Name + " (Netroute_WQ)' Allruncuminflow (mm*basin):             ").c_str(), Allruncuminflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allruncuminflow_mWQ (mm*basin):         ").c_str(), Allruncuminflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allruncumoutflow (mm*basin):            ").c_str(), Allruncumoutflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allruncumoutflow_mWQ (mm*basin):        ").c_str(), Allruncumoutflow_mWQ);
    LogDebug(" ");

    LogMessage(string("'" + Name + " (Netroute_WQ)' AllSdcuminflow (mm*basin):              ").c_str(), AllSdcuminflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' AllSdcuminflow_mWQ (mm*basin):          ").c_str(), AllSdcuminflow_mWQ);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allrechrcuminflow (mm*basin):           ").c_str(), Allrechrcuminflow);
    LogMessage(string("'" + Name + " (Netroute_WQ)' Allrechrcuminflow_mWQ (mm*basin):       ").c_str(), Allrechrcuminflow_mWQ);
    LogDebug(" ");
    LogMessage(string("'" + Name + " (Netroute_WQ)' AllTotal              (mm*basin):       ").c_str(), AllTotal);
    LogDebug(" ");
  } // for Sub


  if(variation == VARIATION_ORG){
    delete hruDelay;
    delete hruDelay_cWQ;
  }
  else if(variation == VARIATION_1){
    delete Clark_hruDelay;
    delete Clark_hruDelay_cWQ;
  }

  delete ssrDelay;
  delete runDelay;
  delete gwDelay;
}

float ClassWQ_Netroute::Function1(float *I, long hh) {
  return runDelay->ChangeLag(I, hh);
}

float ClassWQ_Netroute::Function2(float *X, long hh) {
  return runDelay->ChangeStorage(X, hh);
}

void ClassWQ_Netroute::Reset_Basin_WQ(long hru, float *var, float *var_conc_lay){
  var[hru] = 0.0;
  var_conc_lay[hru] = 0.0;
}

void ClassWQ_Netroute::Reset_WQ(long hru, float *var, float **var_cWQ_lay){
  var[hru] = 0.0;
  for(long Sub = 0; Sub < numsubstances; Sub++){
    var_cWQ_lay[Sub][hru] = 0.0;
  }
}

void ClassWQ_Netroute::Set_WQ(const long hru, float *var, float *var_conc, float Amount, float amount_conc){

  var[hru] = Amount;
  if(Amount > 0.0)
    var_conc[hru] = amount_conc;
  else
    var_conc[hru] = 0.0;
}

ClassWQ_pbsm
* ClassWQ_pbsm::klone(string name) const{
  return new ClassWQ_pbsm(name);
}

void ClassWQ_pbsm::decl(void) {

  Description = "'calculates snow transport and sublimation (Pomeroy and Li, 1999).' \
                 'original version using hru_u,' \
                 'uses hru_Uadjust from walmsley_wind instead of hru_u,' \
                 'using hru_u and a regression to use daily windspeed,' \
                 'uses hru_Uadjust from walmsley_wind instead of hru_u and a regression to use daily windspeed.'";

  variation_set = VARIATION_0 + VARIATION_2;

  declgetvar("*", "hru_u", "(m/s)", &hru_u);


  variation_set = VARIATION_1 + VARIATION_3;

  declgetvar("*", "hru_Uadjust", "(m/s)", &hru_Uadjust);


  variation_set = VARIATION_2 + VARIATION_3;

  declparam("u_D", NHRU, "[1.0]", "0.5", "2.0", "Daily windspeed correction", "()", &u_D);

  declparam("Drift_offset", NHRU, "[0.0]", "0.0", "100.0", "Daily windspeed drift offset correction", "()", &Drift_offset);

  declparam("Drift_slope", NHRU, "[1.0]", "0.5", "2.0", "Daily windspeed drift slope correction", "()", &Drift_slope);

  declparam("Subl_offset", NHRU, "[0.0]", "0.0", "100.0", "Daily windspeed sublimation offset correction", "()", &Subl_offset);

  declparam("Subl_slope", NHRU, "[1.0]", "0.5", "2.0", "Daily windspeed sublimation slope correction", "()", &Subl_slope);


  variation_set = VARIATION_ORG;


  declstatvar("SWE", NHRU, "snow water equivalent", "(mm)", &SWE);

  declstatvar("SWE_max", NHRU, "snow water equivalent seasonal maximum", "(mm)", &SWE_max);

  declstatvar("SWE_conc", NDEFN, "snow water equivalent", "(mg/l)", &SWE_conc, &SWE_conc_lay, numsubstances);

  declvar("Subl", NHRU, "interval sublimation", "(mm/int)", &Subl);

  declvar("Subl_conc", NDEFN, "interval sublimation", "(mm/int)", &Subl_conc, &Subl_conc_lay, numsubstances);

  declvar("Drift_in", NHRU, "interval transport into HRU", "(mm/int)", &Drift_in);

  declvar("Drift_in_conc", NDEFN, "interval transport into HRU", "(mg/l)", &Drift_in_conc, &Drift_in_conc_lay, numsubstances);

  declvar("Drift_out", NHRU, "interval transport out of HRU", "(mm/int)", &Drift_out);

  declvar("Drift_out_conc", NDEFN, "interval transport out of HRU", "(mg/l)", &Drift_out_conc, &Drift_out_conc_lay, numsubstances);

  decldiag("DriftH", NHRU, "interval transport", "(mm/int)", &DriftH);

  decldiag("SublH", NHRU, "interval sublimation", "(mm/int)", &SublH);

  decldiag("BasinSnowLoss", BASIN, "transport out of basin", "(mm/int)", &BasinSnowLoss);

  decldiag("BasinSnowLoss_mWQ", NDEF, "transport out of basin", "(mm/int)", &BasinSnowLoss_mWQ, &BasinSnowLoss_mWQ_lay, numsubstances);

  decldiag("BasinSnowGain", BASIN, "cumulative transport to basin estimated from HRU 1", "(mm/int)", &BasinSnowGain);

  decldiag("BasinSnowGain_mWQ", NDEF, "cumulative transport to basin estimated from HRU 1", "(mm/int)", &BasinSnowGain_mWQ, &BasinSnowGain_mWQ_lay, numsubstances);

  declstatvar("cumSubl", NHRU, "cumulative sublimation", "(mm)", &cumSubl);

  declstatvar("cumSubl_mWQ", NDEFN, "cumulative sublimation solute", "(mm)", &cumSubl_mWQ, &cumSubl_mWQ_lay, numsubstances);

  declstatvar("cumDriftOut", NHRU, "cumulative transport from HRU", "(mm)", &cumDriftOut);

  declstatvar("cumDriftOut_mWQ", NDEFN, "mass solute from HRU", "(mg)", &cumDriftOut_mWQ, &cumDriftOut_mWQ_lay, numsubstances);

  declstatvar("cumBasinSnowLoss", BASIN, "cumulative transport out of basin", "(mm)", &cumBasinSnowLoss);

  declstatvar("cumBasinSnowLoss_mWQ", NDEF, "cumulative mass of solute transport out of basin", "(mg)", &cumBasinSnowLoss_mWQ, &cumBasinSnowLoss_mWQ_lay, numsubstances);

  declstatvar("cumBasinSnowGain", BASIN, "cumulative transport to basin estimated from HRU 1", "(mm)", &cumBasinSnowGain);

  declstatvar("cumBasinSnowGain_mWQ", NDEF, "cumulative mass of solute transport to basin estimated from HRU 1", "(mg)", &cumBasinSnowGain_mWQ, &cumBasinSnowGain_mWQ_lay, numsubstances);

  declstatvar("cumDriftIn", NHRU, "cumulative transport to HRU", "(mm)", &cumDriftIn);

  declstatvar("cumDriftIn_mWQ", NDEFN, "cumulative mass of solute transport to HRU", "(mg)", &cumDriftIn_mWQ, &cumDriftIn_mWQ_lay, numsubstances);

  decllocal("hru_basin", NHRU, "conversion factor", "()", &hru_basin);

  decldiag("DrySnow", NHRU, "DrySnow", "()", &DrySnow);

  declstatvar("SnowAge", NHRU, "SnowAge", "()", &SnowAge);

  declstatvar("cumSno", NHRU, "cumulative snow", "(mm)", &cumSno);

  declstatvar("cumSno_mWQ", NDEFN, "cumulative mass of solute snow", "(mg)", &cumSno_mWQ, &cumSno_mWQ_lay, numsubstances);

  declvar("Prob", NHRU, "Probability", "()", &Prob);

  declvar("snowdepth", NHRU, "depth of snow using Gray/Pomeroy", "(m)", &snowdepth);

  decllocal("SWE_Init", NHRU, "initial SWE", "(mm)", &SWE_Init);

// parameters

  declparam("fetch", NHRU, "[1000.0]", "300.0", "10000.0", "fetch distance", "(m)", &fetch);

  declparam("Ht", NHRU, "[0.1, 0.25, 1.0]", "0.001", "100.0", "vegetation height(m)", "(m)", &Ht);

  declparam("distrib", NHRU, "[0.0, 1.0]", "-10.0", "10.0", "distribution fractions - can sum to 1", "()", &distrib);

  declparam("N_S", NHRU, "[320]", "1", "500", "vegetation number density", "(1/m^2)", &N_S);

  declparam("A_S", NHRU, "[0.003]", "0.0", "2.0", "stalk diameter", "(m)", &A_S);

  declparam("basin_area", BASIN, "3", "1e-6", "1e+09", "total basin area", "(km^2)", &basin_area);

  declparam("hru_area", NHRU, "[1]", "1e-6", "1e+09", "hru area", "(km^2)", &hru_area);

  declparam("inhibit_evap", NHRU, "[0]", "0", "1", "inhibit evaporatation, 1 -> inhibit", "()", &inhibit_evap);

  declparam("inhibit_bs", NHRU, "[0]", "0", "1", "inhibit blowing snow, 1 -> inhibit", "()", &inhibit_bs);

  declparam("inhibit_subl", NHRU, "[0]", "0", "1", "inhibit sublimation, 1 -> inhibit", "()", &inhibit_subl);

  declparam("rain_conc", NDEFN, "0", "0", "1000", "rain solute concentration", "(mg/l)", &rain_conc, &rain_conc_lay, numsubstances);

  declparam("snow_conc", NDEFN, "0", "0", "1000", "snow solute concentration", "(mg/l)", &snow_conc, &snow_conc_lay, numsubstances);

  declparam("Atmos_mWQ", NDEFN, "0", "0", "10", "total basin area", "(mg/int)", &Atmos_mWQ, &Atmos_mWQ_lay, numsubstances);


  decllocal("DrySnow_0", NHRU, "", "", &DrySnow_0);

  decllocal("SnowAge_0", NHRU, "", "", &SnowAge_0);

  decllocal("BasinSnowGain_0", NHRU, "", "", &BasinSnowGain_0);

  decllocal("cumBasinSnowGain_0", NHRU, "", "", &cumBasinSnowGain_0);

  decllocal("BasinSnowLoss_0", NHRU, "", "", &BasinSnowLoss_0);

  decllocal("cumBasinSnowLoss_0", NHRU, "", "", &cumBasinSnowLoss_0);

  decllocal("Subl_0", NHRU, "", "", &Subl_0);

  decllocal("SublH_0", NHRU, "", "", &SublH_0);

  decllocal("cumSubl_0", NHRU, "", "", &cumSubl_0);

  decllocal("Drift_in_0", NHRU, "", "", &Drift_in_0);

  decllocal("cumDriftIn_0", NHRU, "", "", &cumDriftIn_0);

  decllocal("Drift_out_0", NHRU, "", "", &Drift_out_0);

  decllocal("cumDriftOut_0", NHRU, "", "", &cumDriftOut_0);

  decllocal("SWE_0", NHRU, "", "", &SWE_0);

  decllocal("SWE_Init_0", NHRU, "", "", &SWE_Init_0);

  decllocal("cumSno_0", NHRU, "", "", &cumSno_0);

  decllocal("DriftH_0", NHRU, "", "", &DriftH_0);

  decllocal("SublH_0", NHRU, "", "", &SublH_0);

  decllocal("Prob_0", NHRU, "", "", &Prob_0);

  decllocal("rho_0", NHRU, "", "", &rho_0);

  decllocal("z_s_0", NHRU, "", "", &z_s_0);

  declgetvar("*", "hru_t", "(�C)", &hru_t);
  declgetvar("*", "hru_rh", "(%)", &hru_rh);
  declgetvar("*", "hru_newsnow", "()", &hru_newsnow);
  declgetvar("*", "net_snow", "(mm/int)", &net_snow);

}

void ClassWQ_pbsm::init(void) {

  nhru = getdim(NHRU);

  cumBasinSnowLoss[0] = 0.0;
  cumBasinSnowGain[0] = 0.0;
  BasinSnowLoss_mWQ[0] = 0.0;
  BasinSnowGain_mWQ[0] = 0.0;
  cumBasinSnowLoss_mWQ[0] = 0.0;
  cumBasinSnowGain_mWQ[0] = 0.0;


  for (hh = 0; hh < nhru; ++hh) {
    for(long Sub = 0; Sub < numsubstances; ++Sub){
      Reset_WQ(hh, SWE, SWE_conc_lay);
      Reset_WQ(hh, Drift_in, Drift_in_conc_lay);
      Reset_WQ(hh, Drift_out, Drift_out_conc_lay);
      Reset_WQ(hh, cumDriftOut, cumDriftOut_mWQ_lay);
      Reset_WQ(hh, cumDriftIn, cumDriftIn_mWQ_lay);
      Reset_WQ(hh, cumSno, cumSno_mWQ_lay);
      Reset_WQ(hh, cumSubl, cumSubl_mWQ_lay);
    }

    BasinSnowLoss[hh] = 0.0;
    BasinSnowGain[hh] = 0.0;
    cumBasinSnowLoss[hh] = 0.0;
    cumBasinSnowGain[hh] = 0.0;

    SnowAge[hh] = 0.0;
    DrySnow[hh] = 0;
    snowdepth[hh] = 0.0;
    Prob[hh] = 0.0;

    if((hh > 0) && (Ht[hh] < Ht[hh-1]) && distrib[hh-1] > 0){
      CRHMException TExcept(string("'" + Name + " (pbsm)' vegetation heights not in ascending order.").c_str(), WARNING);
      LogError(TExcept);
    }
  }

  for (hh = 0; hh < nhru; ++hh)
    hru_basin[hh] = hru_area[hh]/basin_area[0];
}


void ClassWQ_pbsm::run(void){

  float Znod, Ustar, Ustn, E_StubHt, Lambda, Ut, Uten_Prob;
  float SumDrift, SumDrift_conc, total, transport;
  long Sub = 0;

  for(long Sub = 0; Sub < numsubstances; ++Sub){
    if(getstep() == 1)
      for (hh = 0; chkStruct(); ++hh)
        SWE_Init[hh] = SWE[hh];

    if(Sub == 0) // saves all HRUs
      Save();

    for (hh = 0; chkStruct(); ++hh) {

      if(Sub != 0)
        Restore(hh);

      if(net_snow[hh] > 0.0) {
        SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + Atmos_mWQ_lay[Sub][hh]*net_snow[hh])/(SWE[hh] + net_snow[hh]);
        SWE[hh] = SWE[hh] + net_snow[hh];
        cumSno[hh] = cumSno[hh] + net_snow[hh];
        cumSno_mWQ_lay[Sub][hh] += net_snow[hh]*Atmos_mWQ_lay[Sub][hh];
      }

     if(variation == VARIATION_ORG || variation == VARIATION_2)
       hru_u_ = hru_u[hh];
     else
       hru_u_ = hru_Uadjust[hh];

     if(variation == VARIATION_2 || variation == VARIATION_3)
       hru_u_ = u_D[hh]*hru_u_;

     Reset_WQ(hh, Drift_in, Drift_in_conc_lay);
     Reset_WQ(hh, Drift_out, Drift_out_conc_lay);
     Reset_WQ(hh, Subl, Subl_conc_lay);

     DriftH[hh] = 0.0;
     SublH[hh] = 0.0;
     Prob[hh] = 0.0;

     if(SWE[hh] > 0.0)
       SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + Atmos_mWQ_lay[Sub][hh])/SWE[hh];
     else
       SWE_conc_lay[Sub][hh] = 0.0;
  ;

     if(SWE[hh] > 0.0 && !inhibit_bs[hh]) {

       E_StubHt = Ht[hh] - DepthofSnow(SWE[hh]); // depths(m) SWE(mm)
       if(E_StubHt < 0.0001) E_StubHt = 0.0001;

       Ustar = 0.02264*pow(hru_u_, 1.295f); // Eq. 6.2 rev.,  Ustar over fallow

       if (E_StubHt > 0.01) {
         Znod = (sqr(Ustar)/163.3f)+0.5*N_S[hh]*E_StubHt*A_S[hh]; // Eq. 29, Snowcover Book
         Lambda = N_S[hh]*A_S[hh]*E_StubHt;  // Raupach Eq. 1

         Ustn  = Ustar*sqrt((PBSM_constants::Beta*Lambda)/(1.0+PBSM_constants::Beta*Lambda));

         Uten_Prob = (log(10.0/Znod))/PBSM_constants::KARMAN *min <float> (0.0, Ustar-Ustn);
       }
       else
       {
         Uten_Prob = hru_u_;
       } // end if


       ProbabilityThresholdNew(SWE[hh], hru_t[hh], Uten_Prob, Prob[hh], Ut, hru_newsnow[hh], SnowAge[hh], DrySnow[hh]);

       if (Prob[hh] > 0.001) {
         Ut = Ut * 0.8;

         Pbsm(E_StubHt, Ut, DriftH[hh], SublH[hh], hru_t[hh], hru_u_, hru_rh[hh], fetch[hh], N_S[hh], A_S[hh]);

         if(variation == VARIATION_2 || variation == VARIATION_3){
           DriftH[hh] = Drift_offset[hh] + DriftH[hh]*Drift_slope[hh];
           SublH[hh] = Subl_offset[hh] + SublH[hh]*Subl_slope[hh];
         }

         Drift_out[hh] = DriftH[hh]*Prob[hh]/fetch[hh];
         if(!inhibit_subl[hh])
           Subl[hh] = SublH[hh]*Prob[hh];

  // handle insufficient snow pack

         if(Drift_out[hh] + Subl[hh] > SWE[hh]) {
           Subl[hh] = SWE[hh] * Subl[hh]/(Subl[hh] + Drift_out[hh]);
           Drift_out[hh] = SWE[hh] - Subl[hh];
  //?         SWE[hh] = 0.0;
         } // end if
         else
           SWE[hh] = SWE[hh] - Subl[hh] - Drift_out[hh];

         Drift_out_conc_lay[Sub][hh] = SWE_conc_lay[Sub][hh];
         Subl_conc_lay[Sub][hh] = SWE_conc_lay[Sub][hh];


         cumDriftOut[hh] = cumDriftOut[hh] + Drift_out[hh];
         cumDriftOut_mWQ_lay[Sub][hh] += Drift_out_conc_lay[Sub][hh]*Drift_out[hh];
         cumSubl[hh] = cumSubl[hh] + Subl[hh];
         cumSubl_mWQ_lay[Sub][hh] += Subl_conc_lay[Sub][hh]*Subl[hh];
       }
     } // end if
   } // end for (hh)

   // distribute drift

    if(distrib[0] > 0.0 && Drift_out[0] > 0.0) { // simulate transport entering basin using HRU 1
      float Drft = Drift_out[0]*distrib[0];
      SWE_conc_lay[Sub][0] = (SWE_conc_lay[Sub][0]*SWE[0] + SWE_conc_lay[Sub][0]*Drft)/(SWE[0] + Drft);
      SWE[0] += Drft;
      cumDriftIn[0] += Drft;
      BasinSnowGain[0] = Drft*hru_basin[0]; // **** hru_basin = hru_area/basin_area ****
      BasinSnowGain_mWQ_lay[Sub][0] = SWE_conc_lay[Sub][hh];
      cumBasinSnowGain[0] += BasinSnowGain[0];
      cumDriftIn_mWQ_lay[Sub][0] += Drft*SWE_conc_lay[Sub][0];
      cumBasinSnowGain_mWQ_lay[Sub][0] += BasinSnowGain[0]*SWE_conc_lay[Sub][0];
    }

    BasinSnowLoss[0] = 0.0;
    long LastN = 0;

    if(!inhibit_bs[0]&& nhru == 1){
      BasinSnowLoss[0] = Drift_out[0];
      BasinSnowLoss_mWQ_lay[Sub][0] = BasinSnowLoss[0]*SWE_conc_lay[Sub][0];
      cumBasinSnowLoss[0] += BasinSnowLoss[0];
      cumBasinSnowLoss_mWQ_lay[Sub][0] += BasinSnowLoss_mWQ_lay[Sub][0]; // Not correct - source of snow?
    }

    for (long nn = LastN; chkStruct(nn); ++nn) {
      if(distrib[nn] >= 0.0 && nn+1 < nhru) // skip till last HRU or -ve distribution
        continue;

      SumDrift = 0.0; SumDrift_conc = 0.0;
      for (long hhh = LastN; chkStruct(hhh, nn); ++hhh){ // sum drift over range
  //      if(Drift_out[hhh] > 0.0){
          SumDrift_conc = SumDrift*SumDrift_conc + Drift_out[hhh]*hru_basin[hhh]*Drift_out_conc_lay[Sub][hhh];
          SumDrift += Drift_out[hhh]*hru_basin[hhh];

          if(SumDrift > 0.0)
            SumDrift_conc /= SumDrift;
          else
            SumDrift_conc = 0.0;
  //      }
      }
      if(SumDrift > 0.0){ // drift has occurred!
        for (long hh = LastN + 1; chkStruct(hh, nn+1); ++hh) {
          SWE_max[hh] = SWEfromDepth(Ht[hh]);

          if(hh == nn) { // handle last HRU
            if(distrib[nn] > 0){
              float In = SumDrift/hru_basin[hh]; // remaining drift
              Drift_in_conc_lay[Sub][hh] = SumDrift_conc;
              if(SWE_max[hh] > SWE[hh] + In){ // fill snowpack, remainder leaves basin
                Drift_in[hh] = In;
                SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + SumDrift_conc*In)/(SWE[hh] + In);
                SWE[hh]  += In; // can handle all
                cumDriftIn[hh] += In;
                cumDriftIn_mWQ_lay[Sub][hh] += SumDrift_conc*In;
                transport = 0.0;
              }
              else if(SWE_max[hh] > SWE[hh]){ // cannot handle all
                float used = SWE_max[hh] - SWE[hh];
                Drift_in[hh] = used;
                Drift_in_conc_lay[Sub][hh] = SumDrift_conc;
                cumDriftIn[hh] += Drift_in[hh];
                transport -= (In - used)*hru_basin[hh];
                SWE[hh]  += used; //  has to come last
                SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + SumDrift_conc*used)/(SWE[hh] + used);
                cumDriftIn_mWQ_lay[Sub][hh] += SumDrift_conc*used;
              }
              else // zero or -ve - happens during melt??
                transport = SumDrift;
            }
            else if(distrib[nn] < 0){ // all drift deposited
                float used = SumDrift/hru_basin[hh];
                SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + SumDrift_conc*used)/(SWE[hh] + used);
                Drift_in[hh] = used;
                Drift_in_conc_lay[Sub][hh] = SumDrift_conc;
                SWE[hh]  += SumDrift/hru_basin[hh]; // can handle all
                cumDriftIn[hh] += Drift_in[hh];
                cumDriftIn_mWQ_lay[Sub][hh] += SumDrift/hru_basin[hh]*SumDrift_conc;
                transport = 0.0;
            }
            else // distrib[nn] == 0 -> all excess drift leaves basin
                transport = SumDrift;

            BasinSnowLoss[0] += (transport + Drift_out[hh]*hru_basin[hh]);
            BasinSnowLoss_mWQ_lay[Sub][0] = SumDrift_conc;
            cumBasinSnowLoss[0] += (transport + Drift_out[hh]*hru_basin[hh]);
            cumBasinSnowLoss_mWQ_lay[Sub][0] += (transport*SumDrift_conc + Drift_out[hh]*Drift_out_conc_lay[Sub][hh]*hru_basin[hh]); // check
          }
          else if(SWE_max[hh] > SWE[hh] &&  distrib[hh] > 0.0) {
  // handle intermediate HRUs with available storage and distrib > 0
            total = 0.0;
            for (long jj = hh; chkStruct(jj, nn+1); jj++) // calculate denominator        !!!! nn+1
              total += fabs(distrib[jj]);
  // determine contribution and scale
            transport = SumDrift*fabs(distrib[hh])/total/hru_basin[hh];
            if(SWE_max[hh] > SWE[hh] + transport){ // sufficient capacity
              SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + SumDrift_conc*transport)/(SWE[hh] + transport);
              Drift_in[hh] = transport;
              Drift_in_conc_lay[Sub][hh] = SumDrift_conc;
              cumDriftIn[hh] += Drift_in[hh];
              SWE[hh] += transport;
            }
            else{
              transport = SWE_max[hh] - SWE[hh];  // insufficient capacity
              SWE_conc_lay[Sub][hh] = (SWE_conc_lay[Sub][hh]*SWE[hh] + SumDrift_conc*transport)/(SWE[hh] + transport);
              Drift_in[hh] = transport;
              Drift_in_conc_lay[Sub][hh] = SumDrift_conc;
              cumDriftIn[hh] += Drift_in[hh];
              SWE[hh] = SWE_max[hh];
            }
            SumDrift_conc = SumDrift*SumDrift_conc - transport*hru_basin[hh]*SumDrift_conc;
            SumDrift -= transport*hru_basin[hh]; // remove drift used from total available
            cumDriftIn[hh] += transport;
            cumDriftIn_mWQ_lay[Sub][hh] += transport*SumDrift_conc;
            SumDrift_conc = SumDrift_conc/SumDrift;
          } // end if
        } // end for (hh)
      } // end if
      LastN = nn+1;
    } // end for (nn)

    for (hh = 0; chkStruct(); ++hh) { // snow cover inhibits evaporation

      if(SWE[hh] > 0.0){
        const_cast<long*> (inhibit_evap)[hh] = 1;
        snowdepth[hh] = DepthofSnow(SWE[hh]);
      }
      else{
        const_cast<long*> (inhibit_evap)[hh] = 0;
        snowdepth[hh] = 0.0;
        SWE_conc_lay[Sub][hh] = 0.0;
      }
    } // end for (hh)
  } // end for (Sub)
}

void ClassWQ_pbsm::finish(bool good) {

  if(!good) return;

  float AllcumSubl = 0.0;
  float AllcumCover = cumBasinSnowGain[0] - cumBasinSnowLoss[0];
  long Sub = 0;

  for(hh = 0; chkStruct(); ++hh) {
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumSno     (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSno[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumDriftOut(mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftOut[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumDriftIn (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftIn[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumSubl    (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSubl[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumCover   (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSno[hh]+cumDriftIn[hh]-cumDriftOut[hh]-cumSubl[hh], hru_area[hh], basin_area[0], "*** SWE just before melt");
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' SWE        (mm) (mm*hru) (mm*hru/basin): ").c_str(), SWE[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' SWE_change (mm) (mm*hru) (mm*hru/basin): ").c_str(), SWE[hh] - SWE_Init[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");

    AllcumSubl += cumSubl[hh]*hru_area[hh];
    AllcumCover += (cumSno[hh]+cumDriftIn[hh]-cumDriftOut[hh]-cumSubl[hh])*hru_area[hh];
    LogDebug(" ");

    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumSno_mWQ     (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSno_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumDriftOut_mWQ(mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftOut_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsm)' cumDriftIn_mWQ (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftIn_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");
  }

  LogMessage(string("'" + Name + " (WQ_pbsm)' AllcumSubl  (mm*basin): ").c_str(), AllcumSubl, "*** cumulative sum of all HRUs cumSubl");
  LogMessage(string("'" + Name + " (WQ_pbsm)' AllcumCover (mm*basin): ").c_str(), AllcumCover, "*** SWE just before melt cumulative sum of all HRUs cumCover");

  LogDebug(" ");
  LogMessage("'WQ_pbsm' cumBasinSnowLoss (mm): ", cumBasinSnowLoss[0]);
  LogMessage("'WQ_pbsm' cumBasinSnowGain (mm): ", cumBasinSnowGain[0]);

  LogDebug(" ");
  LogMessage("'WQ_pbsm' cumBasinSnowLoss_mWQ (substance) (mm): ", cumBasinSnowLoss_mWQ_lay[Sub][0]);
  LogMessage("'WQ_pbsm' cumBasinSnowGain_mWQ (substance) (mm): ", cumBasinSnowGain_mWQ_lay[Sub][0]);
  LogDebug(" ");

}

void ClassWQ_pbsm::Reset_WQ(long hru, float *var, float **var_WQ_lay){
  var[hru] = 0.0;
  for(long Sub = 0; Sub < numsubstances; Sub++)
    var_WQ_lay[Sub][hru] = 0.0;
}

void ClassWQ_pbsm::copy_array(float *from, float *to){
  for(hh = 0; chkStruct(); ++hh)
    to[hh] = from[hh];
}

void ClassWQ_pbsm::copy_array(long *from, long *to){
  for(hh = 0; chkStruct(); ++hh)
    to[hh] = from[hh];
}

void ClassWQ_pbsm::copy_basin(float *from, float *to){
  to[0] = from[0];
}

void ClassWQ_pbsm::restore_hru(float *from, float *to, long hh){
  to[hh] = from[hh];
}

void ClassWQ_pbsm::restore_hru(long *from, long *to, long hh){
  to[hh] = from[hh];
}

void ClassWQ_pbsm::Restore(long hh){

  restore_hru(DrySnow_0, DrySnow, hh);
  restore_hru(SnowAge_0, SnowAge, hh);

  restore_hru(BasinSnowGain, BasinSnowGain_0, hh);
  restore_hru(cumBasinSnowGain, cumBasinSnowGain_0, hh);
  restore_hru(BasinSnowLoss, BasinSnowLoss_0, hh);
  restore_hru(cumBasinSnowLoss, cumBasinSnowLoss_0, hh);

  restore_hru(SWE_0, SWE, hh);
  restore_hru(SWE_Init_0, SWE_Init, hh);
  restore_hru(cumSno_0, cumSno, hh);
  restore_hru(DriftH_0, DriftH, hh);
  restore_hru(SublH_0, SublH, hh);
  restore_hru(Prob_0, Prob, hh);
  restore_hru(Drift_in_0, Drift_in, hh);
  restore_hru(cumDriftIn_0, cumDriftIn, hh);
  restore_hru(Drift_out_0, Drift_out, hh);
  restore_hru(cumDriftOut_0, cumDriftOut, hh);
  restore_hru(Subl_0, Subl, hh);
  restore_hru(cumSubl_0, cumSubl, hh);
}

void ClassWQ_pbsm::Save(){

  copy_array(DrySnow, DrySnow_0);
  copy_array(SnowAge, SnowAge_0);

  copy_basin(BasinSnowGain, BasinSnowGain_0);
  copy_basin(cumBasinSnowGain, cumBasinSnowGain_0);
  copy_basin(BasinSnowLoss, BasinSnowLoss_0);
  copy_basin(cumBasinSnowLoss, cumBasinSnowLoss_0);

  copy_array(SWE, SWE_0);
  copy_array(SWE_Init, SWE_Init_0);
  copy_array(cumSno, cumSno_0);
  copy_array(DriftH, DriftH_0);
  copy_array(SublH, SublH_0);
  copy_array(Prob, Prob_0);
  copy_array(Drift_in, Drift_in_0);
  copy_array(cumDriftIn, cumDriftIn_0);
  copy_array(Drift_out, Drift_out_0);
  copy_array(cumDriftOut, cumDriftOut_0);
  copy_array(Subl, Subl_0);
  copy_array(cumSubl, cumSubl_0);
}

void ClassWQ_Netroute_M_D::decl(void) {

// kg/km2 = (mg/l)*mm = mg/m2

  Description = "'Handles the routing of surface runoff, subsurface runoff and HRU routing using the lag and route method.'\
                    'uses Muskingum,' \
                    'uses Clark.'";

  declvar("inflow", NHRU, "inflow from other HRUs", "(mm*km^2/int)", &inflow);

  declvar("inflow_mWQ", NDEFN, "Concentration: inflow from other HRUs", "(g/l)", &inflow_mWQ, &inflow_mWQ_lay, numsubstances);

  declstatvar("cuminflow", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &cuminflow);

  declstatvar("cuminflow_mWQ", NDEFN, "cumulative mass of solute inflow from other HRUs", "(g*km^2)", &cuminflow_mWQ, &cuminflow_mWQ_lay, numsubstances);

  declvar("outflow", NHRU, "HRU outflow", "(mm*km^2/int)", &outflow);

  declvar("outflow_mWQ", NDEFN, "Concentration: HRU outflow", "(g/l)", &outflow_mWQ, &outflow_mWQ_lay, numsubstances);

  declstatvar("cumoutflow", NHRU, "cumulative HRU outflow", "(mm*km^2)", &cumoutflow);

  declstatvar("cumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU outflow", "(g*km^2)", &cumoutflow_mWQ, &cumoutflow_mWQ_lay, numsubstances);

  decldiag("outflow_diverted", NHRU, "HRU outflow diverted to another HRU", "(mm*km^2/int)", &outflow_diverted);

  decldiag("outflow_diverted_conc", NDEFN, "Concentration: HRU outflow diverted to another HRU", "(g/l)", &outflow_diverted_conc, &outflow_diverted_conc_lay, numsubstances);

  declstatvar("cumoutflow_diverted", NHRU, "cumulative HRU outflow diverted to another HRU", "(mm*km^2/int)", &cumoutflow_diverted);

  declstatvar("cumoutflow_diverted_mWQ", NDEFN, "cumulative mass of solute HRU outflow diverted to another HRU", "(g*km^2/int)", &cumoutflow_diverted_mWQ, &cumoutflow_diverted_mWQ_lay, numsubstances);

  declstatvar("cum_to_Sd", NHRU, "cumulative other HRU to depressional storage (Sd) of this HRU", "(mm)", &cum_to_Sd);

  declstatvar("cum_to_Sd_mWQ", NDEFN, "cumulative mass of solute from other HRU to depressional storage (Sd) of this HRU", "(mg)", &cum_to_Sd_mWQ, &cum_to_Sd_mWQ_lay, numsubstances);

  declstatvar("cum_to_soil_rechr", NHRU, "cumulative other HRU to soil_rechr of this HRU", "(mm)", &cum_to_soil_rechr);

  declstatvar("cum_to_soil_rechr_mWQ", NDEFN, "cumulative mass of solute from other HRU to soil_rechr of this HRU", "(mg)", &cum_to_soil_rechr_mWQ, &cum_to_soil_rechr_mWQ_lay, numsubstances);

  declvar("gwinflow", NHRU, "ground water inflow", "(mm*km^2/int)", &gwinflow);

  declvar("gwinflow_mWQ", NDEFN, "Concentration: ground water inflow", "(g/l)", &gwinflow_mWQ, &gwinflow_mWQ_lay, numsubstances);

  declstatvar("gwcuminflow", NHRU, "cumulative gw inflow", "(mm*km^2)", &gwcuminflow);

  declstatvar("gwcuminflow_mWQ", NDEFN, "cumulative mass of solute gw inflow", "(g*km^2)", &gwcuminflow_mWQ, &gwcuminflow_mWQ_lay, numsubstances);

  declvar("gwoutflow", NHRU, "HRU gw outflow", "(mm*km^2/int)", &gwoutflow);

  declvar("gwoutflow_mWQ", NDEFN, "Concentration: HRU gw outflow", "(g/l)", &gwoutflow_mWQ, &gwoutflow_mWQ_lay, numsubstances);

  declstatvar("gwcumoutflow", NHRU, "cumulative HRU gw outflow", "(mm*km^2)", &gwcumoutflow);

  declstatvar("gwcumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU gw outflow", "(g*km^2)", &gwcumoutflow_mWQ, &gwcumoutflow_mWQ_lay, numsubstances);

  decldiag("gwoutflow_diverted", NHRU, "HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwoutflow_diverted);

  decldiag("gwoutflow_diverted_conc", NDEFN, "HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwoutflow_diverted_conc, &gwoutflow_diverted_conc_lay, numsubstances);

  declstatvar("gwcumoutflow_diverted", NHRU, "cumulative HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwcumoutflow_diverted);

  declstatvar("gwcumoutflow_diverted_mWQ", NDEFN, "cumulative mass of solute HRU gw outflow diverted to another HRU", "(mm*km^2/int)", &gwcumoutflow_diverted_mWQ, &gwcumoutflow_diverted_mWQ_lay, numsubstances);

  declvar("ssrinflow", NHRU, "inflow from other HRUs", "(mm*km^2/int)", &ssrinflow);

  declvar("ssrinflow_mWQ", NDEFN, "Concentration: inflow from other HRUs", "(g/l)", &ssrinflow_mWQ, &ssrinflow_mWQ_lay, numsubstances);

  declstatvar("ssrcuminflow", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &ssrcuminflow);

  declstatvar("ssrcuminflow_mWQ", NDEFN, "cumulative mass of solute of inflow from other HRUs", "(g*km^2)", &ssrcuminflow_mWQ, &ssrcuminflow_mWQ_lay, numsubstances);

  declvar("ssroutflow", NHRU, "HRU outflow", "(mm*km^2/int)", &ssroutflow);

  declvar("ssroutflow_mWQ", NDEFN, "Concentration: HRU outflow", "(g/l)", &ssroutflow_mWQ, &ssroutflow_mWQ_lay, numsubstances);

  declstatvar("ssrcumoutflow", NHRU, "cumulative HRU outflow", "(mm*km^2)", &ssrcumoutflow);

  declstatvar("ssrcumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU outflow", "(g*km^2)", &ssrcumoutflow_mWQ, &ssrcumoutflow_mWQ_lay, numsubstances);

  declstatvar("HRU_cumbasinflow", NHRU, "cumulative HRU to basinflow", "(mm*km^2)", &HRU_cumbasinflow);

  declstatvar("HRU_cumbasinflow_mWQ", NDEFN, "cumulative HRU to basinflow", "(mm*km^2)", &HRU_cumbasinflow_mWQ, &HRU_cumbasinflow_mWQ_lay, numsubstances);

  declvar("runinflow", NHRU, "inflow from other HRUs", "(mm*km^2/int)", &runinflow);

  declvar("runinflow_mWQ", NDEFN, "Concentration: inflow from other HRUs", "(g/l)", &runinflow_mWQ, &runinflow_mWQ_lay, numsubstances);

  declstatvar("runcuminflow", NHRU, "cumulative inflow from other HRUs", "(mm*km^2)", &runcuminflow);

  declstatvar("runcuminflow_mWQ", NDEFN, "cumulative mass of solute inflow from other HRUs", "(g*km^2)", &runcuminflow_mWQ, &runcuminflow_mWQ_lay, numsubstances);

  declvar("runoutflow", NHRU, "HRU outflow", "(mm*km^2/int)", &runoutflow);

  declvar("runoutflow_mWQ", NDEFN, "Concentration: HRU outflow", "(g/l)", &runoutflow_mWQ, &runoutflow_mWQ_lay, numsubstances);

  declstatvar("runcumoutflow", NHRU, "cumulative HRU outflow", "(mm*km^2)", &runcumoutflow);

  declstatvar("runcumoutflow_mWQ", NDEFN, "cumulative mass of solute HRU outflow", "(g*km^2)", &runcumoutflow_mWQ, &runcumoutflow_mWQ_lay, numsubstances);

  declstatvar("cum_preferential_flow_to_gw", NHRU, "cumulative other HRU's runoff to gw of this HRU via preferential flow path", "(mm)", &cum_preferential_flow_to_gw);


  declvar("basinflow", BASIN, "basin surface and sub-surface outflow", "(m^3/int)", &basinflow);

  declvar("basinflow_conc", NDEF, "basin surface and sub-surface outflow", "(g/l)", &basinflow_conc, &basinflow_conc_lay, numsubstances);

  declvar("Used", NHRU, "directed to basinbasin surface and sub-surface outflow", "()", &Used);

  declvar("Used_mWQ", NDEFN, "directed to basinbasin surface and sub-surface outflow", "()", &Used_mWQ, &Used_mWQ_lay, numsubstances);

  decldiag("basinflow_s", BASIN, "basin surface and sub-surface outflow", "(m^3/s)", &basinflow_s);

  declstatvar("cumbasinflow", BASIN, "cumulative basin surface and sub-surface outflow", "(m^3)", &cumbasinflow);

  declvar("cumbasinflow_mWQ", NDEF, "cumulative mass of solute basin surface and sub-surface outflow", "(m^3)", &cumbasinflow_mWQ, &cumbasinflow_mWQ_lay, numsubstances);

  declvar("basingw", BASIN, "cumulative basin groundwater outflow", "(m^3/int)", &basingw);

  declvar("basingw_conc", NDEF, "cumulative basin groundwater outflow", "(m^3/int)", &basingw_conc, &basingw_conc_lay, numsubstances);

  decldiag("basingw_s", BASIN, "cumulative basin groundwater outflow", "(m^3/s)", &basingw_s);

  declstatvar("cumbasingw", BASIN, "cumulative basin groundwater outflow", "(m^3)", &cumbasingw);

  declstatvar("cumbasingw_mWQ", NDEF, "cumulative mass of solute basin groundwater outflow", "(m^3)", &cumbasingw_mWQ, &cumbasingw_mWQ_lay, numsubstances);


  decllocal("soil_ssr_Buf", NHRU, "buffer subsurface runoff", "(mm/d)", &soil_ssr_Buf);

  declvar("soil_ssr_Buf_conc", NDEFN, "buffer subsurface runoff", "(mm/d)", &soil_ssr_Buf_conc, &soil_ssr_Buf_conc_lay, numsubstances);

  decllocal("soil_runoff_Buf", NHRU, "buffer rain runoff", "(mm/d)", &soil_runoff_Buf);

  declvar("soil_runoff_Buf_conc", NDEFN, "buffer rain runoff", "(mm/d)", &soil_runoff_Buf_conc, &soil_runoff_Buf_conc_lay, numsubstances);

  decllocal("soil_gw_Buf", NHRU, "buffer rain runoff", "(mm/d)", &soil_gw_Buf);

  declvar("soil_gw_Buf_conc", NDEFN, "buffer soil_gw(gw_flow) runoff", "(mm/d)", &soil_gw_Buf_conc, &soil_gw_Buf_conc_lay, numsubstances);

  decllocal("distrib_sum", NHRU, "HRU distribution sum", "()", &distrib_sum);


  declparam("basin_area", BASIN, "3", "1e-6", "1e09", "Total basin area", "(km^2)", &basin_area);

  declparam("hru_area", NHRU, "[1]", "1e-6", "1e09", "HRU area", "(km^2)", &hru_area);

  declparam("Kstorage", NHRU, "[0.0]", "0.0","200.0", "HRU storage constant", "(d)", &Kstorage);

  declparam("Lag", NHRU, "[0.0]", "0.0","1.0E4.0", "HRU lag delay", "(h)", &Lag);

  declparam("ssrKstorage", NHRU, "[0.0]", "0.0","200.0", "subsurface runoff storage constant", "(d)", &ssrKstorage);

  declparam("ssrLag", NHRU, "[0.0]", "0.0","1.0E4.0", "subsurface runoff lag delay", "(h)", &ssrLag);

  declparam("runKstorage", NHRU, "[0.0]", "0.0","200.0", "runoff storage constant", "(d)", &runKstorage);

  declparam("runLag", NHRU, "[0.0]", "0.0","1.0E4", "runoff lag delay", "(h)", &runLag);

  declparam("gwKstorage", NHRU, "[0.0]", "0.0","200.0", "gw storage constant", "(d)", &gwKstorage);

  declparam("gwLag", NHRU, "[0.0]", "0.0","1.0E4", "gw lag delay", "(h)", &gwLag);

  declparam("gwwhereto", NHRU, "[0]", "-1000", "1000", "send to: 0 - basingw, >0 - other HRU surface input <0 - other abs(-HRU) gw input, or (< -HRUmax or > +HRUmax) - surface basinflow", "()", &gwwhereto);

  declparam("order", NHRU, "[1,2,3,4,5!]", "1","1000", "HRU routing process order", "()", &order);

  declparam("distrib_Route", NDEFN, "[0.0]", "-1.0E6.0", "1.0E6.0", "route this HRU to these HRUs", "()", &distrib, &distrib_hru, nhru);

  declparam("distrib_Basin", NHRU, "[1.0]", "0.0", "100.0", "route this HRU to basin (and other HRU(s) determined by 'distrib_Route')", "()", &distrib_Basin);

  declparam("Sdmax", NHRU, "[0]", "0.0", "1000.0","Maximum depression storage", "(mm)", &Sdmax);

  declparam("soil_rechr_max", NHRU, "[60.0]", "0.0", "350.0", "soil recharge maximum (<= soil_moist_max).", "(mm)", &soil_rechr_max);

  declparam("Sd_ByPass", NHRU, "[0]", "0", "1","0 - normal, 1 - Bypass Pond/Depressional storage (i.e. Sd).", "()", &Sd_ByPass);

  declparam("soil_rechr_ByPass", NHRU, "[1]", "0", "1","0 - normal, 1 - Bypass recharge layer (i.e. soil_rechr).", "()", &soil_rechr_ByPass);


  declparam("Channel_shp", NHRU, "[0]", "0", "2", "rectangular - 0/parabolic - 1/triangular - 2", "()", &route_Cshp);

  declparam("order", NHRU, "[1,2,3,4,5!]", "1","1000", "HRU routing process order", "()", &order);

  declparam("preferential_flow", NHRU, "[0]", "0", "1","0 - no preferential and remain as runoff routing to other HRU, 1 - preferential flow and route runoff to other HRU's gw.", "()", &preferential_flow);


  soil_gwDiv = declgetvar("*", "gw_flow", "(mm/int)", &soil_gw);

  soil_ssrDiv = declgetvar("*", "soil_ssr", "(mm/int)", &soil_ssr);

  soil_runoffDiv = declgetvar("*", "soil_runoff", "(mm/int)", &soil_runoff);

  declgetvar("*", "soil_ssr_conc", "(g)", &soil_ssr_conc, &soil_ssr_conc_lay);

  declgetvar("*", "soil_gw_conc", "(mg)", &soil_gw_conc, &soil_gw_conc_lay);

  declgetvar("*", "soil_runoff_conc", "(mg)", &soil_runoff_conc, &soil_runoff_conc_lay);


  declputvar("*", "Sd", "(mm)", &Sd);

  declputvar("*", "Sd_conc", "(g)", &Sd_conc, &Sd_conc_lay);

  declputvar("*", "soil_moist", "(mm)", &soil_moist);

  declputvar("*", "soil_moist_conc", "(g)", &soil_moist_conc, &soil_moist_conc_lay);

  declputvar("*", "soil_lower", "(mm)", &soil_lower);

  declputvar("*", "conc_bottom", "(g)", &conc_bottom, &conc_bottom_lay);

  declputvar("*", "soil_rechr", "(mm)", &soil_rechr);

  declputvar("*", "redirected_residual", "(mm*km^2/int)", &redirected_residual);

  declputvar("*", "redirected_residual_conc", "(g)", &redirected_residual_conc, &redirected_residual_conc_lay);

  declputvar("*", "cum_redirected_residual", "(mm*km^2/int)", &cum_redirected_residual);

  declputvar("*", "cum_redirected_residual_mWQ", "(g)", &cum_redirected_residual_mWQ, &cum_redirected_residual_mWQ_lay);

  declputvar("*", "gw", "(mm)", &gw);

  declputvar("*", "gw_conc", "(mg/l)", &gw_conc, &gw_conc_lay);

  declputvar("*", "conc_top", "(mg/l)", &conc_top, &conc_top_lay);

  declputvar("*", "conc_bottom", "(mg/l)", &conc_bottom, &conc_bottom_lay);


  variation_set = VARIATION_0;

  decllocal("Ktravel", NHRU, "travel time", "(d)", &Ktravel);

  declparam("route_n", NHRU, "[0.025]", "0.016","0.2", "Manning roughness coefficient", "()", &route_n);

  declparam("route_R", NHRU, "[0.5]", "0.01","1.0E4", "hydraulic radius", "(m)", &route_R);

  declparam("route_S0", NHRU, "[1e-3]", "1e-6","1.0", "longitudinal channel slope", "()", &route_S0);

  declparam("route_L", NHRU, "[200.0]", "0.01","1.0E10", "routing length", "(m)", &route_L);

  declparam("route_X_M", NHRU, "[0.25]", "0.0","0.5", "dimensionless weighting factor", "()", &route_X_M);

  declparam("Channel_shp", NHRU, "[0]", "0", "2", "rectangular - 0/parabolic - 1/triangular - 2", "()", &route_Cshp);


  variation_set = VARIATION_1;

  declparam("Kstorage", NHRU, "[0.0]", "0.0","200.0", "aggregated storage constant", "(d)", &Kstorage);


  variation_set = VARIATION_ORG;

}

void ClassWQ_Netroute_M_D::init(void) {

  nhru = getdim(NHRU);

  try {
    ssrDelay_mWQ = new ClassClark*[numsubstances];
    runDelay_mWQ = new ClassClark*[numsubstances];
    gwDelay_mWQ  = new ClassClark*[numsubstances];

    ssrDelay = new ClassClark(ssrinflow, ssroutflow, ssrKstorage, ssrLag, nhru);
    runDelay = new ClassClark(runinflow, runoutflow, runKstorage, runLag, nhru);
    gwDelay = new ClassClark(gwinflow, gwoutflow, gwKstorage, gwLag, nhru);

    for(Sub = 0; Sub < numsubstances; ++Sub){
      ssrDelay_mWQ[Sub] = new ClassClark(ssrinflow_mWQ_lay[Sub], ssroutflow_mWQ_lay[Sub], ssrKstorage, ssrLag, nhru);
      runDelay_mWQ[Sub] = new ClassClark(runinflow_mWQ_lay[Sub], runoutflow_mWQ_lay[Sub], runKstorage, runLag, nhru);
      gwDelay_mWQ[Sub] = new ClassClark(gwinflow_mWQ_lay[Sub], gwoutflow_mWQ_lay[Sub], gwKstorage, gwLag, nhru);
    }
  }
  catch (std::bad_alloc) {
    CRHMException Except("Could not allocate in module Netroute_M_D." ,TERMINATE);
    LogError(Except);
    throw Except;
  }

  if(soil_ssrDiv > 1){
    String S = "WQ_Netroute_M_D:  \"soil_ssr\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(soil_runoffDiv > 1){
    String S = "WQ_Netroute_M_D:  \"soil_runoff\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(soil_gwDiv > 1){
    String S = "WQ_Netroute_M_D:  \"gw_flow\". Converting to mm/int";
    CRHMException TExcept(S.c_str(), WARNING);
    LogError(TExcept);
  }

  if(variation == VARIATION_ORG){
    const float Vw[3] = {1.67, 1.44, 1.33}; // rectangular - 0/parabolic - 1/triangular - 2

    for(hh = 0; hh < nhru; ++hh){
      float Vavg = (1.0/route_n[hh])*pow(route_R[hh], 2.0/3.0)*pow(route_S0[hh], 0.5f); // (m/s)
      Ktravel[hh] = route_L[hh]/(Vw[route_Cshp[hh]]*Vavg)/86400.0; // (d)
    }

    hruDelay = new ClassMuskingum(inflow, outflow, Ktravel, route_X_M, Lag, nhru);

    hruDelay_mWQ = new ClassMuskingum*[numsubstances];

    for(Sub = 0; Sub < numsubstances; ++Sub)
      hruDelay_mWQ[Sub] = new ClassMuskingum(inflow_mWQ_lay[Sub], outflow_mWQ_lay[Sub], Ktravel, route_X_M, Lag, nhru);

    for(hh = 0; hh < nhru; ++hh){
      if(Ktravel[hh] >= (Global::Interval/(2.0*route_X_M[hh]))){
        String S = string("'" + Name + " (Netroute_M_D) Muskingum coefficient negative in HRU ").c_str() + IntToStr(hh+1);
        CRHMException TExcept(S.c_str(), WARNING);
        LogError(TExcept);
      }

      if(Ktravel[hh] < (Global::Interval/(2.0*(1.0-route_X_M[hh])))){ //    if(hruDelay->c0[hh] < 0.0)
        hruDelay->c0[hh] = 0.0;
        hruDelay->c1[hh] = 1.0;
        hruDelay->c2[hh] = 0.0;

        for(Sub = 0; Sub < numsubstances; ++Sub){
          hruDelay_mWQ[Sub]->c0[hh] = 0.0;
          hruDelay_mWQ[Sub]->c1[hh] = 1.0;
          hruDelay_mWQ[Sub]->c2[hh] = 0.0;
        }
      }
    }
  }
  else if(variation == VARIATION_1){
    Clark_hruDelay = new ClassClark(inflow, outflow, Kstorage, Lag, nhru);

    Clark_hruDelay_mWQ = new ClassClark*[numsubstances];
    for(Sub = 0; Sub < numsubstances; ++Sub)
      Clark_hruDelay_mWQ[Sub] = new ClassClark(inflow_mWQ_lay[Sub], outflow_mWQ_lay[Sub], Kstorage, Lag, nhru);

  }

  Reset_Basin_WQ(0, basinflow, basinflow_conc);
  Reset_Basin_WQ(0, cumbasinflow, cumbasinflow_mWQ);
  Reset_Basin_WQ(0, basingw, basingw_conc);
  Reset_Basin_WQ(0, cumbasingw, cumbasingw_mWQ);

  basinflow_s[0] = 0.0;
  basingw_s[0] = 0.0;

  for(hh = 0; hh < nhru; ++hh) {
    Reset_WQ(hh, inflow, inflow_mWQ_lay);
    Reset_WQ(hh, cuminflow, cuminflow_mWQ_lay);

    Reset_WQ(hh, outflow, outflow_mWQ_lay);
    Reset_WQ(hh, cumoutflow, cumoutflow_mWQ_lay);

    Reset_WQ(hh, gwinflow, gwinflow_mWQ_lay); ;
    Reset_WQ(hh, gwcuminflow, gwcuminflow_mWQ_lay); ;

    Reset_WQ(hh, gwoutflow, gwoutflow_mWQ_lay);
    Reset_WQ(hh, gwcumoutflow, gwcumoutflow_mWQ_lay);

    Reset_WQ(hh, ssrinflow, ssrinflow_mWQ_lay);
    Reset_WQ(hh, ssrcuminflow, ssrcuminflow_mWQ_lay);

    Reset_WQ(hh, ssroutflow, ssroutflow_mWQ_lay);
    Reset_WQ(hh, ssrcumoutflow, ssrcumoutflow_mWQ_lay);

    Reset_WQ(hh, runinflow, runinflow_mWQ_lay);
    Reset_WQ(hh, runcuminflow, runcuminflow_mWQ_lay);

    Reset_WQ(hh, runoutflow, runoutflow_mWQ_lay);
    Reset_WQ(hh, runcumoutflow, runcumoutflow_mWQ_lay);

    Reset_WQ(hh, outflow_diverted, outflow_diverted_conc_lay);
    Reset_WQ(hh, cumoutflow_diverted, cumoutflow_diverted_mWQ_lay);

    Reset_WQ(hh, gwoutflow_diverted, gwoutflow_diverted_conc_lay);
    Reset_WQ(hh, gwcumoutflow_diverted, gwcumoutflow_diverted_mWQ_lay);

    Reset_WQ(hh, cum_to_Sd, cum_to_Sd_mWQ_lay);

    Reset_WQ(hh, cum_to_soil_rechr, cum_to_soil_rechr_mWQ_lay);

    Reset_WQ(hh, HRU_cumbasinflow, HRU_cumbasinflow_mWQ_lay);

    cum_preferential_flow_to_gw[hh] = 0.0;
    soil_ssr_Buf[hh] = 0.0;
    soil_runoff_Buf[hh] = 0.0;
    soil_gw_Buf[hh] = 0.0;


    boolean OK = false;
    for(long jj = 0; chkStruct(jj); ++jj)
      if(order[jj] - 1 == hh){
        OK = true;
        break;
      }

    if(!OK){
        string SS = string("'" + Name + " (Netroute)' the 'order' parameter does not have a unique value for each HRU");
        CRHMException Except(SS.c_str() ,ERR);
        LogError(Except);
        throw Except;
    }
  } // for hh
}

void ClassWQ_Netroute_M_D::run(void) {

  long step = getstep();
  long nstep = step%Global::Freq;

  float gw_Amount = 0.0;
  float gw_Amount_mWQ = 0.0;

  float Amount = 0.0;
  float Amount_mWQ = 0.0;

  Reset_Basin_WQ(0, basinflow, basinflow_conc);
  Reset_Basin_WQ(0, basingw, basingw_conc);

  for(hh = 0; chkStruct(hh); ++hh) { // do HRUs in sequence.
    if(nstep == 1){
      distrib_sum[hh] = 0.0;

      for(long hhh = 0; chkStruct(hhh); ++hhh) { // do HRUs in sequence
        if(distrib_hru[hh][hhh] < 0.0)
          const_cast<float **> (distrib_hru) [hh][hhh] = -distrib_hru[hh][hhh]*hru_area[hh];
        distrib_sum[hh] += distrib_hru[hh][hhh];
      }

      if(distrib_sum[hh] <= 0 && distrib_Basin[hh] <= 0.0){
        const_cast<float *> (distrib_Basin) [hh] = 1;
      }

      distrib_sum[hh] += distrib_Basin[hh];
    }

    for(long Sub = 0; Sub < numsubstances; Sub++){

      if(soil_gwDiv == 1){ // interval value
         soil_gw_Buf[hh] = soil_gw[hh];
         soil_gw_Buf_conc_lay[Sub][hh] = soil_gw_conc_lay[Sub][hh];
      }
      if(soil_ssrDiv == 1){ // interval value
         soil_ssr_Buf[hh] = soil_ssr[hh];
         soil_ssr_Buf_conc_lay[Sub][hh] = soil_ssr_conc_lay[Sub][hh];
      }

      if(soil_runoffDiv == 1){ // interval value
         soil_runoff_Buf[hh] = soil_runoff[hh];
         soil_runoff_Buf_conc_lay[Sub][hh] = soil_runoff_conc_lay[Sub][hh];
      }
    } // Sub
  } // hh

  for(long Sub = 0; Sub < numsubstances; Sub++){
    for(long jj = 0; chkStruct(jj); ++jj){ // HRUs not in sequence

      for(hh = 0; chkStruct(hh); ++hh)
        if(order[hh] - 1 == jj)
          break;

      gw_Amount = 0.0;
      gw_Amount_mWQ = 0.0;

      Amount = 0.0;
      Amount_mWQ = 0.0;

      gwinflow[hh] = soil_gw_Buf[hh]*hru_area[hh];

      gwoutflow_diverted[hh] = 0.0;

      gwinflow_mWQ[hh] = soil_gw_conc_lay[Sub][hh]*soil_gw_Buf[hh]*hru_area[hh];
      gwoutflow_diverted_conc_lay[Sub][hh] = 0.0;

      for(long hhh = 0; chkStruct(hhh); ++hhh) {
        if(gwoutflow[hhh] > 0.0 && gwwhereto[hhh] && (abs(gwwhereto[hhh])-1 == hh || abs(gwwhereto[hhh]) > nhru)){ // handles "gwwhereto" <> 0
          gwoutflow_diverted[hhh] = gwoutflow[hhh]; // gwoutflow_diverted[hh] = gwoutflow[hhh];

          gw_Amount = gwoutflow[hhh]/hru_area[hh]; // units (mm*km^2/int)
          gw_Amount_mWQ = gwoutflow_mWQ_lay[Sub][hhh]/hru_area[hh]; // units (mm*km^2/int)

          gwcumoutflow_diverted[hhh] += gw_Amount;
          gwcumoutflow_diverted_mWQ_lay[Sub][hhh] += gw_Amount_mWQ;

          if(abs(gwwhereto[hhh]) <= nhru){
            if(gwwhereto[hhh] > 0){ // direct to HRU surface
              float free = soil_rechr_max[hh] - soil_rechr[hh];
              float free_mWQ = Amount_mWQ*free/gw_Amount;

              if(free > 0.0 && !soil_rechr_ByPass[hh]){
                if(free > gw_Amount){ // units (mm*km^2/int)

                  conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + Amount_mWQ;
                  soil_moist_conc_lay[Sub][hh] = soil_moist_conc_lay[Sub][hh]*soil_moist[hh] + Amount_mWQ;

                  conc_top_lay[Sub][hh] /= (soil_rechr[hh] + gw_Amount);
                  soil_moist_conc_lay[Sub][hh] /= (soil_moist[hh] + gw_Amount);

                  cum_to_soil_rechr_mWQ_lay[Sub][hh] += gw_Amount_mWQ;

                  if(Sub == numsubstances-1){
                    soil_rechr[hh] += gw_Amount;
                    soil_moist[hh] += gw_Amount;
                    cum_to_soil_rechr[hh] += gw_Amount;
                  }

                  gw_Amount = 0.0;
                  gw_Amount_mWQ = 0.0;
                }
                else {
                  gw_Amount_mWQ = gw_Amount_mWQ*(gw_Amount - free)/Amount;
                  gw_Amount = (gw_Amount - free);

                  conc_top_lay[Sub][hh] = conc_top_lay[Sub][hh]*soil_rechr[hh] + free_mWQ;
                  soil_moist_conc_lay[Sub][hh] = soil_moist_conc_lay[Sub][hh]*soil_moist[hh] + free_mWQ;

                  conc_top_lay[Sub][hh] /= soil_rechr_max[hh];
                  soil_moist_conc_lay[Sub][hh] /= (soil_moist[hh] + gw_Amount);

                  cum_to_soil_rechr_mWQ_lay[Sub][hh] += free_mWQ;

                  if(Sub == numsubstances-1){
                    soil_rechr[hh] = soil_rechr_max[hh];
                    soil_moist[hh] = soil_moist[hh] + free_mWQ;
                    cum_to_soil_rechr[hh] += free;
                  }
                }
              }

              free = Sdmax[hh] - Sd[hh];
              free_mWQ = Amount_mWQ*free/gw_Amount;

              if(free > 0.0 && !Sd_ByPass[hh] && gw_Amount > 0.0){
                if(free > gw_Amount){ // units (mm*km^2/int)
                  Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + Amount_mWQ;

                  Sd_conc_lay[Sub][hh] /= (Sd[hh] + gw_Amount);
                  cum_to_Sd_mWQ_lay[Sub][hh] += Amount_mWQ;

                  if(Sub == numsubstances-1){
                    Sd[hh] += gw_Amount;
                    cum_to_Sd[hh] += gw_Amount;
                  }

                  gw_Amount = 0.0;
                  gw_Amount_mWQ = 0.0;
                }
                else {
                  gw_Amount_mWQ = gw_Amount_mWQ*(gw_Amount - free)/Amount;

                  Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]*Sd[hh] + free_mWQ;
                  cum_to_Sd_mWQ_lay[Sub][hh] += free_mWQ;
                  Sd_conc_lay[Sub][hh] = Sd_conc_lay[Sub][hh]/Sdmax[hh];

                  if(Sub == numsubstances-1){
                    gw_Amount = (gw_Amount - free);
                    Sd[hh] = Sdmax[hh];
                    cum_to_Sd[hh] += free;
                  }
                }
              }
            } // hh > 0
            else{ // hh < 0
              gw_conc_lay[Sub][hh] = gw_conc_lay[Sub][hh]*gw[hh] + gw_Amount_mWQ*gw_Amount;
              gw_conc_lay[Sub][hh] /= (gw[hh] + gw_Amount);

              if(Sub == numsubstances-1)
                gw[hh] += gw_Amount;

              gw_Amount = 0.0;
              gw_Amount_mWQ = 0.0;
            }
          }
          else if(gwwhereto[hh] == 0){ // move to basin gw
            basingw_conc_lay[Sub][0] = basingw_conc_lay[Sub][0]*basingw[0] + gwoutflow_mWQ_lay[Sub][hh]*1000;
            float basingw0 = basingw[0] + gw_Amount*1000; // (m3) end of every day

            if(basingw0 > 0.0)
              basingw_conc_lay[Sub][0] /= basingw0;
            else
              basingw_conc_lay[Sub][0] = 0.0;

            gwcumoutflow_mWQ_lay[Sub][hh] += gw_Amount_mWQ;

            if(Sub == numsubstances-1){
              basingw[0] = basingw0;
              gwcumoutflow[hh] += gw_Amount;
            }

            gw_Amount = 0.0;
            gw_Amount_mWQ = 0.0;
          }
          else{
            if(Sub == numsubstances-1){
              basinflow[0] += gw_Amount*hru_area[hh]*1000; // (m3)
              cumoutflow[hh] += gw_Amount*hru_area[hh];
              HRU_cumbasinflow[hh] +=  gw_Amount;
            }

            gw_Amount = 0.0;
            gw_Amount_mWQ = 0.0;
          } // is HRU in range
        } // handles "gwwhereto" <> 0
      } // for hhh

      gwcuminflow[hh] += gwinflow[hh];

      runinflow[hh] = soil_runoff_Buf[hh]*hru_area[hh];
      ssrinflow[hh] = soil_ssr_Buf[hh]*hru_area[hh];
      inflow[hh] = runoutflow[hh] + ssroutflow[hh] + gw_Amount; // add this HRU runoff and subsurface flow

      runcuminflow[hh] += runinflow[hh];
      runcumoutflow[hh] += runoutflow[hh];
      ssrcuminflow[hh] += ssrinflow[hh];
      ssrcumoutflow[hh] += ssroutflow[hh];
      cuminflow[hh] += inflow[hh];

      runinflow_mWQ_lay[Sub][hh] = soil_runoff_Buf_conc_lay[Sub][hh]*soil_runoff_Buf[hh]*hru_area[hh];
      runcuminflow_mWQ_lay[Sub][hh] += runinflow_mWQ_lay[Sub][hh];
      runcumoutflow_mWQ_lay[Sub][hh] += runoutflow_mWQ_lay[Sub][hh];

      ssrinflow_mWQ_lay[Sub][hh] = soil_ssr_Buf_conc_lay[Sub][hh]*soil_ssr_Buf[hh]*hru_area[hh];
      ssrcuminflow_mWQ_lay[Sub][hh] += ssrinflow_mWQ_lay[Sub][hh];
      ssrcumoutflow_mWQ_lay[Sub][hh] += ssroutflow_mWQ_lay[Sub][hh];

      inflow_mWQ_lay[Sub][hh] = runoutflow_mWQ_lay[Sub][hh] + ssroutflow_mWQ_lay[Sub][hh] + Amount_mWQ; // add this HRU runoff and subsurface flow

      if(outflow[hh] > 0.0){
        Amount = outflow[hh]/hru_area[hh]; // unit area

        Reset_WQ(hh, Used, Used_mWQ_lay);

        if(distrib_Basin[hh] > 0.0){ // direct to basin
          if(Amount > 0.0){
            Used[hh] = Amount*distrib_Basin[hh]/distrib_sum[hh]; // amount over basin area
            Used_mWQ_lay[Sub][hh] = outflow_mWQ_lay[Sub][hh]*distrib_Basin[hh]/distrib_sum[hh]; // fraction used

            basinflow_conc_lay[Sub][0] = Used_mWQ_lay[Sub][hh]/Used[hh];

            cumoutflow_mWQ_lay[Sub][hh] += Used_mWQ_lay[Sub][hh];

            HRU_cumbasinflow_mWQ_lay[Sub][hh] +=  Used_mWQ_lay[Sub][hh];

            if(Sub == numsubstances-1){
              basinflow[0] += Used[hh]*1000; // (m3)
              cumoutflow[hh] += Used[hh];
              HRU_cumbasinflow[hh] +=  Used[hh];
            }
          }
        }

        for(long To = 0; chkStruct(To); ++To) {

          if(hh != To && distrib_hru[hh][To] > 0.0){
              Amount = (outflow[hh] - Used[hh])/hru_area[To]*distrib_hru[hh][To]/(distrib_sum[hh]-distrib_Basin[hh]); // outflow (mm*km^2/int)
              Amount_mWQ = (outflow_mWQ_lay[Sub][hh] - Used_mWQ_lay[Sub][hh])/hru_area[To]*distrib_hru[hh][To]/(distrib_sum[hh] - distrib_Basin[hh]);

              if(preferential_flow[To]) {
                gw_conc_lay[Sub][To] = gw_conc_lay[Sub][To]*gw[To];

                gw_conc_lay[Sub][To] += Amount_mWQ;
                gw_conc_lay[Sub][To] /= gw[To];

                if(Sub == numsubstances-1){
                  gw[To] += Amount;
                  cum_preferential_flow_to_gw[To] += Amount;
                  Amount = 0.0;
                  Amount_mWQ = 0.0;
                }

                Amount = 0.0;
                Amount_mWQ = 0.0;
              }
              else{
                if(!soil_rechr_ByPass[To]){
                  if(soil_rechr[To] + Amount >= soil_rechr_max[To]){ // units (mm*km^2/int)
                    float Excess = soil_rechr[To] + Amount - soil_rechr_max[To];
                    float Free = Amount - Excess;

                    conc_top_lay[Sub][To] = conc_top_lay[Sub][To]*soil_rechr[To] + Amount_mWQ;
                    conc_top_lay[Sub][To] /= (soil_rechr[To] + Amount);

                    soil_moist_conc_lay[Sub][To] = (conc_bottom_lay[Sub][To]*soil_lower[To] + conc_top_lay[Sub][To]*soil_rechr[To] + Amount_mWQ*Free)/(soil_lower[To] + soil_rechr[To] + Free); // present mQW

                    if(Sub == numsubstances-1){
                      soil_rechr[To] += Free;
                      soil_moist[To] = soil_lower[To] + soil_rechr[To];
                    }
                    Amount = Excess;
                    Amount_mWQ = conc_top_lay[Sub][To]*Excess/Amount;
                  }
                  else{
                    conc_top_lay[Sub][To] = conc_top_lay[Sub][To]*soil_rechr[To] + Amount_mWQ;
                    if(soil_rechr[To] + Amount > 0.0)
                      conc_top_lay[Sub][To] /= (soil_rechr[To] + Amount);
                    else
                      conc_top_lay[Sub][To] = 0.0;

                    soil_moist_conc_lay[Sub][To] = (conc_bottom_lay[Sub][To]*soil_lower[To] + conc_top_lay[Sub][To]*(soil_rechr[To] + Amount_mWQ))/(soil_lower[To] + soil_rechr[To] + Amount); // amount used

                    if(Sub == numsubstances-1){
                      soil_rechr[To] = soil_rechr[To] + Amount;
                      soil_moist[To] = soil_lower[To] + soil_rechr[To];
                    }

                    Amount = 0.0;
                    Amount_mWQ = 0.0;
                  }
                }

                if(Amount > 0.0 && !Sd_ByPass[To]){
                  if(Sd[To] + Amount >= Sdmax[To]){ // units (mm*km^2/int)
                    float Excess = Sd[To] + Amount - Sdmax[To];
                    float Free = Amount - Excess;

                    Sd_conc_lay[Sub][To] = Sd_conc_lay[Sub][To]*Sd[To] + Amount_mWQ;
                    Sd_conc_lay[Sub][To] /= (Sd[To] + Amount);

                    cum_to_Sd_mWQ_lay[Sub][To] += Amount_mWQ;

                    if(Sub == numsubstances-1){
                      Sd[To] += Free;
                      cum_to_Sd[To] += Amount;
                    }

                    Amount = Excess;
                    Amount_mWQ = Sd_conc_lay[Sub][To]*Excess;
                  }
                  else{
                    Sd_conc_lay[Sub][To] = Sd_conc_lay[Sub][To]*Sd[To] + Amount_mWQ;

                    if(Sd[To] + Amount > 0.0)
                      Sd_conc_lay[Sub][To] /= (Sd[To] + Amount);
                    else
                      Sd_conc_lay[Sub][To] = 0.0;

                    if(Sub == numsubstances-1){
                      Sd[To] = Sd[To] + Amount;
                      cum_to_Sd[To] += Amount;
                    }

                    Amount = 0.0;
                    Amount_mWQ = 0.0;
                  } // else
                } // if !Sd_ByPass
              } // else

              if(Amount > 0.0){ // handle excess
                  redirected_residual_conc_lay[Sub][To] = Amount_mWQ*hru_area[To]/(redirected_residual[To] + Amount*hru_area[To]);

                if(Sub == numsubstances-1)
                  redirected_residual[To] += Amount*hru_area[To];
              }
          } // contribute to this array
        } // distribute outflow over HRUs (for)
      } // if outflow > 0.0

      if(nstep == 0){ // end of every day
        if(soil_ssrDiv > 1) // daily value - ready for next day
           soil_ssr_Buf[hh] = soil_ssr[hh]/soil_ssrDiv;

        if(soil_runoffDiv > 1) // daily value - ready for next day
           soil_runoff_Buf[hh] = soil_runoff[hh]/soil_runoffDiv;

        if(soil_gwDiv > 1) // daily value - ready for next day
           soil_gw_Buf[hh] = soil_gw[hh]/soil_gwDiv;
      } // end if

      if(variation == VARIATION_ORG){
        hruDelay->DoMuskingum(hh); // need to update for later HRUs
        for(long Sub = 0; Sub < numsubstances; Sub++)
          hruDelay_mWQ[Sub]->DoMuskingum(hh); // need to update for later HRUs
      }
      else if(variation == VARIATION_1){
        Clark_hruDelay->DoClark(hh); // need to update for later HRUs
        for(long Sub = 0; Sub < numsubstances; Sub++)
          Clark_hruDelay_mWQ[Sub]->DoClark(hh); // need to update for later HRUs
      }

    } // for jj accessing hh

    for(hh = 0; chkStruct(hh); ++hh){
      cuminflow[hh] += inflow[hh];

      outflow_diverted[hh] = 0.0;
      if(distrib_sum[hh] > 0.0){ // does not apply to last HRU
        if(Sub == numsubstances-1)
          for(long hhh = 0; chkStruct(hhh); ++hhh){
            outflow_diverted[hh] += outflow[hh]*distrib_hru[hh][hhh]/distrib_sum[hh];
            cumoutflow_diverted[hh] += outflow_diverted[hh];
          }
      }

        cumoutflow_diverted_mWQ_lay[Sub][hh] += outflow_diverted_conc[hh];

        cumbasinflow_mWQ_lay[Sub][0] += basinflow[0]*outflow_mWQ_lay[Sub][hh];
        cumbasingw_mWQ_lay[Sub][0] += basingw[0]*gwoutflow_mWQ_lay[Sub][hh];
    } // end for
  } // Sub
/*
  if(variation == VARIATION_ORG){
    hruDelay->DoMuskingum(); // need to update for later HRUs
    for(long Sub = 0; Sub < numsubstances; Sub++)
      hruDelay_cWQ[Sub]->DoMuskingum(); // need to update for later HRUs
  }
  else if(variation == VARIATION_1){
    Clark_hruDelay->DoClark(); // need to update for later HRUs
    for(long Sub = 0; Sub < numsubstances; Sub++)
      Clark_hruDelay_cWQ[Sub]->DoClark(); // need to update for later HRUs
  }
*/
  runDelay->DoClark();
  ssrDelay->DoClark();
  gwDelay->DoClark();

  for(long Sub = 0; Sub < numsubstances; Sub++){
    runDelay_mWQ[Sub]->DoClark();
    ssrDelay_mWQ[Sub]->DoClark();
    gwDelay_mWQ[Sub]->DoClark();
  }
  basinflow_s[0] = basinflow[0]*Global::Freq/86400.0;
  basingw_s[0] = basingw[0]*Global::Freq/86400.0;

  cumbasinflow[0] += basinflow[0];
  cumbasingw[0] += basingw[0];
} // run

void ClassWQ_Netroute_M_D::finish(bool good){

  float Allcuminflow = 0.0;
  float Allcumoutflow = 0.0;
  float Allcumoutflowdiverted = 0.0;

  float Allgwcuminflow = 0.0;
  float Allgwcumoutflow = 0.0;
  float Allgwcumoutflowdiverted = 0.0;

  float Allssrcuminflow = 0.0;
  float Allssrcumoutflow = 0.0;
  float Allruncuminflow = 0.0;
  float Allruncumoutflow = 0.0;

  float Allcuminflow_mWQ = 0.0;
  float Allcumoutflow_mWQ = 0.0;
  float Allcumoutflowdiverted_mWQ = 0.0;

  float Allgwcuminflow_mWQ = 0.0;
  float Allgwcumoutflow_mWQ = 0.0;
  float Allgwcumoutflowdiverted_mWQ = 0.0;

  float Allssrcuminflow_mWQ = 0.0;
  float Allssrcumoutflow_mWQ = 0.0;
  float Allruncuminflow_mWQ = 0.0;
  float Allruncumoutflow_mWQ = 0.0;

  float AllSdcuminflow = 0.0;
  float Allrechrcuminflow = 0.0;

  float AllSdcuminflow_mWQ = 0.0;
  float Allrechrcuminflow_mWQ = 0.0;
  float AllTotal = 0.0;
  float Total;

  String S = String("H2O");
  LogDebug(S.c_str());
  LogMessage(" ");

  for(hh = 0; chkStruct(); ++hh) {
    LogMessageA(hh, string("'" + Name + " (Netroute_M_D_WQ)' cuminflow                   (mm) (mm*km^2) (mm*basin): ").c_str(), cuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumoutflow                  (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumoutflow_diverted         (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);

    if(variation == VARIATION_ORG)
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' hruDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), hruDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);
    else if(variation == VARIATION_1)
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' hruDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), Clark_hruDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrcuminflow                (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrcumoutflow               (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), ssrDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runoffcuminflow             (mm) (mm*km^2) (mm*basin): ").c_str(), runcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runoffcumoutflow            (mm) (mm*km^2) (mm*basin): ").c_str(), runcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), runDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcuminflow                 (mm) (mm*km^2) (mm*basin): ").c_str(), gwcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcumoutflow                (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcumoutflow_diverted       (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwDelay_in_storage          (mm) (mm*km^2) (mm*basin): ").c_str(), gwDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_to_Sd                   (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_Sd[hh], hru_area[hh], basin_area[0], " *** Added to this HRU Sd");
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_to_soil_rechr           (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_soil_rechr[hh], hru_area[hh], basin_area[0], " *** Added to this HRU recharge");
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_redirected_residual     (mm) (mm*km^2) (mm*basin): ").c_str(), cum_redirected_residual[hh], hru_area[hh], basin_area[0]);
//    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumbasinflow                (mm) (mm*km^2) (mm*basin): ").c_str(), cumbasinflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' HRU_cumbasinflow            (mm) (mm*km^2) (mm*basin): ").c_str(), HRU_cumbasinflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");

//    Total = cumoutflow[hh] + gwcumoutflow[hh] - cumbasinflow[hh] - cum_to_soil_rechr[hh] - cum_to_Sd[hh] - gwcumoutflow[hh]
//            + cumoutflow_diverted[hh] + gwcumoutflow_diverted[hh] - cum_redirected_residual[hh];
//    AllTotal += Total;

//    LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' Total                           (mm) (mm*hru) (mm*hru/basin): ").c_str(), Total, hru_area[hh], basin_area[0], " *** HRU mass balance");
//    LogDebug(" ");

    Allcuminflow += cuminflow[hh];
    Allcumoutflow += cumoutflow[hh];

    Allgwcuminflow += gwcuminflow[hh];
    Allgwcumoutflow += gwcumoutflow[hh];

    Allssrcumoutflow += ssrcumoutflow[hh];
    Allssrcuminflow += ssrcuminflow[hh];
    Allruncuminflow += runcuminflow[hh];
    Allruncumoutflow += runcumoutflow[hh];

    AllSdcuminflow += cum_to_Sd[hh];
    Allrechrcuminflow += cum_to_soil_rechr[hh];

    Allcumoutflowdiverted += cumoutflow_diverted[hh];

    Allgwcumoutflowdiverted += gwcumoutflow_diverted[hh];

    LogMessage(" ");
  }  // for hh

  LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' cumbasinflow (m^3):                   ").c_str(), cumbasinflow[0]);
  LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' cumbasingw   (m^3):                   ").c_str(), cumbasingw[0]);
  LogMessage(" ");

  for(long Sub = 0; Sub < numsubstances; ++Sub){

    String S = String("Substance# ") + (Sub+1);
    LogDebug(S.c_str());
    LogMessage(" ");

    Allcuminflow = 0.0;
    Allcumoutflow = 0.0;
    Allcumoutflowdiverted = 0.0;

    Allgwcuminflow = 0.0;
    Allgwcumoutflow = 0.0;
    Allgwcumoutflowdiverted = 0.0;


    Allssrcuminflow = 0.0;
    Allssrcumoutflow = 0.0;
    Allruncuminflow = 0.0;
    Allruncumoutflow = 0.0;

    Allcuminflow_mWQ = 0.0;
    Allcumoutflow_mWQ = 0.0;
    Allcumoutflowdiverted_mWQ = 0.0;

    Allgwcuminflow_mWQ = 0.0;
    Allgwcumoutflow_mWQ = 0.0;
    Allgwcumoutflowdiverted_mWQ = 0.0;

    Allssrcuminflow_mWQ = 0.0;
    Allssrcumoutflow_mWQ = 0.0;
    Allruncuminflow_mWQ = 0.0;
    Allruncumoutflow_mWQ = 0.0;

    AllSdcuminflow = 0.0;
    Allrechrcuminflow = 0.0;

    AllSdcuminflow_mWQ = 0.0;
    Allrechrcuminflow_mWQ = 0.0;
    AllTotal = 0.0;

    for(hh = 0; chkStruct(); ++hh) {
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cuminflow                   (mm) (mm*km^2) (mm*basin): ").c_str(), cuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cuminflow_mWQ               (mm) (mm*km^2) (mm*basin): ").c_str(), cuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumoutflow                  (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumoutflow_mWQ              (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumoutflow_diverted         (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumoutflow_diverted_mWQ     (mm) (mm*km^2) (mm*basin): ").c_str(), cumoutflow_diverted_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);

      if(variation == VARIATION_ORG)
        LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' hruDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), hruDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);
      else if(variation == VARIATION_1)
        LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' hruDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), Clark_hruDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrcuminflow                (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrcuminflow_mWQ            (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrcumoutflow               (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrcumoutflow_mWQ           (mm) (mm*km^2) (mm*basin): ").c_str(), ssrcumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' ssrDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), ssrDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runoffcuminflow             (mm) (mm*km^2) (mm*basin): ").c_str(), runcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runoffcuminflow_mWQ         (mm) (mm*km^2) (mm*basin): ").c_str(), runcuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runoffcumoutflow            (mm) (mm*km^2) (mm*basin): ").c_str(), runcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runoffcumoutflow_mWQ        (mm) (mm*km^2) (mm*basin): ").c_str(), runcumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' runDelay_in_storage         (mm) (mm*km^2) (mm*basin): ").c_str(), runDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcuminflow                 (mm) (mm*km^2) (mm*basin): ").c_str(), gwcuminflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcuminflow_mWQ             (mm) (mm*km^2) (mm*basin): ").c_str(), gwcuminflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcumoutflow                (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcumoutflow_mWQ            (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcumoutflow_diverted       (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_diverted[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwcumoutflow_diverted_mWQ   (mm) (mm*km^2) (mm*basin): ").c_str(), gwcumoutflow_diverted_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' gwDelay_in_storage          (mm) (mm*km^2) (mm*basin): ").c_str(), gwDelay->Left(hh)/hru_area[hh], hru_area[hh], basin_area[0]);

      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_to_Sd                   (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_Sd[hh], hru_area[hh], basin_area[0], " *** Added to this HRU Sd");
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_to_Sd_mWQ               (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_Sd_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0], " *** Added to this HRU Sd");
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_to_soil_rechr           (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_soil_rechr[hh], hru_area[hh], basin_area[0], " *** Added to this HRU recharge");
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_to_soil_rechr_mWQ       (mm) (mm*km^2) (mm*basin): ").c_str(), cum_to_soil_rechr_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0], " *** Added to this HRU recharge");
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_redirected_residual     (mm) (mm*km^2) (mm*basin): ").c_str(), cum_redirected_residual[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cum_redirected_residual_mWQ (mm) (mm*km^2) (mm*basin): ").c_str(), cum_redirected_residual_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumbasinflow                (mm) (mm*km^2) (mm*basin): ").c_str(), cumbasinflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' cumbasinflow_conc           (mm) (mm*km^2) (mm*basin): ").c_str(), cumbasinflow_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' HRU_cumbasinflow            (mm) (mm*km^2) (mm*basin): ").c_str(), HRU_cumbasinflow[hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' HRU_cumbasinflow_mwq        (mm) (mm*km^2) (mm*basin): ").c_str(), HRU_cumbasinflow_mWQ_lay[Sub][hh]/hru_area[hh], hru_area[hh], basin_area[0]);
      LogDebug(" ");

      Total = cumoutflow[hh] + gwcumoutflow[hh] - cumbasinflow[hh] - cum_to_Sd_mWQ_lay[Sub][hh] - cum_to_soil_rechr[hh] - gwcumoutflow[hh]
             + cumoutflow_diverted[hh] + gwcumoutflow_diverted[hh] - cum_redirected_residual[hh];
      AllTotal += Total;

      LogMessageA(hh, string("'" + Name + " (WQ_Netroute_M_D)' Total                           (mm) (mm*hru) (mm*hru/basin): ").c_str(), Total, hru_area[hh], basin_area[0], " *** HRU mass balance");
      LogDebug(" ");

      Allgwcumoutflowdiverted += gwcumoutflow_diverted[hh];
      Allcuminflow_mWQ += cuminflow_mWQ_lay[Sub][hh];
      Allcumoutflow_mWQ += cumoutflow_mWQ_lay[Sub][hh];
      Allcumoutflowdiverted_mWQ += cumoutflow_diverted_mWQ_lay[Sub][hh];

      Allgwcuminflow_mWQ += gwcuminflow_mWQ_lay[Sub][hh];
      Allgwcumoutflow_mWQ += gwcumoutflow_mWQ_lay[Sub][hh];
      Allgwcumoutflowdiverted_mWQ += gwcumoutflow_diverted_mWQ_lay[Sub][hh];

      Allssrcumoutflow_mWQ += ssrcumoutflow_mWQ_lay[Sub][hh];
      Allssrcuminflow_mWQ += ssrcuminflow_mWQ_lay[Sub][hh];
      Allruncuminflow_mWQ += runcuminflow_mWQ_lay[Sub][hh];
      Allruncumoutflow_mWQ += runcumoutflow_mWQ_lay[Sub][hh];

      AllSdcuminflow_mWQ += cum_to_Sd_mWQ_lay[Sub][hh];
      Allrechrcuminflow_mWQ += cum_to_soil_rechr_mWQ_lay[Sub][hh];
      LogMessage(" ");

      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' cumbasinflow (m^3):                   ").c_str(), cumbasinflow[0]);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' cumbasinflow_mWQ_lay (mg):                  ").c_str(), cumbasinflow_mWQ_lay[Sub][0]);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' cumbasingw   (m^3):                   ").c_str(), cumbasingw[0]);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' cumbasingw_mWQ_lay (mg):                    ").c_str(), cumbasingw_mWQ_lay[Sub][0]);
      LogMessage(" ");

      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allgwcuminflow (mm*basin):              ").c_str(), Allgwcuminflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allgwcuminflow_mWQ (mm*basin):          ").c_str(), Allgwcuminflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allgwcumoutflow (mm*basin):             ").c_str(), Allgwcumoutflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allgwcumoutflow_mWQ (mm*basin):         ").c_str(), Allgwcumoutflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allgwcumoutflowdiverted (mm*basin):     ").c_str(), Allgwcumoutflowdiverted);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allgwcumoutflowdiverted_mWQ (mm*basin): ").c_str(), Allgwcumoutflowdiverted_mWQ);
      LogDebug(" ");

      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allcuminflow (mm*basin):                ").c_str(), Allcuminflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allcuminflow_mWQ (mm*basin):            ").c_str(), Allcuminflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allcumoutflow (mm*basin):               ").c_str(), Allcumoutflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allcumoutflow_mWQ (mm*basin):           ").c_str(), Allcumoutflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allcumoutflowdiverted (mm*basin):       ").c_str(), Allcumoutflowdiverted);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allcumoutflowdiverted_mWQ (mm*basin):   ").c_str(), Allcumoutflowdiverted_mWQ);
      LogDebug(" ");

      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allssrcuminflow (mm*basin):             ").c_str(), Allssrcuminflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allssrcuminflow_mWQ (mm*basin):         ").c_str(), Allssrcuminflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allssrcumoutflow (mm*basin):            ").c_str(), Allssrcumoutflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allssrcumoutflow_mWQ (mm*basin):        ").c_str(), Allssrcumoutflow_mWQ);
      LogDebug(" ");

      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allruncuminflow (mm*basin):             ").c_str(), Allruncuminflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allruncuminflow_mWQ (mm*basin):         ").c_str(), Allruncuminflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allruncumoutflow (mm*basin):            ").c_str(), Allruncumoutflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allruncumoutflow_mWQ (mm*basin):        ").c_str(), Allruncumoutflow_mWQ);
      LogDebug(" ");

      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' AllSdcuminflow (mm*basin):              ").c_str(), AllSdcuminflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' AllSdcuminflow_mWQ (mm*basin):          ").c_str(), AllSdcuminflow_mWQ);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allrechrcuminflow (mm*basin):           ").c_str(), Allrechrcuminflow);
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' Allrechrcuminflow_mWQ (mm*basin):       ").c_str(), Allrechrcuminflow_mWQ);
      LogDebug(" ");
      LogMessage(string("'" + Name + " (WQ_Netroute_M_D)' AllTotal              (mm*basin):       ").c_str(), AllTotal);
      LogDebug(" ");
    }  // for hh
  } // for Sub

  LogDebug(" ");
  if(variation == VARIATION_ORG){
//    for(long ii=0; ii < numsubstances; ++ii)
//      delete[] hruDelay_mWQ[ii];

    delete[] hruDelay_mWQ;
  }
  if(variation == VARIATION_1){
//    for(long ii=0; ii < numsubstances; ++ii)
//      delete[] Clark_hruDelay_mWQ[ii];

    delete[] Clark_hruDelay_mWQ;
  }
//  for(long ii=0; ii < numsubstances; ++ii){
//    delete[] ssrDelay_mWQ[ii];
//    delete[] runDelay_mWQ[ii];
//    delete[] gwDelay_mWQ[ii];
//  }

  delete[] ssrDelay_mWQ;
  delete[] runDelay_mWQ;
  delete[] gwDelay_mWQ;

  delete ssrDelay;
  delete runDelay;
  delete gwDelay;
}

void ClassWQ_Netroute_M_D::Reset_Basin_WQ(long hru, float *var, float *var_conc){
  var[hru] = 0.0;
  var_conc[hru] = 0.0;
}

void ClassWQ_Netroute_M_D::Reset_WQ(long hru, float *var, float **var_WQ_lay){
  var[hru] = 0.0;
  for(long Sub = 0; Sub < numsubstances; Sub++){
    var_WQ_lay[Sub][hru] = 0.0;
  }
}

void ClassWQ_Netroute_M_D::Set_WQ(const long hru, float *var, float *var_conc, float Amount, float amount_conc){

  var[hru] = Amount;
  if(Amount > 0.0)
    var_conc[hru] = amount_conc;
  else
    var_conc[hru] = 0.0;
}

ClassWQ_Netroute_M_D* ClassWQ_Netroute_M_D::klone(string name) const{
  return new ClassWQ_Netroute_M_D(name);
}

float ClassWQ_Netroute_M_D::Function1(float *I, long hh) {
  return runDelay->ChangeLag(I, hh);
}

float ClassWQ_Netroute_M_D::Function2(float *X, long hh) {
  return runDelay->ChangeStorage(X, hh);
}

ClassWQ_Test_Hype* ClassWQ_Test_Hype::klone(string name) const{
  return new ClassWQ_Test_Hype(name);
}

void ClassWQ_Test_Hype::decl(void) {

  Description = "'Provides necessary inputs to test the operation of CRHM and Hype separately.'\
                 'Provides necessary inputs to test WQ_soil, WQ_Netroute and Hype.'\
                 'Provides necessary inputs to test Hype alone.'";

  variation_set = VARIATION_1;

  declvar("soil_moist", NHRU, "moisture content of soil of the HRU.", "(mm)", &soil_moist);
  declvar("soil_rechr", NHRU, "moisture content of soil recharge zone of the HRU.", "(mm)", &soil_rechr);
  declparam("soil_moist_0", NHRU, "[60]", "0","1000", "initial soil moisture.", "()", &soil_moist_0);
  declparam("soil_rechr_0", NHRU, "[30]", "0","1000", "initial soil recharge.", "()", &soil_rechr_0);
  declparam("soil_moist_max", NHRU, "[375.0]", "0.0", "5000.0", "Maximum available water holding capacity of soil profile.", "(mm)", &soil_moist_max);
  declparam("soil_rechr_max", NHRU, "[375.0]", "0.0", "5000.0", "soil recharge maximum (<= soil_moist_max).", "(mm)", &soil_rechr_max);

  variation_set = VARIATION_ORG;

  declvar("infil", NHRU, "infil", "(mm/int)", &infil);
  declvar("snowinfil", NHRU, "snowinfil", "(mm/int)", &snowinfil);
  declvar("runoff", NHRU, "runoff", "(mm/int)", &runoff);
  declvar("meltrunoff", NHRU, "meltrunoff", "(mm/int)", &meltrunoff);
  declvar("hru_evap", NHRU, "hru_evap", "(mm/int)", &hru_evap);
  declvar("hru_cum_evap", NHRU, "hru_cum_evap", "(mm)", &hru_cum_evap);
  declvar("hru_cum_actet", NHRU, "hru_cum_actet", "(mm)", &hru_cum_actet);
  declvar("hru_actet", NHRU, "hru_actet", "(mm/int)", &hru_actet);
  declvar("net_rain", NHRU, "net_rain", "(mm/int)", &net_rain);
  declvar("SWE", NHRU, "SWE", "(mm)", &SWE);
  declvar("SWE_max", NHRU, "maximum seasonal SWE", "(mm)", &SWE_max);
  declvar("hru_t", NHRU, "hru_t", "(mm/int)", &hru_t);
  declvar("SWE_conc", NDEFN, "SWE_conc", "(mg/l)", &SWE_conc, &SWE_conc_lay, numsubstances);


  declparam("runoff_0", NHRU, "[0.0]", "0.0", "100.0", "runoff_0", "(mm/int)", &runoff_0);
  declparam("infil_0", NHRU, "[0.0]", "0.0", "100.0", "infil_0", "(mm/int)", &infil_0);
  declparam("snowinfil_0", NHRU, "[0.0]", "0.0", "100.0", "snowinfil_0", "(mm/int)", &snowinfil_0);
  declparam("runoff_0", NHRU, "[0.0]", "0.0", "100.0", "runoff_0", "(mm/int)", &runoff_0);
  declparam("meltrunoff_0", NHRU, "[0.0]", "0.0", "100.0", "meltrunoff_0", "(mm/int)", &meltrunoff_0);
  declparam("hru_evap_0", NHRU, "[0.0]", "0.0", "100.0", "hru_evap_0", "(mm/int)", &hru_evap_0);
  declparam("hru_actet_0", NHRU, "[0.0]", "0.0", "100.0", "hru_actet_0", "(mm/int)", &hru_actet_0);
  declparam("hru_cum_evap_0", NHRU, "[0.0]", "0.0", "100.0", "hru_cum_evap_0", "(mm/int)", &hru_cum_evap_0);
  declparam("hru_cum_actet_0", NHRU, "[0.0]", "0.0", "100.0", "hru_cum_actet_0", "(mm/int)", &hru_cum_evap_0);
  declparam("net_rain_0", NHRU, "[0.0]", "0.0", "100.0", "net rain", "(mm/int)", &net_rain_0);
  declparam("SWE_0", NHRU, "[0.0]", "0.0", "500.0", "SWE", "(mm)", &SWE_0);
  declparam("hru_t_0", NHRU, "[20.0]", "-50.0", "100.0", "hru_t_0", "(mm/int)", &hru_t_0);
  declparam("SWE_conc_0", NDEFN, "[0.0]", "0.0", "100.0", "SWE_conc_0", "(mg/l)", &SWE_conc_0, &SWE_conc_lay_0, numsubstances);
  declparam("Julian_start", NHRU, "[30]", "0","366", "enable input.", "()", &Julian_start);
  declparam("Julian_end", NHRU, "[30]", "0","366", "disable input.", "()", &Julian_end);

}

void ClassWQ_Test_Hype::init(void) {

  nhru = getdim(NHRU);
  for(hh = 0; hh < nhru; ++hh) {
    infil[hh] = 0.0;
    snowinfil[hh] = 0.0;
    runoff[hh] = 0.0;
    meltrunoff[hh] = 0.0;
    hru_evap[hh] = 0.0;
    hru_cum_evap[hh] = 0.0;
    hru_cum_actet[hh] = 0.0;
    hru_actet[hh] = 0.0;
    net_rain[hh] = 0.0;
    SWE[hh] = 0.0;
    SWE_max[hh] = 0.0;
    hru_t[hh] = 0.0;
    SWE_conc[hh] = 0.0;

    if(variation == VARIATION_1){
      soil_moist[hh] = 0.0;
      soil_rechr[hh] = 0.0;
    }

    for(long Sub = 0; Sub < numsubstances; ++Sub)
      SWE_conc_lay[Sub][hh] = 0.0;
  }
}

void ClassWQ_Test_Hype::run(void) {

  long step = getstep();
  long nstep = step%Global::Freq;

  for(hh = 0; chkStruct(); ++hh) {
    if(nstep == 1){
      if(Julian_start[hh] <= julian("now") && julian("now") < Julian_end[hh]){
        infil[hh] = infil_0[hh];
        snowinfil[hh] = snowinfil_0[hh];
        runoff[hh] = runoff_0[hh];
        meltrunoff[hh] = meltrunoff_0[hh];
        hru_t[hh] = hru_t_0[hh];
        SWE[hh] = SWE_0[hh];
        SWE_max[hh] = SWE[hh];

        for(long Sub = 0; Sub < numsubstances; ++Sub)
          SWE_conc_lay[Sub][hh] = SWE_conc_lay_0[Sub][hh];

        hru_evap[hh] = hru_evap_0[hh];
        hru_actet[hh] = hru_actet_0[hh];
        net_rain[hh] = net_rain_0[hh];
        if(variation == VARIATION_1){
          soil_moist[hh] = soil_moist_0[hh];
          soil_rechr[hh] = soil_moist_0[hh];
        }
      }
      else{
        infil[hh] = 0.0;
        snowinfil[hh] = 0.0;
        runoff[hh] = 0.0;
        meltrunoff[hh] = 0.0;
        hru_evap[hh] = 0.0;
        hru_cum_evap[hh] = 0.0;
        hru_cum_actet[hh] = 0.0;
        hru_actet[hh] = 0.0;
        net_rain[hh] = 0.0;
//        SWE[hh] = 0.0;
//        hru_t[hh] = 0.0;
//        SWE_conc[hh] = 0.0;
//        if(variation == VARIATION_1){
//          soil_moist[hh] = 0.0;
//          soil_rechr[hh] = 0.0;
//        }
      }
    }
  }
}

void ClassWQ_Test_Hype::finish(bool good){

  for(hh = 0; chkStruct(); ++hh) {
//    LogMessageA(hh, string("'" + Name + " (intcp)'  cumnetrain  (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumnet_rain[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");
  }
  LogDebug(" ");
}

ClassWQ_pbsmSnobal* ClassWQ_pbsmSnobal::klone(string name) const{
  return new ClassWQ_pbsmSnobal(name);
}

void ClassWQ_pbsmSnobal::decl(void) {

  Description = "'Special \"pbsm\" module compatible with \"snobal\".' \
                 'original version using hru_u,' \
                 'uses hru_Uadjust from walmsley_wind instead of hru_u,' \
                 'using hru_u and a regression to use daily windspeed,' \
                 'uses hru_Uadjust from walmsley_wind instead of hru_u and a regression to use daily windspeed.'";

  variation_set = VARIATION_0 + VARIATION_2;

  declgetvar("*", "hru_u", "(m/s)", &hru_u);


  variation_set = VARIATION_1 + VARIATION_3;

  declgetvar("*", "hru_Uadjust", "(m/s)", &hru_Uadjust);


  variation_set = VARIATION_2 + VARIATION_3;

  declparam("u_D", NHRU, "[1.0]", "0.5", "2.0", "Daily windspeed correction", "()", &u_D);

  declparam("Drift_offset", NHRU, "[0.0]", "-100.0", "100.0", "Daily windspeed drift offset correction", "()", &Drift_offset);

  declparam("Drift_slope", NHRU, "[1.0]", "0.5", "2.0", "Daily windspeed drift slope correction", "()", &Drift_slope);

  declparam("Subl_offset", NHRU, "[0.0]", "-100.0", "100.0", "Daily windspeed sublimation offset correction", "()", &Subl_offset);

  declparam("Subl_slope", NHRU, "[1.0]", "0.5", "2.0", "Daily windspeed sublimation slope correction", "()", &Subl_slope);


  variation_set = VARIATION_ORG;

  declstatvar("SWE_max", NHRU, "snow water equivalent seasonal maximum", "(mm)", &SWE_max);

  declvar("SWE_conc", NDEFN, "snow water equivalent", "(mg/l)", &SWE_conc, &SWE_conc_lay, numsubstances);

  declvar("Subl", NHRU, "interval sublimation", "(mm/int)", &Subl);

  declvar("Subl_conc", NDEFN, "interval sublimation", "(mm/int)", &Subl_conc, &Subl_conc_lay, numsubstances);

  declvar("cumDriftOut", NHRU, "cumulative transport from HRU", "(mm)", &cumDriftOut);

  declvar("cumDriftOut_mWQ", NDEFN, "mass solute from HRU", "(mg)", &cumDriftOut_mWQ, &cumDriftOut_mWQ_lay, numsubstances);

  declvar("Drift_out", NHRU, "interval transport out of HRU", "(mm/int)", &Drift_out);

  declvar("Drift_out_conc", NDEFN, "interval transport out of HRU", "(mg/l)", &Drift_out_conc, &Drift_out_conc_lay, numsubstances);

  declvar("hru_subl", NHRU, "interval sublimation", "(mm/int)", &Subl);

  declvar("hru_drift", NHRU, "interval composite transport", "(mm/int)", &Drift);

  declvar("Drift_in", NHRU, "interval transport into HRU", "(mm/int)", &Drift_in);

  declvar("Drift_in_conc", NDEFN, "interval transport into HRU", "(mg/l)", &Drift_in_conc, &Drift_in_conc_lay, numsubstances);

  decldiag("DriftH", NHRU, "interval transport", "(mm/int)", &DriftH);

  decldiag("SublH", NHRU, "interval sublimation", "(mm/int)", &SublH);

  decldiag("BasinSnowLoss", BASIN, "transport out of basin", "(mm/int)", &BasinSnowLoss);

  decldiag("BasinSnowLoss_mWQ", NDEF, "transport out of basin", "(mm/int)", &BasinSnowLoss_mWQ, &BasinSnowLoss_mWQ_lay, numsubstances);

  decldiag("BasinSnowGain", BASIN, "cumulative transport to basin estimated from HRU 1", "(mm/int)", &BasinSnowGain);

  decldiag("BasinSnowGain_mWQ", NDEF, "cumulative transport to basin estimated from HRU 1", "(mm/int)", &BasinSnowGain_mWQ, &BasinSnowGain_mWQ_lay, numsubstances);

  declvar("BasinSnowLoss", BASIN, "transport out of basin", "(mm/int)", &BasinSnowLoss);

  declvar("cumSubl", NHRU, "cumulative sublimation", "(mm)", &cumSubl);

  declvar("cumSubl_mWQ", NDEFN, "cumulative sublimation solute", "(mm)", &cumSubl_mWQ, &cumSubl_mWQ_lay, numsubstances);

  declvar("cumDrift", NHRU, "cumulative transport from HRU", "(mm)", &cumDrift);

  declvar("cumBasinSnowLoss", BASIN, "cumulative transport out of basin", "(mm)", &cumBasinSnowLoss);

  declvar("cumBasinSnowLoss_mWQ", NDEF, "cumulative mass of solute transport out of basin", "(mg)", &cumBasinSnowLoss_mWQ, &cumBasinSnowLoss_mWQ_lay, numsubstances);

  declvar("cumBasinSnowGain", BASIN, "cumulative transport to basin estimated from HRU 1", "(mm)", &cumBasinSnowGain);

  declvar("cumBasinSnowGain_mWQ", NDEF, "cumulative mass of solute transport to basin estimated from HRU 1", "(mg)", &cumBasinSnowGain_mWQ, &cumBasinSnowGain_mWQ_lay, numsubstances);

  declvar("cumDriftIn", NHRU, "cumulative transport to HRU", "(mm)", &cumDriftIn);

  declvar("cumDriftIn_mWQ", NDEFN, "cumulative mass of solute transport to HRU", "(mg)", &cumDriftIn_mWQ, &cumDriftIn_mWQ_lay, numsubstances);

  decllocal("hru_basin", NHRU, "conversion factor", "()", &hru_basin);

  decldiag("DrySnow", NHRU, "DrySnow", "()", &DrySnow);

  declvar("SnowAge", NHRU, "SnowAge", "()", &SnowAge);

  declvar("cumSno", NHRU, "cumulative snow", "(mm)", &cumSno);

  declvar("cumSno_mWQ", NDEFN, "cumulative mass of solute snow", "(mg)", &cumSno_mWQ, &cumSno_mWQ_lay, numsubstances);

  decldiag("Prob", NHRU, "Probability", "()", &Prob);

  decldiag("snowdepth", NHRU, "depth of snow using Gray/Pomeroy", "(m)", &snowdepth);

  decllocal("SWE_Init", NHRU, "initial SWE", "(mm)", &SWE_Init);

// parameters

  declparam("fetch", NHRU, "[1000.0]", "300.0", "10000.0", "fetch distance", "(m)", &fetch);

  declparam("Ht", NHRU, "[0.1, 0.25, 1.0]", "0.001", "100.0", "vegetation height(m)", "(m)", &Ht);

  declparam("distrib", NHRU, "[0.0, 1.0]", "-10.0", "10.0", "distribution fractions - can sum to 1", "()", &distrib);

  declparam("N_S", NHRU, "[320]", "1", "500", "vegetation number density", "(1/m^2)", &N_S);

  declparam("A_S", NHRU, "[0.003]", "0.0", "2.0", "stalk diameter", "(m)", &A_S);

  declparam("basin_area", BASIN, "3", "1e-6", "1e+09", "total basin area", "(km^2)", &basin_area);

  declparam("hru_area", NHRU, "[1]", "1e-6", "1e+09", "hru area", "(km^2)", &hru_area);

  declparam("inhibit_evap", NHRU, "[0]", "0", "1", "inhibit evaporatation, 1 -> inhibit", "()", &inhibit_evap);

  declparam("inhibit_bs", NHRU, "[0]", "0", "1", "inhibit blowing snow, 1 -> inhibit", "()", &inhibit_bs);

  declparam("inhibit_subl", NHRU, "[0]", "0", "1", "inhibit sublimation, 1 -> inhibit", "()", &inhibit_subl);

  declparam("rain_conc", NDEFN, "0", "0", "1000", "rain solute concentration", "(mg/l)", &rain_conc, &rain_conc_lay, numsubstances);

  declparam("snow_conc", NDEFN, "0", "0", "1000", "snow solute concentration", "(mg/l)", &snow_conc, &snow_conc_lay, numsubstances);

  declparam("Atmos_mWQ", NDEFN, "0", "0", "10", "total basin area", "(mg/int)", &Atmos_mWQ, &Atmos_mWQ_lay, numsubstances);


  decllocal("DrySnow_0", NHRU, "", "", &DrySnow_0);

  decllocal("SnowAge_0", NHRU, "", "", &SnowAge_0);

  decllocal("BasinSnowGain_0", NHRU, "", "", &BasinSnowGain_0);

  decllocal("cumBasinSnowGain_0", NHRU, "", "", &cumBasinSnowGain_0);

  decllocal("BasinSnowLoss_0", NHRU, "", "", &BasinSnowLoss_0);

  decllocal("cumBasinSnowLoss_0", NHRU, "", "", &cumBasinSnowLoss_0);

  decllocal("Subl_0", NHRU, "", "", &Subl_0);

  decllocal("Subl_0", NHRU, "", "", &Subl_0);

  decllocal("cumSubl_0", NHRU, "", "", &cumSubl_0);

  decllocal("Drift_in_0", NHRU, "", "", &Drift_in_0);

  decllocal("cumDriftIn_0", NHRU, "", "", &cumDriftIn_0);

  decllocal("Drift_out_0", NHRU, "", "", &Drift_out_0);

  decllocal("cumDriftOut_0", NHRU, "", "", &cumDriftOut_0);

  decllocal("SWE_0", NHRU, "", "", &SWE_0);

  decllocal("SWE_Init_0", NHRU, "", "", &SWE_Init_0);

  decllocal("cumSno_0", NHRU, "", "", &cumSno_0);

  declputvar("*", "SWE", "(kg/m^2)", &SWE);
  declgetvar("*", "z_s", "(m)", &z_s);
  declgetvar("*", "rho", "(kg/m^3)", &rho);

  declgetvar("*", "hru_t", "(�C)", &hru_t);
  declgetvar("*", "hru_ea", "(kPa)", &hru_ea);
  declgetvar("*", "hru_newsnow", "()", &hru_newsnow);
  declgetvar("*", "net_snow", "(mm/int)", &net_snow);

}

void ClassWQ_pbsmSnobal::init(void) {

  nhru = getdim(NHRU);

  Reset_WQ(0, cumBasinSnowLoss, cumBasinSnowLoss_mWQ_lay);
  Reset_WQ(0, cumBasinSnowGain, cumBasinSnowGain_mWQ_lay);
  Reset_WQ(0, BasinSnowLoss, BasinSnowLoss_mWQ_lay);
  Reset_WQ(0, BasinSnowGain, BasinSnowGain_mWQ_lay);

  for (hh = 0; hh < nhru; ++hh) {

    Reset_WQ(hh, Drift_in, Drift_in_conc_lay);
    Reset_WQ(hh, Drift_out, Drift_out_conc_lay);
    Reset_WQ(hh, cumDriftOut, cumDriftOut_mWQ_lay);
    Reset_WQ(hh, cumDriftIn, cumDriftIn_mWQ_lay);
    Reset_WQ(hh, cumSno, cumSno_mWQ_lay);
    Reset_WQ(hh, cumSubl, cumSubl_mWQ_lay);

    for(long Sub = 0; Sub < numsubstances; ++Sub)
      SWE_conc_lay[Sub][hh];


    SnowAge[hh] = 0.0;
    DrySnow[hh] = 0;
    snowdepth[hh] = 0.0;
    DriftH[hh] = 0.0;
    SublH[hh] = 0.0;
    Prob[hh] = 0.0;
    hru_basin[hh] = hru_area[hh]/basin_area[0];

    if((hh > 0) && (Ht[hh] < Ht[hh-1]) && distrib[hh-1] > 0){
      CRHMException TExcept(string("'" + Name + " (pbsmSnobal)' vegetation heights not in ascending order.").c_str(), WARNING);
      LogError(TExcept);
    }
  }
}

void ClassWQ_pbsmSnobal::run(void){

  float Znod, Ustar, Ustn, E_StubHt, Lambda, Ut, Uten_Prob;
  float SumDrift, SumDrift_mWQ, total, transport, transport_mWQ;

  for(long Sub = 0; Sub < numsubstances; ++Sub){
    if(getstep() == 1)
      for (hh = 0; chkStruct(); ++hh)
        SWE_Init[hh] = SWE[hh];

    if(Sub == 0) // saves all HRUs
      Save();

    for (hh = 0; chkStruct(); ++hh) {

      if(Sub != 0)
        Restore(hh);

     if(variation == VARIATION_ORG || variation == VARIATION_2)
       hru_u_ = hru_u[hh];
     else
       hru_u_ = hru_Uadjust[hh];

     if(variation == VARIATION_2 || variation == VARIATION_3)
       hru_u_ = u_D[hh]*hru_u_;

     Reset_WQ(hh, Drift_in, Drift_in_conc_lay);
     Reset_WQ(hh, Drift_out, Drift_out_conc_lay);
     Reset_WQ(hh, Subl, Subl_conc_lay);

     Drift[hh] = 0.0;
     DriftH[hh] = 0.0;
     SublH[hh] = 0.0;
     Prob[hh] = 0.0;

     if(SWE[hh] > 0.0)
       SWE_conc_lay[Sub][hh] = ((SWE[hh] - net_snow[hh])*SWE_conc_lay[Sub][hh] + snow_conc_lay[Sub][hh]*net_snow[hh] + Atmos_mWQ_lay[Sub][hh])/SWE[hh]; // netsnow already added in Snobal
     else
       SWE_conc_lay[Sub][hh] = 0.0;

     if(SWE[hh] > 0.0 && !inhibit_bs[hh]) {

       E_StubHt = Ht[hh] - z_s[hh];                      // depths(m)
       if(E_StubHt < 0.0001) E_StubHt = 0.0001;

       Ustar = 0.02264*pow(hru_u_, 1.295f);            // Eq. 6.2 rev.,  Ustar over fallow

       if (E_StubHt > 0.01) {
         Znod = (sqr(Ustar)/163.3f)+0.5*N_S[hh]*E_StubHt*A_S[hh]; // Eq. 29, Snowcover Book
         Lambda = N_S[hh]*A_S[hh]*E_StubHt;                      // Raupach Eq. 1

         Ustn  = Ustar*sqrt((PBSM_constants::Beta*Lambda)/(1.0+PBSM_constants::Beta*Lambda));

         Uten_Prob = (log(10.0/Znod))/PBSM_constants::KARMAN *sqrt(Ustar-Ustn);
       }
       else
       {
         Uten_Prob = hru_u_;
       } // end if

       bool newsnow = net_snow[hh]; // attributed with snow_conc_lay

       ProbabilityThresholdNew(SWE[hh], hru_t[hh], Uten_Prob, Prob[hh], Ut, newsnow, SnowAge[hh], DrySnow[hh]);  // (mm)

       if (Prob[hh] > 0.001) {
         Ut = Ut * 0.8;

         float RH = hru_ea[hh]/estar(hru_t[hh]); // Snobal uses Pascals

         Pbsm(E_StubHt, Ut, DriftH[hh], SublH[hh], hru_t[hh], hru_u_, RH, fetch[hh], N_S[hh], A_S[hh]);

         if(variation == VARIATION_2 || variation == VARIATION_3){
           DriftH[hh] = Drift_offset[hh] + DriftH[hh]*Drift_slope[hh];
           SublH[hh] = Subl_offset[hh] + SublH[hh]*Subl_slope[hh];
         }

         Drift_out[hh] = DriftH[hh]*Prob[hh]/fetch[hh];
         Drift_out_conc_lay[Sub][hh] = SWE_conc_lay[Sub][hh];

         if(!inhibit_subl[hh]){
           Subl[hh] = SublH[hh]*Prob[hh];
// assume WQ evaporates slso
         }

// handle insufficient snow pack

         if(Drift_out[hh] + Subl[hh] > SWE[hh]){
           if(Drift_out[hh] >= SWE[hh]){
             Drift_out[hh] = SWE[hh];
             Subl[hh] = 0.0;
           }
           else{
             Subl[hh] = SWE[hh] - Drift_out[hh];
// assume WQ evaporates slso
           }
         }
         cumDriftOut_mWQ_lay[Sub][hh] += Drift_out_conc_lay[Sub][hh]*Drift_out[hh];
         cumDriftOut[hh] += Drift_out[hh];
         cumSubl[hh] += Subl[hh];
         cumSubl_mWQ_lay[Sub][hh] += Subl[hh]; // stays behind or sublimates
         }
       } // end if
     } // end for (hh)

// distribute drift

    if(distrib[0] > 0.0) { // simulate transport entering basin using HRU 1
      float Drft = Drift_out[0]*distrib[0];
      SWE_conc_lay[Sub][0] = SWE_conc_lay[Sub][0]*SWE[0] + Drift_out_conc_lay[Sub][0]*Drft;
      SWE[0] += Drft;
      SWE_conc_lay[Sub][0] /= SWE[0];
      cumDriftIn[0] += Drft;
      cumDriftIn_mWQ_lay[Sub][0] += Drift_out_conc_lay[Sub][0]*Drft;
      cumBasinSnowGain[0] += Drft*hru_basin[0];  // **** hru_basin = hru_area/basin_area ****
    }

    BasinSnowLoss[0] = 0.0;
    long LastN = 0;

    if(!inhibit_bs[0]&& nhru == 1){
      BasinSnowLoss[0] = Drift_out[0];
      BasinSnowLoss_mWQ_lay[Sub][0] += Drift_out_conc_lay[Sub][0]*Drift_out[0];
      cumBasinSnowLoss[0] += BasinSnowLoss[0];
      cumBasinSnowLoss_mWQ_lay[Sub][0] += BasinSnowLoss_mWQ_lay[Sub][0];
    }

    for (long nn = LastN; chkStruct(nn); ++nn) {
      if(distrib[nn] >= 0.0 && nn+1 < nhru) // skip till last HRU or -ve distribution
        continue;

      SumDrift = 0.0; SumDrift_mWQ = 0.0;
      for (long hhh = LastN; chkStruct(hhh, nn); ++hhh){ // sum drift over range
          SumDrift += Drift_out[hhh]*hru_basin[hhh];
          SumDrift_mWQ += Drift_out[hhh]*hru_basin[hhh]*Drift_out_conc_lay[Sub][hhh];
      }

      if(SumDrift > 0.0){ // drift has occurred!
        for (long hh = LastN + 1; chkStruct(hh, nn+1); ++hh) {

          if(Ht[hh] > z_s[hh])
            SWE_max[hh] = SWE[hh] + rho[hh]*(Ht[hh]-z_s[hh]); // not filled
          else
            SWE_max[hh] = SWE[hh]; // filled or over filled. Wait for snow transport

          if(SWE_max[hh] <= 0.0)
             SWE_max[hh] = Ht[hh];

          if(hh == nn) { // handle last HRU
            if(distrib[nn] > 0){
              float In = SumDrift/hru_basin[hh]; // remaining drift
              if(SWE_max[hh] > SWE[hh] + In){ // fill snowpack, remainder leaves basin
                Drift_in[hh] = In; // can handle all
                Drift_in_conc_lay[Sub][hh] = SumDrift_mWQ/hru_basin[hh]/In;
                cumDriftIn[hh] += Drift_in[hh];
                cumDriftIn_mWQ_lay[Sub][hh] += SumDrift_mWQ/hru_basin[hh];
                transport = 0.0;
                transport_mWQ = 0.0;
              }
              else if(SWE_max[hh] > SWE[hh]){ // cannot handle all
                Drift_in[hh] = SWE_max[hh] - SWE[hh];
                float used = Drift_in[hh]/In;
                Drift_in_conc_lay[Sub][hh] = SumDrift_mWQ*used/Drift_in[hh];
                cumDriftIn[hh] += Drift_in[hh];
                cumDriftIn_mWQ_lay[Sub][hh] += SumDrift_mWQ*used;
                transport -= (In -(SWE_max[hh] - SWE[hh]))*hru_basin[hh];
                transport_mWQ -= SumDrift_mWQ*(1.0 - used)*hru_basin[hh];
              }
              else{ // zero or -ve - happens during melt??
                transport = SumDrift;
                transport_mWQ = SumDrift_mWQ;
              }
            }
            else if(distrib[nn] < 0){ // all drift deposited
                float In = SumDrift/hru_basin[hh]; // remaining drift
                Drift_in[hh] = SumDrift/hru_basin[hh]; // can handle all
                Drift_in_conc_lay[Sub][hh] = SumDrift_mWQ/hru_basin[hh]/In;
                cumDriftIn[hh] += Drift_in[hh];
                cumDriftIn_mWQ_lay[Sub][hh] += SumDrift_mWQ/hru_basin[hh];
                transport = 0.0;
                transport_mWQ = 0.0;
            }
            else{ // distrib[nn] == 0 -> all excess drift leaves basin
                transport = SumDrift;
                transport_mWQ = SumDrift_mWQ;
            }

            BasinSnowLoss[0] += (transport + Drift_out[hh]*hru_basin[hh]);
            BasinSnowLoss_mWQ_lay[Sub][0] += (transport_mWQ + Drift_out_conc_lay[Sub][hh]*Drift_out[hh]*hru_basin[hh]);
            cumBasinSnowLoss[0] += BasinSnowLoss[0];
            cumBasinSnowLoss_mWQ_lay[Sub][0] += BasinSnowLoss_mWQ_lay[Sub][0];
          }
          else if(SWE_max[hh] > SWE[hh] &&  distrib[hh] > 0.0) {
// handle intermediate HRUs with available storage and distrib > 0
            total = 0.0;
            for (long jj = hh; chkStruct(jj, nn+1); jj++) // calculate denominator
              total += fabs(distrib[jj]);
// determine contribution and scale
            transport = SumDrift*fabs(distrib[hh])/total/hru_basin[hh];
            transport_mWQ = SumDrift_mWQ*fabs(distrib[hh])/total/hru_basin[hh];
            if(SWE_max[hh] > SWE[hh] + transport){ // sufficient capacity
              Drift_in[hh] += transport;
            }
            else if(SWE_max[hh] > SWE[hh]){
              transport = SWE_max[hh] - SWE[hh];  // insufficient capacity
              Drift_in[hh] += transport;
            }
            else{
              transport = 0.0;
              transport_mWQ = 0.0;
            }

            SumDrift -= transport*hru_basin[hh]; // remove drift used from total available
            SumDrift_mWQ -= transport_mWQ*hru_basin[hh]; // remove drift used from total available
            cumDriftIn[hh] += transport;
            cumDriftIn_mWQ_lay[Sub][hh] += transport_mWQ;
          } // end if
        } // end for (hh)
        LastN = nn+1;
      } // end if
    } // end for (nn)

    for (hh = 0; chkStruct(); ++hh) { // snow cover inhibits evaporation
      Drift[hh] = Drift_in[hh] - Drift_out[hh];
      if(SWE[hh] > 0.0){
        const_cast<long*> (inhibit_evap)[hh] = 1;
        snowdepth[hh] = DepthofSnow(SWE[hh]);
      }
      else{
        const_cast<long*> (inhibit_evap)[hh] = 0;
        snowdepth[hh] = 0.0;
      }
    } // end for (hh)
  } // end for (Sub)
}

void ClassWQ_pbsmSnobal::finish(bool good) {

  if(!good) return;

  float AllcumSubl = 0.0;
  float AllcumCover = cumBasinSnowGain[0] - cumBasinSnowLoss[0];
  long Sub = 0;

  for(hh = 0; chkStruct(); ++hh) {
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumSno     (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSno[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumDriftOut(mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftOut[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumDriftIn (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftIn[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumSubl    (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSubl[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumCover   (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSno[hh]+cumDriftIn[hh]-cumDriftOut[hh]-cumSubl[hh], hru_area[hh], basin_area[0], "*** SWE just before melt");
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' SWE        (mm) (mm*hru) (mm*hru/basin): ").c_str(), SWE[hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' SWE_change (mm) (mm*hru) (mm*hru/basin): ").c_str(), SWE[hh] - SWE_Init[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");

    AllcumSubl += cumSubl[hh]*hru_area[hh];
    AllcumCover += (cumSno[hh]+cumDriftIn[hh]-cumDriftOut[hh]-cumSubl[hh])*hru_area[hh];
    LogDebug(" ");

    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumSno_mWQ     (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumSno_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumDriftOut_mWQ(mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftOut_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
    LogMessageA(hh, string("'" + Name + " (WQ_pbsmSnobal)' cumDriftIn_mWQ (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumDriftIn_mWQ_lay[Sub][hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");
  }

  LogMessage(string("'" + Name + " (WQ_pbsmSnobal)' AllcumSubl  (mm*basin): ").c_str(), AllcumSubl, "*** cumulative sum of all HRUs cumSubl");
  LogMessage(string("'" + Name + " (WQ_pbsmSnobal)' AllcumCover (mm*basin): ").c_str(), AllcumCover, "*** SWE just before melt cumulative sum of all HRUs cumCover");

  LogDebug(" ");
  LogMessage("'WQ_pbsm' cumBasinSnowLoss (mm): ", cumBasinSnowLoss[0]);
  LogMessage("'WQ_pbsm' cumBasinSnowGain (mm): ", cumBasinSnowGain[0]);

  LogDebug(" ");
  LogMessage("'WQ_pbsm' cumBasinSnowLoss_mWQ (substance) (mm): ", cumBasinSnowLoss_mWQ_lay[Sub][0]);
  LogMessage("'WQ_pbsm' cumBasinSnowGain_mWQ (substance) (mm): ", cumBasinSnowGain_mWQ_lay[Sub][0]);
  LogDebug(" ");
}

void ClassWQ_pbsmSnobal::Reset_WQ(long hru, float *var, float **var_WQ_lay){
  var[hru] = 0.0;
  for(long Sub = 0; Sub < numsubstances; ++Sub)
    var_WQ_lay[Sub][hru] = 0.0;
}

void ClassWQ_pbsmSnobal::copy_array(float *from, float *to){
  for(hh = 0; chkStruct(); ++hh)
    to[hh] = from[hh];
}

void ClassWQ_pbsmSnobal::copy_array(long *from, long *to){
  for(hh = 0; chkStruct(); ++hh)
    to[hh] = from[hh];
}

void ClassWQ_pbsmSnobal::copy_basin(float *from, float *to){
  to[0] = from[0];
}

void ClassWQ_pbsmSnobal::restore_hru(float *from, float *to, long hh){
  to[hh] = from[hh];
}

void ClassWQ_pbsmSnobal::restore_hru(long *from, long *to, long hh){
  to[hh] = from[hh];
}

void ClassWQ_pbsmSnobal::Restore(long hh){

  restore_hru(DrySnow_0, DrySnow, hh);
  restore_hru(SnowAge_0, SnowAge, hh);

  restore_hru(BasinSnowGain_0, BasinSnowGain, hh);
  restore_hru(cumBasinSnowGain_0, cumBasinSnowGain, hh);
  restore_hru(BasinSnowLoss_0, BasinSnowLoss, hh);
  restore_hru(cumBasinSnowLoss_0, cumBasinSnowLoss, hh);

  restore_hru(SWE_0, SWE, hh);
  restore_hru(SWE_Init_0, SWE_Init, hh);
  restore_hru(cumSno_0, cumSno, hh);
  restore_hru(Drift_in_0, Drift_in, hh);
  restore_hru(cumDriftIn_0, cumDriftIn, hh);
  restore_hru(Drift_out_0, Drift_out, hh); // not used in run
  restore_hru(cumDriftOut_0, cumDriftOut, hh);
  restore_hru(Subl_0, Subl, hh);
  restore_hru(cumSubl_0, cumSubl, hh);
}

void ClassWQ_pbsmSnobal::Save(){

  copy_array(DrySnow, DrySnow_0);
  copy_array(SnowAge, SnowAge_0);

  copy_basin(BasinSnowGain, BasinSnowGain_0);
  copy_basin(cumBasinSnowGain, cumBasinSnowGain_0);
  copy_basin(BasinSnowLoss, BasinSnowLoss_0);
  copy_basin(cumBasinSnowLoss, cumBasinSnowLoss_0);

  copy_array(SWE, SWE_0);
  copy_array(SWE_Init, SWE_Init_0);
  copy_array(cumSno, cumSno_0);
  copy_array(Drift_in, Drift_in_0);
  copy_array(cumDriftIn, cumDriftIn_0);
  copy_array(Drift_out, Drift_out_0);
  copy_array(cumDriftOut, cumDriftOut_0);
  copy_array(Subl, Subl_0);
  copy_array(cumSubl, cumSubl_0);
}

ClassWQ_mass_conc* ClassWQ_mass_conc::klone(string name) const{
  return new ClassWQ_mass_conc(name);
}

void ClassWQ_mass_conc::decl(void) {

  Description = "'converts _mWQ to concentrations'";

  declvar("inflow_conc", NDEFN, "Mass: inflow from other HRUs", "(mg)", &inflow_conc, &inflow_conc_lay, numsubstances);
  declvar("outflow_conc", NDEFN, "Mass: inflow from other HRUs", "(mg)", &outflow_conc, &outflow_conc_lay, numsubstances);
  declvar("runflow_conc", NDEFN, "Mass: inflow from other HRUs", "(mg)", &runoutflow_conc, &runoutflow_conc_lay, numsubstances);
  declvar("ssrflow_conc", NDEFN, "Mass: inflow from other HRUs", "(mg)", &ssroutflow_conc, &ssroutflow_conc_lay, numsubstances);
  declvar("gwoutflow_conc", NDEFN, "Mass: inflow from other HRUs", "(mg)", &gwoutflow_conc, &gwoutflow_conc_lay, numsubstances);


  declgetvar("*", "inflow", "(mm/int)", &inflow);
  declgetvar("*", "inflow_mWQ", "(mg)", &inflow_mWQ, &inflow_mWQ_lay);
  declgetvar("*", "outflow", "(mm/int)", &outflow);
  declgetvar("*", "outflow_mWQ", "(mg)", &outflow_mWQ, &outflow_mWQ_lay);
  declgetvar("*", "runoutflow", "(mm/int)", &runoutflow);
  declgetvar("*", "runoutflow_mWQ", "(mg)", &runoutflow_mWQ, &runoutflow_mWQ_lay);
  declgetvar("*", "ssroutflow", "(mm/int)", &ssroutflow);
  declgetvar("*", "ssroutflow_mWQ", "(mg)", &ssroutflow_mWQ, &ssroutflow_mWQ_lay);
  declgetvar("*", "gwoutflow", "(mm/int)", &gwoutflow);
  declgetvar("*", "gwoutflow_mWQ", "(mg)", &gwoutflow_mWQ, &outflow_mWQ_lay);

}

void ClassWQ_mass_conc::init(void) {

  nhru = getdim(NHRU);
  Reset_WQ(inflow_conc_lay);
  Reset_WQ(outflow_conc_lay);
  Reset_WQ(runoutflow_conc_lay);
  Reset_WQ(ssroutflow_conc_lay);
  Reset_WQ(gwoutflow_conc_lay);
}

void ClassWQ_mass_conc::run(void) {

  mass_to_conc(inflow, inflow_mWQ_lay, inflow_conc_lay);
  mass_to_conc(outflow, outflow_mWQ_lay, outflow_conc_lay);
  mass_to_conc(runoutflow, runoutflow_mWQ_lay, runoutflow_conc_lay);
  mass_to_conc(ssroutflow, ssroutflow_mWQ_lay, ssroutflow_conc_lay);
  mass_to_conc(gwoutflow, gwoutflow_mWQ_lay, gwoutflow_conc_lay);
}

void ClassWQ_mass_conc::mass_to_conc(const float *var, const float **var_mWQ, float **var_conc) {

  for(long Sub = 0; Sub < numsubstances; ++Sub){
    for(long hh = 0; hh < nhru; ++hh) {
      if(var[hh] <= 0.0)
         var_conc[Sub][hh] = 0.0;
      else
        var_conc[Sub][hh] = var_mWQ[Sub][hh]/var[hh];
    }
  }
}

void ClassWQ_mass_conc::Reset_WQ(float **var_lay){
  for(long hh = 0; hh < nhru; ++hh) {
    for(long Sub = 0; Sub < numsubstances; ++Sub)
      var_lay[Sub][hh] = 0.0;
  }
}

ClassWQ_Substitute_Hype* ClassWQ_Substitute_Hype::klone(string name) const{
  return new ClassWQ_Substitute_Hype(name);
}

void ClassWQ_Substitute_Hype::decl(void) {

  Description = "'Provides necessary Hype variables conc_top, conc_bottom and conc_below to run other WQ modules to run without proper Hype module.'";

  declstatvar("conc_top", NDEFN, "concentration of inorganic nitrogen in soil moisture per land-soil", "(mg/l)", &conc_top, &conc_top_lay, numsubstances); //

  declstatvar("conc_bottom", NDEFN, "concentration of organic nitrogen in soil moisture per land-soil", "(mg/l)", &conc_bottom, &conc_bottom_lay, numsubstances);

  declstatvar("conc_below", NDEFN, "concentration of soluble (reactive) phosphorus, i.e. phosphate in soil moisture per land-soil", "(mg/l)", &conc_below, &conc_below_lay, numsubstances);

}

void ClassWQ_Substitute_Hype::init(void) {

  nhru = getdim(NHRU);
  for(hh = 0; hh < nhru; ++hh) {
    for(long Sub = 0; Sub < numsubstances; ++Sub){
      conc_top_lay[Sub][hh] = 0.0;
      conc_bottom_lay[Sub][hh] = 0.0;
      conc_below_lay[Sub][hh] = 0.0;
    }
  }
}

void ClassWQ_Substitute_Hype::run(void) {

  long step = getstep();
  long nstep = step%Global::Freq;

  for(hh = 0; chkStruct(); ++hh) {
  }
}

void ClassWQ_Substitute_Hype::finish(bool good){

  for(hh = 0; chkStruct(); ++hh) {
//    LogMessageA(hh, string("'" + Name + " (intcp)'  cumnetrain  (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumnet_rain[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");
  }
  LogDebug(" ");
}

ClassWQ_Gen_Mass_Bal* ClassWQ_Gen_Mass_Bal::klone(string name) const{
  return new ClassWQ_Gen_Mass_Bal(name);
}

void ClassWQ_Gen_Mass_Bal::decl(void) {

  Description = "'Generates project mass balance. Parameter defines how to handle 'cum*** variables''";

  declvar("WQ_Total_mass", BASIN, "Total_mass of liquid water", "(mm)", &WQ_Total_mass);

  declvar("WQ_Total_mass_in", BASIN, "Total_mass in", "(mm)", &WQ_Total_mass_in);

  declvar("WQ_Total_mass_out", BASIN, "Total_mass out", "(mm)", &WQ_Total_mass_out);

  declvar("WQ_Total_mass_atmos", BASIN, "Total_mass interaction with the atmosphere", "(mm)", &WQ_Total_mass_atmos);

  declvar("WQ_hru_Total_mass", NHRU, "Total_mass", "(mm)", &WQ_hru_Total_mass);

  declvar("WQ_hru_Total_mass_in", NHRU, "Total_mass in", "(mm)", &WQ_hru_Total_mass_in);

  declvar("WQ_hru_Total_mass_out", NHRU, "Total_mass out", "(mm)", &WQ_hru_Total_mass_out);

  declvar("WQ_hru_Total_mass_atmos", NHRU, "Total_mass interaction with the atmosphere", "(mm)", &WQ_hru_Total_mass_atmos);


  declparam("basin_area", BASIN, "3", "1e-6", "1e+09", "total basin area", "(km^2)", &basin_area);

  declparam("hru_area", NHRU, "[1]", "1e-6", "1e+09", "hru area", "(km^2)", &hru_area);

  Atmos_Vars_WQ_pbsm = declGMBparam("Atmos_Vars_WQ_pbsm", TWO, "'WQ_pbsm', 'cumSubl_mWQ'", "Variables used in mass balance.",  Atmos_Vars_WQ_pbsm);
  In_Vars_WQ_pbsm = declGMBparam("In_Vars_WQ_pbsm", THREE, "'WQ_pbsm', 'cumSno_mWQ', 'cumDriftIn_mWQ'", "Variables used in mass balance.",  In_Vars_WQ_pbsm);
  Out_Vars_WQ_pbsm = declGMBparam("Out_Vars_WQ_pbsm", TWO, "'WQ_pbsm', 'cumDriftOut_mWQ'", "Variables used in mass balance.",  Out_Vars_WQ_pbsm);

  Atmos_Vars_WQ_pbsmSnowbal = declGMBparam("Atmos_Vars_WQ_pbsmSnowbal", TWO, "'WQ_pbsmSnobal', 'cumSubl_mWQ'", "Variables used in mass balance.",  Atmos_Vars_WQ_pbsmSnowbal);
  In_Vars_WQ_pbsmSnowbal = declGMBparam("In_Vars_WQ_pbsmSnowbal", THREE, "'WQ_pbsmSnobal', 'cumSno_mWQ', 'cumDriftIn_mWQ'", "Variables used in mass balance.",  In_Vars_WQ_pbsmSnowbal);
  Out_Vars_WQ_pbsmSnowbal = declGMBparam("Out_Vars_WQ_pbsmSnowbal", TWO, "'WQ_pbsmSnobal', 'cumDriftOut_mWQ'", "Variables used in mass balance.",  Out_Vars_WQ_pbsmSnowbal);

  Out_Vars_WQ_Netroute_M_D = declGMBparam("Out_Vars_WQ_Netroute_M_D", THREE, "'WQ_Netroute_M_D', 'cumoutflow_mWQ', 'cumoutflow_diverted_mWQ'", "Variables used in mass balance.",  Out_Vars_WQ_Netroute_M_D);

  In_Vars_WQ_Soil = declGMBparam("In_Vars_WQ_Soil", SIX, "'WQ_Soil', 'soil_moist_change_mWQ', 'soil_rechr_change_mWQ', 'soil_bottom_change_mWQ', 'Sd_change_mWQ', 'gw_change_mWQ'", "change in reservoirs values.", In_Vars_WQ_Soil);
  Out_Vars_WQ_Soil = declGMBparam("Out_Vars_WQ_Soil", SIX, "'WQ_Soil', 'cum_rechr_ssr_mWQ', 'cum_soil_runoff_mWQ', 'cum_runoff_to_Sd_mWQ', 'cum_gw_flow_mWQ', 'cum_redirected_residual_mWQ'", "flows values.", Out_Vars_WQ_Soil);
  Atmos_Vars_WQ_Soil = declGMBparam("Atmos_Vars_WQ_Soil", TWO, "'WQ_Soil', 'cum_Sd_evap'", "Sd evaporation.",  Atmos_Vars_WQ_Soil);
}

void ClassWQ_Gen_Mass_Bal::init(void) {

  nhru = getdim(NHRU);

  WQ_Total_mass[0] = 0;
  WQ_Total_mass_in[0] = 0;
  WQ_Total_mass_out[0] = 0;
  WQ_Total_mass_atmos[0] = 0;

  for(hh = 0; chkStruct(); ++hh) {
    WQ_hru_Total_mass = 0;
    WQ_hru_Total_mass_in[hh] = 0;
    WQ_hru_Total_mass_out[hh] = 0;
    WQ_hru_Total_mass_atmos[hh] = 0;
  }
}

void ClassWQ_Gen_Mass_Bal::run(void) {
  WQ_Total_mass_in[0] = 0;
}

void ClassWQ_Gen_Mass_Bal::finish(bool good) {

  LogDebug("**WQ_Gen_Mass_Bal**");
  LogDebug(" ");

  if(In_Vars_WQ_Soil){
    Add_mass(In_Vars_WQ_Soil);
  }

  if(Out_Vars_WQ_Soil){
    Add_mass(Out_Vars_WQ_Soil);
  }

  if(Atmos_Vars_WQ_Soil){
    Add_mass(Atmos_Vars_WQ_Soil);
  }

  if(In_Vars_WQ_pbsm){
    Add_mass(In_Vars_WQ_pbsm);
  }

  if(In_Vars_WQ_pbsmSnowbal){
    Add_mass(In_Vars_WQ_pbsmSnowbal);
  }

  if(Out_Vars_WQ_Netroute_M_D){
    Transfer_mass(Out_Vars_WQ_Netroute_M_D);
  }

  if(Out_Vars_WQ_pbsm){
    Transfer_mass(Out_Vars_WQ_pbsm);
  }

  if(Out_Vars_WQ_pbsmSnowbal){
    Transfer_mass(Out_Vars_WQ_pbsmSnowbal);
  }

  if(Atmos_Vars_WQ_pbsm){
    mass_to_Atmosphere(Atmos_Vars_WQ_pbsm);
  }

  if(Atmos_Vars_WQ_pbsmSnowbal){
    mass_to_Atmosphere(Atmos_Vars_WQ_pbsmSnowbal);
  }
}

void ClassWQ_Gen_Mass_Bal::Add_mass(TStringList* WQ_In_Vars_obs){

  string Module, Variable, Items;
  float Value;

  for(hh = 0; chkStruct(); ++hh) {
    Module = (WQ_In_Vars_obs->Strings[0]).c_str();
    long Cnt = WQ_In_Vars_obs->Count;
    WQ_hru_Total_mass_in[hh] = 0.0;

    for(long vv = 1; vv < Cnt; ++vv){
      Variable = (WQ_In_Vars_obs->Strings[vv]).c_str();
      Value = ((ClassVar*) WQ_In_Vars_obs->Objects[vv])->values[hh]*hru_area[hh];
      WQ_hru_Total_mass_in[hh] += Value;
      LogMessageA(hh, string("'" + Name + " (WQ_hru_Total_mass_in)' " + Module + " " + Variable + "        (mm) (mm*hru) (mm*hru/basin): ").c_str(), Value, hru_area[hh], basin_area[0]);
      if(hh == 0){
        if(vv > 1)
          Items += (" + ");
        Items += (Variable);
      }
    }
    WQ_Total_mass_in[0] +=  WQ_hru_Total_mass_in[hh];
    WQ_Total_mass[0] +=  WQ_hru_Total_mass_in[hh];
  }

  LogDebug(" ");
  LogMessage(string(string("Add WQ_Total_mass_in     ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass_in[0]);
  LogMessage(string(string("Add WQ_Total_mass(basin) ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass[0]);
  LogDebug(" ");
}

void ClassWQ_Gen_Mass_Bal::Subtract_mass(TStringList* WQ_In_Vars_obs){
  string Variable, Module, Items;
  float Value;

  for(hh = 0; chkStruct(); ++hh) {
    Module = (WQ_In_Vars_obs->Strings[0]).c_str();
    long Cnt = WQ_In_Vars_obs->Count;
    WQ_hru_Total_mass_out[hh] = 0.0;

    for(long vv = 1; vv < Cnt; ++vv){
      Variable = (WQ_In_Vars_obs->Strings[vv]).c_str();
      Value = ((ClassVar*) WQ_In_Vars_obs->Objects[vv])->values[hh]*hru_area[hh];
      WQ_hru_Total_mass_out[hh] += Value;
      LogMessageA(hh, string("'" + Name + " (WQ_hru_Total_mass_in)' " + Module + " " + Variable + "        (mm) (mm*hru) (mm*hru/basin): ").c_str(), Value, hru_area[hh], basin_area[0]);
      if(hh == 0){
        if(vv > 1)
          Items += (" + ");
        Items += (Variable);
      }
    }
    WQ_Total_mass_out[0] +=  WQ_hru_Total_mass_out[hh];
    WQ_Total_mass[0] -=  WQ_hru_Total_mass_out[hh];
  }

  LogDebug(" ");
  LogMessage(string(string("Subtract WQ_Total_mass_out    ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass_out[0]);
  LogMessage(string(string("Subtract WQ_Total_mass(basin) ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass[0]);
  LogDebug(" ");
}

void ClassWQ_Gen_Mass_Bal::mass_to_Atmosphere(TStringList* WQ_In_Vars_obs){

  string Variable, Module, Items;
  float Value;

  for(hh = 0; chkStruct(); ++hh) {
    Module = (WQ_In_Vars_obs->Strings[0]).c_str();
    long Cnt = WQ_In_Vars_obs->Count;
    WQ_hru_Total_mass_out[hh] = 0.0;

    for(long vv = 1; vv < Cnt; ++vv){
      Variable = (WQ_In_Vars_obs->Strings[vv]).c_str();
      Value = ((ClassVar*) WQ_In_Vars_obs->Objects[vv])->values[hh]*hru_area[hh];
      WQ_hru_Total_mass_out[hh] += Value;
      LogMessageA(hh, string("'" + Name + " (WQ_hru_Total_mass_in)' " + Module + " " + Variable + "        (mm) (mm*hru) (mm*hru/basin): ").c_str(), Value, hru_area[hh], basin_area[0]);
      if(hh == 0){
        if(vv > 1)
          Items += (" + ");
        Items += (Variable);
      }
    }
    WQ_Total_mass_out[0] +=  WQ_hru_Total_mass_in[hh];
  }

  WQ_Total_mass[0] -=  WQ_Total_mass_out[0];

  LogDebug(" ");
  LogMessage(string(string("Atmosphere WQ_Total_mass_out    ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass_out[0]);
  LogMessage(string(string("Atmosphere WQ_Total_mass(basin) ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass[0]);
  LogDebug(" ");
}

void ClassWQ_Gen_Mass_Bal::Add_mass_first(TStringList* WQ_In_Vars_obs){

  string Module, Variable, Items;
  float Value;

  for(hh = 0; chkStruct(); ++hh) {
    Module = (WQ_In_Vars_obs->Strings[0]).c_str();
    WQ_hru_Total_mass_in[hh] = 0.0;

    long vv = 1;
    Variable = (WQ_In_Vars_obs->Strings[vv]).c_str();
    Value = ((ClassVar*) WQ_In_Vars_obs->Objects[vv])->values[hh]*hru_area[hh];
    WQ_hru_Total_mass_in[hh] += Value;
    LogMessageA(hh, string("'" + Name + " (WQ_hru_Total_mass_in)' " + Module + " " + Variable + "        (mm) (mm*hru) (mm*hru/basin): ").c_str(), Value, hru_area[hh], basin_area[0]);
    if(hh == 0){
      if(vv > 1)
        Items += (" + ");
      Items += (Variable);
    }
    WQ_Total_mass_in[0] +=  WQ_hru_Total_mass_in[hh];
    WQ_Total_mass[0] +=  WQ_hru_Total_mass_in[hh];
  }

  LogDebug(" ");
  LogMessage(string(string("Add first_WQ_Total_mass_in     ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass_in[0]);
  LogMessage(string(string("Add first_WQ_Total_mass(basin) ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass[0]);
  LogDebug(" ");
}

void ClassWQ_Gen_Mass_Bal::Subtract_mass_last(TStringList* WQ_In_Vars_obs){
  string Variable, Module, Items;
  float Value;

  long Cnt = WQ_In_Vars_obs->Count;
  long vv = Cnt-1;

  for(hh = 0; chkStruct(); ++hh) {
    WQ_hru_Total_mass_out[hh] = 0.0;
    Module = (WQ_In_Vars_obs->Strings[0]).c_str();
    WQ_hru_Total_mass_out[hh] = 0.0;

    Variable = (WQ_In_Vars_obs->Strings[vv]).c_str();
    Value = ((ClassVar*) WQ_In_Vars_obs->Objects[vv])->values[hh]*hru_area[hh];
    WQ_hru_Total_mass_out[hh] -= Value;
    LogMessageA(hh, string("'" + Name + " (WQ_hru_Total_mass_in)' " + Module + " " + Variable + "        (mm) (mm*hru) (mm*hru/basin): ").c_str(), Value, hru_area[hh], basin_area[0]);
    if(hh == 0){
      if(vv > 1)
        Items += (" + ");
      Items += (Variable);
    }
    WQ_Total_mass_out[0] +=  WQ_hru_Total_mass_out[hh];
    WQ_Total_mass[0] -=  WQ_hru_Total_mass_out[hh];
  }

  LogDebug(" ");
  LogMessage(string(string("Subtract last_WQ_Total_mass_out    ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass_out[0]);
  LogMessage(string(string("Subtract WQ_Total_mass(basin)      ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass[0]);
  LogDebug(" ");
}

void ClassWQ_Gen_Mass_Bal::Transfer_mass(TStringList* WQ_In_Vars_obs){

  string Variable, Module, Items;
  float Value;

  for(hh = 0; chkStruct(); ++hh) {
    Module = (WQ_In_Vars_obs->Strings[0]).c_str();
    long Cnt = WQ_In_Vars_obs->Count;
    WQ_hru_Total_mass_out[hh] = 0.0;

    for(long vv = 1; vv < Cnt; ++vv){
      Variable = (WQ_In_Vars_obs->Strings[vv]).c_str();
      Value = ((ClassVar*) WQ_In_Vars_obs->Objects[vv])->values[hh]*hru_area[hh];
      WQ_hru_Total_mass_out[hh] += Value;
      LogMessageA(hh, string("'" + Name + " (WQ_hru_Total_mass_in)' " + Module + " " + Variable + "        (mm) (mm*hru) (mm*hru/basin): ").c_str(), Value, hru_area[hh], basin_area[0]);
      if(hh == 0){
        if(vv > 1)
          Items += (" + ");
        Items += (Variable);
      }
    }
    WQ_Total_mass_out[0] +=  WQ_hru_Total_mass_out[hh];
  }

  WQ_Total_mass[0] -=  WQ_Total_mass_out[0];

  LogDebug(" ");
  LogMessage(string(string("Transfer WQ_Total_mass_out    ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass_out[0]);
  LogMessage(string(string("Transfer WQ_Total_mass(basin) ") + WQ_In_Vars_obs->Strings[0].c_str() + " " + Items.c_str() + "  (mm*basin): ").c_str(), WQ_Total_mass[0]);
  LogDebug(" ");
}

ClassGrow_crops_annually* ClassGrow_crops_annually::klone(string name) const{
  return new ClassGrow_crops_annually(name);
}

void ClassGrow_crops_annually::decl(void) {

  Description = "'Misc modules.' \
                 'nothing,' \
                 'Hype fertilizer.'\
                 'Grow_crops.'\
                 'Both.'";

// Hype parameters not affected fertdown1, fertdowd2, mandown1, mandown2, part, resfast and resdown

  variation_set = VARIATION_1 + VARIATION_3;

  declvar("Fert_N_amount1", NHRU, "current hru N fertilizer amount", "(kg/km^2)", &Fert_N_amount1);

  declvar("Fert_P_amount1", NHRU, "current hru P fertilizer amount", "(kg/km^2)", &Fert_P_amount1);

  declvar("Man_N_amount1", NHRU, "current hru N manure amount", "(kg/km^2)", &Man_N_amount1);

  declvar("Man_P_amount1", NHRU, "current hru P manure amount", "(kg/km^2)", &Man_P_amount1);

  declvar("Res_N_amount", NHRU, "current hru N residue amount", "(kg/km^2)", &Res_N_amount);

  declvar("Res_P_amount", NHRU, "current hru P residue amount", "(kg/km^2)", &Res_P_amount);

  declvar("Fert_N_amount2", NHRU, "current hru N fertilizer amount", "(kg/km^2)", &Fert_N_amount2);

  declvar("Fert_P_amount2", NHRU, "current hru P fertilizer amount", "(kg/km^2)", &Fert_P_amount2);

  declvar("Man_N_amount2", NHRU, "current hru N manure amount", "(kg/km^2)", &Man_N_amount2);

  declvar("Man_P_amount2", NHRU, "current hru P manure amount", "(kg/km^2)", &Man_P_amount2);

  declvar("Fertperiod", NHRU, "current period for spreading fertilizer and manure amount", "()", &Fertperiod);

  declvar("Litterperiod", NHRU, "current period for residue amount", "()", &Litterperiod);

  declvar("Fertday", NHRU, "day to apply fertilizer", "()", &Fertday);

  declvar("Manday", NHRU, "day to apply manure", "()", &Manday);

  declvar("Resdayno", NHRU, "day to apply residue", "()", &Resdayno);

  declvar("LockOut", NHRU, "prevents changes of fertilizer and manure amounts on multi day periods", "()", &LockOut);

  declvar("SecondDown_fert", NHRU, "second down", "()", &SecondDown_fert);

  declvar("SecondDown_man", NHRU, "second down", "()", &SecondDown_man);

  ObsCnt_N = declreadobs("Fert_N", NHRU, "annual N fertilizer dates and amount", "(kg/km^2)", &Fert_N, HRU_OBS_Q, true);

  ObsCnt_P = declreadobs("Fert_P", NHRU, "annual P fertilizer dates and amount", "(kg/km^2)", &Fert_P, HRU_OBS_Q, true);

  ObsCntMan_N = declreadobs("Man_N", NHRU, "annual N Manure dates and amount", "(kg/km^2)", &Man_N, HRU_OBS_Q, true);

  ObsCntMan_P = declreadobs("Man_P", NHRU, "annual P Manure dates and amount", "(kg/km^2)", &Man_P, HRU_OBS_Q, true);

  ObsCntRes_N = declreadobs("Res_N", NHRU, "annual N Residues dates and amount", "(kg/km^2)", &Res_N, HRU_OBS_Q, true);

  ObsCntRes_P = declreadobs("Res_P", NHRU, "annual P Residues dates and amount", "(kg/km^2)", &Res_P, HRU_OBS_Q, true);

  ObsCnt_fertperiod = declreadobs("Fert_period", NHRU, "spreading period for feritilzer and manure", "(d)", &Fert_period, HRU_OBS_Q, true);


  declparam("Ag_YearStart", NHRU, "[0]", "0", "10", " suggestions for northern hemisphere - 0, southern hemisphere - 183", "()", &Ag_YearStart);

  declputparam("*", "fertNamount1", "(kg/km^2)", &fertNamount1);

  declputparam("*", "fertNamount2", "(kg/km^2)", &fertNamount2);

  declputparam("*", "fertPamount1", "(kg/km^2)", &fertPamount1);

  declputparam("*", "fertPamount2", "(kg/km^2)", &fertPamount2);

  declputparam("*", "manNamount1", "(kg/km^2)", &manNamount1);

  declputparam("*", "manNamount2", "(kg/km^2)", &manNamount2);

  declputparam("*", "manPamount1", "(kg/km^2)", &manPamount1);

  declputparam("*", "manPamount2", "(kg/km^2)", &manPamount2);

  declputparam("*", "resNamount", "(kg/km^2)", &resNamount);

  declputparam("*", "resPamount", "(kg/km^2)", &resPamount);

  declputparam("*", "fertday1", "()", &fertday1);

  declputparam("*", "fertday2", "()", &fertday2);

  declputparam("*", "manday1", "()", &manday1);

  declputparam("*", "manday2", "()", &manday2);

  declputparam("*", "resdayno", "()", &resdayno);

  declputparam("*", "fertperiod", "()", &fertperiod);

  declputparam("*", "litterperiod", "()", &litterperiod);

  variation_set = VARIATION_2 + VARIATION_3;

  declparam("Htmax", NHRU, "[0.1, 0.25, 1.0]", "0.001", "100.0", "maximum vegetation height", "(m)", &Htmax);

  declparam("Init_Crop_Ht_1", NHRU, "[0.1]", "0.001", "100.0", "initial crop height (1)", "(m)", &Init_Crop_Ht_1);

  declparam("Crop_Grow_Rate_1", NHRU, "[0.8]", "0.0", "1.0", "crop growth rate (1)", "(m/d)", &Crop_Grow_Rate_1);

  declparam("JCrop_Start_1", NHRU, "[250]", "0", "366", "start Julian day (1); JCrop_Start_1 = 0 if no crop", "()", &JCrop_Start_1);

  declparam("JCrop_Harvest_1", NHRU, "[228]", "0", "366", "harvest Julian day (1); JCrop_Harvest_1 = 0 if no crop", "()", &JCrop_Harvest_1);

  declparam("Crop_Htmax_1", NHRU, "[0.1, 0.25, 1.0]", "0.001", "100.0", "maximum vegetation height (1)", "(m)", &Crop_Htmax_1);

  declparam("Init_Crop_Ht_2", NHRU, "[0.1]", "0.001", "100.0", "initial crop height (2)", "(m)", &Init_Crop_Ht_2);

  declparam("Crop_Grow_Rate_2", NHRU, "[0.8]", "0.0", "1.0", "crop growth rate (2)", "(m/d)", &Crop_Grow_Rate_2);

  declparam("JCrop_Start_2", NHRU, "[250]", "0", "366", "start Julian day (2); JCrop_Start_2 = 0 if no crop", "()", &JCrop_Start_2);

  declparam("JCrop_Harvest_2", NHRU, "[228]" , "0", "366", "harvest Julian day (2); JCrop_Harvest_2 = 0 if no crop", "()", &JCrop_Harvest_2);

  declparam("Crop_Htmax_2", NHRU, "[0.1, 0.25, 1.0]", "0.001", "100.0", "maximum vegetation height (2)", "(m)", &Crop_Htmax_2);

  declputparam("*", "Ht", "(m)", &Ht);

  variation_set = VARIATION_ORG;
}

void ClassGrow_crops_annually::init(void) {

  nhru = getdim(NHRU);

  if(variation == VARIATION_1 || VARIATION_3){
    if(ObsCnt_N > -1){
      CRHMException TExcept("Handling N fertilizer from (Fert_N) observation.", WARNING);
      LogError(TExcept);
    }

    if(ObsCnt_P > -1){
      CRHMException TExcept("Handling P fertilizer from (Fert_P) observation.", WARNING);
      LogError(TExcept);
    }

    if(ObsCntMan_N > -1){
      CRHMException TExcept("Handling N manure from (Man_N) observation.", WARNING);
      LogError(TExcept);
    }

    if(ObsCntMan_P > -1){
      CRHMException TExcept("Handling P manure from (Man_P) observation.", WARNING);
      LogError(TExcept);
    }

    if(ObsCntRes_N > -1){
      CRHMException TExcept("Handling N residues from (Res_N) observation.", WARNING);
      LogError(TExcept);
    }

    if(ObsCntRes_P > -1){
      CRHMException TExcept("Handling P residues from (Res_P) observation.", WARNING);
      LogError(TExcept);
    }

    if(ObsCnt_fertperiod > -1){
      CRHMException TExcept("Handling fertilizer and manure period (Fertperiod) observation.", WARNING);
      LogError(TExcept);
    }

    for(hh = 0; hh < nhru; ++hh) {
       const_cast<float *> (Ht)[hh] = Init_Crop_Ht_1[hh];
       Fert_N_amount1[hh] = 0.0;
       Fert_P_amount1[hh] = 0.0;
       Man_N_amount1[hh] = 0.0;
       Man_P_amount1[hh] = 0.0;
       Fert_N_amount2[hh] = 0.0;
       Fert_P_amount2[hh] = 0.0;
       Man_N_amount2[hh] = 0.0;
       Man_P_amount2[hh] = 0.0;
       Res_N_amount[hh] = 0.0;
       Res_P_amount[hh] = 0.0;
       Fertday[hh] = 0;
       Manday[hh] = 0;
       Resdayno[hh] = 0;
       Fertperiod[hh] = 1;
       Litterperiod[hh] = 1;
//       const_cast<float *> (Fert_period)[hh] = 0.0;
       LockOut[hh] = 0;
    }
  } // VARIATION_1

  if(variation == VARIATION_2 || VARIATION_3){

    if(Good_Dates(JCrop_Start_1)){
      CRHMException TExcept("JCrop_Start_1 dates out of range!", TERMINATE);
      LogError(TExcept);
    }

    if(Good_Dates(JCrop_Harvest_1)){
      CRHMException TExcept("JCrop_Harvest_1 dates out of range!", TERMINATE);
      LogError(TExcept);
    }

    if(Good_Dates(JCrop_Start_2)){
      CRHMException TExcept("JCrop_Start_2 dates out of range!", TERMINATE);
      LogError(TExcept);
    }

    if(Good_Dates(JCrop_Harvest_2)){
      CRHMException TExcept("JCrop_Harvest_2 dates out of range!", TERMINATE);
      LogError(TExcept);
    }
    for(hh = 0; hh < nhru; ++hh)
       const_cast<float *> (Ht)[hh] = Init_Crop_Ht_1[hh];
  } // VARIATION_2
}

void ClassGrow_crops_annually::run(void) {

  long step = getstep();
  long nstep = step%Global::Freq;
  long today = julian("now");
  bool Good_N, Good_P;

  if(step == 1){ // first step of run
    for(hh = 0; chkStruct(); ++hh){
      if(variation == VARIATION_1 || VARIATION_3){
        if(ObsCnt_N >= hh){ // file open
          declputparam("*", "fertNamount1", "(kg/km^2)", &fertNamount1);
          const_cast<float *> (fertNamount1)[hh] = 0.0; // set by module
        }
        if(ObsCnt_N >= hh){ // file open
          declputparam("*", "fertNamount2", "(kg/km^2)", &fertNamount2);
          const_cast<float *> (fertNamount2)[hh] = 0.0; // not used
        }
        if(ObsCnt_P >= hh){ // file open
          declputparam("*", "fertPamount1", "(kg/km^2)", &fertPamount1);
          const_cast<float *> (fertPamount1)[hh] = 0.0; // set by module
        }
        if(ObsCnt_P >= hh){ // file open
          declputparam("*", "fertPamount2", "(kg/km^2)", &fertPamount2);
          const_cast<float *> (fertPamount2)[hh] = 0.0; // not used
        }
        if(ObsCntMan_N >= hh){ // file open
          declputparam("*", "manNamount1", "(kg/km^2)", &manNamount1);
          const_cast<float *> (manNamount1)[hh] = 0.0; // set by module
        }
        if(ObsCntMan_P >= hh){ // file open
          declputparam("*", "manPamount1", "(kg/km^2)", &manPamount1);
          const_cast<float *> (manPamount1)[hh] = 0.0; // set by module
        }
        if(ObsCntMan_N >= hh){ // file open
          declputparam("*", "manNamount2", "(kg/km^2)", &manNamount2);
          const_cast<float *> (manNamount2)[hh] = 0.0; // not used
        }
        if(ObsCntMan_P >= hh){ // file open
          declputparam("*", "manPamount2", "(kg/km^2)", &manPamount2);
          const_cast<float *> (manPamount2)[hh] = 0.0; // not used
        }
        if(ObsCntRes_N >= hh){ // file open
          declputparam("*", "resNamount", "(kg/km^2)", &resNamount);
          const_cast<float *> (resNamount)[hh] = 0.0;  // set by module
        }
        if(ObsCntRes_P >= hh){ // file open
          declputparam("*", "resPamount", "(kg/km^2)", &resPamount);
          const_cast<float *> (resPamount)[hh] = 0.0;  // set by module
        }
        if(ObsCnt_N >= hh || ObsCnt_P >= hh){ // file open
          declputparam("*", "fertday1", "(d)", &fertday1);
          const_cast<float *> (fertday1)[hh] = 0; // set by module

          declputparam("*", "fertday2", "(d)", &fertday2);
          const_cast<float *> (fertday2)[hh] = 0; // not used
        }
        if(ObsCntMan_N >= hh || ObsCntMan_P >= hh){ // file open
          declputparam("*", "manday1", "(d)", &manday1);
          const_cast<float *> (manday1)[hh] = 0; // set by module

          declputparam("*", "manday2", "(d)", &manday2);
          const_cast<float *> (manday2)[hh] = 0; // not used
        }
        if(ObsCnt_fertperiod >= hh){ // file open
          declputparam("*", "fertperiod", "(d)", &fertperiod);
          const_cast<float *> (fertperiod)[hh] = 0; // set by module
        }
//        if(ObsCnt_litterperiod >= hh){ // file open  check ????
          declputparam("*", "litterperiod", "(d)", &litterperiod);
//          const_cast<float *> (litterperiod)[hh] = 0; // set by module
//        }
      } // VARIATION_1 or VARIATION_3
    } // for hh
  } // first step of run

  if(nstep == 1){ // beginning of every day
    for(hh = 0; chkStruct(); ++hh){
      if(variation == VARIATION_1 || VARIATION_3){
        if(Ag_YearStart[hh] == today){
          SecondDown_fert[hh] = 0;
          SecondDown_man[hh] = 0;
          fertday1[hh] = 0;
          fertday2[hh] = 0;
          manday1[hh] = 0;
          manday2[hh] = 0;
          resdayno[hh] = 0;
        }
        --LockOut[hh];
        if(LockOut[hh] <= 0){
          if(ObsCnt_fertperiod > 0){ // file open
            if(!(Fert_period[hh] >= fLimit || !Fert_period[hh])){
               fertperiod[hh] = 1;
               if(Fert_period[hh] > 1.0 && Fert_period[hh] < 366){
                fertperiod[hh] = Fert_period[hh];
               }
               Fertperiod[hh] = Fert_period[hh];
              }
          } // file open

          Good_N = TRUE, Good_P = TRUE;
          if(ObsCnt_N >= hh || ObsCnt_P >= hh){ // either file open
            if(ObsCnt_N < 0 || Fert_N[hh] >= fLimit || !Fert_N[hh] || Fert_N[hh] <= 0.0){
              Fert_N_amount1[hh] = 0.0; // set to all zeros
              Fert_N_amount2[hh] = 0.0; // set to all zeros
              fertNamount1[hh] = 0.0; // set to all zeros
              fertNamount2[hh] = 0.0; // set to all zeros
              Good_N = FALSE;
            }
            if(ObsCnt_P < 0 || Fert_P[hh] >= fLimit || !Fert_P[hh] || Fert_P[hh] <= 0.0){
              Fert_P_amount1[hh] = 0.0; // set to all zeros
              Fert_P_amount2[hh] = 0.0; // set to all zeros
              fertPamount1[hh] = 0.0; // set to all zeros
              fertPamount2[hh] = 0.0; // set to all zeros
              Good_P = FALSE;
            }
            if(Good_N || Good_P){
                LockOut[hh] = fertperiod[hh];
                if(!SecondDown_fert[hh]){
                  if(Good_N && Fert_N[hh] > 0.0){
                    fertNamount1[hh] = Fert_N[hh]/Global::Freq; // set to interval values given
                    Fert_N_amount1[hh] = fertNamount1[hh];
                  }
                  if(Good_P && Fert_P[hh] > 0.0){
                    fertPamount1[hh] = Fert_P[hh]/Global::Freq; // set to interval values given
                    Fert_P_amount1[hh] = fertPamount1[hh];
                  }
                  fertday1[hh] = today;
                  Fertday[hh] = today;
                  SecondDown_fert[hh] = 1;
                }
                else{
                  fertNamount1[hh] = 0.0;
                  fertPamount1[hh] = 0.0;
                  Fert_N_amount1[hh] = 0.0; // set to all zeros
                  Fert_P_amount1[hh] = 0.0; // set to all zeros
                  if(Good_N && Fert_N[hh] > 0.0){
                    fertNamount2[hh] = Fert_N[hh]/Global::Freq; // set to interval values given
                    Fert_N_amount2[hh] = fertNamount2[hh];
                  }
                  if(Good_P && Fert_P[hh] > 0.0){
                    fertPamount2[hh] = Fert_P[hh]/Global::Freq; // set to interval values given
                    Fert_P_amount2[hh] = fertPamount2[hh];
                  }
                  fertday1[hh] = 0;
                  fertday2[hh] = today;
                  Fertday[hh] = today;
                  SecondDown_fert[hh] = 0;
                }
              } // fertN or P > 0.0
              else{
                fertday1[hh] = 0;
                fertday2[hh] = 0;
                Fertday[hh] = 0;
              }
          } // file open

          Good_N = TRUE, Good_P = TRUE;
          if(ObsCntMan_N >= hh || ObsCntMan_P >= hh){ // file open
            if(ObsCntMan_N < 0 || Man_N[hh] >= fLimit || !Man_N[hh] || Man_N[hh] <= 0.0){
              Man_N_amount1[hh] = 0.0; // set to all zeros
              Man_N_amount2[hh] = 0.0; // set to all zeros
              manNamount1[hh] = 0.0; // set to all zeros
              manNamount2[hh] = 0.0; // set to all zeros
              Good_N = FALSE;
            }
            if(ObsCntMan_P < 0 || Man_P[hh] >= fLimit || !Man_P[hh] || Man_P[hh] <= 0.0){
              Man_P_amount1[hh] = 0.0; // set to all zeros
              Man_P_amount2[hh] = 0.0; // set to all zeros
              manPamount1[hh] = 0.0; // set to all zeros
              manPamount2[hh] = 0.0; // set to all zeros
              Good_P = FALSE;
            }
            if(Good_N || Good_P){
              LockOut[hh] = fertperiod[hh];
              if(!SecondDown_man[hh]){
                if(Good_N && Man_N[hh] > 0.0){
                  manNamount1[hh] = Man_N[hh]/Global::Freq; // set to interval values given
                  Man_N_amount1[hh] = manNamount1[hh];
                }
                if(Good_P && Man_P[hh] > 0.0){
                  manPamount1[hh] = Man_P[hh]/Global::Freq; // set to interval values given
                  Man_P_amount1[hh] = manPamount1[hh];
                }
                  manday1[hh] = today;
                  Manday[hh] = today;
                  SecondDown_man[hh] = 1;
                }
                else{
                  manNamount1[hh] = 0.0; // set to all zeros
                  manPamount1[hh] = 0.0; // set to all zeros
                  Man_N_amount1[hh] = 0.0; // set to all zeros
                  Man_P_amount1[hh] = 0.0; // set to all zeros
                  if(Good_N && Man_N[hh]){
                    manNamount2[hh] = Man_N[hh]/Global::Freq; // set to interval values given
                    Man_N_amount2[hh] = manNamount2[hh];
                  }
                  if(Good_P && Man_P[hh]){
                    manPamount2[hh] = Man_P[hh]/Global::Freq; // set to interval values given
                    Man_P_amount2[hh] = manPamount2[hh];
                  }
                  manday1[hh] = 0;
                  manday2[hh] = today;
                  Manday[hh] = today;
                  SecondDown_man[hh] = 0;
               }
              } // manNamount1 > 0.0
              else{
                manday1[hh] = 0;
                manday2[hh] = 0;
                Manday[hh] = 0;
              }
            } // else (read amount from open file)


          Good_N = TRUE, Good_P = TRUE;
          if(ObsCntRes_N >= hh && ObsCntRes_P >= hh){ // file open
            if(ObsCntRes_N < 0 || Res_N[hh] >= fLimit || !Res_N[hh] || Res_N[hh] <= 0.0){
              Res_N_amount[hh] = 0.0; // set to all zeros
              resNamount[hh] = 0.0;
              Good_N = False;
            }
            if(ObsCntRes_P < 0 || Res_P[hh] >= fLimit || !Res_P[hh] || Res_P[hh] <= 0.0){
              Res_P_amount[hh] = 0.0; // set to all zeros
              resPamount[hh] = 0.0;
              Good_P = False;
            }
            if(Good_N || Good_P){
              LockOut[hh] = fertperiod[hh];
              if(Good_N && Man_N[hh]){
                resNamount[hh] = Res_N[hh]/Global::Freq; // set to interval values given
                Res_N_amount[hh] = Res_N[hh]; // set to interval values given
              }
              if(Good_P && Man_P[hh]){
                resPamount[hh] = Res_P[hh]/Global::Freq; // set to interval values given
                Res_P_amount[hh] = Res_P[hh]; // set to interval values given
              }
              resdayno[hh] = today;
              Resdayno[hh] = today;
            } // resNamount1 > 0.0
            else{
              resdayno[hh] = 0;
              Resdayno[hh] = 0;
            }
          } // file open
        } // lockout
      } // VARIATION_1 or VARIATION_3

      if(variation == VARIATION_2 || VARIATION_3){
        if(today == JCrop_Start_1[hh])
          const_cast<float *> (Ht)[hh] = Init_Crop_Ht_1[hh];

        if(JCrop_Start_1[hh] != 0 && today >= JCrop_Start_1[hh] && today <= JCrop_Harvest_1[hh])
          const_cast<float *> (Htmax)[hh] = Crop_Htmax_1[hh];

        if(Ht[hh] < Crop_Htmax_1[hh])
          const_cast<float *> (Ht)[hh] =  Ht[hh] + Crop_Grow_Rate_1[hh];
        else
          const_cast<float *> (Ht)[hh] = Htmax[hh];

        if(today == JCrop_Harvest_1[hh])
          const_cast<float *> (Ht)[hh] =  Init_Crop_Ht_1[hh];

        if(today == JCrop_Start_2[hh])
          const_cast<float *> (Ht)[hh] = Init_Crop_Ht_2[hh];

        if(Ht[hh] < Crop_Htmax_2[hh])
          const_cast<float *> (Ht)[hh] =  Ht[hh] + Crop_Grow_Rate_2[hh];
        else
          const_cast<float *> (Ht)[hh] = Htmax[hh];

        if(today == JCrop_Harvest_2[hh])
          const_cast<float *> (Ht)[hh] =  Init_Crop_Ht_2[hh];
      } // VARIATION_2 or VARIATION_3
    } // for hh
  } // start of day
}

void ClassGrow_crops_annually::finish(bool good) {

  for(hh = 0; chkStruct(); ++hh) {
//    LogMessageA(hh, string("'" + Name + " (Grow_crops_annually)'  cumnetrain  (mm) (mm*hru) (mm*hru/basin): ").c_str(), cumnet_rain[hh], hru_area[hh], basin_area[0]);
    LogDebug(" ");
  }
  LogDebug(" ");
}

bool ClassGrow_crops_annually::Good_Dates(const float* dates) {

  for(hh = 0; hh < nhru; ++hh) {
    if(dates[hh] > 366 || dates[hh] < 0);
    return FALSE; // error
  }
  return TRUE;
}

