#pragma once 

namespace fasst {
namespace utils {

// Dense function - calculates density based on temperature
double CalcDensity(double temp1, double wind1, int flag);

// Derivative of density with respect to temperature  
double CalcDensityDerivative(double temp1, int flag);

// Calculate specific heat
double CalcSpecificHeat(double temp1, int flag);

// Calculate thermal conductivity
double CalcThermalConductivity(double temp1, int flag);

// Calculate met date
double CalcMetDate(double year, double doy, double hr, double minute);

// Calculate pressure head
double CalcHead(int i, double smt);

// Calculate soil humidity
double CalcSoilHumidity(int i, double ph, double sms, double st);

// Calculate vapor pressure 
double CalcVaporPressure(int i, double rh, double ap);

// Calculate max infiltration
double CalcMaxInfiltration(int i, double smi);

// Map USDA soil type to USCS
int MapUSDAToUSCS(int lis); 

// USCS sand/silt/clay classification
int USCSSandSiltClay(double sand, double silt, double clay,
                     double p200, double plimit, double llimit,
                     double orgf);

} // namespace utils
} // namespace fasst
