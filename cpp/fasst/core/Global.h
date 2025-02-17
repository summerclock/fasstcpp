#pragma once

#include <array>
#include <memory>
#include <string>
#include <cmath>
#include <limits>
#include "fasst/utils/Array.h"
#include "fasst/fasst_export.h"

namespace fasst {

// System limits
struct SystemLimits
{
  static constexpr int kMaxLayers = 10;          // maxl: maximum number of layers
  static constexpr int kMaxNodes = 100;          // maxn: maximum number of nodes
  static constexpr int kMaxParams = 25;          // maxp: maximum number of soil parameters
  static constexpr int kMaxColumns = 35;         // maxcol: max columns of met data
  static constexpr int kMaxOverlap = 15;         // moverlap: max lines of metfile overlap
  static constexpr int kCanopyLayers = 3;        // nclayers: number of canopy layers
  static constexpr double kMissingValue = 999.0; //: spflag: missing soil, veg parameter flag
};

// Physical constants
struct PhysicalConstants
{
  static constexpr double kPi = 3.141592654;
  static constexpr double kSigma = 5.669e-8;                                 // Stefan-Boltzmann const (W/m^2*K^4)
  static constexpr double kKelvinRef = 273.15;                                    // Reference temperature (K)
  static constexpr double kGravity = 9.81;                                   // grav: gravitational acceleration (m/s^2)
  static constexpr double kVonKarman = 0.4;                                  // vK: von Karman constant
  static constexpr double kEpsilon = std::numeric_limits<double>::epsilon(); // eps: tolerance limit for equality
  static constexpr double kGasConstWater = 8.3143 / 1.8015e-2;               // Rv: gas constant for water vapor (J/kg*K)
  static constexpr double kGasConstAir = 8.3143 / 2.8964e-2;                 // Rd: gas constant for dry air (J/kg*K)
  static constexpr double kVegThermalCond = 0.38;                            // kveg: vegetation thermal conductivity (W/m*K)
};

// Snow parameters
struct SnowParams
{
  static constexpr double kAlbedoNew = 0.8;         // snalbedo: new snow albedo
  static constexpr double kAlbedoOld = 0.5;         // soalbedo: old snow albedo
  static constexpr double kEmissivity = 0.92;       // semis: snow emissivity
  static constexpr double kThermalDiff = 2e-7;      // sthdiff: thermal diffusivity (m^2/s)
  static constexpr double kThermalCond = 0.3492;    // sthcond: thermal conductivity (W/m*K)
  static constexpr double kDensityWet = 550.0;      // sdensw: wet snow density (kg/m^3)
  static constexpr double kDensityDry = 50.0;       // sdensd: dry snow density (kg/m^3)
  static constexpr double kLatentHeatFus = 3.335e5; // lhfus: latent heat of fusion (J/kg)
  static constexpr double kLatentHeatSub = 2.838e6; // lhsub: latent heat of sublimation (J/kg)
};

// Ice parameters
struct IceParams
{
  static constexpr double kAlbedo = 0.7;           // ialbedo: ice albedo
  static constexpr double kEmissivity = 0.9;       // iemis: ice emissivity
  static constexpr double kThermalDiff = 1.167e-6; // ithdiff: thermal diffusivity (m^2/s)
  static constexpr double kThermalCond = 2.1648;   // ithcond: thermal conductivity (W/m*K)
  static constexpr double kDensity = 916.5;        // idens: ice density (kg/m^3)
  static constexpr double kLimitIce = 0.998;       // ilim: closeness of ice to max water content
  static constexpr double kLimitWater = 0.998;     // wlim: closeness of water vapor to porosity
};

// Met file column pointers
struct MetFileColumns
{
  static constexpr int kYear = 1;      // ip_year
  static constexpr int kDayOfYear = 2; // ip_doy
  static constexpr int kHour = 3;      // ip_hr
  static constexpr int kMinute = 4;    // ip_min
  static constexpr int kPressure = 5;  // ip_ap: air pressure (mbar)
  static constexpr int kTemp = 6;      // ip_tmp: air temperature (C)
  static constexpr int kHumidity = 7;  // ip_rh: relative humidity (%)
  static constexpr int kWindSpeed = 8; // ip_ws: wind speed (m/s)
  static constexpr int kWindDir = 9;   // ip_wdir: wind direction
  static constexpr int kPrecip = 10;   // ip_prec: precipitation rate (mm/hr)
  // ...existing columns...
  static constexpr int kVisibility = 34; // ip_vis: visibility (km)
  static constexpr int kAerosol = 35;    // ip_aer: aerosol type
};

class DLL_PUBLIC GlobalData
{
public:
  GlobalData();
  ~GlobalData();

  // Met data and parameters
  struct MetData
  {
    int max_lines{ 0 };                                                                                                                                                 // max_lines
    int num_columns{ 0 };                                                                                                                                               // ncols
    int end_index{ 0 };                                                                                                                                                 // iend
    int measurement_count{ 0 };                                                                                                                                         // met_count
    int single_multi_flag{ 0 };                                                                                                                                         // single_multi_flag
    int water_flag{ 0 };                                                                                                                                                // water_flag
    double missing_data_flag{ 0.0 };                                                                                                                                    // m_flag
    double time_offset{ 0.0 };                                                                                                                                          // time_offset
    double elevation{ 0.0 };                                                                                                                                            // elevation
    double site_elevation{ 0.0 };                                                                                                                                       // elev
    double latitude{ 0.0 };                                                                                                                                             // lat
    double longitude{ 0.0 };                                                                                                                                            // mlong
    double timestep{ 0.0 };                                                                                                                                             // timestep
    double slope{ 0.0 };                                                                                                                                                // slope
    double albedo{ 0.0 };                                                                                                                                               // albedo
    double emissivity{ 0.0 };                                                                                                                                           // emis
    double instrument_height{ 0.0 };                                                                                                                                    // iheight
    double aspect{ 0.0 };                                                                                                                                               // aspect
    double slope_radians{ 0.0 };                                                                                                                                        // sloper
    utils::Array<double, SystemLimits::kMaxColumns, SystemLimits::kMaxLayers> data = utils::Array<double, SystemLimits::kMaxColumns, SystemLimits::kMaxLayers>::Zero(); // met data
  } met_data;

  // Soil and node data
  struct SoilData
  {
    int num_layers{ 0 };                                      // nlayers
    int num_nodes{ 0 };                                       // nnodes
    int ref_node{ 0 };                                        // refn
    utils::ArrayX<int, SystemLimits::kMaxLayers> soil_type{};   // soil_type (icourse)
    utils::ArrayX<int, SystemLimits::kMaxNodes> node_type{};    // node_type
    utils::ArrayX<int, SystemLimits::kMaxNodes> coarse_index{}; // coarse_index

    // Temperature and depth arrays
    utils::Array<double, SystemLimits::kMaxNodes, 1> temperature = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();   // tm
    utils::Array<double, SystemLimits::kMaxNodes, 1> depth_temprature = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero(); // zti
    utils::Array<double, SystemLimits::kMaxNodes, 1> depth_moisture = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();         // zm
    utils::Array<double, SystemLimits::kMaxNodes, 1> moisture = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();      // sm
    utils::Array<double, SystemLimits::kMaxNodes, 1> soil_moisture = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero(); // soil_moist

    // Soil properties arrays
    utils::Array<double, SystemLimits::kMaxNodes, 1> thermal_conductivity = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero(); // grthcond
    utils::Array<double, SystemLimits::kMaxNodes, 1> specific_heat = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();        // grspheat
    utils::Array<double, SystemLimits::kMaxNodes, 1> node_spacing = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();         // nz
    utils::Array<double, SystemLimits::kMaxNodes, 1> delta_z = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();              // delzs

    // Extended arrays (nodes+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> temperature_state = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero(); // stt
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> ice_content = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();       // ice
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> pressure_head = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();     // phead
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> vapor_content = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();     // wvc

    // Layer properties
    utils::Array<double, SystemLimits::kMaxLayers, 1> layer_thickness = utils::Array<double, SystemLimits::kMaxLayers, 1>::Zero(); // lthick
    utils::Array<double, SystemLimits::kMaxLayers, 1> density_factor = utils::Array<double, SystemLimits::kMaxLayers, 1>::Zero();  // rho_fac

    // Soil parameters matrix
    utils::Array<double, SystemLimits::kMaxLayers, SystemLimits::kMaxParams> soil_params_initial = utils::Array<double, SystemLimits::kMaxLayers, SystemLimits::kMaxParams>::Zero(); // isoilp
    utils::Array<double, SystemLimits::kMaxLayers, SystemLimits::kMaxParams> soil_params = utils::Array<double, SystemLimits::kMaxLayers, SystemLimits::kMaxParams>::Zero();         // soilp
    utils::Array<double, SystemLimits::kMaxNodes, SystemLimits::kMaxParams> node_soil_params = utils::Array<double, SystemLimits::kMaxNodes, SystemLimits::kMaxParams>::Zero();      // nsoilp

    // Additional properties
    double ground_albedo{ 0.0 };     // sgralbedo
    double ground_emissivity{ 0.0 }; // sgremis
    double groundwater_level{ 0.0 }; // gwl
    double roughness{ 0.0 };         // rough

    // Node properties
    utils::Array<double, SystemLimits::kMaxNodes, 1> node_pos_initial = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();     // nzi
    utils::Array<double, SystemLimits::kMaxNodes, 1> node_spacing_initial = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero(); // delzsi

    // Type descriptions
    utils::ArrayX<std::string, SystemLimits::kMaxLayers> soil_type_desc{};     // stype
    utils::ArrayX<std::string, SystemLimits::kMaxLayers> soil_class_desc{};    // sclass
    utils::ArrayX<std::string, SystemLimits::kMaxNodes + 2> node_class_desc{}; // nclass
  } soil_data;

  // Snow and ice data
  struct SnowIceData
  {
    utils::ArrayX<int, 2> sn_istat{};                                                    // TODO: sn_istat meaning? snow_ice_status?
    double init_snow_depth{ 0.0 };                                                     //: hsaccum
    double init_ice_thickness{ 0.0 };                                                  //: hi
    double snow_depth{ 0.0 };                                                          //: dsnow
    double init_snow_water_equivalent{ 0.0 };                                          //: iswe
    double refreeze_amount{ 0.0 };                                                     // refreeze
    double new_snow_deposition{ 0.0 };                                                 //: newsd
    double atopf{ 0.0 };                                                               // TODO: atopf meaning? atop freezing point?
    double thermal_conductivity{ 0.0 };                                                //: km
    double specific_heat{ 0.0 };                                                       //: sphm
    double ice_refreezing{ 0.0 };                                                      //: refreezei
    utils::Array<double, 15, 1> snow_statistics = utils::Array<double, 15, 1>::Zero(); //: sn_stat
    double snow_melt_volume{ 0.0 };                                                    //: vsmelt
    double ice_melt_volume{ 0.0 };                                                     //: vimelt
    double mixed_depth{ 0.0 };                                                         //: hm
  } snow_ice_data;

  // Canopy/vegetation data
  struct VegetationData
  {
    // Layer arrays
    utils::Array<double, SystemLimits::kCanopyLayers, 1> leaf_area_index = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero(); // ilai/laif
    utils::Array<double, SystemLimits::kCanopyLayers, 1> clumping = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();        // iclump
    utils::Array<double, SystemLimits::kCanopyLayers, 1> reflectance = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();     // irho
    utils::Array<double, SystemLimits::kCanopyLayers, 1> transmittance = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();   // itau
    utils::Array<double, SystemLimits::kCanopyLayers, 1> alpha = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();           // ialp
    utils::Array<double, SystemLimits::kCanopyLayers, 1> emissivity = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();      // ieps
    utils::Array<double, SystemLimits::kCanopyLayers, 1> dz = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();              // dzveg/idzveg
    utils::Array<double, SystemLimits::kCanopyLayers, 1> storage_liquid = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();  // storcl
    utils::Array<double, SystemLimits::kCanopyLayers, 1> storage_snow = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();    // storcs
    utils::Array<double, SystemLimits::kCanopyLayers, 1> foliage_type = utils::Array<double, SystemLimits::kCanopyLayers, 1>::Zero();    // ifoliage_type
    utils::Array<double, 5, SystemLimits::kCanopyLayers> avect = utils::Array<double, 5, SystemLimits::kCanopyLayers>::Zero();           // avect

    // Vegetation parameters matrix
    utils::Array<double, 18, 17> veg_properties = utils::Array<double, 18, 17>::Zero(); // veg_prp

    // Root parameters
    utils::Array<double, 18, SystemLimits::kMaxNodes> root_k = utils::Array<double, 18, SystemLimits::kMaxNodes>::Zero();           // rk
    utils::Array<double, SystemLimits::kMaxNodes, 1> root_fraction_low = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();  // frl
    utils::Array<double, SystemLimits::kMaxNodes, 1> root_fraction_high = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero(); // frh
    utils::Array<double, SystemLimits::kMaxNodes, 1> root_sink = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();          // sinkr

    // Vegetation types and flags
    int low_veg_type{ 0 };      //: vegl_type
    int veg_type_high{ 0 };     //: vegh_type
    int low_veg_flag{ 0 };      //: veg_flagl
    int veg_flag_high{ 0 };     //: veg_flagh
    int season{ 0 };            // iseason
    int biome_source{ 0 };      //: biome_source
    int new_veg_type_low{ 0 };  //: new_vtl
    int new_veg_type_high{ 0 }; //: new_vth

    // Physical parameters
    double low_foliage_density{ 0.0 };      //: isigfl, initial foliage density
    double low_foliage_emissivity{ 0.0 };   //: iepf, initial foliage emissivity
    double low_foliage_absorptivity{ 0.0 }; //: ifola, initial foliage absorptivity (1 - albedo)
    double low_foliage_height{ 0.0 };       //: ihfol, initial foliage height(cm)
    double high_foliage_density{ 0.0 };     //: isigfh, canopy density
    double high_foliage_height{ 0.0 };      //: izh/zh, canopy height
    double height_n{ 0.0 };                 // iheightn
    double zero_plane{ 0.0 };               // Zd
    double albedo_f{ 0.0 };                 // albf
    double foliage_temp{ 0.0 };             // ftemp
    double roughness_length{ 0.0 };         // z0l

    // Additional parameters
    double leaf_area_index_low{ 0.0 }; // ilail/lail
    double sigma_f_low_val{ 0.0 };     // sigfl
    double state{ 0.0 };
    double foliage_area_val{ 0.0 };   // fola
    double emissivity_f_val{ 0.0 };   // epf
    double height_foliage_val{ 0.0 }; // hfol
    double storage_liquid_low{ 0.0 }; // stll
    double storage_snow_low{ 0.0 };   // stls
  } veg_data;

  // Profile data using fixed-size matrices
  struct ProfileData
  {
    int total_nodes{ 0 };                                                                                                                    // ntot
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> conductivity_upper = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();  // khu(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> conductivity_lower = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();  // khl(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> velocity = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();            // vin(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> flow_upper = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();          // flowu(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> flow_lower = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();          // flowl(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> flux = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();                // fv1(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> source = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();              // source(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> sink = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();                // sink(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> moisture_head_deriv = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero(); // dsmdh(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> temp_old = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();            // too(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> moisture_old = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();        // smoo(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> vapor_old = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();           // woo(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> ice_old = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();             // ioo(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes + 2, 1> pressure_head_old = utils::Array<double, SystemLimits::kMaxNodes + 2, 1>::Zero();   // phoo(maxn+2)
    utils::Array<double, SystemLimits::kMaxNodes, 1> back_fill_temp = utils::Array<double, SystemLimits::kMaxNodes, 1>::Zero();              // bftm(maxn)
  } profile_data;

private:
  // void Initialize();
};

template<typename T>
DLL_PUBLIC inline bool ApproxEqual(const T a, const T b)
{
  return std::abs(static_cast<int>((a - b) * 1e5)) / 1e5 <= PhysicalConstants::kEpsilon;
}

template<typename T>
DLL_PUBLIC inline bool ApproxZero(const T value)
{
  return std::abs(value) <= PhysicalConstants::kEpsilon;
}

template<typename T>
DLL_PUBLIC inline T Clamp(const T value, const T min, const T max)
{
  return std::max(min, std::min(value, max));
}

// Global singleton instance
extern DLL_PUBLIC std::unique_ptr<GlobalData> g_data;

} // namespace fasst
