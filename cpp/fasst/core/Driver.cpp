#include "Driver.h"

#include <fstream>
#include <chrono>
#include <iostream>
#include <cmath>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <glog/logging.h>

#include "fasst/core/Global.h"
#include "fasst/utils/Array.h"
#include "fasst/utils/Functions.h"

namespace fasst {

namespace fs = boost::filesystem;

auto& met = g_data->met_data;
auto& snow_ice = g_data->snow_ice_data;
auto& soil = g_data->soil_data;
auto& veg = g_data->veg_data;
auto& profile = g_data->profile_data;

Driver::Driver() = default;

void Driver::Run(const std::string& input_file)
{
  auto start = std::chrono::steady_clock::now();

  ReadControlFile(input_file);
  ReadMeteorologicalData();
  ReadOldData();
  InitializeSoilProfile();
  OutputData();

  CleanupFiles();

  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  LOG(INFO) << "Total run time: " << duration.count() << "s\n";
}

void Driver::ReadControlFile(const std::string& input_file)
{
  // fasst_driver.F90:197

  std::ifstream input(input_file);
  if (!input) {
    LOG(FATAL) << "Error opening FASST input file: " << input_file;
  }

  LOG(INFO) << "Reading control file: " << input_file;

  // Read met file name
  std::getline(input, line_);
  runtime_.met_file_name = TrimComment(line_);

  // Read infer test parameters
  ReadLine(input, runtime_.infer_test, runtime_.single_multi_flag, runtime_.mgap, runtime_.keep);

  if (runtime_.infer_test != 0 && runtime_.infer_test != 1) {
    LOG(FATAL) << "The infer test flag must be 0 or 1";
  }

  // Read surface parameters
  ReadLine(input, runtime_.gwl, runtime_.vegint, runtime_.nprint, runtime_.flprint, runtime_.sprint, runtime_.frozen, runtime_.mstflag);

  if (runtime_.single_multi_flag == 0) {
    ParseSinglePointMode(input);
  }
  else if (runtime_.single_multi_flag == 1) {
    ParseMultiPointMode(input);
  }
  else {
    LOG(FATAL) << "Invalid single_multi_flag: " << runtime_.single_multi_flag;
  }

  if (runtime_.infer_test != 0 && runtime_.infer_test != 1) {
    LOG(FATAL) << "The infer test flag must be 0 or 1";
  }

  LOG(INFO) << "Successfully read control file";
}

std::string Driver::TrimComment(const std::string& line)
{
  size_t pos = line.find('!');
  std::string result = (pos != std::string::npos) ? line.substr(0, pos) : line;
  return boost::trim_copy(result);
}

void Driver::ParseSinglePointMode(std::ifstream& input)
{
  // Read VITD index
  ReadLine(input, runtime_.vitd_index);

  // Read output file names
  std::getline(input, line_);
  runtime_.output1 = TrimComment(line_); // Met and surface data
  std::getline(input, line_);
  runtime_.output2 = TrimComment(line_); // Vertical profile data
  std::getline(input, line_);
  runtime_.output3 = TrimComment(line_); // Surface flux data
  std::getline(input, line_);
  runtime_.output4 = TrimComment(line_); // Vegetation temperature data
  std::getline(input, line_);
  runtime_.output5 = TrimComment(line_); // Snow processes data

  // Read slope and aspect
  ReadLine(input, met.slope, met.aspect);

  if (met.slope > 90.0) {
    LOG(FATAL) << "Slope cannot be greater than 90 degrees";
  }
  if (std::abs(met.aspect) > 360.0) {
    LOG(FATAL) << "Aspect cannot be greater than +/-360 degrees";
  }

  met.slope_radians = std::min(1.57, met.slope * M_PI / 180.0);
  met.slope_radians = std::round(met.slope_radians * 1e20) / 1e20;

  // Read surface roughness
  ReadLine(input, soil.roughness);

  // Read water flag
  ReadLine(input, met.water_flag);

  if (met.water_flag != 0 && met.water_flag != 1) {
    LOG(FATAL) << "Water flag must be 0 or 1";
  }

  // Handle water parameters
  if (met.water_flag == 1) {
    ReadLine(input, runtime_.water_type); // 1 = lakes/ponds, 2 = rivers
    ReadLine(input, runtime_.water_depth);
    ReadLine(input, runtime_.water_velocity);

    if (runtime_.water_type == 2) {
      met.water_flag = 2;
    }

    // Setup water layers
    soil.num_layers = 3;
    soil.soil_type(0) = 26; // Water
    soil.soil_type(1) = 26; // Water
    soil.soil_type(2) = 7;  // SM

    constexpr double kWaterLayerThickness = 0.02;
    soil.layer_thickness(0) = kWaterLayerThickness;

    if (!ApproxEqual(runtime_.water_depth, SystemLimits::kMissingValue)) {
      soil.layer_thickness(1) = runtime_.water_depth - kWaterLayerThickness;
    }
    else {
      soil.layer_thickness(1) = (runtime_.water_type == 1) ? 10.0 - kWaterLayerThickness : 3.0 - kWaterLayerThickness;
    }
    soil.layer_thickness(2) = 1.0;
  }

  // Process non-inferred parameters
  if (runtime_.infer_test == 0) {
    // Read snow parameters
    ReadLine(input, snow_ice.init_snow_depth, snow_ice.init_snow_water_equivalent);

    if (snow_ice.init_snow_depth < 0.0) {
      LOG(FATAL) << "Initial snow depth must be >= 0";
    }
    if (snow_ice.init_snow_water_equivalent < 0.0) {
      LOG(FATAL) << "Initial snow water equivalent must be >= 0";
    }

    // Read ice thickness
    ReadLine(input, snow_ice.init_ice_thickness);

    if (snow_ice.init_ice_thickness < 0.0) {
      LOG(FATAL) << "Ice thickness must be >= 0";
    }

    // Handle non-water surfaces
    if (met.water_flag == 0) {
      ParseVegetationParameters(input);
      ParseSoilParameters(input);
    }
  }
}

void Driver::ParseVegetationParameters(std::ifstream& input)
{
  // fasst_driver.F90:355

  // Read vegetation parameters

  // Low vegetation type
  ReadLine(input, veg.low_veg_type);

  if (veg.low_veg_type != 0) {
    veg.low_veg_flag = 1;
    ProcessVegetationType(veg.low_veg_type, veg.low_veg_type, veg.low_veg_flag);

    if (veg.low_veg_type >= 13 && veg.low_veg_type <= 15) {
      runtime_.water_type = 1;
      runtime_.water_velocity = 0.0;

      if (veg.low_veg_type == 13) {
        met.water_flag = 1;
        runtime_.water_depth = 1.5;
        runtime_.num_layers = 3;
        runtime_.soil_type(0) = 26;
        runtime_.layer_thickness(0) = 0.02;
        runtime_.soil_type(1) = 26;
        runtime_.layer_thickness(1) = runtime_.water_depth - 0.02;
        runtime_.soil_type(2) = 7;
        runtime_.layer_thickness(2) = 1.0;
      }
      else if (veg.low_veg_type == 14) {
        met.water_flag = 1;
        veg.low_veg_flag = 0;
        veg.low_veg_type = 0;
        runtime_.water_depth = 10.0;
        runtime_.num_layers = 3;
        runtime_.soil_type(0) = 26;
        runtime_.layer_thickness(0) = 0.02;
        runtime_.soil_type(1) = 26;
        runtime_.layer_thickness(1) = runtime_.water_depth - 0.02;
        runtime_.soil_type(2) = 7;
        runtime_.layer_thickness(2) = 1.0;
      }
      else {
        met.water_flag = 3;
        veg.low_veg_type = 0;
        veg.low_veg_flag = 0;
      }
    }

    // Read low vegetation parameters
    ReadLine(input, veg.low_foliage_density);

    if (veg.low_foliage_density > 0.98 && !ApproxEqual(veg.low_foliage_density, SystemLimits::kMissingValue)) {
      veg.low_foliage_density = 0.98;
    }

    // Read other low veg parameters
    ReadLine(input, veg.low_foliage_emissivity);
    ReadLine(input, veg.low_foliage_absorptivity);
    ReadLine(input, veg.low_foliage_height);

    if (ApproxZero(veg.low_foliage_density)) {
      veg.low_veg_type = 0;
      veg.low_veg_flag = 0;
    }

    if (veg.low_veg_type == 12) {
      veg.low_veg_flag = 0;
    }
  }

  // fasst_driver.F90:435
  // high vegetation type (trees)
  ReadLine(input, veg.veg_type_high);

  if (veg.veg_type_high != 0) {
    veg.veg_flag_high = 1;
    ProcessVegetationType(veg.veg_type_high, veg.veg_type_high, veg.veg_flag_high);

    // Check for water bodies
    if (veg.veg_type_high == 14 || veg.veg_type_high == 15) {
      runtime_.water_type = 1;
      veg.veg_type_high = 0;
      veg.veg_flag_high = 0;

      if (veg.veg_type_high == 14) { // Inland water
        met.water_flag = 1;
        runtime_.water_depth = 10.0;
        runtime_.num_layers = 3;
        runtime_.soil_type(0) = 26; // Water
        runtime_.layer_thickness(0) = 0.02;
        runtime_.soil_type(1) = 26; // Water
        runtime_.layer_thickness(1) = runtime_.water_depth - 0.02;
        runtime_.soil_type(2) = 7; // SM
        runtime_.layer_thickness(2) = 1.0;
      }
      else { // Ocean
        met.water_flag = 3;
      }
      runtime_.water_depth = 3.0;
      runtime_.water_velocity = 0.0;
    }

    // Read high vegetation parameters
    ReadLine(input, veg.high_foliage_height, veg.high_foliage_density);

    if (veg.high_foliage_density > 0.98 && !ApproxEqual(veg.high_foliage_density, SystemLimits::kMissingValue)) {
      veg.high_foliage_density = 0.98;
    }

    // Read canopy layer parameters
    for (int j = 0; j < SystemLimits::kCanopyLayers; j++) {
      ReadLine(input,
               veg.foliage_type(j),    // Foliage type for layer
               veg.leaf_area_index(j), // LAI
               veg.clumping(j),        // Markov clumping (typically 1.0)
               veg.dz(j),              // Layer thickness
               veg.reflectance(j),     // Surface reflectance
               veg.transmittance(j),   // Surface transmittance
               veg.alpha(j),           // Shortwave absorption (1-albedo)
               veg.emissivity(j)       // Longwave emissivity
      );

      // Set defaults if needed
      if (ApproxZero(veg.leaf_area_index(j))) {
        veg.leaf_area_index(j) = 0.5;
      }
      if (ApproxZero(veg.clumping(j))) {
        veg.clumping(j) = 0.01;
      }
    }

    if (ApproxZero(veg.high_foliage_density)) {
      veg.veg_type_high = 0;
      veg.veg_flag_high = 0;
    }

    if (veg.veg_type_high == 12) { // Ice cap/permanent snow
      veg.veg_flag_high = 0;
    }
  }
}

void Driver::ProcessVegetationType(int raw_type, int& type, int& flag)
{
  // Process vegetation type codes
  // Returns adjusted vegetation type after handling MODIS/UMD codes

  if (raw_type == 0) {
    return;
  }

  flag = 1;

  if (raw_type >= 1000 && raw_type < 2000) {
    veg.biome_source = 1000; // MODIS/IGBP
  }
  else if (raw_type >= 2000 && raw_type < 3000) {
    veg.biome_source = 2000; // UMD/LDAS
  }

  veg.new_veg_type_low = raw_type - veg.biome_source;
  if (veg.new_veg_type_low > 0) {
    if (veg.biome_source == 1000) { // MODIS
      if (veg.new_veg_type_low != 13 && veg.new_veg_type_low != 16) {
        type = modis_type_(veg.new_veg_type_low);
      }
      else {
        type = 0;
        flag = 0;
      }
    }
    else if (veg.biome_source == 2000) { // UMD
      if (veg.new_veg_type_low != 13) {
        type = umd_type_(veg.new_veg_type_low);
      }
      else {
        type = 0;
        flag = 0;
      }
    }
  }
}

void Driver::ParseSoilParameters(std::ifstream& input)
{
  // fasst_driver.F90:502

  // soil type, number of layers, layer thickness

  ReadLine(input, soil.num_layers);

  bool has_custom_soil = false;

  // Process each layer
  for (int i = 0; i < soil.num_layers; i++) {
    double layer_thickness;
    // Changed from soil.soil_type(i) to soil.soil_type(i)
    ReadLine(input, soil.soil_type(i), layer_thickness);
    soil.layer_thickness(i) = layer_thickness;
    if (std::abs(soil.soil_type(i)) < 100) {
      soil.soil_class_desc(i) = "USCS";
    }
    else if (std::abs(soil.soil_type(i)) >= 100) {
      soil.soil_class_desc(i) = "USDA";
      if (soil.soil_type(i) > 0) {
        int usda_type = soil.soil_type(i) - 100;
        soil.soil_type(i) = utils::MapUSDAToUSCS(usda_type);
      }
      else {
        soil.soil_type(i) = -(std::abs(soil.soil_type(i)) - 100);
      }
    }
    if (soil.soil_type(i) == -1) {
      has_custom_soil = true;
      runtime_.num_custom_soil_layers = i;
      ReadLine(input, soil.soil_type_desc(i));
    }
    if (soil.soil_type(i) == 26) {
      met.water_flag = 1;
      runtime_.water_flag = 1;
      runtime_.water_type = 1; // Default to lake/pond
      runtime_.water_velocity = 0.0;
    }
  }

  // Handle special water layer cases
  if (runtime_.water_flag == 1 && soil.num_layers <= 1) {
    soil.num_layers = 3;
    soil.soil_type(1) = soil.soil_type(0);
    soil.layer_thickness(1) = soil.layer_thickness(0) - 0.02;
    soil.soil_type(0) = 26; // Water
    soil.layer_thickness(0) = 0.02;
    soil.soil_type(2) = 7; // SM
    soil.layer_thickness(2) = 1.0;
  }
  else if (runtime_.num_layers > 0) {
    int total_layers = runtime_.num_layers + soil.num_layers; // kk
    int j = 0;
    utils::ArrayX<double, SystemLimits::kMaxNodes + 3> merged_layer_thickness{};
    utils::ArrayX<int, SystemLimits::kMaxNodes + 3> merged_soil_type{};
    double d1 = 0.0;
    double d2 = 0.0;

    for (int i = 0; i <= total_layers; ++i) {
      if (i <= runtime_.num_layers && runtime_.soil_type(i) == 26) {
        merged_layer_thickness(i) = d1 + runtime_.layer_thickness(i);
        merged_soil_type(i) = runtime_.soil_type(i);
      }
      else if (i > runtime_.num_layers &&
               soil.soil_type(i - runtime_.num_layers) == 26) {
        merged_layer_thickness(i) = d2 + soil.layer_thickness(i - runtime_.num_layers);
        merged_soil_type(i) = soil.soil_type(i - runtime_.num_layers);
      }
      if (merged_layer_thickness(i) > PhysicalConstants::kEpsilon) {
        j = i;
      }
    }
    // fasst_driver.F90:554
    int k;
    utils::ArrayX<double, SystemLimits::kMaxNodes + 3> layer_thickness3{};
    utils::ArrayX<int, SystemLimits::kMaxNodes + 3> soil_type3{};
    Sort2(j, merged_layer_thickness, merged_soil_type, k, layer_thickness3, soil_type3);

    d1 = 0.0;
    for (int i = 0; i < k; i++) {
      layer_thickness3(i) -= d1;
      d1 += layer_thickness3(i);
    }
    d2 = d1;

    // utils::ArrayX<double, SystemLimits::kMaxNodes + 3> layer_thickness2{};
    // utils::ArrayX<int, SystemLimits::kMaxNodes + 3> soil_type2{};

    for (int i = 0; i < total_layers; i++) {
      k = i - j;
      if (k > 0) {
        if (i <= runtime_.num_layers && runtime_.soil_type(i) != 26) {
          merged_layer_thickness(k) = d1 + runtime_.layer_thickness(i);
          merged_soil_type(k) = runtime_.soil_type(i);
        }
        else if (i > runtime_.num_layers && soil.soil_type(i - runtime_.num_layers) != 26) {
          merged_layer_thickness(k) = d2 + soil.layer_thickness(i - runtime_.num_layers);
          merged_soil_type(k) = soil.soil_type(i - runtime_.num_layers);
        }
      }
    }

    int jj;
    utils::ArrayX<double, SystemLimits::kMaxNodes + 3> layer_thickness4{};
    utils::ArrayX<int, SystemLimits::kMaxNodes + 3> soil_type4{};
    Sort2(k, merged_layer_thickness, merged_soil_type, jj, layer_thickness4, soil_type4);

    for (int i = 0; i < jj; i++) {
      layer_thickness4(i) -= d1;
      d1 += layer_thickness3(i);
    }

    soil.num_layers = k + jj;
    for (int i = 0; i < soil.num_layers; i++) {
      if (i < k) {
        soil.layer_thickness(i) = layer_thickness3(i);
        soil.soil_type(i) = soil_type3(i);
      }
      else {
        soil.layer_thickness(i) = layer_thickness4(i - k);
        soil.soil_type(i) = soil_type4(i - k);
      }
    }
  }

  // fasst_driver.F90:594
  if (runtime_.num_custom_soil_layers != 0) {
    ReadLine(input, runtime_.fasst_user_soil);
  }
  runtime_.fasst_user_soil = TrimComment(runtime_.fasst_user_soil);
  ReadLine(input, runtime_.init_temprature_flag);

  // Handle initial soil temperature data
  ReadLine(input, runtime_.init_temprature_flag);
  if (runtime_.init_temprature_flag == "y" || runtime_.init_temprature_flag == "Y") {
    ReadLine(input, runtime_.num_init_soil_temprature_layers);
    for (int j = 0; j < runtime_.num_init_soil_temprature_layers; j++) {
      ReadLine(input, soil.depth_temprature(j), soil.temperature(j));
      if (std::abs(soil.depth_temprature(j)) > 100.0) {
        LOG(FATAL) << "Measured soil temperatures out of range";
      }
      soil.temperature(j) += PhysicalConstants::kKelvinRef;
    }
  }

  // Handle initial soil moisture data
  ReadLine(input, runtime_.init_soil_moisture_flag);
  if (runtime_.init_soil_moisture_flag == "y" || runtime_.init_soil_moisture_flag == "Y") {
    ReadLine(input, runtime_.num_init_soil_moisture_layers);
    for (int j = 0; j < runtime_.num_init_soil_moisture_layers; j++) {
      ReadLine(input, soil.depth_moisture(j), soil.moisture(j));
      if (soil.depth_moisture(j) > 1.0) {
        LOG(FATAL) << "Measured soil moistures out of range";
      }
    }
  }
}

void Driver::ParseMultiPointMode(std::ifstream& input)
{
  // Read input files
  ReadLine(input, runtime_.vitd_input);
  runtime_.vitd_input = TrimComment(runtime_.vitd_input);

  ReadLine(input, runtime_.fasst_user_soil);
  runtime_.fasst_user_soil = TrimComment(runtime_.fasst_user_soil);
}

void Driver::ReadMeteorologicalData()
{
  // Read meteorological data similar to STEP 2
  // ... implementation
}

void Driver::ReadOldData()
{
  // Read old data if it exists similar to STEP 3
  // ... implementation
}

void Driver::InitializeSoilProfile()
{
  // Initialize soil profile similar to STEP 4
  // ... implementation
}

void Driver::OutputData()
{
  // Output data similar to STEP 5
  // ... implementation
}

void Driver::ProcessSinglePoint()
{
  // Process single point calculation
  // ... implementation
}

void Driver::ProcessMultiPoint()
{
  // Process multi point calculation
  // ... implementation
}

void Driver::CleanupFiles() {
  // Close all output files
  // ... implementation
}

void Sort2(int n,
           utils::ArrayX<double, SystemLimits::kMaxNodes + 3>& arr,
           utils::ArrayX<int, SystemLimits::kMaxNodes + 3>& arr1,
           int& ncount,
           utils::ArrayX<double, SystemLimits::kMaxNodes + 3>& b,
           utils::ArrayX<int, SystemLimits::kMaxNodes + 3>& b1)
{
  b.fill(0.0);
  b1.fill(0);

  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  // 根据 arr 排序下标，用以伴随排序
  std::sort(indices.begin(), indices.end(), [&](int A, int B) { return arr(A) < arr(B); });

  for (int i = 0; i < n; ++i) {
    arr(i) = arr(indices[i]);
    arr1(i) = arr1(indices[i]);
  }

  ncount = 1;
  b(0) = arr(0);
  b1(0) = arr1(0);

  for (int j = 1; j < n; ++j) {
    double a = arr(j) * 1000;
    if (std::fabs(a - arr(j - 1) * 1000) > 1e-4) {
      b(ncount) = arr(j);
      b1(ncount) = arr1(j);
      ++ncount;
    }
  }
}

} // namespace fasst
