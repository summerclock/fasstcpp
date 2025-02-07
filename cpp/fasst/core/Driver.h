#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

#include <glog/logging.h>

#include "fasst/core/Global.h"
#include "fasst/utils/Array.h"

namespace fasst {

class Driver
{
public:
  Driver();
  ~Driver() = default;

  // Disallow copy
  Driver(const Driver&) = delete;
  Driver& operator=(const Driver&) = delete;

  // Main entry point
  void Run(const std::string& input_file);

private:
  // Step functions corresponding to original Fortran program sections
  void ReadControlFile(const std::string& input_file);
  void ReadMeteorologicalData();
  void ReadOldData();
  void InitializeSoilProfile();
  void OutputData();

  // Helper methods
  template<typename T>
  void StreamToValue(std::stringstream& ss, T& value)
  {
    ss >> value;
    DCHECK(!ss.fail()) << "Failed to read value from stream";
  }
  template<typename T, typename... Args>
  void StreamToValue(std::stringstream& ss, T& value, Args&... args)
  {
    ss >> value;
    DCHECK(!ss.fail()) << "Failed to read value from stream";
    StreamToValue(ss, args...);
  }
  template<typename T, typename... Args>
  void ReadLine(std::ifstream& input, T& value, Args&... args)
  {
    std::string line;
    std::getline(input, line);
    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line);
    StreamToValue(ss, value, args...);
  }
  void ParseSinglePointMode(std::ifstream& input);
  void ParseMultiPointMode(std::ifstream& input);
  std::string TrimComment(const std::string& line);
  void ParseVegetationParameters(std::ifstream& input);
  void ProcessVegetationType(int raw_type, int& type, int& flag);
  void ParseSoilParameters(std::ifstream& input);
  void ProcessSinglePoint();
  void ProcessMultiPoint();
  void CleanupFiles();

  // Vegetation type conversions
  // convert modis and umd vegetation types to fasst defaults
  // biome_sourece: FASST default = 0, MODIS/IGBP = 1000, UMD/LDAS = 2000
  //                      1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  utils::ArrayX<int, 20> modis_type_{ { 3, 6, 4, 5, 18, 16, 16, 18, 7, 2, 13, 10, 0, 1, 12, 8, 14, 9, 9, 9 } };
  utils::ArrayX<int, 14> umd_type_{ { 3, 6, 5, 4, 18, 3, 2, 11, 11, 2, 1, 8, 0, 15 } };

  // Runtime parameters
  struct RuntimeParams
  {
    int infer_test{ 0 };
    int single_multi_flag{ 0 };
    double mgap{ 0 };
    int keep{ 0 };
    double gwl{ 0.0 };
    double vegint{ 0.0 };
    int nprint{ 0 };
    int flprint{ 0 };
    int sprint{ 0 };
    int frozen{ 100 };
    int mstflag{ 0 };

    int error_code{ 0 };                                                 // error_code: main error flag
    int error_type{ 0 };                                                 // error_type: type of error
    int step{ 0 };                                                       // step: current step counter
    int step_i{ 0 };                                                     // stepi: initial step
    int met_pos{ 0 };                                                    // mpos: met file position
    int dstart{ 0 };                                                     // dstart: data start position
    int wstart{ 0 };                                                     // wstart: water calculations start
    int freq_id{ 0 };                                                    // freq_id: frequency identifier
    int vitd_index{ 0 };                                                 // vitd_index: vegetation index
    int water_type{ 0 };                                                 // wtype: 1=lakes/ponds, 2=rivers
    int error_count{ 0 };                                                // ecount: error counter
    int first_loop{ 0 };                                                 // first_loop: initialization flag
    int num_custom_soil_layers{ 0 };                                     // lcount: custom soil layer count
    int water_flag{ 0 };                                                 // wflag2
    int num_layers{ 0 };                                                 // nlayers1
    int num_init_soil_temprature_layers{ 0 };                            // nt0: number of initial soil temperature layers
    int num_init_soil_moisture_layers{ 0 };                              // nm0: number of initial soil moisture layers
    double water_velocity{ 0 };                                          // wvel: water velocity
    double water_depth{ 0 };                                             // wdepth: water depth
    utils::Array<int, 3, 1> soil_type = utils::Array<int, 3, 1>::Zero(); // soiltype1
    utils::ArrayX<double, 3> layer_thickness{};                          // lthick1: layer thickness
    std::string melt_flag;                                               // meltfl: melting flag
    std::string met_file_name;
    std::string vitd_input;
    std::string fasst_user_soil{ "fasstusersoil.inp" }; // fasstusersoil: user defined soil parameters file
    std::string init_temprature_flag;                   // ttest: initial temperature flag, "y" = have initial soil temps in C
    std::string init_soil_moisture_flag;                // mtest: initial soil moisture flag, "y" = have initial soil moisture (fraction)
    std::string output1;                                // output1: met and surface data
    std::string output2;                                // output2: vertical profile data
    std::string output3;                                // output3: surface flux data
    std::string output4;                                // output4: vegetation temperature data
    std::string output5;                                // output5: snow processes data
  } runtime_;

  std::string line_;

  // Output file streams
  std::ofstream output1_; // met and surface data
  std::ofstream output2_; // vertical profile data
  std::ofstream output3_; // surface flux data
  std::ofstream output4_; // vegetation temperature data
  std::ofstream output5_; // snow processes data
};

void Sort2(int n,
           utils::ArrayX<double, SystemLimits::kMaxNodes + 3>& arr,
           utils::ArrayX<int, SystemLimits::kMaxNodes + 3>& arr1,
           int& ncount,
           utils::ArrayX<double, SystemLimits::kMaxNodes + 3>& b,
           utils::ArrayX<int, SystemLimits::kMaxNodes + 3>& b1);

} // namespace fasst
