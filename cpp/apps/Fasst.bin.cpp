#include "fasst/core/Driver.h"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  try {
    std::string input_file;
    if (argc > 1) {
      input_file = argv[1];
    } else {
      input_file = "gr1_zip.inp";  // Default input file
    }

    fasst::Driver driver;
    driver.Run(input_file);
    return 0;
    
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
