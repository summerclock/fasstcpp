# Create fasst library
add_library(fasst SHARED "")

# Add subdirectories
add_subdirectory(core)
add_subdirectory(modules)
add_subdirectory(utils)

# # Configure library settings
# target_include_directories(fasst
#   PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#     $<INSTALL_INTERFACE:include>
#   PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}
# )

target_link_libraries(fasst
  PUBLIC
    Boost::filesystem
  PRIVATE
    glog::glog
    Boost::program_options
)

# Set properties
set_target_properties(fasst PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)
