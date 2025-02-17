# Create fasst library
add_library(fasstlib SHARED "")

# Add subdirectories
add_subdirectory(core)
add_subdirectory(modules)
add_subdirectory(utils)

# # Configure library settings
# target_include_directories(fasstlib
#   PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#     $<INSTALL_INTERFACE:include>
#   PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}
# )

target_link_libraries(fasstlib
  PUBLIC
    Boost::filesystem
  PRIVATE
    glog::glog
    Boost::program_options
)
