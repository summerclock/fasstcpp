target_sources(fasst
  PRIVATE
    Global.cpp
    Driver.cpp
)

target_link_libraries(fasst
  PRIVATE
  Boost::filesystem
  glog::glog
)