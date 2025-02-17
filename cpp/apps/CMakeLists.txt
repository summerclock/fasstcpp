project(fasst)

add_executable(fasst Fasst.bin.cpp)

target_link_libraries(fasst PRIVATE fasstlib)

install(TARGETS fasst DESTINATION bin)
