project(FasstApp)

add_executable(Fasst Fasst.bin.cpp)

target_link_libraries(Fasst PRIVATE fasst)

install(TARGETS Fasst DESTINATION bin)
