cmake_minimum_required(VERSION 3.10)
project(fasst_utils)

target_sources(fasstlib
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Functions.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Functions.h
)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

file(GLOB TEST_SOURCES "*.test.cpp")

add_executable(matrix_tests ${TEST_SOURCES})
target_link_libraries(matrix_tests GTest::GTest GTest::Main)

# add_test(NAME matrix_tests COMMAND matrix_tests)
gtest_add_tests(TARGET matrix_tests)
