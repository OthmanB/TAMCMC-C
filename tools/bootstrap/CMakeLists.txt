cmake_minimum_required(VERSION 3.12)
project(BootstrapTests)

# Find Boost
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )

# Add the executable for the main program
add_executable(bootstrap_tests bootstrap_algos.cpp bootstrap_algos_unit_test.cpp)

# Link the Boost libraries
target_include_directories(bootstrap_tests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(bootstrap_tests PRIVATE ${Boost_LIBRARIES})

# Add a target for running the tests
add_test(NAME BootstrapTests COMMAND bootstrap_tests)

# Enable testing in CMake
enable_testing()
