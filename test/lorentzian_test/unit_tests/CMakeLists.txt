cmake_minimum_required(VERSION 3.5)
project(ModesTests)

set(CMAKE_CXX_STANDARD 17)

# Set compiler optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)

# Find Boost packages
#find_package(Boost COMPONENTS system filesystem iostreams program_options REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem iostreams unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
	message(FATAL_ERROR "GSL not found. GSL is required for the interpolation")
endif()
set(GSL_LIB GSL::gsl)

# Find Eigen and OpenMP package depending on OS
if(APPLE)
    find_package(OpenMP REQUIRED COMPONENTS C)
    find_package(Eigen3 REQUIRED)
else()
    find_package(OpenMP REQUIRED)
endif()
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )
include_directories(../../../tamcmc/headers)
include_directories(../../../external/Alm/Alm_cpp)

# Bring the core program sources into the project
set(SOURCES_TESTS_BUILD_L
    ../../../tamcmc/sources/build_lorentzian.cpp
    ../../../tamcmc/sources/models.cpp
    ../../../tamcmc/sources/noise_models.cpp
    ../../../tamcmc/sources/function_rot.cpp
    ../../../tamcmc/sources/string_handler.cpp
    ../../../tamcmc/sources/acoefs.cpp
    ../../../tamcmc/sources/interpol.cpp
    ../../../tamcmc/sources/linfit.cpp
    ../../../tamcmc/sources/derivatives_handler.cpp
    ../../../external/Alm/Alm_cpp/activity.cpp
    ../../../external/Alm/Alm_cpp/Alm_interpol.cpp
    ../../../external/Alm/Alm_cpp/bilinear_interpol.cpp
    ../../../external/ARMM/configure_make_star.cpp
    ../../../external/ARMM/readparams_job.cpp
    ../../../external/ARMM/bump_DP.cpp
    ../../../external/ARMM/solver_mm.cpp
    ref_code/build_lorentzian_ref.cpp
    ref_code/models_ref.cpp
    ref_code/models_ref_parallel_RGBonly.cpp
    test_build_l_mode.cpp
)


# Add executable
add_executable(BuildModeTests ${SOURCES_TESTS_BUILD_L})

# Link Eigen and Boost libraries
target_link_libraries(BuildModeTests Eigen3::Eigen ${Boost_LIBRARIES} ${GSL_LIB})

