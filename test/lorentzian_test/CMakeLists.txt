cmake_minimum_required(VERSION 3.5)

project(lorentzian_test VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
set(CMAKE_CXX_FLAGS "-I/usr/include/  -O2")
set(CMAKE_C_FLAGS "-I/usr/include/  -O2")

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
    lorentzian_test.cpp
    ../../tamcmc/sources/acoefs.cpp
    ../../tamcmc/sources/function_rot.cpp
    ../../tamcmc/sources/build_lorentzian.cpp
    ../../tamcmc/sources/string_handler.cpp
    ../../external/Alm/Alm_cpp/bilinear_interpol.cpp
    ../../external/Alm/Alm_cpp/activity.cpp
    ../../external/Alm/Alm_cpp/Alm_interpol.cpp
    ../../external/Alm/Alm_cpp/linspace.cpp
    )

find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
	message(FATAL_ERROR "GSL not found. ajAlm function requires it.")
endif()
set(GSL_LIB GSL::gsl)

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}")
include_directories(../../tamcmc/headers)
include_directories(../../external/Alm/Alm_cpp)

find_package(Boost COMPONENTS system filesystem iostreams program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES_MAIN})
target_link_libraries(${CMAKE_PROJECT_NAME} ${Boost_LIBRARIES} ${GSL_LIB})
