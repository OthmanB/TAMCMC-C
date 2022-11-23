cmake_minimum_required(VERSION 3.4.1)

project(acoefs_test VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
set(CMAKE_CXX_FLAGS "-I/usr/include/  -O2")
set(CMAKE_C_FLAGS "-I/usr/include/  -O2")

#Bring the headers into the project
set(HEADERS  #_COMMON
    ../../tamcmc/headers/acoefs.h	
   )

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
    acoefs_test.cpp
    ../../tamcmc/sources/acoefs.cpp
    )

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )

add_executable(${CMAKE_PROJECT_NAME} ${HEADERS} ${SOURCES_MAIN})
target_link_libraries(${CMAKE_PROJECT_NAME})
