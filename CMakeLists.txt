cmake_minimum_required(VERSION 3.4.1)

project(cpptamcmc VERSION 1.3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)

option(WITH_GSL "GSL Library" ON)
option(WITH_OPENMP "OpenMP Library" ON)
option(BUILD_ON_DALMA "Building for NYU DALMA Supercomputer" OFF)

if(NOT BUILD_ON_DALMA)
	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "-O3")
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		set(CMAKE_CXX_FLAGS "-O3")
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
		set(CMAKE_CXX_FLAGS "-O3")
	endif()
else()
	message(STATUS "Building for NYU DALMA Supercomputer...")
	message(WARNING "   - Remember to load the modules required for DALMA before running cmake")
	message(STATUS  "   - Adding optimisation flags for DALMA...")
	# SSE / AVX switch
	option(USE_AVX "Build with AVX support" ON)
	if(USE_AVX)
	  set(CMAKE_CXX_FLAGS "-fopenmp -O3 -mavx2 -mfma -march=haswell")
	  #set(CMAKE_CXX_FLAGS "-O3 -mavx2 -mfma -march=haswell")
	else()
	  set(CMAKE_CXX_FLAGS "-O3 -msse4.2 -march=westmere")
	endif()
endif()

#Bring the headers into the project
set(HEADERS  #_COMMON
    tamcmc/headers/MALA.h
    tamcmc/headers/data.h
    tamcmc/headers/gnuplot-iostream.h
    tamcmc/headers/matrices.h
    tamcmc/headers/outputs.h
    tamcmc/headers/build_lorentzian.h
    tamcmc/headers/derivatives_handler.h
    tamcmc/headers/interpol.h
    tamcmc/headers/model_def.h
    tamcmc/headers/priors_calc.h
    tamcmc/headers/config.h
    tamcmc/headers/diagnostics.h
    tamcmc/headers/io_ms_global.h
    tamcmc/headers/io_local.h
    tamcmc/headers/models.h
    tamcmc/headers/stats_dictionary.h
    tamcmc/headers/config_presets.h
    tamcmc/headers/function_rot.h
	tamcmc/headers/likelihoods.h
    tamcmc/headers/noise_models.h
    tamcmc/headers/string_handler.h
    tamcmc/headers/io_models.h
   )

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
    tamcmc/sources/main.cpp
    tamcmc/sources/MALA.cpp
    tamcmc/sources/matrices.cpp
    tamcmc/sources/outputs.cpp
    tamcmc/sources/build_lorentzian.cpp
    tamcmc/sources/derivatives_handler.cpp
    tamcmc/sources/interpol.cpp
    tamcmc/sources/model_def.cpp
    tamcmc/sources/priors_calc.cpp
    tamcmc/sources/config.cpp
    tamcmc/sources/diagnostics.cpp
    tamcmc/sources/io_ms_global.cpp
    tamcmc/sources/io_local.cpp
    tamcmc/sources/models.cpp
    tamcmc/sources/stats_dictionary.cpp
    tamcmc/sources/config_presets.cpp
    tamcmc/sources/function_rot.cpp
    tamcmc/sources/likelihoods.cpp
    tamcmc/sources/noise_models.cpp
    tamcmc/sources/string_handler.cpp
    tamcmc/sources/ben_timer.cpp
    tamcmc/sources/multivariate_cholesky.cpp
    tamcmc/sources/random_JB.cpp
    tamcmc/sources/io_models.cpp
    )
 
#Define the sources for the bin2txt program
set(SOURCES_TOOLS_BIN2TXT
	tools/bin2txt_params.cpp
	tamcmc/sources/interpol.cpp
	tamcmc/sources/diagnostics.cpp
	tamcmc/sources/string_handler.cpp
	)
#Define the sources for the getstats program
set(SOURCES_TOOLS_GETSTATS
	tools/read_stats.cpp
	tamcmc/sources/string_handler.cpp
	)

#Define the sources for the getmodel program
set(SOURCES_TOOLS_GETMODEL
	tools/getmodel.cpp
	tamcmc/sources/models.cpp
	tamcmc/sources/model_def.cpp
	tamcmc/sources/string_handler.cpp
	tamcmc/sources/config.cpp
	tamcmc/sources/function_rot.cpp
	tamcmc/sources/likelihoods.cpp
	tamcmc/sources/build_lorentzian.cpp
	tamcmc/sources/interpol.cpp
	tamcmc/sources/matrices.cpp
	tamcmc/sources/io_ms_global.cpp
	tamcmc/sources/io_local.cpp
	tamcmc/sources/noise_models.cpp
	tamcmc/sources/priors_calc.cpp
	tamcmc/sources/stats_dictionary.cpp
	tamcmc/sources/derivatives_handler.cpp
	tamcmc/sources/io_models.cpp
	)
	
# Check and Load dependencies
find_package(Gnuplot REQUIRED)
if(NOT GNUPLOT_FOUND)
	message(FATAL_ERROR "Gnuplot is required for rendering. Cannot pursue. Please install gnuplot")
endif()

if (WITH_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
	    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	else()
		message(FATAL_ERROR "OpenMP not found. Use -DWITH_OPENMP=OFF to build without OpenMP")
	endif()
endif()

if(NOT WITH_GSL)
	#set(DIAG_SRC tamcmc/sources/nogsl/diagnostics.cpp tamcmc/sources/nogsl/diagnostics.h)
else() 
	add_definitions(/DTAMCMC_WITH_GSL)
	find_package(GSL REQUIRED)
	if(NOT GSL_FOUND)
		message(FATAL_ERROR "GSL not found. Use -DWITH_GSL=OFF to build without GSL")
	endif()
endif()

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )
include_directories(tamcmc/headers)

#Dependencies that are not portable (BOOST)
set(Boost_USE_MULTITHREADED=OFF)
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(... ${Boost_INCLUDE_DIRS})
link_directories(... ${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${HEADERS} ${SOURCES_MAIN} ${DIAG_SRC})
target_link_libraries(${CMAKE_PROJECT_NAME}  ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE} ${GSL_CBLAS_LIBRARY} ${GSL_LIBRARY})
#target_link_libraries(${CMAKE_PROJECT_NAME} sharedlib)

add_executable(bin2txt ${HEADERS} ${SOURCES_TOOLS_BIN2TXT} ${DIAG_SRC})
target_link_libraries(bin2txt ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE} ${GSL_CBLAS_LIBRARY} ${GSL_LIBRARY})

add_executable(getstats ${HEADERS} ${SOURCES_TOOLS_GETSTATS})
target_link_libraries(getstats ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE} ${GSL_CBLAS_LIBRARY} ${GSL_LIBRARY})

add_executable(getmodel ${HEADERS} ${SOURCES_TOOLS_GETMODEL})
target_link_libraries(getmodel ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE} ${GSL_CBLAS_LIBRARY} ${GSL_LIBRARY})




