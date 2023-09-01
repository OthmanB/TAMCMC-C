cmake_minimum_required(VERSION 3.4.1)

project(cpptamcmc VERSION 1.7.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

option(WITH_STATIC "Link Libraries Statically. Useful to port from one system to another" OFF)
option(WITH_GSL "GSL Library" ON)
option(WITH_OPENMP "OpenMP Library" ON)
option(BUILD_ON_DALMA "Building for NYU DALMA Supercomputer" OFF)
option(BUILD_ON_ZEN2 " Building on AMD Zen2 generation computers " OFF)
if(NOT BUILD_ON_DALMA)
	if(CMAKE_COMPILER_IS_GNUCXX)
		if(NOT BUILD_ON_ZEN2)		
			set(CMAKE_CXX_FLAGS "-O3 -fopenmp ")
		else()
			set(CMAKE_CXX_FLAGS "-O3 -fopenmp -march=znver2")
		endif()
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		if(NOT WITH_OPENMP)
			set(CMAKE_CXX_FLAGS "-O3")
		else()
			set(CMAKE_CXX_FLAGS "-O3 -fopenmp=libomp")
		endif()
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
		if(NOT WITH_OPENMP)
			set(CMAKE_CXX_FLAGS "-O3")
		else()
			set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
		endif()
	endif()
else()
	message(STATUS "Building for NYU DALMA Supercomputer...")
	message(WARNING "   - Remember to load the modules req	uired for DALMA before running cmake")
	message(STATUS  "   - Adding optimisation flags for DALMA...")
	# SSE / AVX switch
	option(USE_AVX "Build with AVX support" ON)
	if(USE_AVX)
	  set(CMAKE_CXX_FLAGS "-fopenmp -O3 -mavx2 -mfma -march=haswell")
	  #set(CMAKE_CXX_FLAGS "-O3 -mavx2 -mfma -march=haswell")
	else()
	  set(CMAKE_CXX_FLAGS "-O3 -fopenmp -msse4.2 -march=westmere")
	endif()
endif()

# Add the following line to include filesystem library
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lstdc++fs")

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
    tamcmc/sources/io_asymptotic.cpp
    tamcmc/sources/io_ajfit.cpp
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
    tamcmc/sources/linspace.cpp
    tamcmc/sources/acoefs.cpp
    tamcmc/sources/polyfit.cpp
    external/ARMM/solver_mm.cpp
    external/ARMM/bump_DP.cpp
    external/ARMM/linfit.cpp
	external/Alm/Alm_cpp/activity.cpp
	external/Alm/Alm_cpp/Alm_interpol.cpp
	external/Alm/Alm_cpp/bilinear_interpol.cpp
    )
 
#Define the sources for the bin2txt program
set(SOURCES_TOOLS_BIN2TXT
	tools/bin2txt_params.cpp
	tools/quick_samples_stats.cpp
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
	tamcmc/sources/io_asymptotic.cpp
	tamcmc/sources/io_ajfit.cpp
	tamcmc/sources/noise_models.cpp
	tamcmc/sources/priors_calc.cpp
	tamcmc/sources/stats_dictionary.cpp
	tamcmc/sources/derivatives_handler.cpp
	tamcmc/sources/io_models.cpp
	tamcmc/sources/linfit.cpp
    tamcmc/sources/polyfit.cpp
	tamcmc/sources/linspace.cpp
	tamcmc/sources/acoefs.cpp
	external/ARMM/solver_mm.cpp
	external/ARMM/bump_DP.cpp
	external/Alm/Alm_cpp/activity.cpp
	external/Alm/Alm_cpp/Alm_interpol.cpp
	external/Alm/Alm_cpp/bilinear_interpol.cpp
        )
	
# Check and Load dependencies
find_package(Gnuplot REQUIRED)
if(NOT GNUPLOT_FOUND)
	message(FATAL_ERROR "Gnuplot is required for rendering. Cannot pursue. Please install gnuplot")
endif()

if (WITH_OPENMP)
	find_package(OpenMP)
	#find_package(OPENMP)
	if (OPENMP_FOUND)
	    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	else()
	      find_package(OPENMP)
	      if (OPENMP_FOUND)
           	 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
           	 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
           	 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	      else()	
	      	message(FATAL_ERROR "OpenMP not found. Use -DWITH_OPENMP=OFF to build without OpenMP")
	      endif()
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
	if(WITH_STATIC)
    	set(GSL_LIB /usr/local/lib/libgsl.a)
    else()
    	set(GSL_LIB GSL::gsl)
    endif()
endif()

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )
include_directories(tamcmc/headers)

#Dependencies that are not portable (BOOST)
#set(Boost_USE_MULTITHREADED=OFF)
find_package(Boost COMPONENTS system filesystem iostreams program_options REQUIRED)
#find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES_MAIN} ${DIAG_SRC})
add_executable(bin2txt ${SOURCES_TOOLS_BIN2TXT} ${DIAG_SRC})
add_executable(getstats  ${SOURCES_TOOLS_GETSTATS})
add_executable(getmodel ${SOURCES_TOOLS_GETMODEL})

#target_link_libraries(bin2txt ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
if (WITH_STATIC)
	target_link_libraries(${CMAKE_PROJECT_NAME}  -static
		${Boost_LIBRARY_DIRS}/libboost_system.a
		${Boost_LIBRARY_DIRS}/libboost_filesystem.a
		${Boost_LIBRARY_DIRS}/libboost_iostreams.a
		${Boost_LIBRARY_DIRS}/libboost_program_options.a
		/usr/local/lib/libgsl.a
		-lz)
	target_link_libraries(bin2txt -static
  		${Boost_LIBRARY_DIRS}/libboost_system.a
  		${Boost_LIBRARY_DIRS}/libboost_filesystem.a
  		${Boost_LIBRARY_DIRS}/libboost_iostreams.a
  		${Boost_LIBRARY_DIRS}/libboost_program_options.a
  		/usr/local/lib/libgsl.a)
	target_link_libraries(getstats -static
		${Boost_LIBRARY_DIRS}/libboost_system.a
		${Boost_LIBRARY_DIRS}/libboost_filesystem.a
		${Boost_LIBRARY_DIRS}/libboost_iostreams.a
		${Boost_LIBRARY_DIRS}/libboost_program_options.a
		/usr/local/lib/libgsl.a)
	target_link_libraries(getmodel -static
		${Boost_LIBRARY_DIRS}/libboost_system.a
		${Boost_LIBRARY_DIRS}/libboost_filesystem.a
		${Boost_LIBRARY_DIRS}/libboost_iostreams.a
		${Boost_LIBRARY_DIRS}/libboost_program_options.a
		/usr/local/lib/libgsl.a
		-lz)
else()
	target_link_libraries(${CMAKE_PROJECT_NAME} ${Boost_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}) # Boost::boost)
	target_link_libraries(bin2txt ${Boost_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
	target_link_libraries(getstats ${Boost_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})# Boost::boost)
	target_link_libraries(getmodel  ${Boost_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}) # Boost::boost)
endif()



