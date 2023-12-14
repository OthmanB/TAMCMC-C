cmake_minimum_required(VERSION 3.5)
set(CMAKE_BUILD_TYPE Release)
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
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -lstdc++fs")
        if(WITH_OPENMP)
            string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fopenmp")
        endif()
        if(BUILD_ON_ZEN2)
            string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=znver2")
        endif()
    endif()

    if(APPLE)
        if(WITH_STATIC)
            message(FATAL_ERROR "Cannot compile statically on MacOS. Use it only for Linux systems. Turn it off by -DWITH_STATIC=OFF")
        endif()
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
        if(WITH_OPENMP)
            string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fopenmp")
        endif()
    endif()

else()
    message(STATUS "Building for NYU DALMA Supercomputer...")
    message(WARNING " - Remember to load the modules req uired for DALMA before running cmake")
    message(STATUS " - Adding optimisation flags for DALMA...")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fopenmp -march=westmere")
    option(USE_AVX "Build with AVX support" ON)
    if(USE_AVX)
        string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mavx2 -mfma -march=haswell")
    else()
        string(APPEND CMAKE_CXX_FLAGS_RELEASE " -msse4.2")
    endif()
endif()

# Check and Load dependencies
find_package(Gnuplot REQUIRED)
if(NOT GNUPLOT_FOUND)
	message(FATAL_ERROR "Gnuplot is required for rendering. Cannot pursue. Please install gnuplot")
endif()

# The check for OpenMP is only done once
if (WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()	
	message(FATAL_ERROR "OpenMP not found. Use -DWITH_OPENMP=OFF to build without OpenMP")
endif()

if(WITH_GSL)
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
#find_package(ZLIB REQUIRED)
# Set LibArchive_INCLUDE_DIR and LibArchive_LIBRARY in your bashrc
include_directories($ENV{LibArchive_INCLUDE_DIR})
find_library(LibArchive_LIBRARY NAMES archive)
find_path(LibArchive_INCLUDE_DIR archive.h)


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
	tamcmc/sources/acoefs.cpp
	external/ARMM/solver_mm.cpp
	external/ARMM/bump_DP.cpp
	external/Alm/Alm_cpp/activity.cpp
	external/Alm/Alm_cpp/Alm_interpol.cpp
	external/Alm/Alm_cpp/bilinear_interpol.cpp
        )
	
#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )
include_directories(tamcmc/headers)

#Dependencies that are not portable (BOOST)
find_package(Boost COMPONENTS system filesystem iostreams program_options REQUIRED)
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
	target_link_libraries(${CMAKE_PROJECT_NAME} ${Boost_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${LibArchive_LIBRARY}) # ${ZLIB_LIBRARIES}
	target_link_libraries(bin2txt ${Boost_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${LibArchive_LIBRARY})  # ${ZLIB_LIBRARIES}
	target_link_libraries(getstats ${Boost_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
	target_link_libraries(getmodel  ${Boost_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}) 
endif()