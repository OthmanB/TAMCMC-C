#define BOOST_TEST_MODULE Build_l_mode_TESTS
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Dense>
#include "colors.hpp"
#include <chrono>
#include <random>
#include <iostream>
#include "../../../external/ARMM/configure_make_star.h"
#include "../../../external/ARMM/readparams_job.h"
#include "../../../tamcmc/headers/function_rot.h"
#include "../../../tamcmc/headers/string_handler.h"
#include "../../../tamcmc/headers/build_lorentzian.h"
#include "../../..//external/Alm/Alm_cpp/Alm_interpol.h"
#include "../../../tamcmc/headers/models.h"
#include "ref_code/build_lorentzian_ref.h"
#include "ref_code/models_ref.h"
#include "ref_code/models_ref_parallel_RGBonly.h"

// For ARMM models
#include "../../../external/ARMM/data_solver.h"

using Eigen::VectorXd;
using namespace std::chrono;;

struct Setup{
    VectorXd params;
    VectorXi params_length;
};


Setup make_params_RGB_model(Setup inputs_ref, const int lmax, const int Nfreqs, const int Nspline_ref, const int n0, 
        const Cfg_synthetic_star cfg_star, const bool is_old, const std::string width_type , const bool fix_compare=false , const bool exclude_nobias_test=false);
Setup make_params_aj_model(const int lmax, const int Nfreqs, const std::string model_name, 
        const long double Dnu, const long double epsilon, const long double d0l);
external_data load_Alm_tables(const std::string);

/*
BOOST_AUTO_TEST_CASE(build_l_mode_a1l_etaa3_test)
{  
    std::cout << colors::cyan << "Testing output difference between reference code and current TAMCMC build_l_mode_a1l_etaa3..." << colors::white << std::endl;

    const double tolerance = 1e-6; // in percent
    const int numTests = 1000;
    const double resol=1e6/(4.*365.*86400.);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform01(0, 1);
    std::uniform_int_distribution<int> Uniform03(0, 3);
    
    long duration_total1=0, duration_total2=0;
    for (int i = 0; i < numTests; i++) { 
        // p modes
        int l = Uniform03(gen);
        long double fc_l=std::uniform_real_distribution<>(950., 1050.)(gen); 
        long double H_l=1;
        long double gamma_l=1;
        long double f_s1=std::uniform_real_distribution<>(0.5, 5.)(gen); 
        long double f_s2=std::uniform_real_distribution<>(0.5, 5.)(gen); 
        long double eta0=0;
        long double a3=f_s1/100;
        int sw=Uniform01(gen);
        long double asym;
        if (sw ==1){ // This is a switch to test asym = 0 or asym !=0 cases in the build_l_mode
            asym=std::uniform_real_distribution<>(-100, 100.)(gen);
        } else{
            asym=0;
        }
        long double inc=std::uniform_real_distribution<>(0., 90.)(gen); 
        VectorXd V=amplitude_ratio(l, inc);
        long double fmin=fc_l - 500;
        long double fmax=fc_l + 500;
        
        int Ndata=std::ceil((fmax - fmin)/resol);
        Eigen::VectorXd x_l = Eigen::VectorXd::LinSpaced(Ndata, fmin, fmax);
        // --------------------------------------------
        
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new=build_l_mode_a1l_etaa3(x_l, H_l, fc_l, f_s1, f_s2, eta0, a3, asym, gamma_l, l, V);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in seconds
        auto start2 = high_resolution_clock::now();
        VectorXd sols_ref=build_l_mode_a1l_etaa3_ref(x_l, H_l, fc_l, f_s1, f_s2, eta0, a3, asym, gamma_l, l, V);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in seconds
        // Check if the solutions are the same within tolerance
        if ((sols_new - sols_ref).norm() > tolerance/100){
            std::cout << "(sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
        }
    // Compare execution times
    std::cout << "      build_l_mode_a1l_etaa3 NEW: total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "      build_l_mode_a1l_etaa3 REF: total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
}
*/

BOOST_AUTO_TEST_CASE(model_MS_Global_aj_HarveyLike_TEST)
{  
    std::cout << colors::cyan << "Testing output difference between reference code and current TAMCMC model_MS_Global_aj_HarveyLike..." << colors::white << std::endl;
 
    const double tolerance = 1e-6; // in percent
    const int numTests_Acc = 10;
    const int numTests_Time = 50;
    const double resol=1e6/(4.*365.*86400.);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform23(2, 3);
    const bool output_param=false;
    const int Nfreqs=5;
    long duration_total1=0, duration_total2=0;
    
    std::cout << colors::blue << "    ACCURACY TEST" << colors::white<< std::endl;
    long double fmin=0;
    long double fmax=Nfreqs*130 + 250; // 130 is the approximate Dnu randomized within make_params_aj_model()
    for (int i = 0; i < numTests_Acc; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        int lmax = Uniform23(gen);
        long double Dnu=std::uniform_real_distribution<>(129, 130.)(gen);
        long double epsilon=std::uniform_real_distribution<>(0., 0.05)(gen);
        long double d0l=std::uniform_real_distribution<>(-Dnu*0.02, 0)(gen);
        Setup set1=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike", Dnu, epsilon, d0l);
        //Setup set2=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike");
        int Ndata=std::ceil((fmax - fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, fmin, fmax);
        // --------------------------------------------
        VectorXd sols_new=model_MS_Global_aj_HarveyLike(set1.params, set1.params_length, x, output_param);
        VectorXd sols_ref=model_MS_Global_aj_HarveyLike_ref(set1.params, set1.params_length, x, output_param);
        // Check if the solutions are the same within tolerance
        std::cout << " (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
        //std::cout << "(sols_new - sols_ref).sum()  = " << (sols_new - sols_ref).sum()  << std::endl;
        if ((sols_new - sols_ref).norm() > tolerance/100){
            std::cout << "(sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
        }
    std::cout << colors::blue << "    SPEED TEST " << colors::white <<  std::endl;
    fmin=0;
    fmax=Nfreqs*130 + 2500; // 130 is the approximate Dnu randomized within make_params_aj_model()
    for (int i = 0; i < numTests_Time; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        long double Dnu=std::uniform_real_distribution<>(129, 130.)(gen);
        long double epsilon=std::uniform_real_distribution<>(0., 0.05)(gen);
        long double d0l=std::uniform_real_distribution<>(-Dnu*0.02, 0)(gen);
        int lmax = Uniform23(gen);
        Setup set1=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike", Dnu, epsilon, d0l);
        //Setup set2=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike");

        int Ndata=std::ceil((fmax - fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, fmin, fmax);
        // --------------------------------------------
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new=model_MS_Global_aj_HarveyLike(set1.params, set1.params_length, x, output_param);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in seconds
        auto start2 = high_resolution_clock::now();
        VectorXd sols_ref=model_MS_Global_aj_HarveyLike_ref(set1.params, set1.params_length, x, output_param);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in seconds
    }
    // Compare execution times
    std::cout << std::endl;
    std::cout << "      model_MS_Global_aj_HarveyLike NEW: total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "      model_MS_Global_aj_HarveyLike REF: total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
    std::cout << " --- " << std::endl;
}

BOOST_AUTO_TEST_CASE(model_MS_Global_ajAlm_HarveyLike_TEST)
{  
    std::cout << colors::cyan << "Testing output difference between reference code and current TAMCMC model_MS_Global_ajAlm_HarveyLike..." << colors::white << std::endl;

    const double tolerance = 1e-6; // in percent
    const int numTests_Acc = 10;
    const int numTests_Time = 50;
    const double resol=1e6/(4.*365.*86400.);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform23(2, 3);
    std::uniform_int_distribution<int> Uniform02(0, 2);
    std::uniform_int_distribution<int> Uniform_m12(-1, 2);
    
    const bool output_param=false;
    const int Nfreqs=8;
    long duration_total1=0, duration_total2=0;
    
    std::string Alm_tables_dir="../../../../external/Alm/data/Alm_grids_CPP/1deg_grids/";
    external_data Alm_data=load_Alm_tables(Alm_tables_dir);

    std::cout << colors::blue << "    ACCURACY TEST" << colors::white<< std::endl;
    long double fmin=0;
    long double fmax=Nfreqs*130 + 250; // 130 is the approximate Dnu randomized within make_params_aj_model()
    for (int i = 0; i < numTests_Acc; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        int lmax = Uniform23(gen);
        int decompose_Alm=Uniform_m12(gen);
        int filter_code=Uniform02(gen);
        while (filter_code == 1){ // filter_code must be 0, or 2 because case 1 is Gaussian and unprecise with the grid method
            filter_code=Uniform02(gen); // Either Gate or Triangle case
        }
        long double Dnu=std::uniform_real_distribution<>(129, 130.)(gen);
        long double epsilon=std::uniform_real_distribution<>(0., 0.05)(gen);
        long double d0l=std::uniform_real_distribution<>(-Dnu*0.02, 0)(gen);
        Setup set1=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_ajAlm_HarveyLike", Dnu, epsilon, d0l);
        set1.params.conservativeResize(set1.params.size() + 1 ); // We need to add decompose_Alm and filter_code
        set1.params[set1.params.size()-2]=decompose_Alm;
        set1.params[set1.params.size()-1]=filter_code;
        set1.params_length[set1.params_length.size()-1] = set1.params_length[set1.params_length.size()-1] + 2; // We update the parameters length accordingly
        //Setup set2=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_ajAlm_HarveyLike");
        int Ndata=std::ceil((fmax - fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, fmin, fmax);
        // --------------------------------------------
        VectorXd sols_new=model_MS_Global_ajAlm_HarveyLike(set1.params, set1.params_length, x, output_param, Alm_data);
        VectorXd sols_ref=model_MS_Global_ajAlm_HarveyLike_ref(set1.params, set1.params_length, x, output_param, Alm_data);
        // Check if the solutions are the same within tolerance
        std::cout << " (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
        //std::cout << "(sols_new - sols_ref).sum()  = " << (sols_new - sols_ref).sum()  << std::endl;
        if ((sols_new - sols_ref).norm() > tolerance/100){
            //std::cout << "     (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
        }
    std::cout << colors::blue << "    SPEED TEST " << colors::white <<  std::endl;
    fmin=0;
    fmax=Nfreqs*130 + 2500; // 130 is the approximate Dnu randomized within make_params_aj_model()
    for (int i = 0; i < numTests_Time; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        long double Dnu=std::uniform_real_distribution<>(129, 130.)(gen);
        long double epsilon=std::uniform_real_distribution<>(0., 0.05)(gen);
        long double d0l=std::uniform_real_distribution<>(-Dnu*0.02, 0)(gen);
        int lmax = Uniform23(gen);
        int decompose_Alm=Uniform_m12(gen);
        int filter_code=Uniform02(gen);
        while (filter_code == 1){ // filter_code must be 0, or 2 because case 1 is Gaussian and unprecise with the grid method
            filter_code=Uniform02(gen); // Either Gate or Triangle case
        }
        Setup set1=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_ajAlm_HarveyLike", Dnu, epsilon, d0l);
        set1.params.conservativeResize(set1.params.size() + 1 ); // We need to add decompose_Alm and filter_code
        set1.params[set1.params.size()-2]=decompose_Alm;
        set1.params[set1.params.size()-1]=filter_code;
        set1.params_length[set1.params_length.size()-1] = set1.params_length[set1.params_length.size()-1] + 2; // We update the parameters length accordingly
        //Setup set2=make_params_aj_model(lmax, Nfreqs,"model_MS_Global_ajAlm_HarveyLike");

        int Ndata=std::ceil((fmax - fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, fmin, fmax);
        // --------------------------------------------
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new=model_MS_Global_ajAlm_HarveyLike(set1.params, set1.params_length, x, output_param, Alm_data);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in seconds
        auto start2 = high_resolution_clock::now();
        VectorXd sols_ref=model_MS_Global_ajAlm_HarveyLike_ref(set1.params, set1.params_length, x, output_param, Alm_data);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in seconds
    }
    // Compare execution times
    std::cout << std::endl;
    std::cout << "      model_MS_Global_aj_HarveyLike NEW: total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "      model_MS_Global_aj_HarveyLike REF: total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
    std::cout << " --- " << std::endl;
}


BOOST_AUTO_TEST_CASE(model_RGB_asympt_aj_AppWidth_HarveyLike_v4_PARELLELISATION_TEST)
{  
    // This test takes the reference function for the RGB star and introduce its parallelisation 
    // at the model level. Please note that this does not test parallelisation of the ARMM part.
    // This ARMM parallelisation is already included here and is separately tested in unit tests 
    // at github.com/OthmanB/ARMM-solver in the unit_test v1.0--v1.1 branch.
    // Thus, performance gains that are noted here are just for the model paralleliation.
    // An estimate of x1.5 - 2 can be achieved here while the ARMM parallel gives x2 - 5.
    // This leads to an overall total gain of x3 - x10

    const double tolerance = 1e-6; // in percent
    const int numTests_Acc = 10;
    const int numPoints = 20;
    const double resol=1e6/(4.*365.*86400.);
    const std::string cfg_file="../config_ARMM/make_star_test.cfg";
    
    const bool output_param=false;
    const int Nfreqs=4;
    const int Nspline_ref=4;
    const int n0=5;
    
    Eigen::VectorXd nu(numPoints);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform23(2, 3);

    // Load a reference configuration for l=1 mixed modes
    std::unordered_map<std::string, std::string> input_params=readParameterFile(cfg_file);
    Cfg_synthetic_star cfg_star=configure_make_star(input_params); 
    
    long duration_total1=0, duration_total2=0;
    std::cout << colors::cyan << "Testing output difference between make_star with and without all improvments (Fully randomized)..." << colors::white << std::endl;
    for (int i = 0; i < numTests_Acc; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        int lmax = Uniform23(gen);
        // p modes
        cfg_star.Dnu_star=std::uniform_real_distribution<>(10., 40.)(gen); //60
        // g modes
        cfg_star.DPl_star=std::uniform_real_distribution<>(100., 400.)(gen); // 400
        cfg_star.q_star=std::uniform_real_distribution<>(0.05, 1.)(gen); // 0.15
        
        bool is_old=true;
        std::string width_type = "Appourchaux";
        // First make a classical MS model... note the n0 to shift frequencies far enough from (to avoid infinite g mode density)
        Setup inputs_MS_aj=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike", 
            cfg_star.Dnu_star, n0 + cfg_star.epsilon_star, cfg_star.delta0l_percent_star*cfg_star.Dnu_star/100.);
        Setup inputs=make_params_RGB_model(inputs_MS_aj, lmax, Nfreqs, Nspline_ref, n0, cfg_star, is_old, width_type);
        // --------------------------------------------
        // Define ranges 
        cfg_star.fmin=(n0-1)*cfg_star.Dnu_star; // WARNING: GETTING CLOSE TO fmin=0 WILL LEAD TO INFINITE COMPUTATION TIME. PUT n0=4 at lowest
        cfg_star.fmax=1000;
        int Ndata=std::ceil((cfg_star.fmax - cfg_star.fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, cfg_star.fmin, cfg_star.fmax);
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new= model_RGB_asympt_aj_AppWidth_HarveyLike_v4_parallel(inputs.params, inputs.params_length, x, output_param);
        //VectorXd sols_new= model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4_ref(inputs.params, inputs.params_length, x, output_param);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in microseconds
        auto start2 = high_resolution_clock::now();
        VectorXd  sols_ref= model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4_ref(inputs.params, inputs.params_length, x, output_param);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in microseconds
        // Check if the solutions are the same within tolerance
        std::cout << " (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
        //std::cout << "(sols_new - sols_ref).sum()  = " << (sols_new - sols_ref).sum()  << std::endl;
        if ((sols_new - sols_ref).norm() > tolerance/100){
            //std::cout << "     (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
    }
 
    // Compare execution times
    std::cout << "New model_RGB_asympt_aj_AppWidth_HarveyLike_v4 total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "Old model_RGB_asympt_aj_AppWidth_HarveyLike_v4 total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
}


BOOST_AUTO_TEST_CASE(model_RGB_asympt_aj_CteWidth_HarveyLike_v4_PARALLELISATION_TEST)
{  
    // This test takes the reference function for the RGB star and introduce its parallelisation 
    // at the model level. Please note that this does not test parallelisation of the ARMM part.
    // This ARMM parallelisation is already included here and is separately tested in unit tests 
    // at github.com/OthmanB/ARMM-solver in the unit_test v1.0--v1.1 branch.
    // Thus, performance gains that are noted here are just for the model paralleliation.
    // An estimate of x1.5 - 2 can be achieved here while the ARMM parallel gives x2 - 5.
    // This leads to an overall total gain of x3 - x10
    const double tolerance = 1e-6; // in percent
    const int numTests_Acc = 10;
    const int numPoints = 20;
    const double resol=1e6/(4.*365.*86400.);
    const std::string cfg_file="../config_ARMM/make_star_test.cfg";
    
    const bool output_param=false;
    const int Nfreqs=4;
    const int Nspline_ref=4;
    const int n0=5;
    
    Eigen::VectorXd nu(numPoints);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform23(2, 3);

    // Load a reference configuration for l=1 mixed modes
    std::unordered_map<std::string, std::string> input_params=readParameterFile(cfg_file);
    Cfg_synthetic_star cfg_star=configure_make_star(input_params); 
    
    std::cout << " ---- " << std::endl;
    std::cout << " model_RGB_asympt_aj_CteWidth_HarveyLike_v4_PARALLELISATION_TEST " << std::endl;
    std::cout << " ---- " << std::endl;

    long duration_total1=0, duration_total2=0;
    std::cout << colors::cyan << "Testing output difference between make_star with and without all Parallelisation of the models (Fully randomized)..." << colors::white << std::endl;
    for (int i = 0; i < numTests_Acc; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        int lmax = Uniform23(gen);
        // p modes
        cfg_star.Dnu_star=std::uniform_real_distribution<>(10., 40.)(gen); //60
        // g modes
        cfg_star.DPl_star=std::uniform_real_distribution<>(100., 400.)(gen); // 400
        cfg_star.q_star=std::uniform_real_distribution<>(0.05, 1.)(gen); // 0.15
        
        bool is_old=true;
        std::string width_type = "Constant";
        // First make a classical MS model... note the n0 to shift frequencies far enough from (to avoid infinite g mode density)
        Setup inputs_MS_aj=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike", 
            cfg_star.Dnu_star, n0 + cfg_star.epsilon_star, cfg_star.delta0l_percent_star*cfg_star.Dnu_star/100.);

        Setup inputs=make_params_RGB_model(inputs_MS_aj, lmax, Nfreqs, Nspline_ref, n0, cfg_star, is_old, width_type, false, true);
        // --------------------------------------------
        // Define ranges 
        cfg_star.fmin=(n0-1)*cfg_star.Dnu_star; // WARNING: GETTING CLOSE TO fmin=0 WILL LEAD TO INFINITE COMPUTATION TIME. PUT n0=4 at lowest
        cfg_star.fmax=1000;
        int Ndata=std::ceil((cfg_star.fmax - cfg_star.fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, cfg_star.fmin, cfg_star.fmax);
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new= model_RGB_asympt_aj_CteWidth_HarveyLike_v4_parallel(inputs.params, inputs.params_length, x, output_param);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in microseconds
        auto start2 = high_resolution_clock::now();
        VectorXd  sols_ref= model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4_ref(inputs.params, inputs.params_length, x, output_param);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in microseconds
        // Check if the solutions are the same within tolerance
        std::cout << " (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
        //std::cout << "(sols_new - sols_ref).sum()  = " << (sols_new - sols_ref).sum()  << std::endl;
        if ((sols_new - sols_ref).norm() > tolerance/100){
            //std::cout << "     (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
    }
 
    // Compare execution times
    std::cout << "model_RGB_asympt_aj_CteWidth_HarveyLike_v4 total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "Old model_RGB_asympt_aj_CteWidth_HarveyLike_v4 total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
}


BOOST_AUTO_TEST_CASE(model_RGB_asympt_aj_CteWidth_HarveyLike_v4_FINAL_TEST)
{  
    // This test takes the reference function for the RGB star and introduce its parallelisation 
    // at the model level. Please note that this does not test parallelisation of the ARMM part.
    // This ARMM parallelisation is already included here and is separately tested in unit tests 
    // at github.com/OthmanB/ARMM-solver in the unit_test v1.0--v1.1 branch.
    // Thus, performance gains that are noted here are just for the model paralleliation.
    // An estimate of x1.5 - 2 can be achieved here while the ARMM parallel gives x2 - 5.
    // This leads to an overall total gain of x3 - x10
    const double tolerance = 1e-6; // in percent
    const int numTests_Acc = 10;
    const int numPoints = 20;
    const double resol=1e6/(4.*365.*86400.);
    const std::string cfg_file="../config_ARMM/make_star_test.cfg";
    
    const bool output_param=false;
    const int Nfreqs=4;
    const int Nspline_ref=4;
    const int n0=5;
    
    Eigen::VectorXd nu(numPoints);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform23(2, 3);

    // Load a reference configuration for l=1 mixed modes
    std::unordered_map<std::string, std::string> input_params=readParameterFile(cfg_file);
    Cfg_synthetic_star cfg_star=configure_make_star(input_params); 
    long duration_total1=0, duration_total2=0;
    std::cout << colors::cyan << "Testing output difference between make_star with and without all Parallelisation of the models (FINAL CODE)..." << colors::white << std::endl;
    for (int i = 0; i < numTests_Acc; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        int lmax = Uniform23(gen);
        // p modes
        cfg_star.Dnu_star=std::uniform_real_distribution<>(10., 40.)(gen); //60
        // g modes
        cfg_star.DPl_star=std::uniform_real_distribution<>(100., 400.)(gen); // 400
        cfg_star.q_star=std::uniform_real_distribution<>(0.05, 1.)(gen); // 0.15
        cfg_star.rot_core_input=0.75; // NEED TO INVESTIGATE WHY THIS HAS TO BE FORCED

        // First make a classical MS model... note the n0 to shift frequencies far enough from (to avoid infinite g mode density)
        Setup inputs_MS_aj=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike", 
            cfg_star.Dnu_star, n0 + cfg_star.epsilon_star, cfg_star.delta0l_percent_star*cfg_star.Dnu_star/100.);

        std::string width_type = "Constant";
        Setup inputs_aj=make_params_RGB_model(inputs_MS_aj, lmax, Nfreqs, Nspline_ref, n0, cfg_star, false, width_type, true, true);
        Setup inputs_ref=make_params_RGB_model(inputs_MS_aj, lmax, Nfreqs, Nspline_ref, n0, cfg_star, true, width_type, true, true);

        // --------------------------------------------
        // Define ranges 
        cfg_star.fmin=(n0-1)*cfg_star.Dnu_star; // WARNING: GETTING CLOSE TO fmin=0 WILL LEAD TO INFINITE COMPUTATION TIME. PUT n0=4 at lowest
        cfg_star.fmax=1000;
        int Ndata=std::ceil((cfg_star.fmax - cfg_star.fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, cfg_star.fmin, cfg_star.fmax);
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new= model_RGB_asympt_aj_CteWidth_HarveyLike_v4(inputs_aj.params, inputs_aj.params_length, x, output_param);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in microseconds
        auto start2 = high_resolution_clock::now();
        VectorXd  sols_ref= model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4_ref(inputs_ref.params, inputs_ref.params_length, x, output_param);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in microseconds
        // Check if the solutions are the same within tolerance
        std::cout << " (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
        //std::cout << "(sols_new - sols_ref).sum()  = " << (sols_new - sols_ref).sum()  << std::endl;
        if ((sols_new - sols_ref).norm() > tolerance/100){
            //std::cout << "     (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
    }
 
    // Compare execution times
    std::cout << "ARMM parallelisation + Model Parallelisation  model_RGB_asympt_aj_CteWidth_HarveyLike_v4 total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "ARMM parallelisation + Model Parallelisation  model_RGB_asympt_aj_CteWidth_HarveyLike_v4 total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
}


BOOST_AUTO_TEST_CASE(model_RGB_asympt_aj_AppWidth_HarveyLike_v4_FINAL_TEST)
{  
    // This test takes the reference function for the RGB star and introduce its parallelisation 
    // at the model level. Please note that this does not test parallelisation of the ARMM part.
    // This ARMM parallelisation is already included here and is separately tested in unit tests 
    // at github.com/OthmanB/ARMM-solver in the unit_test v1.0--v1.1 branch.
    // Thus, performance gains that are noted here are just for the model paralleliation.
    // An estimate of x1.5 - 2 can be achieved here while the ARMM parallel gives x2 - 5.
    // This leads to an overall total gain of x3 - x10
    const double tolerance = 1e-5; // in percent
    const int numTests_Acc = 10;
    const int numPoints = 20;
    const double resol=1e6/(4.*365.*86400.);
    const std::string cfg_file="../config_ARMM/make_star_test.cfg";
    
    const bool output_param=false;
    const int Nfreqs=4;
    const int Nspline_ref=4;
    const int n0=5;
    
    Eigen::VectorXd nu(numPoints);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform23(2, 3);

    // Load a reference configuration for l=1 mixed modes
    std::unordered_map<std::string, std::string> input_params=readParameterFile(cfg_file);
    Cfg_synthetic_star cfg_star=configure_make_star(input_params); 
    long duration_total1=0, duration_total2=0;
    std::cout << colors::cyan << "Testing output difference between make_star with and without all Parallelisation of the models (FINAL CODE)..." << colors::white <<  std::endl;
    for (int i = 0; i < numTests_Acc; i++) { 
        std::cout << colors::yellow << "[" << i << "]" << std::flush;
        int lmax = Uniform23(gen);
        // p modes
        cfg_star.Dnu_star=std::uniform_real_distribution<>(10., 40.)(gen); //60
        // g modes
        cfg_star.DPl_star=std::uniform_real_distribution<>(100., 400.)(gen); // 400
        cfg_star.q_star=std::uniform_real_distribution<>(0.05, 1.)(gen); // 0.15
        cfg_star.rot_core_input=0.75; // NEED TO INVESTIGATE WHY THIS HAS TO BE FORCED
        // First make a classical MS model... note the n0 to shift frequencies far enough from (to avoid infinite g mode density)
        Setup inputs_MS_aj=make_params_aj_model(lmax, Nfreqs, "model_MS_Global_aj_HarveyLike", 
            cfg_star.Dnu_star, n0 + cfg_star.epsilon_star, cfg_star.delta0l_percent_star*cfg_star.Dnu_star/100.);

        std::string width_type = "Appourchaux";
        Setup inputs_aj=make_params_RGB_model(inputs_MS_aj, lmax, Nfreqs, Nspline_ref, n0, cfg_star, false, width_type, true, true);
        Setup inputs_ref=make_params_RGB_model(inputs_MS_aj, lmax, Nfreqs, Nspline_ref, n0, cfg_star, true, width_type, true, true);

        // --------------------------------------------
        // Define ranges 
        cfg_star.fmin=(n0-1)*cfg_star.Dnu_star; // WARNING: GETTING CLOSE TO fmin=0 WILL LEAD TO INFINITE COMPUTATION TIME. PUT n0=4 at lowest
        cfg_star.fmax=1000;
        int Ndata=std::ceil((cfg_star.fmax - cfg_star.fmin)/resol);
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ndata, cfg_star.fmin, cfg_star.fmax);
        auto start1 = high_resolution_clock::now();
        VectorXd sols_new= model_RGB_asympt_aj_AppWidth_HarveyLike_v4(inputs_aj.params, inputs_aj.params_length, x, output_param);
        auto end1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(end1 - start1).count();
        duration_total1=duration_total1 + static_cast<long>(duration1); // result in microseconds
        auto start2 = high_resolution_clock::now();
        VectorXd  sols_ref= model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4_ref(inputs_ref.params, inputs_ref.params_length, x, output_param);
        auto end2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end2 - start2).count();
        duration_total2=duration_total2 + static_cast<long>(duration2); // result in microseconds
        // Check if the solutions are the same within tolerance
        std::cout << " (sols_new - sols_ref).norm()  = " << (sols_new - sols_ref).norm()  << std::endl;
        if ((sols_new - sols_ref).norm() > tolerance/100){
            BOOST_FAIL("Error: delta greater than the specified tolerance="+dbl_to_str(tolerance)+" % ");
        }
    }
 
    // Compare execution times
    std::cout << "ARMM parallelisation + Model Parallelisation model_RGB_asympt_aj_AppWidth_HarveyLike_v4 total execution time: " << duration_total1/1e6 << " seconds" << std::endl;
    std::cout << "ARMM parallelisation model_RGB_asympt_aj_AppWidth_HarveyLike_v4 total execution time: " << duration_total2/1e6 << " seconds" << std::endl;
}

Setup make_params_RGB_model(Setup inputs_ref, const int lmax, const int Nfreqs, const int Nspline_ref, const int n0,
        const Cfg_synthetic_star cfg_star, const bool is_old, const std::string width_type, const bool fix_compare, const bool exclude_nobias_test){
    // Function designed to generate a parameter asnd param_length vector that is suitable for evolved models
    // lmax : maximum l
    // Nfreqs: number of l=0 frequencies in the test
    // Nspline_ref : number of spline nodes in the test
    // n0 : first radial order. Make sure to have it at least n0>4 to avoid extremely large number of g modes (even infinite for n0=0)
    // cfg_star: A reference configuration file to fix many of the parameters.
    // is_old: If true, it will create a parameter vector suitable for models a1etaa3 insteaod of aj
    // width_type: Either Appourchaux Widths (5 parameters), Constant (1 parameter) or Free (NOT YET IMPLEMENTED)
    // fix_compare: if set to true, it will not generate random values but use a fix set of values. Important to make comparisons
    //              between models of different kinds but expected to give the same results
    //              Whenever possible, the fix value will be taken from cfg_star


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uniform01(0, 1);
    std::uniform_int_distribution<int> Uniform02(0, 2);
    std::uniform_int_distribution<int> Uniform12(1, 2);
    
    int Ncopy0, Ncopy1, Ncopy0_ref, Ncopy1_ref;

    //std::cout << colors::red << "inputs_ref.params :" << inputs_ref.params.transpose() << std::endl;
    //std::cout << "inputs_ref.params_length :" << inputs_ref.params_length.transpose() << colors::white << std::endl;    

    // Generate model vectors with appropriate sizes and set params_length
    Setup inputs;
    int Nsplit;
    if (is_old == true){
        Nsplit=5;
    } else{
        Nsplit=10;
    }
    //    Nmax , lmax , Nfl0 ,8 + Nfl1bias , Nfl2 , Nfl3 , Nsplit= 5 || 1 + 12 + 2,  Nwidth, Nnoise, Ninc, Ncfg
    //const int New_size=inputs_ref.params_length[0] + inputs_ref.params_length[1] + inputs_ref.params_length[2] + 
    //    8 + Nspline_ref + inputs_ref.params_length[4] + inputs_ref.params_length[5] + Nsplit + inputs_ref.params_length[7] +
    //    inputs_ref.params_length[8] + inputs_ref.params_length[9] + inputs_ref.params_length[10] + 4;
    inputs.params_length=inputs_ref.params_length;
    inputs.params_length[3]=8 + Nspline_ref*2;// d0l, DP, alpha_g, q, Hspread, sigma_g, Wfact, Hfact + Nfl1bias_node + Nfl1bias_err
    inputs.params_length[6]=Nsplit; // if is_old == 1 : 5, else 1 + 9 (rot_env, rot_core, a2_env, a2_core, a3_env, a4_env, a5_env, a6_env, eta_switch, asym)
    if (width_type == "Appourchaux" || width_type == "Constant"){
        if (width_type == "Appourchaux"){ // Appourchaux relation for the width
            inputs.params_length[7]=6; // Appourchaux width models have 6 parameters
        }
        if (width_type == "Constant"){
            inputs.params_length[7]=1;
        }
    } else{ // Models with as many width as l=0 modes IF AND ONLY IF we have "Free"
        if (width_type != "Free"){
            std::cerr << "Error: The width_type must be either: Appourchaux, Constant or Free" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    inputs.params_length[10]=inputs.params_length[10] + 4;
    inputs.params.resize(inputs.params_length.sum());
    inputs.params.setConstant(-9999);
    // Use data for Nmax, lmax, Nfl0 from the aj model...
    Ncopy0=0;
    Ncopy1= inputs.params_length[0] + inputs.params_length[1] + inputs.params_length[2];
    inputs.params.segment(Ncopy0, Ncopy1) = inputs_ref.params.segment(Ncopy0, Ncopy1);
    // Then replace the l=1 frequencies by the l=1 mixed modes parameters
    inputs.params[Ncopy1]= cfg_star.delta0l_percent_star*cfg_star.Dnu_star/100; //params[Nmax + lmax + Nfl0];
    inputs.params[Ncopy1+1]=cfg_star.DPl_star;
    inputs.params[Ncopy1+2]=cfg_star.alpha_g_star;
    inputs.params[Ncopy1+3]=cfg_star.q_star;
    inputs.params[Ncopy1+4]=cfg_star.H0_spread; //SOME MODELS DO NOT HAVE IT. CAREFULL
    inputs.params[Ncopy1+5]=0; // sigma_g_l1 is OBSELETE
    inputs.params[Ncopy1+6]=cfg_star.Wfactor;
    inputs.params[Ncopy1+7]=cfg_star.Hfactor; //std::abs(params[Nmax + lmax + Nfl0 + 7]);
    // Bias spline function nodes are set at l=1 p modes frequencies (Not optimal distrinution obviously but this is a TEST function)
    // + randomized error terms
    int el=1;
    VectorXd fref(Nspline_ref), ferr(Nspline_ref);
    for (int en=0; en<Nspline_ref; en++){ //params[Nmax+ lmax + Nfl0 + 8 + i]; 
        inputs.params[Ncopy1+8+en]=(n0 + en +  cfg_star.epsilon_star + el/2) * cfg_star.Dnu_star + cfg_star.delta0l_percent_star*cfg_star.Dnu_star*el*(el+1)/100;
        if (fix_compare == false){
            inputs.params[Ncopy1+8+ Nspline_ref + en]=std::uniform_real_distribution<>(0., cfg_star.Dnu_star/100)(gen); //params[Nmax+ lmax + Nfl0 + 8 + Nferr + i];
        } else{
            inputs.params[Ncopy1+8+ Nspline_ref + en]=0;
        }
    }
    // Use data for Nfl2, Nfl3... Here, we need to differenciate the starting-ending of inputs from inputs_ref
    Ncopy0= inputs.params_length.segment(0, 4).sum();
    Ncopy1= inputs.params_length.segment(0, 6).sum();
    Ncopy0_ref= inputs_ref.params_length.segment(0, 4).sum();
    Ncopy1_ref= inputs_ref.params_length.segment(0, 6).sum();
   inputs.params.segment(Ncopy0, Ncopy1-Ncopy0) = inputs_ref.params.segment(Ncopy0_ref, Ncopy1_ref-Ncopy0_ref);

    // Now, taking care of the rotation parameters...
    // These models always have the assumption that l>1 probe only the envelope rotation, while l=1 are mixed
    // Depending if this is a old or a new model, we modify the aj parameters
    // The old format is: 2 parameters for rot_core and rot_renv + 1 eta + 1 a3  ( + 1 asym )
    // Instead of: 12 parameters for a1,...,a6 + eta ( + 1 asym )
    Ncopy0=inputs.params_length.segment(0,6).sum();
    Ncopy0_ref= inputs_ref.params_length.segment(0, 6).sum();
    if (is_old == true){ // omega env and a3 same as aj
        if (fix_compare == false){
            inputs.params[Ncopy0]=inputs_ref.params[Ncopy0_ref]; // Omega env //std::abs(params[Nmax + lmax + Nf]);
            inputs.params[Ncopy0+1]=std::uniform_real_distribution<>(0., 1.5)(gen); //  Omega core  //std::abs(params[Nmax + lmax + Nf+1]);
            inputs.params[Ncopy0+2]= inputs_ref.params[Ncopy0_ref+12]; ; // eta
            inputs.params[Ncopy0+3]=inputs_ref.params[Ncopy0_ref+4]; // a3
            inputs.params[Ncopy0+4]=inputs_ref.params[Ncopy0_ref+13]; // asym
        } else{
            inputs.params[Ncopy0]=cfg_star.rot_env_input;
            inputs.params[Ncopy0+1]=cfg_star.rot_core_input;
            inputs.params[Ncopy0+2]=0; // eta
            inputs.params[Ncopy0+3]=0.1; // a3
            inputs.params[Ncopy0+4]=-100; // asym
        }

    } else{ // We just add Omega Core in front of the aj parameters model
        if (fix_compare == false){
            Ncopy1=inputs.params_length.segment(0,6).sum();
            Ncopy1_ref=inputs.params_length.segment(0,6).sum();
            inputs.params[Ncopy0]=std::uniform_real_distribution<>(0., 1.5)(gen); // Omega core //std::abs(params[Nmax + lmax + Nf]);
            inputs.params.segment(Ncopy0+1, Ncopy1-Ncopy0)=inputs_ref.params.segment(Ncopy0_ref, Ncopy1_ref-Ncopy0_ref);
            std::cout << "fix_compare == false with is_old == false NOT TESTED" << std::endl;
            exit(EXIT_SUCCESS);
        } else{
            inputs.params[Ncopy0]=cfg_star.rot_env_input;
            inputs.params[Ncopy0+1]=cfg_star.rot_core_input;
            inputs.params[Ncopy0+2]=0; // a2_env
            inputs.params[Ncopy0+3]=0; // a2_core
            inputs.params[Ncopy0+4]=0.1; // a3_env
            inputs.params[Ncopy0+5]=0; // a4_env
            inputs.params[Ncopy0+6]=0; // a5_env
            inputs.params[Ncopy0+7]=0; // a6_env
            inputs.params[Ncopy0+8]=1; // eta_switch
            inputs.params[Ncopy0+9]=-100; // asym        
        }
    }
    if (width_type == "Appourchaux"){ // Appourchaux relation for the width
        Ncopy0=inputs.params_length.segment(0,7).sum();
     	double numax=inputs.params.segment(0,inputs.params_length[0]).mean();
        inputs.params[Ncopy0]=numax; // Appourchaux width models have 6 parameters
        inputs.params[Ncopy0+1]=numax; // Appourchaux width models have 6 parameters
        inputs.params[Ncopy0+2]=std::abs(4./2150.*numax + (1. - 1000.*4./2150.)); // alpha
        inputs.params[Ncopy0+3]=std::abs(0.8/2150.*numax + (4.5 - 1000.*0.8/2150.)); // Gamma_alpha.
        inputs.params[Ncopy0+4]=std::abs(3400./2150.*numax + (1000. - 1000.*3400./2150.)); // Wdip
        inputs.params[Ncopy0+5]=std::abs(2.8/2200.*numax + (1. - 2.8/2200. * 1.)); //DeltaGammadip
    }
    if (width_type == "Constant"){
        Ncopy0=inputs.params_length.segment(0,7).sum();
        if (fix_compare == false){
            inputs.params[Ncopy0]=std::uniform_real_distribution<>(0.05, 5.)(gen); // This is the width
        } else{
            inputs.params[Ncopy0]=0.5;
        }
    }
    if (width_type == "Free"){
        std::cerr << "width_type == Free not yet set" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Copy the remaining vector with information on Width, noise, inclination, configuration
    Ncopy0=inputs.params_length.segment(0,8).sum();
    Ncopy1=inputs.params_length.sum()-4;
    Ncopy0_ref=inputs_ref.params_length.segment(0,8).sum();
    Ncopy1_ref=inputs_ref.params_length.sum();
    inputs.params.segment(Ncopy0, Ncopy1-Ncopy0)=inputs_ref.params.segment(Ncopy0_ref, Ncopy1_ref-Ncopy0_ref);
    // Add extra parameters at the end
    const double sigma_limit=0; //params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];
    inputs.params[inputs.params.size()-1]=Nspline_ref;
    if (fix_compare == false){
        if (exclude_nobias_test == false){
            inputs.params[inputs.params.size()-2]=Uniform02(gen); // Test the "No bias", "Cubic" or "Hermite" bias function
        } else{
            inputs.params[inputs.params.size()-2]=Uniform12(gen); // Test the "Cubic" or "Hermite" bias function. Use this if the "ref" function does not support "No bias"
        }
        inputs.params[inputs.params.size()-3]=Uniform01(gen); // Test the asymptotic ( == 0 ) and the shifting of l=0 for l=1 p modes ( == 1)
    } else{
        inputs.params[inputs.params.size()-2]=1;
        inputs.params[inputs.params.size()-3]=0;
    }
    inputs.params[inputs.params.size()-4]=sigma_limit;

    //std::cout << " --- " << std::endl; 
    //std::cout << "inputs.params :" << inputs.params.transpose() << std::endl;
    //std::cout << "inputs.params_length :" << inputs.params_length.transpose() << std::endl;   
    //std::cout << " --- " << std::endl; 
    //exit(EXIT_SUCCESS);
    
    return inputs;
}

Setup make_params_aj_model(const int lmax, const int Nfreqs, const std::string model_name, 
        const long double Dnu, const long double epsilon, const long double d0l){
        const double trunc=50;
        const int do_amp=0;
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> Uniform01(0, 1);
        
        // p modes
        VectorXd fl(Nfreqs*(lmax+1));
        VectorXi Nfl(4);
        Nfl.setZero();
        int cpt=0;
        for (int el=0; el<=lmax;el++){
            Nfl[el]=Nfreqs;
            for (int en=0; en<Nfreqs; en++){
                long double scatter_freq=std::uniform_real_distribution<>(-Dnu/100, Dnu/100)(gen);
                fl[cpt]=(en + epsilon + el/2) * Dnu + d0l*el*(el+1) + scatter_freq;
                cpt=cpt+1;
            }
        }
        VectorXd aj;
        if (model_name == "model_MS_Global_aj_HarveyLike"){
            aj.resize(13); // 2 parameters for each aj. j={1,2,3,4,5,6} >> 12 aj parameters + eta
            aj[0]=std::uniform_real_distribution<>(0.1, 5)(gen); // a1
            aj[1]=0;
            aj[2]=std::uniform_real_distribution<>(-0.1*aj[0], 0.1*aj[0])(gen); // a2
            aj[3]=0;
            aj[4]=std::uniform_real_distribution<>(-0.025*aj[0], 0.025*aj[0])(gen); // a3
            aj[5]=0;
            aj[6]=std::uniform_real_distribution<>(-0.025*aj[0], 0.025*aj[0])(gen); // a4
            aj[7]=0;
            aj[8]=std::uniform_real_distribution<>(-0.01*aj[0], 0.01*aj[0])(gen); // a5
            aj[9]=0;
            aj[10]=std::uniform_real_distribution<>(-0.005*aj[0], 0.005*aj[0])(gen); // a6
            aj[11]=0;
            aj[12]=0; // We set eta to 0
        }
        if (model_name == "model_MS_Global_ajAlm_HarveyLike"){
            aj.resize(11); // 2 parameters for each aj. j={1, 3, 5 } + epsilon + theta0 + delta + eta >> 6 aj parameters + 4
            aj[0]=std::uniform_real_distribution<>(0.1, 5)(gen); // a1
            aj[1]=0;
            aj[2]=std::uniform_real_distribution<>(-0.025*aj[0], 0.025*aj[0])(gen); // a3
            aj[3]=0;
            aj[4]=std::uniform_real_distribution<>(-0.01*aj[0], 0.01*aj[0])(gen); // a5
            aj[5]=0;
            aj[6]=std::uniform_real_distribution<>(0.001, 0.005)(gen); // epsilon_nl[0]
            aj[7]=0; // epsilon_nl[1] (slope)
            aj[8]=std::uniform_real_distribution<>(0.0, M_PI/2)(gen); // theta0   
            aj[9]=std::uniform_real_distribution<>(0.0, M_PI/4)(gen); // delta  
            aj[10]=1; // We set eta to 1 
                      
        }

        int sw=Uniform01(gen);
        long double asym;
        if (sw ==1){ // This is a switch to test asym = 0 or asym !=0 cases in the build_l_mode
            asym=std::uniform_real_distribution<>(-100, 100.)(gen);
        } else{
            asym=0;
        }
        VectorXd Visibilities(lmax);
        if (lmax == 2){
            Visibilities << 1.5, 0.53;
        } else{
            if (lmax == 3){
                Visibilities << 1.5, 0.53, 0.07;
            } else{
                std::cout << "Incorrect lmax" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        VectorXd Height(Nfreqs);
        for (int en=0; en<Nfreqs; en++){
            Height[en]=std::uniform_real_distribution<>(10, 20)(gen);
        } 
        VectorXd Width(Nfreqs);
        for (int en=0; en<Nfreqs; en++){
            Width[en]=std::uniform_real_distribution<>(0.5, 2)(gen);
        }   
        VectorXd noise(7);
        noise << 0, 1, 1,   0, 1, 1,  0.1; // White noise only
        long double inc=std::uniform_real_distribution<>(0., 90.)(gen); 
        // FINAL VECTOR        
        VectorXd params(Nfreqs + lmax + Nfreqs*(lmax+1) + aj.size() + 1 + Width.size() + noise.size() + 1 + 3);
        params.setZero();
        params.segment(0,Nfreqs) = Height;
        params.segment(Nfreqs , lmax) = Visibilities;
        params.segment(Nfreqs+lmax , fl.size())=fl;
        params.segment(Nfreqs+lmax+ fl.size() , aj.size())=aj;
        params[Nfreqs+lmax+ fl.size()+ aj.size()]=asym;
        params.segment(Nfreqs + lmax + fl.size() + aj.size() + 1 , Width.size())=Width;
        params.segment(Nfreqs + lmax + fl.size() + aj.size() + 1 + Width.size() , noise.size())=noise;
        params[Nfreqs + lmax + fl.size() + aj.size() + 1 + Width.size() + noise.size()] = inc;
        params[Nfreqs + lmax + fl.size() + aj.size() + 1 + Width.size() + noise.size() + 1]=trunc;
        params[Nfreqs + lmax + fl.size() + aj.size() + 1 + Width.size() + noise.size() + 2]=do_amp;
        //std::cout << "    params = " << params.transpose() << std::endl;
        VectorXi params_length(11);
        params_length << Nfreqs, lmax, Nfl[0], Nfl[1], Nfl[2], Nfl[3], aj.size() + 1, Width.size(), noise.size(), 1, 2;
    
        Setup inputs;
        inputs.params=params;
        inputs.params_length=params_length;
    return inputs;       
}

external_data load_Alm_tables(const std::string grid_dir){
    	// Handling the GSL interpolator pointer
	//const std::string grid_dir="external/Alm/data/Alm_grids_CPP/1deg_grids/"; // Using the provided 1 degree grid
	GridData_Alm_fast grids;
	gsl_funcs funcs_data;
    external_data extra_data;
    try{
		grids=loadAllData(grid_dir, "gate");
		// Pre-initialisation of the grid into gsl : Flattening + gsl init
		funcs_data.flat_grid_A10=flatten_grid(grids.A10);
		funcs_data.flat_grid_A11=flatten_grid(grids.A11);
		funcs_data.flat_grid_A20=flatten_grid(grids.A20);
		funcs_data.flat_grid_A21=flatten_grid(grids.A21);
		funcs_data.flat_grid_A22=flatten_grid(grids.A22);
		funcs_data.flat_grid_A30=flatten_grid(grids.A30);
		funcs_data.flat_grid_A31=flatten_grid(grids.A31);
		funcs_data.flat_grid_A32=flatten_grid(grids.A32);
		funcs_data.flat_grid_A33=flatten_grid(grids.A33);
		funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
		funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
		funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
		funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
		funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
		funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
		funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
		funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
		funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
        extra_data.Alm_interp_gate=funcs_data;
		//
		grids=loadAllData(grid_dir, "triangle");
		// Pre-initialisation of the grid into gsl : Flattening + gsl init
		funcs_data.flat_grid_A10=flatten_grid(grids.A10);
		funcs_data.flat_grid_A11=flatten_grid(grids.A11);
		funcs_data.flat_grid_A20=flatten_grid(grids.A20);
		funcs_data.flat_grid_A21=flatten_grid(grids.A21);
		funcs_data.flat_grid_A22=flatten_grid(grids.A22);
		funcs_data.flat_grid_A30=flatten_grid(grids.A30);
		funcs_data.flat_grid_A31=flatten_grid(grids.A31);
		funcs_data.flat_grid_A32=flatten_grid(grids.A32);
		funcs_data.flat_grid_A33=flatten_grid(grids.A33);
		funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
		funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
		funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
		funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
		funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
		funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
		funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
		funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
		funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
        extra_data.Alm_interp_triangle=funcs_data;
	}
	catch (exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		std::cerr << "       It is possible that the grid files do not exist in the pointed directory" << std::endl;
		std::cerr << "       grid_dir is set to (relative path):" << grid_dir << std::endl;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		std::cerr << "Exception of unknown type in io_ms_global.cpp when attempting to load the Alm grid!\n";
		exit(EXIT_FAILURE);
	}
    return extra_data;
}