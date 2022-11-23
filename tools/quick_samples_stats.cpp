#include <Eigen/Dense>
#include <iostream>

double mean_fct(Eigen::VectorXd samples){
   double sum = 0;
   for (int i = 0; i < samples.size(); i++)
      sum = sum + samples[i];
   return sum/samples.size();
}

//calculate median
double median_fct(Eigen::VectorXd samples){
   std::sort(samples.begin(), samples.end());
   if (samples.size() % 2 != 0){
      return samples[samples.size()/2];
   }
   return (samples[(samples.size()-1)/2] + samples[samples.size()/2])/2.0;
}

double stddev_fct(Eigen::VectorXd samples){
    double sum = 0.0;
    const double mean=mean_fct(samples);

    for(int i=0;i<samples.size();i++)
    {
        sum = sum + (samples[i]-mean)*(samples[i]-mean);
    }

    double variance = sum/samples.size();
    return std::sqrt(variance);
}

// For internal test purpose only
// Execute with : g++ quick_samples_stars.cpp -I $EIGEN3_INCLUDE_DIR
/*
int main(int argc, char* argv[]){

   Eigen::VectorXd test(4);
   test << 1.1, 1.2, 1.3, 1.4;

   std::cout << "mean: " << mean_fct(test) << std::endl;
   std::cout << "median: " << median_fct(test) << std::endl;
   std::cout << "stddev: " << stddev_fct(test) << std::endl;
     
}
*/