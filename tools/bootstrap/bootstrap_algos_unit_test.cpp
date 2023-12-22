#define BOOST_TEST_MODULE BootstrapTests
#include "bootstrap_algos.h"
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(BlockBootstrapTest) {
    // Test case for the block bootstrap function
    VectorXd originalSamples(10);
    originalSamples << 1.2, 3.4, 5.6, 7.8, 9.0, 11.2, 13.4, 15.6, 17.8, 19.0;
    int blockSize = 2;
    // Perform block bootstrap resampling
    VectorXd resampledData = blockBootstrap(originalSamples, blockSize);

    // Ensure that the size of the resampled data is correct
    std::cout << "    -------> " << "Block Bootstrap Test" << std::endl;
    BOOST_TEST(resampledData.size() == originalSamples.size());
    BOOST_TEST(resampledData != originalSamples);

    std::cout << "originalSamples : " << originalSamples.transpose() << std::endl;
    std::cout << "resampledData   : " << resampledData.transpose() << std::endl;
    
}

BOOST_AUTO_TEST_CASE(ClassicalBootstrapTest) {
    // Test case for the classical bootstrap function
    VectorXd originalSamples(10);
    originalSamples << 1.2, 3.4, 5.6, 7.8, 9.0, 11.2, 13.4, 15.6, 17.8, 19.0;
    // Perform classical bootstrap resampling
    VectorXd resampledData = Bootstrap(originalSamples);

    // Ensure that the size of the resampled data is correct
    std::cout << "    -------> " << "Classical Bootstrap Test" << std::endl;
    BOOST_TEST(resampledData.size() == originalSamples.size());
    BOOST_TEST(resampledData != originalSamples);
    std::cout << "originalSamples : " << originalSamples.transpose() << std::endl;
    std::cout << "resampledData   : " << resampledData.transpose() << std::endl;

}

