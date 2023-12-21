#define BOOST_TEST_MODULE BootstrapTests
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(BlockBootstrapTest) {
    // Test case for the block bootstrap function

    std::vector<double> originalSamples = {1.2, 3.4, 5.6, 7.8, 9.0, 11.2, 13.4, 15.6, 17.8, 19.0};
    int blockSize = 2;
    // Perform block bootstrap resampling
    std::vector<double> resampledData = blockBootstrap(originalSamples, blockSize);

    // Ensure that the size of the resampled data is correct
    BOOST_TEST(resampledData.size() == originalSamples.size());

    BOOST_TEST(resampledData != originalSamples);
}

BOOST_AUTO_TEST_CASE(ClassicalBootstrapTest) {
    // Test case for the classical bootstrap function

    std::vector<double> originalSamples = {1.2, 3.4, 5.6, 7.8, 9.0, 11.2, 13.4, 15.6, 17.8, 19.0};

    // Perform classical bootstrap resampling
    std::vector<double> resampledData = classicalBootstrap(originalSamples);

    // Ensure that the size of the resampled data is correct
    BOOST_TEST(resampledData.size() == originalSamples.size());

    BOOST_TEST(resampledData != originalSamples);
}

// Entry point for the unit test executable
int main(int argc, char* argv[]) {
    return boost::unit_test::unit_test_main(&init_unit_test, argc, argv);
}
