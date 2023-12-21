#include <iostream>
#include <vector>
#include <cstdlib> // For rand and srand functions
#include <ctime>   // For time function

using namespace std;

// Function to perform block bootstrap resampling
vector<double> blockBootstrap(const vector<double>& samples, int blockSize) {
    // Seed the random number generator with the current time
    srand(static_cast<unsigned int>(time(0)));

    // Get the size of the original samples
    int originalSize = samples.size();

    // Create a vector to store the resampled data
    vector<double> resampledData;

    // Perform block bootstrap resampling
    for (int i = 0; i < originalSize; i += blockSize) {
        // Generate a random index to select a block
        int blockIndex = rand() % (originalSize - blockSize + 1);

        // Append the selected block to the resampled data
        for (int j = 0; j < blockSize; ++j) {
            resampledData.push_back(samples[blockIndex + j]);
        }
    }

    return resampledData;
}

// Function to perform classical bootstrap resampling
std::vector<double> Bootstrap(const std::vector<double>& samples) {
    srand(static_cast<unsigned int>(time(0)));
    int originalSize = samples.size();
    std::vector<double> resampledData;

    // Perform classical bootstrap resampling
    for (int i = 0; i < originalSize; ++i) {
        // Generate a random index to select an element (with replacement)
        int randomIndex = rand() % originalSize;
        resampledData.push_back(samples[randomIndex]);
    }

    return resampledData;
}
