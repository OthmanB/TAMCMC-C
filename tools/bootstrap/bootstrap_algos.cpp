#include "bootstrap_algos.h"

using Eigen::VectorXd;
using namespace std;

// Function to perform block bootstrap resampling
VectorXd blockBootstrap(const VectorXd& samples, int blockSize) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, samples.size() - 1);
    // Seed the random number generator with the current time
   // srand(static_cast<unsigned int>(time(0)));
    int originalSize = samples.size();
    int numBlocks = originalSize / blockSize; // full blocks
    int remainingSamples = originalSize % blockSize; // To deal with the case of originalSize not a multiple of blockSize
    VectorXd resampledData(originalSize);

    // Perform block bootstrap resampling for full blocks
    for (int i = 0; i < numBlocks; ++i) {
        // Select a block
        int blockIndex = rand() % (originalSize - blockSize + 1);
        resampledData.segment(i * blockSize, blockSize) = samples.segment(blockIndex, blockSize);
    }
    // Handling the remaining samples
    if (remainingSamples > 0) {
        // Select a block
        int blockIndex = rand() % (originalSize - blockSize + 1);
        resampledData.tail(remainingSamples) = samples.segment(blockIndex, remainingSamples);
    }
    return resampledData;
}

// Function to perform classical bootstrap resampling
VectorXd Bootstrap(const VectorXd& samples) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, samples.size() - 1);
    //srand(static_cast<unsigned int>(time(0)));
    int originalSize = samples.size();
    VectorXd resampledData(originalSize);

    // Perform classical bootstrap resampling
    for (int i = 0; i < originalSize; ++i) {
        // Generate a random index to select an element (with replacement)
        int randomIndex = rand() % originalSize;
        resampledData[i]=samples[randomIndex];
    }
    return resampledData;
}
