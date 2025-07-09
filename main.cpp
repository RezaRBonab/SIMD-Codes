#include "Tests.h"
#include <limits>

int main() {
    std::vector<std::string> sequences = readFasta("BRCA1.fasta");
    double totalSpeedup = 0.0;
    double minSpeedup = std::numeric_limits<double>::max();
    double maxSpeedup = 0.0;
    int numTests = 0;
    
    std::cout << "Starting SW-SCAN-SIMD performance tests..." << std::endl;
    std::cout << "Number of sequences loaded: " << sequences.size() << std::endl;
    std::cout << "=======================================" << std::endl;

    for (int i = 0; i < sequences.size(); ++i) {
        for (int j = i+1; j < sequences.size(); ++j) {
            std::cout << "\nRunning test " << (numTests + 1) << ": sequences " << i << " and " << j << std::endl;
            double speedup = 0.0;
            Tests test(sequences[i], sequences[j], speedup);
            test.run();
            
            // Update statistics
            totalSpeedup += speedup;
            minSpeedup = std::min(minSpeedup, speedup);
            maxSpeedup = std::max(maxSpeedup, speedup);
            numTests++;
            
            std::cout << "Test " << numTests << " completed with speedup: " << speedup << "x" << std::endl;
        }
    }
    
    std::cout << "\n=======================================" << std::endl;
    std::cout << "All tests completed!" << std::endl;
    std::cout << "=======================================" << std::endl;
    std::cout << "Total number of tests: " << numTests << std::endl;
    std::cout << "Total number of sequences: " << sequences.size() << std::endl;
    
    if (numTests > 0) {
        std::cout << "Average Speedup: " << (totalSpeedup / numTests) << "x" << std::endl;
        std::cout << "Minimum Speedup: " << minSpeedup << "x" << std::endl;
        std::cout << "Maximum Speedup: " << maxSpeedup << "x" << std::endl;
    } else {
        std::cout << "No tests were run." << std::endl;
    }

    return 0;
}