#include "Tests.h"

Tests::Tests(std::string seq1, std::string seq2, double& speedup)
    : seq1(std::move(seq1)), seq2(std::move(seq2)), speedup(speedup)
{ 
    // Initialize variables
    maxScoreSerial = -1;
    maxIndexSerialI = -1;
    maxIndexSerialJ = -1;
    maxScoreScan = -1;
    maxIndexScanI = -1;
    maxIndexScanJ = -1;
    cigarScoreSerial = -1;
    cigarScoreScan = -1;
}

void Tests::time() {
    
    SWAlgos algo;

    auto start1 = std::chrono::high_resolution_clock::now();
    serialCigar = algo.SWSerialAffine(seq1, seq2, maxScoreSerial, maxIndexSerialI, maxIndexSerialJ);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start1;
    std::cout << "SWSerial time: " << duration1.count() << " seconds\n";

    auto start2 = std::chrono::high_resolution_clock::now();
    scanCigar = algo.SWScanAffine32(seq1, seq2, maxScoreScan, maxIndexScanI, maxIndexScanJ);
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = end2 - start2;
    std::cout << "SWScan time:   " << duration2.count() << " seconds\n";
    std::cout << "Speedup: " << duration1.count() / duration2.count() << "x\n";
    speedup = duration1.count() / duration2.count();
    
}

bool Tests::validate() {
    
    bool valid = true;
    if (serialCigar != scanCigar) {
        std::cout << "CIGAR string mismatch:\n";
        std::cout << "Serial: " << serialCigar << "\n";
        std::cout << "Scan:   " << scanCigar << "\n";
        std::cout << "Max score Serial: " << maxScoreSerial << " at (" << maxIndexSerialI << "," << maxIndexSerialJ << ")\n";
        std::cout << "Max score Scan:   " << maxScoreScan << " at (" << maxIndexScanI << "," << maxIndexScanJ << ")\n";
        valid = false;
        if(maxIndexScanI != maxIndexSerialI || maxIndexScanJ != maxIndexSerialJ) {
            std::cout << "** CIGAR STRINGS DO NOT MATCH BUT MAX INDEXES DO NOT MATCH **\n";
        }
        else{

        exit(0); 
        }
    }
    else {
        std::cout << "CIGAR strings match " << scanCigar << "\n";
      std::cout << "Max score Serial: " << maxScoreSerial << " at (" << maxIndexSerialI << "," << maxIndexSerialJ << ")\n";
       std::cout << "Max score Scan:   " << maxScoreScan << " at (" << maxIndexScanI << "," << maxIndexScanJ << ")\n";
    }
    /*
    std::ifstream file1("tracebackTableScan.txt");
    std::ifstream file2("tracebackTableSerial.txt");

    if (!file1 || !file2) {
        std::cerr << "Error opening one of the files.\n";
        return 1;
    }
    std::string line1, line2;
    int row = 0;
    bool identical = true;

    while (std::getline(file1, line1) && std::getline(file2, line2)) {
        std::istringstream stream1(line1);
        std::istringstream stream2(line2);
        int val1, val2;
        int col = 0;

        while (stream1 >> val1 && stream2 >> val2) {
            if (val1 != val2) {
                std::cout << "(scan != serial) Difference at (" << row << ", " << col << "): "
                          << val1 << " != " << val2 << '\n';
                identical = false;
            }
            col++;
        }

        // Handle unequal row lengths
        if ((stream1 >> val1) || (stream2 >> val2)) {
            std::cout << "Row length mismatch at row " << row << '\n';
            identical = false;
        }

        row++;
    }
    file1.close();
    file2.close();

    if (identical) {
        std::cout << "Matrices are identical.\n";
    }

    file1.open("HMatrixScan.txt");
    file2.open("HMatrixSerial.txt");
    if (!file1 || !file2) {
        std::cerr << "Error opening one of the HMatrix files.\n";
        return 1;
    }
    valid = true;
    row = 0;
    while (std::getline(file1, line1) && std::getline(file2, line2)) {
        std::istringstream stream1(line1);
        std::istringstream stream2(line2);
        int val1, val2;
        int col = 0;

        while (stream1 >> val1 && stream2 >> val2) {
            if (val1 != val2) {
                std::cout << "(scan != serial) HMatrix Difference at (" << row << ", " << col << "): "
                          << val1 << " != " << val2 << '\n';
                valid = false;
            }
            col++;
        }

        // Handle unequal row lengths
        if ((stream1 >> val1) || (stream2 >> val2)) {
            std::cout << "HMatrix Row length mismatch at row " << row << '\n';
            valid = false;
        }

        row++;
    }
    file1.close();
    file2.close();
    if (valid) {
        std::cout << "HMatrixes are identical.\n";
    } else {
        std::cout << "HMatrixes are not identical.\n";
    }
    */

    return valid;
}
void Tests::run() {
    time();
   validate();
}