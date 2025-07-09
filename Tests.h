#ifndef Tests_h_
#define Tests_h_

#include "config.h"
#include "SWAlgos.h"

class Tests{
    public:
    Tests( std::string seq1, std::string seq2, double& speedup);
    void run();

    private:
    void time();
    bool validate();

    std::string seq1, seq2;
    int maxScoreSerial, maxIndexSerialI, maxIndexSerialJ, maxScoreScan, maxIndexScanI, maxIndexScanJ, cigarScoreSerial, cigarScoreScan;
    std::vector<std::vector<int>> serialMatrix;
    int8_t* scanMatrix;
    int serialMaxScore, scanMaxScore;
    double& speedup;
    std::string serialCigar, scanCigar;
};

#endif