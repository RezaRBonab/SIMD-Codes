#ifndef SWAlgos_h_
#define SWAlgos_h_

#include "config.h"

enum TraceFlag : uint8_t {
    TB_ZERO    = 0x01,  // H == 0 → stop
    TB_H_DIAG  = 0x02,  // H from diag (match/mismatch)
    TB_H_FROME = 0x04,  // H from E (insertion)
    TB_H_FROMF = 0x08,  // H from F (deletion)
    TB_E_H     = 0x10,  // E from H (gap open)
    TB_E_EXT   = 0x20,  // E from E (gap extend)
    TB_F_H     = 0x40,  // F from H (gap open)
    TB_F_EXT   = 0x80   // F from F (gap extend)
};

class SWAlgos {
 public:
    std::string SWSerialAffine(  std::string& seq1, std::string& seq2, int& maxScore, int& maxScoreI, int& maxScoreJ);
    std::string SWSerialLinear( std::string& seq1, std::string& seq2, int& maxScore, int& maxScoreI, int& maxScoreJ);
    std::string SWScanAffine32( std::string& seq1,  std::string& seq2, int& maxScore, int& maxScoreI, int& maxScoreJ);
    std::string SWScanLinear32(  std::string& seq1,  std::string& seq2, int& maxScore, int& maxScoreI, int& maxScoreJ);
};

#endif 