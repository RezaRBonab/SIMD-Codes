#include "SWAlgos.h"
#include "alignerInfo.hpp"

std::string SWAlgos::SWSerialAffine(std::string &seq1, std::string &seq2,
                              int &maxScore, int &maxScoreI, int &maxScoreJ){

    int m = seq1.size();
    int n = seq2.size();
    
    // DP matrices
    std::vector<std::vector<int>> H(m+1, std::vector<int>(n+1, 0));
    std::vector<std::vector<int>> E(m+1, std::vector<int>(n+1, 0));
    std::vector<std::vector<int>> F(m+1, std::vector<int>(n+1, 0));
    // Single trace-bit table
    std::vector<std::vector<uint8_t>> TB(m+1, std::vector<uint8_t>(n+1, 0));

    int maxI = 0, maxJ = 0;

    // Fill DP + trace flags
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int16_t mm = (seq1[i-1] == seq2[j-1] ? matchScore : misMatchScore);

            // Compute E (insertion)
            int e_open = H[i][j-1] + gapOpen;
            int e_ext  = E[i][j-1] + gapExtend;
            E[i][j] = std::max(e_open, e_ext);
            
            // Compute F (deletion)
            int f_open = H[i-1][j] + gapOpen;
            int f_ext  = F[i-1][j] + gapExtend;
            F[i][j] = std::max(f_open, f_ext);

            // Compute H
            int diag = H[i-1][j-1] + mm;
            H[i][j] = std::max({0, diag, E[i][j], F[i][j]});

            // Build trace flag mask
            uint8_t mask = 0;
            if (H[i][j] == 0) {
                mask |= TB_ZERO;
            } else {
                if (H[i][j] == diag)    mask |= TB_H_DIAG;
                if (H[i][j] == E[i][j]) mask |= TB_H_FROME;
                if (H[i][j] == F[i][j]) mask |= TB_H_FROMF;
            }
            if (E[i][j] == e_open) mask |= TB_E_H;
            if (E[i][j] == e_ext)  mask |= TB_E_EXT;
            if (F[i][j] == f_open) mask |= TB_F_H;
            if (F[i][j] == f_ext)  mask |= TB_F_EXT;

            TB[i][j] = mask;

            if (H[i][j] > maxScore) {
                maxScore = H[i][j];
                maxI = i; maxJ = j;
            }
        }
    }

    // Traceback
    alignerInfo aln;
    int sc_end = m - maxI;
    if (sc_end > 0) aln.newOp('S', sc_end);

    std::cout << "Serial max score: " << maxScore
         << " at (" << maxI << "," << maxJ << ")\n";

    int i = maxI, j = maxJ;
    maxScoreI = i;
    maxScoreJ = j;
    maxScore = H[i][j];

    uint8_t bits = TB[i][j];
    int state = (bits & TB_H_DIAG  ? TB_H_DIAG  :
                 bits & TB_H_FROME ? TB_H_FROME :
                 bits & TB_H_FROMF ? TB_H_FROMF : TB_H_DIAG);

    while (i > 0 && j > 0) {
        bits = TB[i][j];
        if (bits & TB_ZERO) break;

        if (state == TB_H_DIAG) {
            if (bits & TB_H_DIAG) {
                if (seq1[i-1] == seq2[j-1]) aln.newOp('=',1);
                else                         aln.newOp('X',1);
                --i; --j;
            }
            else if (bits & TB_H_FROME) state = TB_H_FROME;
            else if (bits & TB_H_FROMF) state = TB_H_FROMF;
            else break;
        }
        else if (state == TB_H_FROME) {
            if      (bits & TB_E_EXT) { aln.newOp('I',1); --j; }
            else if (bits & TB_E_H  ) { aln.newOp('I',1); --j; state = TB_H_DIAG; }
            else break;
        }
        else {
            if      (bits & TB_F_EXT) { aln.newOp('D',1); --i; }
            else if (bits & TB_F_H  ) { aln.newOp('D',1); --i; state = TB_H_DIAG; }
            else break;
        }
    }
    if (i > 0) aln.newOp('S', i);

    std::string cigar = aln.str();

    // Recompute score correctly
    int recomputed = 0;
    int run = 0;
    for (char c : cigar) {
        if (isdigit(c)) {
            run = run * 10 + (c - '0');
        } else {
            switch (c) {
                case '=': recomputed += run * matchScore;        break;
                case 'X': recomputed += run * misMatchScore;     break;
                case 'I':
                case 'D': recomputed += gapOpen + (run-1)*gapExtend; break;
                case 'S': /* no score */                   break;
            }
            run = 0;
        }
    }

    return cigar;
}