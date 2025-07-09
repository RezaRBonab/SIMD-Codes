#include "SWAlgos.h"
#include "alignerInfo.hpp"

std::string SWAlgos::SWSerialLinear(std::string &seq1, std::string &seq2,
                              int &maxScore, int &maxScoreI, int &maxScoreJ)
{
    std::vector<std::vector<int>> H = std::vector<std::vector<int>>(seq1.length() + 1, std::vector<int>(seq2.length() + 1, 0));
    int gapScore = gapOpen;
    std::vector<std::vector<int>> tracebackTable = std::vector<std::vector<int>>(seq1.length() + 1, std::vector<int>(seq2.length() + 1, SW_ZERO));
    int m = seq1.length();
    int n = seq2.length();

    maxScore = 0;
    maxScoreI = maxScoreJ = 0;

    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            int matchMismatchScore = (seq1[i - 1] == seq2[j - 1]) ? matchScore : misMatchScore;

            // Main recurrence
            int diag = H[i - 1][j - 1] + matchMismatchScore;
            int E = H[i][j - 1] + gapScore;
            int F = H[i - 1][j] + gapScore;
            H[i][j] = std::max({0, diag, E, F});
            
            // Update traceback table based on which value was selected
            if (H[i][j] == 0)
            {
                tracebackTable[i][j] = SW_ZERO; // No alignment
            }
            else if (H[i][j] == diag)
            {
                tracebackTable[i][j] = SW_DIAG; // Diagonal move
            }
            else if (H[i][j] ==  E)
            {
                tracebackTable[i][j] = SW_INS; // Gap in seq1 (horizontal)
            }
            else if (H[i][j] == F)
            {
                tracebackTable[i][j] = SW_DEL; // Gap in seq2 (vertical)
            }
            else
            {
                tracebackTable[i][j] = -1; // No move (should not happen in valid alignment)
            }
            // Track max score
            if (H[i][j] > maxScore)
            {
                maxScore = H[i][j];
                maxScoreI = i;
                maxScoreJ = j;
            }
        }
    }

    std::cout << "Traceback Table Serial:" << std::endl;
    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            std::cout << tracebackTable[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Traceback to find the aligned sequences
    alignerInfo alignment;

    int i = maxScoreI;
    int j = maxScoreJ;
    maxScore = H[i][j];
    
    int softClipEnd = seq1.length() - i;
    if (softClipEnd > 0)
        alignment.newOp('S', softClipEnd);

    int GAP = gapOpen + gapExtend; 
    while (i > 0 && j > 0)
    {
        int tracebackValue = tracebackTable[i][j];
        //Print for debuging

        
        if (tracebackValue == SW_ZERO)
        {
            break;
        }

        if (tracebackValue == SW_DIAG)
        {
            if (seq1[i - 1] == seq2[j - 1])
            {
                alignment.newOp('=', 1);
            }
            else
            {
                alignment.newOp('X', 1);
            }
            --i;
            --j;
                }
        else if (tracebackValue == SW_INS)
        {
            alignment.newOp('I', 1);  // Insertion in seq1 
            --j;  // Move horizontally
        }
        else if (tracebackValue == SW_DEL)
        {
            alignment.newOp('D', 1);  // Deletion in seq1
            --i;  // Move vertically
        }
        else
        {
            // If we can't determine the path, break
            std::cout << "Warning: Could not determine traceback path at i=" << i << ", j=" << j << std::endl;
            break;
        }
    }
    int softClipStart = i;
    if (softClipStart > 0)
        alignment.newOp('S', softClipStart);
    
    std::string cigar = alignment.str();
    return cigar; // Return the CIGAR string as the result of the alignment
}