#include "SWAlgos.h"
#include "alignerInfo.hpp"

std::string SWAlgos::SWScanLinear32(std::string &seq1, std::string &seq2, int &maxScore, int &maxScoreI, int &maxScoreJ)
{
    // Allocate traceback table
    int8_t* tracebackTable = (int8_t*)_mm_malloc((seq1.length() + 1) * (seq2.length() + 1) * sizeof(int8_t), 16);
    if (!tracebackTable) {
        std::cerr << "Memory allocation failed for traceback table!" << std::endl;
        exit(1);
    }
    
    // Initialize traceback table to zero
    for (int i = 0; i <= seq1.length(); ++i) {
        for (int j = 0; j <= seq2.length(); ++j) {
            tracebackTable[i * (seq2.length() + 1) + j] = 0;
        }
    }
    
    const int segWidth = 4; // 32 bit integer width
    const int segLen = (seq1.length() + segWidth - 1) / segWidth;
    int8_t ***queryProfile = makeQueryProfile(seq1, segLen); // Create the query profile for the first sequence

    int m = seq1.length();
    int n = seq2.length();
    int gapScore = gapOpen;

    __m128i vNegLimit = _mm_set1_epi32(INT32_MIN / 2);
    __m128i vZero = _mm_setzero_si128();
    __m128i vMaxH = vZero; // Initialize max H to negative limit
    int32_t score = 0;
    __m128i vGap = _mm_set1_epi32(gapScore); 
    __m128i vNegInfFront = vZero;
    vNegInfFront = _mm_insert_epi32(vNegInfFront, INT32_MIN / 2, 0);                                         // Set first element to negative infinity
    __m128i vSegLenXGap = _mm_add_epi32(vNegInfFront, _mm_slli_si128(_mm_set1_epi32(segLen * -1), 4)); // This is used to calculate the gap extension in the F vector shifted left by 4 bytes

    // what should these values be?
    __m128i vTZero = _mm_set1_epi32(SW_ZERO);                      // Zero vector for T
    __m128i vTDiag = _mm_set1_epi32(SW_DIAG);                      // T for diagonal
    __m128i vTInsE = _mm_set1_epi32(SW_INS);               // T for insertion (linear gap)
    __m128i vTDelF = _mm_set1_epi32(SW_DEL);    // T for deletion (linear gap)

    // Working variables for comparisons
    __m128i case1;
    
    // heap variables
    __m128i *pvE = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16);
    __m128i *pvF = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16);
    __m128i *pvHt = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16);
    __m128i *pvH = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16);
    __m128i *pvHMax = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16);
    __m128i *pvGapper = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16);
    __m128i *pvET = (__m128i *)_mm_malloc(segLen * sizeof(__m128i), 16); // Store E traceback values
    if (!pvE || !pvF || !pvHt || !pvH || !pvHMax || !pvGapper || !pvET)
    {
        std::cerr << "Memory allocation failed!" << std::endl;
        return ""; // Handle memory allocation failure
    }

    for (int i = 0; i < segLen; i++)
    {
        _mm_store_si128(pvE + i, vZero);     // Initialize E to zero
        _mm_store_si128(pvF + i, vZero);     // Initialize F to zero 
        _mm_store_si128(pvH + i, vZero);
        _mm_store_si128(pvHt + i, vZero);
        _mm_store_si128(pvHMax + i, vZero);
        _mm_store_si128(pvET + i, vTZero);  // Initialize E traceback to diagonal opening
    }

    __m128i vGapper = _mm_add_epi32(vZero, vGap);
    for (int i = segLen - 1; i >= 0; i--)
    {
        _mm_store_si128(pvGapper + i, vGapper);
        vGapper = _mm_add_epi32(vGapper, vGap); // Update gap extension for the next segment
    }

    for (int j = 0; j < n; j++)
    {
        __m128i vE;
        __m128i vE_gap;
        __m128i vHt;
        __m128i vF;
        __m128i vF_gap;
        __m128i vH;
        __m128i vHp;
        __m128i vGapper;
        __m128i vT;
        __m128i vET;
        __m128i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm_load_si128(pvH + (segLen - 1));
        vHp = _mm_slli_si128(vHp, 4);
        vHt = vZero;
        vF = vNegLimit;
        for (int i = 0; i < segLen; i++)
        {
            // load the profile for the current segment
            int rdx = getResidueIndex(seq2[j]);
            int segIdx = i;

            // this function sets things in reverse for some reason
            __m128i vProfile = _mm_set_epi32(
                queryProfile[rdx][segIdx][3],
                queryProfile[rdx][segIdx][2],
                queryProfile[rdx][segIdx][1],
                queryProfile[rdx][segIdx][0]);

            // load previous vector values
            vH = _mm_load_si128(pvH + i);
            vE = _mm_load_si128(pvE + i);
            vGapper = _mm_load_si128(pvGapper + i);

            // For linear gap model: E[i][j] = H[i][j-1] + gapOpen
            vE = _mm_add_epi32(vH, vGap); // H + gapOpen (linear gap)
            vE = _mm_max_epi32(vE, vZero); // Ensure no negative values
            
            // Store E traceback - for linear gaps, E traceback is always SW_INS
            _mm_store_si128(pvET + i, vTInsE);
            //vE = vE_opn;
            // calculate F values
            vGapper = _mm_add_epi32(vHt, vGapper);

            vF = _mm_max_epi32(vGapper, vF);

            // this looks right
            // calculate H values
            vHp = _mm_add_epi32(vHp, vProfile); // Add the profile value to the previous H value
            vHp = _mm_max_epi32(vHp, vZero);    // Ensure no negative values

            vHt = _mm_max_epi32(vHp, vE);

            vHt = _mm_max_epi32(vHt, vF); // Get the maximum H value for the current column

            // save computations in previous vectors
            _mm_store_si128((__m128i *)pvH + i, vHp);
            _mm_store_si128((__m128i *)pvHt + i, vHt);
            _mm_store_si128((__m128i *)pvE + i, vE);
            _mm_store_si128((__m128i *)pvF + i, vF); // Store F for use in second pass
            vHp = vH;                                // Update vHp for the next iteration
        }

        // prefix scan on F and H
        vHt = _mm_slli_si128(vHt, 4); // Shift left to make space for the new column
        vGapper = _mm_load_si128(pvGapper);
        vGapper = _mm_add_epi32(vGapper, vHt);
        vF = _mm_max_epi32(vGapper, vF);
        for (int i = 0; i < segWidth - 2; i++)
        {
            __m128i vFt = _mm_slli_si128(vF, 4);
            vFt = _mm_add_epi32(vFt, vSegLenXGap);
            vF = _mm_max_epi32(vFt, vF);
        }

        // calculate max H value after first pass
        vF = _mm_slli_si128(vF, 4); // Shift left to make space for the new column
        vF = _mm_add_epi32(vF, vNegInfFront);
        vH = _mm_max_epi32(vHt, vF);   // this is the max H value for the current column
        vH = _mm_max_epi32(vH, vZero); // Ensure no negative values

        // now do the second pass down the column
        for (int i = 0; i < segLen; i++)
        {
            vHp = _mm_load_si128((__m128i *)pvH + i);  // This is the diagonal value (H[i-1][j-1] + profile)
            vE = _mm_load_si128((__m128i *)pvE + i);   // Load E value from first pass
            
            // For linear gaps, F[i][j] = H[i-1][j] + gapOpen 
            // The F value from prefix scan contains this correctly
            
            // Compute final H following serial logic: max(0, diagonal, E, F)
            // First, determine which path would be taken before clamping to zero
            __m128i vH_unclamped = _mm_max_epi32(vHp, vE);   // max(diagonal, E)
            vH_unclamped = _mm_max_epi32(vH_unclamped, vF);  // max(max(diagonal, E), F)
            
            // Now clamp to zero for the actual H value
            vH = _mm_max_epi32(vH_unclamped, vZero);
            
            // Check if the final H value is zero (this is what determines SW_ZERO in serial version)
            __m128i case0 = _mm_cmpeq_epi32(vH, vZero);  // True if final H == 0
            
            // Determine traceback based on which value was selected for H_unclamped
            // Following the same priority order as the serial version
            __m128i case1 = _mm_cmpeq_epi32(vH_unclamped, vHp);    // Check if H came from diagonal  
            __m128i case2 = _mm_cmpeq_epi32(vH_unclamped, vE);     // Check if H came from E
            __m128i case3 = _mm_cmpeq_epi32(vH_unclamped, vF);     // Check if H came from F

            // Set traceback value based on which path H took, following serial logic priority
            vT = vTZero;  // Start with zero
            
            // Only set non-zero traceback if final H is not zero
            __m128i vT_nonzero = vTZero;
            
            // For linear gaps, use simple traceback values
            vT_nonzero = _mm_blendv_epi8(vT_nonzero, vTDelF, case3);      // Use SW_DEL if H came from F
            vT_nonzero = _mm_blendv_epi8(vT_nonzero, vTInsE, case2);      // Use SW_INS if H came from E  
            vT_nonzero = _mm_blendv_epi8(vT_nonzero, vTDiag, case1);      // Use SW_DIAG if H came from diagonal
            
            // Use the non-zero traceback only if final H is not zero
            vT = _mm_blendv_epi8(vT_nonzero, vTZero, case0);
            
            // Keep zero if H is zero (but this check is redundant since we started with zero)
            
            // Store the traceback value in the traceback table
            int scores[segWidth];
            _mm_storeu_si128((__m128i *)scores, vT);
            for (int lane = 0; lane < segWidth; ++lane) {
                int row = lane * segLen + i + 1; 
                if (row <= m) {
                    tracebackTable[(row*(n+1))+j + 1] = scores[lane];
                }
            }

            _mm_store_si128((__m128i *)pvH + i, vH); // Store the updated H value
            vMaxH = _mm_max_epi32(vMaxH, vH);        // Update the maximum H value found so far

            // SCAN: Only update maxScore and its position, do not fill the scoreMatrix
            int h_values[segWidth];
            _mm_storeu_si128((__m128i*)h_values, vH);
            for (int lane = 0; lane < segWidth; ++lane) {
                int row = lane * segLen + i + 1;  // Correct striped mapping
                if (row <= m) {
                    int score = h_values[lane];
                    if (score > maxScore) {
                        maxScore = score;
                        maxScoreI = row; 
                        maxScoreJ = j + 1;
                    }
                }
            } 
        }
    }


    // Clean up memory
    freeQueryProfile(queryProfile, 26, segLen);
    _mm_free(pvE);
    _mm_free(pvF);
    _mm_free(pvHt);
    _mm_free(pvH);
    _mm_free(pvHMax);
    _mm_free(pvGapper);
    _mm_free(pvET);

    // Perform traceback and CIGAR string generation
    alignerInfo alignment;
    int i = maxScoreI;
    int j = maxScoreJ;
    int softClipEnd = seq1.length() - i;
    if (softClipEnd > 0)
        alignment.newOp('S', softClipEnd);

    std::cout << "Traceback Table Scan:" << std::endl;
    for (int row = 0; row <= m; ++row)
    {
        for (int col = 0; col <= n; ++col)
        {
            std::cout << (int)tracebackTable[(row * (n + 1)) + col] << " ";
        }
        std::cout << std::endl;
    }

    while (i > 0 && j > 0)
    {
        int tracebackValue = tracebackTable[(i * (n + 1)) + j]; // Access the traceback value for the current position

        // If we hit a zero score, we stop the traceback
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
    
    // Free traceback table
    _mm_free(tracebackTable);
    
    return cigar; // Return a message indicating completion
}
