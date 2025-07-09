#include "SWAlgos.h"
#include "alignerInfo.hpp"

std::string SWAlgos::SWScanAffine32(std::string &seq1, std::string &seq2, int &maxScore, int &maxScoreI, int &maxScoreJ)
{
     int8_t* tracebackTable = (int8_t*)_mm_malloc((seq1.length() + 1) * (seq2.length() + 1) * sizeof(int8_t), 16);
    if (!tracebackTable) {
        std::cerr << "Memory allocation failed for scan traceback table!" << std::endl;
        exit(0);
    }
        for (int i = 0; i <= seq1.length(); ++i) {
        for (int j = 0; j <= seq2.length(); ++j) {
            tracebackTable[i * (seq2.length()+ 1) + j] = 0; // Initialize to zero
        }
    }
    const int segWidth = 4; // 32 bit integer width
    const int segLen = (seq1.length() + segWidth - 1) / segWidth;
    int8_t ***queryProfile = makeQueryProfile(seq1, segLen); // Create the query profile for the first sequence

    int m = seq1.length();
    int n = seq2.length();

    __m128i vNegLimit = _mm_set1_epi32(INT32_MIN / 2);
    __m128i vZero = _mm_setzero_si128();
    __m128i vMaxH = vZero; // Initialize max H to negative limit
    __m128i vGapO = _mm_set1_epi32(gapOpen); // Gap open penalty
    __m128i vGapE = _mm_set1_epi32(gapExtend);            // Gap extend penalty
    __m128i vNegInfFront = vZero;
    vNegInfFront = _mm_insert_epi32(vNegInfFront, INT32_MIN / 2, 0);                                         // Set first element to negative infinity
    __m128i vSegLenXGap = _mm_add_epi32(vNegInfFront, _mm_slli_si128(_mm_set1_epi32(segLen * -1), 4)); // This is used to calculate the gap extension in the F vector shifted left by 4 bytes

    // what should these values be?
    __m128i vTZero = _mm_set1_epi32(TB_ZERO);            
    __m128i vTHDiag = _mm_set1_epi32(TB_H_DIAG);                    
    __m128i vTHFromE = _mm_set1_epi32(TB_H_FROME);
    __m128i vTHFromF = _mm_set1_epi32(TB_H_FROMF);
    __m128i vTEH = _mm_set1_epi32(TB_E_H); // E from H (gap open)
    __m128i vTEExt = _mm_set1_epi32(TB_E_EXT); // E from E (gap extend)
    __m128i vTFH = _mm_set1_epi32(TB_F_H); // F from H (gap open)
    __m128i vTFExt = _mm_set1_epi32(TB_F_EXT); // F from F (gap extend)
    
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
        _mm_store_si128(pvET + i, vTZero);  // Initialize E traceback to zero
    }

    __m128i vGapper = _mm_add_epi32(vZero, vGapO);
    for (int i = segLen - 1; i >= 0; i--)
    {
        _mm_store_si128(pvGapper + i, vGapper);
        vGapper = _mm_add_epi32(vGapper, vGapE); // Update gap extension for the next segment
    }

    for (int j = 0; j < n; j++)
    {
        __m128i vE;
        __m128i vE_ext;
        __m128i vE_opn;
        __m128i vHt;
        __m128i vF;
        __m128i vF_ext;
        __m128i vF_opn;
        __m128i vH;
        __m128i vHp;
        __m128i vGapper;
        __m128i vET;
        __m128i vFT;

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

            // calculate new extension and opening values
            vE_opn = _mm_add_epi32(vH, vGapO); // H + gapOpen
            vE_ext = _mm_add_epi32(vE, vGapE);                       // E + gapExtend

            // find E values
            vE = _mm_max_epi32(vE_opn, vE_ext);
            
            // Set E traceback flags based on which values equal the final E
            __m128i e_from_open = _mm_cmpeq_epi32(vE, vE_opn);
            __m128i e_from_ext = _mm_cmpeq_epi32(vE, vE_ext);
            
            vET = vZero;
            vET = _mm_or_si128(vET, _mm_and_si128(e_from_open, vTEH));
            vET = _mm_or_si128(vET, _mm_and_si128(e_from_ext, vTEExt));
            
            // Store E traceback for later use
            _mm_store_si128(pvET + i, vET);

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
            vE = _mm_load_si128((__m128i *)pvE + i);   // Load E value
            vET = _mm_load_si128((__m128i *)pvET + i); // Load E traceback from first pass

            // Recompute F for this position with updated values from prefix scan
            vF_opn = _mm_add_epi32(vH, vGapO); // H + gapOpen
            vF_ext = _mm_add_epi32(vF, vGapE); // F + gapExtend  
            vF = _mm_max_epi32(vF_opn, vF_ext); // Get the maximum F value

            // Set F traceback flags based on which values equal the final F
            __m128i f_from_open = _mm_cmpeq_epi32(vF, vF_opn);
            __m128i f_from_ext = _mm_cmpeq_epi32(vF, vF_ext);
            
            vFT = vZero;
            vFT = _mm_or_si128(vFT, _mm_and_si128(f_from_open, vTFH));
            vFT = _mm_or_si128(vFT, _mm_and_si128(f_from_ext, vTFExt));

            // Compute final H following serial logic: max(0, diagonal, E, F)
            __m128i vH_unclamped = _mm_max_epi32(vHp, vE);   // max(diagonal, E)
            vH_unclamped = _mm_max_epi32(vH_unclamped, vF);  // max(max(diagonal, E), F)
            vH = _mm_max_epi32(vH_unclamped, vZero);        // Clamp to zero for the actual H value
            
            // Build trace flag mask following serial logic
            __m128i vT = vZero;
            
            // Check if H is zero
            __m128i is_zero = _mm_cmpeq_epi32(vH, vZero);
            vT = _mm_blendv_epi8(vT, vTZero, is_zero);
            
            // If H is not zero, add the appropriate H traceback flags
            __m128i not_zero = _mm_xor_si128(is_zero, _mm_set1_epi32(0xFFFFFFFF));
            
            // Check which source H came from
            __m128i h_from_diag = _mm_cmpeq_epi32(vH, vHp);
            __m128i h_from_e = _mm_cmpeq_epi32(vH, vE);
            __m128i h_from_f = _mm_cmpeq_epi32(vH, vF);
            
            // Add H source flags only if H is not zero
            __m128i diag_flag = _mm_and_si128(not_zero, _mm_and_si128(h_from_diag, vTHDiag));
            __m128i e_flag = _mm_and_si128(not_zero, _mm_and_si128(h_from_e, vTHFromE));
            __m128i f_flag = _mm_and_si128(not_zero, _mm_and_si128(h_from_f, vTHFromF));
            
            vT = _mm_or_si128(vT, diag_flag);
            vT = _mm_or_si128(vT, e_flag);
            vT = _mm_or_si128(vT, f_flag);
            
            // Always add E and F traceback information (like in serial version)
            vT = _mm_or_si128(vT, vET);
            vT = _mm_or_si128(vT, vFT);

            // Store the traceback value in the traceback table
            int scores[segWidth];
            _mm_storeu_si128((__m128i *)scores, vT);
            for (int lane = 0; lane < segWidth; ++lane) {
                int row = lane * segLen + i + 1; 
                if (row <= m) {
                    tracebackTable[(row*(n+1))+j + 1] = (int8_t)(scores[lane] & 0xFF);
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

    alignerInfo aln;
    int sc_end = m - maxScoreI;
    if (sc_end > 0) aln.newOp('S', sc_end);

    uint8_t bits = tracebackTable[(maxScoreI * (n + 1)) + maxScoreJ];
    int state = (bits & TB_H_DIAG  ? TB_H_DIAG  :
                 bits & TB_H_FROME ? TB_H_FROME :
                 bits & TB_H_FROMF ? TB_H_FROMF : TB_H_DIAG);


    int i = maxScoreI;
    int j = maxScoreJ;

    while (i > 0 && j > 0) {
        bits = tracebackTable[(i * (n + 1)) + j];
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
    std::cout << "Scan recomputed score: " << recomputed << std::endl;

    return cigar;
}