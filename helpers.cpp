#include "config.h"

const int gapOpen = -2;
const int gapExtend = -2;
const int matchScore = 2;
const int misMatchScore = -1;
const int SIMD_WIDTH = 4;

int8_t*** makeQueryProfile(const std::string& query, int segLen)
{
    int8_t*** profile = new int8_t**[26];
    int n = query.length();
    for (int r = 0; r < 26; r++)
    {
        profile[r] = new int8_t*[segLen];
        for (int s = 0; s < segLen; s++)
        {
            profile[r][s] = (int8_t*)_mm_malloc(SIMD_WIDTH * sizeof(int8_t), 16);

            for (int k = 0; k < SIMD_WIDTH; k++)
            {
                int query_idx = s + segLen * k;
                if (query_idx < n)
                {
                    char query_char = query[query_idx];
                    char db_char = 'A' + r;
                    profile[r][s][k] = query_char == db_char ? matchScore : misMatchScore;
                }
                else
                {
                    profile[r][s][k] = 0;
                }
            }

        }
    }
    return profile;
}

void freeQueryProfile(int8_t*** profile, int numResidues, int segLen)
{
    for (int r = 0; r < numResidues; r++)
    {
        for (int s = 0; s < segLen; s++)
        {
            _mm_free(profile[r][s]);
        }
        delete[] profile[r];
    }
    delete[] profile;
}

int getResidueIndex(char c)
{
    return c - 'A';
}

std::vector<std::string> readFasta(std::string filename)
{
    // a vector to store the sequences
    std::vector<std::string> sequences;

    // Read Multiple sequences from the input file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        return {}; // Return an empty vector if file can't be opened
    }
    std::string line;
    std::string seq;
    while (std::getline(inputFile, line))
    {
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line[0] == '>')
        {
            if (!seq.empty())
            {
                sequences.push_back(seq);
                seq.clear();
            }
        }
        else
        {
            seq += line;
        }
    }
    sequences.push_back(seq);
    inputFile.close();

    return sequences;
}