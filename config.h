#ifndef CONFIG_H_
#define CONFIG_H_

#define SW_ZERO 0
#define SW_INS 1
#define SW_DEL 2
#define SW_DIAG 3
#define SW_DIAG_E 4
#define SW_INS_E 5
#define SW_DIAG_F 6
#define SW_DEL_F 7

#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <string>
#include <algorithm>
#include <smmintrin.h>
#include <stdint.h>
#include <cmath>
#include <sstream>

extern const int gapOpen;
extern const int gapExtend;
extern const int matchScore;
extern const int misMatchScore;
extern const int SIMD_WIDTH;

int8_t*** makeQueryProfile(const std::string& query, int segLen);
void freeQueryProfile(int8_t*** profile, int numResidues, int segLen);
int getResidueIndex(char c);
std::vector<std::string> readFasta(std::string filename);
#endif