#ifndef ALIGNER_INFO_HPP
#define ALIGNER_INFO_HPP

#include "config.h"

class alignerInfo
{
public:
    alignerInfo();

    // Flush current run then start (or extend) a new one
    void newOp(char op, long len);

    // After all newOp calls, call this to get the final CIGAR string
    std::string str();

private:
    void push();

    char op_;
    int count_;
    std::vector<char> alignment_;
};

#endif // ALIGNER_INFO_HPP
