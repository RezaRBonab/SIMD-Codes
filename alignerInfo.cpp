#include "alignerInfo.hpp"

alignerInfo::alignerInfo()
    : op_('\0'), count_(0)
{
}

void alignerInfo::newOp(char op, long len)
{
    if (op == op_)
    {
        count_ += len;
    }
    else
    {
        push();
        op_ = op;
        count_ = len;
    }
}


std::string alignerInfo::str()
{
    push(); // flush last run
    return std::string(alignment_.rbegin(), alignment_.rend());
}

void alignerInfo::push()
{
    if (count_ <= 0)
        return;

    std::string countStr(std::to_string(count_));

    alignment_.push_back(op_);

    for (auto it = countStr.rbegin(); it != countStr.rend(); ++it)
        alignment_.push_back(*it);
}