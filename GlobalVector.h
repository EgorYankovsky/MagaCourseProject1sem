#pragma once

#include "Logger/Logger.h"
#include "Vector.h"

#include <vector>

class GlobalVector : public Vector
{
private:
    std::vector<double> _values{};
    size_t getSize() const override { return _values.size(); }

public:

    double operator() (size_t i) const override;
    double& operator() (size_t i) override;
    __declspec(property(get = getSize)) size_t Size;

    GlobalVector();
    GlobalVector(size_t size);
    ~GlobalVector() {}
};

