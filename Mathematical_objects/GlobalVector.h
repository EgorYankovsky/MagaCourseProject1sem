#pragma once

#include "..\Logger\Logger.h"
#include "Vector.h"
#include "LocalVector.h"

#include <vector>
#include <array>

class GlobalVector : public Vector
{
private:
    std::vector<double> _values{};
    size_t getSize() const override { return _values.size(); }
    void addLocalVectorValues(const std::array<size_t, 12> localRibs, const LocalVector& b);

public:

    double operator() (size_t i) const override;
    double& operator() (size_t i) override;
    __declspec(property(get = getSize)) size_t Size;

    GlobalVector();
    GlobalVector(size_t size);
    ~GlobalVector() {}

    void Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points,
        std::vector<std::pair<size_t, size_t>> generatedRibs);
    void CommitBoundaryConditions(std::vector<std::array<size_t, 6>> borderRibs, std::vector<std::array<double, 3>> points, std::vector<std::pair<size_t, size_t>> generatedRibs);

};

