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
    void addLocalVectorValues(const std::array<size_t, 12> localRibs, const LocalVector& b);

public:
    size_t getSize() const override { return _values.size(); }

    double operator() (size_t i) const override;
    double& operator() (size_t i) override;
    __declspec(property(get = getSize)) size_t Size;

    GlobalVector();
    GlobalVector(size_t size);
    ~GlobalVector() {}

    void Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points,
        std::vector<std::pair<size_t, size_t>> generatedRibs);
    void CommitBoundaryConditions(std::vector<std::array<size_t, 6>> borderRibs, std::vector<std::array<double, 3>> points, std::vector<std::pair<size_t, size_t>> generatedRibs);

    double Norma() const;

    friend double operator* (const GlobalVector v1, const GlobalVector v2);
    friend GlobalVector operator* (const double a, const GlobalVector v);
    friend GlobalVector operator+ (const GlobalVector v1, const GlobalVector v2);
    friend GlobalVector operator- (const GlobalVector v1, const GlobalVector v2);

};

