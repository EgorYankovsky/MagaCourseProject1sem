#pragma once

#include <vector>
#include <array>
#include <algorithm>

#include "Logger/Logger.h"
#include "Matrix.h"
#include "LocalMatrix.h"

class GlobalMatrix : public Matrix {
private:
    size_t _size = 0;

    bool _isPortraitGenerated = false;
    inline bool isPortraitGenerated() const { return _isPortraitGenerated; }

    std::vector<double> _al{};
    std::vector<double> _au{};
    std::vector<double> _di{};

    std::vector<size_t> _ig{};
    std::vector<size_t> _jg{};

    void addLocalMatrixValues(const std::array<size_t, 12> localRibs, const LocalMatrix& G, const LocalMatrix& M);

    double getAlValue(size_t i, size_t j) const;
    double getAuValue(size_t i, size_t j) const;
    

public:
    inline size_t getSize() const override { return _ig.size(); }

    void GeneratePortrait(std::vector<std::array<size_t, 13>> areas, size_t ribsAmount);
    void Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points,
        std::vector<std::pair<size_t, size_t>> generatedRibs, std::vector<std::pair<size_t, std::pair<double, double>>> areasInfo);

    void CommitBoundaryConditions(std::vector<std::array<size_t, 6>> borderRibs);

    __declspec(property(get = isPortraitGenerated)) bool IsPortraitGenerated;
    __declspec(property(get = getSize)) size_t Size;
    
    double getValue(size_t i, size_t j);

    double operator() (size_t i, size_t j) const override;   // getter.
    double& operator() (size_t i, size_t j) override;        // setter.

    GlobalMatrix();
    ~GlobalMatrix();
};