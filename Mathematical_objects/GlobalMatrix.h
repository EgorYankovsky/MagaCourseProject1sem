#pragma once

#include "..\Logger\Logger.h"
#include "Matrix.h"
#include "LocalMatrix.h"
#include "GlobalVector.h"

#include <vector>
#include <array>
#include <algorithm>


class GlobalMatrix : public Matrix {
private:
    size_t _size = 0;

    bool _isPortraitGenerated = false;
    inline bool isPortraitGenerated() const { return _isPortraitGenerated; }

    void addLocalMatrixValues(const std::array<size_t, 12> localRibs, const LocalMatrix& G, const LocalMatrix& M);

    double getAlValue(size_t i, size_t j) const;
    double getAuValue(size_t i, size_t j) const;
    

public:
    inline size_t getSize() const override { return _di.size(); }

    void GeneratePortrait(std::vector<std::array<size_t, 13>> areas, size_t ribsAmount);
    void Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points,
        std::vector<std::pair<size_t, size_t>> generatedRibs, std::vector<std::pair<size_t, std::pair<double, double>>> areasInfo);

    void CommitBoundaryConditions(std::vector<std::array<size_t, 6>> borderRibs);

    std::vector<double> _al{};
    std::vector<double> _au{};
    std::vector<double> _di{};

    std::vector<size_t> _ig{};
    std::vector<size_t> _jg{};

    inline double AL(size_t i) const { return i < _al.size() ? _al[i] : throw "_al argument out of range."; }
    inline double AU(size_t i) const { return i < _au.size() ? _au[i] : throw "_au argument out of range."; }
    inline double DI(size_t i) const { return i < _di.size() ? _di[i] : throw "_di argument out of range."; }

    inline size_t IG(size_t i) const { return i < _ig.size() ? _ig[i] : throw "_ig argument out of range."; }
    inline size_t JG(size_t i) const { return i < _jg.size() ? _jg[i] : throw "_jg argument out of range."; }

    __declspec(property(get = isPortraitGenerated)) bool IsPortraitGenerated;
    __declspec(property(get = getSize)) size_t Size;
    
    double getValue(size_t i, size_t j);

    double operator() (size_t i, size_t j) const override;   // getter.
    double& operator() (size_t i, size_t j) override;        // setter.

    GlobalMatrix();
    ~GlobalMatrix();

    friend GlobalVector operator*(const GlobalMatrix A, const GlobalVector b);
};