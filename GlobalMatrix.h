#pragma once

#include <vector>
#include <array>
#include <algorithm>

#include "Logger/Logger.h"
#include "Matrix.h"

class GlobalMatrix : public Matrix {
private:
    size_t _size = 0;
    inline size_t getSize() const override { return _size; }

    bool _isPortraitGenerated = false;
    inline bool isPortraitGenerated() const { return _isPortraitGenerated; }

    std::vector<double> _ai{};
    std::vector<double> _aj{};
    std::vector<double> _di{};

    std::vector<size_t> _ig{};
    std::vector<size_t> _jg{};
    

public:
    void GeneratePortrait(std::vector<std::array<size_t, 13>> areas, size_t ribsAmount);
    void Fill();

    __declspec(property(get = isPortraitGenerated)) bool IsPortraitGenerated;
    __declspec(property(get = getSize)) size_t Size;
    
    double operator() (size_t i, size_t j) const override;   // getter.
    double& operator() (size_t i, size_t j) override;        // setter.

    GlobalMatrix();
    ~GlobalMatrix();
};

