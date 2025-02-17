#pragma once

#include "Vector.h"
#include "LocalMatrix.h"
#include "..\Functions\Function.h"

class LocalVector : public Vector {
public:
    LocalVector() {};
    LocalVector(std::array<double, 8> x, std::array<double, 8> y, std::array<double, 8> z) : _x(x), _y(y), _z(z) {
        M = new LocalMatrix(1.0, _x, _y, _z, LMType::Mass);
        generate();
    }
    ~LocalVector() {};
    double operator() (size_t i) const override { return _values[i]; }
    double& operator() (size_t i) override { return _values[i]; }
    size_t getSize() const override { return 12; }
private:
    std::array<double, 8> _x{};
    std::array<double, 8> _y{};
    std::array<double, 8> _z{};
    std::array<double, 12> _values{};
    void generate();
    LocalMatrix* M = nullptr;
};