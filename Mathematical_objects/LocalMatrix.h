#pragma once

#include "Matrix.h"
#include "JacobiMatrix.h"
#include "..\Integration\Integration.h"
#include "..\Functions\BasisFunction.h"

typedef JacobiMatrix J;

enum class LMType {
    Stiffness,
    Mass,
    NotStated
};

class LocalMatrix : public Matrix {
private:
    double _koef;
    const size_t _localMatrixSize = 12;
    LMType _matrixType = LMType::NotStated;
    std::array<double, 8> _x{};
    std::array<double, 8> _y{};
    std::array<double, 8> _z{};

    std::array<std::array<double, 12>, 12> _values{};

    void generate();
    void generateG();
    void generateM();

public:
    LocalMatrix() { _koef = 0.0; }
    LocalMatrix(double koef, std::array<double, 8> x, std::array<double, 8> y, std::array<double, 8> z, LMType matrixType) :
        _koef(koef), _matrixType(matrixType), _x(x), _y(y), _z(z) {
        generate();
    }
    ~LocalMatrix() {}
    double operator() (size_t i, size_t j) const override { return _values[i][j]; };
    double& operator() (size_t i, size_t j) override { return _values[i][j]; };
    inline LMType GetMatrixType() const { return _matrixType; }
    size_t getSize() const override { return _localMatrixSize; };
};

std::function<double(double, double, double)> operator* (std::function<double(double, double, double)> f1,
                                                         std::function<double(double, double, double)> f2);

std::function<double(double, double, double)> operator+ (std::function<double(double, double, double)> f1,
                                                         std::function<double(double, double, double)> f2);

std::function<double(double, double, double)> operator/ (std::function<double(double, double, double)> f1,
                                                         std::function<double(double, double, double)> f2);