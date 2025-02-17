#pragma once

#include "../Mathematical_objects/MathematicalHeader.h"

class Solver {
public:
    Solver() : _eps(1e-15), _maxIters(10'000), _A(nullptr), _b(nullptr) {};
    Solver(double eps, size_t maxIters) : _eps(eps), _maxIters(maxIters), _A(nullptr), _b(nullptr) {}
    ~Solver() {};
    virtual GlobalVector Solve(const GlobalMatrix& A, const GlobalVector& b) const = 0;
protected:
    GlobalMatrix* _A;
    GlobalVector* _b;
    double _eps;
    size_t _maxIters;
};