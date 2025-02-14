#pragma once
class Solver {
public:
    Solver() : _eps(1e-15), _maxIters(10'000) {};
    Solver(double eps, size_t maxIters) : _eps(eps), _maxIters(maxIters) {}
    ~Solver() {};
    virtual void Solve() = 0;
protected:
    double _eps;
    size_t _maxIters;
};

