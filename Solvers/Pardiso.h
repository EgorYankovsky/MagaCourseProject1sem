#pragma once

#include "Solver.h"

class Pardiso : public Solver {
public:
    Pardiso() : Solver() {};
    Pardiso(double eps, size_t maxIters) : Solver(eps, maxIters) {};
    ~Pardiso() {};
    GlobalVector Solve(const GlobalMatrix& A, const GlobalVector& b) const override;
};

