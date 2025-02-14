#pragma once

#include "Solver.h"

class LOS : public Solver {
public:
    LOS() : Solver() {};
    LOS(double eps, size_t maxIters) : Solver(eps, maxIters) {};
    ~LOS() {};
    GlobalVector* Solve(const GlobalMatrix& A, const GlobalVector& b) override;
};

