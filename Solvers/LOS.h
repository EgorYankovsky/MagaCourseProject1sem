#pragma once

#include "Solver.h"

class LOS : public Solver {
public:
    LOS() : Solver() {};
    LOS(double eps, size_t maxIters) : Solver(eps, maxIters) {};
    ~LOS() {};
    void Solve() override {};
};

