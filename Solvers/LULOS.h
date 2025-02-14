#pragma once

#include "Solver.h"

class LULOS : public Solver {
public:    
    LULOS() : Solver() {};
    LULOS(double eps, size_t maxIters) : Solver(eps, maxIters) {};
    ~LULOS() {};
    void Solve() override {};
};

