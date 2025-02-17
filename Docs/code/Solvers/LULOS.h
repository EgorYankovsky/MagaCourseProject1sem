#pragma once

#include "Solver.h"

class LULOS : public Solver {
private:
    void PartialLU(GlobalMatrix& M) const;
    GlobalVector Forward(const GlobalMatrix& M, const GlobalVector& v) const;
    GlobalVector Backward(const GlobalMatrix& M, const GlobalVector& v) const;
public:    
    LULOS() : Solver() {};
    LULOS(double eps, size_t maxIters) : Solver(eps, maxIters) {};
    ~LULOS() {};
    GlobalVector Solve(const GlobalMatrix& A, const GlobalVector& b) const override;
};

