#include "LOS.h"

GlobalVector LOS::Solve(const GlobalMatrix& A, const GlobalVector& b) const {
    GlobalVector x(b.Size); GlobalVector _x(b.Size);
    GlobalVector r(b.Size); GlobalVector _r(b.Size);
    GlobalVector z(b.Size); GlobalVector _z(b.Size);
    GlobalVector p(b.Size); GlobalVector _p(b.Size);
    double alph(0.0); double beta(0.0);
    r = b - A * x; z = r; p = A * r; size_t iter(0);
    do {
        _x = x; _z = z; _r = r; _p = p;
        alph = (_p * _r) / (_p * _p);
        x = _x + alph * _z; r = _r - alph * _p;
        beta = -1.0 * (_p * (A * r)) / (_p * _p);
        z = r + beta * _z; p = A * r + beta * _p;
        if (iter % 4 == 0) std::cout << iter << ". " << std::scientific << r.Norma() / b.Norma() << std::endl;
        ++iter;
    } while (iter < _maxIters and r.Norma() / b.Norma() >= _eps);
    return x;
}