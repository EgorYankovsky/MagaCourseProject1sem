#include "LULOS.h"

void LULOS::PartialLU(GlobalMatrix& M) const {
    for (size_t i(0); i < M.getSize(); ++i) {
        for (size_t j(M.IG(i)); j < M.IG(i + 1); ++j) {
            int jCol = M._jg[j];
            int jk = M._ig[jCol];
            int k = M._ig[i];

            int sdvig = M._jg[M._ig[i]] - M._jg[M._ig[jCol]];
            if (sdvig > 0) jk += sdvig;
            else k -= sdvig;
            double sumL = 0.0;
            double sumU = 0.0;

            int iii(0);
            for (; k < j && jk < M._ig[jCol + 1]; k++, jk++) {
                iii++;
                std::cout << iii << " " << k << " " << jk << std::endl;
                sumL += M._al[k] * M._au[jk];
                sumU += M._au[k] * M._al[jk];
            }
            M._al[j] -= sumL;
            M._au[j] -= sumU;
            M._au[j] /= M._di[jCol];
        }
        double sumD = 0.0;
        for (int j = M._ig[i]; j < M._ig[i + 1]; j++)
            sumD += M._al[j] * M._au[j];
        M._di[i] -= sumD;
    }
}

GlobalVector LULOS::Forward(const GlobalMatrix& M, const GlobalVector& v) const {
    GlobalVector res = v;
    for (size_t i(0); i < M.Size; ++i) {
        for (size_t j(M.IG(i)); j < M.IG(i + 1); ++j)
            res(i) = res(i) - M.AL(j) * res(M.JG(j));
        res(i) = res(i) / M.DI(i);
    }
    return res;
}

GlobalVector LULOS::Backward(const GlobalMatrix& M, const GlobalVector& v) const {
    GlobalVector res = v;
    for (long long i(M.Size - 1); i >= 0; --i)
        for (long long j(M.IG(i + 1) - 1); j >= long long(M.IG(i)); --j)
            res(M.JG(j)) = res(M.JG(j)) - M.AU(j) * res(i);
    return res;
}

GlobalVector LULOS::Solve(const GlobalMatrix& A, const GlobalVector& b) const {
    GlobalVector x(b.Size);
    GlobalVector x_(b.Size);
    
    GlobalMatrix LU = A;

    PartialLU(LU);
    
    bool _a = LU._al == A._al;
    bool _b = LU._au == A._au;
    bool _c = LU._di == A._di;

    GlobalVector r = Forward(LU, b - A * x);
    GlobalVector r0 = r;
    GlobalVector r_(b.Size);

    GlobalVector z = Backward(LU, r);
    GlobalVector z_(b.Size);

    GlobalVector p = Forward(LU, A * z);
    GlobalVector p_(b.Size);

    GlobalVector tmp(b.Size);

    double alph(0.0);
    double beta(0.0);
    size_t iter(0);

    do
    {
        x_ = x;
        z_ = z;
        r_ = r;
        p_ = p;
        alph = (p_ * r_) / (p_ * p_);
        x = x_ + alph * z_;
        r = r_ - alph * p_;
        tmp = Forward(LU, A * Backward(LU, r));
        beta = -1.0 * (p_ * tmp) / (p_ * p_);
        z = Backward(LU, r) + beta * z_;
        p = tmp + beta * p_;
        ++iter;
        std::cout << iter << ". " << (r.Norma() * r.Norma()) / (r0.Norma() * r0.Norma()) << std::endl;
    } while (iter < _maxIters && (r.Norma() * r.Norma()) / (r0.Norma() * r0.Norma()) >= _eps * _eps);

    return x;
}
