#include "JacobiMatrix.h"

double JacobiMatrix::dxde(double eps, double eta, double zeta) {
    return 0.5 * (W_(eta) * W_(zeta) * (_x[1] - _x[0]) +
                   W(eta) * W_(zeta) * (_x[3] - _x[2]) +
                  W_(eta) *  W(zeta) * (_x[5] - _x[4]) +
                   W(eta) *  W(zeta) * (_x[7] - _x[6]));
}

inline double JacobiMatrix::dyde(double eps, double eta, double zeta) {
    return 0.5 * (W_(eps) * W_(zeta) * (_y[1] - _y[0]) +
                   W(eps) * W_(zeta) * (_y[3] - _y[2]) +
                  W_(eps) *  W(zeta) * (_y[5] - _y[4]) +
                   W(eps) *  W(zeta) * (_y[7] - _y[6]));
}

inline double JacobiMatrix::dzde(double eps, double eta, double zeta) {
    return 0.5 * (W_(eps) * W_(eta) * (_z[1] - _z[0]) +
                   W(eps) * W_(eta) * (_z[3] - _z[2]) +
                  W_(eps) *  W(eta) * (_z[5] - _z[4]) +
                   W(eps) *  W(eta) * (_z[7] - _z[6]));
}

inline double JacobiMatrix::dxdn(double eps, double eta, double zeta) {
    return 0.5 * (W_(eta) * W_(zeta) * (_x[2] - _x[0]) +
                   W(eta) * W_(zeta) * (_x[3] - _x[1]) +
                  W_(eta) *  W(zeta) * (_x[6] - _x[4]) +
                   W(eta) *  W(zeta) * (_x[7] - _x[5]));
}

inline double JacobiMatrix::dydn(double eps, double eta, double zeta) {
    return 0.5 * (W_(eps) * W_(zeta) * (_y[2] - _y[0]) +
                   W(eps) * W_(zeta) * (_y[3] - _y[1]) +
                  W_(eps) *  W(zeta) * (_y[6] - _y[4]) +
                   W(eps) *  W(zeta) * (_y[7] - _y[5]));
}

inline double JacobiMatrix::dzdn(double eps, double eta, double zeta) {
    return 0.5 * (W_(eps) * W_(eta) * (_z[2] - _z[0]) +
                   W(eps) * W_(eta) * (_z[3] - _z[1]) +
                  W_(eps) *  W(eta) * (_z[6] - _z[4]) +
                   W(eps) *  W(eta) * (_z[7] - _z[5]));
}

inline double JacobiMatrix::dxdc(double eps, double eta, double zeta) {
    return 0.5 * (W_(eta) * W_(zeta) * (_x[4] - _x[0]) +
                   W(eta) * W_(zeta) * (_x[5] - _x[1]) +
                  W_(eta) *  W(zeta) * (_x[6] - _x[2]) +
                   W(eta) *  W(zeta) * (_x[7] - _x[3]));
}

inline double JacobiMatrix::dydc(double eps, double eta, double zeta) {
    return 0.5 * (W_(eps) * W_(zeta) * (_y[4] - _y[0]) +
                   W(eps) * W_(zeta) * (_y[5] - _y[1]) +
                  W_(eps) *  W(zeta) * (_y[6] - _y[2]) +
                   W(eps) *  W(zeta) * (_y[7] - _y[3]));
}

inline double JacobiMatrix::dzdc(double eps, double eta, double zeta) {
    return 0.5 * (W_(eps) * W_(eta) * (_z[4] - _z[0]) +
                   W(eps) * W_(eta) * (_z[5] - _z[1]) +
                  W_(eps) *  W(eta) * (_z[6] - _z[2]) +
                   W(eps) *  W(eta) * (_z[7] - _z[3]));
}

JacobiMatrix::~JacobiMatrix() {}

std::function<double(double, double, double)> const JacobiMatrix::J(size_t i, size_t j) {
    std::function<double(double, double, double)> f;
    f = dxde;
    
    switch (i, j)
    {
    case (0, 0): //f = dxde;
    default:
        break;
    }
    return f;
}


