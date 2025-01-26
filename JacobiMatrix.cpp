#include "JacobiMatrix.h"

std::array<double, 8> JacobiMatrix::_x = {};
std::array<double, 8> JacobiMatrix::_y = {};
std::array<double, 8> JacobiMatrix::_z = {};

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

void JacobiMatrix::SetValues(std::array<double, 8> x, std::array<double, 8> y, std::array<double, 8> z) {
    _x = x;
    _y = y;
    _z = z;
}

std::function<double(double, double, double)> const JacobiMatrix::GetValueAt(size_t i, size_t j) {
    if (i == 0 and j == 0) return dxde;
    if (i == 0 and j == 1) return dyde;
    if (i == 0 and j == 2) return dzde;

    if (i == 1 and j == 0) return dxdn;
    if (i == 1 and j == 1) return dydn;
    if (i == 1 and j == 2) return dzdn;

    if (i == 2 and j == 0) return dxdc;
    if (i == 2 and j == 1) return dydc;
    if (i == 2 and j == 2) return dzdc;
}

std::function<double(double, double, double)> const JacobiMatrix::GetValueAtTransposed(size_t i, size_t j) { return GetValueAt(j, i); }

std::function<double(double, double, double)> const JacobiMatrix::GetDeterminant() {
    return dxde;
}
