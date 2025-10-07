#include "LocalVector.h"

typedef Function F;

#include <iostream>

void LocalVector::generate() {

    auto findLen = [](vector v) { return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); };
    auto scMult = [](vector v1, vector v2) { return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]; };

    std::array<vector, 12> diffVec = {
        vector{_x[1] - _x[0], _y[1] - _y[0], _z[1] - _z[0]},
        vector{_x[3] - _x[2], _y[3] - _y[2], _z[3] - _z[2]},
        vector{_x[5] - _x[4], _y[5] - _y[4], _z[5] - _z[4]},
        vector{_x[7] - _x[6], _y[7] - _y[6], _z[7] - _z[6]},

        vector{_x[2] - _x[0], _y[2] - _y[0], _z[2] - _z[0]},
        vector{_x[3] - _x[1], _y[3] - _y[1], _z[3] - _z[1]},
        vector{_x[6] - _x[4], _y[6] - _y[4], _z[6] - _z[4]},
        vector{_x[7] - _x[5], _y[7] - _y[5], _z[7] - _z[5]},

        vector{_x[4] - _x[0], _y[4] - _y[0], _z[4] - _z[0]},
        vector{_x[5] - _x[1], _y[5] - _y[1], _z[5] - _z[1]},
        vector{_x[6] - _x[2], _y[6] - _y[2], _z[6] - _z[2]},
        vector{_x[7] - _x[3], _y[7] - _y[3], _z[7] - _z[3]}
    };

    std::array<double, 12> f{};
#pragma omp parallel sections num_threads(12)
    {
#pragma omp section
        { f[0] = scMult(F::TestF0(_x[0] + 0.5 * diffVec[0][0], _y[0] + 0.5 * diffVec[0][1], _z[0] + 0.5 * diffVec[0][2], 0), diffVec[0]) / findLen(diffVec[0]); }
#pragma omp section
        { f[1] = scMult(F::TestF0(_x[2] + 0.5 * diffVec[1][0], _y[2] + 0.5 * diffVec[1][1], _z[2] + 0.5 * diffVec[1][2], 0), diffVec[1]) / findLen(diffVec[1]); }
#pragma omp section
        { f[2] = scMult(F::TestF0(_x[4] + 0.5 * diffVec[2][0], _y[4] + 0.5 * diffVec[2][1], _z[4] + 0.5 * diffVec[2][2], 0), diffVec[2]) / findLen(diffVec[2]); }
#pragma omp section
        { f[3] = scMult(F::TestF0(_x[6] + 0.5 * diffVec[3][0], _y[6] + 0.5 * diffVec[3][1], _z[6] + 0.5 * diffVec[3][2], 0), diffVec[3]) / findLen(diffVec[3]); }

#pragma omp section
        { f[4] = scMult(F::TestF0(_x[0] + 0.5 * diffVec[4][0], _y[0] + 0.5 * diffVec[4][1], _z[0] + 0.5 * diffVec[4][2], 0), diffVec[4]) / findLen(diffVec[4]); }
#pragma omp section
        { f[5] = scMult(F::TestF0(_x[1] + 0.5 * diffVec[5][0], _y[1] + 0.5 * diffVec[5][1], _z[1] + 0.5 * diffVec[5][2], 0), diffVec[5]) / findLen(diffVec[5]); }
#pragma omp section
        { f[6] = scMult(F::TestF0(_x[4] + 0.5 * diffVec[6][0], _y[4] + 0.5 * diffVec[6][1], _z[4] + 0.5 * diffVec[6][2], 0), diffVec[6]) / findLen(diffVec[6]); }
#pragma omp section
        { f[7] = scMult(F::TestF0(_x[5] + 0.5 * diffVec[7][0], _y[5] + 0.5 * diffVec[7][1], _z[5] + 0.5 * diffVec[7][2], 0), diffVec[7]) / findLen(diffVec[7]); }

#pragma omp section
        { f[8] = scMult(F::TestF0(_x[0] + 0.5 * diffVec[8][0], _y[0] + 0.5 * diffVec[8][1], _z[0] + 0.5 * diffVec[8][2], 0), diffVec[8]) / findLen(diffVec[8]); }
#pragma omp section
        { f[9] = scMult(F::TestF0(_x[1] + 0.5 * diffVec[9][0], _y[1] + 0.5 * diffVec[9][1], _z[1] + 0.5 * diffVec[9][2], 0), diffVec[9]) / findLen(diffVec[9]); }
#pragma omp section
        { f[10] = scMult(F::TestF0(_x[2] + 0.5 * diffVec[10][0], _y[2] + 0.5 * diffVec[10][1], _z[2] + 0.5 * diffVec[10][2], 0), diffVec[10]) / findLen(diffVec[10]); }
#pragma omp section
        { f[11] = scMult(F::TestF0(_x[3] + 0.5 * diffVec[11][0], _y[3] + 0.5 * diffVec[11][1], _z[3] + 0.5 * diffVec[11][2], 0), diffVec[11]) / findLen(diffVec[11]); }
    }

    for (size_t i(0); i < 12; ++i)
        for (size_t j(0); j < 12; ++j)
            _values[i] += M->operator()(i, j) * f[j];
}

void LocalVector::generateNew() {
    J::SetValues(_x, _y, _z);
    for (int i(0); i < 12; ++i) {
        auto f = Function::TestFf(_x, _y, _z);
        _values[i] = Integration::Gauss3(BasisFunction::getAt(i) * (J::GetValueAtInverseNoDet(0, i / 4) * f[0] +
                                                                    J::GetValueAtInverseNoDet(1, i / 4) * f[1] + 
                                                                    J::GetValueAtInverseNoDet(2, i / 4) * f[2]));
    }
}