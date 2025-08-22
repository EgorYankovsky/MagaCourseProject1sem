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

    std::array<vector, 12> middle = {
        vector{_x[0] + 0.5 *  diffVec[0][0], _y[0] + 0.5 *  diffVec[0][1], _z[0] + 0.5 *  diffVec[0][2]},
        vector{_x[2] + 0.5 *  diffVec[1][0], _y[2] + 0.5 *  diffVec[1][1], _z[2] + 0.5 *  diffVec[1][2]},
        vector{_x[4] + 0.5 *  diffVec[2][0], _y[4] + 0.5 *  diffVec[2][1], _z[4] + 0.5 *  diffVec[2][2]},
        vector{_x[6] + 0.5 *  diffVec[3][0], _y[6] + 0.5 *  diffVec[3][1], _z[6] + 0.5 *  diffVec[3][2]},

        vector{_x[0] + 0.5 *  diffVec[4][0], _y[0] + 0.5 *  diffVec[4][1], _z[0] + 0.5 *  diffVec[4][2]},
        vector{_x[1] + 0.5 *  diffVec[5][0], _y[1] + 0.5 *  diffVec[5][1], _z[1] + 0.5 *  diffVec[5][2]},
        vector{_x[4] + 0.5 *  diffVec[6][0], _y[4] + 0.5 *  diffVec[6][1], _z[4] + 0.5 *  diffVec[6][2]},
        vector{_x[5] + 0.5 *  diffVec[7][0], _y[5] + 0.5 *  diffVec[7][1], _z[5] + 0.5 *  diffVec[7][2]},

        vector{_x[0] + 0.5 *  diffVec[8][0], _y[0] + 0.5 *  diffVec[8][1], _z[0] + 0.5 *  diffVec[8][2]},
        vector{_x[1] + 0.5 *  diffVec[9][0], _y[1] + 0.5 *  diffVec[9][1], _z[1] + 0.5 *  diffVec[9][2]},
        vector{_x[2] + 0.5 * diffVec[10][0], _y[2] + 0.5 * diffVec[10][1], _z[2] + 0.5 * diffVec[10][2]},
        vector{_x[3] + 0.5 * diffVec[11][0], _y[3] + 0.5 * diffVec[11][1], _z[3] + 0.5 * diffVec[11][2]},
    };
    
    std::array<vector, 12> middleValues = {
        F::TestF1(middle[0][0], middle[0][1], middle[0][2], 0),
        F::TestF1(middle[1][0], middle[1][1], middle[1][2], 0),
        F::TestF1(middle[2][0], middle[2][1], middle[2][2], 0),
        F::TestF1(middle[3][0], middle[3][1], middle[3][2], 0),

        F::TestF1(middle[4][0], middle[4][1], middle[4][2], 0),
        F::TestF1(middle[5][0], middle[5][1], middle[5][2], 0),
        F::TestF1(middle[6][0], middle[6][1], middle[6][2], 0),
        F::TestF1(middle[7][0], middle[7][1], middle[7][2], 0),

        F::TestF1( middle[8][0],  middle[8][1],  middle[8][2], 0),
        F::TestF1( middle[9][0],  middle[9][1],  middle[9][2], 0),
        F::TestF1(middle[10][0], middle[10][1], middle[10][2], 0),
        F::TestF1(middle[11][0], middle[11][1], middle[11][2], 0),
    };

    std::array<double, 12> f = {
        scMult( middleValues[0],  diffVec[0]) / findLen( diffVec[0]),
        scMult( middleValues[1],  diffVec[1]) / findLen( diffVec[1]),
        scMult( middleValues[2],  diffVec[2]) / findLen( diffVec[2]),
        scMult( middleValues[3],  diffVec[3]) / findLen( diffVec[3]),

        scMult( middleValues[4],  diffVec[4]) / findLen( diffVec[4]),
        scMult( middleValues[5],  diffVec[5]) / findLen( diffVec[5]),
        scMult( middleValues[6],  diffVec[6]) / findLen( diffVec[6]),
        scMult( middleValues[7],  diffVec[7]) / findLen( diffVec[7]),

        scMult( middleValues[8],  diffVec[8]) / findLen( diffVec[8]),
        scMult( middleValues[9],  diffVec[9]) / findLen( diffVec[9]),
        scMult(middleValues[10], diffVec[10]) / findLen(diffVec[10]),
        scMult(middleValues[11], diffVec[11]) / findLen(diffVec[11])
    };

    /*
    std::array<double, 12> vf = {
        scMult(Function::TestF1(_x[0] + 0.5 * diffVec[0][0], _y[0] + 0.5 * diffVec[0][1], _z[0] + 0.5 * diffVec[0][2], 0.0), diffVec[0]) / findLen(diffVec[0]),
        scMult(Function::TestF1(_x[2] + 0.5 * diffVec[2][0], _y[2] + 0.5 * diffVec[2][1], _z[2] + 0.5 * diffVec[2][2], 0.0), diffVec[2]) / findLen(diffVec[2]),
        scMult(Function::TestF1(_x[4] + 0.5 * diffVec[4][0], _y[4] + 0.5 * diffVec[4][1], _z[4] + 0.5 * diffVec[4][2], 0.0), diffVec[4]) / findLen(diffVec[4]),
        scMult(Function::TestF1(_x[6] + 0.5 * diffVec[6][0], _y[6] + 0.5 * diffVec[6][1], _z[6] + 0.5 * diffVec[6][2], 0.0), diffVec[6]) / findLen(diffVec[6]),

        scMult(Function::TestF1(_x[0] + 0.5 * diffVec[0][0], _y[0] + 0.5 * diffVec[0][1], _z[0] + 0.5 * diffVec[0][2], 0.0), diffVec[0]) / findLen(diffVec[0]),
        scMult(Function::TestF1(_x[1] + 0.5 * diffVec[1][0], _y[1] + 0.5 * diffVec[1][1], _z[1] + 0.5 * diffVec[1][2], 0.0), diffVec[1]) / findLen(diffVec[1]),
        scMult(Function::TestF1(_x[4] + 0.5 * diffVec[4][0], _y[4] + 0.5 * diffVec[4][1], _z[4] + 0.5 * diffVec[4][2], 0.0), diffVec[4]) / findLen(diffVec[4]),
        scMult(Function::TestF1(_x[5] + 0.5 * diffVec[5][0], _y[5] + 0.5 * diffVec[5][1], _z[5] + 0.5 * diffVec[5][2], 0.0), diffVec[5]) / findLen(diffVec[5]),

        scMult(Function::TestF1(_x[0] + 0.5 * diffVec[0][0], _y[0] + 0.5 * diffVec[0][1], _z[0] + 0.5 * diffVec[0][2], 0.0), diffVec[0]) / findLen(diffVec[0]),
        scMult(Function::TestF1(_x[1] + 0.5 * diffVec[1][0], _y[1] + 0.5 * diffVec[1][1], _z[2] + 0.5 * diffVec[1][2], 0.0), diffVec[1]) / findLen(diffVec[1]),
        scMult(Function::TestF1(_x[2] + 0.5 * diffVec[2][0], _y[2] + 0.5 * diffVec[2][1], _z[4] + 0.5 * diffVec[2][2], 0.0), diffVec[2]) / findLen(diffVec[2]),
        scMult(Function::TestF1(_x[3] + 0.5 * diffVec[3][0], _y[3] + 0.5 * diffVec[3][1], _z[6] + 0.5 * diffVec[3][2], 0.0), diffVec[3]) / findLen(diffVec[3]),
    };
    */

    for (size_t i(0); i < 12; ++i)
        for (size_t j(0); j < 12; ++j)
            _values[i] += M->operator()(i, j) * f[j];
}
