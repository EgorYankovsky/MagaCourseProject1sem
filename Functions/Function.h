#pragma once

#include <array>

typedef std::array<double, 3> vector;

class Function {
public:
    Function() = delete;

    // F = (2.0, 2.0, 2.0).
    static vector TestF0(double t0, double t1, double t2, double time) { return vector{ 2.0, 2.0, 2.0 }; }

    // F = (y, z, x).
    static vector TestF1(double t0, double t1, double t2, double time) { return vector{ t1, t2, t0 }; }
};

