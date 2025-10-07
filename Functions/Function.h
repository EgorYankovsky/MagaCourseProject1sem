#pragma once

#include <array>
#include <functional>

typedef std::array<double, 3> vector;
typedef std::array<std::function<double(double, double, double)>, 3> vectorf;

class Function {
public:
    Function() = delete;

    // F = (2.0, 2.0, 2.0).
    static vector TestF0(double t0, double t1, double t2, double time) {
        return vector{ 1.0,
                       1.0,
                       1.0 };
    }

    // F = (y, z, x).
    static vector TestF1(double t0, double t1, double t2, double time) { return vector{ t1, t2, t0 }; }

    static vectorf TestFf(std::array<double, 8> x,
        std::array<double, 8> y,
        std::array<double, 8> z) {
        return vectorf{ Constant(3.0),
                        Constant(2.0),
                        Constant(1.0) };
    }

    static std::function<double(double, double, double)> Constant(double value) {
        std::function<double(double, double, double)> f = [value](double eps, double eta, double zeta) { return value; };
        return f;
    }

    static std::function<double(double, double, double)> X(std::array<double, 8> x) {
        std::function<double(double, double, double)> f = [x](double eps, double eta, double zeta) 
            { return 0.125 * ((1.0 - eps) * (1.0 - eta) * (1.0 - zeta) * x[0] +
                              (1.0 + eps) * (1.0 - eta) * (1.0 - zeta) * x[1] + 
                              (1.0 - eps) * (1.0 + eta) * (1.0 - zeta) * x[2] + 
                              (1.0 + eps) * (1.0 + eta) * (1.0 - zeta) * x[3] + 
                              (1.0 - eps) * (1.0 - eta) * (1.0 + zeta) * x[4] + 
                              (1.0 + eps) * (1.0 - eta) * (1.0 + zeta) * x[5] + 
                              (1.0 - eps) * (1.0 + eta) * (1.0 + zeta) * x[6] + 
                              (1.0 + eps) * (1.0 + eta) * (1.0 + zeta) * x[7]); };
        return f;
    }

    static std::function<double(double, double, double)> Y(std::array<double, 8> y) {
        std::function<double(double, double, double)> f = [y](double eps, double eta, double zeta) 
            { return 0.125 * ((1.0 - eps) * (1.0 - eta) * (1.0 - zeta) * y[0] +
                              (1.0 + eps) * (1.0 - eta) * (1.0 - zeta) * y[1] + 
                              (1.0 - eps) * (1.0 + eta) * (1.0 - zeta) * y[2] + 
                              (1.0 + eps) * (1.0 + eta) * (1.0 - zeta) * y[3] + 
                              (1.0 - eps) * (1.0 - eta) * (1.0 + zeta) * y[4] + 
                              (1.0 + eps) * (1.0 - eta) * (1.0 + zeta) * y[5] + 
                              (1.0 - eps) * (1.0 + eta) * (1.0 + zeta) * y[6] + 
                              (1.0 + eps) * (1.0 + eta) * (1.0 + zeta) * y[7]); };
        return f;
    }

    static std::function<double(double, double, double)> Z(std::array<double, 8> z) {
        std::function<double(double, double, double)> f = [z](double eps, double eta, double zeta) 
            { return 0.125 * ((1.0 - eps) * (1.0 - eta) * (1.0 - zeta) * z[0] +
                              (1.0 + eps) * (1.0 - eta) * (1.0 - zeta) * z[1] + 
                              (1.0 - eps) * (1.0 + eta) * (1.0 - zeta) * z[2] + 
                              (1.0 + eps) * (1.0 + eta) * (1.0 - zeta) * z[3] + 
                              (1.0 - eps) * (1.0 - eta) * (1.0 + zeta) * z[4] + 
                              (1.0 + eps) * (1.0 - eta) * (1.0 + zeta) * z[5] + 
                              (1.0 - eps) * (1.0 + eta) * (1.0 + zeta) * z[6] + 
                              (1.0 + eps) * (1.0 + eta) * (1.0 + zeta) * z[7]); };
        return f;
    }

    // A = (y, z, x).
    static vector TestA(double t0, double t1, double t2, double time) { 
        return vector{ 3.0, 
                       2.0, 
                       1.0 }; 
    }
};