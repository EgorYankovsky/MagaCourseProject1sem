#pragma once

#include <functional>
#include <array>

class JacobiMatrix {
private:

    // Points arrays.
    std::array<double, 8> _x{};
    std::array<double, 8> _y{};
    std::array<double, 8> _z{};

    // Template functions.
    inline double const W_(double t) { return 0.5 * (1 - t); }
    inline double const W(double t)  { return 0.5 * (1 + t); }

    // Local derivatives.
    double dxde(double eps, double eta, double zeta); 
    inline double dyde(double eps, double eta, double zeta);
    inline double dzde(double eps, double eta, double zeta);
    
    inline double dxdn(double eps, double eta, double zeta);
    inline double dydn(double eps, double eta, double zeta); 
    inline double dzdn(double eps, double eta, double zeta);
    
    inline double dxdc(double eps, double eta, double zeta); 
    inline double dydc(double eps, double eta, double zeta); 
    inline double dzdc(double eps, double eta, double zeta);

public:

    JacobiMatrix(std::array<double, 8> x,
                 std::array<double, 8> y,
                 std::array<double, 8> z) : _x(x), _y(y), _z(z) {}
    ~JacobiMatrix();

    std::function<double(double, double, double)> const J(size_t i, size_t j);
    std::function<double(double, double, double)> const J_1(size_t i, size_t j);
    std::function<double(double, double, double)> const Jt(size_t i, size_t j);
};

