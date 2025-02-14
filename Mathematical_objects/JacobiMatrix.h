#pragma once

#include <functional>
#include <array>

class JacobiMatrix {
private:

    // Points arrays.
    static std::array<double, 8> _x;
    static std::array<double, 8> _y;
    static std::array<double, 8> _z;

    // Template functions.
    static inline double const W_(double t) { return 0.5 * (1.0 - t); }
    static inline double const W(double t)  { return 0.5 * (1.0 + t); }

    // Local derivatives. 
    /*
    static inline double dxde(double eps, double eta, double zeta); 
    static inline double dyde(double eps, double eta, double zeta);
    static inline double dzde(double eps, double eta, double zeta);
    
    static inline double dxdn(double eps, double eta, double zeta);
    static inline double dydn(double eps, double eta, double zeta); 
    static inline double dzdn(double eps, double eta, double zeta);
    
    static inline double dxdc(double eps, double eta, double zeta); 
    static inline double dydc(double eps, double eta, double zeta); 
    static inline double dzdc(double eps, double eta, double zeta);
    */

public:

    JacobiMatrix() = delete;
    static void SetValues(std::array<double, 8> x, std::array<double, 8> y, std::array<double, 8> z);
    
    static inline double dxde(double eps, double eta, double zeta);
    static inline double dyde(double eps, double eta, double zeta);
    static inline double dzde(double eps, double eta, double zeta);

    static inline double dxdn(double eps, double eta, double zeta);
    static inline double dydn(double eps, double eta, double zeta);
    static inline double dzdn(double eps, double eta, double zeta);

    static inline double dxdc(double eps, double eta, double zeta);
    static inline double dydc(double eps, double eta, double zeta);
    static inline double dzdc(double eps, double eta, double zeta);

    static std::function<double(double, double, double)> const GetValueAt(size_t i, size_t j);
    static std::function<double(double, double, double)> const GetValueAtInverse(size_t i, size_t j);
    static std::function<double(double, double, double)> const GetValueAtTransposed(size_t i, size_t j);
    static std::function<double(double, double, double)> const GetDeterminant();
};

