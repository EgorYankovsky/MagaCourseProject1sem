#pragma once

#include <array>

typedef double(*function)(double, double, double);

class Integration {
private:

    static const std::array<double, 2> t2;
    static const std::array<double, 2> tau2;
    
    static const std::array<double, 3> t3;
    static const std::array<double, 3> tau3;
    
    static const std::array<double, 4> t4;
    static const std::array<double, 4> tau4;
    
    static const std::array<double, 5> t5;
    static const std::array<double, 5> tau5;
public:
    Integration() = delete;
    static double Gauss2(function f, double x0, double x1, double y0, double y1, double z0, double z1);
    static double Gauss3(function f, double x0, double x1, double y0, double y1, double z0, double z1);
    static double Gauss4(function f, double x0, double x1, double y0, double y1, double z0, double z1);
    static double Gauss5(function f, double x0, double x1, double y0, double y1, double z0, double z1);
};

