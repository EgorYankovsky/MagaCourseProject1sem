#include "Integration.h"

const std::array<double, 2> Integration::t2 = { 0.577'350'2692,
                                               -0.577'350'2692 };

const std::array<double, 2> Integration::tau2 = { 1.0,
                                                  1.0 };

const std::array<double, 3> Integration::t3{ 0.0,
                                             0.774'596'669'24,
                                            -0.774'596'669'24 };
const std::array<double, 3> Integration::tau3{ 8.0 / 9.0,
                                               5.0 / 9.0,
                                               5.0 / 9.0 };

const std::array<double, 4> Integration::t4{ -0.861'136'311'6,
                                             -0.339'981'043'6,
                                              0.339'981'043'6,
                                              0.861'136'311'6 };
const std::array<double, 4> Integration::tau4{ 0.347'854'845'1,
                                               0.652'145'154'9,
                                               0.652'145'154'9,
                                               0.347'854'845'1 };

const std::array<double, 5> Integration::t5{ -0.906'179'845'9,
                                             -0.538'469'310'1,
                                              0.0,
                                              0.538'469'310'1,
                                              0.906'179'845'9 };
const std::array<double, 5> Integration::tau5{ 0.236'926'885'1,
                                               0.478'628'670'5,
                                               0.568'888'888'9,
                                               0.478'628'670'5,
                                               0.236'926'885'1 };

double Integration::Gauss2(function f, double x0, double x1, double y0, double y1, double z0, double z1) {
    double ans(0.0);
    for (size_t k(0); k < 2; ++k) 
        for (size_t j(0); j < 2; ++j) 
            for (size_t i(0); i < 2; ++i) 
                ans += tau2[k] * tau2[j] * tau2[i] * f(t2[k], t2[j], t2[i]);
    return ans;
}

double Integration::Gauss3(function f, double x0, double x1, double y0, double y1, double z0, double z1)
{
    double ans(0.0);
    for (size_t k(0); k < 3; ++k)
        for (size_t j(0); j < 3; ++j)
            for (size_t i(0); i < 3; ++i)
                ans += tau3[k] * tau3[j] * tau3[i] * f(t3[k], t3[j], t3[i]);
    return ans;
}

double Integration::Gauss4(function f, double x0, double x1, double y0, double y1, double z0, double z1)
{
    double ans(0.0);
    for (size_t k(0); k < 4; ++k)
        for (size_t j(0); j < 4; ++j)
            for (size_t i(0); i < 4; ++i)
                ans += tau4[k] * tau4[j] * tau4[i] * f(t4[k], t4[j], t4[i]);
    return ans;
}

double Integration::Gauss5(function f, double x0, double x1, double y0, double y1, double z0, double z1)
{
    double ans(0.0);
    for (size_t k(0); k < 5; ++k)
        for (size_t j(0); j < 5; ++j)
            for (size_t i(0); i < 5; ++i)
                ans += tau5[k] * tau5[j] * tau5[i] * f(t5[k], t5[j], t5[i]);
    return ans;
}
