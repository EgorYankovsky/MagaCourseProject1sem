#include "BasisFunction.h"


std::function<double(double, double, double)> BasisFunction::_phi_0 = [](double t0, double t1, double t2) {return W_(t1) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_1 = [](double t0, double t1, double t2) {return W (t1) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_2 = [](double t0, double t1, double t2) {return W_(t1) * W (t2); };
std::function<double(double, double, double)> BasisFunction::_phi_3 = [](double t0, double t1, double t2) {return W (t1) * W (t2); };

std::function<double(double, double, double)> BasisFunction::_phi_4 = [](double t0, double t1, double t2) {return W_(t0) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_5 = [](double t0, double t1, double t2) {return W (t0) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_6 = [](double t0, double t1, double t2) {return W_(t0) * W (t2); };
std::function<double(double, double, double)> BasisFunction::_phi_7 = [](double t0, double t1, double t2) {return W (t0) * W (t2); };

std::function<double(double, double, double)> BasisFunction::_phi_8 = [](double t0, double t1, double t2)  {return W_(t0) * W_(t1); };
std::function<double(double, double, double)> BasisFunction::_phi_9 = [](double t0, double t1, double t2)  {return W (t0) * W_(t1); };
std::function<double(double, double, double)> BasisFunction::_phi_10 = [](double t0, double t1, double t2) {return W_(t0) * W (t1); };
std::function<double(double, double, double)> BasisFunction::_phi_11 = [](double t0, double t1, double t2) {return W (t0) * W (t1); };




std::array<double, 3> BasisFunction::getVectorF(double t0, double t1, double t2, double time, std::array<double, 12> weights) {
    const std::array<size_t, 12> switchV{
                                 0, 3, 8, 11,
                                 1, 2, 9, 10,
                                 4, 5, 6, 7 };
    std::array<double, 3> ans {};
    for (size_t i(0); i < 12; ++i) {
        //ans[0] += 
    }

    return ans;
}