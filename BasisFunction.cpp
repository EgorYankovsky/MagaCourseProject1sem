#include "BasisFunction.h"


std::function<double(double, double, double)> BasisFunction::_phi_0 = [](double t0, double t1, double t2) {return W_(t1) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_1 = [](double t0, double t1, double t2) {return W(t1) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_2 = [](double t0, double t1, double t2) {return W_(t1) * W(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_3 = [](double t0, double t1, double t2) {return W(t1) * W(t2); };

std::function<double(double, double, double)> BasisFunction::_phi_4 = [](double t0, double t1, double t2) {return W_(t0) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_5 = [](double t0, double t1, double t2) {return W(t0) * W_(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_6 = [](double t0, double t1, double t2) {return W_(t0) * W(t2); };
std::function<double(double, double, double)> BasisFunction::_phi_7 = [](double t0, double t1, double t2) {return W(t0) * W(t2); };

std::function<double(double, double, double)> BasisFunction::_phi_8 = [](double t0, double t1, double t2) {return W_(t0) * W_(t1); };
std::function<double(double, double, double)> BasisFunction::_phi_9 = [](double t0, double t1, double t2) {return W(t0) * W_(t1); };
std::function<double(double, double, double)> BasisFunction::_phi_10 = [](double t0, double t1, double t2) {return W_(t0) * W(t1); };
std::function<double(double, double, double)> BasisFunction::_phi_11 = [](double t0, double t1, double t2) {return W(t0) * W(t1); };