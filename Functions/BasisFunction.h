#pragma once

#include <functional>
#include <array>

class BasisFunction {

    static double W_(double t) { return 0.5 * (1.0 - t); }
    static double W(double t) { return 0.5 * (1.0 + t); }

    static std::function<double(double, double, double)> _phi_0; 
    static std::function<double(double, double, double)> _phi_1; 
    static std::function<double(double, double, double)> _phi_2; 
    static std::function<double(double, double, double)> _phi_3; 

    static std::function<double(double, double, double)> _phi_4; 
    static std::function<double(double, double, double)> _phi_5; 
    static std::function<double(double, double, double)> _phi_6; 
    static std::function<double(double, double, double)> _phi_7; 
    
    static std::function<double(double, double, double)> _phi_8; 
    static std::function<double(double, double, double)> _phi_9; 
    static std::function<double(double, double, double)> _phi_10;
    static std::function<double(double, double, double)> _phi_11;

public:
    BasisFunction() = delete;
    
    static std::array<std::function<double(double, double, double)>, 3> getRotAt(size_t i) {
        std::array<std::function<double(double, double, double)>, 3> a;
        switch (i)
        {
        case 0: return { [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return -0.5 * W_(t1); } ,
                         [](double t0, double t1, double t2) {return  0.5 * W_(t2); } };

        case 1: return { [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return -0.5 * W (t1); } ,
                         [](double t0, double t1, double t2) {return -0.5 * W_(t2); } };
        
        case 2: return { [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return  0.5 * W_(t1); } ,
                         [](double t0, double t1, double t2) {return  0.5 * W (t2); } };
        
        case 3: return { [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return  0.5 * W (t1); } ,
                         [](double t0, double t1, double t2) {return -0.5 * W (t2); } };



        case 4: return { [](double t0, double t1, double t2) {return  0.5 * W_(t0); } ,
                         [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return -0.5 * W_(t2); } };

        case 5: return { [](double t0, double t1, double t2) {return  0.5 * W (t0); } ,
                         [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return  0.5 * W_(t2); } };

        case 6: return { [](double t0, double t1, double t2) {return -0.5 * W_(t0); } ,
                         [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return -0.5 * W (t2); } };

        case 7: return { [](double t0, double t1, double t2) {return -0.5 * W (t0); } ,
                         [](double t0, double t1, double t2) {return  0.0; } ,
                         [](double t0, double t1, double t2) {return  0.5 * W (t2); } };



        case 8: return { [](double t0, double t1, double t2) {return -0.5 * W_(t0); } ,
                         [](double t0, double t1, double t2) {return  0.5 * W_(t1); } ,
                         [](double t0, double t1, double t2) {return  0.0; } };

        case 9: return { [](double t0, double t1, double t2) {return -0.5 * W (t0); } ,
                         [](double t0, double t1, double t2) {return -0.5 * W_(t1); } ,
                         [](double t0, double t1, double t2) {return  0.0; } };

        case 10: return { [](double t0, double t1, double t2) {return  0.5 * W_(t0); } ,
                          [](double t0, double t1, double t2) {return  0.5 * W (t1); } ,
                          [](double t0, double t1, double t2) {return  0.0; } };

        case 11: return { [](double t0, double t1, double t2) {return  0.5 * W (t0); } ,
                          [](double t0, double t1, double t2) {return -0.5 * W (t1); } ,
                          [](double t0, double t1, double t2) {return  0.0; } };
        
        default: throw std::exception("Out of index value");
        }
    }

    static std::array < std::function<double(double, double, double)>, 3> get_at(size_t i) {
        switch (i)
        {
        case 0: return std::array<std::function<double(double, double, double)>, 3> {_phi_0, [](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; } };
        case 1: return std::array<std::function<double(double, double, double)>, 3> {_phi_1, [](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; } };
        case 2: return std::array<std::function<double(double, double, double)>, 3> {_phi_2, [](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; } };
        case 3: return std::array<std::function<double(double, double, double)>, 3> {_phi_3, [](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; } };

        case 4: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, _phi_4, [](double, double, double) { return 0.0; } };
        case 5: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, _phi_5, [](double, double, double) { return 0.0; } };
        case 6: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, _phi_6, [](double, double, double) { return 0.0; } };
        case 7: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, _phi_7, [](double, double, double) { return 0.0; } };

        case 8: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; }, _phi_8 };
        case 9: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; }, _phi_9 };
        case 10: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; }, _phi_10 };
        case 11: return std::array<std::function<double(double, double, double)>, 3> {[](double, double, double) { return 0.0; }, [](double, double, double) { return 0.0; }, _phi_11 };
        default: throw std::exception("Out of index value");
        }
    }

    static std::function<double(double, double, double)> getAt (size_t i) {
        switch (i)
        {
        case 0: return _phi_0;
        case 1: return _phi_1;
        case 2: return _phi_2;
        case 3: return _phi_3;

        case 4: return _phi_4;
        case 5: return _phi_5;
        case 6: return _phi_6;
        case 7: return _phi_7;

        case 8: return _phi_8;
        case 9: return _phi_9;
        case 10: return _phi_10;
        case 11: return _phi_11;
        default: throw std::exception("Out of index value");
        }
    }

    static std::array<double, 3> getVectorF(double t0, double t1, double t2, double time, std::array<double, 12> w);
};