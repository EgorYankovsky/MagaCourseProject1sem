#pragma once

#include <array>
#include <cfloat>

//#include "Function.h"

typedef std::array<double, 3> vector;
typedef std::array<std::array<double, 3>, 3> matrix;
typedef std::array<double, 8> cube;

namespace coordinates_converter {

#pragma region 3DimFunctions
    inline static double_t dxde(double_t var_eps, double_t nu, double_t khi, cube x);
    inline static double_t dxdn(double_t var_eps, double_t nu, double_t khi, cube x);
    inline static double_t dxdc(double_t var_eps, double_t nu, double_t khi, cube x);


    inline static double_t dyde(double_t var_eps, double_t nu, double_t khi, cube y);
    inline static double_t dydn(double_t var_eps, double_t nu, double_t khi, cube y);
    inline static double_t dydc(double_t var_eps, double_t nu, double_t khi, cube y);


    inline static double_t dzde(double_t var_eps, double_t nu, double_t khi, cube z);
    inline static double_t dzdn(double_t var_eps, double_t nu, double_t khi, cube z);
    inline static double_t dzdc(double_t var_eps, double_t nu, double_t khi, cube z);

    inline static double_t F1(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x);

    inline static double_t dF1deps(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x);
    inline static double_t dF1deta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x);
    inline static double_t dF1dzeta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x);

    inline static double_t F2(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y);

    inline static double_t dF2deps(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y);
    inline static double_t dF2deta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y);
    inline static double_t dF2dzeta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y);

    inline static double_t F3(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z);

    inline static double_t dF3deps(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z);
    inline static double_t dF3deta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z);
    inline static double_t dF3dzeta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z);
#pragma endregion

    class Jacobian {
    public:
        static matrix find_inverse(double_t var_eps, double_t nu, double_t khi,
                                   cube cfs_x, cube cfs_y, cube cfs_z);
    };

    static void solve_gauss(matrix& A, vector& b, vector& solution);

    vector convert_from_xyz(double x0, double y0, double z0,
                            const cube& x, const cube& y, const cube& z,
                            double_t _eps = 1e-6, double_t max_iter = 1'000);
}	// namespace coordinates_converter.