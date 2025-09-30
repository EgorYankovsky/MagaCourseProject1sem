#include "CoordinatesConverter.h"

inline static double_t coordinates_converter::dxde(double_t var_eps, double_t nu, double_t khi, cube x) {
    return 0.125 * ((1.0 - nu) * (1.0 - khi) * (x[1] - x[0]) + (1.0 + nu) * (1.0 - khi) * (x[3] - x[2]) +
                    (1.0 - nu) * (1.0 + khi) * (x[5] - x[4]) + (1.0 + nu) * (1.0 + khi) * (x[7] - x[6]));}

inline static double_t coordinates_converter::dxdn(double_t var_eps, double_t nu, double_t khi, cube x) {
    return 0.125 * ((1.0 - var_eps) * (1.0 - khi) * (x[2] - x[0]) + (1.0 + var_eps) * (1.0 - khi) * (x[3] - x[1]) +
                    (1.0 - var_eps) * (1.0 + khi) * (x[6] - x[4]) + (1.0 + var_eps) * (1.0 + khi) * (x[7] - x[5])); }

inline static double_t coordinates_converter::dxdc(double_t var_eps, double_t nu, double_t khi, cube x) {
    return 0.125 * ((1.0 - var_eps) * (1.0 - nu) * (x[4] - x[0]) + (1.0 + var_eps) * (1.0 - nu) * (x[5] - x[1]) +
                    (1.0 - var_eps) * (1.0 + nu) * (x[6] - x[2]) + (1.0 + var_eps) * (1.0 + nu) * (x[7] - x[3])); }



inline static double_t coordinates_converter::dyde(double_t var_eps, double_t nu, double_t khi, cube y) {
    return 0.125 * ((1.0 - nu) * (1.0 - khi) * (y[1] - y[0]) + (1.0 + nu) * (1.0 - khi) * (y[3] - y[2]) +
                    (1.0 - nu) * (1.0 + khi) * (y[5] - y[4]) + (1.0 + nu) * (1.0 + khi) * (y[7] - y[6])); }

inline static double_t coordinates_converter::dydn(double_t var_eps, double_t nu, double_t khi, cube y) {
    return 0.125 * ((1.0 - var_eps) * (1.0 - khi) * (y[2] - y[0]) + (1.0 + var_eps) * (1.0 - khi) * (y[3] - y[1]) +
                    (1.0 - var_eps) * (1.0 + khi) * (y[6] - y[4]) + (1.0 + var_eps) * (1.0 + khi) * (y[7] - y[5]));}

inline static double_t coordinates_converter::dydc(double_t var_eps, double_t nu, double_t khi, cube y) {
    return 0.125 * ((1.0 - var_eps) * (1.0 - nu) * (y[4] - y[0]) + (1.0 + var_eps) * (1.0 - nu) * (y[5] - y[1]) +
                    (1.0 - var_eps) * (1.0 + nu) * (y[6] - y[2]) + (1.0 + var_eps) * (1.0 + nu) * (y[7] - y[3])); }



inline static double_t coordinates_converter::dzde(double_t var_eps, double_t nu, double_t khi, cube z) {
    return 0.125 * ((1.0 - nu) * (1.0 - khi) * (z[1] - z[0]) + (1.0 + nu) * (1.0 - khi) * (z[3] - z[2]) +
                    (1.0 - nu) * (1.0 + khi) * (z[5] - z[4]) + (1.0 + nu) * (1.0 + khi) * (z[7] - z[6])); }

inline static double_t coordinates_converter::dzdn(double_t var_eps, double_t nu, double_t khi, cube z) {
    return 0.125 * ((1.0 - var_eps) * (1.0 - khi) * (z[2] - z[0]) + (1.0 + var_eps) * (1.0 - khi) * (z[3] - z[1]) +
                    (1.0 - var_eps) * (1.0 + khi) * (z[6] - z[4]) + (1.0 + var_eps) * (1.0 + khi) * (z[7] - z[5])); }

inline static double_t coordinates_converter::dzdc(double_t var_eps, double_t nu, double_t khi, cube z) {
    return 0.125 * ((1.0 - var_eps) * (1.0 - nu) * (z[4] - z[0]) + (1.0 + var_eps) * (1.0 - nu) * (z[5] - z[1]) +
                    (1.0 - var_eps) * (1.0 + nu) * (z[6] - z[2]) + (1.0 + var_eps) * (1.0 + nu) * (z[7] - z[3])); }

inline static double_t coordinates_converter::F1(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x) {
    return cfs[0] + cfs[1] * var_eps + cfs[2] * nu + cfs[3] * khi +
           cfs[4] * var_eps * nu + cfs[5] * var_eps * khi + cfs[6] * nu * khi +
               cfs[7] * var_eps * nu * khi - 8.0 * _x; }

inline static double_t coordinates_converter::dF1deps(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x) {
    return cfs[1] + cfs[4] * nu + cfs[5] * khi + cfs[7] * nu * khi; }

inline static double_t coordinates_converter::dF1deta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x) {
    return cfs[2] + cfs[4] * var_eps + cfs[6] * khi + cfs[7] * var_eps * khi; }

inline static double_t coordinates_converter::dF1dzeta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _x) {
    return cfs[3] + cfs[5] * var_eps + cfs[6] * nu + cfs[7] * var_eps * nu; }

inline static double_t coordinates_converter::F2(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y) {
    return cfs[0] + cfs[1] * var_eps + cfs[2] * nu + cfs[3] * khi +
           cfs[4] * var_eps * nu + cfs[5] * var_eps * khi + cfs[6] * nu * khi +
           cfs[7] * var_eps * nu * khi - 8.0 * _y; }

inline static double_t coordinates_converter::dF2deps(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y) {
    return cfs[1] + cfs[4] * nu + cfs[5] * khi + cfs[7] * nu * khi; }

inline static double_t coordinates_converter::dF2deta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y) {
    return cfs[2] + cfs[4] * var_eps + cfs[6] * khi + cfs[7] * var_eps * khi; }

inline static double_t coordinates_converter::dF2dzeta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _y) {
    return cfs[3] + cfs[5] * var_eps + cfs[6] * nu + cfs[7] * var_eps * nu; }

inline static double_t coordinates_converter::F3(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z) {
    return cfs[0] + cfs[1] * var_eps + cfs[2] * nu + cfs[3] * khi +
           cfs[4] * var_eps * nu + cfs[5] * var_eps * khi + cfs[6] * nu * khi +
           cfs[7] * var_eps * nu * khi - 8.0 * _z; }

inline static double_t coordinates_converter::dF3deps(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z) {
    return cfs[1] + cfs[4] * nu + cfs[5] * khi + cfs[7] * nu * khi; }

inline static double_t coordinates_converter::dF3deta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z) {
    return cfs[2] + cfs[4] * var_eps + cfs[6] * khi + cfs[7] * var_eps * khi; }

inline static double_t coordinates_converter::dF3dzeta(double_t var_eps, double_t nu, double_t khi, cube cfs, double_t _z) {
    return cfs[3] + cfs[5] * var_eps + cfs[6] * nu + cfs[7] * var_eps * nu; }

matrix coordinates_converter::Jacobian::find_inverse(double_t var_eps, double_t nu, double_t khi,
    cube cfs_x, cube cfs_y, cube cfs_z) {
    auto a = dxde(var_eps, nu, khi, cfs_x);
    auto b = dxdn(var_eps, nu, khi, cfs_x);
    auto c = dxdc(var_eps, nu, khi, cfs_x);

    auto d = dyde(var_eps, nu, khi, cfs_y);
    auto e = dydn(var_eps, nu, khi, cfs_y);
    auto f = dydc(var_eps, nu, khi, cfs_y);

    auto g = dzde(var_eps, nu, khi, cfs_z);
    auto h = dzdn(var_eps, nu, khi, cfs_z);
    auto i = dzdc(var_eps, nu, khi, cfs_z);

    auto detJ = a * e * i + b * f * g + c * d * h -
        c * e * g - b * d * i - a * f * h;

    return matrix{ {{(e * i - h * f) / detJ, (h * c - b * i) / detJ, (b * f - c * e) / detJ},
                    {(g * f - d * i) / detJ, (a * i - c * g) / detJ, (c * d - a * f) / detJ},
                    {(d * h - g * e) / detJ, (b * g - a * h) / detJ, (a * e - b * d) / detJ}} };
}



void coordinates_converter::solve_gauss(matrix& A, vector& b, vector& solution) {
    const int n = 3;

    for (int k = 0; k < n; ++k) {
        int max_row = k;
        double max_val = std::abs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > max_val) {
                max_val = std::abs(A[i][k]);
                max_row = i;
            }
        }

        if (max_row != k) {
            std::swap(A[k], A[max_row]);
            std::swap(b[k], b[max_row]);
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j)
                A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        solution[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            solution[i] -= A[i][j] * solution[j];
        solution[i] /= A[i][i];
    }
}

vector coordinates_converter::convert_from_xyz(double x0, double y0, double z0,
					                           const cube& x, const cube& y, const cube& z,
					                           double_t _eps, double_t max_iter) {
	double eps(0.0), eta(0.0), zeta(0.0);

       cube cfs_x{ x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7], // * Constant
                  -x[0] + x[1] - x[2] + x[3] - x[4] + x[5] - x[6] + x[7], // * eps
                  -x[0] - x[1] + x[2] + x[3] - x[4] - x[5] + x[6] + x[7], // * eta
                  -x[0] - x[1] - x[2] - x[3] + x[4] + x[5] + x[6] + x[7], // * zeta
                   x[0] - x[1] - x[2] + x[3] + x[4] - x[5] - x[6] + x[7], // * eps eta
                   x[0] - x[1] + x[2] - x[3] - x[4] + x[5] - x[6] + x[7], // * eps zeta
                   x[0] + x[1] - x[2] - x[3] - x[4] - x[5] + x[6] + x[7], // * eta zeta
                  -x[0] + x[1] + x[2] - x[3] + x[4] - x[5] - x[6] + x[7]  // * eps eta zeta
       };

       cube cfs_y{ y[0] + y[1] + y[2] + y[3] + y[4] + y[5] + y[6] + y[7],
                  -y[0] + y[1] - y[2] + y[3] - y[4] + y[5] - y[6] + y[7],
                  -y[0] - y[1] + y[2] + y[3] - y[4] - y[5] + y[6] + y[7],
                  -y[0] - y[1] - y[2] - y[3] + y[4] + y[5] + y[6] + y[7],
                   y[0] - y[1] - y[2] + y[3] + y[4] - y[5] - y[6] + y[7],
                   y[0] - y[1] + y[2] - y[3] - y[4] + y[5] - y[6] + y[7],
                   y[0] + y[1] - y[2] - y[3] - y[4] - y[5] + y[6] + y[7],
                  -y[0] + y[1] + y[2] - y[3] + y[4] - y[5] - y[6] + y[7] };

       cube cfs_z{ z[0] + z[1] + z[2] + z[3] + z[4] + z[5] + z[6] + z[7],
                  -z[0] + z[1] - z[2] + z[3] - z[4] + z[5] - z[6] + z[7],
                  -z[0] - z[1] + z[2] + z[3] - z[4] - z[5] + z[6] + z[7],
                  -z[0] - z[1] - z[2] - z[3] + z[4] + z[5] + z[6] + z[7],
                   z[0] - z[1] - z[2] + z[3] + z[4] - z[5] - z[6] + z[7],
                   z[0] - z[1] + z[2] - z[3] - z[4] + z[5] - z[6] + z[7],
                   z[0] + z[1] - z[2] - z[3] - z[4] - z[5] + z[6] + z[7],
                  -z[0] + z[1] + z[2] - z[3] + z[4] - z[5] - z[6] + z[7] };

       size_t current_iter = 0; double_t current_eps = DBL_MAX;


       matrix A; 
       vector b;

       const double_t v0_1 = F1(eps, eta, zeta, cfs_x, x0);
       const double_t v0_2 = F2(eps, eta, zeta, cfs_y, y0);
       const double_t v0_3 = F3(eps, eta, zeta, cfs_z, z0);

       double_t vi_1 = v0_1;
       double_t vi_2 = v0_2;
       double_t vi_3 = v0_3;

       // Solving process.
       while (current_eps > _eps and current_iter < max_iter) {
           A[0][0] = dF1deps(eps, eta, zeta, cfs_x, x0); A[0][1] = dF1deta(eps, eta, zeta, cfs_x, x0); A[0][2] = dF1dzeta(eps, eta, zeta, cfs_x, x0); b[0] = -F1(eps, eta, zeta, cfs_x, x0);
           A[1][0] = dF2deps(eps, eta, zeta, cfs_y, y0); A[1][1] = dF2deta(eps, eta, zeta, cfs_y, y0); A[1][2] = dF2dzeta(eps, eta, zeta, cfs_y, y0); b[1] = -F2(eps, eta, zeta, cfs_y, y0);
           A[2][0] = dF3deps(eps, eta, zeta, cfs_z, z0); A[2][1] = dF3deta(eps, eta, zeta, cfs_z, z0); A[2][2] = dF3dzeta(eps, eta, zeta, cfs_z, z0); b[2] = -F3(eps, eta, zeta, cfs_z, z0);


           std::array<double_t, 3> delta_x;
           solve_gauss(A, b, delta_x);

           double_t v_1 = F1(eps + delta_x[0], eta + delta_x[1], zeta + delta_x[2], cfs_x, x0);
           double_t v_2 = F2(eps + delta_x[0], eta + delta_x[1], zeta + delta_x[2], cfs_y, y0);
           double_t v_3 = F3(eps + delta_x[0], eta + delta_x[1], zeta + delta_x[2], cfs_z, z0);
           eps += delta_x[0]; eta += delta_x[1]; zeta += delta_x[2];
           current_eps = sqrt(v_1 * v_1 + v_2 * v_2 + v_3 * v_3) / sqrt(v0_1 * v0_1 + v0_2 * v0_2 + v0_3 * v0_3);
           current_iter++;
       }
	return vector {eps, eta, zeta};
}