#include "LocalMatrix.h"

std::function<double(double, double, double)> operator* (std::function<double(double, double, double)> f1,
														 std::function<double(double, double, double)> f2) {
	return [f1, f2](double t0, double t1, double t2) { return f1(t0, t1, t2) * f2(t0, t1, t2); };
}

std::function<double(double, double, double)> operator+(std::function<double(double, double, double)> f1,
														std::function<double(double, double, double)> f2) {
	return [f1, f2](double t0, double t1, double t2) { return f1(t0, t1, t2) + f2(t0, t1, t2); };
}

std::function<double(double, double, double)> operator/(std::function<double(double, double, double)> f1, 
														std::function<double(double, double, double)> f2) {
	return [f1, f2](double t0, double t1, double t2) { return f1(t0, t1, t2) / f2(t0, t1, t2); };
}

void LocalMatrix::generate() {
	switch (_matrixType) {
	case LMType::Stiffness:
		generateG();
		break;
	case LMType::Mass:
		generateM();
		break;
	case LMType::NotStated:
		break;
	default:
		break;
	}
}

void LocalMatrix::generateG() {
	for (size_t i(0); i < _localMatrixSize; ++i) {
		for (size_t j(0); j < _localMatrixSize; ++j) {

			J::SetValues(_x, _y, _z);

			auto rotPhi_i = BasisFunction::getRotAt(i);
			auto rotPhi_j = BasisFunction::getRotAt(j);

			std::array < std::function<double(double, double, double)>, 3> v1{
				J::GetValueAtTransposed(0, 0) * rotPhi_i[0] + J::GetValueAtTransposed(0, 1) * rotPhi_i[1] + J::GetValueAtTransposed(0, 2) * rotPhi_i[2],
				J::GetValueAtTransposed(1, 0) * rotPhi_i[0] + J::GetValueAtTransposed(1, 1) * rotPhi_i[1] + J::GetValueAtTransposed(1, 2) * rotPhi_i[2],
				J::GetValueAtTransposed(2, 0) * rotPhi_i[0] + J::GetValueAtTransposed(2, 1) * rotPhi_i[1] + J::GetValueAtTransposed(2, 2) * rotPhi_i[2],
			};
			
			std::array < std::function<double(double, double, double)>, 3> v2{
				J::GetValueAtTransposed(0, 0) * rotPhi_j[0] + J::GetValueAtTransposed(0, 1) * rotPhi_j[1] + J::GetValueAtTransposed(0, 2) * rotPhi_j[2],
				J::GetValueAtTransposed(1, 0) * rotPhi_j[0] + J::GetValueAtTransposed(1, 1) * rotPhi_j[1] + J::GetValueAtTransposed(1, 2) * rotPhi_j[2],
				J::GetValueAtTransposed(2, 0) * rotPhi_j[0] + J::GetValueAtTransposed(2, 1) * rotPhi_j[1] + J::GetValueAtTransposed(2, 2) * rotPhi_j[2],
			};
			
			for (size_t ii(0); ii < 3; ++ii) {
				for (size_t jj(0); jj < 3; ++jj) {
					v1[ii] = v1[ii] + J::GetValueAtTransposed(ii, jj) * rotPhi_i[jj];
					v2[ii] = v2[ii] + J::GetValueAtTransposed(ii, jj) * rotPhi_j[jj];
				}
			}

			for (size_t k(0); k < 3; ++k) _values[i][j] += Integration::Gauss5((v1[k] * v2[k]) / J::GetDeterminant());
			_values[i][j] /= _koef;
		}
	}
}

void LocalMatrix::generateM() {
	for (size_t i(0); i < _localMatrixSize; ++i) {
		for (size_t j(0); j < _localMatrixSize; ++j) {
			for (size_t k(0); k < 3; ++k) {
				J::SetValues(_x, _y, _z);
				_values[i][j] += Integration::Gauss5(J::GetValueAtInverse(k, i / 4) * BasisFunction::getAt(i) *
													 J::GetValueAtInverse(k, j / 4) * BasisFunction::getAt(j) *
													 J::GetDeterminant());
			}
			_values[i][j] *= _koef;
		}
	}
}
