#include "LocalMatrix.h"

typedef Integration I;

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

std::function<double(double, double, double)> operator/(double v, std::function<double(double, double, double)> f) {
	return [v, f](double t0, double t1, double t2) { return v / f(t0, t1, t2); };
}

void LocalMatrix::generate() {
	switch (_matrixType) {
	case LMType::Stiffness:
		generateG();
		break;
	case LMType::Mass:
		generateM1();
		break;
	case LMType::NotStated:
	default:
		break;
	}
}

void LocalMatrix::generateG() {
	J::SetValues(_x, _y, _z);
	for (size_t i(0); i < _localMatrixSize; ++i) {
		for (size_t j(0); j < _localMatrixSize; ++j) {

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
			
			//for (size_t ii(0); ii < 3; ++ii) {
			//	for (size_t jj(0); jj < 3; ++jj) {
			//		v1[ii] = v1[ii] + J::GetValueAtTransposed(ii, jj) * rotPhi_i[jj];
			//		v2[ii] = v2[ii] + J::GetValueAtTransposed(ii, jj) * rotPhi_j[jj];
			//	}
			//}

			for (size_t k(0); k < 3; ++k) _values[i][j] += Integration::Gauss5((v1[k] * v2[k]) / J::GetDeterminant());
			_values[i][j] /= _koef;
		}
	}
}

void LocalMatrix::generateM() {
	J::SetValues(_x, _y, _z);
	for (size_t i(0); i < _localMatrixSize; ++i) {
		for (size_t j(0); j < _localMatrixSize; ++j) {
			std::function<double(double, double, double)> f = [](double, double, double) { return 0.0; };
			for (size_t k(0); k < 3; ++k)
				f = f + J::GetValueAtInverse(k, i / 4) * J::GetValueAtInverse(k, j / 4);
			_values[i][j] = _koef * Integration::Gauss3(BasisFunction::getAt(i) * BasisFunction::getAt(j) * f);
			//_values[i][j] = _koef * Integration::Gauss3(J::GetDeterminant() * BasisFunction::getAt(i) * 
			//																  BasisFunction::getAt(j) * f);
		}
	}
}

void LocalMatrix::generateM1() {
	J::SetValues(_x, _y, _z);
	for (size_t i(0); i < _localMatrixSize; ++i) {
		for (size_t j(0); j < _localMatrixSize; ++j) {
//			_values[i][j] = _koef * I::Gauss3((BasisFunction::getAt(i) * BasisFunction::getAt(j) *
//											   (J::GetValueAtInverseNoDet(0, i / 4) * J::GetValueAtInverseNoDet(0, j / 4) +
//											    J::GetValueAtInverseNoDet(1, i / 4) * J::GetValueAtInverseNoDet(1, j / 4) +
//											    J::GetValueAtInverseNoDet(2, i / 4) * J::GetValueAtInverseNoDet(2, j / 4))) / J::GetDeterminant());
			_values[i][j] = _koef * I::Gauss3((BasisFunction::getAt(i) * BasisFunction::getAt(j) *
											   (J::GetValueAtInverseNoDet(0, i / 4) * J::GetValueAtInverseNoDet(0, j / 4) +
											    J::GetValueAtInverseNoDet(1, i / 4) * J::GetValueAtInverseNoDet(1, j / 4) +
											    J::GetValueAtInverseNoDet(2, i / 4) * J::GetValueAtInverseNoDet(2, j / 4))) / J::GetDeterminant());
		}
	}
}