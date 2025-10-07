#include "MainHeader.h"

static auto SelectTest() -> std::string;

/* Mesh checking
def check_correctness(_x, _y, _z, _indexes):
    a = _x[_indexes[1]] - _x[_indexes[0]]
    b = _y[_indexes[1]] - _y[_indexes[0]]
    c = _z[_indexes[1]] - _z[_indexes[0]]

    d = _x[_indexes[2]] - _x[_indexes[0]]
    e = _y[_indexes[2]] - _y[_indexes[0]]
    f = _z[_indexes[2]] - _z[_indexes[0]]

    def fun(x, y, z):
        return (b * f - c * e) * (x - _x[_indexes[0]]) - (a * f - c * d) * (y - _y[_indexes[0]]) + (a * e - b * d) * (z - _z[_indexes[0]])
    return fun(_x[_indexes[3]], _y[_indexes[3]], _z[_indexes[3]]) == 0
*/

int main() {

    //TestLocalMatrixAndVector();
    //TestNewVector();
    //return 0;
    
    auto inputPath = SelectTest();
    Mesh myMesh;

    MeshFileStreamer::Read(myMesh, inputPath);
    if (!myMesh.CheckData()) {
        Logger::ConsoleOutput("Error during data checking.", NotificationColor::Alert);
        return -1;
    }
    
    MeshGenerator::Generate3DMesh(myMesh);
    MeshFileStreamer::Write(&myMesh, FileExtension::Txt);

    Drawer::DrawMesh(PictureOutput::ShowOnDesktop);
    //return 0;
    FEM myFEM3D;
    
    myFEM3D.GetMeshData(&myMesh);
    myFEM3D.Type = EquationType::Elliptical;
    myFEM3D.BuildMatrixAndVector();
    myFEM3D.SetSolver(new LOS());
    myFEM3D.Solve();
    myFEM3D.WriteAnswer();
    myFEM3D.ConsoleTestOutput();
    return 0;
}

static auto SelectTest() -> std::string {
    std::cout << "Select test num:" << std::endl;
    std::cout << "(0) Standard cubic mesh." << std::endl;
    std::cout << "(1) Pointed mesh." << std::endl;
    std::cout << "(2) Diagonal mesh." << std::endl;
    std::cout << "-> ";
    size_t input(0);
    std::cin >> input;
    switch (input)
    {
    case 0: return standardCubicMesh;
    case 1: return PointedMesh;
    case 2: return DiagonalMesh;
    default:
        system("cls");
        return SelectTest();
    }
}

void CheckMass(const std::array<double, 8>& x, 
               const std::array<double, 8>& y, 
               const std::array<double, 8>& z,
               const LocalMatrix& m) {
    double hx = x[1] - x[0];
    double hy = y[2] - y[0];
    double hz = z[4] - z[0];

    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            if (i / 4 == j / 4)  {
                auto fi = BasisFunction::getAt(i);
                auto fj = BasisFunction::getAt(j);
                auto integrated = Integration::Gauss3(fi * fj);
                if (i / 4 == 0)
                    integrated *= ((hy * hz) / (2.0 * hx));
                else if (i / 4 == 1)
                    integrated *= ((hx * hz) / (2.0 * hy));
                else if (i / 4 == 2)
                    integrated *= ((hy * hx) / (2.0 * hz));
                std::cout << std::scientific << std::setprecision(6) << integrated << " ";
            }
            else
                std::cout << std::scientific << std::setprecision(6) << 0.0 << " ";
            if (j % 4 == 3) std::cout << "\t";
        }
        std::cout << std::endl;
        if (i % 4 == 3) std::cout << std::endl;
    }
}

void CheckStiffness(const std::array<double, 8>& x,
                    const std::array<double, 8>& y,
                    const std::array<double, 8>& z,
                    const LocalMatrix& m) {
    double hx = x[1] - x[0];
    double hy = y[2] - y[0];
    double hz = z[4] - z[0];

    std::array<std::array<double, 12>, 12 > test{};

    for (int i = 0; i < 12; ++i)
    {
        for (int j  = 0; j < 12; ++j)
        {
            auto rotPhi_i = BasisFunction::getRotAt(i);
            auto rotPhi_j = BasisFunction::getRotAt(j);            
            double ans = 0.0;
            ans += (2.0 * hx) / (hy * hz) * Integration::Gauss3(rotPhi_i[0] * rotPhi_j[0]);
            ans += (2.0 * hy) / (hx * hz) * Integration::Gauss3(rotPhi_i[1] * rotPhi_j[1]);
            ans += (2.0 * hz) / (hx * hy) * Integration::Gauss3(rotPhi_i[2] * rotPhi_j[2]);
            test[i][j] = ans;
            
            std::cout << std::scientific << std::setprecision(6) << ans << " ";
            if (j % 4 == 3) std::cout << "\t";
        }
        std::cout << std::endl;
        if (i % 4 == 3) std::cout << std::endl;
    }

    for (int i = 0; i < 12 ; ++i) {
        for (int j = 0; j < 12; ++j) {
            if (std::abs(test[i][j] - m(i, j)) > 1e-6) {
                std::cout << "Error at (" << i << ", " << j << "): expected " << test[i][j] << ", got " << m(i, j) << std::endl;
            }
        }
    }
}

void TestLocalMatrixAndVector() {
    std::array<double, 8> x{ -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, 1.0 };
    std::array<double, 8> y{ -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0 };
    std::array<double, 8> z{ -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0, 1.0 };

    //std::array<double, 8> x{ 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
    //std::array<double, 8> y{ 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };

    //std::array<double, 8> x{ 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
    //std::array<double, 8> y{ 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 4.0, 4.0 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0 };

    //std::array<double, 8> x{ 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5 };
    //std::array<double, 8> y{ 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5 };

    //std::array<double, 8> x{ 0.0, 8.0, 0.0, 8.0, 0.0, 8.0, 0.0, 8.0 };
    //std::array<double, 8> y{ 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 8.0, 8.0 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0 };


    //std::array<double, 8> x{ 0.0, 2.0,  1.0, 2.5,  0.0, 1.0, -1.0, 4.0 };
    //std::array<double, 8> y{ 0.0, 0.0,  2.0, 2.0, -1.0, 0.0,  3.0, 3.0 };
    //std::array<double, 8> z{ 0.0, 1.0, -1.0, 0.0,  1.0, 2.0,  1.0, 3.0 };


    JacobiMatrix::SetValues(x, y, z);
    std::cout << "|J| = " << JacobiMatrix::GetDeterminant()(0.0, 0.0, 0.0) << std::endl;
    std::cout << "|J| = " << JacobiMatrix::GetDeterminant()(0.5, 0.5, 0.5) << std::endl;
    std::cout << "|J| = " << JacobiMatrix::GetDeterminant()(1.0, 1.0, 1.0) << std::endl;
    std::cout << "S |J| de dn dc = " << Integration::Gauss3(JacobiMatrix::GetDeterminant()) << std::endl;
    //return 0;
    std::cout << std::endl << "J => " << std::endl;
    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j) {
            std::cout << JacobiMatrix::GetValueAt(i, j)(0.0, 0.0, 0.0) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "J^(-1) => " << std::endl;
    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j) {
            std::cout << JacobiMatrix::GetValueAtInverse(i, j)(0.0, 0.0, 0.0) << " ";
        }
        std::cout << std::endl;
    }

    //return 0;


    LocalMatrix LM(1.0, x, y, z, LMType::Stiffness);

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    CheckStiffness(x, y, z, LM);

    //LocalVector lv(x, y, z);
    //for (int i = 0; i < 12; ++i) {
    //    if (i % 4 == 0) std::cout << std::endl;
    //    std::cout << std::scientific << std::setprecision(6) << lv(i) << std::endl;
    //}
}

void TestNewVector() {
    //std::array<double, 8> x{ -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, 1.0 };
    //std::array<double, 8> y{ -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0 };
    //std::array<double, 8> z{ -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0, 1.0 };

    std::array<double, 8> x{ 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
    std::array<double, 8> y{ 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
    std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };

    //std::array<double, 8> x{ 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
    //std::array<double, 8> y{ 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 4.0, 4.0 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0 };

    //std::array<double, 8> x{ 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5 };
    //std::array<double, 8> y{ 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5 };

    //std::array<double, 8> x{ 0.0, 8.0, 0.0, 8.0, 0.0, 8.0, 0.0, 8.0 };
    //std::array<double, 8> y{ 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 8.0, 8.0 };
    //std::array<double, 8> z{ 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0 };


    //std::array<double, 8> x{ 0.0, 2.0,  1.0, 2.5,  0.0, 1.0, -1.0, 4.0 };
    //std::array<double, 8> y{ 0.0, 0.0,  2.0, 2.0, -1.0, 0.0,  3.0, 3.0 };
    //std::array<double, 8> z{ 0.0, 1.0, -1.0, 0.0,  1.0, 2.0,  1.0, 3.0 };


    LocalVector lv1(x, y, z);
    //lv1.generateNew();
    for (int i = 0; i < 12; ++i) {
        if (i % 4 == 0) std::cout << std::endl;
        std::cout << std::scientific << std::setprecision(6) << lv1(i) << std::endl;
    }

}
