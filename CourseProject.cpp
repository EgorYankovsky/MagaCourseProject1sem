#include "MainHeader.h"

#include "Drawer\\Drawer.h"

#include "LocalMatrix.h"
#include "LocalVector.h"

static auto SelectTest() -> std::string;

int main() {
    
    //std::array<double, 8> x{ -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5 };
    //std::array<double, 8> y{ -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5 };
    //std::array<double, 8> z{ -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5 };

    std::array<double, 8> x{ -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, 1.0 };
    std::array<double, 8> y{ -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0 };
    std::array<double, 8> z{ -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0, 1.0 };

    
    //LocalVector lv(x, y, z);

    //JacobiMatrix::SetValues(x, y, z);

    // For matrix of stiffness.

    /*
    for (size_t i(0); i < 12; ++i) {
        auto f = BasisFunction::getRotAt(i);
        std::array < std::function<double(double, double, double)>, 3> _v1{
            J::GetValueAtTransposed(0, 0) * f[0] + J::GetValueAtTransposed(0, 1) * f[1] + J::GetValueAtTransposed(0, 2) * f[2],
            J::GetValueAtTransposed(1, 0) * f[0] + J::GetValueAtTransposed(1, 1) * f[1] + J::GetValueAtTransposed(1, 2) * f[2],
            J::GetValueAtTransposed(2, 0) * f[0] + J::GetValueAtTransposed(2, 1) * f[1] + J::GetValueAtTransposed(2, 2) * f[2],
        };

        for (size_t ii(0); ii < 3; ++ii) {
            std::cout << std::scientific << std::setprecision(6) << _v1[ii](0.5, 0.5, 0.5) << std::endl;
        }

        std::cout << std::endl;
    }
    */

    // For matrix of mass.
    // 
    //std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(BasisFunction::getAt(0)) << std::endl;
    /*
    std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(JacobiMatrix::dydn) << std::endl;
    std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(JacobiMatrix::dzdc) << std::endl;
    std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(JacobiMatrix::dydc) << std::endl;
    std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(JacobiMatrix::dzdn) << std::endl;
    std::cout << std::scientific << std::setprecision(6) << JacobiMatrix::GetDeterminant()(0, 0, 1) << std::endl;


    std::cout << std::scientific << std::setprecision(6) << JacobiMatrix::GetValueAtInverse(0, 0)(0, 0, 1) << std::endl;
    
    std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(JacobiMatrix::GetValueAtInverse(0, 0)) << std::endl;
    std::cout << std::scientific << std::setprecision(6) << Integration::Gauss5(JacobiMatrix::GetValueAtInverse(0, 0) * BasisFunction::getAt(0)) << std::endl;
    */
    //std::cout << Integration::Gauss5(JacobiMatrix::GetValueAtInverse(k, i / 4) * BasisFunction::getAt(i) *
    //                                 JacobiMatrix::GetValueAtInverse(k, j / 4) * BasisFunction::getAt(j) *
    //                                 JacobiMatrix::GetDeterminant());

    //return 0;
    

    //std::array<double, 8> x{ 0, 1, 0, 1, 0, 1, 0, 1 };
    //std::array<double, 8> y{ 0, 0, 1, 1, 0, 0, 1, 1 };
    //std::array<double, 8> z{ 0, 0, 0, 0, 1, 1, 1, 1 };

    
    JacobiMatrix::SetValues(x, y, z);
    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j) {
            std::cout << JacobiMatrix::GetValueAt(i, j)(0.0, 0.0, 0.0) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j) {
            std::cout << JacobiMatrix::GetValueAtTransposed(i, j)(0.0, 0.0, 0.0) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j) {
            std::cout << JacobiMatrix::GetValueAtInverse(i, j)(0.0, 0.0, 0.0) << " ";
        }
        std::cout << std::endl;
    }
    
    //return 0;
    

    LocalMatrix LM(1.0, x, y, z, LMType::Stiffness);

    for (size_t i(0); i < 12; ++i) {
        for (size_t j(0); j < 12; ++j) {
            std::cout << std::scientific << std::setprecision(6) << LM(i, j) << " ";
            if (j % 4 == 3) std::cout << "\t";
        }
        std::cout << std::endl;
        if (i % 4 == 3) std::cout << std::endl;
    }

    return 0;
    auto inputPath = SelectTest();
    Mesh myMesh;
    MeshFileStreamer::Read(myMesh, inputPath);
    if (!myMesh.CheckData()) {
        Logger::ConsoleOutput("Error during data checking.", NotificationColor::Alert);
        return -1;
    }
    MeshGenerator::Generate3DMesh(myMesh);
    MeshFileStreamer::Write(&myMesh, FileExtension::Txt);

    Drawer::DrawMesh(PictureOutput::SaveAsFile);

    FEM myFEM3D;
    myFEM3D.GetMeshData(&myMesh);
    myFEM3D.Type = EquationType::Elliptical;
    myFEM3D.StartSolution();
    


    return 0;
}

static auto SelectTest() -> std::string {
    std::cout << "Select test num:" << std::endl;
    std::cout << "(0) Standard cubic mesh." << std::endl;
    std::cout << "(1) Emerald mesh." << std::endl;
    std::cout << "(2) Beveled pyramid mesh." << std::endl;
    std::cout << "(3) Hourglass-shaped mesh." << std::endl;
    std::cout << "(4) Bath mesh." << std::endl;
    std::cout << "(5) Detailed emerald mesh." << std::endl;
    std::cout << "(6) Random figure mesh." << std::endl;
    std::cout << "-> ";
    size_t input(0);
    std::cin >> input;
    switch (input)
    {
    case 0: return standartInputPath;
    case 1: return inputPath1;
    case 2: return inputPath2;
    case 3: return inputPath3;
    case 4: return inputPath4;
    case 5: return inputPath5;
    case 6: return testInputPath;
    default:
        system("cls");
        return SelectTest();
    }
}