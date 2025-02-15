#include "MainHeader.h"

static auto SelectTest() -> std::string;

int main() {

    /*
    std::array<double, 8> x{ -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, 1.0 };
    std::array<double, 8> y{ -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0 };
    std::array<double, 8> z{ -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0, 1.0 };
    //std::array<double, 8> x{ 0.0, 2.0,  1.0, 2.5,  0.0, 1.0, -1.0, 4.0 };
    //std::array<double, 8> y{ 0.0, 0.0,  2.0, 2.0, -1.0, 0.0,  3.0, 3.0 };
    //std::array<double, 8> z{ 0.0, 1.0, -1.0, 0.0,  1.0, 2.0,  1.0, 3.0 };

    
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
    
    return 0;
    

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
    */
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
    myFEM3D.BuildMatrixAndVector();
    myFEM3D.SetSolver(new LOS());
    myFEM3D.Solve();
    myFEM3D.WriteAnswer();
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