#include "MainHeader.h"

#include "Drawer\\Drawer.h"
#include "JacobiMatrix.h"

static auto SelectTest() -> std::string;

int main() {
    std::array<double, 8> x{ 0, 1, 0, 1, 0, 1, 0, 1 };
    std::array<double, 8> y{ 0, 0, 1, 1, 0, 0, 1, 1 };
    std::array<double, 8> z{ 0, 0, 0, 0, 1, 1, 1, 1 };
    JacobiMatrix::SetValues(x, y, z);

    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j)
            std::cout << std::scientific << std::setprecision(3) << (JacobiMatrix::GetValueAt(i, j))(0.85, 0.65, 0.5) << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (size_t i(0); i < 3; ++i) {
        for (size_t j(0); j < 3; ++j)
            std::cout << std::scientific << std::setprecision(3) << (JacobiMatrix::GetValueAtTransposed(i, j))(0.85, 0.65, 0.5) << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << (JacobiMatrix::GetDeterminant())(0.85, 0.65, 0.5) << std::endl;
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