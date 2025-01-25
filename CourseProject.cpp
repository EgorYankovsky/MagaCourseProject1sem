#include "MainHeader.h"

#include <functional>

#include "Drawer\\Drawer.h"

static auto SelectTest() -> std::string;

int subtract(int a, int b) { return a - b; }

int main() {
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