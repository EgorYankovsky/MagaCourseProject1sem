#include "MainHeader.h"

#include "Drawer\\Drawer.h"
#include "Integration\\Integration.h"

static auto SelectTest() -> std::string;

int main() {
    auto inputPath = SelectTest();
    Mesh myMesh;
    MeshFileStreamer::Read(myMesh, inputPath);
    if (!myMesh.CheckData()) return -1;
    MeshGenerator::Generate3DMesh(myMesh);
    MeshFileStreamer::Write(&myMesh, FileExtension::Txt);
    return 0;
    myMesh.Generate();
    Drawer::DrawMesh();
    return 0;
}

static auto SelectTest() -> std::string {
    std::cout << "Select test num:" << std::endl;
    std::cout << "(0) Standard cubic mesh." << std::endl;
    std::cout << "(1) Emerald mesh." << std::endl;
    std::cout << "(2) Beveled pyramid mesh." << std::endl;
    std::cout << "(3) Hourglass-shaped mesh." << std::endl;
    std::cout << "(4) Car mesh." << std::endl;
    std::cout << "(5) C*ck mesh." << std::endl;
    std::cout << "(6) Test mesh." << std::endl;
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