#include "Mesh.h"
#include "Drawer.h"

int main() {
    Mesh myMesh;
    ReadData(myMesh, standartInputPath);
    auto isDataCorrect = myMesh.CheckData();
    if (!isDataCorrect) return -1;
    myMesh.Generate();
    myMesh.FileWriteGeneratedPoints();
    Drawer::DrawMesh();
    return 0;
}