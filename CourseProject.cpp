#include "Mesh.h"
#include "Drawer.h"

int main() {
    Mesh myMesh;
    ReadData(myMesh, standartInputPath);
    if (!myMesh.CheckData()) return -1;
    myMesh.Generate();
    myMesh.FileWriteGeneratedPoints();
    myMesh.FileWriteGeneratedRibs();
    myMesh.FileWriteGeneratedAreas();
    Drawer::DrawMesh();
    return 0;
}