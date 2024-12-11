#include "Mesh.h"
#include "Drawer.h"

int main() {
    Mesh myMesh;
    ReadData(myMesh, standartInputPath);
    myMesh.Generate();
    myMesh.FileWriteGeneratedPoints();
    Drawer::DrawMesh();
    return 0;
}