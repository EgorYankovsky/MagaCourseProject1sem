#pragma once

#include "Mesh.h"

class MeshGenerator {
private:
    static void GenerateListOfPointsAboveX(Mesh& mesh);
    static void GenerateListOfPointsAboveY(Mesh& mesh);
    static void GenerateListOfPointsAboveZ(Mesh& mesh);
    static void GenerateListOfAreas(Mesh& mesh);
    static void GenerateListOfRibs(Mesh& mesh);
    static void GenerateListOfBorders(Mesh& mesh);
public:
    MeshGenerator() = delete;
    static void Generate3DMesh(Mesh& mesh);
};