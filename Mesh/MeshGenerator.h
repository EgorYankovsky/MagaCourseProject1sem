#pragma once

#include "Mesh.h"

#include <cassert>
#include <vector>

class MeshGenerator {
private:
    static void GenerateListOfPoints(Mesh& mesh);
    static void GenerateListOfAreas(Mesh& mesh);
    static void GenerateListOfRibs(Mesh& mesh);
    static void GenerateListOfBorders(Mesh& mesh);
public:
    MeshGenerator() = delete;
    static void Generate3DMesh(Mesh& mesh);
    friend void Sort(std::vector<Point>& arr);
};