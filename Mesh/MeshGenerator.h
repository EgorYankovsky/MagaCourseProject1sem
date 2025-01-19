#pragma once

#include "Mesh.h"

#include <cassert>
#include <vector>

typedef std::vector<std::vector<std::vector<Point>>> area3D;
typedef std::vector<std::vector<Point>> square2D;
typedef std::vector<Point> line1D;

class MeshGenerator {
private:
    static int SelectAreaNum(Mesh& mesh, std::array<size_t, 12> arr);
    static void GenerateListOfPoints(Mesh& mesh);
    static void GenerateListOfAreas(Mesh& mesh);
    static void GenerateListOfRibs(Mesh& mesh);
    static void GenerateListOfBorders(Mesh& mesh);
public:
    MeshGenerator() = delete;
    static void Generate3DMesh(Mesh& mesh);
    friend void Sort(std::vector<Point>& arr);
};