#pragma once

#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <iomanip>

#include "Mesh/MeshHeader.h"
#include "FEM/FEM.h"
#include "Solvers/MainSolverHeader.h"
#include "Drawer/Drawer.h"

#pragma region badMeshes
const std::string emeraldMesh = "Data\\Input\\emeraldMesh.txt";
const std::string skewedPyramidMesh = "Data\\Input\\skewedPyramidMesh.txt";
const std::string bathMesh = "Data\\Input\\bathMesh.txt";
const std::string detailedEmeraldMesh = "Data\\Input\\detailedEmeraldMesh.txt";
const std::string randomFigureMesh = "Data\\Input\\randomFigureMesh.txt";
const std::string summerTestInputPath = "Data\\Input\\summerTestInputPath.txt";
#pragma endregion

const std::string standardCubicMesh = "Data\\Input\\StandardCubicMesh.txt";
const std::string PointedMesh = "Data\\Input\\PointedMesh.txt";
const std::string DiagonalMesh = "Data\\Input\\DiagonalMesh.txt";

void CheckMass(const std::array<double, 8>& x, 
               const std::array<double, 8>& y, 
               const std::array<double, 8>& z,
               const LocalMatrix& m);

void CheckStiffness(const std::array<double, 8>& x, 
                    const std::array<double, 8>& y, 
                    const std::array<double, 8>& z,
                    const LocalMatrix& m);

void TestLocalMatrixAndVector();

void TestNewVector();

//const std::string defaultOutputPointsPath = "Data\\Generated\\generatedPoints.txt";
//const std::string defaultOutputRibsPath = "Data\\Generated\\generatedRibs.txt";
//const std::string defaultOutputAreasPath = "Data\\Generated\\generatedAreas.txt";