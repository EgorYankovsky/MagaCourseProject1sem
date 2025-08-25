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

const std::string standartInputPath = "Data\\Input\\input.txt";
const std::string testInputPath = "Data\\Input\\testInput.txt";
const std::string summerTestInputPath = "Data\\Input\\test_straight.txt";
const std::string inputPath1 = "Data\\Input\\input1.txt";
const std::string inputPath2 = "Data\\Input\\input2.txt";
const std::string inputPath3 = "Data\\Input\\input3.txt";
const std::string inputPath4 = "Data\\Input\\input4.txt";
const std::string inputPath5 = "Data\\Input\\input5.txt";

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