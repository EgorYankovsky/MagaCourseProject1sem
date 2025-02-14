#pragma once

#include <vector>
#include <array>

#include "../Mathematical_objects/MathematicalHeader.h"
#include "../Logger/Logger.h"
#include "../Mesh/Mesh.h"
#include "../Solvers/MainSolverHeader.h"

typedef std::array<double, 3> vector;

enum class EquationType {
    Hyperbolical,
    Parabolical,
    Elliptical,
    NotStated
};

enum class InputExtension {
    Txt,
    Bin
};

class FEM {

    bool _isDataCommited = false;
    bool isDataCommited() const { return _isDataCommited; }

    bool _isSolved = false;
    bool isSolved() const { return _isSolved; }

    // Based.
    std::vector<std::array<size_t, 13>> _areas{};
    std::vector<std::array<size_t, 3>> _borderRibs{};
    std::vector<std::array<size_t, 6>> _newBorderRibs{};
    std::vector<std::array<double, 3>> _points{};
    std::vector<std::pair<size_t, size_t>> _generatedRibs{};
    std::vector<std::pair<size_t, std::pair<double, double>>> _areasInfo{};
    std::vector<double> _time{};

    void SolveElliptical();
    void SolveParabolical();
    void SolveHyperbolical();

    GlobalMatrix* A;
    GlobalVector* b;
    GlobalVector* x;

    Solver* _s;

public:

    void ReadMeshData(InputExtension ie);
    void GetMeshData(const Mesh* mesh);
    void BuildMatrixAndVector();
    vector GetSolutionAtPoint(double x, double y, double z);
    void SetSolver(Solver* s);

    EquationType Type = EquationType::NotStated;

    __declspec(property(get = isDataCommited)) bool IsDataCommited;
    __declspec(property(get = isSolved)) bool IsSolved;

    FEM();
    FEM(std::vector<std::array<size_t, 13>> areas,
        std::vector<std::array<size_t, 3>> borderRibs,
        std::vector<std::array<double, 3>> points,
        std::vector<std::pair<size_t, size_t>> generatedRibs) :
        _areas(areas), _borderRibs(borderRibs),
        _points(points), _generatedRibs(generatedRibs) {}
    ~FEM();
};

