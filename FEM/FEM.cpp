#include "FEM.h"

void FEM::SolveElliptical() {
    A = new GlobalMatrix();
    b = new GlobalVector(_generatedRibs.size());
    x = GlobalVector(_generatedRibs.size());

    A->GeneratePortrait(_areas, _generatedRibs.size());
    A->Fill(_areas, _points, _generatedRibs, _areasInfo); // To slow.
    b->Fill(_areas, _points, _generatedRibs);             // To slow.
    A->CommitBoundaryConditions(_newBorderRibs);
    b->CommitBoundaryConditions(_newBorderRibs, _points, _generatedRibs);
}

void FEM::SolveParabolical() {
    Logger::ConsoleOutput("Couldn't read Parabolic problem", NotificationColor::Alert);
    exit(-1);
}

void FEM::SolveHyperbolical() {
    Logger::ConsoleOutput("Couldn't read Hyperbolic problem", NotificationColor::Alert);
    exit(-1);
}

void FEM::ReadMeshData(InputExtension ie) {
    Logger::ConsoleOutput("Couldn't read from file", NotificationColor::Alert);
    exit(-1);
}

void FEM::GetMeshData(const Mesh* mesh) {
    for (const auto& point : mesh->Points) {
        std::array<double, 3> _point = { point.x, point.y, point.z };
        _points.push_back(_point);
    }

    for (const auto& area : mesh->getAreasAsRibs()) {
        std::array<size_t, 13> _area = { area.subdomainNum_,
            area.refs_[0], area.refs_[1], area.refs_[2], area.refs_[3],
            area.refs_[4], area.refs_[5], area.refs_[6], area.refs_[7],
            area.refs_[8], area.refs_[9], area.refs_[10], area.refs_[11] };
        _areas.push_back(_area);
    }

    for (const auto& rib : mesh->getRibsRefs()) {
        std::pair<size_t, size_t> _rib{ rib.p1, rib.p2 };
        _generatedRibs.push_back(_rib);
    }

    //for (const auto& borderRib : mesh->getBorderRibs()) {
    //    std::array<size_t, 3> _borderRib{ borderRib.type_, borderRib.formulaNum_, borderRib.ribRef_ };
    //    _borderRibs.push_back(_borderRib);
    //}

    for (const auto& border : mesh->getNewBorderRibs()) {
        std::array<size_t, 6> rwe{border[0], border[1], border[2], border[3], border[4], border[5]};
        _newBorderRibs.push_back(rwe);
    }

    for (const auto& areasInfos : mesh->AreasInfo) {
        std::pair<size_t, std::pair<double, double>> _area{ areasInfos.subdomainNum_, {areasInfos.mu_, areasInfos.sigma_} };
        _areasInfo.push_back(_area);
    }

    _isDataCommited = true;
}

void FEM::BuildMatrixAndVector() {
    switch (Type) {
    case EquationType::Hyperbolical:
        SolveHyperbolical();
        break;
    case EquationType::Parabolical:
        SolveParabolical();
        break;
    case EquationType::Elliptical:
        SolveElliptical();
        break;
    case EquationType::NotStated:
    default:
        Logger::ConsoleOutput("Equation type didn't stated. Exit program.", NotificationColor::Alert);
        exit(-1);
        break;
    }
}

vector FEM::GetSolutionAtPoint(double x, double y, double z)
{
    Logger::ConsoleOutput("Can't get solution at point.", NotificationColor::Alert);
    exit(-1);
    return vector();
}

void FEM::SetSolver(Solver* s) {
    _s = s;
}

void FEM::Solve() {
    x = _s->Solve(*A, *b);
}

void FEM::WriteAnswer() {
    std::ofstream fout("Data/Output/solution.txt");
    for (size_t i(0); i < x.getSize(); ++i) fout << i << ". " << std::setprecision(15) << std::scientific << x(i) << std::endl;
    fout.close();
}

FEM::FEM() {
    Logger::ConsoleOutput("FEM declared, but it's empty", NotificationColor::Warning);
}

FEM::~FEM() {}
