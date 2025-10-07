#include "FEM.h"

void FEM::SolveElliptical() {
    A = new GlobalMatrix();
    b = new GlobalVector(_generatedRibs.size());
    x = GlobalVector(_generatedRibs.size());

    A->GeneratePortrait(_areas, _generatedRibs.size());
    A->Fill(_areas, _points, _generatedRibs, _areasInfo); // To slow.
    b->Fill(_areas, _points, _generatedRibs);             // To slow.
    A->CommitBoundaryConditions(_newBorderRibs);
    b->CommitBoundaryConditions(_areas, _points, _generatedRibs);
    //b->CommitBoundaryConditions(_newBorderRibs, _points, _generatedRibs);
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

void FEM::ConsoleTestOutput() {

    auto switcher = [](int& num) {
        if (num == 0) num = 0;
        else if (num == 1) num = 3;
        else if (num == 2) num = 8;
        else if (num == 3) num = 11;
        else if (num == 4) num = 1;
        else if (num == 5) num = 2;
        else if (num == 6) num = 9;
        else if (num == 7) num = 10;
        else if (num == 8) num = 4;
        else if (num == 9) num = 5;
        else if (num == 10) num = 6;
        else if (num == 11) num = 7;
        else throw std::exception("Conversation error.");
        };

    auto mv_multiplication = [](const matrix& m, const vector& v) -> vector {
        return vector{ m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
                       m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
                       m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] };
        };

    for (const auto& area : _areas) {
        double center_x(0.0), center_y(0.0), center_z(0.0);

        
        std::array<double, 8> x_array{ _points[_generatedRibs[area[5]].first][0], _points[_generatedRibs[area[6]].first][0],
                                       _points[_generatedRibs[area[7]].first][0], _points[_generatedRibs[area[8]].first][0],
                                       _points[_generatedRibs[area[5]].second][0], _points[_generatedRibs[area[6]].second][0],
                                       _points[_generatedRibs[area[7]].second][0], _points[_generatedRibs[area[8]].second][0], };

        std::array<double, 8> y_array{ _points[_generatedRibs[area[5]].first][1], _points[_generatedRibs[area[6]].first][1],
                                       _points[_generatedRibs[area[7]].first][1], _points[_generatedRibs[area[8]].first][1],
                                       _points[_generatedRibs[area[5]].second][1], _points[_generatedRibs[area[6]].second][1],
                                       _points[_generatedRibs[area[7]].second][1], _points[_generatedRibs[area[8]].second][1], };

        std::array<double, 8> z_array{ _points[_generatedRibs[area[5]].first][2], _points[_generatedRibs[area[6]].first][2],
                                       _points[_generatedRibs[area[7]].first][2], _points[_generatedRibs[area[8]].first][2],
                                       _points[_generatedRibs[area[5]].second][2], _points[_generatedRibs[area[6]].second][2],
                                       _points[_generatedRibs[area[7]].second][2], _points[_generatedRibs[area[8]].second][2], };

        for (const auto& val : x_array) center_x += val; center_x /= 8.0;
        for (const auto& val : y_array) center_y += val; center_y /= 8.0;
        for (const auto& val : z_array) center_z += val; center_z /= 8.0;

        decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(center_x, center_y, center_z,
                                                                                   x_array, y_array, z_array);
        decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
            x_array, y_array, z_array);

        vector ans{ 0.0, 0.0, 0.0 };
        for (int i = 0; i < 12; ++i) {
            int index = i; switcher(index);
            decltype(auto) phi = vector{ BasisFunction::get_at(i)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                         BasisFunction::get_at(i)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                         BasisFunction::get_at(i)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2]) };
            decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);
            decltype(auto) q_i = x(area[index + 1]);
            ans[0] += q_i * psiVector[0]; ans[1] += q_i * psiVector[1]; ans[2] += q_i * psiVector[2];
        }
        std::cout << std::scientific << std::setprecision(6) << center_x << " " << center_y << " " << center_z << " | ";
        std::cout << std::scientific << std::setprecision(6) << ans[0] << " " << ans[1] << " " << ans[2] << std::endl;
        // Account geometrical center of current area.
        // Find x, y, z.
        //auto value = GetSolutionAtPoint(center_x, center_y, center_z);
    }
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
