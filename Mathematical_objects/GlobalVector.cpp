#include "GlobalVector.h"

void GlobalVector::addLocalVectorValues(const std::array<size_t, 12> localRibs, const LocalVector& b) {
    //const std::array<size_t, 12> switchV{
    //    0, 3, 8, 11,
    //    1, 2, 9, 10,
    //    4, 5, 6, 7 };
    for (size_t i(0); i < b.getSize(); ++i)
        _values[localRibs[i]] += b(i);
}

double GlobalVector::operator()(size_t i) const {
    if (i >= Size) {
        Logger::ConsoleOutput("Index run out of Vector range.", NotificationColor::Alert);
        exit(-1);
    }
    return i < Size ? _values[i] : throw "Index run out of Vector range.";
}

double& GlobalVector::operator()(size_t i)
{
    if (i >= Size) {
        Logger::ConsoleOutput("Index run out of Vector range.", NotificationColor::Alert);
        exit(-1);
    }
    return _values[i];
}

GlobalVector::GlobalVector() {
    Logger::ConsoleOutput("Global vector initialized, but it's empty", NotificationColor::Warning);
}

GlobalVector::GlobalVector(size_t size) {
    _values.resize(size);
}

void GlobalVector::Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points, 
                        std::vector<std::pair<size_t, size_t>> generatedRibs) {
    for (const auto& area : areas) {
        std::array<size_t, 12> localArea{ area[1], area[2], area[3], area[4],
                                          area[5], area[6], area[7], area[8],
                                          area[9], area[10], area[11], area[12] };

        std::array<size_t, 12> swiftArea{
            localArea[0], localArea[3], localArea[8], localArea[11],
            localArea[1], localArea[2], localArea[9], localArea[10],
            localArea[4], localArea[5], localArea[6], localArea[7]
        };

        std::array<double, 8> xPoints = { points[generatedRibs[area[1]].first][0], points[generatedRibs[area[1]].second][0],
                                          points[generatedRibs[area[4]].first][0], points[generatedRibs[area[4]].second][0],

                                          points[generatedRibs[area[9]].first][0], points[generatedRibs[area[9]].second][0],
                                          points[generatedRibs[area[12]].first][0], points[generatedRibs[area[12]].second][0] };

        std::array<double, 8> yPoints = { points[generatedRibs[area[1]].first][1], points[generatedRibs[area[1]].second][1],
                                          points[generatedRibs[area[4]].first][1], points[generatedRibs[area[4]].second][1],

                                          points[generatedRibs[area[9]].first][1], points[generatedRibs[area[9]].second][1],
                                          points[generatedRibs[area[12]].first][1], points[generatedRibs[area[12]].second][1] };

        std::array<double, 8> zPoints = { points[generatedRibs[area[1]].first][2], points[generatedRibs[area[1]].second][2],
                                          points[generatedRibs[area[4]].first][2], points[generatedRibs[area[4]].second][2],

                                          points[generatedRibs[area[9]].first][2], points[generatedRibs[area[9]].second][2],
                                          points[generatedRibs[area[12]].first][2], points[generatedRibs[area[12]].second][2] };

        LocalVector b(xPoints, yPoints, zPoints);
        addLocalVectorValues(swiftArea, b);
    }
}

void GlobalVector::CommitBoundaryConditions(std::vector<std::array<size_t, 6>> border_ribs, 
                                            std::vector<std::array<double, 3>> points, 
                                            std::vector<std::pair<size_t, size_t>> generated_ribs) {
    for (const auto& square : border_ribs) {
        //size_t r0 = square[2];
        //size_t r1 = square[3];
        //std::array<double, 3> _x = { points[generatedRibs[r0].first][0], points[generatedRibs[r0].second][0], points[generatedRibs[r1].second][0] };
        //std::array<double, 3> _y = { points[generatedRibs[r0].first][1], points[generatedRibs[r0].second][1], points[generatedRibs[r1].second][1] };
        //std::array<double, 3> _z = { points[generatedRibs[r0].first][2], points[generatedRibs[r0].second][2], points[generatedRibs[r1].second][2] };

        //auto getNormal = [_x, _y, _z]() -> vector {
        //    auto v = vector{ (_y[1] - _y[0]) * (_z[2] - _z[0]) - (_z[1] - _z[0]) * (_y[2] - _y[0]),
        //             -1.0 * ((_x[1] - _x[0]) * (_z[2] - _z[0]) - (_z[1] - _z[0]) * (_x[2] - _x[0])),
        //                     (_x[1] - _x[0]) * (_y[2] - _y[0]) - (_y[1] - _y[0]) * (_x[2] - _x[0]) };
        //    double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        //    return vector{ v[0] / len, v[1] / len, v[2] / len };
        //    };

        //auto normal = getNormal();
        switch (square[0]) {
        case 2:
        case 3:
            Logger::ConsoleOutput("Can't commit boundary conditions of 2nd or 3rd type.", NotificationColor::Alert);
            exit(-1);
            break;
        case 1:
            for (size_t ii(2); ii < 6; ++ii) {
                
                auto get_normal = [points](std::pair<size_t, size_t> vector_points) -> vector {
                    double x0 = points[vector_points.first][0];
                    double y0 = points[vector_points.first][1];
                    double z0 = points[vector_points.first][2];
                    double x1 = points[vector_points.second][0];
                    double y1 = points[vector_points.second][1];
                    double z1 = points[vector_points.second][2];
                    auto v = vector{ x1 - x0, y1 - y0, z1 - z0 };
                    double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                    return vector{ v[0] / len, v[1] / len, v[2] / len };
                };

                auto normal = get_normal(generated_ribs[square[ii]]);
                
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[square[ii]].first][0] + points[generated_ribs[square[ii]].second][0]),
                                                    0.5 * (points[generated_ribs[square[ii]].first][1] + points[generated_ribs[square[ii]].second][1]), 
                                                    0.5 * (points[generated_ribs[square[ii]].first][2] + points[generated_ribs[square[ii]].second][2]), };

                auto fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                auto fValue = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];
                
                
                
                
                _values[square[ii]] = fValue;
            }
            break;
        default:
            break;
        }
    }
}

void GlobalVector::CommitBoundaryConditions(const std::vector<std::array<size_t, 13>>& areas,
    const std::vector<std::array<double, 3>>& points,
    const std::vector<std::pair<size_t, size_t>>& generated_ribs) {
    // nx, ny, nz - amount of points along axes.
    
    std::set<size_t> committed_ribs{};

    int nx = areas[0][2] + 1;
    int ny = (areas[0][5] + nx) / (2 * nx - 1);
    int nz = (areas[areas.size() - 1][12] + 1 + nx * ny) / (3 * nx * ny - nx - ny);


    auto get_normal = [points](const std::pair<size_t, size_t>& vector_points) -> vector {
        double x0 = points[vector_points.first][0]; double y0 = points[vector_points.first][1];
        double z0 = points[vector_points.first][2]; double x1 = points[vector_points.second][0];
        double y1 = points[vector_points.second][1]; double z1 = points[vector_points.second][2];
        auto v = vector{ x1 - x0, y1 - y0, z1 - z0 };
        double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        return vector{ v[0] / len, v[1] / len, v[2] / len }; };

    
    auto mv_multiplication = [](const matrix & m, const vector & v) -> vector {
        return vector{ m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
                       m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
                       m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] };
    };

    auto switcher = [](ptrdiff_t& num) {
        if (num == 0) num = num;
        else if (num == 1) num = 4;
        else if (num == 2) num = 5;
        else if (num == 3) num = 1;
        else if (num == 4) num = 8;
        else if (num == 5) num = 9;
        else if (num == 6) num = 10;
        else if (num == 7) num = 11;
        else if (num == 8) num = 2;
        else if (num == 9) num = 6;
        else if (num == 10) num = 7;
        else if (num == 11) num = 3;
        else throw std::exception("Conversation error.");
        };

    // oXY0
    for (int i(0); i < ny - 1; ++i) {
        for (int j(0); j < nx - 1; ++j) {
            decltype(auto) selected_area = areas[i * (nx - 1) + j]; // Dirty moment. Int and size_t multiplication.

            std::array<double, 8> x_array{ points[generated_ribs[selected_area[5]].first][0], points[generated_ribs[selected_area[6]].first][0],
                                           points[generated_ribs[selected_area[7]].first][0], points[generated_ribs[selected_area[8]].first][0],
                                           points[generated_ribs[selected_area[5]].second][0], points[generated_ribs[selected_area[6]].second][0],
                                           points[generated_ribs[selected_area[7]].second][0], points[generated_ribs[selected_area[8]].second][0], };

            std::array<double, 8> y_array{ points[generated_ribs[selected_area[5]].first][1], points[generated_ribs[selected_area[6]].first][1],
                                           points[generated_ribs[selected_area[7]].first][1], points[generated_ribs[selected_area[8]].first][1],
                                           points[generated_ribs[selected_area[5]].second][1], points[generated_ribs[selected_area[6]].second][1],
                                           points[generated_ribs[selected_area[7]].second][1], points[generated_ribs[selected_area[8]].second][1], };

            std::array<double, 8> z_array{ points[generated_ribs[selected_area[5]].first][2], points[generated_ribs[selected_area[6]].first][2],
                                           points[generated_ribs[selected_area[7]].first][2], points[generated_ribs[selected_area[8]].first][2],
                                           points[generated_ribs[selected_area[5]].second][2], points[generated_ribs[selected_area[6]].second][2],
                                           points[generated_ribs[selected_area[7]].second][2], points[generated_ribs[selected_area[8]].second][2], };

            std::array<size_t, 4> boundary_flat{ selected_area[1], selected_area[2], selected_area[3], selected_area[4] };
            for (const auto& rib : boundary_flat) {
                if (committed_ribs.find(rib) != committed_ribs.end()) continue;
                committed_ribs.insert(rib);
                auto normal = get_normal(generated_ribs[rib]);
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[rib].first][0] + points[generated_ribs[rib].second][0]),
                                                    0.5 * (points[generated_ribs[rib].first][1] + points[generated_ribs[rib].second][1]), 
                                                    0.5 * (points[generated_ribs[rib].first][2] + points[generated_ribs[rib].second][2]), };
                
                decltype(auto) fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                decltype(auto) numerator = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];

                decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(middle_point[0], middle_point[1], middle_point[2], 
                                                                                           x_array, y_array, z_array); // Account basis function psi at the current rib's middle point.
                
                decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
                                                                                                x_array, y_array, z_array);

                decltype(auto) current_local_rib_index = std::find(selected_area.begin() + 1, selected_area.end(), rib) - (selected_area.begin() + 1);
                
                switcher(current_local_rib_index);
                
                decltype(auto) phi = vector{ BasisFunction::get_at(current_local_rib_index)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]), 
                                             BasisFunction::get_at(current_local_rib_index)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2])};
                
                decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);

                decltype(auto) denominator = psiVector[0] * normal[0] + psiVector[1] * normal[1] + psiVector[2] * normal[2];

                _values[rib] = numerator / denominator;
            }
        }
    }


    // oXY1
    for (int i(0); i < ny - 1; ++i) {
        for (int j(0); j < nx - 1; ++j) {
            decltype(auto) selected_area = areas[(nz - 2) * (nx - 1) * (ny - 1) + i * (nx - 1) + j]; // Dirty moment. Int and size_t multiplication.

            std::array<double, 8> x_array{ points[generated_ribs[selected_area[5]].first][0], points[generated_ribs[selected_area[6]].first][0],
                                           points[generated_ribs[selected_area[7]].first][0], points[generated_ribs[selected_area[8]].first][0],
                                           points[generated_ribs[selected_area[5]].second][0], points[generated_ribs[selected_area[6]].second][0],
                                           points[generated_ribs[selected_area[7]].second][0], points[generated_ribs[selected_area[8]].second][0], };

            std::array<double, 8> y_array{ points[generated_ribs[selected_area[5]].first][1], points[generated_ribs[selected_area[6]].first][1],
                                           points[generated_ribs[selected_area[7]].first][1], points[generated_ribs[selected_area[8]].first][1],
                                           points[generated_ribs[selected_area[5]].second][1], points[generated_ribs[selected_area[6]].second][1],
                                           points[generated_ribs[selected_area[7]].second][1], points[generated_ribs[selected_area[8]].second][1], };

            std::array<double, 8> z_array{ points[generated_ribs[selected_area[5]].first][2], points[generated_ribs[selected_area[6]].first][2],
                                           points[generated_ribs[selected_area[7]].first][2], points[generated_ribs[selected_area[8]].first][2],
                                           points[generated_ribs[selected_area[5]].second][2], points[generated_ribs[selected_area[6]].second][2],
                                           points[generated_ribs[selected_area[7]].second][2], points[generated_ribs[selected_area[8]].second][2], };

            std::array<size_t, 4> boundary_flat{ selected_area[9], selected_area[10], selected_area[11], selected_area[12] };
            for (const auto& rib : boundary_flat) {
                if (committed_ribs.find(rib) != committed_ribs.end()) continue;
                committed_ribs.insert(rib);
                auto normal = get_normal(generated_ribs[rib]);
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[rib].first][0] + points[generated_ribs[rib].second][0]),
                                                    0.5 * (points[generated_ribs[rib].first][1] + points[generated_ribs[rib].second][1]),
                                                    0.5 * (points[generated_ribs[rib].first][2] + points[generated_ribs[rib].second][2]), };

                decltype(auto) fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                decltype(auto) numerator = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];

                decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(middle_point[0], middle_point[1], middle_point[2],
                    x_array, y_array, z_array); // Account basis function psi at the current rib's middle point.

                decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
                    x_array, y_array, z_array);

                decltype(auto) current_local_rib_index = std::find(selected_area.begin() + 1, selected_area.end(), rib) - (selected_area.begin() + 1);

                switcher(current_local_rib_index);

                decltype(auto) phi = vector{ BasisFunction::get_at(current_local_rib_index)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2]) };

                decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);

                decltype(auto) denominator = psiVector[0] * normal[0] + psiVector[1] * normal[1] + psiVector[2] * normal[2];

                _values[rib] = numerator / denominator;
            }
        }
    }


    // oX0Z
    for (int i(0); i < nz - 1; ++i) {
        for (int j(0); j < nx - 1; ++j) {
            decltype(auto) selected_area = areas[(nx - 1) * (ny - 1) * i + (nx - 1) * j]; // Dirty moment. Int and size_t multiplication.

            std::array<double, 8> x_array{ points[generated_ribs[selected_area[5]].first][0], points[generated_ribs[selected_area[6]].first][0],
                                           points[generated_ribs[selected_area[7]].first][0], points[generated_ribs[selected_area[8]].first][0],
                                           points[generated_ribs[selected_area[5]].second][0], points[generated_ribs[selected_area[6]].second][0],
                                           points[generated_ribs[selected_area[7]].second][0], points[generated_ribs[selected_area[8]].second][0], };

            std::array<double, 8> y_array{ points[generated_ribs[selected_area[5]].first][1], points[generated_ribs[selected_area[6]].first][1],
                                           points[generated_ribs[selected_area[7]].first][1], points[generated_ribs[selected_area[8]].first][1],
                                           points[generated_ribs[selected_area[5]].second][1], points[generated_ribs[selected_area[6]].second][1],
                                           points[generated_ribs[selected_area[7]].second][1], points[generated_ribs[selected_area[8]].second][1], };

            std::array<double, 8> z_array{ points[generated_ribs[selected_area[5]].first][2], points[generated_ribs[selected_area[6]].first][2],
                                           points[generated_ribs[selected_area[7]].first][2], points[generated_ribs[selected_area[8]].first][2],
                                           points[generated_ribs[selected_area[5]].second][2], points[generated_ribs[selected_area[6]].second][2],
                                           points[generated_ribs[selected_area[7]].second][2], points[generated_ribs[selected_area[8]].second][2], };

            std::array<size_t, 4> boundary_flat{ selected_area[2], selected_area[5], selected_area[7], selected_area[10] };
            for (const auto& rib : boundary_flat) {
                if (committed_ribs.find(rib) != committed_ribs.end()) continue;
                committed_ribs.insert(rib);
                auto normal = get_normal(generated_ribs[rib]);
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[rib].first][0] + points[generated_ribs[rib].second][0]),
                                                    0.5 * (points[generated_ribs[rib].first][1] + points[generated_ribs[rib].second][1]),
                                                    0.5 * (points[generated_ribs[rib].first][2] + points[generated_ribs[rib].second][2]), };

                decltype(auto) fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                decltype(auto) numerator = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];

                decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(middle_point[0], middle_point[1], middle_point[2],
                    x_array, y_array, z_array); // Account basis function psi at the current rib's middle point.

                decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
                    x_array, y_array, z_array);

                decltype(auto) current_local_rib_index = std::find(selected_area.begin() + 1, selected_area.end(), rib) - (selected_area.begin() + 1);

                switcher(current_local_rib_index);

                decltype(auto) phi = vector{ BasisFunction::get_at(current_local_rib_index)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2]) };

                decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);

                decltype(auto) denominator = psiVector[0] * normal[0] + psiVector[1] * normal[1] + psiVector[2] * normal[2];

                _values[rib] = numerator / denominator;
            }
        }
    }


    // oX1Z
    for (int i(0); i < nz - 1; ++i) {
        for (int j(0); j < nx - 1; ++j) {
            decltype(auto) selected_area = areas[(nx - 1) * (ny - 1) * i + (nx - 1) * (j + 1) - 1]; // Dirty moment. Int and size_t multiplication.

            std::array<double, 8> x_array{ points[generated_ribs[selected_area[5]].first][0], points[generated_ribs[selected_area[6]].first][0],
                                           points[generated_ribs[selected_area[7]].first][0], points[generated_ribs[selected_area[8]].first][0],
                                           points[generated_ribs[selected_area[5]].second][0], points[generated_ribs[selected_area[6]].second][0],
                                           points[generated_ribs[selected_area[7]].second][0], points[generated_ribs[selected_area[8]].second][0], };

            std::array<double, 8> y_array{ points[generated_ribs[selected_area[5]].first][1], points[generated_ribs[selected_area[6]].first][1],
                                           points[generated_ribs[selected_area[7]].first][1], points[generated_ribs[selected_area[8]].first][1],
                                           points[generated_ribs[selected_area[5]].second][1], points[generated_ribs[selected_area[6]].second][1],
                                           points[generated_ribs[selected_area[7]].second][1], points[generated_ribs[selected_area[8]].second][1], };

            std::array<double, 8> z_array{ points[generated_ribs[selected_area[5]].first][2], points[generated_ribs[selected_area[6]].first][2],
                                           points[generated_ribs[selected_area[7]].first][2], points[generated_ribs[selected_area[8]].first][2],
                                           points[generated_ribs[selected_area[5]].second][2], points[generated_ribs[selected_area[6]].second][2],
                                           points[generated_ribs[selected_area[7]].second][2], points[generated_ribs[selected_area[8]].second][2], };

            std::array<size_t, 4> boundary_flat{ selected_area[3], selected_area[6], selected_area[8], selected_area[11] };
            for (const auto& rib : boundary_flat) {
                if (committed_ribs.find(rib) != committed_ribs.end()) continue;
                committed_ribs.insert(rib);
                auto normal = get_normal(generated_ribs[rib]);
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[rib].first][0] + points[generated_ribs[rib].second][0]),
                                                    0.5 * (points[generated_ribs[rib].first][1] + points[generated_ribs[rib].second][1]),
                                                    0.5 * (points[generated_ribs[rib].first][2] + points[generated_ribs[rib].second][2]), };

                decltype(auto) fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                decltype(auto) numerator = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];

                decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(middle_point[0], middle_point[1], middle_point[2],
                    x_array, y_array, z_array); // Account basis function psi at the current rib's middle point.

                decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
                    x_array, y_array, z_array);

                decltype(auto) current_local_rib_index = std::find(selected_area.begin() + 1, selected_area.end(), rib) - (selected_area.begin() + 1);

                switcher(current_local_rib_index);

                decltype(auto) phi = vector{ BasisFunction::get_at(current_local_rib_index)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2]) };

                decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);

                decltype(auto) denominator = psiVector[0] * normal[0] + psiVector[1] * normal[1] + psiVector[2] * normal[2];

                _values[rib] = numerator / denominator;
            }
        }
    }

    // o0YZ
    for (int i(0); i < nz - 1; ++i) {
        for (int j(0); j < nx - 1; ++j) {
            decltype(auto) selected_area = areas[i * (nx - 1) * (ny - 1) + j]; // Dirty moment. Int and size_t multiplication.

            std::array<double, 8> x_array{ points[generated_ribs[selected_area[5]].first][0], points[generated_ribs[selected_area[6]].first][0],
                                           points[generated_ribs[selected_area[7]].first][0], points[generated_ribs[selected_area[8]].first][0],
                                           points[generated_ribs[selected_area[5]].second][0], points[generated_ribs[selected_area[6]].second][0],
                                           points[generated_ribs[selected_area[7]].second][0], points[generated_ribs[selected_area[8]].second][0], };

            std::array<double, 8> y_array{ points[generated_ribs[selected_area[5]].first][1], points[generated_ribs[selected_area[6]].first][1],
                                           points[generated_ribs[selected_area[7]].first][1], points[generated_ribs[selected_area[8]].first][1],
                                           points[generated_ribs[selected_area[5]].second][1], points[generated_ribs[selected_area[6]].second][1],
                                           points[generated_ribs[selected_area[7]].second][1], points[generated_ribs[selected_area[8]].second][1], };

            std::array<double, 8> z_array{ points[generated_ribs[selected_area[5]].first][2], points[generated_ribs[selected_area[6]].first][2],
                                           points[generated_ribs[selected_area[7]].first][2], points[generated_ribs[selected_area[8]].first][2],
                                           points[generated_ribs[selected_area[5]].second][2], points[generated_ribs[selected_area[6]].second][2],
                                           points[generated_ribs[selected_area[7]].second][2], points[generated_ribs[selected_area[8]].second][2], };

            std::array<size_t, 4> boundary_flat{ selected_area[1], selected_area[5], selected_area[6], selected_area[9] };
            for (const auto& rib : boundary_flat) {
                if (committed_ribs.find(rib) != committed_ribs.end()) continue;
                committed_ribs.insert(rib);
                auto normal = get_normal(generated_ribs[rib]);
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[rib].first][0] + points[generated_ribs[rib].second][0]),
                                                    0.5 * (points[generated_ribs[rib].first][1] + points[generated_ribs[rib].second][1]),
                                                    0.5 * (points[generated_ribs[rib].first][2] + points[generated_ribs[rib].second][2]), };

                decltype(auto) fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                decltype(auto) numerator = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];

                decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(middle_point[0], middle_point[1], middle_point[2],
                    x_array, y_array, z_array); // Account basis function psi at the current rib's middle point.

                decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
                    x_array, y_array, z_array);

                decltype(auto) current_local_rib_index = std::find(selected_area.begin() + 1, selected_area.end(), rib) - (selected_area.begin() + 1);

                switcher(current_local_rib_index);

                decltype(auto) phi = vector{ BasisFunction::get_at(current_local_rib_index)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2]) };

                decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);

                decltype(auto) denominator = psiVector[0] * normal[0] + psiVector[1] * normal[1] + psiVector[2] * normal[2];

                _values[rib] = numerator / denominator;
            }
        }
    }

    // o1YZ
    for (int i(0); i < nz - 1; ++i) {
        for (int j(0); j < nx - 1; ++j) {
            decltype(auto) selected_area = areas[i * (nx - 1) * (ny - 1) + j + (ny - 2) * (nx - 1)]; // Dirty moment. Int and size_t multiplication.

            std::array<double, 8> x_array{ points[generated_ribs[selected_area[5]].first][0], points[generated_ribs[selected_area[6]].first][0],
                                           points[generated_ribs[selected_area[7]].first][0], points[generated_ribs[selected_area[8]].first][0],
                                           points[generated_ribs[selected_area[5]].second][0], points[generated_ribs[selected_area[6]].second][0],
                                           points[generated_ribs[selected_area[7]].second][0], points[generated_ribs[selected_area[8]].second][0], };

            std::array<double, 8> y_array{ points[generated_ribs[selected_area[5]].first][1], points[generated_ribs[selected_area[6]].first][1],
                                           points[generated_ribs[selected_area[7]].first][1], points[generated_ribs[selected_area[8]].first][1],
                                           points[generated_ribs[selected_area[5]].second][1], points[generated_ribs[selected_area[6]].second][1],
                                           points[generated_ribs[selected_area[7]].second][1], points[generated_ribs[selected_area[8]].second][1], };

            std::array<double, 8> z_array{ points[generated_ribs[selected_area[5]].first][2], points[generated_ribs[selected_area[6]].first][2],
                                           points[generated_ribs[selected_area[7]].first][2], points[generated_ribs[selected_area[8]].first][2],
                                           points[generated_ribs[selected_area[5]].second][2], points[generated_ribs[selected_area[6]].second][2],
                                           points[generated_ribs[selected_area[7]].second][2], points[generated_ribs[selected_area[8]].second][2], };

            std::array<size_t, 4> boundary_flat{ selected_area[4], selected_area[7], selected_area[8], selected_area[12] };
            for (const auto& rib : boundary_flat) {
                if (committed_ribs.find(rib) != committed_ribs.end()) continue;
                committed_ribs.insert(rib);
                auto normal = get_normal(generated_ribs[rib]);
                std::array<double, 3> middle_point{ 0.5 * (points[generated_ribs[rib].first][0] + points[generated_ribs[rib].second][0]),
                                                    0.5 * (points[generated_ribs[rib].first][1] + points[generated_ribs[rib].second][1]),
                                                    0.5 * (points[generated_ribs[rib].first][2] + points[generated_ribs[rib].second][2]), };

                decltype(auto) fVector = Function::TestA(middle_point[0], middle_point[1], middle_point[2], 0.0);
                decltype(auto) numerator = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];

                decltype(auto) local_coordinates = coordinates_converter::convert_from_xyz(middle_point[0], middle_point[1], middle_point[2],
                    x_array, y_array, z_array); // Account basis function psi at the current rib's middle point.

                decltype(auto) inverse_Jacobian = coordinates_converter::Jacobian::find_inverse(local_coordinates[0], local_coordinates[1], local_coordinates[2],
                    x_array, y_array, z_array);

                decltype(auto) current_local_rib_index = std::find(selected_area.begin() + 1, selected_area.end(), rib) - (selected_area.begin() + 1);

                switcher(current_local_rib_index);

                decltype(auto) phi = vector{ BasisFunction::get_at(current_local_rib_index)[0](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[1](local_coordinates[0], local_coordinates[1], local_coordinates[2]),
                                             BasisFunction::get_at(current_local_rib_index)[2](local_coordinates[0], local_coordinates[1], local_coordinates[2]) };

                decltype(auto) psiVector = mv_multiplication(inverse_Jacobian, phi);

                decltype(auto) denominator = psiVector[0] * normal[0] + psiVector[1] * normal[1] + psiVector[2] * normal[2];

                _values[rib] = numerator / denominator;
            }
        }
    }
}

double GlobalVector::Norma() const {
    double sum(0.0);
    for (const auto& value : _values) sum += value * value;
    return sqrt(sum);
}

double operator*(const GlobalVector v1, const GlobalVector v2) {
    if (v1.Size != v2.Size) Logger::ConsoleOutput("During vector multiplication vectors have different size", NotificationColor::Alert);
    double sum(0.0);
    for (size_t i(0); i < v1.Size; ++i) sum += v1(i) * v2(i);
    return sum;
}

GlobalVector operator*(const double a, const GlobalVector v) {
    GlobalVector result(v.Size);
    for (size_t i(0); i < v.Size; ++i) result(i) = a * v(i);
    return result;
}

GlobalVector operator+(const GlobalVector v1, const GlobalVector v2) {
    if (v1.Size != v2.Size) Logger::ConsoleOutput("During vector multiplication vectors have different size", NotificationColor::Alert);
    GlobalVector result(v1.Size);
    for (size_t i(0); i < v1.Size; ++i) result(i) = v1(i) + v2(i);
    return result;
}

GlobalVector operator-(const GlobalVector v1, const GlobalVector v2) {
    if (v1.Size != v2.Size) Logger::ConsoleOutput("During vector multiplication vectors have different size", NotificationColor::Alert);
    GlobalVector result(v1.Size);
    for (size_t i(0); i < v1.Size; ++i) result(i) = v1(i) - v2(i);
    return result;
}
