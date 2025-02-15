#include "GlobalVector.h"

void GlobalVector::addLocalVectorValues(const std::array<size_t, 12> localRibs, const LocalVector& b) {
    const std::array<size_t, 12> switchV{
        0, 3, 8, 11,
        1, 2, 9, 10,
        4, 5, 6, 7 };
    //const std::array<size_t, 12> switchV{
    //    0, 4, 5, 1,
    //    8, 9, 10, 11,
    //    2, 6, 7, 3 };
    for (size_t i(0); i < b.getSize(); ++i)
        _values[localRibs[i]] += b(switchV[i]);
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
        std::array<size_t, 12> localArea{ area[1], area[4], area[9], area[12],
                                  area[2], area[3], area[10], area[11],
                                  area[5], area[6], area[7], area[8] };

        std::array<double, 8> xPoints = { points[generatedRibs[localArea[0]].first][0], points[generatedRibs[localArea[0]].second][0],
                                          points[generatedRibs[localArea[1]].first][0], points[generatedRibs[localArea[1]].second][0],

                                          points[generatedRibs[localArea[2]].first][0], points[generatedRibs[localArea[2]].second][0],
                                          points[generatedRibs[localArea[3]].first][0], points[generatedRibs[localArea[3]].second][0] };

        std::array<double, 8> yPoints = { points[generatedRibs[localArea[0]].first][1], points[generatedRibs[localArea[0]].second][1],
                                          points[generatedRibs[localArea[1]].first][1], points[generatedRibs[localArea[1]].second][1],

                                          points[generatedRibs[localArea[2]].first][1], points[generatedRibs[localArea[2]].second][1],
                                          points[generatedRibs[localArea[3]].first][1], points[generatedRibs[localArea[3]].second][1] };

        std::array<double, 8> zPoints = { points[generatedRibs[localArea[0]].first][2], points[generatedRibs[localArea[0]].second][2],
                                          points[generatedRibs[localArea[1]].first][2], points[generatedRibs[localArea[1]].second][2],

                                          points[generatedRibs[localArea[2]].first][2], points[generatedRibs[localArea[2]].second][2],
                                          points[generatedRibs[localArea[3]].first][2], points[generatedRibs[localArea[3]].second][2] };

        LocalVector b(xPoints, yPoints, zPoints);

    }
}

void GlobalVector::CommitBoundaryConditions(std::vector<std::array<size_t, 6>> borderRibs, std::vector<std::array<double, 3>> points, std::vector<std::pair<size_t, size_t>> generatedRibs) {
    for (const auto& square : borderRibs) {
        size_t r0 = square[2];
        size_t r1 = square[3];
        std::array<double, 3> _x = { points[generatedRibs[r0].first][0], points[generatedRibs[r0].second][0], points[generatedRibs[r1].second][0] };
        std::array<double, 3> _y = { points[generatedRibs[r0].first][1], points[generatedRibs[r0].second][1], points[generatedRibs[r1].second][1] };
        std::array<double, 3> _z = { points[generatedRibs[r0].first][2], points[generatedRibs[r0].second][2], points[generatedRibs[r1].second][2] };

        auto getNormal = [_x, _y, _z]() -> vector {
            auto v = vector{ (_y[1] - _y[0]) * (_z[2] - _z[0]) - (_z[1] - _z[0]) * (_y[2] - _y[0]),
                     -1.0 * ((_x[1] - _x[0]) * (_z[2] - _z[0]) - (_z[1] - _z[0]) * (_x[2] - _x[0])),
                             (_x[1] - _x[0]) * (_y[2] - _y[0]) - (_y[1] - _y[0]) * (_x[2] - _x[0]) };
            double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return vector{ v[0] / len, v[1] / len, v[2] / len };
            };

        auto normal = getNormal();
        switch (square[0]) {
        case 2:
        case 3:
            Logger::ConsoleOutput("Can't commit boundary conditions of 2nd or 3rd type.", NotificationColor::Alert);
            exit(-1);
            break;
        case 1:
            for (size_t ii(2); ii < 6; ++ii) {
                std::array<double, 3> middlePoint{ 0.5 * (points[generatedRibs[square[ii]].first][0] + points[generatedRibs[square[ii]].second][0]),
                                                   0.5 * (points[generatedRibs[square[ii]].first][1] + points[generatedRibs[square[ii]].second][1]), 
                                                   0.5 * (points[generatedRibs[square[ii]].first][2] + points[generatedRibs[square[ii]].second][2]), };
                auto fVector = Function::TestA(middlePoint[0], middlePoint[1], middlePoint[2], 0.0);
                auto fValue = fVector[0] * normal[0] + fVector[1] * normal[1] + fVector[2] * normal[2];
                _values[square[ii]] = fValue;
            }
            break;
        default:
            break;
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
