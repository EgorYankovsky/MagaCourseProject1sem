#include "GlobalVector.h"

void GlobalVector::addLocalVectorValues(const std::array<size_t, 12> localRibs, const LocalVector& b) {
    const std::array<size_t, 12> switchV{
        0, 3, 8, 11,
        1, 2, 9, 10,
        4, 5, 6, 7 };
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

void GlobalVector::CommitBoundaryConditions(std::vector<std::array<size_t, 3>> borderRibs) {
    for (const auto& rib : borderRibs) {
        switch (rib[0]) {
        case 2:
        case 3:
            Logger::ConsoleOutput("Can't commit boundary conditions of 2nd or 3rd type.", NotificationColor::Alert);
            exit(-1);
            break;
        case 1:

        default:
            break;
        }
    }
}
