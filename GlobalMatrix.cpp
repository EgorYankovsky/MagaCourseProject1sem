#include "GlobalMatrix.h"

void GlobalMatrix::GeneratePortrait(std::vector<std::array<size_t, 13>> areas, size_t ribsAmount) {
    _ig.resize(ribsAmount + 1);

    std::vector<std::vector<size_t>> additionalVector(ribsAmount);
    
    for (const auto& area : areas)
        for (size_t i(1); i < 13; ++i)
            for (size_t j(1); j < 13; ++j)
                if (area[i] < area[j] and std::find(additionalVector[area[i]].begin(), additionalVector[area[i]].end(), area[j]) == additionalVector[area[i]].end()) {
                    additionalVector[area[i]].push_back(area[j]);
                    std::sort(additionalVector[area[i]].begin(), additionalVector[area[i]].end());
                }

    _ig[0] = 0;
    for (size_t i(0); i < ribsAmount; ++i) {
        _ig[i + 1] = _ig[i] + additionalVector[i].size();
        _jg.insert(_jg.end(), additionalVector[i].begin(), additionalVector[i].end());
    }
}

void GlobalMatrix::Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points,
    std::vector<std::pair<size_t, size_t>> generatedRibs, std::vector<std::pair<size_t, std::pair<double, double>>> areasInfo) {
    for (const auto& area : areas) {

    }
}

double GlobalMatrix::operator()(size_t i, size_t j) const {
    Logger::ConsoleOutput("Can't get value from global matrix.", NotificationColor::Alert);
    exit(-1);
    return 0.0;
}

double& GlobalMatrix::operator()(size_t i, size_t j)
{
    Logger::ConsoleOutput("Can't set value for global matrix.", NotificationColor::Alert);
    exit(-1);
}

GlobalMatrix::GlobalMatrix() {
}

GlobalMatrix::~GlobalMatrix() {
}
