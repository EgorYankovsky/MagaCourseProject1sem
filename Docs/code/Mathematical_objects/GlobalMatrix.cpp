#include "GlobalMatrix.h"

void GlobalMatrix::addLocalMatrixValues(const std::array<size_t, 12> localRibs, const LocalMatrix& G, const LocalMatrix& M) {
    const std::array<size_t, 12> switchV{
        0, 3, 8, 11,
        1, 2, 9, 10,
        4, 5, 6, 7 };
    int ii(0);
    for (const auto& i : localRibs) {
        int jj(0);
        for (const auto& j : localRibs) {
            int ind(0);
            if (i - j == 0)
                _di[i] += G(switchV[ii], switchV[jj]) + M(switchV[ii], switchV[jj]);
            else if (i - j < 0) {
                ind = _ig[j];
                for (; ind <= _ig[j + 1] - 1; ind++) if (_jg[ind] == i) break;
                _au[ind] += G(switchV[ii], switchV[jj]) + M(switchV[ii], switchV[jj]);
            }
            else if (i - j > 0) {
                ind = _ig[i];
                for (; ind <= _ig[i + 1] - 1; ind++) if (_jg[ind] == j) break;
                _al[ind] += G(switchV[ii], switchV[jj]) + M(switchV[ii], switchV[jj]);
            } ++jj; } ++ii; }
}

double GlobalMatrix::getAlValue(size_t i, size_t j) const {
    for (size_t ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
        if (_jg[_ig[i] + ii] == j) return _al[_ig[i] + ii];    
    return 0.0;
}

double GlobalMatrix::getAuValue(size_t i, size_t j) const {
    for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
        if (_jg[_ig[i] + ii] == j) return _au[_ig[i] + ii];
    return 0.0;
}

void GlobalMatrix::GeneratePortrait(std::vector<std::array<size_t, 13>> areas, size_t ribsAmount) {
    _ig.resize(ribsAmount + 1);
    std::vector<std::vector<size_t>> additionalVector(ribsAmount);
    for (const auto& area : areas)
        for (size_t i(1); i < 13; ++i)
            for (size_t j(1); j < 13; ++j)
                if (area[j] < area[i] and std::find(additionalVector[area[i]].begin(), additionalVector[area[i]].end(), area[j]) == additionalVector[area[i]].end()) {
                    additionalVector[area[i]].push_back(area[j]);
                    std::sort(additionalVector[area[i]].begin(), additionalVector[area[i]].end()); }
    _ig[0] = 0;
    for (size_t i(0); i < ribsAmount; i++) {
        _ig[i + 1] = _ig[i] + additionalVector[i].size();
        _jg.insert(_jg.end(), additionalVector[i].begin(), additionalVector[i].end()); }
    _di.resize(ribsAmount); _al.resize(_jg.size()); _au.resize(_jg.size());
}

void GlobalMatrix::Fill(std::vector<std::array<size_t, 13>> areas, std::vector<std::array<double, 3>> points,
    std::vector<std::pair<size_t, size_t>> generatedRibs, std::vector<std::pair<size_t, std::pair<double, double>>> areasInfo) {
    for (const auto& area : areas) {
        auto mu = [area, areasInfo](){
            for (const auto& info : areasInfo)
                if (area[0] == info.first) return info.second.first; };
        auto sigma = [area, areasInfo]() {
            for (const auto& info : areasInfo)
                if (area[0] == info.first) return info.second.second; };
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
        LocalMatrix localG(mu(), xPoints, yPoints, zPoints, LMType::Stiffness);
        LocalMatrix localM(sigma(), xPoints, yPoints, zPoints, LMType::Mass);
        addLocalMatrixValues(localArea, localG, localM); }
}

void GlobalMatrix::CommitBoundaryConditions(std::vector<std::array<size_t, 6>> borderRibs) {
    for (const auto& rib : borderRibs) {
        switch (rib[0]) {
        case 2: case 3:
            Logger::ConsoleOutput("Can't commit boundary conditions of 2nd or 3rd type.", NotificationColor::Alert);
            exit(-1); break;
        case 1:
            for (size_t i(2); i < 6; ++i) {
                for (size_t j(_ig[rib[i]]); j < _ig[rib[i] + 1]; ++j) _al[j] = 0;
                _di[rib[i]] = 1.0;
                for (size_t j = 0; j < _jg.size(); ++j) if (_jg[j] == rib[i]) _au[j] = 0;
            } break;
        default: break; } }
}

double GlobalMatrix::getValue(size_t i, size_t j) {
    if (i == j) return _di[i];
    else if (i - j < 0) getAuValue(i, j);
    else if (i - j > 0) getAlValue(j, i);
}

double GlobalMatrix::operator()(size_t i, size_t j) const {
    if (i == j) return _di[i];
    else if (i - j < 0) getAuValue(i, j);
    else if (i - j > 0) getAlValue(i, j);
}

double& GlobalMatrix::operator()(size_t i, size_t j) {
    Logger::ConsoleOutput("Can't set value for global matrix.", NotificationColor::Alert);
    exit(-1); }

GlobalVector operator*(const GlobalMatrix A, const GlobalVector b) {
    if (A.Size != b.Size) Logger::ConsoleOutput("Matrix and vector have different sizes during multiplication", NotificationColor::Alert);
    GlobalVector ans(b.Size);
    for (size_t i(0); i < b.Size; ++i) {
        for (size_t j(0); j < A._ig[i + 1] - A._ig[i]; ++j) {
            ans(i) += A._al[A._ig[i] + j] * b(A._jg[A._ig[i] + j]);
            ans(A._jg[A._ig[i] + j]) += A._au[A._ig[i] + j] * b(i); }
        ans(i) += A._di[i] * b(i); }
    return ans;
}