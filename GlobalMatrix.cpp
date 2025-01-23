#include "GlobalMatrix.h"

void GlobalMatrix::GeneratePortrait(std::vector<std::array<size_t, 13>> areas, size_t ribsAmount) {
    _ig.resize(ribsAmount + 1);

    std::vector<std::vector<size_t>> additionalVector(ribsAmount);
    
    for (const auto& area : areas) {
        for (size_t i(1); i < 13; ++i) {
            for (size_t j(1); j < 13; ++j) {
                if (area[i] < area[j] and std::find(additionalVector[area[i]].begin(), additionalVector[area[i]].end(), area[j]) == additionalVector[area[i]].end()) {
                    additionalVector[area[i]].push_back(area[j]);
                    std::sort(additionalVector[area[i]].begin(), additionalVector[area[i]].end());
                }
            }
        }
    }

    _ig[0] = 0;
    for (size_t i(0); i < ribsAmount; ++i) {
        _ig[i + 1] = _ig[i] + additionalVector[i].size();
        _jg.insert(_jg.end(), additionalVector[i].begin(), additionalVector[i].end());
    }

    /*
    List<List<int>> arr = [];

    // ! Дерьмодристный момент.
    for (int i = 0; i < arrPtLen; i++)
        arr.Add(new List<int>());

    foreach(var _elem in arrEl)
        foreach(var point in _elem)
        foreach(var pnt in _elem)
        if (pnt < point && Array.IndexOf(arr[point].ToArray(), pnt) == -1)
        {
            arr[point].Add(pnt);
            arr[point].Sort();
        }

    m._ig[0] = 0;
    for (int i = 0; i < arrPtLen; i++)
    {
        m._ig[i + 1] = m._ig[i] + arr[i].Count;
        m._jg.AddRange(arr[i]);
    }
    */

}

void GlobalMatrix::Fill() {
    Logger::ConsoleOutput("Can't fill matrix", NotificationColor::Alert);
    exit(-1);
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
