#pragma once

#include "..\\DataTypes.h"
#include "..\\Logger\\Logger.h"

#include <cassert>
#include <cmath>
#include <string>
#include <algorithm>

class Mesh {
private:
    // Based.
    bool isGenerated_ = false;
    bool isDeclarated_ = false;
    
    
    size_t linesAmountX_ = 0;
    size_t linesAmountY_ = 0;
    size_t linesAmountZ_ = 0;
    
    std::vector<Point> points_{};

    size_t subdomainsAmount_ = 0;
    std::vector<std::array<size_t, 7>> subdomains_{};
    std::vector<AreaInfo> areasInfo_{};
    std::vector<AreaRibs> areasRibs_{};

    std::vector<std::pair<size_t, double_t>> delimitersX_{};
    std::vector<std::pair<size_t, double_t>> delimitersY_{};
    std::vector<std::pair<size_t, double_t>> delimitersZ_{};

    size_t bordersAmount_ = 0;
    std::vector<Border> borders_{};

    std::vector<BorderLine> borderRibs_{};

    std::vector<std::array<size_t, 6>> newBorders_{};

    std::vector<AreaPoints> areasPoints_{};
    std::vector<RibRef> referableRibs_{};

    // Additional.    
    std::vector<Point> immutablePoints_{};
    std::vector<std::array<size_t, 7>> immutableSubdomains_{};
    std::vector<Border> immutableBorders_{};
    std::vector<size_t> numRefsOfLinesAboveX{};
    std::vector<size_t> numRefsOfLinesAboveY{};
    std::vector<size_t> numRefsOfLinesAboveZ{};


    void organizeBorders();

public:
    Mesh() { Logger::ConsoleOutput("Mesh declared, but it's empty.", NotificationColor::Warning); };
    ~Mesh() {};

    bool CheckData();
    void CommitData(std::vector<std::string>* data);

    inline bool isGenerated() const { return isGenerated_; }
    inline bool isDeclarated() const { return isDeclarated_; }

    inline size_t getLinesAmountX() const { return linesAmountX_; }
    inline size_t getLinesAmountY() const { return linesAmountY_; }
    inline size_t getLinesAmountZ() const { return linesAmountZ_; }

    inline std::vector<Point> getPoints() const { return points_; }
    inline std::vector<AreaRibs> getAreasAsRibs() const { return areasRibs_; }
    inline std::vector<RibRef> getRibsRefs() const { return referableRibs_; }
    inline std::vector<Border> getBorders() const { return borders_; }
    inline std::vector<BorderLine> getBorderRibs() const { return borderRibs_; }
    inline std::vector<std::array<size_t, 6>> getNewBorderRibs() const { return newBorders_; }
    inline std::vector<AreaInfo> getAreaInfo() const { return areasInfo_; }

    __declspec(property(get = getLinesAmountX)) size_t LinesAmountX;
    __declspec(property(get = getLinesAmountY)) size_t LinesAmountY;
    __declspec(property(get = getLinesAmountZ)) size_t LinesAmountZ;
    __declspec(property(get = isGeneraed)) bool IsGenerated;
    __declspec(property(get = isDeclarated)) bool IsDeclarated;
    __declspec(property(get = getPoints)) std::vector<Point> Points;
    __declspec(property(get = getAreaInfo)) std::vector<AreaInfo> AreasInfo;

    friend class MeshGenerator;
};