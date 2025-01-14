#pragma once

#include "..\\DataTypes.h"
#include "..\\Logger\\Logger.h"

#include <cassert>
#include <cmath>
#include <string>

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

    std::vector<std::pair<size_t, double_t>> delimetersX_{};
    std::vector<std::pair<size_t, double_t>> delimetersY_{};
    std::vector<std::pair<size_t, double_t>> delimetersZ_{};

    size_t bordersAmount_ = 0;
    std::vector<Border> borders_{};


    std::vector<AreaPoints> areasPoints_{};
    std::vector<RibRef> referableRibs_{};

    // Additional.    
    std::vector<Point> immutablePoints_{};
    std::vector<Border> immutableBorders_{};

    void generateAboveX();
    void generateAboveY();

    void generatePoints();
    void generateRibsArray();
    void generateAreasArray();
    void generateBorderArray();

public:
    Mesh() { Logger::ConsoleOutput("Mesh declared, but empty.", NotificationColor::Warning); };
    ~Mesh() {};

    void Generate();
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

    __declspec(property(get = getLinesAmountX)) size_t LinesAmountX;
    __declspec(property(get = getLinesAmountY)) size_t LinesAmountY;
    __declspec(property(get = getLinesAmountZ)) size_t LinesAmountZ;
    __declspec(property(get = isGeneraed)) bool IsGenerated;
    __declspec(property(get = isDeclarated)) bool IsDeclarated;
    __declspec(property(get = getPoints)) std::vector<Point> Points;

    friend class MeshGenerator;
};