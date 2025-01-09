#pragma once

#include "..\\DataTypes.h"
#include "..\\Logger\\Logger.h"

#include <cassert>
#include <cmath>

class Mesh {
private:
    // Based.
    bool isGenerated_ = false;
    bool isDeclarated_ = false;
    size_t linesAmountX_ = 0;
    size_t linesAmountY_ = 0;
    size_t linesAmountZ_ = 0;
    std::vector<AreaPoints> areasPoints_{};
    std::vector<AreaRibs> areasRibs_{};
    std::vector<Border> borders_{};
    std::vector<Point> points_{};
    std::vector<Point> immutablePoints_{};
    std::vector<RibRef> referableRibs_{};
    std::vector<AreaInfo> areasInfo_{};

    // Additional.
    size_t subdomainsAmount_ = 0;
    std::vector<std::array<size_t, 7>> subdomains_{};
    std::vector<std::pair<size_t, double_t>> delimetersX_{};
    std::vector<std::pair<size_t, double_t>> delimetersY_{};
    std::vector<std::pair<size_t, double_t>> delimetersZ_{};
    
    static std::string defaultOutputPointsPath;// = "Data\\Generated\\generatedPoints.txt";
    static std::string defaultOutputRibsPath;// = "Data\\Generated\\generatedRibs.txt";
    static std::string defaultOutputAreasPath;// = "Data\\Generated\\generatedAreas.txt";
    
    void generateAboveX();
    void generateAboveY();
    void generateAboveZ();
    void generatePoints();
    void generateRibsArray();
    void generateAreasArray();
    void generateBorderArray();

public:
    Mesh() { Logger::ConsoleOutput("Mesh declared, but empty.", NotificationColor::Warning); };
    ~Mesh() {};

    void Generate();
    bool CheckData();

    inline bool isGenerated() const { return isGenerated_; }
    inline bool isDeclarated() const { return isDeclarated_; }

    void FileWriteGeneratedPoints(std::string fileName = defaultOutputPointsPath);
    void FileWriteGeneratedRibs(std::string fileName = defaultOutputRibsPath);
    void FileWriteGeneratedAreas(std::string fileName = defaultOutputAreasPath);

    inline size_t getLinesAmountX() const { return linesAmountX_; }
    inline size_t getLinesAmountY() const { return linesAmountY_; }
    inline size_t getLinesAmountZ() const { return linesAmountZ_; }
    inline std::vector<Point> getPoints() const { return points_; }

    __declspec(property(get = getLinesAmountX)) size_t LinesAmountX;
    __declspec(property(get = getLinesAmountY)) size_t LinesAmountY;
    __declspec(property(get = getLinesAmountZ)) size_t LinesAmountZ;
    __declspec(property(get = isGeneraed)) bool IsGenerated;
    __declspec(property(get = isDeclarated)) bool IsDeclarated;
    __declspec(property(get = getPoints)) std::vector<Point> Points;

    friend void Sort(std::vector<Point>& arr);
    friend void ReadData(Mesh& _mesh, std::string inputData);
};