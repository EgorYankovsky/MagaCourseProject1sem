#pragma once

#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <iomanip>

const std::string standartInputPath = "Data\\input.txt";
const std::string inputPath1 = "Data\\input1.txt";
const std::string inputPath2 = "Data\\input2.txt";
const std::string inputPath3 = "Data\\input3.txt";
const std::string inputPath4 = "Data\\input4.txt";
const std::string inputPath5 = "Data\\input5.txt";
const std::string defaultOutputPointsPath = "Data\\generatedPoints.txt";
const std::string defaultOutputRibsPath = "Data\\generatedRibs.txt";
const std::string defaultOutputAreasPath = "Data\\generatedAreas.txt";

struct Point {
    double_t x = 0.0;
    double_t y = 0.0;
    double_t z = 0.0;
    Point() {};
    Point(double_t x_, double_t y_, double_t z_) : x(x_), y(y_), z(z_) {}
};

struct AreaPoints {
    uint32_t subdomainNum;
    std::array<uint32_t, 8> refs_{};
};

struct AreaRibs {
    size_t subdomainNum_ = 0;
    std::array<size_t, 12> refs_{};
    AreaRibs() {};
    AreaRibs(size_t subdomainNum, std::array<size_t, 12> refs) : subdomainNum_(subdomainNum), refs_(refs) {};
};

struct Rib {
    Point p1;
    Point p2;
    Rib() {};
    Rib(Point p1_, Point p2_) : p1(p1_), p2(p2_) {}
};

struct RibRef {
    size_t p1 = 0;
    size_t p2 = 0;
    RibRef() {};
    RibRef(size_t p1_, size_t p2_) : p1(p1_), p2(p2_) {}
};

struct AreaInfo {
    size_t subdomainNum_ = 0;
    double_t mu_ = 0.0;
    double_t sigma_ = 0.0;
    AreaInfo() {};
    AreaInfo(double_t mu, double_t sigma) : mu_(mu), sigma_(sigma) {};
};

struct Border {
    size_t type = 0;
    size_t formulaNum = 0;
    std::array<size_t, 6> refs{};
};