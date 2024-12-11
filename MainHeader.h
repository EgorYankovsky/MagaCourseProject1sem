#pragma once

#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <iomanip>

const std::string standartInputPath = "Data\\input.txt";
const std::string defaultOutputPointsPath = "Data\\generatedPoints.txt";


struct Point {
    double_t x = 0.0;
    double_t y = 0.0;
    double_t z = 0.0;
    Point() {};
    Point(double_t x_, double_t y_, double_t z_) : x(x_), y(y_), z(z_) {}
};

struct Area {
    uint32_t subdomainNum;
    uint32_t p1;
    uint32_t p2;
    uint32_t p3;
    uint32_t p4;
    uint32_t p5;
    uint32_t p6;
    uint32_t p7;
    uint32_t p8;
};

struct Border {
    uint32_t borderNum;
    uint32_t p1;
    uint32_t p2;
    uint32_t p3;
    uint32_t p4;
};