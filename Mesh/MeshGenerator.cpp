#include "MeshGenerator.h"

void Sort(std::vector<Point>& arr) {
    Point temp;
    bool isSorted(false);
    while (!isSorted) {
        isSorted = true;
        for (size_t counter(0); counter < arr.size() - 1; ++counter) {
            if (arr[counter].z > arr[counter + 1].z) {
                temp = arr[counter];
                arr[counter] = arr[counter + 1];
                arr[counter + 1] = temp;
                isSorted = false;
            }
            else {
                if (arr[counter].z == arr[counter + 1].z) {
                    if (arr[counter].y > arr[counter + 1].y) {
                        temp = arr[counter];
                        arr[counter] = arr[counter + 1];
                        arr[counter + 1] = temp;
                        isSorted = false;
                    }
                    else {
                        if (arr[counter].y == arr[counter + 1].y) {
                            if (arr[counter].x > arr[counter + 1].x) {
                                temp = arr[counter];
                                arr[counter] = arr[counter + 1];
                                arr[counter + 1] = temp;
                                isSorted = false;
                            }
                        }
                    }
                }
            }
        }
    }
}

void MeshGenerator::Generate3DMesh(Mesh& mesh) {
    assert(mesh.IsDeclarated);
    GenerateListOfPoints(mesh);
    GenerateListOfRibs(mesh);
    GenerateListOfAreas(mesh);
    GenerateListOfBorders(mesh);
}

// Try to optimize memory.
void MeshGenerator::GenerateListOfPoints(Mesh& mesh) {
    
    // Construct 3D area.
    std::vector<std::vector<std::vector<Point>>> figure{};
    figure.resize(mesh.LinesAmountZ);
    for (auto& square : figure) {
        square.resize(mesh.LinesAmountY);
        for (auto& line : square) line.resize(mesh.LinesAmountX);
    }

    // Fill 3D area.
    auto sxy = mesh.LinesAmountX * mesh.LinesAmountY;
    auto lx = mesh.LinesAmountX;
    for (size_t k(0); k < mesh.LinesAmountZ; ++k)
        for (size_t j(0); j < mesh.LinesAmountY; ++j)
            for (size_t i(0); i < mesh.LinesAmountX; ++i)
                figure[k][j][i] = mesh.immutablePoints_[k * sxy + j * lx + i];


    // Generation above X-axis.
    for (size_t k(0); k < mesh.LinesAmountZ; ++k) {
        for (size_t j(0); j < mesh.LinesAmountY; ++j) {
            std::vector<Point> lineToBuild{};
            lineToBuild = figure[k][j];
            size_t shift = mesh.LinesAmountX - 1;
            for (const auto& info : mesh.delimitersX_) {
                auto rightBorderIter = lineToBuild.end() - shift;
                auto amountOfDelimiters = info.first;
                auto coefficientOfDelimiter = info.second;

                double denum = 0.0;
                for (size_t ii(0); ii < amountOfDelimiters; ++ii)
                    denum += pow(coefficientOfDelimiter, ii);

                double x0 = (*(rightBorderIter - 1)).x;
                double x1 = (*rightBorderIter).x;

                double y0 = (*(rightBorderIter - 1)).y;
                double y1 = (*rightBorderIter).y;
                
                double z0 = (*(rightBorderIter - 1)).z;
                double z1 = (*rightBorderIter).z;

                double deltX = x1 - x0;
                double deltY = y1 - y0;
                double deltZ = z1 - z0;

                double xh = deltX / denum;
                double yh = deltY / denum;
                double zh = deltZ / denum;

                double multiplier = 0.0;
                for (size_t ii(0); ii < amountOfDelimiters - 1; ++ii) {
                    multiplier += pow(coefficientOfDelimiter, ii);
                    auto pointToInsert = Point(x0 + xh * multiplier, 
                                               y0 + yh * multiplier, 
                                               z0 + zh * multiplier);
                    lineToBuild.insert(lineToBuild.end() - shift, pointToInsert);
                }
                shift--;
            }
            figure[k][j] = lineToBuild;
        }
    }
    mesh.linesAmountX_ = figure[0][0].size();

    // Generation above Y-axis.
    for (size_t k(0); k < mesh.LinesAmountZ; ++k) {
        std::vector<std::vector<Point>> squareToBuild{};
        squareToBuild = figure[k];
        size_t shift = mesh.LinesAmountY - 1;
        for (const auto info : mesh.delimitersY_) {
            auto rightBorderIter = squareToBuild.end() - shift;
            auto amountOfDelimiters = info.first;
            auto coefficientOfDelimiter = info.second;

            std::vector<Point> v0(mesh.LinesAmountX);
            std::vector<Point> v1(mesh.LinesAmountX);
            std::copy((*(rightBorderIter - 1)).begin(), (*(rightBorderIter - 1)).end(), v0.begin());
            std::copy((*(rightBorderIter)).begin(), (*(rightBorderIter)).end(), v1.begin());

            std::vector<std::vector<Point>> subSquareToBuild(amountOfDelimiters - 1);
            for (auto& line : subSquareToBuild) line.resize(mesh.LinesAmountX);

            double denum = 0.0;
            for (size_t ii(0); ii < amountOfDelimiters; ++ii)
                denum += pow(coefficientOfDelimiter, ii);


            for (size_t i(0); i < mesh.LinesAmountX; ++i) {
                double x0 = v0[i].x;
                double x1 = v1[i].x;

                double y0 = v0[i].y;
                double y1 = v1[i].y;

                double z0 = v0[i].z;
                double z1 = v1[i].z;

                double deltX = x1 - x0;
                double deltY = y1 - y0;
                double deltZ = z1 - z0;

                double xh = deltX / denum;
                double yh = deltY / denum;
                double zh = deltZ / denum;

                double multiplier = 0.0;
                for (size_t ii(0); ii < amountOfDelimiters - 1; ++ii) {
                    multiplier += pow(coefficientOfDelimiter, ii);
                    auto pointToInsert = Point(x0 + xh * multiplier,
                                               y0 + yh * multiplier,
                                               z0 + zh * multiplier);

                    subSquareToBuild[ii][i] = pointToInsert;
                }
            }

            for (auto line : subSquareToBuild)
                squareToBuild.insert(squareToBuild.end() - shift, line);
            shift--;
        }
        figure[k] = squareToBuild;
    }
    mesh.linesAmountY_ = figure[0].size();

    // Generate above Z-axis.
    std::vector<std::vector<std::vector<Point>>> areaToBuild{};
    areaToBuild = figure;
    size_t shift = mesh.LinesAmountZ - 1;
    for (const auto info : mesh.delimitersZ_) {
        auto rightBorderIter = areaToBuild.end() - shift;
        auto amountOfDelimiters = info.first;
        auto coefficientOfDelimiter = info.second;

        std::vector<std::vector<Point>> s0(mesh.LinesAmountY); for (auto& line : s0) line.resize(mesh.LinesAmountX);
        std::vector<std::vector<Point>> s1(mesh.LinesAmountY); for (auto& line : s1) line.resize(mesh.LinesAmountX);
        std::copy((*(rightBorderIter - 1)).begin(), (*(rightBorderIter - 1)).end(), s0.begin());
        std::copy((*(rightBorderIter)).begin(), (*(rightBorderIter)).end(), s1.begin());

        std::vector<std::vector<std::vector<Point>>> subAreaToBuild(amountOfDelimiters - 1);
        for (auto& square : subAreaToBuild) {
            square.resize(mesh.LinesAmountY);
            for (auto& line : square)
                line.resize(mesh.LinesAmountX);
        }

        double denum = 0.0;
        for (size_t ii(0); ii < amountOfDelimiters; ++ii)
            denum += pow(coefficientOfDelimiter, ii);

        for (size_t j(0); j < mesh.LinesAmountY; ++j) {
            for (size_t i(0); i < mesh.LinesAmountX; ++i) {
                double x0 = s0[j][i].x;
                double x1 = s1[j][i].x;

                double y0 = s0[j][i].y;
                double y1 = s1[j][i].y;

                double z0 = s0[j][i].z;
                double z1 = s1[j][i].z;

                double deltX = x1 - x0;
                double deltY = y1 - y0;
                double deltZ = z1 - z0;

                double xh = deltX / denum;
                double yh = deltY / denum;
                double zh = deltZ / denum;

                double multiplier = 0.0;
                for (size_t ii(0); ii < amountOfDelimiters - 1; ++ii) {
                    multiplier += pow(coefficientOfDelimiter, ii);
                    auto pointToInsert = Point(x0 + xh * multiplier,
                        y0 + yh * multiplier,
                        z0 + zh * multiplier);

                    subAreaToBuild[ii][j][i] = pointToInsert;
                }
            }
        }
        for (auto square : subAreaToBuild)
            areaToBuild.insert(areaToBuild.end() - shift, square);
        shift--;
    }
    figure = areaToBuild;
    mesh.linesAmountZ_ = figure.size();


    // Convert to line-format.
    mesh.linesAmountZ_ = figure.size();
    mesh.points_.resize(mesh.linesAmountX_ * mesh.linesAmountY_ * mesh.linesAmountZ_);
    sxy = mesh.LinesAmountX * mesh.LinesAmountY;
    lx = mesh.LinesAmountX;
    for (size_t k(0); k < mesh.LinesAmountZ; ++k)
        for (size_t j(0); j < mesh.LinesAmountY; ++j)
            for (size_t i(0); i < mesh.LinesAmountX; ++i)
                mesh.points_[k * sxy + j * lx + i] = figure[k][j][i];
}

void MeshGenerator::GenerateListOfAreas(Mesh& mesh) {
    Logger::ConsoleOutput("Couldn't generate list of areas", NotificationColor::Warning);
}

void MeshGenerator::GenerateListOfRibs(Mesh& mesh) {
    auto nx = mesh.LinesAmountX;
    auto ny = mesh.LinesAmountY;
    auto nz = mesh.LinesAmountZ;
    auto nxny = mesh.LinesAmountX * mesh.LinesAmountY;
    for (int k = 0; k < mesh.LinesAmountZ; k++)
    {
        for (int j = 0; j < mesh.LinesAmountY; j++)
        {
            for (int i = 0; i < mesh.LinesAmountX - 1; i++)
                mesh.referableRibs_.emplace_back(k * nxny + nx * j + i, k * nxny + nx * j + i + 1);
            if (j != mesh.LinesAmountY - 1)
                for (int i = 0; i < nx; i++)
                    mesh.referableRibs_.emplace_back(k * nxny + nx * j + i, k * nxny + nx * (j + 1) + i);
        }
        if (k != mesh.LinesAmountZ - 1)
            for (int j = 0; j < mesh.LinesAmountY; j++)
                for (int i = 0; i < mesh.LinesAmountX; i++)
                    mesh.referableRibs_.emplace_back(k * nxny + nx * j + i, (k + 1) * nxny + nx * j + i);
    }
    Logger::ConsoleOutput("Ribs array generated.", NotificationColor::Passed);
}

void MeshGenerator::GenerateListOfBorders(Mesh& mesh) {
    Logger::ConsoleOutput("Couldn't generate list of borders", NotificationColor::Warning);
}
