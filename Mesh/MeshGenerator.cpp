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
    for (size_t i(0); i < mesh.LinesAmountX; ++i) mesh.numRefsOfLinesAboveX.push_back(i);
    for (size_t i(0); i < mesh.LinesAmountY; ++i) mesh.numRefsOfLinesAboveY.push_back(i);
    for (size_t i(0); i < mesh.LinesAmountZ; ++i) mesh.numRefsOfLinesAboveZ.push_back(i);
    GenerateListOfPoints(mesh);
    GenerateListOfRibs(mesh);
    GenerateListOfAreas(mesh);
    GenerateListOfBorders(mesh);
}

int MeshGenerator::SelectAreaNum(Mesh& mesh, std::array<size_t, 12> arr) {
    auto sxy = mesh.LinesAmountX * mesh.LinesAmountY;

    auto p1 = mesh.referableRibs_[arr[0]].p1;
    auto p2 = mesh.referableRibs_[arr[11]].p2;

    auto lx0 = p1 % mesh.LinesAmountX;
    auto lx1 = p2 % mesh.LinesAmountX;
    
    auto ly0 = (p1 % sxy) / mesh.LinesAmountX;
    auto ly1 = (p2 % sxy) / mesh.LinesAmountX;
    
    auto lz0 = p1 / sxy;
    auto lz1 = p2 / sxy;

    //std::cout << lx0 << " " << lx1 << " " << ly0 << " " << ly1 << " " << lz0 << " " << lz1 << std::endl;

    for (const auto& area : mesh.subdomains_) 
        if (mesh.numRefsOfLinesAboveX[area[1]] <= lx0 and lx1 <= mesh.numRefsOfLinesAboveX[area[2]] and
            mesh.numRefsOfLinesAboveY[area[3]] <= ly0 and ly1 <= mesh.numRefsOfLinesAboveY[area[4]] and
            mesh.numRefsOfLinesAboveZ[area[5]] <= lz0 and lz1 <= mesh.numRefsOfLinesAboveZ[area[6]])
            return area[0];

    Logger::ConsoleOutput("Error during selection of area num", NotificationColor::Alert);
    return NAN;
}

// Try to optimize memory.
void MeshGenerator::GenerateListOfPoints(Mesh& mesh) {
    
    // Construct 3D area.
    area3D figure{};
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
            line1D lineToBuild{};
            lineToBuild = figure[k][j];
            size_t shift = mesh.LinesAmountX - 1;
            
            auto iterOnXRefs = mesh.numRefsOfLinesAboveX.begin() + 1;
            for (const auto& info : mesh.delimitersX_) {
                auto rightBorderIter = lineToBuild.end() - shift;
                auto amountOfDelimiters = info.first;
                auto coefficientOfDelimiter = info.second;

                *iterOnXRefs = *(iterOnXRefs - 1) + amountOfDelimiters;
                iterOnXRefs++;

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
        square2D squareToBuild{};
        squareToBuild = figure[k];
        size_t shift = mesh.LinesAmountY - 1;

        auto iterOnYRefs = mesh.numRefsOfLinesAboveY.begin() + 1;
        for (const auto info : mesh.delimitersY_) {
            auto rightBorderIter = squareToBuild.end() - shift;
            auto amountOfDelimiters = info.first;
            auto coefficientOfDelimiter = info.second;

            *iterOnYRefs = *(iterOnYRefs - 1) + amountOfDelimiters;
            iterOnYRefs++;

            line1D v0(mesh.LinesAmountX);
            line1D v1(mesh.LinesAmountX);
            std::copy((*(rightBorderIter - 1)).begin(), (*(rightBorderIter - 1)).end(), v0.begin());
            std::copy((*(rightBorderIter)).begin(), (*(rightBorderIter)).end(), v1.begin());

            square2D subSquareToBuild(amountOfDelimiters - 1);
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
    area3D areaToBuild{};
    areaToBuild = figure;
    size_t shift = mesh.LinesAmountZ - 1;

    auto iterOnZRefs = mesh.numRefsOfLinesAboveZ.begin() + 1;
    for (const auto info : mesh.delimitersZ_) {
        auto rightBorderIter = areaToBuild.end() - shift;
        auto amountOfDelimiters = info.first;
        auto coefficientOfDelimiter = info.second;

        *iterOnZRefs = *(iterOnZRefs - 1) + amountOfDelimiters;
        iterOnZRefs++;

        square2D s0(mesh.LinesAmountY); for (auto& line : s0) line.resize(mesh.LinesAmountX);
        square2D s1(mesh.LinesAmountY); for (auto& line : s1) line.resize(mesh.LinesAmountX);
        std::copy((*(rightBorderIter - 1)).begin(), (*(rightBorderIter - 1)).end(), s0.begin());
        std::copy((*(rightBorderIter)).begin(), (*(rightBorderIter)).end(), s1.begin());

        area3D subAreaToBuild(amountOfDelimiters - 1);
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

    size_t rx = mesh.LinesAmountX - 1;
    size_t ry = mesh.LinesAmountY - 1;

    size_t nx = mesh.LinesAmountX;
    size_t ny = mesh.LinesAmountY;

    size_t rxy = rx * ny + ry * nx;
    size_t nxy = nx * ny;
    size_t nz = mesh.LinesAmountZ;

    for (size_t k(0); k < nz - 1; ++k)
        for (size_t j(0); j < ny - 1; ++j)
            for (size_t i(0); i < nx - 1; ++i)
            {
                size_t curr = i + j * (nx + rx) + k * (rxy + nxy);
                std::array<size_t, 12> arrI = {
                    curr, curr + rx, curr + rx + 1, curr + rx + nx,
                    curr + rxy - j * rx, curr + rxy + 1 - j * rx, curr + rxy + nx - j * rx, curr + rxy + nx + 1 - j * rx,
                    curr + rxy + nxy, curr + rxy + nxy + rx, curr + rxy + nxy + rx + 1, curr + rxy + nxy + rx + nx };

                int areaNumber = SelectAreaNum(mesh, arrI);
                mesh.areasRibs_.emplace_back(areaNumber, arrI);
            }
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
