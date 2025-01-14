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

}

void MeshGenerator::GenerateListOfAreas(Mesh& mesh) {
    Logger::ConsoleOutput("Couldn't generate list of areas", NotificationColor::Warning);
}

void MeshGenerator::GenerateListOfRibs(Mesh& mesh) {
    Logger::ConsoleOutput("Couldn't generate list of ribs", NotificationColor::Warning);
}

void MeshGenerator::GenerateListOfBorders(Mesh& mesh) {
    Logger::ConsoleOutput("Couldn't generate list of borders", NotificationColor::Warning);
}
