#include "Mesh.h"

std::string Mesh::defaultOutputPointsPath = "Data\\Generated\\generatedPoints.txt";
std::string Mesh::defaultOutputRibsPath = "Data\\Generated\\generatedRibs.txt";
std::string Mesh::defaultOutputAreasPath = "Data\\Generated\\generatedAreas.txt";

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

// Fully works, but not optimal.
void Mesh::generateAboveX() {
    auto lxy = linesAmountX_ * linesAmountY_;
    for (size_t i(0); i < linesAmountZ_; ++i) {
        for (size_t j(0); j < linesAmountY_; ++j) {
            for (size_t k(0); k < linesAmountX_ - 1; ++k) {
                std::vector<Point> pointsToInsert{};
                auto& pnt1 = points_[i * lxy + j * linesAmountX_ + k];
                auto& pnt2 = points_[i * lxy + j * linesAmountX_ + k + 1];

                // Find difference above axis.
                double_t dx = pnt2.x - pnt1.x;
                double_t dy = pnt2.y - pnt1.y;
                double_t dz = pnt2.z - pnt1.z;

                // Get delimiters amount and it's coefficient.
                auto delimAmount = delimetersX_[k].first;
                auto delimCoef = delimetersX_[k].second;

                double_t denominator(0.0);
                for (size_t ii(0); ii < delimAmount; ++ii)
                    denominator += pow(delimCoef, ii);

                double_t hx0 = dx / denominator;
                double_t hy0 = dy / denominator;
                double_t hz0 = dz / denominator;

                double_t mnoz(0.0);
                for (size_t ii(0); ii < delimAmount - 1; ++ii) {
                    mnoz += pow(delimCoef, ii);
                    pointsToInsert.push_back(Point(
                        pnt1.x + hx0 * mnoz,
                        pnt1.y + hy0 * mnoz,
                        pnt1.z + hz0 * mnoz
                    ));
                }
                if (pointsToInsert.size() > 0) {
                    Sort(pointsToInsert);
                    points_.insert(points_.end(), pointsToInsert.begin(), pointsToInsert.end());
                }
            }
        }
    }
    Sort(points_);
    linesAmountX_ = points_.size() / (linesAmountY_ * linesAmountZ_);
    Logger::ConsoleOutput("Generated above X.", NotificationColor::Warning);
}

// Fully works, but not optimal.
void Mesh::generateAboveY() {
    auto lxy = linesAmountX_ * linesAmountY_;
    for (size_t i(0); i < linesAmountZ_; ++i) {
        for (size_t j(0); j < linesAmountY_ - 1; ++j) {
            std::vector<Point> pointsToInsert{};
            for (size_t k(0); k < linesAmountX_; ++k) {
                auto& pnt1 = points_[i * lxy + j * linesAmountX_ + k];
                auto& pnt2 = points_[i * lxy + (j + 1) * linesAmountX_ + k];

                // Find difference above axis.
                double_t dx = pnt2.x - pnt1.x;
                double_t dy = pnt2.y - pnt1.y;
                double_t dz = pnt2.z - pnt1.z;

                // Get delimiters amount and it's coefficient.
                auto delimAmount = delimetersY_[j].first;
                auto delimCoef = delimetersY_[j].second;

                double_t denominator(0.0);
                for (size_t ii(0); ii < delimAmount; ++ii)
                    denominator += pow(delimCoef, ii);

                double_t hx0 = dx / denominator;
                double_t hy0 = dy / denominator;
                double_t hz0 = dz / denominator;

                double_t mnoz(0.0);
                for (size_t ii(0); ii < delimAmount - 1; ++ii) {
                    mnoz += pow(delimCoef, ii);
                    pointsToInsert.push_back(Point(
                        pnt1.x + hx0 * mnoz,
                        pnt1.y + hy0 * mnoz,
                        pnt1.z + hz0 * mnoz
                    ));
                }
            }
            if (pointsToInsert.size() > 0) {
                Sort(pointsToInsert);
                points_.insert(points_.end(), pointsToInsert.begin(), pointsToInsert.end());
            }
        }
    }
    Sort(points_);
    linesAmountY_ = points_.size() / (linesAmountX_ * linesAmountZ_);
    Logger::ConsoleOutput("Generated above Y.", NotificationColor::Warning);
    Logger::ConsoleOutput("Possible mistakes during generation above Y.", NotificationColor::Alert);
}

// Fully works.
void Mesh::generateAboveZ() {
    auto lxy = linesAmountX_ * linesAmountY_;
    for (size_t i(0); i < linesAmountZ_ - 1; ++i) {
        std::vector<Point> pointsToInsert{};
        for (size_t j(0); j < linesAmountY_; ++j) {
            for (size_t k(0); k < linesAmountX_; ++k) {
                auto& pnt1 = immutablePoints_[i * lxy + j * linesAmountX_ + k];
                auto& pnt2 = immutablePoints_[(i + 1) * lxy + j * linesAmountX_ + k];
                
                // Find difference above axis.
                double_t dx = pnt2.x - pnt1.x;
                double_t dy = pnt2.y - pnt1.y;
                double_t dz = pnt2.z - pnt1.z;

                // Get delimiters amount and it's coefficient.
                auto delimAmount = delimetersZ_[i].first;
                auto delimCoef = delimetersZ_[i].second;
                
                double_t denominator(0.0);
                for (size_t ii(0); ii < delimAmount; ++ii) 
                    denominator += pow(delimCoef, ii);

                double_t hx0 = dx / denominator;
                double_t hy0 = dy / denominator;
                double_t hz0 = dz / denominator;

                double_t mnoz(0.0);
                for (size_t ii(0); ii < delimAmount - 1; ++ii) {
                    mnoz += pow(delimCoef, ii);
                    pointsToInsert.push_back(Point(
                        pnt1.x + hx0 * mnoz,
                        pnt1.y + hy0 * mnoz,
                        pnt1.z + hz0 * mnoz
                    ));
                }
            }
        }
        if (pointsToInsert.size() > 0) {
            Sort(pointsToInsert);
            points_.insert(points_.begin() + (i + 1) * lxy, pointsToInsert.begin(), pointsToInsert.end());
        }
    }
    linesAmountZ_ = points_.size() / lxy;
    Logger::ConsoleOutput("Generated above Z.", NotificationColor::Passed);
}

// Fully works.
void Mesh::generatePoints() {
    generateAboveZ();
    generateAboveY();
    generateAboveX();
    Logger::ConsoleOutput("Points generated.", NotificationColor::Passed);
}

void Mesh::generateRibsArray() {
    auto nx = linesAmountX_;
    auto ny = linesAmountY_;
    auto nz = linesAmountZ_;
    auto nxny = linesAmountX_ * linesAmountY_;
    for (int k = 0; k < linesAmountZ_; k++)
    {
        for (int j = 0; j < linesAmountY_; j++)
        {
            for (int i = 0; i < linesAmountX_ - 1; i++)
                referableRibs_.emplace_back(k * nxny + nx * j + i, k * nxny + nx * j + i + 1);
            if (j != linesAmountY_ - 1)
                for (int i = 0; i < nx; i++)
                    referableRibs_.emplace_back(k * nxny + nx * j + i, k * nxny + nx * (j + 1) + i);
        }
        if (k != linesAmountZ_ - 1)
            for (int j = 0; j < linesAmountY_; j++)
                for (int i = 0; i < linesAmountX_; i++)
                    referableRibs_.emplace_back(k * nxny + nx * j + i, (k + 1) * nxny + nx * j + i);
    }
    Logger::ConsoleOutput("Ribs array generated.", NotificationColor::Passed);
}

void Mesh::generateAreasArray() {
    
    size_t rx = linesAmountX_ - 1;
    size_t ry = linesAmountY_ - 1;

    size_t nx = linesAmountX_;
    size_t ny = linesAmountY_;

    size_t rxy = rx * ny + ry * nx;
    size_t nxy = nx * ny;
    size_t nz = linesAmountZ_;

    for (size_t k(0); k < nz - 1; ++k)
        for (size_t j(0); j < ny - 1; ++j)
            for (size_t i(0); i < nx - 1; ++i)
            {
                size_t curr = i + j * (nx + rx) + k * (rxy + nxy);
                std::array<size_t, 12> arrI = { 
                    curr, curr + rx, curr + rx + 1, curr + rx + nx,
                    curr + rxy - j * rx, curr + rxy + 1 - j * rx, curr + rxy + nx - j * rx, curr + rxy + nx + 1 - j * rx,
                    curr + rxy + nxy, curr + rxy + nxy + rx, curr + rxy + nxy + rx + 1, curr + rxy + nxy + rx + nx };
                areasRibs_.emplace_back(1, arrI);
            }
    Logger::ConsoleOutput("Area number doesn't commits!", NotificationColor::Alert);
    Logger::ConsoleOutput("Areas array generated.", NotificationColor::Warning);
}

void Mesh::generateBorderArray() {

}

void Mesh::Generate() {
    assert(isDeclarated_);
    generatePoints();
    generateRibsArray();
    generateAreasArray();
}

bool Mesh::CheckData() {
    if (linesAmountX_ * linesAmountY_ * linesAmountZ_ != points_.size()) return false;
    if (linesAmountX_ - 1 != delimetersX_.size()) return false;
    if (linesAmountY_ - 1 != delimetersY_.size()) return false;
    if (linesAmountZ_ - 1 != delimetersZ_.size()) return false;
    size_t maxLineX = 0;
    size_t maxLineY = 0;
    size_t maxLineZ = 0;
    for (const auto& subdomain : subdomains_) {
        maxLineX = maxLineX < subdomain[2] ? subdomain[2] : maxLineX;
        maxLineY = maxLineY < subdomain[4] ? subdomain[4] : maxLineY;
        maxLineZ = maxLineZ < subdomain[6] ? subdomain[6] : maxLineZ;
    }
    if (linesAmountX_ - 1 != maxLineX) return false;
    if (linesAmountY_ - 1 != maxLineY) return false;
    if (linesAmountZ_ - 1 != maxLineZ) return false;
    Logger::ConsoleOutput("Mesh checked and declared.", NotificationColor::Passed);
    return true;
}

void Mesh::FileWriteGeneratedPoints(std::string fileName) {
    std::ofstream fout(fileName);
    fout << points_.size() << std::endl; 
    for (auto& point : points_)
        fout << std::setprecision(15) << std::scientific << point.x << " " << point.y << " " << point.z << std::endl;
    fout.close();
}

void Mesh::FileWriteGeneratedRibs(std::string fileName) {
    std::ofstream fout(fileName);
    fout << referableRibs_.size() << std::endl;
    for (auto& refRib : referableRibs_)
        fout << refRib.p1 << " " << refRib.p2 << std::endl;
    fout.close();
}

void Mesh::FileWriteGeneratedAreas(std::string fileName) {
    std::ofstream fout(fileName);
    fout << areasRibs_.size() << std::endl;
    for (auto& area : areasRibs_) {
        fout << area.subdomainNum_ << " ";
        for (auto& value : area.refs_) fout << value << " ";
        fout << std::endl;
    }
    fout.close();
}

void ReadData(Mesh& _mesh, std::string inputData) {
    // Read nodes (lines) amount above X,Y,Z.
    std::ifstream fin(inputData);
    fin >> _mesh.linesAmountX_ >> _mesh.linesAmountY_ >> _mesh.linesAmountZ_;

    // Read all points.
    _mesh.points_.resize(_mesh.linesAmountX_ * _mesh.linesAmountY_ * _mesh.linesAmountZ_);
    for (auto& point : _mesh.points_) fin >> point.x >> point.y >> point.z;
    _mesh.immutablePoints_ = _mesh.points_;

    // Read all subdomains.
    fin >> _mesh.subdomainsAmount_;
    _mesh.subdomains_.resize(_mesh.subdomainsAmount_);
    for (auto& subdomain : _mesh.subdomains_) fin >> subdomain[0] >> subdomain[1] >> subdomain[2] >>
        subdomain[3] >> subdomain[4] >> subdomain[5] >> subdomain[6];
    _mesh.areasInfo_.resize(_mesh.subdomainsAmount_);
    for (auto& subdomainInfo : _mesh.areasInfo_)
        fin >> subdomainInfo.subdomainNum_ >> subdomainInfo.mu_ >> subdomainInfo.sigma_;

    // Read delimiters above X.
    _mesh.delimetersX_.resize(_mesh.linesAmountX_ - 1);
    for (auto& delimeter : _mesh.delimetersX_) fin >> delimeter.first >> delimeter.second;
    
    // Read delimiters above Y.
    _mesh.delimetersY_.resize(_mesh.linesAmountY_ - 1);
    for (auto& delimeter : _mesh.delimetersY_) fin >> delimeter.first >> delimeter.second;
    
    // Read delimiters above Z.
    _mesh.delimetersZ_.resize(_mesh.linesAmountZ_ - 1);
    for (auto& delimeter : _mesh.delimetersZ_) fin >> delimeter.first >> delimeter.second;

    // Read information about borders.
    size_t bordersAmount(0);
    fin >> bordersAmount;
    _mesh.borders_.resize(bordersAmount);
    for (auto& border : _mesh.borders_)
        fin >> border.type >> border.formulaNum >>
            border.refs[0] >> border.refs[1] >> border.refs[2] >>
            border.refs[3] >> border.refs[4] >> border.refs[5];
    _mesh.isDeclarated_ = true;
    Logger::ConsoleOutput("Mesh completely read.", NotificationColor::Passed);
    fin.close();
}
