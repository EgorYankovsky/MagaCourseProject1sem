#include "Mesh.h"

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

                // Get delimeters amount and it's coefficient.
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

                // Get delimeters amount and it's coefficient.
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

                // Get delimeters amount and it's coefficient.
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
}

void Mesh::Generate() {
    assert(isDeclarated_);
    generateAboveZ();
    generateAboveY();
    generateAboveX();
}

bool Mesh::CheckData() {
    if (linesAmountX_ * linesAmountY_ * linesAmountZ_ != points_.size()) return false;
    if (linesAmountX_ - 1 != delimetersX_.size()) return false;
    if (linesAmountY_ - 1 != delimetersY_.size()) return false;
    if (linesAmountZ_ - 1 != delimetersZ_.size()) return false;
    uint32_t maxLineX = 0;
    uint32_t maxLineY = 0;
    uint32_t maxLineZ = 0;
    for (const auto& subdomain : subdomains_) {
        maxLineX = maxLineX < subdomain[2] ? subdomain[2] : maxLineX;
        maxLineY = maxLineY < subdomain[4] ? subdomain[4] : maxLineY;
        maxLineZ = maxLineZ < subdomain[6] ? subdomain[6] : maxLineZ;
    }
    if (linesAmountX_ - 1 != maxLineX) return false;
    if (linesAmountY_ - 1 != maxLineY) return false;
    if (linesAmountZ_ - 1 != maxLineZ) return false;
    return true;
}

void Mesh::FileWriteGeneratedPoints(std::string fileName) {
    std::ofstream fout(fileName);
    for (auto& point : points_)
        fout << std::setprecision(15) << std::scientific << point.x << " " << point.y << " " << point.z << std::endl;
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

    // Read delimetes above X.
    _mesh.delimetersX_.resize(_mesh.linesAmountX_ - 1);
    for (auto& delimeter : _mesh.delimetersX_) fin >> delimeter.first >> delimeter.second;
    
    // Read delimetes above Y.
    _mesh.delimetersY_.resize(_mesh.linesAmountY_ - 1);
    for (auto& delimeter : _mesh.delimetersY_) fin >> delimeter.first >> delimeter.second;
    
    // Read delimetes above Z.
    _mesh.delimetersZ_.resize(_mesh.linesAmountZ_ - 1);
    for (auto& delimeter : _mesh.delimetersZ_) fin >> delimeter.first >> delimeter.second;


    _mesh.isDeclarated_ = true;
    fin.close();
}
