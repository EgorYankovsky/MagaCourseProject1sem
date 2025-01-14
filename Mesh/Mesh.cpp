#include "Mesh.h"


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
                    //Sort(pointsToInsert);
                    points_.insert(points_.end(), pointsToInsert.begin(), pointsToInsert.end());
                }
            }
        }
    }
    //Sort(points_);
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
                //Sort(pointsToInsert);
                points_.insert(points_.end(), pointsToInsert.begin(), pointsToInsert.end());
            }
        }
    }
    //Sort(points_);
    linesAmountY_ = points_.size() / (linesAmountX_ * linesAmountZ_);
    Logger::ConsoleOutput("Generated above Y.", NotificationColor::Warning);
    Logger::ConsoleOutput("Possible mistakes during generation above Y.", NotificationColor::Alert);
}

// Fully works.
void Mesh::generatePoints() {
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
    isDeclarated_ = true;
    return true;
}

void Mesh::CommitData(std::vector<std::string>* data) {
    // Select first item of vector.
    auto currentItem = data->begin();
    
    // Commit lines amount above X,Y,Z axis.
    linesAmountX_ = std::stoul(*currentItem); currentItem++;
    linesAmountY_ = std::stoul(*currentItem); currentItem++;
    linesAmountZ_ = std::stoul(*currentItem); currentItem++;

    // Commit area description.
    for (size_t i(0); i < linesAmountX_ * linesAmountY_ * linesAmountZ_; ++i) {
        points_.emplace_back(std::stod(*currentItem),           // X 
                             std::stod(*(currentItem + 1)),     // Y
                             std::stod(*(currentItem + 2)));    // Z
        currentItem += 3;
    }
    
    // Commit unique areas description.
    subdomainsAmount_ = std::stoul(*currentItem); currentItem++;
    for (size_t i(0); i < subdomainsAmount_; ++i) {

        std::array<size_t, 7> currentArray = { std::stoul(*currentItem),            // Formula num.
                                               std::stoul(*(currentItem + 1)),      // X0.
                                               std::stoul(*(currentItem + 2)),      // X1.
                                               std::stoul(*(currentItem + 3)),      // Y0.
                                               std::stoul(*(currentItem + 4)),      // Y1.
                                               std::stoul(*(currentItem + 5)),      // Z0.
                                               std::stoul(*(currentItem + 6)) };    // Z1.
        subdomains_.push_back(currentArray);
        currentItem += 7;
    }

    for (size_t i(0); i < subdomainsAmount_; ++i) {
        areasInfo_.emplace_back(std::stoul(*currentItem),           // Area num.
                                std::stod(*(currentItem + 1)),      // Mu_i.
                                std::stod(*(currentItem + 2)));     // Sigma_i.
        currentItem += 3;
    }

    for (size_t i(0); i < linesAmountX_ - 1; ++i) {
        delimetersX_.emplace_back(std::stoul(*currentItem),         // Delimiters amount above X.
                                  std::stod(*(currentItem + 1)));   // Delimiters coefficient above X.
        currentItem += 2;
    }

    for (size_t i(0); i < linesAmountX_ - 1; ++i) {
        delimetersY_.emplace_back(std::stoul(*currentItem),         // Delimiters amount above Y.
            std::stod(*(currentItem + 1)));                         // Delimiters coefficient above Y.
        currentItem += 2;
    }

    for (size_t i(0); i < linesAmountX_ - 1; ++i) {
        delimetersZ_.emplace_back(std::stoul(*currentItem),         // Delimiters amount above Z.
            std::stod(*(currentItem + 1)));                         // Delimiters coefficient above Z.
        currentItem += 2;
    }

    bordersAmount_ = std::stoul(*currentItem); currentItem++;
    for (size_t i(0); i < bordersAmount_; ++i) {
        borders_.emplace_back(std::stoul(*currentItem),             // Border type.
                              std::stoul(*(currentItem + 1)),       // Border formula num.
                              std::stoul(*(currentItem + 2)),       // X0.
                              std::stoul(*(currentItem + 3)),       // X1.
                              std::stoul(*(currentItem + 4)),       // Y0.
                              std::stoul(*(currentItem + 5)),       // Y1.
                              std::stoul(*(currentItem + 6)),       // Z0.
                              std::stoul(*(currentItem + 7)));      // Z1.
        currentItem += 8;
    }
}