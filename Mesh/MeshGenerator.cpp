#include "MeshGenerator.h"

void MeshGenerator::Generate3DMesh(Mesh& mesh) {
    GenerateListOfPointsAboveZ(mesh);
    GenerateListOfPointsAboveY(mesh);
    GenerateListOfPointsAboveX(mesh);
    GenerateListOfRibs(mesh);
    GenerateListOfAreas(mesh);
    GenerateListOfBorders(mesh);
}

void MeshGenerator::GenerateListOfPointsAboveX(Mesh& mesh) {
    Logger::ConsoleOutput("Couldn't generate list of points above X axis", NotificationColor::Warning);
}

void MeshGenerator::GenerateListOfPointsAboveY(Mesh& mesh){
    Logger::ConsoleOutput("Couldn't generate list of points above Y axis", NotificationColor::Warning);
}

void MeshGenerator::GenerateListOfPointsAboveZ(Mesh& mesh) {
    auto lxy = mesh.LinesAmountX * mesh.LinesAmountY;
    for (size_t i(0); i < mesh.LinesAmountZ - 1; ++i) {
        std::vector<Point> pointsToInsert{};
        for (size_t j(0); j < mesh.LinesAmountY; ++j) {
            for (size_t k(0); k < mesh.LinesAmountX; ++k) {
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

    Logger::ConsoleOutput("Couldn't generate list of points above Z axis", NotificationColor::Warning);
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
