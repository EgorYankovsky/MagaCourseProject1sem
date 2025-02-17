#include "Mesh.h"

void Mesh::organizeBorders() {
    size_t bordersIISwifted(0); size_t bordersIIISwifted(0);
    for (auto& border : borders_) {
        if (border.type_ == 2) {
            std::iter_swap(borders_.begin() + bordersIISwifted, &border);
            bordersIISwifted++; }
        if (border.type_ == 3) {
            std::iter_swap(borders_.begin() + bordersIISwifted + bordersIIISwifted, &border);
            bordersIIISwifted++; } }
}

bool Mesh::CheckData() {
    if (linesAmountX_ * linesAmountY_ * linesAmountZ_ != points_.size()) return false;
    if (linesAmountX_ - 1 != delimitersX_.size()) return false;
    if (linesAmountY_ - 1 != delimitersY_.size()) return false;
    if (linesAmountZ_ - 1 != delimitersZ_.size()) return false;
    size_t maxLineX = 0; size_t maxLineY = 0; size_t maxLineZ = 0;
    for (const auto& subdomain : subdomains_) {
        maxLineX = maxLineX < subdomain[2] ? subdomain[2] : maxLineX;
        maxLineY = maxLineY < subdomain[4] ? subdomain[4] : maxLineY;
        maxLineZ = maxLineZ < subdomain[6] ? subdomain[6] : maxLineZ; }
    if (linesAmountX_ - 1 != maxLineX) return false;
    if (linesAmountY_ - 1 != maxLineY) return false;
    if (linesAmountZ_ - 1 != maxLineZ) return false;
    Logger::ConsoleOutput("Mesh checked and declared.", NotificationColor::Passed);
    isDeclarated_ = true; return true;
}

void Mesh::CommitData(std::vector<std::string>* data) {
    auto currentItem = data->begin(); // Select first item of vector.
    // Commit lines amount above X,Y,Z axis.
    linesAmountX_ = std::stoul(*currentItem); currentItem++;
    linesAmountY_ = std::stoul(*currentItem); currentItem++;
    linesAmountZ_ = std::stoul(*currentItem); currentItem++;
    // Commit area description.
    for (size_t i(0); i < linesAmountX_ * linesAmountY_ * linesAmountZ_; ++i) {
        points_.emplace_back(std::stod(*currentItem),           // X 
                             std::stod(*(currentItem + 1)),     // Y
                             std::stod(*(currentItem + 2)));    // Z
        currentItem += 3; }
    immutablePoints_ = points_;
    // Commit unique areas description.
    subdomainsAmount_ = std::stoul(*currentItem); currentItem++;
    for (size_t i(0); i < subdomainsAmount_; ++i) {
        std::array<size_t, 7> currentArray = { std::stoul(*currentItem), // Formula num.
            std::stoul(*(currentItem + 1)),      // X0.
            std::stoul(*(currentItem + 2)),      // X1.
            std::stoul(*(currentItem + 3)),      // Y0.
            std::stoul(*(currentItem + 4)),      // Y1.
            std::stoul(*(currentItem + 5)),      // Z0.
            std::stoul(*(currentItem + 6)) };    // Z1.
        subdomains_.push_back(currentArray);
        currentItem += 7; }
    immutableSubdomains_ = subdomains_;
    // Commit unique areas coefficients description.
    for (size_t i(0); i < subdomainsAmount_; ++i) {
        areasInfo_.emplace_back(std::stoul(*currentItem),           // Area num.
                                std::stod(*(currentItem + 1)),      // Mu_i.
                                std::stod(*(currentItem + 2)));     // Sigma_i.
        currentItem += 3; }
    // Commit delimiters above X description.
    for (size_t i(0); i < linesAmountX_ - 1; ++i) {
        delimitersX_.emplace_back(std::stoul(*currentItem),
            std::stod(*(currentItem + 1)));
        currentItem += 2; }
    // Commit delimiters above Y description.
    for (size_t i(0); i < linesAmountY_ - 1; ++i) {
        delimitersY_.emplace_back(std::stoul(*currentItem),
            std::stod(*(currentItem + 1)));
        currentItem += 2; }
    // Commit delimiters above Z description.
    for (size_t i(0); i < linesAmountZ_ - 1; ++i) {
        delimitersZ_.emplace_back(std::stoul(*currentItem),
            std::stod(*(currentItem + 1)));
        currentItem += 2; }
    // Commit information about borders.
    bordersAmount_ = std::stoul(*currentItem); currentItem++;
    for (size_t i(0); i < bordersAmount_; ++i) {
        borders_.emplace_back(std::stoul(*currentItem), // Border type.
            std::stoul(*(currentItem + 1)),       // Border formula num.
            std::stoul(*(currentItem + 2)),       // X0.
            std::stoul(*(currentItem + 3)),       // X1.
            std::stoul(*(currentItem + 4)),       // Y0.
            std::stoul(*(currentItem + 5)),       // Y1.
            std::stoul(*(currentItem + 6)),       // Z0.
            std::stoul(*(currentItem + 7)));      // Z1.
        currentItem += 8; }
    organizeBorders(); immutableBorders_ = borders_;
}