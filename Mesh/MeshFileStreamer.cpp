#include "MeshFileStreamer.h"

#include <stdlib.h>

std::string MeshFileStreamer::pathToText = "Data/Generated/Text/";
std::string MeshFileStreamer::pathToBin = "Data/Generated/Bin/";

void MeshFileStreamer::Read(Mesh& _mesh, std::string path) {
    std::ifstream file(path, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        Logger::ConsoleOutput("Such file doesn't exist.", NotificationColor::Alert);
        exit(-1);
    }
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        Logger::ConsoleOutput("Error during file fetching", NotificationColor::Alert);
        exit(-1);
    }
    file.close();

    std::vector<std::string> numbers;
    /*
    * -1 -> no reading
    *  0 -> delayed reading
    *  1 -> read word
    *  2 -> read number
    */
    int readMode = -1;
    std::string word;
    std::string mes;
    for (auto character : buffer) {
        switch (character) {
        case 10: // \n
            if (readMode == 2 && !word.empty()) {
                numbers.push_back(word);
                word = std::string();
            }
            break;
        case 13: // \r
            break;
        case 91: // [
            readMode = 1;
            if (!word.empty()) {
                numbers.push_back(word);
                word = std::string();
            }
            break;
        case 93: // ]
            readMode = 2;
            mes = word + " already read!";
            Logger::ConsoleOutput(mes, NotificationColor::Passed);
            word = std::string();
            break;
        case 32: // (space)
            if (readMode == 1) word += character;
            if (readMode == 2 && !word.empty()) {
                numbers.push_back(word);
                word = std::string();
            }
            break;
        default:
            if (readMode == -1) {
                Logger::ConsoleOutput("Failure during reading", NotificationColor::Passed);
                exit(-1);
            }
            if (readMode == 0) {
                mes = "Unexpected character: " + character;
                Logger::ConsoleOutput(mes, NotificationColor::Passed);
            }
            if (readMode == 1 || readMode == 2) word += character;
        }
    }
    numbers.push_back(word);
    _mesh.CommitData(&numbers);
}

void MeshFileStreamer::Write(const Mesh* mesh, FileExtension fe) {
    switch (fe)
    {
    case Txt:
        WriteTxt(mesh);
        break;
    case Bin:
        WriteBin(mesh);
        break;
    default:
        break;
    }
}

void MeshFileStreamer::WriteTxt(const Mesh* mesh) {
    std::ofstream fout(pathToText + "GeneratedPoints.txt");

    // Write array of points.
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of points. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total points ] " << mesh->Points.size() << std::endl;
    for (const auto &pnt : mesh->getPoints())
        fout << std::scientific << std::setprecision(15) << pnt.x << " " << pnt.y << " " << pnt.z << std::endl;
    fout.close();


    // Write array of areas.
    fout.open(pathToText + "GeneratedAreas.txt", std::ios::out);
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of areas. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total areas ] " << mesh->getAreasAsRibs().size() << std::endl;
    for (const auto& area : mesh->getAreasAsRibs())
        fout << area.subdomainNum_ << " " << area.refs_[0] << " " << area.refs_[1] << " "
            << area.refs_[2] << " " << area.refs_[3] << " " << area.refs_[4] << " "
            << area.refs_[5] << " " << area.refs_[6] << " " << area.refs_[7] << " "
            << area.refs_[8] << " " << area.refs_[9] << " " << area.refs_[10] << " "
            << area.refs_[11] << std::endl;
    fout.close();

    
    // Write array of ribs.
    fout.open(pathToText + "GeneratedRibs.txt", std::ios::out);
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of areas. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total ribs ] " << mesh->getRibsRefs().size() << std::endl;
    for (const auto& rib : mesh->getRibsRefs())
        fout << rib.p1 << " " << rib.p2 << std::endl;
    fout.close();


    // Write array of border ribs.
    fout.open(pathToText + "GeneratedBorderRibs.txt", std::ios::out);
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of border ribs. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total border ribs ] " << mesh->getBorderRibs().size() << std::endl;
    for (const auto& borderRib : mesh->getBorderRibs())
        fout << borderRib.type_ << " " << borderRib.formulaNum_ << " " << borderRib.ribRef_ << std::endl;
    fout.close();
}

void MeshFileStreamer::WriteBin(const Mesh* mesh) {
    std::ofstream fout(pathToBin + "GeneratedPoints.bin", std::ios::out | std::ios::binary);
    
    // Write array of points.
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of points. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total points ] " << mesh->Points.size() << std::endl;
    for (const auto& pnt : mesh->getPoints())
        fout << std::scientific << std::setprecision(15) << pnt.x << " " << pnt.y << " " << pnt.z << std::endl;
    fout.close();


    // Write array of areas.
    fout.open(pathToBin + "GeneratedAreas.bin", std::ios::out);
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of areas. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total areas ] " << mesh->getAreasAsRibs().size() << std::endl;
    for (const auto& area : mesh->getAreasAsRibs())
        fout << area.subdomainNum_ << " " << area.refs_[0] << " " << area.refs_[1] << " "
        << area.refs_[2] << " " << area.refs_[3] << " " << area.refs_[4] << " "
        << area.refs_[5] << " " << area.refs_[6] << " " << area.refs_[7] << " "
        << area.refs_[8] << " " << area.refs_[9] << " " << area.refs_[10] << " "
        << area.refs_[11] << std::endl;
    fout.close();


    // Write array of ribs.
    fout.open(pathToBin + "GeneratedRibs.bin", std::ios::out);
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of areas. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total ribs ] " << mesh->getRibsRefs().size() << std::endl;
    for (const auto& rib : mesh->getRibsRefs())
        fout << rib.p1 << " " << rib.p2 << std::endl;
    fout.close();


    // Write array of border ribs.
    fout.open(pathToBin + "GeneratedBorderRibs.bin", std::ios::out);
    if (!fout.is_open()) {
        Logger::ConsoleOutput("Error during file writing array of border ribs. File isn't open", NotificationColor::Alert);
        exit(-1);
    }
    fout << "[ Total border ribs ] " << mesh->getBorderRibs().size() << std::endl;
    for (const auto& borderRib : mesh->getBorderRibs())
        fout << borderRib.type_ << " " << borderRib.formulaNum_ << " " << borderRib.ribRef_ << std::endl;
    fout.close();
}
