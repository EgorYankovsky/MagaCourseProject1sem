#include "MeshFileStreamer.h"

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
            Logger::ConsoleOutput(mes, NotificationColor::Warning);
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
                Logger::ConsoleOutput("Failure during reading", NotificationColor::Alert);
                exit(-1);
            }
            if (readMode == 0) {
                mes = "Unexpected character: " + character;
                Logger::ConsoleOutput(mes, NotificationColor::Passed);
            }
            if (readMode == 1 || readMode == 2) word += character;
        }
    }
}
