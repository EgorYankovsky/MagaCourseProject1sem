#pragma once

#include <fstream>

#include "Mesh.h"
#include "../Logger/Logger.h"

enum FileExtension {
    Txt,
    Bin
};

class MeshFileStreamer {
private:
    static std::string const pathToText;
    static std::string const pathToBin;
    static void WriteTxt(const Mesh* mesh, std::string path = pathToText);
    static void WriteBin(const Mesh* mesh, std::string path = pathToBin);
public:
    MeshFileStreamer() = delete;
    static void Read(Mesh& _mesh, std::string path);
    static void Write(const Mesh* mesh, FileExtension fe);
};

