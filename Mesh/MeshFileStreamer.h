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
    static std::string pathToText;
    static std::string pathToBin;
    static void WriteTxt(const Mesh* mesh);
    static void WriteBin(const Mesh* mesh);
public:
    MeshFileStreamer() = delete;
    static void Read(Mesh& _mesh, std::string path);
    static void Write(const Mesh* mesh, FileExtension fe);
};

