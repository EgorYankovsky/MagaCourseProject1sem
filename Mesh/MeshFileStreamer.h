#pragma once

#include <fstream>

#include "Mesh.h"
#include "../Logger/Logger.h"

class MeshFileStreamer
{
public:
    MeshFileStreamer() = delete;
    static void Read(Mesh& _mesh, std::string path);
    static void WriteTxt(const Mesh* mesh, std::string path);
    static void WriteBin(const Mesh* mesh, std::string path);
};

