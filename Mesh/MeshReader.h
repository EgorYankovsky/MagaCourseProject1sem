#pragma once

#include <string>
#include <fstream>
#include <sstream>

#include "Mesh.h"

class MeshReader
{
public:
    MeshReader() = delete;
    void Read(Mesh& mesh, const std::string path);
};

