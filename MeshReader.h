#pragma once

#include <string>
#include "Mesh.h"

class MeshReader
{
public:
    MeshReader() = delete;
    void Read(std::string path);
};

