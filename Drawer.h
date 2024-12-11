#pragma once

#include <Windows.h>
#include <string>

class Drawer {
private:
    static const char* scriptName;
    static const char* pythonScriptName;
    static const char* fileName;
public:
    Drawer() = delete;
    static void DrawMesh();
    static void DrawSolution(/*const FEM fem*/);
};

