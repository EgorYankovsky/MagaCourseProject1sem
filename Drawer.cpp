#include "Drawer.h"

#include <iostream>

const char* Drawer::scriptName = "python";
const char* Drawer::pythonScriptName = "pythonScripts\\drawer.py";
const char* Drawer::fileName = "Data\\generatedPoints.txt";


void Drawer::DrawMesh() {
    char commandToRun[100];
    strcpy_s(commandToRun, scriptName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, pythonScriptName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, fileName);
    system(commandToRun);
}
