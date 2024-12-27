#include "Drawer.h"

const char* Drawer::scriptName = "python";
const char* Drawer::pythonScriptName = "pythonScripts\\drawer.py";
const char* Drawer::pointsFileName = "Data\\Generated\\generatedPoints.txt";
const char* Drawer::ribsFileName = "Data\\Generated\\generatedRibs.txt";

void Drawer::DrawMesh() {
    char commandToRun[100];
    strcpy_s(commandToRun, scriptName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, pythonScriptName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, pointsFileName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, ribsFileName);
    system(commandToRun);
}
