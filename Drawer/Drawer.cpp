#include "Drawer.h"

const char* Drawer::scriptName = "python";
const char* Drawer::pythonScriptName = "pythonScripts\\drawer.py";
const char* Drawer::pointsFileName = "Data\\Generated\\Text\\GeneratedPoints.txt";
const char* Drawer::ribsFileName = "Data\\Generated\\Text\\GeneratedRibs.txt";

void Drawer::DrawMesh() {
    char commandToRun[_MAX_PATH];
    strcpy_s(commandToRun, scriptName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, pythonScriptName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, pointsFileName);
    strcat_s(commandToRun, " ");
    strcat_s(commandToRun, ribsFileName);
    system(commandToRun);
}
