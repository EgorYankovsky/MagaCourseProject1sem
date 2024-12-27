#include "MeshReader.h"

void MeshReader::Read(Mesh& mesh, const std::string path) {
	std::ifstream fin(path);
	std::stringstream buffer;
	if (fin.is_open()) {
		buffer << fin.rdbuf();
	}
	else {
		Logger::ConsoleOutput("File to read doesn't exists.", NotificationColor::Alert);
		Logger::FileOutput("File to read doesn't exists.", NotificationColor::Alert);
	}
	fin.close();
}
