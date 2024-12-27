#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

enum class NotificationColor {
    Passed,
    Warning,
    Alert
};

class Logger {
private:
    static void SetColor(NotificationColor clr);
    static void ResetColor();
public:
    Logger() = delete;
    static void ConsoleOutput(std::string textMessage, NotificationColor clr);
    static void FileOutput(std::string textMessage, NotificationColor clr);
};