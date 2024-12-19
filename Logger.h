#pragma once

#include <iostream>
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
};