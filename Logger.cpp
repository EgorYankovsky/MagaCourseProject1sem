#include "Logger.h"

void Logger::SetColor(NotificationColor clr) {
    switch (clr)
    {
    case NotificationColor::Passed:
        std::cout << "\033[" << 32 << "m";
        break;
    case NotificationColor::Warning:
        std::cout << "\033[" << 33 << "m";
        break;
    case NotificationColor::Alert:
        std::cout << "\033[" << 31 << "m";
        break;
    }
}

void Logger::ResetColor() { std::cout << "\033[0m"; }

void Logger::ConsoleOutput(std::string textMessage, NotificationColor clr) {
    Logger::SetColor(clr);
    std::cout << textMessage << std::endl;
    Logger::ResetColor();
}
