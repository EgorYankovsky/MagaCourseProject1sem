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
    switch (clr) {
    case NotificationColor::Passed:
#ifdef DEBUG
        std::cout << textMessage << std::endl;
#endif
        break;
    case NotificationColor::Alert:
        std::cout << '\a' << std::endl;
    case NotificationColor::Warning:
        auto stringSize = textMessage.size();
        for (size_t i(0); i < stringSize + 4; ++i) std::cout << "*";
        std::cout << std::endl;
        std::cout << "* " << textMessage << " *" << std::endl;
        for (size_t i(0); i < stringSize + 4; ++i) std::cout << "*";
        std::cout << std::endl;
        break;
    }
    Logger::ResetColor();
}

void Logger::FileOutput(std::string textMessage, NotificationColor clr) {
    std::ofstream fout("log\\log.txt");
    time_t currentTime = time(NULL);
    switch (clr) {
    case NotificationColor::Passed:
#ifdef DEBUG
        fout << "Passed task." << ctime(&currentTime) << std::endl;
        fout << textMessage << std::endl;
#endif
        break;
    case NotificationColor::Warning:
        fout << "Warning! " << ctime(&currentTime);
        fout << textMessage << std::endl;
        break;
    case NotificationColor::Alert:
        fout << "Alert!!! " << ctime(&currentTime);
        fout << textMessage << std::endl;
        break;
    }
    fout.close();
}
