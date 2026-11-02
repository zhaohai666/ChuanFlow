#pragma once

#include <iostream>
#include <sstream>
#include <string>

namespace cf::log
{

template<typename... Args>
std::string join(Args&&... args)
{
    std::ostringstream os;
    (os << ... << args);
    return os.str();
}

template<typename... Args>
void info(Args&&... args)
{
    std::cout << "[chflow] " << join(std::forward<Args>(args)...) << '\n';
}

template<typename... Args>
void warn(Args&&... args)
{
    std::cerr << "[chflow] warn: " << join(std::forward<Args>(args)...) << '\n';
}

template<typename... Args>
void error(Args&&... args)
{
    std::cerr << "[chflow] error: " << join(std::forward<Args>(args)...) << '\n';
}

} // namespace cf::log
