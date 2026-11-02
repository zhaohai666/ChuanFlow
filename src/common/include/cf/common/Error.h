#pragma once

#include <stdexcept>
#include <string>

namespace cf
{

class Error : public std::runtime_error
{
public:
    explicit Error(const std::string& what) : std::runtime_error(what) {}
};

} // namespace cf

#define cfFatal(Msg) throw ::cf::Error(std::string(__FILE__) + ":" + std::to_string(__LINE__) + " " + (Msg))
