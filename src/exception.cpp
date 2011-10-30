#include "interface/exception.hpp"

#include <sstream>
#include <typeinfo>

#include <stdexcept>


using namespace theta;

const char* Exception::what() const throw(){
    std::stringstream ss;
    ss << typeid(*this).name() << ": " << message;
    whatstring = ss.str();
    return whatstring.c_str();
}

Exception::Exception(const std::string & m): runtime_error(m), message(m){}
ConfigurationException::ConfigurationException(const std::string & msg): Exception(msg){}

void fail_assert(const char * filename, int lineno, const char * expression){
    std::stringstream ss;
    ss << "Assertion '" << expression << "' failed in " << filename << ":" << lineno;
    throw std::logic_error(ss.str());
}

DatabaseException::DatabaseException(const std::string & s): Exception(s){}

MinimizationException::MinimizationException(const std::string & s): Exception(s){}

ExitException::ExitException(const std::string & message_): message(message_){}
