#include "interface/exception.hpp"

#include <sstream>
#include <typeinfo>


using namespace theta;

const char* Exception::what() const throw(){
    std::stringstream ss;
    ss << typeid(*this).name() << ": " << message;
    whatstring = ss.str();
    return whatstring.c_str();
}

InvalidArgumentException::InvalidArgumentException(const std::string & m) : FatalException(m) {}
Exception::Exception(const std::string & m):message(m){}
ConfigurationException::ConfigurationException(const std::string & msg): Exception(msg){}
NotFoundException::NotFoundException(const std::string & msg): Exception(msg){}
MathException::MathException(const std::string & m): Exception(m){}

FatalException::FatalException(const Exception & ex){
    message = ex.what();
}

FatalException::FatalException(const std::string & message_): message(message_){
}

void fail_assert(const char * filename, int lineno, const char * expression){
    std::stringstream ss;
    ss << "Assertion '" << expression << "' failed in " << filename << ":" << lineno;
    throw FatalException(ss.str());
}

