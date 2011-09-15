#include "interface/exception.hpp"

using namespace theta;

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
