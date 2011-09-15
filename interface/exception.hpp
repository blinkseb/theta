#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <sstream>
#include <typeinfo>


namespace theta {

/** \brief Base class for all Exceptions used in %theta
 */
class Exception: virtual public std::exception{
public:
    /// The human-readable, english message 
    std::string message;
    
    /// Constructor taking a message intended for the user which will be written to Exception::message
    Exception(const std::string & msg);
    
    /// Declare destructor virtual
    virtual ~Exception() throw(){}
    
    /// override std::exception::what to print out the message
    virtual const char* what() const throw(){
         std::stringstream ss;
         ss << typeid(*this).name() << ": " << message;
         whatstring = ss.str();
         return whatstring.c_str();
    }
protected:
    /// buffer for the const char * returned by what()
    mutable std::string whatstring;
};

/** \brief Thrown during configuration file processing.
 */
class ConfigurationException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    ConfigurationException(const std::string & msg);
};

/** \brief Thrown in case of a request for a non-existing member in a conatiner-like structure is made. 
 */
class NotFoundException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
   NotFoundException(const std::string & message);
};

/** \brief Thrown in case of invalid mathematical constructs like domain errors
 */
class MathException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    MathException(const std::string &);
};

/// \brief Thrown in case of database errors.
class DatabaseException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    DatabaseException(const std::string & s): Exception(s){}
};

/** \brief Thrown in case of minimization errors.
 *
 * \sa theta::Minimizer
 */
class MinimizationException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    MinimizationException(const std::string & s): Exception(s){}
};

/** \brief Exception class to indicate an internal error
 * 
 * An exception of this kind should usually not be caught: it indicates a serious error
 * which prevents further execution, such as violation of pre / postconditions, i.e., logic errors
 * which could be prevented by proper programming.
 * 
 * In order to prevent catching FatalException in a catch(Exception &) statement, FatalException is not part
 * of the usual exception hierarchy of theta.
 */
class FatalException{
public:
    std::string message;
       
    /** \brief Construct from a "usual" Exception
     * 
     * The error message displayed will be taken from ex. 
     */
    explicit FatalException(const Exception & ex);
    
    /// Construct directly from an error message
    explicit FatalException(const std::string & message);
};

/// \brief General exception to indicate that arguments passed to a function are invalid (=do not correspond to documentation)
class InvalidArgumentException: public FatalException{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    InvalidArgumentException(const std::string &);
};


/** \brief Thrown in case an immediate exit of the main program is requested. 
 *
 * As this exception is caught at a different place, it is not part of the usual Exception
 * hierarchy.
 */
class ExitException {
public:
   std::string message;

   /// Constructor taking a message intended for the user which will be written to Exception::message
   ExitException(const std::string & message_): message(message_){}
};

}

#endif
