#ifndef CSMP_EXCEPTION_H
#define CSMP_EXCEPTION_H

#include <string>
#include <exception>
#include "CSMP_global_enumerations.h"

namespace csmp {

/// classification of exception; warning (possible problem), error (actual problem; fatal (need to terminate)
enum CSMP_MESSAGE { CSMP_INFO = 0, CSMP_EXCEPTION = -1, CSMP_WARNING = -2, CSMP_ERROR = -3, CSMP_FATAL_ERROR = -4 };


// Support for old code
#if defined(INFO)
#warning "INFO has already been defined somewhere in an included dependency"
#undef INFO
#endif
#define INFO CSMP_INFO

#if defined(EXCEPTION)
#warning "EXCEPTION has already been defined somewhere in an included dependency"
#undef EXCEPTION
#endif
#define EXCEPTION CSMP_EXCEPTION

#if defined(WARNING)
#warning "WARNING has already been defined somewhere in an included dependency"
#undef WARNING
#endif
#define WARNING CSMP_WARNING

#if defined(ERROR)
#warning "ERROR has already been defined somewhere in an included dependency"
#undef ERROR
#endif
#define ERROR CSMP_ERROR

#if defined(FATAL_ERROR)
#warning "FATAL_ERROR has already been defined somewhere in an included dependency"
#undef FATAL_ERROR
#endif
#define FATAL_ERROR CSMP_FATAL_ERROR

/// convert enum to a printable string
std::string  parseMessage( CSMP_MESSAGE );

/** @brief base class for CSMP-specific exception handling using an exception qualifier (CSMP_MESSAGE),
    an originator (the class and method from where the exception comes from, and a string
    that contains the message that describes what happened locally.
    
     @note no need to use a csmp::Exception if the case is well described by one of the 
     standard exceptions.
     
      @author Stephan Matthai
      @date 1999
*/
class Exception : std::exception {
  public:
    /// use as: throw csmp::Exception( CSMP_WARNING, "class::method where this comes from", messsage );
    Exception( CSMP_MESSAGE, 
               const std::string& originator,
               const std::string& message );

    Exception( CSMP_MESSAGE, 
               const std::string& originator,
               const std::string& parameter,
               const std::string& message );

    virtual ~Exception() throw();
                          
    Exception( const Exception& );
    Exception& operator=( const Exception& );
  
    virtual const char* what() const throw();

    /// returns what() for the underlying standard exception
    std::string   What() const;
  
    /// the class and method (or function) from which the exception was thrown
    std::string   Originator() const;
  
    /// what went wrong in as much as was known at the executation of the program from where the exception was thrown
    CSMP_MESSAGE  Message() const;
  
    /// prints all data assocatiated with the exception (following the convention for all CSMP classes)
    void          Out(std::ostream& os) const;

  private:
    CSMP_MESSAGE csmp_exception_;
    std::string  originator_;
    std::string  message_;
    Exception();
 };

} // csmp

#endif



