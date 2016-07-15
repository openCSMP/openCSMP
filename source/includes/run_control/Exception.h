#ifndef CSMP_EXCEPTION_H
#define CSMP_EXCEPTION_H

#include <string>
#include <exception>
#include "CSMP_global_enumerations.h"

namespace csmp {

/// classification of screen output
enum CSMP_OUTPUT_LEVEL { VERBOSE=100, SILENT=0, IMPORTANT_OUTPUT=1, PROCESS_OUTPUT=2, INFO_OUTPUT=3 };

/// classification of screen output
enum CSMP_MESSAGE { INFO=0, EXCEPTION=-1, WARNING=-2, ERROR=-3, FATAL_ERROR=-4 };

/// convert enum to a printable string
std::string  parseMessage( CSMP_MESSAGE );

/// base class for CSMP-specific exception handling 
class Exception : std::exception {
  public:
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

    std::string   What() const;
    std::string   Originator() const;
    CSMP_MESSAGE  Message() const;
    void          Out() const;

  private:
    CSMP_MESSAGE csmp_exception_;
    std::string  originator_;
    std::string  message_;
    Exception();
 };

} // csmp

#endif



