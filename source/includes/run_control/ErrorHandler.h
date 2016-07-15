#ifndef CSMP_ERROR_HANDLER_H
#define CSMP_ERROR_HANDLER_H

#include "CSMP_number_types.h"
#include "GenericSingleton.h"
#include "Exception.h"
#include "BE_Time.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>

namespace csmp {

class ErrorHandler : public GenericSingleton<ErrorHandler> {

  friend class GenericSingleton<ErrorHandler>;

  public:    
    ~ErrorHandler();

    void   Verbose( bool verb );
    void   Verbose( size_t verb );
    size_t Verbose() const;

    void   Debug( bool debug );
    bool   Debug() const;

    bool operator> ( size_t verbose_level )  const;
    bool operator>=( size_t verbose_level ) const;
    bool operator==( size_t verbose_level ) const;
    bool operator!=( size_t verbose_level ) const;
    
    void WriteErrorsToFile( const char* err_file="ErrorHandler.log" );
    
    void notice( CSMP_MESSAGE err_type, 
                 const std::string& source, 
                 const std::string& message );
                 
    void notice( CSMP_MESSAGE err_type,
                 const std::string& source,
                 const std::string& message1,
                 const std::string& message2);

    void notice( CSMP_MESSAGE err_type,
                 const char* source,
                 const std::string& message );

    void notice( CSMP_MESSAGE err_type,
                 const char* source,
                 const std::string& message1,
                 const std::string& message2);

    void notice( CSMP_MESSAGE err_type,
                 const char* source,
                 const char* message );

    void notice( CSMP_MESSAGE err_type,
                 const char* source,
                 const char* message1,
                 const char* message2);

    void PrintInfos() const;
    void PrintWarnings() const;
    void PrintErrors() const;
    size_t MaxNumberOfErrors() const;
    void Out() const;

  private:
    ErrorHandler();

    void Print( const std::multimap<std::string,std::string>& ) const;

  private:

    BE_Time      timer_;
    size_t       verbose;
    bool         debug;
    long         errors;
    long         warnings;
    long         total;
    const int32  error_limit;

    //       originator & message      
    std::multimap<std::string,std::string>  minfos;
    std::multimap<std::string,std::string>  mwarnings;
    std::multimap<std::string,std::string>  merrors;
    //  time       message
    std::map<std::string,std::string>       error_sequence;
 }; 
 


/**
 
@class ErrorHandler ErrorHandler "main_library/ErrorHandler.h"

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1999

@section motivation Motivation 
 
The ErrorHandler implements a reporting scheme for CSMP designed to
take advantage of platform-specific message capabilities like error
windows or specific options that can be specified for core-dumps etc. 
Thus, the ErrorHandler also tries to support the typical Macintosh error
message boxes called Alert() and future platform-specific capabilities
are planned. The error handler further needs to document the progress
of a run such that it is documented for reference if the results are 
used for publication etc. For long time-dependent runs the aim is
to provide diagnostics of when and under what conditions an
error arose.  

As another important task, the ErrorHandler classifies (by enumeration),
the different messages that may be sent to him as INFOs, WARNINGs, ERRORs, and 
FATAL_ERRORs: 

@code
enum csmp_error{ INFO, EXCEPTION, WARNING, ERROR, FATAL_ERROR };
@endcode

If a FATAL_ERROR is reported, the ErrorHandler will terminate
the run. This role is given to the error handler with its future role as
exception handler in mind. Present FATAL_ERRORs will become recoverable
once exception handling is consistently implemented throughout CSMP.
 
 
@section design Design Intent

The ErrorHandler logs information strings to 'cout' and to the log
file 'ErrorHandler.log'. To be able to trace a message to its source, a string
that consists of the class and method name from which the message originated
must be supplied as first message argument. In case of an overloaded method
additional information is included in braces which identifies the called
method uniquely, for example:

@code
ErrorHandler::Notice(char*, char*):
@endcode

@section messages Messages 

Messages which are written to file, are automatically preceded by a time
string which tells the user when the message originated. Since messages
may also contain numerical data in the description of the error message
the ErrorHandler interface is overloaded to accept std::string arguments
(std::string has overloaded operators that allow the method writer to convert 
numbers to strings). The functionality of the ErrorHandler is planned to be 
extended to handle top level (C++) exceptions. For this purpose the standard
exception base classes exception, logic_error, and runtime_error will be
subclassed for CSMP-specific error types.

A goal of the exception handling is to store the transient state
of a run before it gets aborted such that it can be carried on from this
state at a later time, when the error condition has been properly
diagnosed.
 
 
@section applicability Applicability

The ErrorHandler is used inside high-level CSMP objects like the
Model etc.. Only those which carry enough information that it makes sense
to collect it into user-readable diagnostics on a CSP
run. Atomic generic classes report errors directly to the
'cout' stream to prevent a loss in their generality.
 
 
@section structure Structure

ErrorHandler is a base-level class using STL maps to store the
different types of error messages.
 
 
@section participants Participants

A singleton instance(static) is created upon loading the csmp core library.
 
 
@section collaborations Collaborations

ErrorHandler uses the functionality of STL maps and the std::string class.
The Standard_IO_Handler is used in dialogs with the CSP user and the
object BE_Time is contained to log the time at which errors occur.
 
 
@section consequences Consequences

The consequent usage of the ErrorHandler turns it into an 'Observer'
object the functionality of which may be widely extended in future
implementations. For the time being, the user can use the '*.log' file
and the screen output to obtain a detailed progress report for his/her
CSMP run. Interaction with the run, however, is still very limited.
 
 
@section implementation Implementation

The ErrorHandler stores the message strings in 3 maps corresponding
to the message types: INFO, WARNING, ERROR. Since only a single FATAL_ERROR
can occur before the run is terminated, it is lumped into the ERROR
map. In the log file, all messages are logged sequentially. 

The error log file is written by the destructor of the error handler. Because
the destructor of a global object is called last, this assures that all
messages that arise in the destruction process are logged as well. 
 
 
@section application Application Examples

The ErrorHandler is implemented as a Singleton pattern (GoF).

A reference to the ErrorHanlder instance is acquired as such

@code
#include "ErrorHandler.h"
...
ErrorHandler& csmp_error( ErrorHandler::Instance() );
@endcode

after which one can call the handler

@code
csmp_error.notice( ERROR, "IterativeAlgorithm<dim>::Iterations", "Solution did not converge! " );
@endcode

alternatively, if i.e. just a single call is to be made

@code
ErrorHandler::Instance().notice( ERROR, "IterativeAlgorithm<dim>::Iterations", "Solution did not converge! " );
@endcode

@todo SKM - catch some standard exceptions

@todo SKM - use typename capability to find out from which inherited object the message came

@todo SKM - allow user to intervene and fix some problems

*/


/**
    Returns the maximum number of errors that the ErrorHandler permits before 
    issueing a fatal (program terminating) error.
*/
inline size_t ErrorHandler::MaxNumberOfErrors() const
 { return static_cast<size_t>(error_limit); }



} // csmp
 
 #endif






 
 
