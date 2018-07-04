#ifndef CSMP_ERROR_HANDLER_H
#define CSMP_ERROR_HANDLER_H

#include "CSMP_number_types.h"
#include "GenericSingleton.h"
#include "Exception.h"
#include "BE_Time.h"

namespace csmp {

/// reporting level set for the error handler; default=verbose
enum CSMP_OUTPUT_LEVEL { VERBOSE=100, SILENT=0, IMPORTANT_OUTPUT=1, PROCESS_OUTPUT=2, INFO_OUTPUT=3 };

/**
 
@brief The singleton ErrorHandler implements an issue reporting scheme to support 
platform-specific message capabilities like error
windows or specific options that can be specified for core-dumps etc. 
The error handler further documents the progression of issues that occur during
a run so that warnings etc. leading up to an error can be examined in the 
runtime.log file. Thus, the aim is
to provide diagnostics of when and under what conditions an
error arose.  

@author S.K. Matthai
@author Stephen G. Roberts
@date 1999

@section messages Messages 

The ErrorHandler classifies (by enumeration),
the different messages that are sent to it as INFOs, WARNINGs, ERRORs, and
FATAL_ERRORs: 

@code
enum csmp_error{ INFO, EXCEPTION, WARNING, ERROR, FATAL_ERROR };
@endcode

If a FATAL_ERROR is reported, the ErrorHandler will terminate
the run. This role is given to the error handler with its future role as
exception handler in mind. Present FATAL_ERRORs will become recoverable
once exception handling is consistently implemented throughout CSMP.
 
@section design Design Intent

The ErrorHandler logs information strings to 'cout' / 'cerr' and to the log
file 'ErrorHandler.log'. To be able to trace a message to its source, a string
that consists of the class and method name from which the message originated
must be supplied as first message argument. In case of an overloaded method
additional information is included in braces which identifies the called
method uniquely, for example:

@code
ErrorHandler::Notice(char*, char*):
@endcode

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
Model etc.. Only those which carry enough information so that it makes sense
to collect user-readable diagnostics from a CSMP run are logged. 
Atomic generic classes report errors directly to the
'cerr' stream.
 

@section structure Structure

ErrorHandler is a base-level class using STL maps to store the
different types of error messages.
 
 
@section participants Participants

A singleton instance(static) is created upon loading the csmp core library.
 
 
@section collaborations Collaborations

ErrorHandler uses the functionality of STL maps and the std::string class.
The Standard_IO_Handler is used in dialogs with a console user.
The object BE_Time is contained to log the time at which errors occur.
 
 
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
...
csmp_error.notice( ERROR, "MyClass::ProblematicMethod()", "what actually went wrong" );

@endcode

after which one can call the handler

@code
csmp_error.notice( ERROR, "IterativeAlgorithm<dim>::Iterations", "Solution did not converge! " );
@endcode

alternatively, if i.e. just a single call is to be made

@code
ErrorHandler::Instance().notice( ERROR, "IterativeAlgorithm<dim>::Iterations", "Solution did not converge! " );
@endcode

@todo SKM - resolve how ErrorHandler should interact with exception handling

@todo SKM - use typename capability to find out from which inherited object the message came

@todo SKM - allow user to intervene and fix some problems

*/
class ErrorHandler : public GenericSingleton<ErrorHandler> {

  friend class GenericSingleton<ErrorHandler>;

  public:
    /// writes a file called 'Runtime.log' with all the current messages when the program terminates
    ~ErrorHandler();

    /// sets output reporting level to verbose (reporting of all events to cout and cerr)
    void   Verbose( bool VERBOSE );
  
    /// reports whether the reporting level is set to VERBOSE
    bool   Verbose() const;
  
    /// sets output to corresponding CSMP_OUTPUT_LEVEL
    void   ReportingLevel( CSMP_OUTPUT_LEVEL verb );
  
    /// returns the current setting of the CSMP_OUTPUT_LEVEL
    CSMP_OUTPUT_LEVEL ReportingLevel() const;

    /// use as @code ErrorHandler& csmp_error( ErrorHandler::Instance() ); csmp_error.notice( WARNING, "class::method where this comes from", messsage ); @endcode
    void notice( CSMP_MESSAGE err_type,
                 const std::string& source, 
                 const std::string& message );
  
    /// 4 argument version of notice
    void notice( CSMP_MESSAGE err_type,
                 const std::string& source,
                 const std::string& message1,
                 const std::string& message2 );

    void   PrintInfos() const;
    void   PrintWarnings() const;
    void   PrintErrors() const;
  
    /// returns how many errors are allowed before program termination is triggered
    size_t MaximumNumberOfErrors() const;
  
    /// writes the content of message containers stored by the handler to file
    void WriteErrorsToFile( const char* err_file="ErrorHandler.log" );

    /// writes contents of contained maps to file
    void Out() const;

  private:
    /// sets default to VERBOSE, error limit to 10'000, and counters to zero; private because ErrorHandler is a singleton
    ErrorHandler();

  private:

    BE_Time            timer_;        ///< to measure time level (replace by C++11 time facilities
    CSMP_OUTPUT_LEVEL  verbose_;      ///< CSMP_OUTPUT_LEVEL stored as a short enumeration (see csmp::Exception)
    long               errors_;       ///< errors incurred where something actually went wrong as opposed to a problem with- or without consequences
    long               warnings_;     ///< problems incurred
    long               total_;        ///< total issues logged
    const size_t       error_limit_;  ///< number of errors that are allowed to occur before program terminates

    //       originator & message      
    std::multimap<std::string,std::string>  minfos_;          ///< log of info messages
    std::multimap<std::string,std::string>  mwarnings_;       ///< log of warnings
    std::multimap<std::string,std::string>  merrors_;         ///< log of errors
    //  time       message
    std::map<std::string,std::string>       error_sequence_;  ///< chronological sequence of messages
 }; 
 
} // csmp
 
 #endif






 
 
