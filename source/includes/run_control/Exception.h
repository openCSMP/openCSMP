#ifndef CSMP_EXCEPTION_H
#define CSMP_EXCEPTION_H

#include <cstdint>
#include <string>
#include <exception>

namespace csmp {

/// classification of exception; warning (possible problem), error (actual problem; fatal (need to terminate)
enum CSMP_MESSAGE : std::int8_t { INFO=0, EXCEPTION=-1, WARNING=-2, ERROR=-3, FATAL_ERROR=-4 };

/**
@brief CSMP exception class.

Thrown by CSMP methods to signal error conditions. Inherits publicly
from std::exception so it can be caught as std::exception& as well as
csmp::Exception&.

In debug builds (NDEBUG not defined), the exception prints its content
to stdout when thrown. In release builds no output is produced until
the exception is caught and handled.

@section usage Usage

@code
throw csmp::Exception( ERROR,
                       "MyClass::MyMethod",
                       "something went wrong" );

throw csmp::Exception( ERROR,
                       "MyClass::MyMethod",
                       variable_name,
                       "something went wrong with this variable" );
@endcode
*/
class Exception : public std::exception
{
public:
    /**
    Constructs an exception with a severity, originator, and message.
    In debug builds, prints the exception content to stdout immediately.
    */
    Exception( CSMP_MESSAGE     severity,
               const std::string& originator,
               const std::string& message );

    /**
    Constructs an exception with a severity, originator, parameter
    name, and message. The parameter name is prepended to the message.
    In debug builds, prints the exception content to stdout immediately.
    */
    Exception( CSMP_MESSAGE     severity,
               const std::string& originator,
               const std::string& parameter,
               const std::string& message );

    Exception( const Exception& ) = default;
    Exception& operator=( const Exception& ) = default;

    virtual ~Exception() noexcept = default;

    /**
    Returns the full message string (parameter + message if applicable).
    Overrides std::exception::what().
    */
    virtual const char* what() const noexcept override;

    /** Returns the full message as a std::string. */
    std::string   What()       const noexcept;

    /** Returns the originator string (class::method). */
    std::string   Originator() const noexcept;

    /** Returns the severity level. */
    CSMP_MESSAGE  Message()    const noexcept;

    /** Prints all exception data to stdout. */
    void          Out()        const noexcept;

private:
    CSMP_MESSAGE csmp_exception_;
    std::string  originator_;
    std::string  message_;

    /** Called by constructors in debug builds to print and optionally
        pause. Separated from the constructor so the constructor body
        remains clean and the debug behaviour can be controlled
        independently.
        
        @attention Define CSMP_MANUAL_TESTING, if you want ERROR to stop programme at exception construction
    */
    void DebugOutput() const;

    Exception() = delete;
};

/** Converts a CSMP_MESSAGE severity to a human-readable string. */
std::string parseMessage( CSMP_MESSAGE msg ) noexcept;

} // namespace csmp

#endif // CSMP_EXCEPTION_H

