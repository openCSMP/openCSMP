// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

// Exception.cpp

#include "Exception.h"
#include <iostream>

using namespace std;

namespace csmp {

Exception::Exception( CSMP_MESSAGE     severity,
                      const string&    orig,
                      const string&    msg )
    : csmp_exception_( severity ),
      originator_( orig ),
      message_( msg )
{
#ifndef NDEBUG
    DebugOutput();
#endif
}

Exception::Exception( CSMP_MESSAGE     severity,
                      const string&    orig,
                      const string&    param,
                      const string&    msg )
    : csmp_exception_( severity ),
      originator_( orig ),
      message_( param + "  " + msg )
{
#ifndef NDEBUG
    DebugOutput();
#endif
}

const char* Exception::what() const noexcept
{
    return message_.c_str();
}

string Exception::What() const noexcept
{
    return message_;
}

string Exception::Originator() const noexcept
{
    return originator_;
}

CSMP_MESSAGE Exception::Message() const noexcept
{
    return csmp_exception_;
}

void Exception::Out() const noexcept
{
    cout << "\n" << parseMessage( csmp_exception_ )
         << ": " << originator_ << "\n"
         << message_ << "\n";
}

void Exception::DebugOutput() const
{
    Out();
    // In debug builds, pause only for ERROR and FATAL_ERROR to avoid
    // interrupting normal WARNING flow in automated tests.
    if ( csmp_exception_ == ERROR || csmp_exception_ == FATAL_ERROR )
    {
#ifdef CSMP_MANUAL_TESTING
        cout << "\nHit return to continue." << endl;
        getchar();
#endif
    }
}

string parseMessage( CSMP_MESSAGE msg ) noexcept
{
    switch ( msg )
    {
        case INFO:        return "INFO";
        case WARNING:     return "WARNING";
        case ERROR:       return "ERROR";
        case FATAL_ERROR: return "FATAL_ERROR";
        case EXCEPTION:
        default:          return "EXCEPTION";
    }
}

} // namespace csmp

