// Copyright (c) 2001 by Andrei Alexandrescu
// Copyright (c) 2001. Addison-Wesley.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STATIC_CHECK_INC_
#define STATIC_CHECK_INC_

namespace Loki
{

/**
@namespace Loki
@brief Andrei Alexandrescu
*/

////////////////////////////////////////////////////////////////////////////////
// Helper structure for the STATIC_CHECK macro
////////////////////////////////////////////////////////////////////////////////

    template<int> struct CompileTimeError;
    template<> struct CompileTimeError<true> {};
}

////////////////////////////////////////////////////////////////////////////////
// macro STATIC_CHECK
// Invocation: STATIC_CHECK(expr, id)
// where:
// expr is a compile-time integral or pointer expression
// id is a C++ identifier that does not need to be defined
// If expr is zero, id will appear in a compile-time error message.
////////////////////////////////////////////////////////////////////////////////

#define STATIC_CHECK(expr, msg) \
    { Loki::CompileTimeError<((expr) != 0)> ERROR_##msg; (void)ERROR_##msg; } 


////////////////////////////////////////////////////////////////////////////////
// Change log:
// March 20, 2001: add extra parens to STATIC_CHECK - it looked like a fun 
//     definition
// June 20, 2001: ported by Nick Thurn to gcc 2.95.3. Kudos, Nick!!!
////////////////////////////////////////////////////////////////////////////////

#endif // STATIC_CHECK_INC_
