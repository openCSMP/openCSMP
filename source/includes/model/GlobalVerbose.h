// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GLOBAL_VERBOSE_H
#define GLOBAL_VERBOSE_H

#include "GenericSingleton.h"

namespace csmp {

 /**
 @brief Globally visible switch to turn console output from classes on/off.
 
 @author P. Lang
 @date Aug 2011
 
 To improve global status of the verbose variable. No explicit instantiation required any more.
 
 Usage:
 
 @code
 #include "GlobalVerbose.h"
 ...
 bool& global_verbose( GlobalVerbose::Instance().globalVerbose );
 @endcode
 
 */
class GlobalVerbose  : public GenericSingleton<GlobalVerbose>{

  friend class GenericSingleton<GlobalVerbose>;

  public:
    ~GlobalVerbose();

  public:
    bool globalVerbose;

  private:
    GlobalVerbose() : globalVerbose( true ) {}
 };
 


} // csmp

#endif // GLOBAL_VERBOSE_H
