// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MODELTIME_H
#define MODELTIME_H

#include "GenericSingleton.h"

namespace csmp {

/**

 @brief The only "global" variable in CSMP computations that is used to store simulated
 (as opposed to computation) time so that it can be accessed across classes and processors.
 
 @author P. Lang
 @date Aug 2011
 
 To improve global status of the verbose variable. No explicit instantiation required any more.
 
 Usage:
 
 @code
 #include "ModelTime.h"
 ...
 bool& model_time( ModelTime::Instance().modelTime );
 @endcode
 
 */
class ModelTime  : public GenericSingleton<ModelTime>{

  friend class GenericSingleton<ModelTime>;

  public:
    ~ModelTime();

  public:
    double modelTime;

  private:
    ModelTime() : modelTime( 0. ) {}
 };
 

} // csmp

#endif // MODELTIME_H
