// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MODEL_TIME_TO_INTEGER_H
#define MODEL_TIME_TO_INTEGER_H

#include "CSMP_definitions.h"

namespace csmp {

template<typename fT>
class ModelTimeToInteger {
  public:
    ModelTimeToInteger();
    ~ModelTimeToInteger();
    
   long ModelTimeInMinutes( fT t ); 
   long ModelTimeInHours( fT t ); 
   long ModelTimeInDays( fT t ); 
   long ModelTimeInWeeks( fT t ); 
   long ModelTimeInMonths( fT t ); 
   long ModelTimeInYears( fT t ); 
   long ModelTimeInKiloYears( fT t ); 
  
  private:
    const fT minute;
    const fT hour;
    const fT day;
    const fT week;
    const fT month;
    const fT year;
    const fT kiloyear;
 };

} // end namespace csmp












#endif
