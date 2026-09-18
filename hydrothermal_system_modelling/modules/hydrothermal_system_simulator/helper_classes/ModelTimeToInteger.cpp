// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ModelTimeToInteger.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

template<typename fT>
ModelTimeToInteger<fT>::ModelTimeToInteger()
 : minute(60.00),   		   	// minute in sec
   hour(3600.0),            // hour in sec
   day(86400.0),            // day in sec
   week(604800.0),          // week (7 days) in sec
   month(2419200.0),        // month (4 weeks) in sec
   year(31536000.0),        // year in sec
   kiloyear(31536000000.0)  // 1000 years in sec
 {
 }


template<typename fT>
ModelTimeToInteger<fT>::~ModelTimeToInteger()
  {
  }

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInMinutes( fT t )     
{ return static_cast<long>(rint( t / minute )); } 

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInHours( fT t )     
{ return static_cast<long>(rint( t / hour )); } 

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInDays( fT t )      
{ return static_cast<long>(rint( t / day )); }

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInWeeks( fT t )     
{ return static_cast<long>(rint( t / week )); }

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInMonths( fT t )     
{ return static_cast<long>(rint( t / month )); }

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInYears( fT t )     
{ return static_cast<long>(rint( t / year )); }

template<typename fT>
long ModelTimeToInteger<fT>::ModelTimeInKiloYears( fT t ) 
{ return static_cast<long>(rint( t / kiloyear )); } 


template class ModelTimeToInteger<double>;

} // end namespace csmp
