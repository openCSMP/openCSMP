//
//  license.h
//
//  for the compilation of the binary library to be shipped with a CSMP application
//
//  Created by Stephan Matthai on 8/16/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_LICENSE_H
#define CSMP_LICENSE_H


namespace csmp {

// in days
enum LICENSE_PERIOD { UNLIMITED=0, TRIAL=31, YEAR=365, THREE_YEAR=1095 };

std::string parse( LICENSE_PERIOD );

// stores start date and secs remaining  
struct Expiration {
   Expiration( LICENSE_PERIOD p, int startYear, int startMonth, int startDay );
   LICENSE_PERIOD  duration;
   std::tm         start_time;
   std::time_t     secs_remaining; // remaining seconds of licensing period
};



// ----------------------------------------------------
// here the license details need to be specified
// ----------------------------------------------------
Expiration  aus(YEAR,2018,2,28);



// enum to string
inline std::string parse( LICENSE_PERIOD p )
 {
    if ( p == UNLIMITED )  return "UNLIMITED";
    if ( p == TRIAL )      return "TRIAL";
    if ( p == YEAR )       return "ONE-YEAR";
    if ( p == THREE_YEAR ) return "THREE-YEAR";
    return "'invalid'";
 }



// initialisation
inline Expiration::Expiration( LICENSE_PERIOD p, int startYear, int startMonth, int startDay )
  : duration(p)
 {
      // establishing a time-record for the start date of the license
      // ------------------------------------------------------------
	    start_time.tm_sec=0.;		        /* seconds after the minute [0-60] */
	    start_time.tm_min=0.;		        /* minutes after the hour [0-59] */
	    start_time.tm_hour=0.;	        /* hours since midnight [0-23]  = midnight */
	    start_time.tm_mday=startDay;	  /* day of the month [1-31] */
	    start_time.tm_mon=startMonth;		/* months since January [0-11] */
	    start_time.tm_year=startYear-1900;	/* years since 1900 */
     
      // establishing how many seconds still remain
      std::time_t now = std::time(NULL);

	  #ifdef _MSC_VER
		#if _MSC_VER < 1800   // detecting versions older than 2013
			secs_remaining = static_cast<long>(duration) * 86400L - static_cast<long>(ceil(std::difftime( now, std::mktime(&start_time))) );
        #else
	        secs_remaining = static_cast<long>(duration) * 86400L - lrint(std::difftime( now, std::mktime(&start_time)) );
        #endif
	  #else
			secs_remaining = static_cast<long>(duration) * 86400L - lrint(std::difftime( now, std::mktime(&start_time)) );
      #endif
//      std::cout <<"\nExpiration: remaining time: "<< secs_remaining <<"\n";
 }

} // end csmp

#endif
