#ifndef BRUCE_ECKEL_TIMER_UTILITY_H
#define BRUCE_ECKEL_TIMER_UTILITY_H

#include <ctime>
#include "CSMP_definitions.h"

/** Bruce Eckel's Time class from "Thinking in C++" p. 314
*/
class BE_Time {
  public:
    /// constructor records current system time
    BE_Time() { mark(); }
    ~BE_Time() {}
    void mark() {
      lflag = aflag = 0;
      std::time(&time);
    }
    /// time as a day/month/year string
    const char* ascii() {
      updateAscii();
      return Ascii;
    }   
    /// difference in seconds
    double delta( const BE_Time& dt ) const {
      return std::difftime( time, dt.time );
    }
    int DaylightSavings() {
      updateLocal();
      return local.tm_isdst;
    }
    /// days since January 1
    int DayOfYear() {
      updateLocal();
      return local.tm_yday;
    }
    /// days since Sunday
    int DayOfWeek() {
      updateLocal();
      return local.tm_wday;
    }
    /// years since 1/1/1900
    int Since1900() {
      updateLocal();
      return local.tm_year;
    }
    /// since January
    int Month() {
      updateLocal();
      return local.tm_mon;
    }    
    int DayOfMonth() {
      updateLocal();
      return local.tm_mday;
    }
    /// since midnight (24 hour clock)
    int Hour() {
      updateLocal();
      return local.tm_hour;
    }
    int Minute() {
      updateLocal();
      return local.tm_min;
    }
    int Second() {
      updateLocal();
      return local.tm_sec;
    }
    /// gets ascii hours,minutes, seconds and tenth's of seconds out of the time-string
    void Time( long& yr, long& hr, long& mi, long& sec );
    void Time( long& year, long& month, long& day );

private:
    std::time_t    time;
    std::tm        local;
    char           Ascii[26];
    unsigned char  lflag, aflag;
  
    /// copy time into local time buffer
    void updateLocal() {
      if (!lflag) {
          local = *std::localtime(&time);
          lflag++;
        }
    }
    /// stores asctime results in local variable
    void updateAscii() {
      if (!aflag) {
          updateLocal();
          std::strcpy( Ascii, asctime(&local) );
          aflag++;
       }
    }
 };

#endif
