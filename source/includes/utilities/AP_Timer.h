#ifndef CSMP_TIMER_H
#define CSMP_TIMER_H

#include "CSMP_definitions.h"

/**
  Timer to track program running time (wall clock time).
  Could use std::clock() and CLOCK_PER_SEC w/<ctime>
  and not have to use and ifdefs in the code, may sacrifice
  some precision. see www.boost.org
  Adapted from Boost A Paluszny 2006 
*/

namespace csmp {

class AP_Timer {
public:
  /// Constructor.
  AP_Timer ();
  
  ~AP_Timer() {}

  /// Reset the timer to zero and activate at the same time.
  AP_Timer & Reset();

  /// Stop the timer (and freeze the elapsed time). 
  AP_Timer & Stop();

  /// Continiue the timer from where Stop() had frozen it. 
  AP_Timer & Continue();

  /// Return the currently elapsed time. 
  double64 ElapsedTime();

  /// Set the currently elapsed time
  void SetElapsedTime(const double64 & elapsedTime);

  /// Returns the elapsed time as a string.
  std::string AsString();

  std::string AverageTimeAsString(size_t uSteps);

  static double64 GetSecondsPerClockTick();

protected:
  
  /// read the system time and then update m_elapsedTime accordingly
  void UpdateTimes();


protected:
  /// true if the timer is running
  bool m_bTimerIsRunning;

  /// the absolute time the timer was created or restarted
  double64 m_fStartTime;

  /// the relative time since the timer was created or restarted
  double64 m_fElapsedTime;

  /// precision of the system clock
  static const double64 m_secondsPerTick;

  /// output to stream
  friend std::ostream & operator << (std::ostream & os, AP_Timer & timer);

}; // Timer

extern std::ostream & operator << (std::ostream & os, AP_Timer & timer);

} //end namespace csmp

#endif  
