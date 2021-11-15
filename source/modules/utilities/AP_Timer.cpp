#include "AP_Timer.h"
#ifdef _WIN32
  #include "Windows.h"
#elif __APPLE__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#else
  #include "sys/time.h"
#endif

using namespace std;

namespace csmp {
#ifndef _WIN32
   typedef unsigned long long LARGE_INTEGER;
#endif

const double AP_Timer::m_secondsPerTick(AP_Timer::GetSecondsPerClockTick());

//Constructor.
AP_Timer::AP_Timer ()
{
  Reset();
} // constructor


//Get processor spped from windows registry and set fSecondsPerTick accordingly.
double AP_Timer::GetSecondsPerClockTick()
{
  static double sec_per_tick(0);

  if (sec_per_tick == 0.0f)
  {
#ifdef _WIN32
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    sec_per_tick = static_cast<double>(1.0f / frequency.QuadPart);
#elif __APPLE__
    mach_timebase_info_data_t info;
    mach_timebase_info(&info);
    sec_per_tick = info.numer/info.denom;
#else
    struct timespec ts;
    clock_getres(CLOCK_MONOTONIC, &ts);
    sec_per_tick = static_cast<double>(ts.tv_sec);
#endif
  }
  return sec_per_tick;
} // GetSecondsPerClockTick()


//Reset the timer (elapsed time is reset to zero).
//return: reference to the timer object.
AP_Timer & AP_Timer::Reset()
{
  m_fStartTime = 0;
  UpdateTimes();
  m_fStartTime = m_fElapsedTime; 
  m_bTimerIsRunning = true;
  return(*this);
} // Reset()

 
//Stop the timer (elapsed time is frozen).
//return: reference to the timer object.
AP_Timer & AP_Timer::Stop()
{
  if ( m_bTimerIsRunning )
  {
    UpdateTimes();
    m_bTimerIsRunning = false;
  }
  return(*this);
} // Stop()


//Return the currently elapsed time.
//return elapsed time return as double.
double  AP_Timer::ElapsedTime()
{
  if ( m_bTimerIsRunning )
  {
    UpdateTimes();
  }
  return(fabs(m_fElapsedTime));
} // ElapsedTime


//Continue the timer (elapsed time is continues).
//return reference to the timer object.
AP_Timer & AP_Timer::Continue()
{
  if ( !m_bTimerIsRunning )
  {
    double tmp  = m_fElapsedTime;
    UpdateTimes();
    m_fStartTime += m_fElapsedTime  - tmp;
    m_bTimerIsRunning = true;
  }
  return(*this); 
} // Continue


//Set the currently elapsed time.
//param: elapsedTime
//       the time you want the timer to contiue from.
void AP_Timer::SetElapsedTime(const double & elapsedTime)
{
  if ( m_bTimerIsRunning )
  {
    UpdateTimes();
  }
  m_fStartTime += m_fElapsedTime - elapsedTime;
}


// Returns the elapsed time as a string.
// formatted as Hh MMm SS.SSSs
// return: the string with the elapsed time.
std::string AP_Timer::AsString()
{
  if ( m_bTimerIsRunning )
  {
    UpdateTimes();
  }

  size_t uTotalSeconds(static_cast<size_t>(fabs(m_fElapsedTime)));
  size_t uHours(uTotalSeconds/60/60);
  size_t uMinutes(uTotalSeconds/60%60);
  double fSeconds(fabs(m_fElapsedTime) - uTotalSeconds + uTotalSeconds%60);
  char timeStr[1024];
  sprintf(timeStr, "time: %lu h %02lu m %06.3f s", uHours, uMinutes, fabs(fSeconds));
  return( timeStr );

} // AsString()


//Returns the elapsed time as a string.
//formatted as Hh MMm SS.SSSs
//return: the string with the elapsed time.
std::string AP_Timer::AverageTimeAsString(size_t uSteps)
{
  if ( m_bTimerIsRunning )
  {
    UpdateTimes();
  }
  double fAverageElapsedTime(m_fElapsedTime/static_cast<double>(uSteps));
  size_t uTotalSeconds(static_cast<size_t>(fabs(fAverageElapsedTime)));
  size_t uHours(uTotalSeconds/60/60);
  size_t uMinutes(uTotalSeconds/60%60);
  double fSeconds(fabs(fAverageElapsedTime) - uTotalSeconds + uTotalSeconds%60);
  char timeStr[50];
  sprintf(timeStr, "time: %lu h %02lu m %06.3f s", uHours, uMinutes, fabs(fSeconds));
  return( timeStr );
} // AsString()

//Update the elapsed time from system time call.
void  AP_Timer::UpdateTimes()
{
#ifdef _WIN32
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    double endTime = static_cast<double>(now.QuadPart) * m_secondsPerTick;
    m_fElapsedTime = endTime - m_fStartTime;
#elif __APPLE__
    uint64_t now = mach_absolute_time();
    double endTime = static_cast<double>(now) * m_secondsPerTick;
    m_fElapsedTime = endTime - m_fStartTime;
#else
    timeval now;
    gettimeofday(&now, NULL);
    m_fElapsedTime = now.tv_sec - m_fStartTime;
#endif

} // UpdateTimes

std::ostream & operator << (std::ostream & os, AP_Timer & timer)
{
  return( os<<timer.AsString() );
} // operator <<

} //end namespace

