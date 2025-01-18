#ifndef TIMER_HPP
#define TIMER_HPP

#include <ctime>

namespace csmp
  {

  /// Trivial timer class, returns in ~seconds (PL 2011)
  class Timer
    {
    public:
      Timer() : ticks_(0) {}
      ~Timer() {}
      
      void Start() { ticks_ = clock(); }
      double StopClock() { return (clock() - ticks_); }
      double Stop() { return StopClock() / static_cast<double>(CLOCKS_PER_SEC); }

    private:
      Timer( const Timer& );
      Timer& operator = ( const Timer& );

    private:
      clock_t ticks_;
    };

  } // csmp

#endif
