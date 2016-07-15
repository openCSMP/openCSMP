#ifndef COMPUTATIONAL_SETTINGS_H
#define COMPUTATIONAL_SETTINGS_H

#include "CSMP_definitions.h"

namespace csmp {

enum TIME_STEPPING { BACKWARD_EULER, CRANK_NICHOLSON, FORWARD };
enum TIME_STRATEGY { PEDANTIC, CONSERVATIVE, MODERATE, AGGRESSIVE, DARING, PRESCRIBED , CAREFUL};

TIME_STRATEGY parseTimeStrategy( const char* );
const char*   parseTimeStrategy( TIME_STRATEGY );

template<size_t> class Model;

/// base class for generic configurations
class ComputationalSettings {
  public:
    ComputationalSettings();
    virtual ~ComputationalSettings();
    // for output of results
    void           SetOutputTimes( const std::set<double64>& times, bool overwrite=true );
    void           AddOutputTime( double64 time );
    
    bool           IsOutputTime( double64 time, double64 tolerance ) const;
    
    double64       TimeToNearestOutputTime( double64 current_time ) const;
    double64       NearestOutputTime( double64 current_time ) const;
    size_t         OutputTimePosition(double64 time);

    /// For output of monitoring data

    void           SetMonitorTimes( const std::set<double64>& times, bool overwrite=true );
    void           AddMonitorTime( double64 time );

    bool           IsMonitorTime( double64 time, double64 tolerance ) const;

    double64       TimeToNearestMonitorTime( double64 current_time ) const;
    double64       NearestMonitorTime( double64 current_time ) const;
    size_t         MonitorTimePosition(double64 time);

    /// Time control
    
    void           Duration( double64 duration );
    double64       Duration() const;
    
    void           TimeSteppingApproach( TIME_STRATEGY time_stepping_approach );
    TIME_STRATEGY  TimeSteppingApproach() const;
    
    void           TimeIncrement( double64 );
    double64       TimeIncrement() const;
  
    /// based on the time-stepping strategy, pressure-step and saturation step multipliers are suggested 
    void EstablishMultipliers( double64& pf_multiplier, 
                               double64& adv_multiplier ) const;

    /// algorithm that monitors the flow velocity change from timestep to timestep and suggests a timeincrement on this basis
    template<size_t dim>
    double64 TimeIncrementFromVelocityChange( const Model<dim>& sg, 
                                              double64 log_velocity_change, 
                                              double64 current_delta_t, 
                                              double64 max_delta_t );
    void            RemoveAllOutputtimes(){ output_times_.clear();}
    void            RemoveAllMonitortimes(){ output_times_.clear();}

    /// step through output times
    double64                           PopOutputTime();
    std::set<double64>::const_iterator OutputTimesBegin() const;
    std::set<double64>::const_iterator OutputTimesEnd() const;

    /// step through monitor times
    double64                           PopMonitorTime();
    std::set<double64>::const_iterator MonitorTimesBegin() const;
    std::set<double64>::const_iterator MonitorTimesEnd() const;
    
    virtual void   Out() const;
    virtual void   Out( const char* filename ) const;
  
  private:
    TIME_STRATEGY       time_strategy_;
    double64            run_duration_;
    double64            time_increment_;
    std::set<double64>  output_times_;
    std::set<double64>  monitor_times_;
    std::set<double64>::const_iterator  it_; // points to next output time
    std::set<double64>::const_iterator  itm_; // points to next monitor time
};



inline std::set<double64>::const_iterator ComputationalSettings::OutputTimesBegin() const
 { return output_times_.begin(); }
 
inline std::set<double64>::const_iterator ComputationalSettings::OutputTimesEnd() const
 { return output_times_.end(); }

inline std::set<double64>::const_iterator ComputationalSettings::MonitorTimesBegin() const
 { return monitor_times_.begin(); }

inline std::set<double64>::const_iterator ComputationalSettings::MonitorTimesEnd() const
 { return monitor_times_.end(); }


} // end namespace csmp

/**
 
@class ComputationalSettings ComputationalSettings "utilities/CommandLineParser.h"
@author S.K. Matthaei
@date 2003
 
@section motivation Motivation
 
Record run configuration and parameters for later use to steer the 
simulation as it runs.  
 
 
@section design Design Intent

To have everything that concerns the time-stepping control for a 
run in a single object.  
 
 
@section applicability Applicability

Any CSMP model.
 
 
@section consequences Consequences

You can inherit from this baseclass.
 
 
@section examples Application Examples

Set the object up from an input file and use it to control the 
time stepping and output of data from a run:  

 @code
  ComputationalSettings  run_settings;
  
  double64        duration(run_settings.Duration());
  
  double64 time_until_output = run_settings.NearestOutputTime( model_time );
  
  if ( run_settings.IsOutputTime( model_time, 1. )  ) {
    ...
    }
@endcode
*/

#endif
