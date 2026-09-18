// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COMPUTATIONAL_SETTINGS_H
#define COMPUTATIONAL_SETTINGS_H

#include "CSMP_definitions.h"

namespace csmp {

enum TIME_STEPPING { BACKWARD_EULER, CRANK_NICHOLSON, FORWARD };
enum TIME_STRATEGY { PEDANTIC, CONSERVATIVE, MODERATE, AGGRESSIVE, DARING, PRESCRIBED, CAREFUL, ASYNCHRONOUS};

TIME_STRATEGY parseTimeStrategy( const char* );
const char*   parseTimeStrategy( TIME_STRATEGY );

template<uint32_t> class Model;

/// base class for generic configurations
class ComputationalSettings {
  public:
    ComputationalSettings();
    virtual ~ComputationalSettings();
  
    // for output of results
  
    void           SetOutputTimes( const std::set<double>& times, bool overwrite=true );
    void           AddOutputTime( double time );
    
    bool           IsOutputTime( double time, double tolerance ) const;
    
    double       TimeToNearestOutputTime( double current_time ) const;
    double       NearestOutputTime( double current_time ) const;
    size_t         OutputTimePosition(double time);

    // For output of monitoring data

    void           SetMonitorTimes( const std::set<double>& times, bool overwrite=true );
    void           AddMonitorTime( double time );

    bool           IsMonitorTime( double time, double tolerance ) const;

    double       TimeToNearestMonitorTime( double current_time ) const;
    double       NearestMonitorTime( double current_time ) const;
    size_t         MonitorTimePosition(double time);

    // Time control
    
    void           Duration( double duration );
    double       Duration() const;
    
    void           TimeSteppingApproach( TIME_STRATEGY time_stepping_approach );
    TIME_STRATEGY  TimeSteppingApproach() const;
    
    void           TimeIncrement( double );
    double       TimeIncrement() const;
  
    /// based on the time-stepping strategy, pressure-step and saturation step multipliers are suggested 
    void EstablishMultipliers( double& pf_multiplier, 
                               double& adv_multiplier ) const;

    /// algorithm that monitors the flow velocity change from timestep to timestep and suggests a timeincrement on this basis
    template<uint32_t dim>
    double TimeIncrementFromVelocityChange( const Model<dim>& sg, 
                                              double log_velocity_change, 
                                              double current_delta_t, 
                                              double max_delta_t );
  
    void            RemoveAllOutputtimes() { output_times_.clear(); }
    void            RemoveAllMonitortimes(){ output_times_.clear(); }

    /// step through output times
    double                           PopOutputTime();
    std::set<double>::const_iterator OutputTimesBegin() const;
    std::set<double>::const_iterator OutputTimesEnd() const;

    /// step through monitor times
    double                           PopMonitorTime();
    std::set<double>::const_iterator MonitorTimesBegin() const;
    std::set<double>::const_iterator MonitorTimesEnd() const;
    
    virtual void   Out() const;
    virtual void   Out( const char* filename ) const;
  
  private:
    TIME_STRATEGY       time_strategy_;
    double            run_duration_;
    double            time_increment_;
    std::set<double>  output_times_;
    std::set<double>  monitor_times_;
    std::set<double>::const_iterator  it_; // points to next output time
    std::set<double>::const_iterator  itm_; // points to next monitor time
};

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
  
  double        duration(run_settings.Duration());
  
  double time_until_output = run_settings.NearestOutputTime( model_time );
  
  if ( run_settings.IsOutputTime( model_time, 1. )  ) {
    ...
    }
@endcode
*/

#endif
