#include "ComputationalSettings.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"

using namespace std;

namespace csmp {


/**

The default constructor sets the time stepping strategy to conservative
and the run duration to zero. No output times are entered into the 
output time set.  

The default time increment is 60 seconds.
*/
ComputationalSettings::ComputationalSettings()
    : time_strategy_(CONSERVATIVE),
      time_increment_(60.), // 1 minute = default
      run_duration_(0.)
{
    it_ = output_times_.begin();
    itm_ = monitor_times_.begin();
}


ComputationalSettings::~ComputationalSettings()
{
}




std::set<double64>::const_iterator ComputationalSettings::OutputTimesBegin() const
 { return output_times_.begin(); }
 
std::set<double64>::const_iterator ComputationalSettings::OutputTimesEnd() const
 { return output_times_.end(); }

std::set<double64>::const_iterator ComputationalSettings::MonitorTimesBegin() const
 { return monitor_times_.begin(); }

std::set<double64>::const_iterator ComputationalSettings::MonitorTimesEnd() const
 { return monitor_times_.end(); }



/**

Changes the supplied arguments pf- and advection multiplier to reflect
the user-defined timestepping approach; options are:

pedantic, conservative, moderate=normal, aggressive, daring, prescribed  

@section arguments Input Arguments 

The strategy chosen by the user with the options:

@code
PEDANTIC, CONSERVATIVE, MODERATE, AGGRESSIVE, DARING.
@endcode

@return The fluid pressure solve frequency multiplier and a multiplier that
decides by how much the Courant number shall be overstepped. 
The latter must be <= 1 if the explicit advection scheme is used.  
*/
void ComputationalSettings::EstablishMultipliers( double64& pf_multiplier, 
                                                  double64& adv_multiplier ) const
{
    switch ( time_strategy_ ) {
    case PEDANTIC:  // this will invoke the explicit second-order scheme
        pf_multiplier  = 1.;
        adv_multiplier = 0.999999;
        break;

    case CONSERVATIVE: // IMP-IMS basic version
        pf_multiplier  = 2.;
        adv_multiplier = 1.;
        break;

    case MODERATE:
        pf_multiplier  = 3.;
        adv_multiplier = 2.;
        break;

    case AGGRESSIVE:
        pf_multiplier  = 10.;
        adv_multiplier =  5.;
        break;
    case DARING:
        pf_multiplier  = 20.;
        adv_multiplier =  8.;
        break;
    case PRESCRIBED:
        pf_multiplier  = 1.;
        adv_multiplier = 1.;
        break;
    case CAREFUL:
        pf_multiplier  = 1.;
        adv_multiplier = 0.7;
        break;
    case ASYNCHRONOUS:  // this will invoke the explicit DES scheme
        pf_multiplier  = 1.;
        adv_multiplier = 0.5;
        break; 
    default:
        cout <<"\nComputationalSettings::EstablishMultipliers: Directive not recognized."<< endl;
    }

} // end EstablishMultipliers




/** Permits the user to enter a series of unique output times as supplied
in the argument set.  Either over-writing or adding to existing list is optional.
Overwrite is the default settings.
*/
void ComputationalSettings::SetOutputTimes( const set<double64>& times ,bool overwrite)
{
    if ( times.empty() )
        cout <<"\nComputationalSettings::SetOutputTimes: Supplied set is empty."<< endl;
    
    if (overwrite)
        output_times_ = times;
    else
        for (set<double64>::iterator sit = times.begin(); sit!=times.end();sit++)
            output_times_.insert(*sit);
}

/** Permits the user to enter a series of unique monitoring times as supplied
in the argument set.  Either over-writing or adding to existing list is optional.
Overwrite is the default settings.
*/
void ComputationalSettings::SetMonitorTimes( const set<double64>& times ,bool overwrite)
{
    if ( times.empty() )
        cout <<"\nComputationalSettings::SetMonitorTimes: Supplied set is empty."<< endl;

    if (overwrite)
        monitor_times_ = times;
    else
        for (set<double64>::iterator sit = times.begin(); sit!=times.end();sit++)
            monitor_times_.insert(*sit);
}


/**

Checks whether the nearest time-increment lies within the user-specified
tolerance. If so the method returns true, else false.  
*/
bool ComputationalSettings::IsOutputTime( double64 value, double64 tol ) const
{
    // brute force approach
    //    for ( set<double64>::const_iterator it=output_times_.begin();
    //          it!=output_times_.end(); it++ ) {
    //        if ( *it <= (value + tol) && *it >= (value - tol) ) return true;
    //    }
    set<double64>  time_difference;

    for ( set<double64>::const_iterator it=output_times_.begin();
          it!=output_times_.end(); it++ )
        // if the output time has no occurred yet
        if ( (*it) - value >= 0. )
            time_difference.insert( fabs((*it) - value) );

    if ( time_difference.empty() ) return false;

    return (*time_difference.begin()) < tol;
}

/**

Checks whether the nearest time-increment lies within the user-specified
tolerance. If so the method returns true, else false.
*/
bool ComputationalSettings::IsMonitorTime( double64 value, double64 tol ) const
{
    // brute force approach
    //    for ( set<double64>::const_iterator it=output_times_.begin();
    //          it!=output_times_.end(); it++ ) {
    //        if ( *it <= (value + tol) && *it >= (value - tol) ) return true;
    //    }
    set<double64>  time_difference;

    for ( set<double64>::const_iterator it=monitor_times_.begin();
          it!=monitor_times_.end(); it++ )
        // if the output time has no occurred yet
        if ( (*it) - value >= 0. )
            time_difference.insert( fabs((*it) - value) );

    if ( time_difference.empty() ) return false;

    return (*time_difference.begin()) < tol;
}


/**

Returns time difference to closest output time as currently stored by
the object. If the output time occurs in the past, the time difference 
would be negative but it is not considered. Only output times which
occur in the future are considered and they are all positive.  
*/
double64  ComputationalSettings::TimeToNearestOutputTime( double64 current_time ) const
{
    set<double64>  time_difference;
    
    for ( set<double64>::const_iterator it=output_times_.begin();
          it!=output_times_.end(); it++ )
        // if the output time has no occurred yet
        if ( (*it) - current_time > 0. )
            time_difference.insert( (*it) - current_time );

    if ( time_difference.empty() ) return run_duration_;

    return (*time_difference.begin());
    
} // end NearestOutputTime

/**

Returns time difference to closest monitor time as currently stored by
the object. If the output time occurs in the past, the time difference
would be negative but it is not considered. Only output times which
occur in the future are considered and they are all positive.
*/
double64  ComputationalSettings::TimeToNearestMonitorTime( double64 current_time ) const
{
    set<double64>  time_difference;

    for ( set<double64>::const_iterator it=monitor_times_.begin();
          it!=monitor_times_.end(); it++ )
        // if the output time has no occurred yet
        if ( (*it) - current_time > 0. )
            time_difference.insert( (*it) - current_time );

    if ( time_difference.empty() ) return run_duration_;

    return (*time_difference.begin());

} // end TimeToNearestOutputTime


/**

Returns the nearest registered output time step with respect to the current_time.
 
*/
double64  ComputationalSettings::NearestOutputTime( double64 current_time ) const
{
    set<double64>  nextOutputTimes;
    
    for ( set<double64>::const_iterator it=output_times_.begin();
          it!=output_times_.end(); it++ )
        // if the output time has no occurred yet
        if ( (*it) - current_time > 0. )
            nextOutputTimes.insert( (*it) );

    if ( nextOutputTimes.empty() ) return run_duration_;

    return (*nextOutputTimes.begin());
    
} // end NearestOutputTime

/**

Returns the nearest registered monitoring time step with respect to the current_time.

*/
double64  ComputationalSettings::NearestMonitorTime( double64 current_time ) const
{
    set<double64>  nextOutputTimes;

    for ( set<double64>::const_iterator it=monitor_times_.begin();
          it!=monitor_times_.end(); it++ )
        // if the output time has no occurred yet
        if ( (*it) - current_time > 0. )
            nextOutputTimes.insert( (*it) );

    if ( nextOutputTimes.empty() ) return run_duration_;

    return (*nextOutputTimes.begin());

} // end NearestOutputTime


/** Adds a new output time increment to the current set of output times.
*/
void ComputationalSettings::AddOutputTime( double64 time )
{
    output_times_.insert( time );
}

/** Adds a new output time increment to the current set of output times.
*/
void ComputationalSettings::AddMonitorTime( double64 time )
{
    monitor_times_.insert( time );
}

/** Sets the duration of the current run to the argument value.
*/
void ComputationalSettings::Duration( double64 duration )
{
    if ( duration < 0. ) {
        cout <<"\nComputationalSettings::Duration: Can't assign negative value: ";
        cout << duration << endl;
        return;
    }
    if ( duration > static_cast<double64>(145065600000000000.) ) {
        cout <<"\nComputationalSettings::Duration: Longer than age of earth ?: ";
        cout << duration << endl;
        return;
    }
    run_duration_ = duration;
}


/** Returns the scheduled duration of the current simulation.
*/
double64 ComputationalSettings::Duration() const
{
    return run_duration_;
}



/** Sets the time-stepping strategy parameter.
*/
void ComputationalSettings::TimeSteppingApproach( TIME_STRATEGY time_stepping_approach )
{
    time_strategy_ = time_stepping_approach;
}



/**

Returns the current time stepping approach as captured by one of the 
switches of the enumeration TIME_STRATEGY:

@code
PEDANTIC, CONSERVATIVE, MODERATE, AGGRESSIVE, DARING
@endcode
 */
TIME_STRATEGY ComputationalSettings::TimeSteppingApproach() const
{
    return time_strategy_;
}


size_t ComputationalSettings::OutputTimePosition(double64 time)
{
    return distance(output_times_.begin(),output_times_.find(NearestOutputTime(time)));
}

size_t ComputationalSettings::MonitorTimePosition(double64 time)
{
    return distance(output_times_.begin(),output_times_.find(NearestOutputTime(time)));
}

void ComputationalSettings::TimeIncrement( double64 dt )
{
    time_increment_ = dt;
}

double64 ComputationalSettings::TimeIncrement() const
{
    return time_increment_;
}


/**

Calculates the maximum permissible time increment based on the current change
in velocity and the maximum allowed change in velocity for the next time step. 
If the velocity does not change, the maximum time increment provided by the user
is returned. The variables "volume flux" and "previous volume flux" are compared
and must be computed by the user in the main program. The velocity change is only
computed if the velocity is below a minimum threshold (1.0e-15 m/s).  

@section arguments Input Arguments 

- The allowed velocity change (in log units)

- The current time increment (in sec)

- The maximum time increment (in sec).

*/
template<size_t dim>
double64 ComputationalSettings::TimeIncrementFromVelocityChange( const Model<dim>& sg, 
                                                                 double64 log_velocity_change,
                                                                 double64 current_delta_t,
                                                                 double64 max_delta_t )
{
    csmp::Index  vfc_key(sg.Database().StorageKey("volume flux")),
            vfp_key(sg.Database().StorageKey("previous volume flux"));
    double64           vfc, vfp, rate, time_increment(max_delta_t), local_increment;
    const double64     zero(0.), min_change(0.01), max_change(100.);
    const double64     min_velocity(1.0e-15);
    const double64     log_min_velocity(log10(min_velocity));

    // test if a useful velocity change is provided
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( log_velocity_change <= zero ) {
        csmp_error.notice( ERROR, "ComputationalSettings::TimeIncrementFromVelocityChange",
                           "The desired change in velocity (in log units) is zero or below zero",
                           "Returning user specified maximum time increment");
        return max_delta_t;
    }
    if ( log_velocity_change < min_change ) {
        csmp_error.notice( WARNING, "ComputationalSettings::TimeIncrementFromVelocityChange",
                           "The desired change in velocity (in log units) is low",
                           "Time increment may become excessively small" );
    }
    if ( log_velocity_change > max_change ) {
        csmp_error.notice( WARNING, "ComputationalSettings::TimeIncrementFromVelocityChange",
                           "The desired change in velocity (in log units) is large",
                           "Influence of velocity change may not be reflected in time increment");
    }

    // loop over elements and compare volume fluxes
    typename std::vector<Element<dim>*>::const_iterator firs = sg.Region("Model").ElementsBegin();
    typename std::vector<Element<dim>*>::const_iterator last = sg.Region("Model").ElementsEnd();
    while ( firs != last ) {
        // get current and previous volume flux
        vfc = (*firs)->Read( vfc_key );
        vfp = (*firs)->Read( vfp_key );

        if ( vfc > min_velocity ) vfc = std::log10(vfc);
        else                      vfc = log_min_velocity;
        if ( vfp > min_velocity ) vfp = std::log10(vfp);
        else                      vfp = log_min_velocity;

        // compute local rate
        rate  = std::fabs(vfc-vfp);
        rate /= current_delta_t;

        // compute local time increment if rate of velocity change is not zero and something flows
        if ( rate != zero and vfc > log_min_velocity ) {
            local_increment  = log_velocity_change; // maximum allowed velocity change (in log units)
            local_increment /= rate;                // current rate of velocity change

            // check if local increment is less than current increment
            time_increment = std::min( time_increment, local_increment );
        }
        firs++;
    }

    return time_increment;
}

template double64
ComputationalSettings::TimeIncrementFromVelocityChange<1U>( const Model<1U>&, 
double64, double64, double64 );
template double64
ComputationalSettings::TimeIncrementFromVelocityChange<2U>( const Model<2U>&, 
double64, double64, double64 );
template double64
ComputationalSettings::TimeIncrementFromVelocityChange<3U>( const Model<3U>&, 
double64, double64, double64 );





/**

Increments the time-step iterator to the next stored output time and
returns the corresponding time value.  
*/
double64 ComputationalSettings::PopOutputTime()
{
    if ( it_ == output_times_.end() ) it_ = output_times_.begin();

    double64 output_time = *it_++;

    return output_time;
}



/** Prints the current computational settings to screen.

*/
void ComputationalSettings::Out() const
{
    cout <<"\nComputationalSettings::Out:";
    cout <<"\ntime strategy: "<< parseTimeStrategy( time_strategy_ );
    cout <<"\nrun duration:  "<< run_duration_;
    cout <<"\noutput times:  ";
    for ( std::set<double64>::const_iterator oit=output_times_.begin();
          oit!=output_times_.end(); oit++ )
        cout << (*oit) <<", ";
    cout << endl << endl;

    if (!monitor_times_.empty()){
        cout <<"\nmonitor times:  ";
        for ( std::set<double64>::const_iterator oit=monitor_times_.begin();
              oit!=monitor_times_.end(); oit++ )
            cout << (*oit) <<", ";
        cout << endl << endl;
    }
}



/**

Writes the current run parameters to the ASCII text file that the user
specifies as the first method argument.  
*/
void ComputationalSettings::Out( const char* filename ) const
{
    ofstream  ofs( filename );
    if ( !ofs.is_open() ) {
        cout <<"\nComputationalSettings::Out: Can't open output file: '"<< filename <<"'"<< endl;
        return;
    }
    
    ofs <<"\nComputationalSettings::Out: To '"<< filename <<"'";
    ofs <<"\ntime strategy: "<< parseTimeStrategy( time_strategy_ );
    ofs <<"\nrun duration:  "<< run_duration_;
    ofs <<"\noutput times:  ";
    for ( std::set<double64>::const_iterator oit=output_times_.begin();
          oit!=output_times_.end(); oit++ )
        ofs << (*oit) <<", ";
    ofs << endl << endl;
    
    ofs.close();
}




/**

Translates time stepping strategy as expressed by a character string
into the enumeration TIME_STRATEGY.  
*/
TIME_STRATEGY parseTimeStrategy( const char* time_strategy )
{
    string strategy( time_strategy );
    
    if      ( strategy == "PEDANTIC" )     return PEDANTIC;
    else if ( strategy == "CONSERVATIVE" ) return CONSERVATIVE;
    else if ( strategy == "MODERATE" )     return MODERATE;
    else if ( strategy == "AGGRESSIVE" )   return AGGRESSIVE;
    else if ( strategy == "DARING" )       return DARING;
    else if ( strategy == "pedantic" )     return PEDANTIC;
    else if ( strategy == "conservative" ) return CONSERVATIVE;
    else if ( strategy == "moderate" )     return MODERATE;
    else if ( strategy == "aggressive" )   return AGGRESSIVE;
    else if ( strategy == "daring" )       return DARING;
    else if ( strategy == "prescribed" )   return PRESCRIBED;
    else if ( strategy == "PRESCRIBED" )   return PRESCRIBED;
    else if ( strategy == "careful" )   return CAREFUL;
    else if ( strategy == "CAREFUL" )   return CAREFUL;
    else if ( strategy == "asynchronous" )   return ASYNCHRONOUS;
    else if ( strategy == "ASYNCHRONOUS" )   return ASYNCHRONOUS;    

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    csmp_error.notice( FATAL_ERROR, "parseTimeStrategy(const char*):", "Unable to parse time-stepping strategy", time_strategy );

    return CONSERVATIVE; // the default value
}



/**

Translates time-stepping strategy as expressed by the enumeration
TIME_STRATEGY into a character string. 
 
*/
const char* parseTimeStrategy( TIME_STRATEGY  time_strategy )
{
    if      ( time_strategy == PEDANTIC )     return "PEDANTIC";
    else if ( time_strategy == CONSERVATIVE ) return "CONSERVATIVE";
    else if ( time_strategy == MODERATE )     return "MODERATE";
    else if ( time_strategy == AGGRESSIVE )   return "AGGRESSIVE";
    else if ( time_strategy == DARING )       return "DARING";
    else if ( time_strategy == PRESCRIBED )   return "PRESCRIBED";
    else if ( time_strategy == CAREFUL )      return "CAREFUL";
    else if ( time_strategy == ASYNCHRONOUS ) return "ASYNCHRONOUS";

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    csmp_error.notice( FATAL_ERROR, "parseTimeStrategy(TIME_STRATEGY):","Unable to parse time-stepping strategy." );

    return "conservative";
}

} // end namespace csmp













