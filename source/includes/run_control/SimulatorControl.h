#ifndef SIMULATORCONTROL_H
#define SIMULATORCONTROL_H
#include "ErrorHandler.h"
#include "SimulatorSetup.h"
#include "SimulationSignalHandler.h"
#include "ComputationalSettings.h"
#include "VTU_Interface.h"

namespace csmp{

enum WELL_MODE { SHUTIN, PRODUCER, INJECTOR};

/** @author Julian E. Mindel
 */

template <size_t dim>
class SimulatorControl
{
public:
    SimulatorControl();
    ~SimulatorControl();
    SimulatorControl(SimulatorSetup<dim>* ss);
    void CatchSignals();
    SimulationSignalHandler& GetSignalHandler(){return sig_;}
    void OutputRestartFile(bool manually_triggered=false);
    void Initialize();
    std::string GetName(){return project_name_;}
    SimulatorSetup<dim>* GetSS(){ return this->simulator_setup_;}
    ComputationalSettings& RunSettings();
    void ManageWellRates();
    virtual void Run();

    bool CheckForDirichletEssentialCondition(csmp::Index key);
    void PrintRangeOfVariablesToScreen();
    bool ReadControlOptions();
    virtual bool ReadSimulatorSpecificControlOptions()=0;

    // specific necessary property keys to be accessed from derived classes.
    Index& GetThicknessKey(){return thickness_;}
    Index& GetWellRadiusKey(){return wellRadius_;}
    virtual void LoadKeysAndNames()=0;

    double64 MaxDifferenceScalarNodalProperty(Index &snp1Key, Index& snp2Key );


    // flags and value control on places
    void SetBoundaryFlag(const char* prop_name, std::string boundary, VARIABLE_FLAG flag);
    void SetBoundaryValue(const char* prop_name,std::string boundary,double64 value,VARIABLE_FLAG flag);
    void SetRegionFlag(const char* prop_name,std::string region,VARIABLE_FLAG flag);
    void SetRegionValue(std::string prop_name,std::string region,double64 value,VARIABLE_FLAG flag);
    void SetFlagNearestToPoint(Index key,double64 x, double64 y, double64 z, VARIABLE_FLAG flag);
    void SetValueNearestToPoint(Index key,double64 x, double64 y, double64 z, double64 value,VARIABLE_FLAG flag,SUBDOMAIN_PART sub=COMPLETE);
    void SetPropertyFlagToTopBoundaries(const char* property_name, VARIABLE_FLAG flag);

    // time control
    double64 GetSimulationEndTime(){return simulationEndTime_;}
    double64 GetSimulationStartTime(){return simulationStartTime_;}
    std::vector<std::pair<double64,std::string> >& GetIntervals(){return interval_start_times_;}
    double64 GetIntervalEndTime(){return intervalEndTime_;}
    double64 GetIntervalStartTime(){return intervalStartTime_;}

    size_t GetIntervalNumber(double64 current_time);
    std::string GetIntervalName(double64 current_time);
    void InsertNewTimeInterval(std::string interval_name,double64 starttime);
    double64 GetSimulationTime(){return simulationTime_;}
    //double64 GetGlobalTimestep(){return time_increment_;}

    void SetThicknessFactorFromWellRadii();

    void SetSimulationEndTime(double64 t){simulationEndTime_=t;}
    void SetSimulationStartTime(double64 t){simulationStartTime_=t;}
    void SetIntervalEndTime(double64 t){intervalEndTime_=t;}
    void SetIntervalStartTime(double64 t){intervalStartTime_=t;}
    void SetSimulationTime(double64 t);
    void SetGlobalTimestep(double64 t){time_increment_=t;}
    double64 GetCurrentIntervalEndTime(double64 current_simulation_time);
    double64 GetCurrentIntervalStartTime(double64 current_simulation_time);
    double64 CalculateCurrentRunTime();
    void SetOutputTimes(size_t n, double64 starttime=0.0);
    void SetMonitorTimes(size_t n, double64 starttime=0.0);
    void LoadMonitoringTimesFromRestartedModel();
    virtual void ComputeTimeIncrement()=0;
    void TrimTimeIncrementWithMonitorTime();
    void TrimTimeIncrementWithOutputTime();
    void TrimTimeIncrementWithEndOfInterval();
    double64 TimeIncrement(){return time_increment_;}
    void TimeIncrement(double64 time_increment){this->time_increment_=time_increment;}
    double64 GetInitialSimulatorStartTimeFromIntervals();
    void OutputSimulationTimeToScreen(double64 time);

    void SyncOutputAntMonitoringTimesToModel();
    bool MonitorAllTimesteps() {return this->monitor_all_timesteps_;}
    void MonitorAllTimesteps(bool mon_all_timesteps) {this->monitor_all_timesteps_=mon_all_timesteps;}

    void SetReferencePointValuesModelInterior();
    void SetReferencePointValuesModelPerimeter();
    void BeginInterval();

    // interactive
    bool InteractiveInputNodeValue(csmp::Index key, SUBDOMAIN_PART sub = INTERIOR);
    bool AskUserToContinue();
    bool YesOrNo();
    bool CheckTimeIntervals();
    void InteractivelyManageEssentialConditions();
    double64 GetRealInput(bool success=true);
    int32 GetIntegerInput();
    std::string GetStringInput();
    std::string ParseTrueOrFalse(bool trueorfalse);
    std::string ParseWellMode(WELL_MODE mode);
    WELL_MODE ParseIntegerToWellMode(int32 i);
    std::multimap<std::string,std::pair<std::string ,WELL_MODE> >& GetWellModesPerInterval(){return wellmodes_per_interval_;}
    std::multimap <std::string,std::pair<std::vector<double64>,double64> >& GetReferencePointsAndValues(){return this->refpoints_;}

    // output
    VTU_Interface<dim>* GetVTU_Interface(){return vtu_;}
    void OutputVariableToRegionVTUFiles();
    virtual void OutputSimulatorSpecificControlFileSection()=0;
    void OutputSampleControlFile();

    void UpdateCurrentRunTime();
    void Output();
    void QuietSolvers();
    void Verbose(bool verbose){verbose_=verbose;error_handler_.Verbose( verbose );}
    bool Verbose(){return verbose_;}
    void TriggerSAMGSetup(bool trigger){triggerSAMGsetup_=trigger;}
    bool TriggerSAMGSetup(){return triggerSAMGsetup_;}
    bool& Restart(){return this->restart_;}
    void Restart(bool restart){this->restart_=restart;}


protected:
    clock_t runTimeStart_;
    PropertyDatabase<dim>& pdb_;
    ErrorHandler& error_handler_;
    double64 initial_omp_wtime_;

private:
    bool stop_;
    bool restart_;
    std::string project_name_;
    SimulatorSetup<dim>* simulator_setup_;
    bool monitor_all_timesteps_;

    static size_t outputVTU_counter_;
    // for simulated time
    double64 simulationStartTime_,simulationEndTime_;
    double64 simulationTime_,time_increment_;
    double64 intervalStartTime_,intervalEndTime_;
    std::vector<std::pair<double64,std::string> > interval_start_times_;

    // for actual cpu run time.
    double64 restartFileOuputRunTimeInterval_;
    clock_t startClockTicks_;
    double64 runtimeOnLastRestart_;

    // keys and output
    bool verbose_;
    Index run_time_Key_;
    Index model_time_Key_;
    Index vtu_frames_Key_;
    Index monitor_frames_Key_;
    Index thickness_,wellRadius_;
    VTU_Interface<dim>* vtu_;
    SimulationSignalHandler sig_;

    std::multimap<std::string,std::pair<std::string ,WELL_MODE> > wellmodes_per_interval_;
    std::multimap <std::string,std::pair<std::vector<double64>,double64> > refpoints_;
    bool triggerSAMGsetup_;
};
}

#endif // SIMULATORCONTROL_H
