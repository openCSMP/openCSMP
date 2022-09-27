#include "SimulatorControl.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#include <cassert>

using namespace std;

namespace csmp {

template <uint32_t dim>
size_t SimulatorControl<dim>::outputVTU_counter_=0;

template <uint32_t dim>
SimulatorControl<dim>::SimulatorControl(SimulatorSetup<dim>* ss):
    simulator_setup_(ss),restart_(false),verbose_(ss->Verbose()),
    project_name_(ss->GetProjectName()),pdb_(ss->GetModel()->Database()),
    vtu_frames_Key_(simulator_setup_->GetModel()->Database().StorageKey("vtu frames")),
    monitor_frames_Key_(simulator_setup_->GetModel()->Database().StorageKey("monitor frames")),
    run_time_Key_(simulator_setup_->GetModel()->Database().StorageKey("run time")),
    model_time_Key_(simulator_setup_->GetModel()->Database().StorageKey("model time")),
    thickness_(this->pdb_.StorageKey("thickness")),
    wellRadius_(this->pdb_.StorageKey("well radius")),
    vtu_(NULL),
    restartFileOuputRunTimeInterval_(10e32),
    intervalEndTime_(0.0),intervalStartTime_(0.0),
    simulationStartTime_(0.0),simulationEndTime_(0.0),simulationTime_(0.0),
    time_increment_(0.0),triggerSAMGsetup_(false),
    monitor_all_timesteps_(false),
    error_handler_( ErrorHandler::Instance() )
{
}

template <uint32_t dim>
void SimulatorControl<dim>::Initialize()
{
    cout<<" *****************************************************************"<<endl;
    cout<<" Initializing Base Simulator Control."<<endl;
    this->Restart(this->GetSS()->Restart());

    this->ReadControlOptions();

    ReadSimulatorSpecificControlOptions();

    LoadKeysAndNames();

    this->GetSS()->CheckInputRanges();

    vtu_= new VTU_Interface<dim>(*this->simulator_setup_->GetModel(),this->GetName().c_str());
    vtu_->OmitZeroInFileName(false);

    this->SetThicknessFactorFromWellRadii();

    this->startClockTicks_=clock();
    this->runtimeOnLastRestart_=0.0;
#ifdef _OPENMP
    initial_omp_wtime_=omp_get_wtime(); // open mp wall time in seconds from some reference point.
#endif

    this->GetSS()->LoadMonitor();

    if (!this->Restart())
        this->SetSimulationTime(this->GetInitialSimulatorStartTimeFromIntervals());
    else
        this->SetSimulationTime(this->GetSS()->GetModel()->Read(this->GetSS()->GetModel()->Database().StorageKey("model time")));

    cout<<" Done Initializing Base Simulator Control."<<endl;
    cout<<" *****************************************************************"<<endl;
}


template<uint32_t dim>
double SimulatorControl<dim>::GetInitialSimulatorStartTimeFromIntervals()
{
    double time = 10e52;
    for (auto mit = this->GetIntervals().begin(); mit != this->GetIntervals().end();mit++){
        if (mit->first < time )
            time = mit-> first;
    }
    return time;
}


template<uint32_t dim>
string SimulatorControl<dim>::ParseTrueOrFalse(bool trueorfalse)
{
    if (trueorfalse)
        return "true";
    else
        return "false";
}

template<uint32_t dim>
string SimulatorControl<dim>::ParseWellMode(WELL_MODE mode)
{
    if (mode==SHUTIN)
        return "SHUTOFF";
    else if (mode==PRODUCER)
        return "PRODUCER";
    else if (mode==INJECTOR)
        return "INJECTOR";
    else
        return "UNKNOWN";
}

template<uint32_t dim>
WELL_MODE SimulatorControl<dim>::ParseIntegerToWellMode(int32_t i)
{
    if (i==0)
        return SHUTIN;
    else if (i==1)
        return PRODUCER;
    else if (i==2)
        return INJECTOR;
    else
        return SHUTIN; // for any other mode, the well is simply parsed as SHUT IN.
}

template<uint32_t dim>
size_t SimulatorControl<dim>::GetIntervalNumber(double current_time){
    //    auto ith = std::lower_bound(this->GetIntervals().begin(),this->GetIntervals().end(),make_pair(current_time,string("test_string")),
    //                                [](pair<double,string> lhs, pair<double,string> rhs) -> bool { return lhs.second < rhs.second; });
    auto ith = lower_bound( this->GetIntervals().begin(),this->GetIntervals().end(), make_pair(current_time,string("test_string")),
    
                                [](pair<double,string> lhs, pair<double,string> rhs) -> bool { return lhs.first < rhs.first; } );

    size_t int_number = ith - this->GetIntervals().begin();
    return int_number;
}

template<uint32_t dim>
string SimulatorControl<dim>::GetIntervalName(double current_time){

//    auto ith = std::lower_bound(this->GetIntervals().begin(),this->GetIntervals().end(),make_pair(current_time,string("test_string")),
//                                [](pair<double,string> lhs, pair<double,string> rhs) -> bool { return lhs.first < rhs.first; });
    string int_name = "unknown";
    size_t ith = 0;
    for (auto it = this->GetIntervals().begin(); it !=this->GetIntervals().end(); it++){
        if (it->first >= current_time){
            int_name = it -> second;
            break;
        }
        ith++;
    }
    if (ith == this->GetIntervals().size())
        int_name = this->GetIntervals().back().second;
//    cout<<"interval name :"<<int_name<<endl;
//    cin.get();
    return int_name;
}

template <uint32_t dim>
void SimulatorControl<dim>::QuietSolvers()
{
#ifdef CSMP_WITH_SAMG_SOLVER
    auto samgbeg=this->GetSS()->GetAllSAMGSettings().begin();
    auto samgend=this->GetSS()->GetAllSAMGSettings().end();

    for (auto samgit=samgbeg;samgit!=samgend;samgit++){
        samgit->second->Set_idmp(-1);  // minimal screen output
        samgit->second->Set_iout1(-1); // minimal screen output
        samgit->second->Set_iout2(-1); // minimal screen output
        samgit->second->Set_mode_mess(-3);  // no output from solver at all
    }
#else
    /// add extra functionality for alternative solver if needed
#endif
}


template <uint32_t dim>
void SimulatorControl<dim>::OutputRestartFile(bool manually_triggered)
{
    static size_t restartfilecounter(1);
    /// @todo make an option here, that if the connectivity hasn't changed we do not need to output the connectivity, saving time.
    //static bool output_connectivity_file(false);
    //clock_t currentClockTicks  = clock();
    ScalarVariable current_runtime;
    this->GetSS()->GetModel()->Read(run_time_Key_,current_runtime);

    double runTimeSinceLastRestart=current_runtime() - runtimeOnLastRestart_;

    if ((runTimeSinceLastRestart>=restartFileOuputRunTimeInterval_) || manually_triggered)
    {
        if (restartfilecounter % 2 != 0){
            simulator_setup_->GetModel()->OutputToBinaryFile(string(this->GetName()+"_restartFile0").c_str());
        }
        else if (restartfilecounter % 2 ==0){
            simulator_setup_->GetModel()->OutputToBinaryFile(string(this->GetName()+"_restartFile1").c_str());
        }
        this->runtimeOnLastRestart_=current_runtime();
        restartfilecounter++;
    }
}

template <uint32_t dim>
ComputationalSettings& SimulatorControl<dim>::RunSettings()
{
    return this->simulator_setup_->GetRunSettings();
}

template<uint32_t dim>
void SimulatorControl<dim>::PrintRangeOfVariablesToScreen()
{
    list<SimulatorSetupParameter>::iterator lit;
    list<SimulatorSetupParameter>::iterator litbeg=this->GetSS()->GetParameterList().begin();
    list<SimulatorSetupParameter>::iterator litend=this->GetSS()->GetParameterList().end();
    cout<<"--------------------------------------------------------"<<endl;
    cout<<" Printing property ranges."<<endl;
    cout<<" Input Variables."<<endl;
    double pmin, pmax;
    const PropertyDatabase<dim>& p_ref = this->GetSS()->GetModel()->Database();
    for ( lit = litbeg; lit != litend; lit++){
        string usage=lit->usage;
        std::transform(usage.begin(), usage.end(), usage.begin(), ::tolower);
        if (usage=="input"){
            this->GetSS()->GetModel()->MinMaxOf( (lit->name).c_str(), pmin, pmax );
            cout << scientific << setprecision(5) <<" Range of variable :";
            cout << pmin <<" to "<< pmax <<" name: "<<(lit->name) <<" ["<< p_ref.Unit(lit->name.c_str()) <<"]: '"<<endl;
        }
    }
    cout<<endl;
    cout<<" Computed Variables."<<endl;
    for ( lit = litbeg; lit != litend; lit++){
        string usage=lit->usage;
        std::transform(usage.begin(), usage.end(), usage.begin(), ::tolower);
        if (usage=="computed"){
            this->GetSS()->GetModel()->MinMaxOf( (lit->name).c_str(), pmin, pmax );
            cout << scientific << setprecision(5) <<" Range of variable :";
            cout << pmin <<" to "<< pmax <<" name: "<<(lit->name) <<" ["<< p_ref.Unit(lit->name.c_str()) <<"]: '"<<endl;
        }
    }
}

template <uint32_t dim>
void SimulatorControl<dim>::UpdateCurrentRunTime()
{
#ifdef _OPENMP
    double runtime=(omp_get_wtime()-this->initial_omp_wtime_);
    this->GetSS()->GetModel()->Store(run_time_Key_,makeScalar(PLAIN,runtime));
#else
    unsigned long millisec ((clock() - startClockTicks_) * 1000000 / CLOCKS_PER_SEC);
    double runtime=(double)(millisec)/1000000.;
    //float runtime=float((clock() - startClockTicks_))/float(CLOCKS_PER_SEC);
    //double runtime=double((clock() - startClockTicks_))/double(CLOCKS_PER_SEC);

    this->GetSS()->GetModel()->Store(run_time_Key_,makeScalar(PLAIN,runtime));
#endif
}

template <uint32_t dim>
void SimulatorControl<dim>::Output()
{
    // register the timestep, model time, and run time for possible output.
    this->GetSS()->GetModel()->Store(this->pdb_.StorageKey("timestep"),makeScalar(PLAIN,this->time_increment_));
    this->GetSS()->GetModel()->Store(model_time_Key_,makeScalar(PLAIN,this->GetSimulationTime()));
    this->UpdateCurrentRunTime();

    this->OutputRestartFile(false);

    if (this->RunSettings().IsOutputTime(this->GetSimulationTime(),10e-8))
    {
        if (Verbose()) this->PrintRangeOfVariablesToScreen();
        this->OutputVariableToRegionVTUFiles();
        this->GetSS()->GetSimulatorMonitor()->Monitor();
    } else if (this->RunSettings().IsMonitorTime(this->GetSimulationTime(),10e-8)|| this->MonitorAllTimesteps())
        this->GetSS()->GetSimulatorMonitor()->Monitor();

}

template<uint32_t dim>
void SimulatorControl<dim>::OutputVariableToRegionVTUFiles()
{
    /// @todo put conditional output options here (i.e. OutputTimes)
    string prefix=this->GetName();
    for (list<string>::iterator it = this->GetSS()->GetRegOutList().begin(); it!=this->GetSS()->GetRegOutList().end();it++)
        if (this->GetSS()->GetModel()->ContainsRegion( (*it).c_str()))
            vtu_->OutputDataToVTU(prefix.c_str(),this->GetSS()->GetVTUPropList(), this->GetSS()->GetModel()->Region((*it).c_str()),this->outputVTU_counter_);
        else if (this->GetSS()->GetModel()->ContainsBoundary( static_cast <const string> (*it)))
            vtu_->OutputDataToVTU(prefix.c_str(),this->GetSS()->GetVTUPropList(), this->GetSS()->GetModel()->Boundary(*it),this->outputVTU_counter_);
    this->outputVTU_counter_++;
}

template<uint32_t dim>
double SimulatorControl<dim>::GetCurrentIntervalEndTime(double current_simulation_time)
{
    /** Here we look in the stored time interval start times to see which is the the closest
    to the current simulation time.  If none is found then GetSimulationEndTime() is returned. */

    set<double>  possibleIntervalEndTimes;
    /// check all simulation start times, and choose the closest one.
    for ( auto it=this->GetIntervals().begin();it!=this->GetIntervals().end(); it++ )
        if ( (it->first) - current_simulation_time > 0. )
            possibleIntervalEndTimes.insert( it->first );
    if ( possibleIntervalEndTimes.empty() ){
        return this->GetSimulationEndTime();
    }
    else {
        return *possibleIntervalEndTimes.begin();
    }
}

template<uint32_t dim>
double SimulatorControl<dim>::GetCurrentIntervalStartTime(double current_simulation_time)
{
    /** Here we look in the stored time interval start times to see which is the the closest
    previous to the current simulation time.
    If none is found then GetSimulationStartTime() is returned. */

    set<double>  possibleIntervalStartTimes;
    for ( auto it=this->GetIntervals().begin();it!=this->GetIntervals().end(); it++ )
        if ( current_simulation_time - (it->first)  >= 0. )
            possibleIntervalStartTimes.insert( (it->first) );
    if ( possibleIntervalStartTimes.empty() )
        return this->GetSimulationStartTime();
    else
        return *possibleIntervalStartTimes.rbegin();
}

template<uint32_t dim>
bool SimulatorControl<dim>::YesOrNo()
{
    bool answer(false),yesorno(false);
    while (!answer){
        string choice="";
        getline(cin, choice);
        std::transform(choice.begin(), choice.end(), choice.begin(), ::toupper);
        if (choice.substr(0,1)=="Y"){
            yesorno = true;
            answer=true;
        }
        else if (choice.substr(0,1)=="N"){
            yesorno = false;
            answer=true;
        }
        else
            cout<<" Please write either y or n (one character)."<<endl;
    }
    return yesorno;
}

template<uint32_t dim>
bool SimulatorControl<dim>::CheckForDirichletEssentialCondition(csmp::Index key)
{
    // checks that at least one node of the domain has a pressure condition.
    for ( typename vector<csmp::Node<dim>*>::const_iterator nit=this->GetSS()->GetModel()->Region("Model").NodesBegin();
          nit!=this->GetSS()->GetModel()->Region("Model").NodesEnd(); nit++ )
        if ( (*nit)->Status(key) == DIRICH ) return true;

    cout<<" You do not seem to have any DIRICH flag set for this variable. "<<endl;
    cout<<" If your simulator needs it (say, for fluid pressure conditions in incompressible flow) "<<endl;
    cout<<" Consider setting at least one boundary or interior point to a reference value. "<<endl;
    return false;
}

template<uint32_t dim>
bool SimulatorControl<dim>::InteractiveInputNodeValue(csmp::Index key, SUBDOMAIN_PART sub)
{
    cout<<" Please type in the coordinates of the point where you would like to set an essential condition "<<endl;

    cout<<" x: ";double x=this->GetRealInput();
    cout<<" y: ";double y=this->GetRealInput();
    cout<<" z: ";double z=this->GetRealInput();
    cout<<" Please type the value you would like to assign: ";
    double value=this->GetRealInput();
    cout<<" output :"<<value<<endl;
    this->SetValueNearestToPoint(key,x,y,z,value,DIRICH,sub);
    if (value >= 10e50){
        cout<<" Unreasonable value for assignment. Re-starting the process..."<<endl;
        return false;
    }
    else
        return true;
}


/// This function sets Dirichlet values on points closest to the supplied coordinates.  Note that
/// This will only address those nodes in the interior of the model.
template<uint32_t dim>
void SimulatorControl<dim>::SetReferencePointValuesModelInterior()
{
    for (auto it = this->GetReferencePointsAndValues().begin(); it != this->GetReferencePointsAndValues().end();it++){
        double x=it->second.first[0];
        double y=it->second.first[1];
        double z=it->second.first[2];
        double value=it->second.second;
        this->SetValueNearestToPoint(this->pdb_.StorageKey(it->first.c_str()),x,y,z,value,DIRICH,INTERIOR);
    }
}

/// This function sets Dirichlet values on points closest to the supplied coordinates.  Note that
/// This will only address those nodes in the perimeter of the model.
template<uint32_t dim>
void SimulatorControl<dim>::SetReferencePointValuesModelPerimeter()
{
    for (auto it = this->GetReferencePointsAndValues().begin(); it != this->GetReferencePointsAndValues().end();it++){
        double x=it->second.first[0];
        double y=it->second.first[1];
        double z=it->second.first[2];
        double value=it->second.second;
        this->SetValueNearestToPoint(this->pdb_.StorageKey(it->first.c_str()),x,y,z,value,DIRICH,INTERIOR);
    }
}

template<uint32_t dim>
void SimulatorControl<dim>::SetThicknessFactorFromWellRadii()
{
    ScalarVariable wellradius;
    for (vector<string>::iterator wit=this->GetSS()->GetWells().begin();wit!=this->GetSS()->GetWells().end();wit++){

        Region<dim>& wref=this->GetSS()->GetModel()->Region(wit->c_str());

        if (   !((!containsVolumeElements(wref) &&
                  containsSurfaceElements(wref) &&
                  !containsLineElements(wref))
                 || (!containsVolumeElements(wref) &&
                     !containsSurfaceElements(wref) &&
                     containsLineElements(wref))
                 || (containsVolumeElements(wref) &&
                     !containsSurfaceElements(wref) &&
                     !containsLineElements(wref))) )
            error_handler_.Note( ERROR, "SimulatorControl::SetThicknessFactorFromWellRadii()",
                                   "Wells should contain only elements of the same dimension (all volumes, all surfaces, or all lines)","Check your geometry." );

        if (containsLineElements(wref)){
            wref.Read(this->wellRadius_,wellradius);
            for (auto eit=wref.CellsBegin();eit!=wref.CellsEnd();eit++)
                (*eit)->Store(this->thickness_,makeScalar(PLAIN,wellradius()*wellradius()*PI));
        }
    }
}

template <uint32_t dim>
void SimulatorControl<dim>::InteractivelyManageEssentialConditions()
{
    map<int32_t,string> location_choice;
    location_choice[1]="boundary";
    location_choice[2]="point";
    location_choice[3]="region";
    int32_t choice(0);
    bool finished(false);
    while (!finished){
        cout<<"----------------------------------------------------------------------------------------"<<endl;
        cout<<" The following options will allow you set essential conditions on boundaries (as BCs), "<<endl;
        cout<<" regions, and point locations (based on closest node)."<<endl;
        cout<<" Please select where you would like to set essential conditions."<<endl;
        cout<<" 1. A boundary."<<endl;
        cout<<" 2. A point (closest node to it)"<<endl;
        cout<<" 3. A Region (includes wells)."<<endl;
        cout<<" 4. Quit essential conditions setup."<<endl;
        cout<<" Your choice: ";

        choice=this->GetIntegerInput();
        cout<<" You chose: '"<<choice<<"'"<<endl;
        if (choice==4){
            finished=true;
            break;
        }


        if (location_choice.find(choice)!=location_choice.end()){
            cout<<" You have chosen to set essential conditions on a "<<location_choice.at(choice)<<"."<<endl;
            cout<<"----------------------------------------------------------------------------------------"<<endl;

            if (location_choice.at(choice)=="boundary"){
                cout<<" The following is a list of available boundaries and their flags:"<<endl;
                size_t counter(1);
                map<size_t,string> bnumbertoname;
                for( typename Model<dim>::boundaryConstIterator bit = this->GetSS()->GetModel()->BoundariesBegin() ; bit != this->GetSS()->GetModel()->BoundariesEnd(); ++bit ){
                    string bname=bit->first;
                    cout<<" "<<counter<<". "<<bname<<endl;
                    bnumbertoname[counter]=bname;
                    counter++;
                }
                cout<<" Please select a boundary to which you would like to assign an essential condition: ";
                int32_t bn=this->GetIntegerInput();

                string variable="";
                size_t fcounter(0);
                while (!this->pdb_.IsDefined(variable.c_str()) && fcounter < 5){
                    cout<<" Please type the variable name that you would like to specify at this region. "<<endl;
                    variable=this->GetStringInput();
                    if (this->pdb_.IsDefined(variable.c_str())){
                        cout<<" Property name found."<<endl;
                        cout<<" Please type the value you would like to assign: ";
                        double value=this->GetRealInput();
                        VARIABLE_FLAG vflag;
                        cout<<" Would you like to set a Dirichlet flag on this region? (y/n) ";
                        if (this->YesOrNo())
                            vflag=DIRICH;
                        else
                            vflag=PLAIN;

                        this->SetBoundaryValue(variable.c_str(),bnumbertoname[bn],value,vflag);
                    }
                    else {
                        cout<<" Property name not found please try again (hint: read your -variables.txt file!)"<<endl;
                        fcounter++;
                    }


                }
                if (fcounter >= 5)
                    cout<<" You have attempted to write the variable name at least 5 times. Re-starting the process..."<<endl;
            }
            else if (location_choice.at(choice)=="point"){
                size_t fcounter(0);
                string variable="";
                while (!this->pdb_.IsDefined(variable.c_str()) && fcounter < 5){
                    cout<<" Please type the variable name that you would like to specify near this point. "<<endl;
                    variable=this->GetStringInput();
                    fcounter++;
                }
                if (fcounter< 5)
                    this->InteractiveInputNodeValue(this->pdb_.StorageKey(variable.c_str()));
            }
            else if (location_choice.at(choice)=="region"){
                cout<<" The following is a list of available regions and their flags:"<<endl;
                typename std::map<std::string,csmp::Region<dim> >::const_iterator uregbeg=this->GetSS()->GetModel()->UniqueRegionsBegin();
                typename std::map<std::string,csmp::Region<dim> >::const_iterator uregend=this->GetSS()->GetModel()->UniqueRegionsEnd();
                size_t counter(1);
                map<size_t,string> regnumbertoname;
                for ( typename std::map<std::string,csmp::Region<dim> >::const_iterator rit=uregbeg;rit!=uregend;rit++)
                {
                    string regname=rit->first;
                    cout<<" "<<counter<<". "<<regname<<endl;
                    regnumbertoname[counter]=regname;
                    counter++;
                }
                cout<<" Please select a region to which you would like to assign an essential condition: ";
                int32_t regn=this->GetIntegerInput();

                string variable="";
                size_t fcounter(0);
                while (!this->pdb_.IsDefined(variable.c_str()) && fcounter < 5){
                    cout<<" Please type the variable name that you would like to specify at this region. "<<endl;
                    variable=this->GetStringInput();
                    if (this->pdb_.IsDefined(variable.c_str())){
                        cout<<" Property name found."<<endl;
                        cout<<" Please type the value you would like to assign: ";
                        double value=this->GetRealInput();
                        VARIABLE_FLAG vflag;
                        cout<<" Would you like to set a Dirichlet flag on this region? (y/n) ";
                        if (this->YesOrNo())
                            vflag=DIRICH;
                        else
                            vflag=PLAIN;

                        this->SetRegionValue(variable.c_str(),regnumbertoname[regn],value,vflag);
                    }
                    else {
                        cout<<" Property name not found please try again (hint: read your -variables.txt file!)"<<endl;
                        fcounter++;
                    }


                }
                if (fcounter >= 5)
                    cout<<" You have attempted to write the variable name at least 5 times. Re-starting the process..."<<endl;

            }
        }
    }
}


//template <uint32_t dim>
//void SimulatorControl<dim>::ManageTimeIntervals()
//{
//    cout<<"-----------------------------------------------------------"<<endl;
//    cout<<" These options will allow you divide the duration into time intervals."<<endl;
//    cout<<" Your simulation will simulate "<<this->GetSS()->GetRunSettings().Duration()<<" [s]"<<endl;
//    cout<<" You may now enter the number of intervals in which to subdivide this time."<<endl;
//    cout<<" Note that a value of below 1 will run a single interval."<<endl;
//    cout<<" Number of intervals?: ";
//    cout.flush();
//    std::string choice="";
//    getline(cin, choice);
//    cout<<endl;
//    size_t nintervals=stringToNumber<uint32_t>(choice);
//    cout<<" You have chosen to establish "<<nintervals<<" time intervals."<<endl;
//    cout.flush();

//    if (nintervals>1){
//        bool checked_intervals(false);
//        while (!checked_intervals) {
//            cout<<"-----------------------------------------------------------"<<endl;
//            cout<<" Please enter the interval times. The order is irrelevant."<<endl;
//            cout<<" Values must be below the simulation duration: "<<this->GetSS()->GetRunSettings().Duration()<<" [s]"<<endl;
//            this->interval_start_times_.clear();
//            this->interval_start_times_.insert(make_pair(this->GetSS()->GetRunSettings().Duration()));
//            for (auto i = 0; i < nintervals-1;i++)
//            {
//                cout<<" End time for interval "<<i+1<<"?: ";
//                cout.flush();
//                choice="";
//                getline(cin, choice);
//                cout<<endl;
//                double endtime=stringToNumber<double>(choice);
//                if (endtime < this->GetSS()->GetRunSettings().Duration())
//                    this->interval_start_times_.insert(endtime);
//                else
//                    cout<<" interval endtime is larger than the duration settings on your -configuration file."<<endl;
//            }
//            checked_intervals=this->CheckTimeIntervals();
//        }
//    }



//}

template <uint32_t dim>
void SimulatorControl<dim>::TrimTimeIncrementWithOutputTime()
{
    this->time_increment_=std::min(this->time_increment_,this->RunSettings().NearestOutputTime(this->GetSimulationTime())-this->GetSimulationTime());
}

template <uint32_t dim>
void SimulatorControl<dim>::TrimTimeIncrementWithMonitorTime()
{
    if (!this->monitor_all_timesteps_)
        this->time_increment_=std::min(this->time_increment_,this->RunSettings().NearestMonitorTime(this->GetSimulationTime())-this->GetSimulationTime());
}

template <uint32_t dim>
void SimulatorControl<dim>::TrimTimeIncrementWithEndOfInterval()
{
    this->time_increment_=this->GetIntervalEndTime()-this->GetSimulationTime();
}

template <uint32_t dim>
void SimulatorControl<dim>::ManageWellRates()
{
    cout<<" Would you like to (re)set well rates? (y/n): ";
    bool setupwells=this->YesOrNo();
    if (setupwells){
        bool finished(false);
        while (!finished) {
            map<size_t,string> wellnumbertoname;
            cout<<" Note: Each well is associated to a particular flow rate variable."<<endl;
            cout<<" The following is a list of wells, their current state, and associated flow rate variable:"<<endl;
            size_t counter(1);
            double well_rate_max(0),well_rate_min(0);

            for (vector<string>::iterator lit= this->GetSS()->GetWells().begin();lit!=this->GetSS()->GetWells().end();lit++){
                this->GetSS()->GetModel()->Region( (*lit).c_str() ).MinMaxOf(this->GetSS()->GetWellRateVars().at(*lit).first.c_str() , well_rate_min, well_rate_max );
                cout<<" "<<counter<<". "<<*lit<<" rate var: '"<<this->GetSS()->GetWellRateVars().at(*lit).first<<"' set to : "<<well_rate_min<<endl;
                wellnumbertoname[counter]=*lit;
                counter++;
            }
            cout<<" "<<counter<<". Quit well rate manager."<<endl;
            int32_t regn=this->GetSS()->GetWells().size()+2;
            while ( regn < 1 || regn > this->GetSS()->GetWells().size() ){
                cout<<" Please select a well to which you would like to assign a rate, or "<<this->GetSS()->GetWells().size()+1<<" to quit."<<endl;
                regn=GetIntegerInput();
                if (regn == this->GetSS()->GetWells().size()+1 ){
                    finished=true;
                    break;
                }
                if (regn >=1 && regn <= this->GetSS()->GetWells().size()){
                    cout<<" Please type the rate value you would like to assign to well '"<<wellnumbertoname[regn]<<"' :";
                    double value=this->GetRealInput();
                    this->SetRegionValue((this->GetSS()->GetWellRateVars().at(wellnumbertoname[regn]).first),wellnumbertoname[regn],value,PLAIN);
                }
            }
        }
    }
}

template <uint32_t dim>
double SimulatorControl<dim>::GetRealInput(bool success)
{
    std::string choice="";
    getline(cin, choice);
    if (!choice.empty()){
        return std::stod(choice);
        success=true;
    }
    else{
        return 10e50;
        success=false;
    }
}

template <uint32_t dim>
int32_t SimulatorControl<dim>::GetIntegerInput()
{
    std::string choice="";
    getline(cin, choice);
    if (!choice.empty())
        return std::stoi(choice);
    else
        return 0;
}

template <uint32_t dim>
string SimulatorControl<dim>::GetStringInput()
{
    std::string choice="";
    getline(cin, choice);
    return choice;
}

template <uint32_t dim>
void SimulatorControl<dim>::SetBoundaryFlag(const char* property_name, string boundary,VARIABLE_FLAG flag)
{
    if (this->GetSS()->GetModel()->ContainsBoundary(boundary))
        this->GetSS()->GetModel()->Boundary(boundary).ChangePropertyStatus(property_name, flag);
    else
        error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::SetBoundaryFlag()","Model does not contain boundary: ",boundary.c_str() );
}

template <uint32_t dim>
void SimulatorControl<dim>::SetPropertyFlagToTopBoundaries(const char* property_name, VARIABLE_FLAG flag)
{
    for (auto itb=this->GetSS()->GetTopBoundaries().begin();itb!=this->GetSS()->GetTopBoundaries().end();itb++){
        this->SetBoundaryFlag(property_name,*itb,flag);
    }
}

template <uint32_t dim>
void SimulatorControl<dim>::SetBoundaryValue(const char* prop_name,string boundary,double value,VARIABLE_FLAG flag)
{
    csmp::Index key=this->GetSS()->GetModel()->Database().StorageKey(prop_name);
    if (key.type==SCALAR)
        this->GetSS()->GetModel()->Boundary(boundary.c_str()).InputPropertyValue(prop_name, makeScalar(flag,value));
    else
        error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::SetBoundaryValue()","Only scalars may be set on the boundary","(for now)" );
}

template <uint32_t dim>
void SimulatorControl<dim>::SetRegionValue(string prop_name,string region,double value,VARIABLE_FLAG flag)
{
    csmp::Index key=this->GetSS()->GetModel()->Database().StorageKey(prop_name.c_str());
    if (key.type==SCALAR)
        this->GetSS()->GetModel()->Region(region.c_str()).InputPropertyValue(prop_name.c_str(), makeScalar(flag,value));
    else
        error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::SetRegionValue()","Only scalars may be set","(for now)" );
}

template <uint32_t dim>
void SimulatorControl<dim>::SetRegionFlag(const char* prop_name,string region,VARIABLE_FLAG flag)
{
    if (this->GetSS()->GetModel()->ContainsRegion(region.c_str()))
        this->GetSS()->GetModel()->Region(region.c_str()).ChangePropertyStatus(prop_name, flag);
    else
        error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::SetRegionFlag()","Model does not contain region: ",region.c_str() );
}


template <uint32_t dim>
void SimulatorControl<dim>::SetFlagNearestToPoint(Index key,double x, double y, double z, VARIABLE_FLAG flag)
{

}

template <uint32_t dim>
void SimulatorControl<dim>::SetValueNearestToPoint(Index key,double x, double y, double z, double value,VARIABLE_FLAG flag,SUBDOMAIN_PART sub)
{
    auto nodes_begin=this->GetSS()->GetModel()->Region("Model").NodesBegin();
    typename vector<Node<dim>*>::const_iterator nodes_end;
    if (sub==COMPLETE)
        nodes_end=this->GetSS()->GetModel()->Region("Model").NodesEnd();
    else if (sub==INTERIOR)
        nodes_end=this->GetSS()->GetModel()->Region("Model").PerimeterNodesBegin();
    double mind(10e50),distance(0.);

    Node<dim>* n(nullptr);
    for ( auto npit= nodes_begin; npit!=nodes_end;npit++) {
        distance=pow((pow(((*npit)->x()-x),2.0)+pow(((*npit)->y()-y),2.0)+pow(((*npit)->z()-z),2.0)),0.5);
        if (distance<mind) {
            n=*npit;
            mind=distance;
        }
    }
  
    assert( n != nullptr );
    n->Store(key,makeScalar(flag,value));
}



template <uint32_t dim>
void SimulatorControl<dim>::BeginInterval()
{
    this->SetIntervalStartTime(this->GetCurrentIntervalStartTime(this->GetSimulationTime()));
    this->SetIntervalEndTime(this->GetCurrentIntervalEndTime(this->GetSimulationTime()));
    this->GetSS()->GetModel()->Store(pdb_.StorageKey("interval"),makeScalar(PLAIN,this->GetIntervalNumber(this->GetSimulationTime())));
}

template <uint32_t dim>
bool SimulatorControl<dim>::AskUserToContinue()
{
    cout<<"  Do you want to continue the simulation for another time interval? (y/n):";
    return this->YesOrNo();
}

template <uint32_t dim>
bool SimulatorControl<dim>::CheckTimeIntervals()
{
    cout<<" Simulation Start Time    [s]: "<<this->GetSimulationStartTime()<<endl;
    cout<<" Simulation End Time      [s]: "<<this->GetSimulationEndTime()<<endl;

    cout<<" Checking time interval settings...."<<endl;
    cout<<" Number of time intervals    : "<<this->GetIntervals().size()<<endl;
    cout<<" You have entered the following time intervals "<<endl;

    size_t counter(0);
    double fake_simulation_time(this->GetSimulationStartTime());
    for (auto it = this->GetIntervals().begin(); it != this->GetIntervals().end();it++){
        cout<<"-------------------------------------------------------------------------------------"<<endl;
        cout<<" Interval "<<counter+1<<" named: "<<it->second<<" from :"<<this->GetCurrentIntervalStartTime(fake_simulation_time)
           <<" [s] to: "<<this->GetCurrentIntervalEndTime(fake_simulation_time)<<" [s]. "<<endl;
        fake_simulation_time=GetCurrentIntervalEndTime(fake_simulation_time);
        counter++;
    }
    cout<<"-------------------------------------------------------------------------------------"<<endl;
    if (this->GetSS()->Interactive()){
        cout<<" Is this information correct? (y/n): ";
        return this->YesOrNo();
    }
    else
        return true;

}

template <uint32_t dim>
void SimulatorControl<dim>::CatchSignals()
{
    SimulationSignalHandler sig;
    this->UpdateCurrentRunTime();

    if (sig.OutputVTUSignal())
    {
        simulator_setup_->GetModel()->Store( model_time_Key_, makeScalar(PLAIN,this->GetSimulationTime()));

        string prefix=this->GetName();
        vtu_->SetSuffixText("_ManualTrigger");
        for (list<string>::iterator it = this->GetSS()->GetRegOutList().begin(); it!=this->GetSS()->GetRegOutList().end();it++)
            if (this->GetSS()->GetModel()->ContainsRegion( (*it).c_str()))
                vtu_->OutputDataToVTU(prefix.c_str(),this->GetSS()->GetVTUPropList(), this->GetSS()->GetModel()->Region((*it).c_str()),this->outputVTU_counter_);
            else if (this->GetSS()->GetModel()->ContainsBoundary( static_cast <const string> (*it)))
                vtu_->OutputDataToVTU(prefix.c_str(),this->GetSS()->GetVTUPropList(), this->GetSS()->GetModel()->Boundary(*it),this->outputVTU_counter_);
        vtu_->SetSuffixText("");
        this->outputVTU_counter_++;

        sig.OutputVTUSignal(false);
    }

    if (sig.OutputMonitorDataSignal()) {
        this->GetSS()->GetSimulatorMonitor()->Monitor();
        sig.OutputMonitorDataSignal(false);
    }

    if (sig.OutputRestartFileSignal()) {
        this->OutputRestartFile(true);
        sig.OutputRestartFileSignal(false);
    }

    if (sig.QuitSignal()){
        error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::CatchSignals()","User has chosen to quit the simulation.","Shutting down SimulatorControl." );
        //exit(1);
    }

}

template <uint32_t dim>
void SimulatorControl<dim>::SetSimulationTime(double t)
{
    this->simulationTime_=t;
}

/// returns run-time, in seconds, since an initial reference point
/// normally set by the BeginInterval() method.
template <uint32_t dim>
double SimulatorControl<dim>::CalculateCurrentRunTime()
{

#ifdef _OPENMP
    return (omp_get_wtime()-this->initial_omp_wtime_);
#else
    clock_t actual_runtime=clock();
    return double(actual_runtime-this->runTimeStart_)/CLOCKS_PER_SEC;
#endif

}

template<uint32_t dim>
void SimulatorControl<dim>::SyncOutputAntMonitoringTimesToModel()
{
    cout<<" Syncing Output Times (monitor and vtu frames) to Model..."<<endl;
    if (simulator_setup_->GetModel()->Database().IsDefined("vtu frames"))
    {

        if (!isnan(simulator_setup_->GetModel()->Read( vtu_frames_Key_ ) ) &&
                simulator_setup_->GetModel()->Read( vtu_frames_Key_ )>0.0 )
        {
            this->SetOutputTimes(simulator_setup_->GetModel()->Read( vtu_frames_Key_ ), this->GetSimulationStartTime());
            //            if (double(simulator_setup_->GetModel()->Database().Components("output times"))>simulator_setup_->GetModel()->Read( vtu_frames_Key_) )
            //                this->SetOutputTimes(simulator_setup_->GetModel()->Read( vtu_frames_Key_ ), this->GetSimulationStartTime());
            //            else
            //                throw csmp::Exception (ERROR,"SimulatorControl<dim>::SyncOutputTimesToModel()"
            //                                       ,"vtu frames cannot be larger than the amount of output times to save. Decrease vtu frames or increase the size of output times in the variables file.");
        }
        else{
            cout<<" SimulatorControl<dim>::SyncOutputTimesToModel() Warning.  No number of 'vtu frames' set, using output times from the configuration file or whatever was set in the restart file."<<endl;
            this->SetOutputTimes(1, this->GetSimulationStartTime());
        }

    }
    else
        throw csmp::Exception (ERROR,"SimulatorControl<dim>::SyncOutputTimesToModel()"
                               ,"'vtu frames' needs to be defined variables. Set them to zero to use pre-defined output times from your configuration file as monitor and vtu output times.");

    if (simulator_setup_->GetModel()->Database().IsDefined("monitor frames"))
    {

        if (!isnan(simulator_setup_->GetModel()->Read( monitor_frames_Key_ ) ) &&
                simulator_setup_->GetModel()->Read( monitor_frames_Key_ )>0.0 )
        {
            if (!this->restart_)
                this->SetMonitorTimes(simulator_setup_->GetModel()->Read( monitor_frames_Key_ ),this->GetSimulationStartTime());

        }
        else{
            cout<<" SimulatorControl<dim>::SyncOutputTimesToModel() Warning.  No number of 'monitor frames' set, using output times from the configuration file or whatever was set in the restart file."<<endl;
            cout<<" Please press enter to confirm and continue."<<endl;
            cin.get();
        }

    }
    else
        throw csmp::Exception (ERROR,"SimulatorControl<dim>::SyncOutputTimesToModel()"
                               ,"'monitor frames' needs to be a defined variable. Set them to zero to use pre-defined output times from your configuration file as monitor and vtu output times.");


    ArrayVariable av;

    if (simulator_setup_->GetModel()->Database().IsDefined("output times")){
        Index avkey=this->GetSS()->GetModel()->Database().StorageKey("output times");
        this->GetSS()->GetModel()->Read(avkey,av);
        cout<<" Reading "<<av.Size()<<" Output times..."<<endl;
        for (auto i = 0 ; i < av.Size();i++)
            if (!isnan(av(i))){
                cout<<" "<<i<<" time: "<<av(i)<<endl;
                this->RunSettings().AddOutputTime(av(i));
            }
        this->GetSS()->GetModel()->DeleteProperty("output times");
        cout<<" Done."<<endl;
    }


    if (simulator_setup_->GetModel()->Database().IsDefined("monitor times")){
        Index avkey=this->GetSS()->GetModel()->Database().StorageKey("monitor times");
        this->GetSS()->GetModel()->Read(avkey,av);
        cout<<" Reading "<<av.Size()<<" Monitoring times..."<<endl;
        for (auto i = 0 ; i < av.Size();i++)
            if (!isnan(av(i)))
                this->RunSettings().AddMonitorTime(av(i));
        this->GetSS()->GetModel()->DeleteProperty("monitor times");

        cout<<" Done."<<endl;
    }

    size_t nout = std::distance(this->RunSettings().OutputTimesBegin(),this->RunSettings().OutputTimesEnd());
    size_t nmon = std::distance(this->RunSettings().MonitorTimesBegin(),this->RunSettings().MonitorTimesEnd());

    SimulatorSetupParameter sp;
    sp.name="output times";        sp.notation="OT";   sp.unit="s";    sp.type=ARRAY; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="computed";sp.vsize=nout;
    this->GetSS()->GetModel()->CreateProperty(sp.name.c_str(),sp.notation.c_str(),sp.unit.c_str(),sp.type,sp.placement,sp.vsize,sp.min,sp.max,sp.usage);
    sp.key=this->GetSS()->GetModel()->Database().StorageKey(sp.name.c_str());
    this->GetSS()->GetParameterList().push_back(sp);
    sp.name="monitor times";       sp.notation="MONT"; sp.unit="s";    sp.type=ARRAY; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="computed";sp.vsize=nmon;
    this->GetSS()->GetModel()->CreateProperty(sp.name.c_str(),sp.notation.c_str(),sp.unit.c_str(),sp.type,sp.placement,sp.vsize,sp.min,sp.max,sp.usage);
    sp.key=this->GetSS()->GetModel()->Database().StorageKey(sp.name.c_str());
    this->GetSS()->GetParameterList().push_back(sp);

    cout<<" Done Setting Output Frequency."<<endl;
}

template<uint32_t dim>
void SimulatorControl<dim>::SetOutputTimes(size_t n, double starttime)
{
    /** @todo double check if it isn't convenient to set the output intervals related
      to start time (i.e. (duration-starttime)/n).  At the moment, this does not seem so. --Julian
    */
    double outputtimeinterval=this->RunSettings().Duration()/static_cast<double>(n);

    cout<<" Generating "<<n <<" output times.";
    set<double> ots;
    ots.insert(starttime);

    double outputtime(starttime);
    for (auto i = 0 ; i < n ;i++)
    {
        outputtime+=outputtimeinterval;
        ots.insert(outputtime);
    }
    this->RunSettings().SetOutputTimes(ots,false); // this will not overwrite existing output times.
}

template<uint32_t dim>
void SimulatorControl<dim>::SetMonitorTimes(size_t n, double starttime)
{
    double outputtimeinterval=this->RunSettings().Duration()/static_cast<double>(n);

    cout<<" Generating "<<n <<" monitor times.";
    set<double> ots;
    ots.insert(starttime);

    double outputtime(starttime);
    for (auto i = 0 ; i < n ;i++)
    {
        outputtime+=outputtimeinterval;
        ots.insert(outputtime);
    }
    this->RunSettings().SetMonitorTimes(ots,false); // this will not overwrite existing monitor times.
}


template<uint32_t dim>
void SimulatorControl<dim>::LoadMonitoringTimesFromRestartedModel()
{

    set<double> ots;
    ArrayVariable mtimes;
    this->GetSS()->GetModel()->Read(monitor_frames_Key_,mtimes);
    cout<<" Loading "<<mtimes.Size() <<" monitor times from restarted model.";
    for (auto i = 0 ; i < mtimes.Size() ;i++)
        ots.insert(mtimes(i));
    this->RunSettings().SetMonitorTimes(ots,true); // this will not overwrite existing monitor times.
}

template <uint32_t dim>
bool SimulatorControl<dim>::ReadControlOptions()
{
    std::string  text_line;
    std::string token;
    std::ifstream ifs( string(this->GetSS()->GetProjectName()+"-control.txt").c_str() );

    if ((isInputFileEmpty(ifs) || !ifs.is_open()) && !this->restart_) {
        cout<<" ******************** ATTENTION **********************************"<<endl;
        cout<<" It seems that you are missing a file (or it is empty) named: "<<this->GetSS()->GetProjectName()<<"-control.txt"<<endl;
        cout<<" This file can be generated taking into account general and simulator-specific options. "<<endl;
        cout<<" NOTE: Tensors and Array variables will be skipped by this process. "<<endl;
        cout<<"       All unnecessary lines can be deleted manually if desired."<<endl;
        cout<<" \nWould you like to generate a control file? [y/n]: "<<endl;
        if (this->YesOrNo())
            this->OutputSampleControlFile();
        else
            error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::OutputSampleControlFile", " A sample control file is needed by the simulator",
                                   " Please create it (or re-run the simulator so that it will create a sample one for you) and restart the simulator.");
    }

    std::getline( ifs, text_line );
    if( ifs.eof())
        return true;
    cout<<" Reading control options from *-control.txt file."<<endl;

    while ((!ifs.eof()))
    {
        if( text_line[0] != '#' && (text_line.length() >= 1))
        {
            vector<string> listoftokens;
            std::stringstream iss(text_line);

            while (getline(iss,token,'\t'))
                listoftokens.push_back(token);

            if (this->Verbose()) cout<<" SimulatorControl<dim> Reading control option: '"<<listoftokens[0]<<"'"<<endl;

            string tok1=listoftokens[0];
            std::transform(tok1.begin(), tok1.end(), tok1.begin(), ::toupper);


            // get the simulation endtime
            if (tok1.find("END TIME")!=std::string::npos) {
                if (listoftokens.size()<3)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for simulation duration"," set 'end time' in the -control.txt file.");
                double time_value(0);
                token=listoftokens[2];
                std::transform(token.begin(), token.end(), token.begin(), ::tolower);
                std::string time_unit(token.substr(0,1));
                if (time_unit=="s")
                    time_value=std::stod(listoftokens[1]);
                else if (time_unit=="h")
                    time_value=3600.0 * std::stod(listoftokens[1]);
                else if (time_unit=="d")
                    time_value=86400.0 * std::stod(listoftokens[1]);
                else if (time_unit=="y")
                    time_value=86400.0*365.0 * std::stod(listoftokens[1]);
                else
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," End time value or unit not recognized.","Should be one of seconds, hours, days, years");

                this->GetSS()->RunSettings().Duration(time_value);
                this->SetSimulationEndTime(time_value);
                this->RunSettings().AddOutputTime(time_value);
                this->RunSettings().AddMonitorTime(time_value);
            }

            if (tok1.find("VTU FRAMES")!=std::string::npos) {
                if (listoftokens.size()<2)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for vtu frames"," set the number in the -control.txt file.");
                size_t value(0);
                value=std::stoul(listoftokens[1]);
                // We now modify the model variables accordingly, to fit the current control options,
                // only if restart is not active.
                if (!this->GetSS()->Restart())
                    this->GetSS()->GetModel()->Store(this->vtu_frames_Key_,makeScalar(PLAIN,value));

            }

            if (tok1.find("MONITOR FRAMES")!=std::string::npos) {
                if (listoftokens.size()<2)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for monitor frames");
                size_t value(0);
                value=std::stoul(listoftokens[1]);
                // We now modify the model variables accordingly, to fit the current control options,
                // only if restart is not active.
                if (!this->GetSS()->Restart())
                    this->GetSS()->GetModel()->Store(this->monitor_frames_Key_,makeScalar(PLAIN,value));
            }

            if (tok1.find("OUTPUT TIMES")!=std::string::npos) {
                if (listoftokens.size()<3)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for output times."," Please provide a unit type and at least one output time.");

                double time_value(0);
                token=listoftokens[1];
                std::transform(token.begin(), token.end(), token.begin(), ::tolower);
                std::string time_unit(token.substr(0,1));

                /// Note: These times will output monitoring as well!
                for  (auto i = 2 ; i < listoftokens.size();i++){
                    if (time_unit=="s")
                        time_value=std::stod(listoftokens[i]);
                    else if (time_unit=="h")
                        time_value=3600.0 * std::stod(listoftokens[i]);
                    else if (time_unit=="d")
                        time_value=86400.0 * std::stod(listoftokens[i]);
                    else if (time_unit=="y")
                        time_value=86400.0*365.0* std::stod(listoftokens[i]);
                    else
                        error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Output time value or unit not recognized."," Should be one of seconds, hours, days, years");

                    this->RunSettings().AddOutputTime(time_value);
                    this->RunSettings().AddMonitorTime(time_value);
                }
            }

            if (tok1.find("MONITOR TIMES")!=std::string::npos) {
                if (listoftokens.size()<3)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for monitor times."," Please provide a unit type and at least one output time.");

                double time_value(0);
                token=listoftokens[1];
                std::transform(token.begin(), token.end(), token.begin(), ::tolower);
                std::string time_unit(token.substr(0,1));

                for  (auto i = 2 ; i < listoftokens.size();i++){
                    if (time_unit=="s")
                        time_value=std::stod(listoftokens[i]);
                    else if (time_unit=="h")
                        time_value=3600.0*std::stod(listoftokens[i]);
                    else if (time_unit=="d")
                        time_value=86400.0*std::stod(listoftokens[i]);
                    else if (time_unit=="y")
                        time_value=86400.0*365.0*std::stod(listoftokens[i]);
                    else
                        error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Output time value or unit not recognized."," Should be one of seconds, hours, days, years");

                    this->RunSettings().AddMonitorTime(time_value);
                }
            }


            cout<<" TOKEN :  "<<tok1<<endl;
            if (tok1.find("RESTART FILE")!=std::wstring::npos) {
                if (listoftokens.size()<2)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for restart file output runtime interval."," Please provide a unit type and at least one output time.");
                //                cout<<"reading restart options: "<<listoftokens[0]<<" "<<listoftokens[1]<<endl;
                double value(0);
                value=std::stod(listoftokens[1]);
                //                cout<<" value: "<<value<<endl;
                this->restartFileOuputRunTimeInterval_=value;
                //                cin.get();
            }

            // Get interval start times.
            if (tok1.find("INTERVAL NAME")!=std::string::npos) {
                if (listoftokens.size()<3)
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()"," Not enough parameters provided (or read) for simulation time interval"," check the interval settings of your -control.txt file.");

                double time_value(0);
                token=listoftokens[3];
                std::transform(token.begin(), token.end(), token.begin(), ::tolower);
                string time_unit(token.substr(0,1));
                if (time_unit=="s")
                    time_value=std::stod(listoftokens[2]);
                else if (time_unit=="h")
                    time_value=3600.0*std::stod(listoftokens[2]);
                else if (time_unit=="d")
                    time_value=86400.0*std::stod(listoftokens[2]);
                else if (time_unit=="y")
                    time_value=86400.0*365.0*std::stod(listoftokens[2]);
                else{
                    string errmsg="Interval start time value or unit not recognized: '"+time_unit+"' text_line: '"+text_line+"'";
                    error_handler_.Note(EXCEPTION,"SimulatorControl<dim>::ReadControlOptions()",errmsg.c_str()," Should be one of seconds, hours, days, years");
                }
                this->InsertNewTimeInterval(listoftokens[1],time_value);
                this->RunSettings().AddOutputTime(time_value);
                this->RunSettings().AddMonitorTime(time_value);
            }

            // Read VTU output variables for the full model.
            if (tok1.find("VTU OUTPUT")!=std::string::npos) {
                string vname;
                for  (auto i = 1 ; i < listoftokens.size();i++){
                    vname=listoftokens[i];
                    this->GetSS()->GetVTUPropList().push_back(vname);
                }
            }

            // Read Monitor output variables for the full model.
            if (tok1.find("MONITOR RANGE")!=std::string::npos) {
                string vname;
                for  (auto i = 1 ; i < listoftokens.size();i++){
                    vname=listoftokens[i];
                    if (this->Verbose()) cout<<" monitored range property names: "<<vname;
                    this->GetSS()->GetMonitoredRangePropList().push_back(vname);
                }
                cout<<endl;
            }

            if (tok1.find("MONITOR INTEGRAL")!=std::string::npos) {
                string vname;
                for  (auto i = 1 ; i < listoftokens.size();i++){
                    vname=listoftokens[i];
                    if (this->Verbose()) cout<<" monitored integral property names: "<<vname;
                    this->GetSS()->GetMonitoredIntegralPropList().push_back(vname);
                }
                cout<<endl;
            }

            if (tok1.find("MONITOR VALUE")!=std::string::npos) {
                string vname;
                for  (auto i = 1 ; i < listoftokens.size();i++){
                    vname=listoftokens[i];
                    if (this->Verbose()) {cout<<" monitored value property names: "<<vname;}
                    if (this->GetSS()->GetModel()->Database().IsDefined(vname.c_str())){
                        if (this->GetSS()->GetModel()->Database().StorageKey(vname.c_str()).place == MODEL){
                            this->GetSS()->GetMonitoredValuePropList().push_back(vname);
                        }
                        else{
                            error_handler_.Note(FATAL_ERROR,"SimulatorControl<dim>::ReadControlOptions()"," Monitored singular values need to be variables created on the MODEL.",
                                                  " Double check your code and/or your CreateParameterList method of your derived SimulatorSetup class ");
                        }
                    }

                }
                cout<<endl;
            }

            if (tok1.find("SUBDOMAIN RANGES")!=std::string::npos
                    || tok1.find("SUBDOMAIN INTEGRALS")!=std::string::npos
                    || tok1.find("SUBDOMAIN DIMENSION")!=std::string::npos
                    || tok1.find("SUBDOMAIN PERIMETER")!=std::string::npos){
                set<string> possible_regions;
                for ( auto git=this->GetSS()->GetModel()->UniqueRegionsBegin(); git!=this->GetSS()->GetModel()->UniqueRegionsEnd(); git++ )
                    possible_regions.insert((*git).first);
                for ( auto git=this->GetSS()->GetModel()->RegionsBegin(); git!=this->GetSS()->GetModel()->RegionsEnd(); git++ )
                    possible_regions.insert((*git).first);
                set<string> possible_boundaries;
                for ( auto bit=this->GetSS()->GetModel()->BoundariesBegin(); bit!=this->GetSS()->GetModel()->BoundariesEnd(); bit++ )
                    possible_boundaries.insert((*bit).first);
                for ( auto bit=this->GetSS()->GetModel()->BoundariesBegin(); bit!=this->GetSS()->GetModel()->BoundariesEnd(); bit++ )
                    possible_boundaries.insert((*bit).first);

                string subdomain_name;
                for  (auto i = 1 ; i < listoftokens.size();i++){
                    subdomain_name=listoftokens[i];

                    if ( possible_regions.find(subdomain_name) != possible_regions.end()
                         || possible_boundaries.find(subdomain_name) != possible_boundaries.end()){
                        if (this->Verbose()) cout<<" "<<tok1<<" subdomain names: "<<subdomain_name;
                        if ( tok1.find("SUBDOMAIN RANGES")!=std::string::npos)
                            this->GetSS()->GetMonitoredRegionRangeList().push_back(subdomain_name);
                        if ( tok1.find("SUBDOMAIN INTEGRALS")!=std::string::npos)
                            this->GetSS()->GetMonitoredRegionIntegralList().push_back(subdomain_name);
                        if ( tok1.find("SUBDOMAIN DIMENSION")!=std::string::npos)
                            this->GetSS()->GetMonitoredRegionDimensionList().push_back(subdomain_name);
                        if ( tok1.find("SUBDOMAIN PERIMETER")!=std::string::npos)
                            this->GetSS()->GetMonitoredRegionPerimeterList().push_back(subdomain_name);
                    }
                }
                cout<<endl;
            }


            if (tok1.find("MONITOR ALL TIMESTEPS")!=std::string::npos)
                this->MonitorAllTimesteps(true);
        }
        std::getline( ifs, text_line  );
    }
    ifs.close();

    cout<<" Finished reading the control parameters for this simulator. "<<endl;
    cout<<"-------------------------------------------------------------------------------------"<<endl;

    // Now remove possible duplicate names in the monitored variables and vtu properties list (respectively)
    this->GetSS()->GetVTUPropList().sort();
    this->GetSS()->GetMonitoredPropList().sort();

    this->GetSS()->GetVTUPropList().unique();
    this->GetSS()->GetMonitoredPropList().unique();

    double starttime(10e32);
    for (auto it = this->GetIntervals().begin(); it != this->GetIntervals().end(); ++it ){
        if (it->first < starttime)
            starttime=it->first;
    }

    this->SetSimulationStartTime(starttime);
    this->CheckTimeIntervals();
    this->SyncOutputAntMonitoringTimesToModel();

    if (restart_)
        this->SetSimulationStartTime(this->simulator_setup_->GetModel()->Read(model_time_Key_));

    return true;
}

template <uint32_t dim>
void SimulatorControl<dim>::InsertNewTimeInterval(string interval_name,double starttime)
{
    if (starttime >= this->GetSimulationTime()){
        this->GetIntervals().push_back(make_pair(starttime,interval_name));
    }
    std::sort(this->GetIntervals().begin(),this->GetIntervals().end());
}

/**

@todo  Hualp! -> use standard library

*/
template <uint32_t dim>
void SimulatorControl<dim>::OutputSimulationTimeToScreen(double time)
{
    auto days    = floor(time/86400.0);
    auto hours   = floor(( time - days*86400 )/3600);
    auto minutes = floor((time-days*86400-hours*3600)/60);
    auto seconds = floor(time - days*86400 - hours*3600 - minutes*60);
    string stime = "[      ] Simulated time: "+std::to_string(days) +" days, "+std::to_string(hours)+" hours, "+std::to_string(minutes)+" minutes, "+std::to_string(seconds)+" seconds                    ";
    cout<<"\r"<<stime<<"\r["<<flush;
}

template<uint32_t dim>
void SimulatorControl<dim>::OutputSampleControlFile()
{
    ofstream fout(this->GetSS()->GetProjectName()+"-control.txt",std::ofstream::out | std::ofstream::app);

    cout<<" Outputting to sample file: "<<this->GetSS()->GetProjectName()<<"-control.txt"<<endl;
    fout<<"# "<<this->GetSS()->GetProjectName()<<" control file. "<<endl;
    fout<<"# Do not worry about empty lines, or the vertical order in which parameters are read."<<endl;
    fout<<endl;  // skip one line.
    fout<<"# Simulation End time information (seconds, hours, days, or years)"<<endl;
    fout<<"end time\t1.0\tdays"<<endl;
    fout<<endl;  // skip one line.
    fout<<"# VTU frames output to disk in this simulation"<<endl;
    fout<<"# (WARNING: This number could limit your timestep if it is too large!)"<<endl;
    fout<<"vtu frames\t5"<<endl;
    fout<<endl;  // skip one line.
    fout<<"# Monitor frames output to disk in this simulation"<<endl;
    fout<<"# (WARNING: This number could limit your timestep if it is too large!)"<<endl;
    fout<<"monitor frames\t10000"<<endl;
    fout<<endl;  // skip one line.
    fout<<"# The following option sets up the computing time interval between which a safety restart file is written."<<endl;
    fout<<"#restart file ouput runtime interval\t3600"<<endl;
    fout<<"# Simulation Intervals (each interval can have a name or number assigned)"<<endl;
    fout<<"# Interval\tinterval name\tstart time"<<endl;
    fout<<"interval\t1\t0.0\tdays"<<endl;
    fout<<"interval\t2\t0.5\tdays"<<endl;
    fout<<"#interval\t3\t0.0\tdays"<<endl;
    fout<<endl;
    fout<<"# Specific output times (both monitor and vtu frames will be output for these times)"<<endl;
    fout<<"# output times\tunit(secs,days,hours,years)\ttime1\ttime2\ttime3...."<<endl;
    fout<<"#output times\tdays\t0.9\t1.5"<<endl;
    fout<<"#output times\tdays\t1.9\t2.9"<<endl;
    fout<<endl;
    fout<<"# Specific monitor times (only monitor data will be output/appended for these times)"<<endl;
    fout<<"# monitor times\tunit(secs,days,hours,years)\ttime1\ttime2\ttime3...."<<endl;
    fout<<"#monitor times\tdays\t0.9\t1.5"<<endl;
    fout<<"#monitor times\tdays\t1.9\t2.9"<<endl;
    fout<<endl;
    fout<<"# VTU output variables (separated by tabs). Uncomment to use."<<endl;
    fout<<"# By default, only the full model is output."<<endl;
    fout<<"# Depending on the <vtu> keyword in the regions file,"<<endl;
    fout<<"# these variables will also be output as separate vtu's for those regions"<<endl;
    fout<<"# You may change the order or reduce the list at will."<<endl;
    fout<<"#vtu output\t";
    std::list<SimulatorSetupParameter>::iterator dictBegin ( this->GetSS()->GetParameterList().begin());
    std::list<SimulatorSetupParameter>::iterator dictEnd   ( this->GetSS()->GetParameterList().end());

    size_t counter(0);
    for( std::list<SimulatorSetupParameter>::iterator pit = dictBegin; pit!=dictEnd; pit++){
        fout<<pit->name<<"\t";
        if (counter > 5){
            counter=0;
            fout<<endl;
            fout<<"#vtu output\t";
        }
        counter++;
    }
    fout<<endl;
    fout<<endl;
    fout<<"# Monitored variables (separated by tabs). Uncomment to use."<<endl;
    fout<<"# You may also change the order or reduce the list at will. "<<endl;
    fout<<"# You can request per-model-subdomain (or full model) monitoring via subdomain keywords below."<<endl;
    fout<<endl;
    fout<<"# VERY IMPORTANT: A region must exist (if created on the fly)before the monitor attempts to calculate, or else monitor will crash"<<endl;
    fout<<"# In 3D, monitor dimension measures volume, while monitor perimeter measures surface area."<<endl;
    fout<<"# In 2D, monitor dimension measures area, while monitor perimeter measures contour lengthurface area."<<endl;
    vector<string> monitor_tokens;
    monitor_tokens.push_back("#monitor range\t");
    monitor_tokens.push_back("#monitor integral\t");
    monitor_tokens.push_back("#monitor value\t");
    monitor_tokens.push_back("#subdomain ranges\t");
    monitor_tokens.push_back("#subdomain integrals\t");
    monitor_tokens.push_back("#subdomain dimension\t");
    monitor_tokens.push_back("#subdomain perimeter\t");

    map<string,set<PLACEMENT> > token_to_placements;
    set<PLACEMENT> all_possible_placements;
    variablePlacementSet(all_possible_placements);

    token_to_placements.insert(make_pair("#monitor range\t",all_possible_placements));
    token_to_placements.insert(make_pair("#monitor integral\t",all_possible_placements));
    set<PLACEMENT> only_model;only_model.insert(MODEL);
    token_to_placements.insert(make_pair("#monitor value\t",only_model));

    for (auto itm = monitor_tokens.begin(); itm!=monitor_tokens.end();itm++)
    {
        string token=*itm;
        fout<<token;
        counter=0;
        map<string,set<PLACEMENT> >::iterator has_placement=token_to_placements.find(token);
        if (has_placement!=token_to_placements.end()){
            auto i = 0;
            for( std::list<SimulatorSetupParameter>::iterator pit = dictBegin; pit!=dictEnd; pit++){
                if( has_placement->second.find(this->GetSS()->GetModel()->Database().StorageKey(pit->name.c_str()).place)!= has_placement->second.end()
                        && this->GetSS()->GetInternalVars().find(pit->name)== this->GetSS()->GetInternalVars().end()){
                    fout<<pit->name<<"\t";
                    if (counter > 6 || i == this->GetSS()->GetParameterList().size()-1){
                        counter = 0;
                        fout<<endl;
                        if (i < this->GetSS()->GetParameterList().size()-1)
                            fout<<token;
                    }
                    counter++;
                }
                i++;
            }
            fout<<endl;
            fout<<endl;
        } else {
            vector<string> possible_subdomains;
            for ( typename map<string,Region<dim> >::const_iterator
                  git=this->GetSS()->GetModel()->UniqueRegionsBegin(); git!=this->GetSS()->GetModel()->UniqueRegionsEnd(); git++ )
                possible_subdomains.push_back((*git).first);
            for ( typename map<string,Region<dim> >::const_iterator
                  git=this->GetSS()->GetModel()->RegionsBegin(); git!=this->GetSS()->GetModel()->RegionsEnd(); git++ )
                possible_subdomains.push_back((*git).first);
            for ( auto bit=this->GetSS()->GetModel()->BoundariesBegin(); bit!=this->GetSS()->GetModel()->BoundariesEnd(); bit++ )
                possible_subdomains.push_back((*bit).first);
            for ( auto bit=this->GetSS()->GetModel()->BoundariesBegin(); bit!=this->GetSS()->GetModel()->BoundariesEnd(); bit++ )
                possible_subdomains.push_back((*bit).first);

            for( auto i = 0; i<possible_subdomains.size(); i++){
                fout<<possible_subdomains[i]<<"\t";
                if (counter > 6 || i == possible_subdomains.size()-1){
                    counter = 0;
                    fout<<endl;
                    if (i < possible_subdomains.size()-1)
                        fout<<token;
                }

                counter++;
            }
            fout<<endl;
        }
    }

    fout<<endl;
    fout<<"# Uncomment the option 'monitor all timesteps' to override the monitoring frames and output every timestep"<<endl;
    fout<<"# of ths simulation"<<endl;
    fout<<"#monitor all timesteps"<<endl;
    fout<<endl;

    fout.close();

    // Now, this will output the section specific to each simulator.
    OutputSimulatorSpecificControlFileSection();
    error_handler_.Note( EXCEPTION, "SimulatorControl<dim>::OutputSampleControlFile", "A sample control file has been output",
                           "Please check it and restart the simulator.");
}

template <uint32_t dim>
double SimulatorControl<dim>::MaxDifferenceScalarNodalProperty(Index &snp1Key, Index& snp2Key )
{
    Region<dim>& mref(this->GetSS()->GetModel()->Region("Model"));

    double diff(0.), max(0.);
    const auto nodesEnd( mref.NodesEnd() );
    for( auto it = mref.NodesBegin(); it != nodesEnd; ++it )
    {
        diff = fabs( (*it)->Read(snp1Key) - (*it)->Read(snp2Key) );
        if (diff > max)
            max = diff;
    }
    return diff;
}


template <uint32_t dim>
void SimulatorControl<dim>::Run()
{
    error_handler_.Note( FATAL_ERROR, "SimulatorControl<dim>::Run",
                           " Huh? You are using the Run() method from the base class!" ," Create a run method in your derived SimulatorControl class");
}



template <uint32_t dim>
SimulatorControl<dim>::~SimulatorControl()
{
    delete(vtu_);
}


template class SimulatorControl<1U>;
template class SimulatorControl<2U>;
template class SimulatorControl<3U>;

}
