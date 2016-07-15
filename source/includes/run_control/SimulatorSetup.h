#ifndef SIMULATORSETUP_H
#define SIMULATORSETUP_H

#include "Index.h"
#include "SimulatorSetupParameter.h"
#include "ErrorHandler.h"

#include <list>
#include <string>
#include <algorithm>
#include <iostream>

#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"
#include "PDE_Integrator.h"

#include "SimulatorMonitor.h"
#include "PropertyAtPointVisitor.h"
#include "ComputationalSettings.h"

namespace csmp {

template<size_t dim>
class SimulatorSetup {
public:
    SimulatorSetup(std::string geometry_prefix,
                   std::string project_name,
                   std::vector<bool > options,
                   bool verbose = false);
    ~SimulatorSetup();

    //-----------------------------------------------------------------
    // Initializer methods to be used by subclasses
    void LoadModel();   // Load models (2D and 3D for now)
    void Initialize(bool verbose = false);
    virtual void CreateParameterList()=0;
    virtual void CheckOptions()=0;
    void AddPreExistingModelVariables();
    std::vector<bool>& GetOptions(){return options_;}
    virtual std::vector<bool> LoadDefaultOptions();
    bool CheckVariables();
    void CreateAllProperties();
    void CreatePreviousProperties();
    bool CheckProperties();
    void CheckInputRanges();
    void CheckModel();
    void CheckModelMinMaxCoordinates();
    void OutputSampleVariablesFile(std::string filename = "SimulatorSetupSampleVariables.txt");
    void OutputSampleConfigurationFile();
    void OutputSampleRegionsFile();
    void OutputParameterListCode();
    bool ReadOptionsFromRegionsFile();
    //void ReadMonitoringFilesForRestart();
    void LoadBoundaryAndWellNames();
    void LoadWells();
    void LoadMonitor();
    void LoadControlSimulatorSetupParameters();
    void ModelSetupFromConfigFile();
    ComputationalSettings& RunSettings(){return this->run_settings_;}

    //-----------------------------------------------------------------
    // Auxiliary internal methods (left public on purpose for now!)
    //void AddKey( const std::string &name);
    //void AddKeys( );
    void OutputParameterListToScreen();


    void SetupWellRatesBasedOnRateVariables();
    void TriggerRestart(){this->restart_=true;}
    bool YesOrNo();

    //-----------------------------------------------------------------
    // Access methods
    bool Interactive(){return interactive_;}
    bool Restart(){return restart_;}
    void Restart(bool restart){this->restart_=restart;}
    Model<dim>* GetModel(){return this->model_;}

    bool Verbose(){return verbose_;}
    void Verbose(bool choice){ verbose_=choice;}

    std::string GetGeometryPrefix(){return this->geometry_file_prefix_;}
    std::string GetProjectName(){return this->name_;}

    ComputationalSettings& GetRunSettings(){return this->run_settings_;}
    //std::string& PropName(const std::string& prop);
    //Index& PropKey(const std::string &prop);
    //bool PropNameExists(const std:: std::string & prop);

    std::vector<std::string>& GetWells(){return well_regions_;}
    std::map<std::string,std::pair<std::string,std::string > >& GetWellRateVars(){return well_rate_vars_;}
    std::vector<std::string>& GetTopBoundaries(){return top_boundaries_;}
    std::vector<std::string>& GetLateralBoundaries(){return lateral_boundaries_;}
    std::vector<std::string>& GetBottomBoundaries(){return bottom_boundaries_;}

    // Output lists
    std::list<std::string>& GetVTUPropList(){return vtuOutputProps_;}
    std::list<std::string>& GetRegOutList(){return vtuOutputRegions_;}

    std::list<std::string>& GetMonitoredPropList(){return monitoredProps_;}
    std::list<std::string>& GetMonitoredRegionsList(){return monitoredRegions_;}

    std::list<std::string>& GetMonitoredRangePropList(){return monitoredRangeProps_;}
    std::list<std::string>& GetMonitoredIntegralPropList(){return monitoredIntegralProps_;}
    std::list<std::string>& GetMonitoredValuePropList(){return monitoredValueProps_;}
    std::list<std::string>& GetMonitoredRegionIntegralList() {return monitoredRegionIntegralList_;}
    std::list<std::string>& GetMonitoredRegionRangeList()    {return monitoredRegionRangeList_;}
    std::list<std::string>& GetMonitoredRegionDimensionList(){return monitoredRegionDimensionList_;}
    std::list<std::string>& GetMonitoredRegionPerimeterList(){return monitoredRegionPerimeterList_;}

    //std::list<std::string>& GetIntegralProps(){return integral_properties_;}
    //std::list<std::string>& GetRangeProps(){return range_properties_;}
    std::list<SimulatorSetupParameter>& GetParameterList(){return parameter_list_;}
    SimulatorMonitor<dim>* GetSimulatorMonitor(){return simulator_monitor_;}
    void InsertNewMonitorDataValues(double64 time);

    //-----------------------------------------------------------------
    // Integrator + solver setup management
    /** @todo in principle, integrators will be executed in the sequence they were created.
        This could be done by checking the order of the SAMG_Solvers, to avoid creating wrappers for PDE_Integrator and NCFVT classes.
        This seems easier said than done, as many things need to happen between solver calls.  Perhaps it is better to allow
        the usage of names are keys and call the correct one when desired.  This is left to the user.
    */
    void Add_DESCompatible_Integrator(std::string name);
    void AddLegacy_FE_Integrator(std::string name);
    void DeleteLegacy_FE_Integrator(std::string name); /// J.E.M. Added because we need to remove Legacy FE integrators if we switch variable flags to PLAIN.
    void AddLegacy_FEFV_Integrator(std::string name,std::vector<std::string> variables);

    PropertyDatabase<dim>& Database(){return this->GetModel()->Database();}
    std::set<std::string>& GetInternalVars(){return internal_vars_;}

#ifdef CSMP_WITH_SAMG_SOLVER
    //FEFV_Algorithm<dim>* GetDESCompatibleIntegrator(SAMG_Solver* samg) ;
    PDE_Integrator<dim,Region>* GetLegacy_FE_Integrator(SAMG_Solver* samg);
    NodeCenteredFiniteVolumeTransport<dim>* GetLegacy_FEFV_Integrator(SAMG_Solver* samg);

    //FEFV_Algorithm<dim>* GetDESCompatibleIntegrator(std::string name);
    PDE_Integrator<dim,Region>* GetLegacy_FE_Integrator(std::string name);
    NodeCenteredFiniteVolumeTransport<dim>* GetLegacy_FEFV_Integrator(std::string name);
    void QuietIntegrators();

    //std::map< SAMG_Solver*,FEFV_Algorithm<dim>*  >&  GetDESCompatibleIntegrators(){return  des_compatible_integrators_;}
    std::map< SAMG_Solver*,PDE_Integrator<dim,Region>* >&  GetLegacy_FE_Integrators(){return legacy_FE_integrators_;}
    std::map< SAMG_Solver*,NodeCenteredFiniteVolumeTransport<dim>* >& GetLegacy_FEFV_Integrators() {return legacy_FEFV_integrators_;}

    SAMG_Settings* GetIntegratorSolverSettings(std::string name);
    std::map<std::string,SAMG_Settings*>& GetAllSAMGSettings() {return solver_settings_;}
#else
    /// add extra functionality for laternative solve
#endif

protected:
    std::list<SimulatorSetupParameter> parameter_list_;
    bool created_lists_;
    //std::map<std::string,std::string> dict_;          // dictionary of necessary variables and its names
    //std::map<std::string,csmp::Index> idict_;         // dictionary of necessary variables and its keys
    bool interactive_;           // interactive

private:

    // construction
    std::string name_;                   // Prefix for -configuration.txt, -regions.txt, and -variables.txt
    std::string geometry_file_prefix_;   // Prefix for .asc and .dat files.

    std::vector<bool> options_;
    bool validated_;             // Define if settings have been validated or not.

    bool debug_;                 // Debug Mode
    bool verbose_;               // Verbose Mode
    bool restart_;               // Restart

#ifdef CSMP_WITH_SAMG_SOLVER
    // -------------------------------------
    // integrators and associated solvers.
    //std::map< SAMG_Solver*,FEFV_Algorithm<dim>* >                    des_compatible_integrators_;
    std::map< SAMG_Solver*,PDE_Integrator<dim,Region>* >             legacy_FE_integrators_;
    std::map< SAMG_Solver*,NodeCenteredFiniteVolumeTransport<dim>* > legacy_FEFV_integrators_;

    std::map<std::string,SAMG_Solver*> name_samgsolver_;
    std::map<std::string,SAMG_Settings*> solver_settings_;
#else
    /// add extra functionality for laternative solve
#endif

    // -------------------------------------
    // property and geometry groups.
    Model<dim>* model_;
    ComputationalSettings  run_settings_;
    std::vector<std::string> lateral_boundaries_,top_boundaries_,bottom_boundaries_;
    std::vector<std::string> well_regions_;
    std::map<std::string,std::pair<std::string,std::string > > well_rate_vars_; //map of well names to well_rate variables and distributed fluid volume source variables.
    std::set<std::string> internal_vars_;
    //key is the well name, pair is constructed out of region variable (e.g. volume flow rate)
    //and distributed (element placed) variable (e.g. fluid volume source).

    // -------------------------------------
    // Monitoring and Output
    SimulatorMonitor<dim>* simulator_monitor_;
    std::list<std::string> monitoredProps_;      // this list is created in the CreateParametersList method.
    std::list<std::string> monitoredRegions_;    // this list is generated when reading the -regions file.

    std::list<std::string> vtuOutputProps_;      // this list is created in the CreateParametersList method.
    std::list<std::string> vtuOutputRegions_;    // this list is generated when reading the -regions file. ("Model" is in it by default)

    std::list<std::string> monitoredRangeProps_,monitoredIntegralProps_,monitoredValueProps_; // these lists are generated by the -control.txt file.
    std::list<std::string> monitoredRegionIntegralList_, monitoredRegionRangeList_,monitoredRegionDimensionList_, monitoredRegionPerimeterList_; // these lists are generated by the -control.txt file.

    // Not used yet.
    PropertyAtPointVisitor<dim>* property_at_point_;
    std::map<std::string,std::map<size_t,std::vector<double64> > > points_to_monitor_;
    std::map<size_t,Element<dim>*> point_in_element_in_region_;

};

} // end csmp

#endif // SIMULATORSETUP_H

