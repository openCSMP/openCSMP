#include <ctime>
#include "SimulatorSetup.h"
#include "Exception.h"
#include "PropertyDatabase.h"
#include "SimulationSignalHandler.h"
#include "PL_Utilities.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "InputDataManager.h"
#include "OS_Utilities.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
SimulatorSetup<dim>::SimulatorSetup(std::string geometry_prefix,
                                    std::string project_name,
                                    vector<bool > options,
                                    bool verbose):

    created_lists_(false),
    name_(project_name),
    geometry_file_prefix_(geometry_prefix),
    options_(options),
    debug_(false),
    verbose_(verbose),
    restart_(false),
    model_(NULL),
    simulator_monitor_(NULL),
    property_at_point_(NULL)
{}

template<uint32_t dim>
void SimulatorSetup<dim>::Initialize(bool verbose)
{
    this->Verbose(verbose);

    if (this->geometry_file_prefix_=="restart")
        this->Restart(true);

    CheckOptions();

    this->LoadModel();

    CreateParameterList();

    this->LoadControlSimulatorSetupParameters();

    this->CreateAllProperties();

    this->CheckModel();

    if (!this->Restart()) this->ModelSetupFromConfigFile();

    this->ReadOptionsFromRegionsFile();

}

template<uint32_t dim>
void SimulatorSetup<dim>::LoadControlSimulatorSetupParameters()
{
    //required parameters for control mechanisms to work.
    SimulatorSetupParameter sp;
    // scalars
    sp.name="run time";            sp.notation="RT";   sp.unit="s";    sp.type=SCALAR; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="computed"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);
    sp.name="model time";          sp.notation="MT";   sp.unit="s";    sp.type=SCALAR; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="computed"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);
    sp.name="timestep";            sp.notation="TS";   sp.unit="s";    sp.type=SCALAR; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="computed"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);
    sp.name="vtu frames";          sp.notation="VTU";  sp.unit="1";    sp.type=SCALAR; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="input"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);
    sp.name="monitor frames";      sp.notation="MONF"; sp.unit="1";    sp.type=SCALAR; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="input"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);
    sp.name="interval";            sp.notation="I";    sp.unit="--";    sp.type=SCALAR; sp.min= 0.00E00; sp.max=1.00E+50; sp.placement=MODEL;   sp.usage="computed"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);
    sp.name="well radius";         sp.notation="WR";   sp.unit="m";     sp.type=SCALAR; sp.min=0.00E-50; sp.max=1.00E+00; sp.placement=REGION; sp.usage="input"; this->GetParameterList().push_back(sp);
    internal_vars_.insert(sp.name);

}

template<uint32_t dim>
void SimulatorSetup<dim>::CreateAllProperties()
  {
    // This function creates properties in the model assuming that, if restart is not an active option, the model is completely empty.
    // NEVERTHELESS, for safety, IsDefined and other precautions are left in place with our without restart.

    // first we go through all the parameter list and make sure that vsize is defined!
    for (auto lit = this->GetParameterList().begin(); lit!= this->GetParameterList().end(); lit++){
        SimulatorSetupParameter& p(*lit);
        if (p.type== SCALAR)
            p.vsize=1;
        else if (p.type== VECTOR)
            p.vsize=dim;
        else if (p.type== TENSOR)
            p.vsize=dim*dim;
        else if ((p.type== ARRAY || p.type==FLAGGEDARRAY)){
            if (p.vsize <= 0 && p.vsize > 10e6){
                p.Out();
                throw csmp::Exception(EXCEPTION,"SimulatorSetup<dim>::CreateAllProperties()", "(Flagged)Array variable has unusual size.");
            }
        }
        else {
            p.Out();
            throw csmp::Exception(EXCEPTION,"SimulatorSetup<dim>::CreateAllProperties()", "Failed to create vsize for this variable type.");
        }
    }

    //Go through all existing variables in the model, and load the sizes of arrays and flagged arrays.
    for (auto lit = this->GetParameterList().begin(); lit!= this->GetParameterList().end(); lit++)
      {
        SimulatorSetupParameter& p(*lit);

        if (!this->Database().IsDefined(p.name.c_str())) {
            if (this->Verbose()) cout<<" creating property for :"<<p.name<<" vmin:"<<p.min<<" vmax:"<<p.max<<endl;
            this->GetModel()->CreateProperty(p.name.c_str(),p.notation.c_str(),p.unit.c_str(),p.type,p.placement,p.vsize,p.min,p.max,p.usage);
            p.key=this->Database().StorageKey(p.name.c_str());
        }
        else{
            if (p.usage=="computed" && !this->Restart()) {
                this->GetModel()->Database().DeleteProperty(p.name.c_str());
                this->GetModel()->CreateProperty(p.name.c_str(),p.notation.c_str(),p.unit.c_str(),p.type,p.placement,p.vsize,p.min,p.max,p.usage);
                p.key=this->Database().StorageKey(p.name.c_str());
            }
            else {
                if (this->Verbose()) cout<<" Variable '"<<p.name<<"' already exists. if its size and type is correct, it will be loaded untouched. "<<endl;
                if (this->Database().Parameter(p.name.c_str()).key.type!=p.type)
                    throw csmp::Exception(FATAL_ERROR, "SimulatorSetup<dim>::CreateAllProperties()", " Variable type in model does not correspond to the type in the expected parameter list.");
                if (this->Database().Parameter(p.name.c_str()).usage!=p.usage)
                    throw csmp::Exception(FATAL_ERROR, "SimulatorSetup<dim>::CreateAllProperties()", " Variable usage in model does not correspond to the usage in the expected parameter list.");

            }
        }
    }
}

template<uint32_t dim>
void SimulatorSetup<dim>::CheckModel(){
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if (this->Verbose()) cout<<" SimulatorSetup<dim>::CheckModel()..."<<endl;
    if (this->Verbose()) this->CheckModelMinMaxCoordinates();

    if( !created_lists_ )
        csmp_error.Note(FATAL_ERROR,"\n SimulatorSetup<dim>::",
                          "Something is wrong in the creation of default parameter list. This is a developer issue.",
                          "\n Revise the method CreateParameterList.");

    bool associated_notations_all_variables=this->CheckVariables();
    if( !associated_notations_all_variables)
    {
        this->OutputSampleVariablesFile();
        csmp_error.Note(FATAL_ERROR,"\n SimulatorSetup<dim>::",
                          "Property name or notation is incorrect.",
                          "\n Please specify the notation of your variables according to the template file SimulatorSetupSampleVariables.txt supplied after this message is printed out.");

    }

    bool checked_range_and_placement_of_properties =this->CheckProperties();
    if ( !checked_range_and_placement_of_properties)
    {
        this->OutputSampleVariablesFile();
        csmp_error.Note(FATAL_ERROR,"\n SimulatorSetup<dim>::",
                          "Property placement or range is incorrect.",
                          "\n Please specify the notation of your variables according to the template file SimulatorSetupSampleVariables.txt supplied after this message is printed out.");

    }

    if (this->Verbose()) cout<<" Exiting SimulatorSetup<dim>::CheckModel()"<<endl;
}

/**
 * @brief SimulatorSetup<dim>::AddPreExistingModelVariables this method adds pre-existing model variables to the parameter list of the SimulatorSetup.
 * Typically, these "pre-existing" variables are loaded from the variables file when creating the model.
 * The parameter list is sorted following this addition.
 * Loading pre-existing variables will replace any existing parameters with the same name in the parameter list.
 */
template<uint32_t dim>
void SimulatorSetup<dim>::AddPreExistingModelVariables()
{
    map<string,csmp::Index> pre_existing_props_in_model;
    this->GetModel()->Database().ListVariables(pre_existing_props_in_model);
    for (map<string,csmp::Index>::iterator peim_it = pre_existing_props_in_model.begin(); peim_it != pre_existing_props_in_model.end();peim_it++){
        SimulatorSetupParameter p;
        double min,max;
        p.name=peim_it->first;
        p.notation=this->GetModel()->Database().Parameter(p.name.c_str()).notation;
        p.unit=this->GetModel()->Database().Unit(p.name.c_str());
        p.type=peim_it->second.type;
        this->GetModel()->Database().RangeOf(p.name.c_str(),min,max);
        p.min=min;
        p.max=max;
        p.placement=peim_it->second.place;
        p.usage=this->GetModel()->Database().Usage(p.name.c_str());
        std::list<SimulatorSetupParameter>::iterator findIter = std::find(this->parameter_list_.begin(), this->parameter_list_.end(), p);
        if (findIter==this->parameter_list_.end())
            this->parameter_list_.push_back(p);
        else //replace the parameter list with the version from the original variables file
            *findIter=p;

    }
    this->parameter_list_.sort(compare_setup_parameter_nocase);
}



template<uint32_t dim>
bool SimulatorSetup<dim>::YesOrNo()
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
void SimulatorSetup<dim>::ModelSetupFromConfigFile()
{
    // -----------------------------------------------------------------------
    // 1. Assigning material properties, initial conditions, and boundary
    //    conditions from *-configuration.txt files.
    // -----------------------------------------------------------------------
    InputDataManager<dim>  model_configuration;

    // first test if configuration file exists.
    ifstream  ifs;
    ifs.open( this->GetProjectName()+"-configuration.txt");

    if (!ifs.is_open()) {
        cout<<" ******************** ATTENTION **********************************"<<endl;
        cout<<" It seems that you are missing a file named: "<<this->GetProjectName()<<"-configuration.txt"<<endl;
        cout<<" This file can be generated taking into account your current region's file and variables. "<<endl;
        cout<<" NOTE: Tensors and Array variables will be skipped by this process. "<<endl;
        cout<<"       Regional properties will be assigned to the complete region. "<<endl;
        cout<<"       Configuration block for applying flags on Regions will be empty. "<<endl;
        cout<<"       Boundary conditions will be written for all boundaries, but commented out."<<endl;
        cout<<"       Boundary conditions on vector variables will be omitted (for the time being)."<<endl;
        cout<<"       All unnecessary lines can be deleted manually if desired."<<endl;
        cout<<" \nWould you like to generate a configuration file? [y/n]: "<<endl;
        if (this->YesOrNo())
            this->OutputSampleConfigurationFile();
    }
    else {
        ifs.close();
    }

    if (dim==2)
        model_configuration.ConfigureFromFile( *this->GetModel(),
                                               this->GetProjectName().c_str(),
                                               false,   // groupname from parameter range
                                               true,    // default property values
                                               true,    // regional property values
                                               false,   // boundary conditions for box-shaped model
                                               true,    // essential conditions for regions
                                               true,    // properties & essential conditions for csmp::Boundaries
                                               this->RunSettings() );

    else if (dim==3)
        model_configuration.ConfigureFromFile( *this->GetModel(), this->GetProjectName().c_str(),
                                               false,   // groupname from parameter range
                                               true,    // default property values
                                               true,    // regional property values
                                               false,   // boundary conditions for box-shaped model
                                               true,    // essential conditions for regions
                                               true,    // properties & essential conditions for csmp::Boundaries
                                               this->RunSettings() );
    else
        throw csmp::Exception( ERROR, "SimulatorSetup<dim>::ModelSetupFromConfigFile()",
                               "can only configure 2D and 3D models from file." );


}

template<uint32_t dim>
bool SimulatorSetup<dim>::CheckVariables()
{

    if (this->Verbose()) this->OutputParameterListToScreen();
    //cin.get();
    //================================================================================
    // CHECK MODEL PARAMETER LIST AND LOCAL SETUP CLASS PARAMETER LIST FOR DUPLICATE VARIABLE NAMES
    const PropertyDatabase<dim>& pdb=model_->Database();
    auto modelplistBegin ( pdb.Begin());
    auto modelplistEnd ( pdb.End());

    // ---------------------------------------------------------------------
    // extract the names and notations from the Property Database for comparison.
    // detect duplicate names and/or notations
    // Checking model's variables.
    set<std::string> list_of_names;
    multimap<std::string,std::string> duplicates;

    for(auto pit = modelplistBegin; pit!=modelplistEnd;pit++){
        if (list_of_names.find(pit->second.name)==list_of_names.end())
            list_of_names.insert(pit->second.name);
        else{
            duplicates.insert(std::make_pair(pit->second.notation,pit->first));
        }
    }

    if (!duplicates.empty()){
        cout<<"---------------------------------------------------"<<endl;
        cout<<" Duplicate parameter found in model database."<<endl;
        for (auto mit = duplicates.begin(); mit!=duplicates.end();mit++)
            cout<<" notation: "<<mit->first<<" name: "<<mit->second<<endl;
        string errorMessage(" Detected same notation for at least two variables inside the Property Database. This comes from the used variables file!");
        throw csmp::Exception( FATAL_ERROR, "SimulatorSetup<dim>::CheckVariables():",errorMessage.c_str() );
    }
    // ---------------------------------------------------------------------
    // Checking setup class' variables.
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    list_of_names.clear();
    //We check that there are no duplicates in the desired parameter list as well.
    //extract the names and notations from the Property Database for comparison.
    for(std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin;pit!=dictEnd;pit++){
        if (list_of_names.find(pit->name)==list_of_names.end())
            list_of_names.insert(pit->name);
        else
            duplicates.insert(std::make_pair(pit->notation,pit->name));
    }

    if (!duplicates.empty()){
        cout<<"---------------------------------------------------------------------------------------"<<endl;
        cout<<" Duplicate parameter list from your setup class. (not necessarily your variables file!)"<<endl;
        for (auto mit = duplicates.begin(); mit!=duplicates.end();mit++)
            cout<<" notation: "<<mit->first<<" name: "<<mit->second<<endl;
        string errorMessage(" Detected same notation for at least two variables. This comes from the desired parameter list in your own setup class.");
        throw csmp::Exception( FATAL_ERROR, "SimulatorSetup<dim>::CheckVariables():",errorMessage.c_str() );
    }

    // -----------------------------------------------------------------------
    // detect missing variables.
    // now we construct the dictionary and check for missing variables.
    std::map<std::string,std::string > missing_vars;
    for (auto classpit= this->GetParameterList().begin(); classpit != this->GetParameterList().end();classpit++){
        bool found=false;
        for(auto pit = modelplistBegin;        pit!=modelplistEnd;pit++)
            if (pit->second.name==classpit->name)
                found=true;
        if (!found)
            missing_vars.insert(make_pair(classpit->name,classpit->notation));
    }

    if( missing_vars.size() > 0 ){

        cerr<<"\n SimulatorSetup<dim>::CheckVariables"<<endl;
        cerr<<"\n - "<<missing_vars.size()<<" variables are missing (out of "<<parameter_list_.size()<<" needed) from your configuration file:"<<endl;
        cerr<<"\n These variables are needed for this simulator. Please find the generated text file with the missing variables and add them to your variables file."<<endl;
        cerr<<"\n You are free to change the names and ranges of the variables however notations and placements are strict."<<endl;
        cerr<<"\n Summary:"<<endl;

        for( std::map<std::string,std::string >::const_iterator pit = missing_vars.begin();pit!=missing_vars.end();pit++)
            cerr<<"\n Property = "<<(*pit).first<<" , Notation = "<<(*pit).second;
        cerr<<endl;
        cerr<<"\n End of missing variables summary."<<endl;
        cerr<<endl;
        cerr<<endl;

        //this->AddPreExistingModelVariables();
        //this->OutputParameterListCode();

        return false;
    }
    return true;
}

template<uint32_t dim>
bool SimulatorSetup<dim>::CheckProperties()
{
    // Check placement, type and usage status of the variables
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    const PropertyDatabase<dim>& p=model_->Database();
    double vmin, vmax;
    std::string property_name, notation_name;
    bool no_mistakes(true);

    for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){

        notation_name = pit->notation;
        property_name = pit->name;
        p.RangeOf( property_name.c_str() , vmin, vmax );

        if ( (   p.Type( property_name.c_str() ) != pit->type
                 || p.Placement( property_name.c_str() ) != pit->placement
                 || p.Usage( property_name.c_str() ) != pit->usage ))
        {
            cerr<<" Variable '"<<property_name<<"'( min = "<<vmin<<", max = "<<vmax<< "):";
            cerr<<" must be a '"<<parseType(pit->type)<<" but is a '"<<parseType(p.Type( property_name.c_str() ))<<"'"<<endl;
            cerr<<" placed on the '"<<parsePlacement(pit->placement)<<"' but is placed on '"<<parsePlacement(p.Placement( property_name.c_str() ))<<"'"<<endl;
            cerr<<" should be within range ( min = "<<pit->min<<", max = "<<pit->max<< "):"<<endl;
            no_mistakes = false;
            cin.get();
        }


    }
    return no_mistakes;
}


template <uint32_t dim>
void SimulatorSetup<dim>::CheckModelMinMaxCoordinates(){
    Point<dim> minp,maxp;
    this->GetModel()->MinMaxCoordinates(minp,maxp);
    cout<<" min max coords:"<<endl;
    cout<<" max :"; maxp.Out();
    cout<<" min :"; minp.Out();
}

template <uint32_t dim>
void SimulatorSetup<dim>::CheckInputRanges(){


    // Check placement, type and usage status of the variables
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    const PropertyDatabase<dim>& p=model_->Database();
    double vmin, vmax;
    std::string property_name, usage;
    bool no_mistakes(true);
    for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){
        property_name = pit->name;

        usage=p.Usage( property_name.c_str());
        std::transform(usage.begin(), usage.end(), usage.begin(), ::tolower);
        this->GetModel()->MinMaxOf( property_name.c_str(), vmin, vmax );

        if ((isnan(vmin)|| isnan(vmax)) && usage=="input" ) {
            cout<<" property name (SS): "<<property_name<<" vmin:"<<vmin<<" vmax:"<<vmax<<endl;
            no_mistakes=false;
            cerr<<" Variable '"<<property_name<<"'( min = "<<vmin<<", max = "<<vmax<< "):";
            cerr<<" has usage: '"<<usage<<"'. Yet it has not been initialized properly."<<endl;
        }
    }

    if (!no_mistakes)
    {
        cerr<<" Please input correct static and/or initial conditions for properties listed above, "<<endl;
        cerr<<" or correct their usage to 'computed', if you will be calculating them inside your algorithm. "<<endl;
        string errorMessage(" One or more input variables are incorrect. Check your -configuration.txt file.");
        throw csmp::Exception( FATAL_ERROR, "SimulatorSetup<dim>::CheckInputRanges():", errorMessage );
    }

}



// Read output options as well as region types (whether they are wells, boundaries)

template <uint32_t dim>
bool SimulatorSetup<dim>::ReadOptionsFromRegionsFile()
{
    std::string  text_line;
    std::string token;

    std::ifstream ifs( string(this->GetProjectName()+"-regions.txt").c_str() );

    if ( !ifs.is_open() )
        return false;

    vtuOutputRegions_.push_back("Model"); // by default

    // Skip the two first lines
    std::getline( ifs, text_line );
    std::getline( ifs, text_line );

    std::getline( ifs, text_line );

    /// @todo flag if eof is reached.  Perhaps throw an exception.  Some region should be in there!
    if( ifs.eof())
        return true;

    cout<<"\n Reading options from regions file..."<<endl;
    // Reading specific names
    while ((text_line.length() >= 1) && (!ifs.eof()))
    {
        if( text_line[0] != '#')
        {
            vector<string> listoftokens;
            std::stringstream iss(text_line);

            while (getline(iss,token,'\t'))
                listoftokens.push_back(token);

            if (this->Verbose()) cout<<" SimulatorSetup<dim> Reading options for region/boundary/well name: '"<<listoftokens[0]<<"'"<<endl;
            // name of region cannot exist twice in a list.
            for (vector<string>::iterator lit = ++listoftokens.begin(); lit!=listoftokens.end();lit++){
                if (this->GetModel()->ContainsRegion( listoftokens[0].c_str() )||this->GetModel()->ContainsBoundary(listoftokens[0]) )
                {
                    string tok=*lit;
                    std::transform(tok.begin(), tok.end(), tok.begin(), ::toupper);
                    if (std::find(lateral_boundaries_.begin(), lateral_boundaries_.end(),listoftokens[0] ) == lateral_boundaries_.end()
                            && std::find(bottom_boundaries_.begin(), bottom_boundaries_.end(),listoftokens[0] ) == bottom_boundaries_.end()
                            && std::find(top_boundaries_.begin(), top_boundaries_.end(),listoftokens[0] ) == top_boundaries_.end()
                            && std::find(well_regions_.begin(), well_regions_.end(),listoftokens[0] ) == well_regions_.end())
                    {
                        if (tok.find("LATERAL")!=std::string::npos)
                            lateral_boundaries_.push_back(listoftokens[0]);

                        if (tok.find("BOTTOM")!=std::string::npos)
                            bottom_boundaries_.push_back(listoftokens[0]);

                        if (tok.find("TOP")!=std::string::npos)
                            top_boundaries_.push_back(listoftokens[0]);

                        if (tok.find("WELL")!=std::string::npos)
                            well_regions_.push_back(listoftokens[0]);
                    }

                    if( std::find(vtuOutputRegions_.begin(), vtuOutputRegions_.end(),listoftokens[0] ) == vtuOutputRegions_.end() ){
                        if (tok.find("VTU")!=std::string::npos)
                            vtuOutputRegions_.push_back(listoftokens[0]);
                    }
                }
            }
        }
        std::getline( ifs, text_line  );
    }


    cout<<" Using the following as lateral boundaries: ";
    if (lateral_boundaries_.empty())
        cout<<" None"<<endl;
    else
    {
        for (vector<string>::iterator sit=lateral_boundaries_.begin();sit!=lateral_boundaries_.end();sit++)
            cout<<"  "<<*sit;
        cout<<endl;
    }

    cout<<" Using the following as top boundaries: ";
    if (top_boundaries_.empty())
        cout<<" None"<<endl;
    else
    {
        for (vector<string>::iterator sit=top_boundaries_.begin();sit!=top_boundaries_.end();sit++)
            cout<<"  "<<*sit;
        cout<<endl;
    }

    cout<<" Using the following as bottom boundaries: ";
    if (bottom_boundaries_.empty())
        cout<<" None"<<endl;
    else
        cout<<endl;
    for (vector<string>::iterator sit=bottom_boundaries_.begin();sit!=bottom_boundaries_.end();sit++)
        cout<<"  "<<*sit;
    cout<<endl;

    cout<<" Using the following as well regions: ";
    if (well_regions_.empty())
        cout<<" None"<<endl;
    else
        cout<<endl;
    for (vector<string>::iterator sit=well_regions_.begin();sit!=well_regions_.end();sit++)
        cout<<"  "<<*sit;
    cout<<endl;

    cout<<" Output Regions to VTU : ";
    if (vtuOutputRegions_.empty())
        cout<<" None"<<endl;
    else
        cout<<endl;
    for (list<string>::iterator sit=vtuOutputRegions_.begin();sit!=vtuOutputRegions_.end();sit++)
        cout<<"  "<<*sit;
    cout<<endl;

    cout<<" Output monitoring data for : ";
    if (monitoredRegions_.empty())
        cout<<" None"<<endl;
    else
        cout<<endl;
    for (list<string>::iterator sit=monitoredRegions_.begin();sit!=monitoredRegions_.end();sit++)
        cout<<"  "<<*sit;
    cout<<endl;

    ifs.close();
    return true;
}

template <uint32_t dim>
void SimulatorSetup<dim>::SetupWellRatesBasedOnRateVariables()
{
    /** @todo This should probably use a map of well names vs which field is supposed to be used to read in the rates, as well as which phase, if other than single
        phase flow is intended.
        @todo first thing to do is check that all input information froms from fields defined on the region.
    */
    static bool first_call(true);

    if (well_rate_vars_.empty() && first_call ) {
        cout<<"  WARNING: You are calling SetupWellsBasedOnRateVariables(), while your well_rate_vars_ is empty."<<endl;
        cout<<"  It is likely that your wells will NOT operate properly, or not at all."<<endl;
        cout<<"  Please press enter to confirm and continue."<<endl;
        cin.get();
        first_call=false;
    }

    if (first_call ) {
        Region<dim>& rref(this->GetModel()->Region("Model"));
        /// @todo this needs to be constructed better.  Better logic for other users to create their own simulators.
        for (auto wrvit = this->GetWellRateVars().begin(); wrvit != this->GetWellRateVars().end();wrvit++)
            rref.InputPropertyValue(wrvit->second.second.c_str(),makeScalar(PLAIN,0.0));
        first_call=false;
    }

    if (!well_regions_.empty()) {
        for( auto mit = well_rate_vars_.begin(); mit!=well_rate_vars_.end(); mit++){
            if (this->Verbose()) cout<<" Setting up fluid volume sources from well rates...";
            if ( model_->Database().Placement( mit->second.first.c_str() )==REGION)
            {
                ScalarVariable vflowrate;
                for (vector<std::string>::iterator sit = well_regions_.begin();sit != well_regions_.end();sit++)
                {
                    // get the min and max rates.
                    model_->Region((*sit).c_str()).Read( this->Database().StorageKey(mit->second.first.c_str()), vflowrate );

                    if ( isnan(vflowrate()) )
                        throw csmp::Exception( ERROR, "SimulatorSetup::SetupWellsBasedOnRates()","nan injection rate detected" );

                    // The volume method calculates volume, surface, or length depending on the type of element that composes the well.
                    double wellvolume=model_->Region((*sit).c_str()).Volume();
                    auto ebegin=model_->Region((*sit).c_str()).CellsBegin();
                    auto eend=model_->Region((*sit).c_str()).CellsEnd();
                    // calculate volume flow rate per unit length of well (this assumes it is uniform throughout the well!)
                    /// Special Note: Here, we do not divide by the thickness (to get [m^3/(m^3*s)]) because the PDE_Integrator
                    /// does not multiply by any thickness when performing integrals over line elements. If it did,
                    /// then we would be cancelling our division by thickness here, with the multiplication by thickness when we
                    /// integrate the fluid volume source with FEM at the PDE_Integrator.  As a result, we do not need to divide by thickness.
                    for ( auto eit=ebegin;eit!=eend;eit++)
                        (*eit)->Store(this->Database().StorageKey(mit->second.second.c_str()) ,makeScalar(PLAIN,vflowrate()/wellvolume));

                    if (!containsLineElements(model_->Region((*sit).c_str()))) {
                        if (fabs(vflowrate() - model_->Region((*sit).c_str()).VolumeIntegral( mit->second.second.c_str(),false,this->Verbose())) > 10.e-15) {
                            cout<<" rate: "<<vflowrate<< " well: "<<(*sit)<<" integrated dist. rate: "<<model_->Region((*sit).c_str()).VolumeIntegral(mit->second.second.c_str(),false,false)<<endl;
                            cout<<" rate diff: "<<vflowrate-model_->Region((*sit).c_str()).VolumeIntegral(mit->second.second.c_str(),false,false)<<endl;
                            throw csmp::Exception( ERROR, "SimulatorSetup::SetupWellsBasedOnRates()",
                                                   "distribution of well rate for 2D or 3D  elements to 'fluid volume source' " );
                        }
                    }
                    else{
                        if (fabs(vflowrate() - model_->Region((*sit).c_str()).VolumeIntegral(mit->second.second.c_str(),false,this->Verbose())) > 10.e-15) {
                            cout<<" rate: "<<vflowrate<< " well: "<<(*sit)<<" integrated dist. rate: "<<model_->Region((*sit).c_str()).VolumeIntegral(mit->second.second.c_str(),false,false)<<endl;
                            cout<<" rate diff: "<<vflowrate-model_->Region((*sit).c_str()).VolumeIntegral( mit->second.second.c_str(),false,false)<<endl;
                            throw csmp::Exception( ERROR, "SimulatorSetup::SetupWellsBasedOnRates()",
                                                   "distribution of well rate to 'fluid volume source' does not pass check." );
                        }
                    }
                }

            }
            else
            {
                throw csmp::Exception( ERROR, "SimulatorSetup::SetupWellsBasedOnRates()",
                                       "Well rates may only be placed on the REGION for now." );
            }

        }
    }
}

template<uint32_t dim>
void SimulatorSetup<dim>::OutputSampleVariablesFile(string filename)
{
    ofstream fout(filename);
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    std::string property_name, notation_name;

    cout<<" Outputting sample file to SimulatorSetupSampleVariables.txt"<<endl;
    fout<<"name"<<"\t"<<"MParameter"<<"\t"<<"unit"<<"\t"<<"type"<<"\t"<<"min."<<"\t"<<"max."<<"\t"<<"place"<<"\t"<<"usage"<<endl;
    fout<<endl;
    /// Note: if properties in model match those in parameter list, then keep the original list version. (Julian 1-07-2014)
    for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){

        notation_name = pit->notation;
        property_name = pit->name;

        PLACEMENT place=pit->placement;
        VARIABLE_TYPE type=pit->type;
        if (notation_name=="OT"){
            //cout<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<1000<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<parsePlacement(place)<<"\t"<<pit->usage<<endl;
            fout<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<1000<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<parsePlacement(place)<<"\t"<<pit->usage<<endl;
        }
        else if (notation_name=="MONT"){
            //cout<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<10000<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<parsePlacement(place)<<"\t"<<pit->usage<<endl;
            fout<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<10000<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<parsePlacement(place)<<"\t"<<pit->usage<<endl;
        }
        else {
            //cout<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<type<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<parsePlacement(place)<<"\t"<<pit->usage<<endl;
            fout<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<type<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<parsePlacement(place)<<"\t"<<pit->usage<<endl;
        }
    }
    fout.close();

    // template
}

template<uint32_t dim>
void SimulatorSetup<dim>::OutputSampleConfigurationFile()
{
    ofstream fout(this->GetProjectName()+"-configuration.txt");
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    std::string property_name;
    cout<<" Outputting sample file to "<<this->GetProjectName()<<"-configuration.txt"<<endl;
    fout<<this->GetProjectName()<<endl;
    fout<<endl;  // skip one line.

    // output block 1  (global definitions)
    fout<<"# (1) global (default) material properties"<<endl;
    for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){
        property_name = pit->name;
        if (internal_vars_.find(property_name)==internal_vars_.end()){
            VARIABLE_TYPE type=pit->type;
            if (pit->usage=="input" && type==SCALAR){
                fout<<property_name<<"\t"<<pit->min<<endl;
            }
            else if (pit->usage=="input" && type==VECTOR){
                fout<<property_name<<"\t0.0";
                if (dim==2)
                    fout<<"\t0.0";
                if (dim==3)
                    fout<<"\t0.0"<<"\t0.0";
                fout<<endl;
            }
        }
    }
    fout<<endl; //skip line between blocks

    // output block 2  (regional definitions)
    typename std::map<std::string,csmp::Region<dim> >::iterator rbeg=this->GetModel()->UniqueRegionsBegin();
    typename std::map<std::string,csmp::Region<dim> >::iterator rend=this->GetModel()->UniqueRegionsEnd();
    typename std::map<std::string,csmp::Region<dim> >::iterator rit;
    fout<<"# (2) Regional properties."<<endl;
    for (rit = rbeg; rit!=rend;rit++){
        for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){
            property_name = pit->name;
            if (internal_vars_.find(property_name)==internal_vars_.end()){
                VARIABLE_TYPE type=pit->type;
                if (pit->usage=="input" && type==SCALAR && rit->first!="Model"){
                    fout<<"#"<<rit->first<<"\t"<<"complete"<<"\t";
                    fout<<property_name<<"\t"<<pit->min<<endl;
                }
                else if (pit->usage=="input" && type==VECTOR && rit->first!="Model"){
                    fout<<"#"<<rit->first<<"\t"<<"complete"<<"\t";
                    fout<<property_name<<"\t0.0";
                    if (dim==2)
                        fout<<"\t0.0";
                    if (dim==3)
                        fout<<"\t0.0"<<"\t0.0";
                    fout<<endl;
                }
            }
        }
    }
    fout<<endl; //skip line between blocks

    fout<<"# (3) Essential conditions on regions."<<endl;
    fout<<"# triggers imposition of either default or local regional properties in a region"<<endl;
    fout<<endl;

    typename std::map<std::string,csmp::Boundary<dim> >::iterator bbeg=this->GetModel()->BoundariesBegin();
    typename std::map<std::string,csmp::Boundary<dim> >::iterator bend=this->GetModel()->BoundariesEnd();
    typename std::map<std::string,csmp::Boundary<dim> >::iterator bit;

    fout<<"# (4) Boundary conditions. Apply flag and value to variable on boundary."<<endl;
    for (bit = bbeg; bit!=bend;bit++){
        for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){
            property_name = pit->name;
            if (internal_vars_.find(property_name)==internal_vars_.end()){
                VARIABLE_TYPE type=pit->type;
                if (type==SCALAR){
                    fout<<"#"<<bit->first<<"\t"<<"complete"<<"\t"<<"Dirichlet"<<"\t";
                    fout<<property_name<<"\t"<<pit->min<<endl;
                }
            }
        }
        fout<<"#--------------------------------------------"<<endl;
    }
    fout<<endl;

    fout<<"# (5) Computational Settings"<<endl;
    fout<<"#(this does not affect the simulation, it only exists for the model constructor)"<<endl;
    fout<<"time stepping\tCAREFUL"<<endl;
    fout<<"duration\t0.0"<<endl;
    fout.close();
    throw csmp::Exception( EXCEPTION, "SimulatorSetup<dim>::OutputSampleConfigurationFile",
                           " A sample configuration file has been output. Please restart the simulator." );
}

template<uint32_t dim>
void SimulatorSetup<dim>::OutputSampleRegionsFile()
{
    //for input
    ifstream  ifsext( (this->GetGeometryPrefix()+".asc").c_str() );
    string stext_line;

    //for output
    ofstream fout(this->GetProjectName()+"-regions.txt");
    set<string> regions;
    cout<<" Outputting sample file to "<<this->GetProjectName()<<"-regions.txt"<<endl;
    fout<<this->GetProjectName()<<endl;
    fout<<"no properties"<<endl;

    if (ifsext.is_open())
    {
        //read the first line, containing all the column headers.
        getline(ifsext,stext_line);
        getline(ifsext,stext_line);  // skip two lines

        getline(ifsext,stext_line);
        istringstream istline(stext_line);
        string token;
        getline(istline, token,' ');
        size_t num_regs=atoi(token.c_str()); // get the number of regions to read.

        getline(ifsext,stext_line); // skip one more line
        auto i = 0;
        while (i < num_regs) {
            getline(ifsext,stext_line);
            istringstream istdata(stext_line);
            getline(istdata, token,' ');
            regions.insert(token);
            i++;
        }
        list<string> ordered_list;
        for (set<string>::iterator sit = regions.begin(); sit != regions.end() ; sit++)
            ordered_list.push_back(*sit);
        ordered_list.sort();
        for (list<string>::iterator lit = ordered_list.begin() ; lit != ordered_list.end(); lit++){
            cout<<" Region: "<<*lit<<endl;
            fout<<*lit<<endl;
        }
    }
    else {
        throw csmp::Exception( ERROR, "SimulatorSetup<dim>::OutputSampleRegionsFile()",
                               " Could not open .asc file for reading!" );
    }
    throw csmp::Exception( EXCEPTION, "SimulatorSetup<dim>::OutputSampleRegionsFile()",
                           " A sample regions file has been output. Please restart the simulator." );
}

template<uint32_t dim>
void SimulatorSetup<dim>::OutputParameterListCode()
{
    ofstream fout("SimulatorSetupParameterListCode.txt");
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    std::string property_name, notation_name;

    cout<<" Outputting sample C++ parameter setup code to SimulatorSetupParameterListCode.txt"<<endl;
    fout<<"name"<<"\t"<<"MParameter"<<"\t"<<"unit"<<"\t"<<"type"<<"\t"<<"min."<<"\t"<<"max."<<"\t"<<"place"<<"\t"<<"usage"<<endl;
    fout<<endl;
    /// Note: if properties in model match those in parameter list, then keep the original list version. (Julian 1-07-2014)
    for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){

        notation_name = pit->notation;
        property_name = pit->name;

        PLACEMENT place=pit->placement;
        VARIABLE_TYPE type=pit->type;
        fout<<"p.name='"<<property_name<<
              "';p.notation='"<<notation_name<<
              "';p.unit='"<<pit->unit<<
              "';p.type="<<parseType(type)<<
              ";p.min="<<pit->min<<
              ";p.max="<<pit->max<<
              ";p.placement="<<parsePlacement(place)<<
              ";p.usage='"<<pit->usage<<"';this->parameter_list_.push_back(p);"<<endl;

    }
    fout.close();

    // template
}

template<uint32_t dim>
void SimulatorSetup<dim>::OutputParameterListToScreen()
{
    parameter_list_.sort(compare_setup_parameter_nocase);
    std::list<SimulatorSetupParameter>::const_iterator dictBegin ( parameter_list_.begin());
    std::list<SimulatorSetupParameter>::const_iterator dictEnd   ( parameter_list_.end());
    std::string property_name, notation_name;

    cout<<"name"<<"\t"<<"MParameter"<<"\t"<<"unit"<<"\t"<<"type"<<"\t"<<"min."<<"\t"<<"max."<<"\t"<<"place"<<"\t"<<"usage"<<endl;
    for( std::list<SimulatorSetupParameter>::const_iterator pit = dictBegin; pit!=dictEnd; pit++){

        notation_name = pit->notation;
        property_name = pit->name;

        PLACEMENT place=pit->placement;
        VARIABLE_TYPE type=pit->type;

        cout<<"\n "<<property_name<<"\t"<<notation_name<<"\t"<<pit->unit<<"\t"<<type<<"\t"<<pit->min<<"\t"<<pit->max<<"\t"<<place<<"\t"<<pit->usage;
    }
    cout<<endl;
    cout<<endl;
}

template<uint32_t dim>
vector<bool> SimulatorSetup<dim>::LoadDefaultOptions()
{

    throw csmp::Exception( ERROR, "SimulatorSetup::GetDefaultOptions()",
                           " Your default options should be set in the derived class." );
}


template<uint32_t dim>
void SimulatorSetup<dim>::LoadModel()
{
    //****************************************************************************************************
    //  Creating geometry
    cout <<" Building the Geometrical Model for "<<this->name_<<endl;
    if (!restart_)
    {
        ifstream  ifs;
        ifs.open( this->GetProjectName()+"-regions.txt");

        if (!ifs.is_open()) {
            cout<<" ******************** ATTENTION **********************************"<<endl;
            cout<<" It seems that you are missing a file named: "<<this->GetProjectName()<<"-regions.txt"<<endl;
            cout<<" This file can be generated using the file: "<<geometry_file_prefix_<<".asc"<<endl;
            cout<<" NOTE: All names present in the file will be used as names. "<<endl;
            cout<<"       It is up to you to delete the ones that should not exist in your simulator run. "<<endl;
            cout<<"       You should then add any specific tokens next to each region for which you want specific output"<<endl;
            cout<<"       Possible tokens are: lateral, bottom, top, well, vtu, monitor."<<endl;
            cout<<" \nWould you like to generate a regions file? [y/n]: "<<endl;
            if (this->YesOrNo())
                this->OutputSampleRegionsFile();
        }
        else {
            ifs.close();
        }
        const char* null_variable_file(NULL);
        switch( dim ){
        case 2:
            this->model_ = dynamic_cast<Model<dim>*> ( new ANSYS_Model2D(geometry_file_prefix_.data(),name_.c_str(),null_variable_file));
            break;
        case 3:
            this->model_ = dynamic_cast<Model<dim>*> ( new ANSYS_Model3D(geometry_file_prefix_.data(),name_.c_str(),null_variable_file,true,true));
            break;
        default:
            throw csmp::Exception( ERROR, "SimulatorSetup<dim>::LoadModel()",
                                   "Only 2D and 3D models are supported" );
        }
    }
    else{

        string file0_name(string(this->name_+"_restartFile0"));
        string file1_name(string(this->name_+"_restartFile1"));
        set<string> file_list;

        FILE* pFile;
        pFile = fopen (string(file0_name+".vset").c_str(),"r");
        if (pFile!=NULL){
            file_list.insert(file0_name+".vset");
            fclose(pFile);
        }

        pFile = fopen (string(file1_name+".vset").c_str(),"r");
        if (pFile!=NULL){
            file_list.insert(file1_name+".vset");
            fclose(pFile);
        }

        string filename;
        if (file_list.size() == 2 )
            filename = CompareFileModifiedTimeStamps((file0_name+".vset"),file1_name+".vset");
        else
            filename = (*file_list.begin());
        pFile = fopen (filename.c_str(),"r");


        if (pFile!=NULL) {
            fclose(pFile);
            cout<<" Restart option has been toggled. Reading from saved binary restart file."<<endl;
            cout<<" using filename: "<<filename<<" press <enter>"<<endl;
//            cin.get();
            switch( dim ){
            case 2:
                model_ = new Model<dim>(filename.substr(0, filename.size()-5));
              break;
            case 3:
                model_ = new Model<dim>(filename.substr(0, filename.size()-5));
              break;
            default:
              throw csmp::Exception( ERROR, "SimulatorSetup<dim>::LoadModel:",
                                    "Only the restart of 2D and 3D models is supported." );
            }
        }
        else {
            throw csmp::Exception( ERROR, "SimulatorSetup<dim>::LoadModel:",
                                   "Failed to restart simulation. Are the files missing or corrupted?" );
        }
    }

    this->GetModel()->Verbose(this->Verbose());
    if (this->Restart()) {
        cout <<" Finished loading the model from restart files."<<endl;
//        cin.get();
    }
    cout <<" Finished Building the Geometrical Model...Project Name:"<<this->name_<<endl;

}

template<uint32_t dim>
void SimulatorSetup<dim>::LoadMonitor()
{
    // This is for the legacy output of the region monitor.
    cout<<" Loading Simulation Monitor...";
    list<string> empty_list;
    this->simulator_monitor_=new SimulatorMonitor<dim>(this->GetProjectName(),
                                                       *this->GetModel(),
                                                       this->Restart(),
                                                       this->GetMonitoredIntegralPropList(),
                                                       this->GetMonitoredRangePropList(),
                                                       this->GetMonitoredValuePropList(),
                                                       this->GetMonitoredRegionIntegralList(),
                                                       this->GetMonitoredRegionRangeList(),
                                                       this->GetMonitoredRegionDimensionList(),
                                                       this->GetMonitoredRegionPerimeterList(),
                                                       empty_list, /// not active for now. Until I get time for testing and implementing control file option. Julian 22-08-2014
                                                       empty_list, /// not active for now. Until I get time for testing and implementing control file option. Julian 22-08-2014
                                                       this->Verbose());

    cout<<" Done."<<endl;
}



template<uint32_t dim>
SimulatorSetup<dim>::~SimulatorSetup()
{
#ifdef CSMP_WITH_SAMG_SOLVER
    /// @todo investigate if this way of deleting pointers actually works better than simply calling delete on it->first and it->second.

    for (auto it = solver_settings_.begin();it != solver_settings_.end();it++){
        it->second=NULL;  // Setting pointers to null, as these will be deleted by the solver class when destructed.
    }

    for (auto it = name_samgsolver_.begin();it != name_samgsolver_.end();it++){
        SAMG_Solver* solver=it->second;
        it->second=NULL;
        // Avoid deleting solvers existing in the legacy_FEFV_integrators. Are deleted by the NCFTAlgorithm class. -- Julian 24-09-2013
        if (legacy_FEFV_integrators_.find(solver)==legacy_FEFV_integrators_.end())
            if (solver)
                delete(solver);
    }

    /*
    Roman, 2014. Deactivate FEFV_Algorithm for reconstruction purposes
    for (auto it = des_compatible_integrators_.begin();it != des_compatible_integrators_.end();it++){
        auto integrator = it->second;
        it->second=NULL;
        if (integrator)
            delete(integrator);
    }
    */

    for (auto it = legacy_FE_integrators_.begin();it != legacy_FE_integrators_.end();it++){
        auto integrator = it->second;
        it->second=NULL;
        if (integrator)
            delete(integrator);
    }

    for (auto it = legacy_FEFV_integrators_.begin();it != legacy_FEFV_integrators_.end();it++){
        auto integrator = it->second;
        it->second=NULL;
        if (integrator)
            delete(integrator);
    }
    if (simulator_monitor_)
        delete(simulator_monitor_);
    if (model_)
        delete(model_);

#else
    /// add extra functionality for alternative solver if needed
#endif

}

#ifdef CSMP_WITH_SAMG_SOLVER
/*
Roman, 2014. Deactivate FEFV_Algorithm for reconstruction purposes
template<uint32_t dim>
void SimulatorSetup<dim>::Add_DESCompatible_Integrator(std::string name)
{

    if (solver_settings_.find(name)==solver_settings_.end())
    {
        SAMG_Solver* sol = new SAMG_Solver();
        name_samgsolver_[name]=sol;
        solver_settings_[name]=sol->SolverSettings();
        FEFV_Algorithm<dim>* pint=new FEFV_Algorithm<dim> (sol);
        des_compatible_integrators_[sol]=pint;
    }
}
*/

template<uint32_t dim>
void SimulatorSetup<dim>::AddLegacy_FE_Integrator(std::string name)
{
    if (solver_settings_.find(name)==solver_settings_.end())
    {
        SAMG_Solver* sol = new SAMG_Solver();
        name_samgsolver_[name]=sol;
        solver_settings_[name]=dynamic_cast<SAMG_Settings*>(sol->GetSolverSettings());
        PDE_Integrator<dim,Element>* pint=new PDE_Integrator<dim,Element>(*sol);
        legacy_FE_integrators_[sol]=pint;
    }
}

template<uint32_t dim>
void SimulatorSetup<dim>::DeleteLegacy_FE_Integrator(std::string name)
{
    if (solver_settings_.find(name)!=solver_settings_.end())
    {
        SAMG_Solver* sol = name_samgsolver_[name];
        PDE_Integrator<dim,Element>* pint=GetLegacy_FE_Integrator(name);
        //auto isol=name_samgsolver_.find(name);
        auto iset=solver_settings_.find(name);
        //auto iint=legacy_FE_integrators_.find(name_samgsolver_[name]);
        //delete(settings);
        delete(sol);
        delete(pint);
        //name_samgsolver_.erase(isol);
        solver_settings_.erase(iset);
        //legacy_FE_integrators_.erase(iint);
    }
}

template<uint32_t dim>
void SimulatorSetup<dim>::AddLegacy_FEFV_Integrator(std::string name,vector<string> variables)
{
    bool second_order_in_space(false),second_order_in_time(false);
    string name_upper=name;
    std::transform(name_upper.begin(), name_upper.end(), name_upper.begin(), ::toupper);
    if (solver_settings_.find(name)==solver_settings_.end())
    {
        if (name_upper.find("2ND")!=std::string::npos)
            second_order_in_space=true;

        /// Since Legacy NCFVT classes create their own solver instances,
        /// the method here extracts the associated Solver and Setting pointers to be used.
        /// Julian E. Mindel 24-09-2013


        if (name_upper.find("IMPLICIT")!=std::string::npos){
            NodeCenteredFiniteVolumeTransport<dim>* pint=new NodeCenteredFiniteVolumeTransport<dim> (
                        "Model",
                        *this->GetModel(),
                        variables[0].c_str(),
                    variables[1].c_str(),
                    variables[2].c_str(),
                    variables[3].c_str(),
                    variables[4].c_str(),
                    second_order_in_space,
                    second_order_in_time,
                    variables[5].c_str());

            pint->WithLsmGradientLimiter(*this->GetModel());

            SAMG_Solver* sol = pint->GetSolver();
            name_samgsolver_[name]=sol;
            solver_settings_[name]=dynamic_cast<SAMG_Settings*>(sol->GetSolverSettings());
            legacy_FEFV_integrators_[sol]=pint;
        }
        else {
            if (name_upper.find("MASSBASED")==std::string::npos){
                ExplicitNodeCenteredFiniteVolumeTransport<dim,ExplicitStencilProcessor>* pint=new ExplicitNodeCenteredFiniteVolumeTransport<dim,ExplicitStencilProcessor> (
                            "Model",
                            *this->GetModel(),
                            variables[0].c_str(),
                        variables[1].c_str(),
                        variables[2].c_str(),
                        variables[3].c_str(),
                        variables[4].c_str(),
                        second_order_in_space,//second order in space
                        second_order_in_time,
                        variables[5].c_str());

                SAMG_Solver* sol = pint->GetSolver();
                name_samgsolver_[name]=sol;
                solver_settings_[name]=dynamic_cast<SAMG_Settings*>(sol->GetSolverSettings());
                legacy_FEFV_integrators_[sol]=pint;
            }
            else{
                ExplicitMassBasedTransport<dim, MassBasedStencilProcessor>* pint=new ExplicitMassBasedTransport<dim, MassBasedStencilProcessor> (
                            "Model",
                            *this->GetModel(),
                            variables[0].c_str(),
                        variables[1].c_str(),
                        variables[2].c_str(),
                        variables[3].c_str(),
                        variables[4].c_str(),
                        second_order_in_space,//second order in space
                        second_order_in_time);
                SAMG_Solver* sol = pint->GetSolver();
                name_samgsolver_[name]=sol;
                solver_settings_[name]=dynamic_cast<SAMG_Settings*>(sol->GetSolverSettings());
                legacy_FEFV_integrators_[sol]=pint;
            }
        }
    }
}

template<uint32_t dim>
void SimulatorSetup<dim>::QuietIntegrators()
{
    for (auto iti = this->GetLegacy_FE_Integrators().begin(); iti != this->GetLegacy_FE_Integrators().end();iti++){
        iti->second->Verbose(false);
    }
    for (auto iti = this->GetLegacy_FEFV_Integrators().begin(); iti != this->GetLegacy_FEFV_Integrators().end();iti++){
        iti->second->Verbose(false);
    }
}

template <uint32_t dim>
SAMG_Settings* SimulatorSetup<dim>::GetIntegratorSolverSettings(string name)
{
    if (solver_settings_.find(name)==solver_settings_.end()){
        throw csmp::Exception(FATAL_ERROR,"SimulatorSetup<dim>::GetIntegratorSolverSettings"," Did not find a solver settings object for this name of an integrator.\n Did you add the integrator to the SimulatorSetup?");
    }
    return solver_settings_.at(name);
}

/*
Roman, 2014. Deactivate FEFV_Algorithm for reconstruction purposes
template<uint32_t dim>
FEFV_Algorithm<dim>* SimulatorSetup<dim>::GetDESCompatibleIntegrator(SAMG_Solver* samg)
{
    if (des_compatible_integrators_.find(samg)==des_compatible_integrators_.end()){
        cout<<"SimulatorSetup<dim>::GetDESCompatibleIntegrator. Did not find a solver integrator for this SAMG Solver entry."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    return des_compatible_integrators_.at(samg);
}
*/

template<uint32_t dim>
PDE_Integrator<dim,Element>* SimulatorSetup<dim>::GetLegacy_FE_Integrator(SAMG_Solver* samg)
{
    if (legacy_FE_integrators_.find(samg)==legacy_FE_integrators_.end()){
        cout<<"SimulatorSetup<dim>::GetLegacy_FE_Integrator. Did not find a solver integrator for this SAMG Solver entry."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    return legacy_FE_integrators_.at(samg);
}

template<uint32_t dim>
NodeCenteredFiniteVolumeTransport<dim>* SimulatorSetup<dim>::GetLegacy_FEFV_Integrator(SAMG_Solver* samg)
{
    if (legacy_FEFV_integrators_.find(samg)==legacy_FEFV_integrators_.end()){
        cout<<"SimulatorSetup<dim>::GetLegacy_FEFV_Integrator. Did not find a solver integrator for this SAMG Solver entry."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    return legacy_FEFV_integrators_.at(samg);
}

/*
Roman, 2014. Deactivate FEFV_Algorithm for reconstruction purposes
template<uint32_t dim>
FEFV_Algorithm<dim>* SimulatorSetup<dim>::GetDESCompatibleIntegrator(string name)
{
    if (name_samgsolver_.find(name)==name_samgsolver_.end()){
        cout<<"SimulatorSetup<dim>::GetDESCompatibleIntegrator. Did not find an SAMG_Solver associated to this integrator name: "<<name<<endl;
        if (!name_samgsolver_.empty()){
            cout<<" List of inserted integrators."<<endl;
            for (auto nit=name_samgsolver_.begin();nit!=name_samgsolver_.end();nit++)
                cout<<nit->first<<endl;
        }
        else
            cout<<" Name of integrator does not exist in the list."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    if (des_compatible_integrators_.find(name_samgsolver_.at(name))==des_compatible_integrators_.end()){
        cout<<"SimulatorSetup<dim>::GetDESCompatibleIntegrator. Did not find a solver integrator for this SAMG Solver entry."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    return des_compatible_integrators_.at(name_samgsolver_.at(name));
}
*/

template<uint32_t dim>
PDE_Integrator<dim,Element>* SimulatorSetup<dim>::GetLegacy_FE_Integrator(string name)
{
    if (name_samgsolver_.find(name)==name_samgsolver_.end())
    {
        cout<<"SimulatorSetup<dim>::GetLegacy_FE_Integrator. Did not find an SAMG_Solver associated to this integrator name: "<<name<<endl;
        if (!name_samgsolver_.empty()){
            cout<<" List of inserted integrators."<<endl;
            for (auto nit=name_samgsolver_.begin();nit!=name_samgsolver_.end();nit++)
                cout<<nit->first<<endl;
        }
        else
            cout<<" Name of integrator does not exist in the list."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    if (legacy_FE_integrators_.find(name_samgsolver_.at(name))==legacy_FE_integrators_.end()){
        cout<<"SimulatorSetup<dim>::GetLegacy_FE_Integrator. Did not find a solver integrator for this SAMG Solver entry."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    return legacy_FE_integrators_.at(name_samgsolver_.at(name));
}

template<uint32_t dim>
NodeCenteredFiniteVolumeTransport<dim>* SimulatorSetup<dim>::GetLegacy_FEFV_Integrator(string name)
{
    if (name_samgsolver_.find(name)==name_samgsolver_.end()){
        cout<<"SimulatorSetup<dim>::GetLegacy_FEFV_Integrator. Did not find an SAMG_Solver associated to this integrator name: "<<name<<endl;
        if (!name_samgsolver_.empty()){
            cout<<" List of inserted integrators."<<endl;
            for (auto nit=name_samgsolver_.begin();nit!=name_samgsolver_.end();nit++)
                cout<<nit->first<<endl;
        }
        else
            cout<<" Name of integrator does not exist in the list."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    if (legacy_FEFV_integrators_.find(name_samgsolver_.at(name))==legacy_FEFV_integrators_.end()){
        cout<<"SimulatorSetup<dim>::GetLegacy_FEFV_Integrator. Did not find a solver integrator for this SAMG Solver entry."<<endl;
        cout<<" Did you add the integrator to the SimulatorSetup?"<<endl;
    }
    return legacy_FEFV_integrators_.at(name_samgsolver_.at(name));
}
#else
    /// add extra functionality for alternative solver if needed
#endif

template class SimulatorSetup<1U>;
template class SimulatorSetup<2U>;
template class SimulatorSetup<3U>;
}


