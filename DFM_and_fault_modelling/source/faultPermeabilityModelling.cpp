//
//  faultPermeabilityModelling.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/19/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#include "faultPermeabilityModelling.h"
#include "MechanicalProperties.h"
#include "Model.h"
#include "FaultFlowPropertyCalculator.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "ModelTime.h"
#include "VTU_Interface.h"
#include "SKUA_Interface.h"
#include "Standard_IO_Handler.h"

using namespace std;

namespace csmp {

/**  inputFromFile

     Reads user input from text file '*-permeability-modeling-data.txt' into its argument variables.
     File lists variables sequentially starting with header line, followed by number of faults or fractures,
     fault name strings, and finally stress relevant data.
     
     @attention this function also searches the model for subregions with the same name 
     but an appended number. These regions will also be included. 
     Naming conventions need to be adapted if this undesired.
     
*/
template<uint32_t dim>
void inputFromFile( Model<dim>& model,
                    const char* file_name,
                    set<string>& faults,
                    double& depth, double& rdensity, double& overburden,
                    double& Sv, double& SH, double& Sh, double& trend )
 {
    assert( file_name != NULL );

    ifstream ifs(file_name);
    if ( !ifs.is_open() )
      throw csmp::Exception( ERROR, "inputFromFile",
                             file_name,   "file specifying permeability-model input variables is missing");
    if ( !faults.empty() )
      faults.clear();
   
    set<string> faults_basic_set;

    // reading header line printing it to screen and swallowing empty line thereafter
    char text[256];
    ifs.getline( text, 256 );
    cout <<"\ninputFromFile: file header: "<< text << endl;
    ifs.getline( text, 256 );
   
    int n_faults(0);
    ifs >> n_faults;
    assert( n_faults > 0 );
    assert( n_faults < 10000 );
    string  fault_name;
    for ( int n=0; n<n_faults; n++ ) {
         ifs >> fault_name;
         if ( !fault_name.empty() ) faults_basic_set.insert(fault_name);
         else
         throw csmp::Exception( ERROR, "inputFromFile:", "encountered empty region name.");
      }
   
    bool first_call(true);
    for ( auto git=model.UniqueRegionsBegin(); git!=model.UniqueRegionsEnd(); git++ )
      {
         for ( auto it=faults_basic_set.begin(); it!=faults_basic_set.end(); ++it )
             if ( (*git).first.find(*it) != string::npos )
             {
                 faults.insert( (*git).first );
                 if ( first_call ) {
                      cout << "\n\tregion(s) incorporated into the input list: ";
                      first_call = false;
                   }
                 cout << (*git).first <<" ";
                 break;
             }
      }
    cout <<"\n\n";

    // data required for initializng the stress state
    ifs >> depth >> rdensity >> overburden >> Sv >> SH >> Sh >> trend;
    ifs.close();

    cout <<"\ninputFromFile: region names and stress data stored in '"<< file_name <<"' read successfully."<< endl;
   
 } // end inputFromFile

template void inputFromFile( Model<3>&, const char*, set<string>&, double&, double&, double&,
                             double&, double&, double&, double& );

template void inputFromFile( Model<2>&, const char*, set<string>&, double&, double&, double&,
                             double&, double&, double&, double& );




/// writing fault modeling input data into text file
void outputToFile( const char* file_name,
                   const set<string>& faults,
                   double depth, double rdensity, double overburden,
                   double Sv, double SH, double Sh, double trend )
 {
    assert( file_name != NULL );

    ofstream ofs(file_name);
    if ( !ofs.is_open() )
      throw csmp::Exception( ERROR, "inputToFile",
                             file_name,   "file specifying permeability-model input variables is missing");

    ofs <<"'"<< file_name <<"' fault regions to be included into fault property modeling.\n\n";

    ofs << faults.size() <<" ";
    for ( auto it=faults.begin(); it!=faults.end(); ++it )
      ofs << (*it) <<" ";

    // data required for initializng the stress state
    ofs <<" "<< depth <<" "<< rdensity <<" "<< overburden <<" "<< Sv <<" "<< SH <<" "<< Sh <<" "<< trend;
    ofs <<"\n";
    ofs.close();
   
    cout <<"\ninputToFile: '"<< file_name <<"' written successfully."<< endl;
   
 } // end outputToFile





/**  =========================================================================================

     FAULT PERMEABILITY MODELING

     Driving all options provided:
 
     - contouring distance on fault planes
     - computation of shear and fault normal stress followed by evaluation of failure criteria
     - capped dilatiation of critically stressed parts of the fault planes
     
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    double& model_time( ModelTime::Instance().modelTime );
    model_time = 0.;

    =========================================================================================
    TODO: use multiphysics variable set
*/
void faultPermeabilityModelling( const char* model_name )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    string  report(model_name);
    report +="-";
    report +="fault_permeability_modeling";
    Standard_IO_Handler  stdio( report.c_str() );

    string  fault_permeability_modeling_data(model_name);
    fault_permeability_modeling_data +="-permeability-modeling-data.txt";

    // 1. reading model with assigned properties etc. from CSMP native binary format
    Model<3U>    model(string{model_name});
    set<string>  faults,faults_basic_set;
    string       fault_name, dummy;
    int          n_faults(0);

    if ( !model.Database().IsDefined("displacement") )
      csmp_error.Note( FATAL_ERROR, "faultPermeabilityModelling",
                        "current set of variables does not support geomechanical calculations; use 'ReservoirMultiPhysics-variables.txt' in stead.");

    // 1.1 checking input properties that are relevant to the property modeling
    // porosity
    double var_min, var_max;
    model.Database().RangeOf( "porosity", var_min, var_max );
    cout <<"\nfault_permeability_modeling: testing ranges of input variables..."<< endl;
    if ( printRangeOfVariable( model, "porosity", false ) < var_min ) {
         cout <<"\n\tValue(s) below minimum 'porosity' (phi_min="<< var_min <<") detected\n";
         if ( stdio.YesNo("faultPermeabilityModelling: replace by minimum value from PropertyDatabase") )
           {
              csmp::Index phi_key = model.Database().StorageKey("porosity");
              Region<3U>&  mref=model.Region("Model");
              for ( auto it=mref.CellsBegin(); it!=mref.CellsEnd(); it++ )
                if ( (*it)->Read(phi_key) < var_min )
                  (*it)->Store( phi_key, makeScalar( (*it)->Status(phi_key),var_min) );
           }
      }

    // 2. getting the faults in the flow domain and other necessary data from screen or file
    double depth, rdensity, overburden, Sv, SH, Sh, trend;
    if ( stdio.YesNo("faultPermeabilityModelling: do you have a fault and stress data file (*.-permeability-modeling-data.txt)") )
      inputFromFile( model, fault_permeability_modeling_data.c_str(), faults, depth, rdensity, overburden, Sv, SH, Sh, trend );
    else {
        cout <<"\n\n"<<"faultPermeabilityModelling: Enter number and names of faults to be considered: ";
        cin >> n_faults;
        for ( int i=0; i<n_faults; i++ ) {
             cin >> fault_name;
             faults_basic_set.insert(fault_name);
          }

        for ( map<string,csmp::Region<3U> >::const_iterator
              git=model.UniqueRegionsBegin(); git!=model.UniqueRegionsEnd(); git++ )
          {
             for ( set<string>::const_iterator it=faults_basic_set.begin(); it!=faults_basic_set.end(); ++it )
                 if ( (*git).first.find( *it )!=string::npos)
                 {
                     faults.insert( (*git).first );
                     cout<< "Fault Region( "<< (*git).first <<" ) was identified."<<endl;
                     break;
                 }
          }

        // 3. Establishing the state of stress
        cout <<"\nfaultPermeabilityModelling: Enter information about in situ stress and the location of measurement.\n";
        cout <<"\n\tEnter subsurface depth of measurement, ambient rock density, and litholostatic load: ";
        cin >> depth >> rdensity >> overburden;
        cout <<"\n\tEnter values for Sv, SH, Sh as normalized by litholostatic load: ";
        cin >> Sv >> SH >> Sh;
        cout <<"\n\tEnter the trend of SH max (azimuth=0-180o): ";
        cin >> trend;
     }
   
    outputToFile( fault_permeability_modeling_data.c_str(), faults, depth, rdensity, overburden, Sv, SH, Sh, trend );

    InSituStress stress( depth, rdensity, overburden, Sv, SH, Sh, trend );
   
    cout <<"\n\n\nfaultPermeabilityModelling: stress trend (0-360o): "<< trend << endl;
    string  stress_file_name(model_name);
    stress_file_name +="-stress_tensor";
    Point<3U> xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max );
    const double max_dimension(std::max(xyz_max[0U]-xyz_min[0U],std::max(xyz_max[1U]-xyz_min[1U],xyz_max[2U]-xyz_min[2U])));
    stress.OutStressTensorToVTK( stress_file_name, trend, max_dimension/5. );
   
    // 4. Calculating normalized distances from fault centerlines to fault tips
    model.MergeRegions( faults, "faults" );
    cout <<"\n\n\nfaultPermeabilityModelling: the following faults will be considered in the property modeling:\n\n";
    for ( set<string>::const_iterator i=faults.begin(); i!=faults.end(); i++ ) cout << (*i) <<" "; cout <<"\n\n";
 
    FaultFlowPropertyCalculator<3>  calculator;
    cout <<"\n\n\nfaultPermeabilityModelling: calculating geometrical properties of faults...\n";
    // initial values are needed here else nan will be printed when range is queried since not all of the model consists of faults
    model.InputPropertyValue( "tipline distance", makeScalar(PLAIN,0.) );
    model.InputPropertyValue( "centerline distance", makeScalar(PLAIN,0.) );
    model.InputPropertyValue( "fault size", makeScalar(PLAIN,0.) );
    calculator.CalculateTipAndCenterLineDistances( model, faults );
    printRangeOfVariable( model, "faults", "tipline distance" );
    printRangeOfVariable( model, "faults", "centerline distance" );
    printRangeOfVariable( model, "fault size" );
    printRangeOfVariable( model, "porosity" );
    printRangeOfVariable( model, "faults", "porosity" );
   
    // 5. Estimating fault dilatancy
    cout <<"\n\n\nfaultPermeabilityModelling: fault dilatation calculation in model '" << model_name;
    cout <<"', input variable ranges:\n";
    printRangeOfVariable( model, "faults", "porosity" );
    printRangeOfVariable( model, "faults", "permeability" );
    printRangeOfVariable( model, "faults", "fluid pressure" );
 
    // basement treated as granite:   E       B       Poisson's  Ts     coeff. sliding friction
    MechanicalProperties  rock_props( 30.0e9, 50.0e9, 0.2,       7.7e6, 0.6 );
    rock_props.UCS_   = 5e7; // unconfined compressive strength of shattered granite (low value)
    rock_props.alpha_ = 30.; // angle of internal friction = angle at which Riedel 1 fractures branch from the main fault trace
  
    calculator.EvaluateFailureAndDilatation( model, faults, rock_props, stress );
    cout <<"\n\n\nfaultPermeabilityModelling: 'failure' and 'dilatation' ranges after dilatation calculation:\n";
    printRangeOfVariable( model, "faults", "failure" );
    printRangeOfVariable( model, "faults", "dilatation" );
   
    // 6. Computing fault permeability from dilatancy and distance from fault tip
    double max_aperture = stdio.RecordChoice("\n\n\nfault_permeability_modeling: Enter a limit on fracture aperture: ");
   
    cout <<"\n\n\nfaultPermeabilityModelling: fault permeability calculation in model '" << model_name;  
    cout <<"', input value ranges:\n";
    printRangeOfVariable( model, "faults", "permeability" );
    printRangeOfVariable( model, "faults", "fault size" );
    printRangeOfVariable( model, "faults", "fracture intensity" );
    printRangeOfVariable( model, "faults", "centerline distance" );
    printRangeOfVariable( model, "faults", "dilatation" );
    const bool compute_halo_permeability = stdio.YesNo("fault_permeability_modeling: Do you want to calculate a fault permeability halo array for SKUA");
    const double tortuosity(1.5);
    calculator.PermeabilityFromStressAndDilatation( model, faults, max_aperture, tortuosity, compute_halo_permeability );

    // 6b. giving the user a choice how to proceed with the results from the permeability modeling
    if ( !stdio.YesNo("faultPermeabilityModelling: Continue with the new composite fault permeability? (if not, you can use either the fracture- or breccia k).") ) {
         Region<3>&  fref(model.Region("faults"));
         if ( stdio.YesNo("\tDo you want to replace the composite fault permeability with the breccia permeability") )
           fref.CopyReplace("fault breccia permeability","permeability");
         else if ( stdio.YesNo("\tDo you want to replace the composite fault permeability with the fracture permeability") )
           fref.CopyReplace("fracture permeability","permeability");
      }
      
    cout <<"\n\n\nfaultPermeabilityModelling: 'permeability' range after fault permeability calculation:\n";
    printRangeOfVariable( model, "faults", "permeability" );
   
    cout <<"\n\n\nfaultPermeabilityModelling: completed successfully."<< endl;
   

    // 7. Output of newly computed parameters to VTU
    list<string>  fault_vars;
    fault_vars.push_back("normal stress");
    fault_vars.push_back("shear stress");
    fault_vars.push_back("effective stress");
    fault_vars.push_back("overburden stress");
    fault_vars.push_back("failure");
    fault_vars.push_back("dilatation");
    fault_vars.push_back("permeability");
    fault_vars.push_back("fracture permeability");
    fault_vars.push_back("fault breccia permeability");
    fault_vars.push_back("tipline distance");
    fault_vars.push_back("centerline distance");
    // Region properties cannot be output to VTU
    // fault_vars.push_back("fault size");
    // fault_vars.push_back("slip patch area");
   
    VTU_Interface<3U>  vtk_output( model, "fault-property-modeling");
    vtk_output.OutputDataToVTU( (string(model_name) + "-properties"), fault_vars, string("faults"), 0 );
 
   
    // 8. smooth the result variables if so desired by user
    if ( stdio.YesNo("faultPermeabilityModelling: Do you want to smooth the computed permeability values") ) {
         smoothElementVariable( model, "faults", "permeability", "dummy node", 1 );
         printRangeOfVariable( model, "faults", "permeability" );
         vtk_output.OutputDataToVTU( (string(model_name) + "-properties1_k_smoothed"), fault_vars, string("faults"), 0 );
      }

    // 9. scale the result variable if so desired
    if ( stdio.YesNo("faultPermeabilityModelling: Do you want to scale the computed permeability values") ) {
         double scale_factor = stdio.RecordChoice("\t\tEnter scale factor: ");
         csmp::Index k_key     = model.Database().StorageKey("permeability");
         Region<3U>& fref      = model.Region("faults");
         ScalarVariable  sc;
         for ( auto it=fref.CellsBegin(); it!=fref.CellsEnd(); ++it ) {
              (*it)->Read( k_key, sc );
              (*it)->Store( k_key, sc *= scale_factor );
           }
         printRangeOfVariable( model, "faults", "permeability" );
         vtk_output.OutputDataToVTU( (string(model_name) + "-properties2_k_scaled"), fault_vars, string("faults"), 0 );
      }
   
    // 10. scale the result variable if so desired
    printRangeOfVariable( model, "faults", "permeability" );
    if ( stdio.YesNo("faultPermeabilityModelling: Do you want to manually impose upper- and lower limits on fault permeability") ) {
         double k_max    = stdio.RecordChoice("\t\tEnter desired upper pemeability limit (m2): ");
         double k_min    = stdio.RecordChoice("\t\tEnter desired lower pemeability limit (m2): ");
         csmp::Index k_key = model.Database().StorageKey("permeability");
         Region<3U>& fref  = model.Region("faults");
         ScalarVariable  sc;
         for ( auto it=fref.CellsBegin(); it!=fref.CellsEnd(); ++it ) {
              (*it)->Read( k_key, sc );
              if      ( sc() < k_min ) sc = k_min;
              else if ( sc() > k_max ) sc = k_max;
              (*it)->Store( k_key, sc );
           }
         printRangeOfVariable( model, "faults", "permeability" );
         vtk_output.OutputDataToVTU( (string(model_name) + "-properties3_k_limited"), fault_vars, string("faults"), 0 );
      }

    // 11. Output of newly computed permeability data to SKUA
    if ( stdio.YesNo("faultPermeabilityModelling: Do you want to output permeability values as point data to SKUA") ) {
         cout <<"\n\n\nfaultPermeabilityModelling: output of permeability data to SKUA point cloud."<< endl;
         // attaching to element array variable with 6 entries
         // PropertyHandle<3U>  fault_zone_perm( model, "fault zone permeability", ARRAY, ELEMENT );
         SKUA_Interface      SKUA_output;
         set<string>  scalar_elmt_vars;
         scalar_elmt_vars.insert("permeability");
         scalar_elmt_vars.insert("shear stress");
         scalar_elmt_vars.insert("effective stress");
         scalar_elmt_vars.insert("failure");
         SKUA_output.VariablesToPointCloud( model, faults, (string(model_name) + "-SKUA_point_cloud.txt").c_str(), scalar_elmt_vars );
         SKUA_output.SurfaceArrayVariableToPointCloud( model, faults, (string(model_name) +
                                                     "-perm-SKUA_laterally_extended_point_cloud.txt").c_str(), "fault zone permeability" );
      }
   
    // 12. saving model to disk
    if ( stdio.YesNo("faultPermeabilityModelling: Do you want to overwrite 'permeability' etc. with the newly computed values?") )
      model.OutputToBinaryFile( model_name );
   
    stdio.Out();
    cout <<"\n\nfaultPermeabilityModelling: That's it!"<< endl;
   
 } // emd fault_permeability_modeling

} // end csmp
