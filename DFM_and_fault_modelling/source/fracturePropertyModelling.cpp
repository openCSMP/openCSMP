//
//  fracturePropertyModelling.cpp
//
//  Created by Stephan Matthai on 1/11/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "MechanicalProperties.h"
#include "FaultFlowPropertyCalculator.h"
#include "FFSA_Aperture.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "ModelTime.h"
#include "VTU_Interface.h"
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
#include "fracturePropertyModelling.h"
#include "compareFloats.h"
#include "geometricCalculations.h"

using namespace std;

/* USER INTERFACE DESIGN CONSIDERATIONS

Implemented approach
--------------------

1. First guess of minimum aperture is provided, set-by-set, via the configuration file

2. For a given state of stress and fluid pressure, an analysis is performed of which
   fractures are critically stressed and what their failure mode is

3. Based on the failure mode, user can choose 3 different flavours of the Cruikshank
   model to compute apertures for the critically stressed fracture segments.
 
   - these apertures get assigned only if they are greater than the original apertures 
     which are taken as lower bounds
   
4. A further option allows to reduce the aperture of fracture segments that experienced
   compressive failure
   
5. After examination of generated aperture distributions, user can impose upper and lower
   bounds on the computed apertures before saving the model to disk.


Revision
--------

- user has no idea what the right aperture should be 
- the in situ stress state may be uncertain
  
  -> User wants to experiment with the tool, creating different aperture
     patterns, keeping the ones that seem geologically plausible

- missing: the Olson & Laubach models which impose linear- or sublinear
  relationships between fracture length and aperture

*/

namespace csmp {

/// reading fault modeling input data
template<uint32_t dim>
void inputFromFile( Model<dim>& model,
                    const char* file_name,
                    set<string>& faults,
                    double& depth, double& rdensity, double& overburden,
                    double& Sv, double& SH, double& Sh, double& trend );
  
/// writing fault modeling input data into text file
void outputToFile( const char* file_name,
                   const set<string>& faults,
                   double depth, double rdensity, double overburden,
                   double Sv, double SH, double Sh, double trend );

template<uint32_t dim>
void imposeUpperLimit( Model<dim>& model, const char* fracture_region, const char* property, double limit );

template<uint32_t dim>
void parallelPlatePermeabilityFromAperture( Model<dim>& model, const char* fracture_region );

/** 
    restricts value range of scalar element property.
*/
template<uint32_t dim>
void imposeUpperLimit( Model<dim>& model, const char* fracture_region, const char* property, double limit )
 {
     Region<dim>&  gref(model.Region(fracture_region));

     const csmp::Index  key(model.Database().StorageKey(property));
     assert( key.place == ELEMENT );
     assert( key.type == SCALAR );
   
     for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
       if ( (*it)->Read(key) > limit )
         (*it)->Store( key, makeScalar( (*it)->Status(key), limit ) );
   
 } // end imposeUpperLimit
 
template void imposeUpperLimit( Model<3>&, const char*, const char*, double );
template void imposeUpperLimit( Model<2>&, const char*, const char*, double );





/**
     FRACTURE PERMEABILITY MODELING

     Driving all options provided:
     
     - Cruishank et al. (2001) aperture model for critically stressed fractures
     
     - dilatation of fractures in shear
     
     - fractures treated like small faults
     
     - compaction of closed fractures due to compressive failure
 
Names of discretized variables used in the fracture modeling.
 
thickness                     h	m	1	1.00E-06	1.00E+03	element
porosity                      PHI	X	1	1.00E-05	1.00E+00	element
failure					FF	X	1	-1.		4.		element   - integer -1 to 4, 0= stable / intact
normal stress				SN	Pa m2	1	-1.00E+08	1.00E+08	element
shear stress				SS	Pa m2	1	-1.00E+08	1.00E+08	element
overburden stress			OS	Pa m2	1	0.		1.00E+08	element
tipline distance			DTL	m	1	0.		1.00E+05	node
centerline distance			DCL	m	1	0.		1.00E+05	node

parameters needed when fractures are treated as faults:

fault size				FS	m	1	1.		1.00E+05	region
slip patch area				AS	m2	1	0.1		1.00E+10	region
effective stress			SSE	Pa m2	1	-1.00E+08	1.00E+08	element
dilatation				DIL	X	1	-1.		1.		element
fracture permeability         KFR	m2	1	1.00E-25	1.00E+00	element
fault breccia permeability		KBR	m2	1	1.00E-25	1.00E+00	element

upscaling parameters obtained at this step (that should be output)
P32
fracture matrix interface area		FMA	m2	1	1.00E-20	1.00E+08	element
fracture porosity			PHIF	X	1	1.00E-20	1.00E+08	element
fracture intensity			FI	1/m	1	0.		1000.		element
block radius				BR	m	1	1.00E-20	1.00E+08	element
fracture matrix flux ratio		QFQM	X	1	1.00E-20	1.00E+08	element

*/
template<uint32_t dim>
void fracturePropertyModelling( const char* model_name )
 {
    // 0. Model initialization from CSMP binary file
    // ______________________________________________________________________________________________

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    double& model_time( ModelTime::Instance().modelTime );
    model_time = 0.;

    // reading model with assigned properties etc. from CSMP native binary format
    string  report = string(model_name) + "-user-choices";
    report +="-";
    report +="fracturePropertyModelling";
    Standard_IO_Handler  stdio( report.c_str() );

    string  fracture_modeling_data(model_name);
    fracture_modeling_data +="-permeability-modeling-data";

    Model<dim>   model(string{model_name});
    set<string>  fractures,fractures_basic_set; // non-unique region for property modeling
    string       fault_name, dummy;
    int          n_fracture_sets(0);


    // 1. checking input properties that are relevant to the property modeling
    // ______________________________________________________________________________________________

    // porosity
    double var_min, var_max;
    model.Database().RangeOf( "porosity", var_min, var_max );
    cout <<"\nfracturePropertyModelling: testing ranges of input variables..."<< endl;
    if ( printRangeOfVariable( model, "porosity", false ) < var_min ) {
         cout <<"\n\tValue(s) below minimum 'porosity' (phi_min="<< var_min <<") detected\n";
         if ( stdio.YesNo("fracturePropertyModelling: replace by minimum value from PropertyDatabase") )
           {
              csmp::Index phi_key = model.Database().StorageKey("porosity");
              Region<dim>&  mref=model.Region("Model");
              for ( auto it=mref.CellsBegin(); it!=mref.CellsEnd(); it++ )
                if ( (*it)->Read(phi_key) < var_min )
                  (*it)->Store( phi_key, makeScalar( (*it)->Status(phi_key),var_min) );
           }
      }
    printRangeOfVariable( model, "shear modulus", false );
    printRangeOfVariable( model, "tensile strength", false );
    printRangeOfVariable( model, "unconfined compressive strength", false );
    printRangeOfVariable( model, "cohesion", false );
    printRangeOfVariable( model, "joint roughness coefficient", false );


    // 2. reading region data, in situ stress, and fluid pressure either from file or from screen
    // __________________________________________________________________________________________

    double depth, rdensity, overburden, Sv, SH, Sh, trend;
    // 2.1 fracture set and stress state input from file
    if ( stdio.YesNo("fracturePropertyModelling: do you have a fracture-set and a stress data file (*-permeability-modeling-data.txt)?") )
      inputFromFile( model, std::string(fracture_modeling_data+".txt").c_str(), fractures, depth, rdensity, overburden, Sv, SH, Sh, trend );
    else {
        cout <<"\n\nfracturePropertyModelling: Enter number and names of fracture/fracture sets to be considered: ";
        cin >> n_fracture_sets;
        for ( int i=0; i<n_fracture_sets; i++ ) {
             cin >> fault_name;
             fractures_basic_set.insert(fault_name);
          }

        for ( auto git=model.UniqueRegionsBegin(); git!=model.UniqueRegionsEnd(); git++ )
           for ( auto it=fractures_basic_set.begin(); it!=fractures_basic_set.end(); ++it )
               if ( (*git).first.find(*it) != string::npos ) {
                    fractures.insert( (*git).first );
                    cout<< "\tregion '"<< (*git).first <<"' was identified."<<endl;
                    break;
                 }

        // 2.2 input the state of stress from the console
        cout <<"\nfracturePropertyModelling: Enter information about in situ stress and the location of measurement.\n";
        cout <<"\n\tEnter subsurface depth of measurement, ambient rock density, and litho-static load: ";
        cin >> depth >> rdensity >> overburden;
        cout <<"\n\tEnter values for Sv, SH, Sh as normalized by litholostatic load: ";
        cin >> Sv >> SH >> Sh;
        cout <<"\n\tEnter the trend of SH max (azimuth=0-180o): ";
        cin >> trend;
     }

    // 2.3 Splitting input regions if they are discontigouos
    bool discontiguous_regions(false);
    for ( auto it=fractures.begin(); it!=fractures.end(); ++it )
      if ( !model.Region((*it).c_str()).IsContiguous() ) {
           cerr <<"\n\tfracturePropertyModelling: discovered discontiguous region: "<< (*it);
           discontiguous_regions = true;
        }
    if ( discontiguous_regions )
      if ( !stdio.YesNo("fracturePropertyModelling: discontiguous regions may compromise aperture estimates; do you want to continue?") )
        throw Exception( ERROR, "fracturePropertyModelling", "user chose not to proceed with discontiguous fracture regions.");

    // 2.4 forming a joint fracture reagion for visualisation
    model.MergeRegions( fractures, "fractures" );
    cout <<"\n\n\nfracturePropertyModelling: the following fracture(s) / fracture sets will be considered in the property modeling:\n\n";
    for ( set<string>::const_iterator i=fractures.begin(); i!=fractures.end(); i++ ) cout << (*i) <<" "; cout <<"\n\n";
   
    outputToFile( std::string( fracture_modeling_data + "-contiguous-regions-revised.txt").c_str(), fractures, depth, rdensity, overburden, Sv, SH, Sh, trend );

    // 2.5 in situ stress state initialization and output to file
    InSituStress stress( depth, rdensity, overburden, Sv, SH, Sh, 0. );
    const double XZ_plane_max_stress_azimuth{ trend }; // default is N-S=0., trend rotation is clockwise in degrees
    if constexpr ( dim == 3 ) {
           stress.StressState().Rotate( 'y', trend );
           stress.Trend( XZ_plane_max_stress_azimuth );
      }
    // in a 2D model, we still want to rotate the azimuth that lies in the XZ plane, but this is achieved by how the 3D stresses are assigned to 2D
    else if constexpr ( dim == 2 ) {
           stress.StressState().Rotate( 'y', trend );
           stress.Trend( XZ_plane_max_stress_azimuth );
      }
    cout <<"\n\n\nfracturePropertyModelling: stress trend (0-360o clockwise looking down): "<< trend << endl;
    string  stress_file_name(model_name);
    stress_file_name +="-stress_tensor";
    // setting the tensor location to lower front corner of the model (x-max, y-min, z-max)
    Point<dim> xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max );
    double max_dimension{1.};
    if constexpr ( dim == 3 ) {
         max_dimension = (std::max(xyz_max[0U]-xyz_min[0U],std::max(xyz_max[1U]-xyz_min[1U],xyz_max[2U]-xyz_min[2U])));
         // deliberately ignoring the depth of stress measurement stored in InSituStress::subsurface_depth_.
         stress.SampleLocation( xyz_max[0U], xyz_min[1U], xyz_max[2U] );
      }
    else if constexpr ( dim == 2 ) {
         max_dimension = ( std::max(xyz_max[0U]-xyz_min[0U], xyz_max[1U]-xyz_min[1U]) );
         stress.SampleLocation( (xyz_max[0U]-xyz_min[0U])/2., depth, (xyz_max[1U]-xyz_min[1U])/2. ); // Y set to 0 for model in XZ plane
      }
    stress.OutStressTensorToVTK( stress_file_name, trend, max_dimension/5. );
   
   
    // 3. Calculating geometrical fracture characteristics
    // ______________________________________________________________________________________________

    // 'tipline distance', 'centerline distance', 'fault size'

    // TODO: Add option to calculate tipline distance by solving Laplace equation with homogeneous boundary conditions

    FaultFlowPropertyCalculator<dim>  calculator;
    cout <<"\n\n\nfracturePropertyModelling: calculating geometrical properties of fractures...\n";
    // initial values are needed here else nan will be printed when range is queried since not all of the model consists of faults
    model.InputPropertyValue( "tipline distance", makeScalar(PLAIN,0.) );
    model.InputPropertyValue( "centerline distance", makeScalar(PLAIN,0.) );
    // using the term fault size here because it is more general than 'fracture'
    // (a fracture that has beeen sheared is considered as a fault)
    model.InputPropertyValue( "fault size", makeScalar(PLAIN,0.) );
    calculator.CalculateTipAndCenterLineDistances( model, fractures );

    // reporting
    list<string>  fracture_vars;
    fracture_vars.push_back("tipline distance");
    fracture_vars.push_back("centerline distance");
    printRangeOfVariable( model, "fractures", "tipline distance" );
    printRangeOfVariable( model, "fractures", "centerline distance" );
    const double max_fault_size = printRangeOfVariable( model, "fault size", true );
    if ( max_fault_size > printModelDimensions( model, false ) / 2. ) {
          csmp_error.Note( WARNING, "fracturePropertyModelling:",
                            "largest fracture is greater than 1/2 of the long axis of model.");
          if ( !stdio.YesNo("fracturePropertyModelling: calculations with models smaller than the REV may lead to erratic values; do you want to proceed?") )
             throw Exception( ERROR, "fracturePropertyModelling", "attempted computation with a model that is too small to capture fracture statistics.");
      }
    list<string>  result_vars;
    result_vars.push_back("thickness");
    result_vars.push_back("permeability");
   
   
    // 4. Evaluating failure criteria for fracture aperture modeling
    // ______________________________________________________________________________________________

    cout <<"\n\n\nfracturePropertyModelling: Entering fracture-aperture calculation in model '" << model_name;
    cout <<"', input variable ranges:\n";
    printRangeOfVariable( model, "fractures", "fluid pressure" );
    if ( dim == 3 && stdio.YesNo("fracturePropertyModelling: Do you want to replace fluid pressure by a depth-based hydrostatic pressure gradient?") ) {
          cout <<"To compute vertical variation of pressure values; enter fluid pressure at model top (Pa) and average fluid density (kg/m3): ";
          double pf_model_top(100325.), fluid_density(1000.);
          cin >> pf_model_top >> fluid_density;
          calculator.CalculateLinearFluidPressure( model, pf_model_top, fluid_density );
          cout <<"\n\n\tRecomputed 'fluid pressure' assuming it is hydrostatic (rho_f=1000kg/m3).\n";
          printRangeOfVariable( model,"fluid pressure", true );
      }
    //  criteria: 1 to 5, for a value of zero, the fracture is stable
    calculator.EvaluateStressesAndFailurePotential( model, fractures, stress );
    // allowing the user to increase fluid pressure to trigger failure
    {
        double fmin{0.}, fmax{0.};
        model.MinMaxOf("failure", fmin, fmax );
        // failure did not occur yet
        if ( approximatelyEqual(fmin,fmax) && approximatelyEqual(fmax,0.) ) {
             if ( stdio.YesNo("fracturePropertyModelling: no failure thus far, do you want to increase fluid pressure?") ) {
                    cout <<"\n"<<"fracturePropertyModelling: Enter pressure (Pa) to add to existing values: ";
                    const csmp::Index pkey = model.Database().StorageKey("fluid pressure");
                    double p_increment{0.};
                    cin >> p_increment;
                    for ( auto& nit : model.Region("Model").NodeVector() )
                      nit->Store( pkey, makeScalar( nit->Status(pkey), nit->Read(pkey) + p_increment ) );
                    cout <<"\n\n\tRecomputed 'fluid pressure' after adding "<< p_increment <<" Pa.\n";
                    printRangeOfVariable( model,"fluid pressure", true );
               }
             calculator.EvaluateStressesAndFailurePotential( model, fractures, stress );
          }
    }

    // TODO: allow incremental rotation of the stress tensor to calculate an ensemble of potential aperture patterns

    cout <<"\n\n\nfracturePropertyModelling: 'failure' and 'effective stress' ranges after dilatation calculation:\n";
    fracture_vars.push_back("normal stress");
    fracture_vars.push_back("shear stress");
    fracture_vars.push_back("effective stress");
    fracture_vars.push_back("overburden stress"); 
    fracture_vars.push_back("failure");
    printRangeOfVariable( model, "fractures", "failure" );
    printRangeOfVariable( model, "fractures", "effective stress" );

    VTU_Interface<dim>  vtk_output( model, "DFM_modeling");
    cout <<"\n\nfracturePropertyModelling: in situ 'stress' is output as 3 vectors representing the principal stresses; scale and invert these for display.\n";
    vtk_output.OutputDataToVTU( (string(model_name) + "-mechanical_characteristics").c_str(), fracture_vars, "fractures", static_cast<int>(0) );

   
    // 5. Computing fracture permeability from aperture and distance from fault tip
    // ______________________________________________________________________________________________

    cout <<"\n\nfracturePropertyModelling: fracture aperture (='thickness') and permeability modelling loop:";
    cout <<"\n\t(suggestion: calculate the aperture of the critically stressed fractures first; then deal with the closure apertures)\n";
      {
         double max_aperture(0.);
         double param_min, param_max;
         model.Region("fractures").MinMaxOf( "failure", param_min, param_max );
         long choice(0);
         while ( (choice=stdio.RecordIntChoice("\tChoose option:\
                                               \n\t\t(1) Elliptic fracture (Cruikshank 91') effective stress- and fracture radius-dependent single-valued maximum aperture\
                                               \n\t\t(2) Elliptic fracture (Cruikshank 91') with aperture variations across fracture plane\
                                               \n\t\t(3) Rectangular fracture (Matthai 06') for elliptic/rectangular fractures open along their side boundaries\
                                               \n\t\t(4) Sheared fracture (Barton 85') shear-displacement and dilation from JRC_mobilised\
                                               \n\t\t(5) Assign lower limit of aperture to fractures undergoing compressional failure\
                                               \n\t\t(6) Apply upper limit of aperture to all fractures\
                                               \n\t\t(7) Treat fractures like faults with internal structure (use to compute k of sheared fractures)\
                                               \n\t\t(8) Re-read configuration file (with fracture and matrix data)\
                                               \n\t\t(9) exit from loop.\n")) < 8 )
           {
              switch( choice ) {
                 case 1: // CRUIKSHANK CONSTANT APERTURE
                     // --------------------------------
                     // aperture is reported via the variable 'thickness'
                     cout <<"\ncalculating fracture-wise constant maximum apertures using an elliptical crack model...\n";
                     cout <<"\ninput variables:\n";
                     printRangeOfVariable( model, "fractures", "failure" );
                     printRangeOfVariable( model, "Youngs modulus", false );
                     printRangeOfVariable( model, "Poissons ratio", false );
                     printRangeOfVariable( model, "fluid pressure", false );
                     cout <<"\naperture distribution BEFORE calculation:\n";
                     printRangeOfVariable( model, "fractures", "thickness" );
                     if ( param_max < 1 ) {
                          csmp_error.Note( INFO, "fracturePropertyModelling:", "none of the fractures was critically stressed; nothing was done" );
                          parallelPlatePermeabilityFromAperture( model, "fractures" );
                       }
                     else {
                          calculator.ConstantFractureApertureFromCruikshankModel( model, fractures );
                          cout <<"\nAperture distribution AFTER calculation:\n";
                          max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                       }
                     if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU?") )
                       vtk_output.OutputDataToVTU( (string(model_name) + "-Cruikshank_constant_aperture").c_str(), result_vars, "fractures", static_cast<int>(0) );
                   break;
                  
                 case 2: // CRUIKSHANK VARIABLE APERTURE
                     // --------------------------------
                     // aperture varies across fracture as a function of distance from tipline (a=0)
                     cout <<"\ncalculating apertures varying across the fracture planes using an elliptical crack model...\n";
                     cout <<"\ninput variables:\n";
                     printRangeOfVariable( model, "fractures", "failure" );
                     printRangeOfVariable( model, "shear modulus", false );
                     printRangeOfVariable( model, "Poissons ratio", false );
                     printRangeOfVariable( model, "fluid pressure", false );
                     cout <<"\nAperture distribution BEFORE calculation:\n";
                     max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                     if ( param_max < 1 ) {
                          csmp_error.Note( INFO, "fracturePropertyModelling:", "none of the fractures was critically stressed; nothing was done" );
                          parallelPlatePermeabilityFromAperture( model, "fractures" );
                       }
                     else {
                          calculator.VariableFractureApertureFromCruikshankModel( model, fractures );
                          cout <<"\naperture distribution AFTER calculation:\n";
                          max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                       }
                     if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU?") )
                       vtk_output.OutputDataToVTU( (string(model_name) + "-Cruikshank_variable_aperture").c_str(), result_vars, "fractures", static_cast<int>(0) );
                   break;

                 case 3: // ELONGATED FRACTURES WITH VARIABLE APERTURE
                     // ----------------------------------------------
                     // aperture varies across fracture as a function of distance from centerline (a=a_max)
                     cout <<"\ncalculating apertures varying across the planes of elongated fractures...\n";
                     cout <<"\ninput variables:\n";
                     printRangeOfVariable( model, "fractures", "failure" );
                     printRangeOfVariable( model, "shear modulus", false );
                     printRangeOfVariable( model, "Poissons ratio", false );
                     printRangeOfVariable( model, "fluid pressure", false );
                     cout <<"\naperture distribution BEFORE calculation:\n";
                     max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                     if ( param_max < 1 ) {
                          csmp_error.Note( INFO, "fracturePropertyModelling:", "none of the fractures was critically stressed; nothing was done" );
                          parallelPlatePermeabilityFromAperture( model, "fractures" );
                       }
                     else {
                          calculator.VariableApertureForElongatedFractures( model, fractures );
                          cout <<"\naperture distribution AFTER calculation:\n";
                          max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                       }
                     if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU?") )
                       vtk_output.OutputDataToVTU( (string(model_name) + "-elongated_fracture_aperture").c_str(), result_vars, "fractures", static_cast<int>(0) );
                   break;
                 case 4: // DILATION DURING SHEAR using Barton et al. 1985
                     // --------------------------------------------------
                     // aperture is uniform across fracture
                     cout <<"\ncalculating fracture aperture inferring shear displacement from far-field stress...\n";
                     cout <<"\ninput variables:\n";
                     printRangeOfVariable( model, "fractures", "failure" );
                     printRangeOfVariable( model, "shear modulus", false );
                     printRangeOfVariable( model, "Poissons ratio", false );
                     printRangeOfVariable( model, "fluid pressure", false );
                     cout <<"\naperture distribution BEFORE calculation:\n";
                     max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                     
                     if constexpr ( dim == 2 ) {
                          // demo_FFSA();
                          fractureApertureFromShearAndNormalStress( model, stress );
                          cout <<"\naperture distribution AFTER calculation:\n";
                          max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                       }
                     else throw csmp::Exception( ERROR, "fracturePropertyModelling", "shear-dilation not implemented in 3D yet" );

                     if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU?") )
                       vtk_output.OutputDataToVTU( (string(model_name) + "-elongated_fracture_aperture").c_str(), result_vars, "fractures", static_cast<int>(0) );
                   break;

                 case 5: // CLOSURE APERTURE TODO: test this compression option
                     // --------------------
                     cout <<"\ncalculating apertures of fracture segments failed under compression (seff >= UCS)...\n";
                     cout <<"\naperture distribution BEFORE calculation:\n";
                     printRangeOfVariable( model, "fractures", "thickness" );
                     if ( param_max >= 0 )
                       csmp_error.Note( INFO, "fracturePropertyModelling:", "none of the fractures failed in compression. Nothing was done." );
                     else {
                         const double minimum_aperture = stdio.RecordChoice("Enter desired aperture for closed fractures (failure=0)");
                         calculator.CompressiveFailureAperture( model, fractures, minimum_aperture );
                         cout <<"\naperture distribution AFTER calculation:\n";
                         max_aperture = printRangeOfVariable( model, "fractures", "thickness" );
                         if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU?") )
                           vtk_output.OutputDataToVTU( (string(model_name) + "-compacted_fracture_aperture").c_str(), result_vars, "fractures", static_cast<int>(0) );
                      }
                   break;
                  
                 case 6: // UPPER LIMIT IMPOSED ON FRACTURE APERTURE
                      {  // ----------------------------------------
                         cout <<"\n\n\nfracturePropertyModelling: current maximum fracture aperture, a = "<< max_aperture;
                         cout <<" m; imposing an upper limit on fracture aperture...\n";
                         const double maximum_aperture = stdio.RecordChoice("\tEnter desired maximum fracture aperture:");
                         imposeUpperLimit( model, "fractures", "thickness", maximum_aperture );
                         parallelPlatePermeabilityFromAperture( model, "fractures" );
                         if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU?") )
                           vtk_output.OutputDataToVTU( (string(model_name) + "-capped_aperture").c_str(), result_vars, "fractures", static_cast<int>(0) );
                      }
                   break;
                  
                 case 7: // FAULTS
                      {  // --------------------------------
                         // fractures are treated like faults
                         cout <<"\n\n\nfracturePropertyModelling: treating fractures as dilating faults ";
                         cout <<"using average mechanical properties of model (only permeability is modeled)...\n";
                         const Region<dim>&  mref = model.Region("Model");
                         assert( model.Database().Placement("Youngs modulus") == REGION );
                         const csmp::Index Y_key  = model.Database().StorageKey("Youngs modulus");
                         const csmp::Index mu_key = model.Database().StorageKey("Poissons ratio");
                         const csmp::Index Ts_key = model.Database().StorageKey("tensile strength");
                         const double  E  = mref.Read(Y_key);
                         const double  mu = mref.Read(mu_key);
                         const double  B  = E / (3. * (1. - mu)); // see Birch (1961)
                         const double  Ts = mref.Read(Ts_key);
                         const double  cfs(0.5); // coefficient of frictional sliding

                         MechanicalProperties  rock_props( E, B, mu, Ts, cfs );

                         const bool recompute_hydrostatic_pf(false);
                         PropertyHandle<dim>  fracture_intensity( model, "fracture intensity", SCALAR, ELEMENT );
                         // assigning constant fracture intensity fractures per meter to fault zone
                         const double f_intensity = stdio.RecordChoice("\tEnter fracture intensity (fracs-per-m) within fault zones:");
                         model.Region("fractures").InputPropertyValue( "fracture intensity", makeScalar(PLAIN,f_intensity) );
                         PropertyHandle<dim>  dilation( model, "dilatation", SCALAR, ELEMENT );
                   
                         calculator.EvaluateFailureAndDilatation( model, fractures, rock_props, stress, recompute_hydrostatic_pf );
                   
                         cout <<"', input parameter ranges:\n";
                         printRangeOfVariable( model, "fractures", "fault size" );
                         printRangeOfVariable( model, "fractures", "fracture intensity" );
                         printRangeOfVariable( model, "fractures", "centerline distance" );
                         printRangeOfVariable( model, "fractures", "dilatation" );
                         cout <<"\npermeability distribution BEFORE calculation:\n";
                         printRangeOfVariable( model, "fractures", "permeability" );
                   
                         PropertyHandle<dim>  frac_perm( model, "fracture permeability", SCALAR, ELEMENT );
                         PropertyHandle<dim>  brec_perm( model, "fault breccia permeability", SCALAR, ELEMENT );
                   
                         const double tortuosity(1.5);
                         const bool compute_halo_permeability(false);
 
                         calculator.PermeabilityFromStressAndDilatation( model, fractures, max_aperture, tortuosity, compute_halo_permeability );

                         cout <<"\npermeability distribution AFTER calculation:\n";
                         printRangeOfVariable( model, "fractures", "permeability" );
                   
                         result_vars.push_back("fracture permeability");
                         result_vars.push_back("fault breccia permeability");
                   
                         if ( stdio.YesNo("fracturePropertyModelling: Do you want to output results to VTU") )
                           vtk_output.OutputDataToVTU( (string(model_name) + "-fractures_treated_as_faults").c_str(), result_vars, "fractures", static_cast<int>(0) );

                         // fault permeability options
                         if ( !stdio.YesNo("fracturePropertyModelling: Continue with new fault permeability? (if not, you can use fracture- or breccia k)?") )
                           {
                              Region<dim>&  fref(model.Region("fractures"));
                              if ( stdio.YesNo("\tDo you want to replace the composite fault permeability with 'fault breccia permeability'?") )
                                fref.CopyReplace("fault breccia permeability","permeability");
                              else if ( stdio.YesNo("\tDo you want to replace the composite fault permeability with 'fracture permeability'?") )
                                fref.CopyReplace("fracture permeability","permeability");
                           }
                      }
                   break;
                  
                 case 8: // RE-READ CONFIGURATION FILE
                   // --------------------------------
                     cout <<"\n\nfracturePropertyModelling: If you partitioned fracture sets their new names in config file may no longer be recognized.";
                     if ( stdio.YesNo("\tDo you still want to re-read the configuration file") )
                        InputDataManager<dim>().ConfigureFromFile( model, model_name, false,
                                                                  true,    // 2) default prop.values
                                                                  true,    // 3) group prop.values
                                                                  false,   // 4) essential box-boundary conditions
                                                                  false ); // 5) essential flags
                   break;
                
                 default:
                   csmp_error.Note( ERROR, "fracturePropertyModelling:", "fracture aperture option was not recognized; only fracture permeability is computed");
                   parallelPlatePermeabilityFromAperture( model, "fractures" );
               }
           }
        
      } // end fracture property modeling
   

    // 6. smooth the result variables if so desired by user
    // ______________________________________________________________________________________________

    if ( stdio.YesNo("fracturePropertyModelling: Do you want to smooth the computed permeability values?") ) {
         // using fluid pressure1 as a temporary variable
         smoothElementVariable( model, "fractures", "permeability", "fluid pressure1", 1 );
         printRangeOfVariable( model, "fractures", "permeability" );
         vtk_output.OutputDataToVTU( (string(model_name) + "-properties1_k_smoothed").c_str(), result_vars, "fractures", static_cast<int>(0) );
      }

    // 7. scale the result variable if so desired
    // ______________________________________________________________________________________________

    if ( stdio.YesNo("fracturePropertyModelling: Do you want to scale the computed permeability values?") ) {
         double scale_factor = stdio.RecordChoice("\t\tEnter scale factor: ");
         csmp::Index k_key     = model.Database().StorageKey("permeability");
         Region<dim>& fref      = model.Region("fractures");
         ScalarVariable  sc;
         for ( auto it=fref.CellsBegin(); it!=fref.CellsEnd(); ++it ) {
              (*it)->Read( k_key, sc );
              (*it)->Store( k_key, sc *= scale_factor );
           }
         printRangeOfVariable( model, "fractures", "permeability" );
         vtk_output.OutputDataToVTU( (string(model_name) + "-properties2_k_scaled").c_str(), result_vars, "fractures", static_cast<int>(0) );
      }

   
    // 8. bracket the result variable if so desired
    // ______________________________________________________________________________________________

    printRangeOfVariable( model, "fractures", "permeability" );
    if ( stdio.YesNo("fracturePropertyModelling: Do you want to manually impose upper- and lower limits on fracture permeability?") ) {
         double k_max    = stdio.RecordChoice("\t\tEnter desired upper pemeability limit (m2): ");
         double k_min    = stdio.RecordChoice("\t\tEnter desired lower pemeability limit (m2): ");
         csmp::Index k_key = model.Database().StorageKey("permeability");
         Region<dim>& fref  = model.Region("fractures");
         ScalarVariable  sc;
         for ( auto it=fref.CellsBegin(); it!=fref.CellsEnd(); ++it ) {
              (*it)->Read( k_key, sc );
              if      ( sc() < k_min ) sc = k_min;
              else if ( sc() > k_max ) sc = k_max;
              (*it)->Store( k_key, sc );
           }
         printRangeOfVariable( model, "fractures", "permeability" );
         vtk_output.OutputDataToVTU( (string(model_name) + "-properties3_k_limited").c_str(), result_vars, "fractures", static_cast<int>(0) );
      }

   
    // 9. saving model to disk
    // ______________________________________________________________________________________________

    if ( stdio.YesNo("fracturePropertyModelling: Do you want to save the computed 'thickness'and 'permeability' values overwriting the former ones?") )
      model.OutputToBinaryFile( model_name );
   
    stdio.Out();
    cout <<"\n\nfracturePropertyModelling: That's it!"<< endl;
   
 } // end fracturePropertyModelling
 
template void fracturePropertyModelling<3>( const char* );
template void fracturePropertyModelling<2>( const char* );










/** 
    for all surface elements, the parallel plate permeability is calculated from thickness.
*/
template<uint32_t dim>
void parallelPlatePermeabilityFromAperture( Model<dim>& model, const char* fracture_region )
 {
     Region<dim>&  gref(model.Region(fracture_region));

    const csmp::Index  thi_key(model.Database().StorageKey("thickness"));
    const csmp::Index  k_key(model.Database().StorageKey("permeability"));
   
     for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
       if ( (*it)->IsSurface() ) {
            double aperture = (*it)->Read( thi_key );
            (*it)->Store( k_key, makeScalar( (*it)->Status(k_key),  (aperture*aperture)/12. ) );
         }
   
 } // end parallelPlatePermeabilityFromSurfaceThickness

template void parallelPlatePermeabilityFromAperture( Model<3>&, const char* );
template void parallelPlatePermeabilityFromAperture( Model<2>&, const char* );





/// helper function: that returns smallest angle of line to X-axis [0,90]
static double lineAngleToXAxisFromNormal(const Point<2>& n)
{
    double alpha = std::abs(std::atan2(n[0], n[1])) * 180. / std::numbers::pi;
    // range folding
    if ( alpha > 90. ) alpha = 180. - alpha;
    
    return alpha;
}

/**
      Fracture-by-fracture calculation of input parameters to FFSA algorithm of M. Liem, based on Barton's model/
      The computed fracture parameters are collected into FFSA_InputParameters vector.
*/
void fractureApertureFromShearAndNormalStress( Model<2>& model, const InSituStress& stress )
 {
     // 0. selecting unique fracture regions identified by name-string and numbers attached to them
     // -------------------------------------------------------------------------------------------
     set<string> fracture_regions;
     for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit )
       // distinguish fracture regions (=1D line element check is omitted)
       if ( (contains((*rit).first, "frac") ||
             contains((*rit).first, "FRAC") ||
             contains((*rit).first, "Frac")) && containsDigit((*rit).first) )
        {
           fracture_regions.insert( (*rit).first );
        }

     // 1. collecting data from unique fracture regions identified by name-string and numbers attached to them
     // ------------------------------------------------------------------------------------------------------
     FFSA_InputParameters ffsa_params( fracture_regions.size() );

     long N_fracture{0}; // fracture number
     for ( const auto& region : fracture_regions ) {
           Region<2>& frac = model.Region( region );

           // 1. Find unit normal to average fracture plane and calculating angle between
           //    fracture and x-axis
           vector<Point<2>> fracture_pts; fracture_pts.reserve( frac.Nodes() );
           for ( const auto& node : frac.NodeVector() ) fracture_pts.push_back( node->Coordinate() );
           const Point<2> unrml = normalToAveragePlaneThroughPointCloud( fracture_pts );
           frac.Store( model.Database().StorageKey("fracture normal"), makeVector(ANY,ANY,unrml[0],unrml[1]) );
           // angle between fracture plane normal and x-axis
           ffsa_params.alpha[N_fracture] = lineAngleToXAxisFromNormal( unrml );
 
           // 2. initialising far-field stress analysis parameters
           ffsa_params.L[N_fracture]         = frac.Volume(); ///< = length of a line element region
           ffsa_params.p_f[N_fracture]       = frac.Average("fluid pressure") * 1.0e-6; // converting into MPa
           ffsa_params.JRC[N_fracture]       = frac.Read( model.Database().StorageKey("joint roughness coefficient") );     ///< JRC 0-20 [-]
           ffsa_params.sigma_c[N_fracture]   = frac.Read( model.Database().StorageKey("unconfined compressive strength") ) * 1.0e-6; ///< UCS [MPa]
           // NB: Barton and co-workers treat JCS as a degraded or scale-modified UCS to account for:
           // weathering of joint walls, microcracking, gouge infill, asperity damage scale effects; here I set it to 0.9 UCS
           ffsa_params.JCS[N_fracture]       = 0.9 * ffsa_params.sigma_c[N_fracture]; ///< Joint wall compressive strength, JCS [MPa]
           // slope of the stress closure curve at very low stress (~joint compressibility = elastic compliance of asperities)
           // (TODO: correlate K_ni with asperity height, hr: kn0 proportional to UCS / hr)
           ffsa_params.K_ni[N_fracture]      = 1000.;  ///< Initial normal stiffness [MPa/mm] guestimate
           // delta_n,max approx (0.6-1.0) hr ) = closure aperture (converted to mm); substituted by default aperture(=thickness) of line elements
           ffsa_params.vm_factor[N_fracture] = 0.7 * frac.Average("thickness") * 1.0e3;   ///< Factor for maximum possible joint closure [-]
           // Residual friction angle (of the fractured rock) should be estimater from from direct shear tests (with roughness removed)
           // hard crystalline rocks (30-35°), sandstone (28-32°), limestone(25-30°), shale and other weak rocks (20-25°)
           ffsa_params.phi_r[N_fracture]     = 25.;       ///< Residual friction angle [°]
           ffsa_params.E_mod[N_fracture]     = frac.Read( model.Database().StorageKey("Youngs modulus") ) * 1.0e-6; ///< E-modulus of fractured rock [MPa]
           ffsa_params.nu[N_fracture]        = frac.Read( model.Database().StorageKey("Poissons ratio") );          ///< Poisson's ratio of fractured rock [-]
           // should be proportional to JRC increasing with roughness up to ~4
           ffsa_params.C_g[N_fracture]       = 1.1; ///< Proportionality between displacement and shear stress [-]
            
           N_fracture++;
        }

     // 2. computing fracture apertures
     // -------------------------------
     // 2.1 assigning far-field parameters (stresses in MPa)
     ffsa_params.sigma_H = stress.SH() * stress.P_conf() * 1.0e-6;
     ffsa_params.sigma_h = stress.Sh() * stress.P_conf() * 1.0e-6;
     ffsa_params.beta    = stress.Trend();
     // 2.2 assigning method qualifiers
     ffsa_params.method_displacement = "mob"; // Iterative search
     ffsa_params.method_dilation = "integrate_pos"; // Integration, positive increments only
     ffsa_params.method_peak_displacement = "Asadollahi";
     ffsa_params.method_sigma_eff = "end"; // sigma_EFF = sigma_eff
     // optional M parameter (testing the optional field logic)
     ffsa_params.M.resize(ffsa_params.N_fractures); ffsa_params.M.fill(numeric_limits<double>::quiet_NaN());
     ffsa_params.M(0) = 0.8; // Specify M for frac 1

     FFSA_FractureAperture calculator;
     auto [results_matrix, debug_data] = calculator.ApertureFromFarFieldStress( ffsa_params );

     
     // 3. writing aperture back to the model using the parameter 'thickness'
     // ---------------------------------------------------------------------
     cout << "=================================================================\n";
     cout << "FFSA Aperture Calculation Results (N_frac = " << fracture_regions.size() << ")\n";
     cout << "=================================================================\n";
     cout << "Units: L [m], Stress [MPa], Displacement [mm]\n\n";
     // (a negative mobilised JRC indicates compaction failure)
     for ( long i = 0; i < static_cast<long>(fracture_regions.size()); ++i ) {
           cout << "--- Test Fracture " << i + 1 << " ---\n";
           cout << "  Calculated Stresses: sigma_n = " << debug_data.sigma_n(i)
                     << " MPa, sigma_s = " << debug_data.sigma_s(i)
                     << " MPa, sigma_eff = " << debug_data.sigma_eff(i) << " MPa\n";
           cout << "  Initial Aperture (a_0): " << debug_data.a_0(i) << " mm\n";
           cout << "  Normal Closure (delta_n): " << debug_data.delta_n(i) << " mm\n";
           cout << "  Shear Displacement (delta_s): " << debug_data.delta_s(i) << " mm\n";
           cout << "  Shear Dilation (delta_d): " << debug_data.delta_d(i) << " mm\n";
           cout << "  Final Aperture (a): " << debug_data.a(i) << " mm\n";
           cout << "  Mobilized JRC: " << debug_data.JRC_mob(i) << " (-)\n";
           cout << "  Mobilized Dilation Angle: " << debug_data.phi_d_mob(i) << " deg\n";
           cout << "  Permeability Factor (k): " << debug_data.k(i) << " m^2\n\n";
        }
        
     // TODO: (maybe use the deviation of element from mean fracture direction to adjust dilation/opening)
     N_fracture = 0ul; // fracture number
     const csmp::Index thi_key  = model.Database().StorageKey("thickness");
     const csmp::Index fail_key = model.Database().StorageKey("failure");
     for ( const auto& region : fracture_regions ) {
           // only if a closure, shear displacement or dilatation occurred, results are written to model
           if ( debug_data.delta_n(N_fracture) != 0. || debug_data.delta_s(N_fracture) != 0. || debug_data.delta_d(N_fracture) != 0. ) {
                Region<2>& frac = model.Region( region );
                // fracture aperture (converting from mm to m)
                auto aperture = debug_data.a(N_fracture) * 1.0e-3;
                // assigned to all elements in the region
                for ( auto& it : frac.CellVector() ) it->Store( thi_key, makeScalar(it->Status(thi_key),aperture) );
                // recording the failure regime
                // compression
                if ( debug_data.delta_n(N_fracture) > 0. && debug_data.phi_d_mob(N_fracture) < 0. )
                  for ( auto& it : frac.CellVector() ) it->Store( fail_key, makeScalar(it->Status(fail_key),-1.) );
                // shear reactivitation (frictional sliding)
                if ( debug_data.delta_s(N_fracture) > 0. && debug_data.delta_d(N_fracture) > 0. )
                  for ( auto& it : frac.CellVector() ) it->Store( fail_key, makeScalar(it->Status(fail_key),1.) );
                
             }
           N_fracture++;
        }
     
 } // end fractureApertureFromShearAndNormalStress





} // end csmp
