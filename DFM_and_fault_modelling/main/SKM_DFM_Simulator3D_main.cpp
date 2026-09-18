// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

// ANALYSIS OF DISCRETE FRACTURE AND (ROCK)MATRIX MODELS
// SKM, Octobre 2013
#include "FracMan_Interface.h"
#include "ANSYS_Model3D.h"
#include "ModelTime.h"
#include "InputDataManager.h"
#include "Standard_IO_Handler.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "ComputationalSettings.h"
#include "GlobalVerbose.h"
#include "convertANSYS_to_CSMP_model.h"
#include "fracturePropertyModelling.h"
#include "pointPropertyMapping.h"
#include "DFM_PermeabilityAndFluxRatioAnalyzer.h"
#include "VTU_Interface.h"

using namespace std;
using namespace csmp;

namespace csmp {
      /// Siroos' implemtation of the velocity-volume averaging technique
      template<uint32_t> void fullTensorPermeability_Durlofsky( const char* model_name );
   }

/**
     2D flow-based upscaling tools including determination of relative permeabilities
     
     @todo:  ANSYS model conversation: add option to deal with different configuration files 
 
     @todo:  integrate Shaho's simulator to carry our relative permeability measurements.
     
     @attention fractures should be significantly smaller than the bouding box of the model; give user some feedback
     
*/
int main()
 {
    cout <<"____________________________________________________________________________________________________\n\n";
    cout <<"\tCSMP++ Discrete Fracture & Matrix (DFM) Simulator, vs. 2025 alpha.\n";
    cout <<"____________________________________________________________________________________________________\n\n";
#ifdef NDEBUG
    cout <<"\n\trunning in RELEASE mode."<< endl;
#else
    cout <<"\n\trunning in DEBUG mode."<< endl;
#endif
    // CURRENT PHYSICAL VARIABLES FILE
    const string  simulation_variable_dataset("CSMP_k_analysis_variables.txt");
   
    cout <<"\nmain: What would you like to do? - simulation options: "<< endl
//         <<"\t(0) Read FRACMAN fracture file and eliminate proximal and overlapping fractures."<< endl
         <<"\t(1) COMPULSORY: read and convert ANSYS model into CSMP native binary format (all other options require this as input)."<< endl
         <<"\t(2) map properties from region-specific point clouds to CSMP model."<< endl

         /** Fracture permeability modeling:
         
            - preprocessing step: failure criteria from in situ stress and single-valued pore pressure (config file)
            - constant apertures and fracture permeabilities from config file
            - opening-mode fracture: constant maximum aperture
            - Cruikshank model for variable apertuer of planar critically stressed fractures
            - Closure aperture (new param, read from config file), but also compute compaction failure (exclude from computational domain?)
                
            - TODO: Olson model of "frozen" fracture aperture, linear and sublinear scaling
            - TODO: shear failure model using JRC data, as developed by Siroos
            - fractures like faults: dilatancy and poro-perm / grainsize relationship
          
            - TODO: create concentric permeability rings for each fracture
         */
         <<"\t(3) calculate failure criteria, fracture apertures, and permeabilities for a given in situ stress, mechanical rock properties and fluid pressure."<< endl
   
         <<"\t(4) analyse diagonal-only tensor k_equiv, qfqm, and velocity spectra."<< endl

         <<"\t(5) analyse full-tensor k_equiv."<< endl // TODO: finalize full tensor permeability option

         <<"\t(6) output 'thickness' and 'permeability' of selected regions to VTU file."<< endl;

// TODO: move multiphase flow methods to a new relative permeability analysis tool
//         <<"\t(6) analyse imbibition relative permeability in 2D using experimentally determined matrix properties (Shaho's simulator)."<< endl
//         <<"\t(7) analyse imbibition relative permeability in 3D with Brooks-Corey matrix properties."<< endl;
   
    int choice(0);
    cin >> choice;
    
    try {
          cout <<"\nCSMP++ DFM_Simulator (3D): Enter name of three-dimensional simulation model: ";
          string model_name("undefined");
          cin >> model_name;
          // TODO: add check for whether file exists and prompt user for different name if not
      
          // Read FRACMAN .fab file and eliminate proximal and overlapping fractures
          if ( choice == 0 ) {
               // TODO: FRACMAN_Interface( model_name.c_str(), false );
               cerr <<"\nmain: Option currently not available.\n";
               return 0;
            }
          // ANSYS input file set + regions + configuration file
          if ( choice == 1 ) {
               return !convertANSYS3D_to_CSMP_model( model_name.c_str(), simulation_variable_dataset.c_str() );
            }
          // point clouds to CSMP model
          else if ( choice == 2 ) {
               cout <<"\nCSMP++ DFM_Simulator: mapping property values from points in specific regions to CSMP binary. ";
               pointPropertyMapping( model_name.c_str() );
               return 0;
            }
          // calculate failure criteria for the fractures
          else if ( choice == 3 ) {
               // descend into calculation option submenu
               fracturePropertyModelling<3>( model_name.c_str() );
               return 0;
            }
          // equivalent permeability and fracture matrix flux ratio
          else if ( choice == 4 ) {
               const bool with_vtk_output(true);
               DFM_PermeabilityAndFluxRatioAnalyzer<3U>  analyzer( model_name.c_str(), with_vtk_output );
               analyzer.DiagonalTensorAnalysis( model_name.c_str(), with_vtk_output );
               return 0;
            }
          // equivalent full tensor permeability and fracture matrix flux ratio
          else if ( choice == 5 ) {
               fullTensorPermeability_Durlofsky<3U>( model_name.c_str() );
               // DFM_PermeabilityAndFluxRatioAnalyzer<3U>  analyzer( model_name.c_str(), verbose );
               // const bool with_vtk_output(true);
               // analyzer.FullTensorAnalysis( model_name.c_str(), with_vtk_output );
               return 0;
            }
          // 'thickness' and 'permeability' output of selected regions
          else if ( choice == 6 ) {
               list<string>  result_vars;
               result_vars.push_back("thickness");
               result_vars.push_back("permeability");
               Model<3U>          model(model_name);
               VTU_Interface<3U>  vtk_output( model, "DFM_modeling");
               cout <<"\nmain: Enter number and names of region where 'thickness' and 'permeability' shall be output to VTK files: ";
               size_t n_regions;
               cin >> n_regions;
               for ( size_t i=0U; i<n_regions; ++i ) {
                    string region;
                    cin >> region;
                    if ( !region.empty() )
                      vtk_output.OutputDataToVTU( (string(model_name) + "-mechanical_characteristics").c_str(), result_vars, region.c_str(), static_cast<int>(0) );
                 }
               return 0;
            }
          else
          throw csmp::Exception( FATAL_ERROR, "main", "choice of model not recognized.");
      }

      // catching all possible standard and csmp::Exceptions
      catch( bad_alloc& ba ) {
           cout <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
        }
      catch( bad_cast& ba ) {
           cout <<"\nbad_cast: Type casting error caused by: "<< ba.what() << endl;
        }
      catch( bad_exception& ba ) {
           cout <<"\nbad_exception: Exception error caused by: "<< ba.what() << endl;
        }
      catch( bad_typeid& ba ) {
           cout <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << endl;
        }
      catch( ios_base::failure& ba ) {
           cout <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
        }
      // standard logic errors
      catch( domain_error& ba ) {
           cout <<"\ndomain_error: Logic error caused by: "<< ba.what() << endl;
        }
      catch( invalid_argument& ba ) {
           cout <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << endl;
        }
      catch( length_error& ba ) {
           cout <<"\nlength_error: Logic error caused by: "<< ba.what() << endl;
         }
      catch( out_of_range& ba ) {
           cout <<"\nout_of_range: Logic error caused by: "<< ba.what() << endl;
        }
      // runtime errors
      catch( overflow_error& ba ) {
           cout <<"\noverflow_error: Runtime error caused by: "<< ba.what() << endl;
        }
      catch( range_error& ba ) {
           cout <<"\nrange_error: Runtime error caused by: "<< ba.what() << endl;
        }
      catch( underflow_error& ba ) {
           cout <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << endl;
        }
      catch( Exception& ba ) {
           cout <<"\nException: Exception raised by: "<< ba.What() << endl;
           cout <<"\nDiagnostics:"<< endl;
           ba.Out();
        }

    cout <<"\nmain: That's it."<< endl;
   
    return 0;
    
 } // end main
