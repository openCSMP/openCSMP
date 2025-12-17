// ANALYSIS OF DISCRETE FRACTURE AND (ROCK)MATRIX MODELS
// SKM, December 2025
#include "ANSYS_Model2D.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "ComputationalSettings.h"
#include "GlobalVerbose.h"
#include "ModelTime.h"
#include "Standard_IO_Handler.h"
#include "CSMP_highLevelUtilities.h"
#include "InputDataManager.h"
#include "fracturePropertyModelling.h"
#include "pointPropertyMapping.h"
#include "DFM_PermeabilityAndFluxRatioAnalyzer.h"
#include "VTU_Interface.h"

using namespace std;
using namespace csmp;

namespace csmp {
      /// Siroos' implementation of Durlofsky's velocity-volume averaging technique
      template<uint32_t> void fullTensorPermeability_Durlofsky( const char* model_name );
      
      /// check
      bool checkNeighborConnectivityOfElementsInFractureRegion( const Model<2>& model );
      // TODO: write method that checks modelling assumptions (dim-1 fracs), meaningful input prop ranges,
      // TODO: write method that checks that binary model is there and intact
   }

/**
     2D flow-based upscaling tools including determination of relative permeabilities
     
     @todo:  ANSYS model conversation: add option to deal with different configuration files 
 
     @todo:  integrate Shaho's simulator to carry our relative permeability measurements.
     
     @attention fractures should be significantly smaller than the bouding box of the model; give user some feedback
     
*/
int main()
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    cout <<"____________________________________________________________________________________________________\n\n";
    cout <<"\tCSMP++ Discrete Fracture & Matrix (DFM) Simulator, 2D version 2025 alpha.\n";
    cout <<"____________________________________________________________________________________________________\n\n";
#ifdef NDEBUG
    cout <<"\n\trunning in RELEASE mode."<< endl;
#else
    cout <<"\n\trunning in DEBUG mode."<< endl;
#endif
    // CURRENT PHYSICAL VARIABLES FILE
    const string  simulation_variable_dataset("CSMP_k_analysis_variables.txt");
   
    cout <<"\nmain: What would you like to do? - simulation options: "<< endl
         <<"\t(1) COMPULSORY: read and convert ANSYS model into CSMP native binary format (all other options require this as input)."<< endl

         <<"\t(2) map properties from region-specific point clouds to CSMP model."<< endl

         <<"\t(3) calculate failure criteria, fracture apertures, and permeabilities for a given in situ stress, mechanical rock properties and fluid pressure."<< endl
   
         <<"\t(4) analyse diagonal-only tensor k_equiv, qfqm, and velocity spectra."<< endl

         <<"\t(5) analyse full-tensor k_equiv."<< endl

         <<"\t(6) output 'thickness' and 'permeability' of selected regions to VTU file."<< endl;

//         <<"\t(6) analyse imbibition relative permeability in 2D using experimentally determined matrix properties (Shaho's simulator)."<< endl
   
    int choice(0);
    cin >> choice;
    
    // TODO: let CSMP detect whether the CSMP binary model is already there
    try {
          // tested with 'DFM_modelling_tests/UnitSquareFracs_irregular' as test model, with s1 in NE-SW orientation
          cout <<"\nCSMP++ DFM_Simulator2D: Enter name of simulation two-dimensional model: ";
          string model_name("undefined");
          cin >> model_name;

          // 1. Input mesh from ICEMCFD-ANSYS: reading ANSYS-CSP input file set + regions + configuration file,
          //    assign permeability and porosity and setup the model the data from the configuration file.
          // =============================================================================================
          //    (this step is needed if there is not already a 2D DFM model in CSMP native format)
          if ( choice == 1 ) {
              ANSYS_Model2D ansys_model( model_name.c_str(), simulation_variable_dataset.c_str() );
              Standard_IO_Handler  stdio( model_name.c_str() );
              printModelDimensions( ansys_model );
              ansys_model.RegionsOut();
              ansys_model.BoundariesOut();
              checkNeighborConnectivityOfElementsInFractureRegion( ansys_model );

              // 1.1 reading in model configuration
              if ( ansys_model.Database().IsDefined("run duration") ) {
                  ComputationalSettings  run_settings;
                  InputDataManager<2>().ConfigureFromFile( ansys_model,
                                                           model_name.c_str(),
                                                           false,
                                                           true,   // 2) default prop.values
                                                           true,   // 3) group prop.values
                                                           false,  // 4) essential box-boundary conditions
                                                           true,   // 5) essential flags
                                                           true,   // 6) boundary conditions
                                                           run_settings );

                  // saving time-stepping parameters to model
                  csmp::Index duration_key = ansys_model.Database().StorageKey("run duration");
                  ansys_model.Store(duration_key,makeScalar(PLAIN,run_settings.Duration()));
                  csmp::Index dt_key = ansys_model.Database().StorageKey("default time increment");
                  ansys_model.Store(dt_key,makeScalar(PLAIN,864000.));
                  csmp::Index tout_key = ansys_model.Database().StorageKey("output times");
                  ArrayVariable  output_times(tout_key); // array is automatically initialized to size of array in variables file
                  if ( output_times.Size() < distance( run_settings.OutputTimesBegin(), run_settings.OutputTimesEnd() ) )
                    csmp_error.Note( ERROR, "convertANSYS_to_CSMP_model", "number of 'output times' exceeds maximum defined in '-variables.txt' file.");
                  else {
                      uint32_t counter(0U);
                      for ( auto it=run_settings.OutputTimesBegin(); it!=run_settings.OutputTimesEnd(); it++ )
                        output_times(counter++) = (*it);
                      ansys_model.Store( tout_key, output_times );
                    }
                }
              // without time stepping information
              else
                InputDataManager<2>().ConfigureFromFile( ansys_model,
                                                         model_name.c_str(),
                                                         false,
                                                         true,   // 2) default prop.values
                                                         true,   // 3) group prop.values
                                                         false,  // 4) essential box-boundary conditions
                                                         false ); // 5) essential flags

              // 1.2 Partitioning fracture regions into contiguous fractures
              if ( stdio.YesNo("convertANSYS_to_CSMP_model: Do you want to partition disconnected regions into contiguous subregions") ) {
                   set<string>  original_region_names;
                   for ( auto it=ansys_model.UniqueRegionsBegin(); it!=ansys_model.UniqueRegionsEnd(); ++it )
                     if ( (*it).first != "Model" ) original_region_names.insert( (*it).first.c_str() );
                
                   // partitioning regions without revisiting new partitions that can inserted into region map
                   for ( auto it=original_region_names.begin(); it!=original_region_names.end(); ++it )
                     ansys_model.PartitionRegionIntoContiguousSubRegions( (*it).c_str() );
                }
                
               // 1.3 Checking that all fractures consist of line elements with valid neighbors
               checkNeighborConnectivityOfElementsInFractureRegion( ansys_model );

               // 1.4 writing the CSMP binary
               ansys_model.OutputToBinaryFile( model_name.c_str() );
 
               return 0;
            }

          // 2. Assign permeability and porosity data from a point cloud to the CSMP native model
          // =============================================================================================
          // (use this method to model an inhomgeneous or even heterogeneous matrix)
          else if ( choice == 2 ) {
               cout <<"\nCSMP++ DFM_Simulator2D: mapping property values from points in specific regions to CSMP binary. ";
               pointPropertyMapping( model_name.c_str() );
               return 0;
            }
            
          // 3. Calculate effective stress and evaluate failure criteria for the fractures
          // =============================================================================================
          else if ( choice == 3 ) {
               // descends into calculation option submenu
               fracturePropertyModelling<2>( model_name.c_str() );
               return 0;
            }
            
          // 4. Compute Cartesian equivalent permeability & fracture matrix flux ratio (using bfluxes)
          // =============================================================================================
          else if ( choice == 4 ) {
               const bool with_vtk_output(true);
               DFM_PermeabilityAndFluxRatioAnalyzer<2>  analyzer( model_name.c_str(), with_vtk_output );
               analyzer.DiagonalTensorAnalysis( model_name.c_str(), with_vtk_output );
               return 0;
            }

          // 5. Compute the equivalent permeability & fracture matrix flux ratio (using Durlofsky's method)
          // =============================================================================================
          //    (the output will be a full permeability tensor and directional fracture matrix flux ratio
          else if ( choice == 5 ) {
               fullTensorPermeability_Durlofsky<2>( model_name.c_str() );
               // DFM_PermeabilityAndFluxRatioAnalyzer<3U>  analyzer( model_name.c_str(), verbose );
               // const bool with_vtk_output(true);
               // analyzer.FullTensorAnalysis( model_name.c_str(), with_vtk_output );
               return 0;
            }

          // 6. Output of 'thickness' (=fracture aperture) and 'permeability for selected regions lower-dim fractures
          // ========================================================================================================
          else if ( choice == 6 ) {
               list<string>  result_vars;
               result_vars.push_back("thickness"); // volume modifier for lower dimensional fractures
               result_vars.push_back("permeability");
               Model<2>          model(model_name);
               VTU_Interface<2>  vtk_output( model, "DFM_modeling");
               cout <<"\nCSMP++ DFM_Simulator2D: Enter number and names of region where 'thickness' and 'permeability' shall be output to VTK files: ";
               size_t n_regions;
               cin >> n_regions;
               for ( size_t i=0; i<n_regions; ++i ) {
                    string region;
                    cin >> region;
                    if ( !region.empty() )
                      vtk_output.OutputDataToVTU( (string(model_name) + "-mechanical_characteristics").c_str(),
                                                   result_vars, region.c_str(), static_cast<int>(0) );
                 }
               return 0;
            }
          else
          throw csmp::Exception( FATAL_ERROR, "CSMP++ DFM_Simulator2D", "choice of model not recognized.");
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
    
 } // end main (2D version)






namespace csmp {

/**
        Any region that contains FRAC or 'frac' in its name
 */
bool checkNeighborConnectivityOfElementsInFractureRegion( const Model<2>& model )
 {
    bool everything_OK{true};
    
    cout <<"\ncheckNeighborConnectivityOfElementsInFractureRegion: checking model '"<< model.Name() <<"'"<< endl;
    for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit )
      if ( contains((*rit).first, "FRAC") || contains((*rit).first, "frac") || contains((*rit).first, "FAULT") || contains((*rit).first, "fault") )
        {
           long errors{0};
           cout <<"\n\t"<<"checking region '"<< (*rit).first <<"'";
           for ( const auto& eit : (*rit).second.CellVector() ) {
                 assert( eit != nullptr );
                 // are the frac elements line elements? (must be for this DFM approach to work)
                 if ( !eit->IsLine() ) { cout <<"  elmt "<< eit->Idx() <<": is a "<< parseFiniteElementType( eit->FE_Type() ) <<", "; errors++; }
                 // are the neighbors of each fracture element of the same dimensionality?
                 for ( uint32_t i{0u}; i<eit->Neighbors(); ++i )
                   if ( eit->Neighbor(i) != nullptr && !eit->Neighbor(i)->IsLine() ) {
                        cout <<"  elmt "<< eit->Idx() <<": neighbor "<< i <<": is a "<< parseFiniteElementType( eit->Neighbor(i)->FE_Type() ) <<", ";
                        errors++;
                     }
             }
          if ( errors > 0 ) {
               everything_OK = false;
               cout <<" failed.";
            }
          else cout <<" OK.";
       }
 
    if ( !everything_OK ) return false;
    return true;
 
 } // end


} // end csmp
