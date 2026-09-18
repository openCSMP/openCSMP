// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  convertANSYS_to_CSMP_model.cpp
//  CSMP_DFM_Upscaling
//
//

#include "convertANSYS_to_CSMP_model.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "RegionInterface.h"
#include "ModelTime.h"
#include "InputDataManager.h"
#include "Standard_IO_Handler.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "ComputationalSettings.h"
#include "GlobalVerbose.h"
#include "Node.h"
#include "VTK_Interface.h"

/** read_and_configure_ANSYS_model

    1. Reads an ANSYS model that was output by ANSYS in CSMP format,
      (an enormous number of checks are performed here to make sure 
       that the input mesh is fit for purpose, also, the user has the 
       option to rotate coordinate axis to make the model match CSMP,
       which assumes that Y axis points upward, X from W to E, and 
       Z towards us, from N to S),
    
    2. using the information provided by the "-regions.txt' file,
       selects the mesh regions that shall be kept for the simulation
      (ANSYS outputs every curve, surface and volume that was part of the original CAD model,
       but many of these elements are not needed for the CSMP-based simulations)
       
    3. assigns material properties to regions, fluid properties to nodes, and 
       essential and initial conditions from "-configuration.txt" file
      (these can later be overwritten)
      
    4. gets a preliminary set-up for the simulation run, regarding duration,
       output time steps and time-stepping strategy. This is stored as ComputationalSettings object.
       
    @attention If the variable 'run duration' is defined in the var_file, this info is stored in the model
 
    @test not tested yet   
*/

//#define convertANSYS_to_CSMP_model_TESTING

using namespace std;

namespace csmp {

void restrictLowerDimensionalElementConnectivityByRegion( Region<2U>& );
void restrictLowerDimensionalElementConnectivityByRegion( Region<3U>& );

/**
    Converts ANSYS -> CSP input deck into the CSMP native file format and writes this to disk.
    There are 2 options.
    
    1) the model can be rotated to match the coordinate system that is standard to the oil industry
    where y points downward.
    
    2) Each region of the model can be partitioned into contiguous subregions.
    If this results into multiple subregions per region, these are distinguished by the extensions
    1..n, as appended to the region name.

    @attention neighbor ambiguity of lower-dimensional elements (which may have multiple neighbors
    per face, but only one neighbor pointer per face) is removed by considering only elements of 
    the same region as potential neighbors. This requires the user to avoid such manifolds
    by appropriate partitioning of the regions in the ANSYS model.
*/
bool convertANSYS3D_to_CSMP_model( const char* ansys_model_name, const char* var_file_name )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      Standard_IO_Handler  stdio( ansys_model_name );
      int choice(-1);
    
      if ( stdio.YesNo("convertANSYS_to_CSMP_model: Do you want to rotate the coordinate system") ) {
           cout <<"\n\tRotate from ANSYS (1) or SKUA (2) to CSMP coordinate system ? ";
           cin >> choice;
        }

      ANSYS_Model3D  ansys_model( ansys_model_name, var_file_name );
      Region<3U>&    mref(ansys_model.Region("Model"));
     
      printModelDimensions( ansys_model );
      cout <<"\nread_and_configure_ANSYS_model: model volume:       "<< mref.Volume() << endl;
      cout <<"\nread_and_configure_ANSYS_model: model surface area: "<< mref.SurfaceArea() << endl;

      // 1. reading only selected blocks
      // ------------------------------------------------------
      if ( ansys_model.Database().IsDefined("run duration") ) {
          ComputationalSettings  run_settings;
          InputDataManager<3U>().ConfigureFromFile( ansys_model, 
                                               ansys_model_name,
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
        InputDataManager<3U>().ConfigureFromFile( ansys_model,
                                                  ansys_model_name,
                                                  false,
                                                  true,   // 2) default prop.values
                                                  true,   // 3) group prop.values
                                                  false,  // 4) essential box-boundary conditions
                                                  false ); // 5) essential flags
    

      // 2. rotating model if so required by user
      // ------------------------------------------------------
      if ( choice == 1 ) { // the y-coordinate needs to become z and z=-y
          cout <<"\n\tTransforming ANSYS coordinate system to CSMP: turning y-coordinate into z and z into -y..."<< endl;
          for ( auto nit=mref.NodesBegin(); nit!=mref.NodesEnd(); ++nit ) {
               double old_y = (*nit)->y();
               (*nit)->y( (*nit)->z() );
               (*nit)->z( -old_y );
            }
        }
      else if ( choice == 2 ) {
          cout <<"\n\tTransforming SKUA coordinate system to CSMP: Nothing is done because this operation would violate righthand rule..."<< endl;
        }
      else cout <<"\nconvertANSYS_to_CSMP_model: rotation choice not recognized, nothing was done."<< endl;
 
          
      // 3. partitioning regions into contiguous subregions
      // --------------------------------------------------
      if ( stdio.YesNo("convertANSYS_to_CSMP_model: Do you want to partition disconnected regions into contiguous subregions") ) {
           set<string>  original_region_names;
           for ( auto it=ansys_model.UniqueRegionsBegin(); it!=ansys_model.UniqueRegionsEnd(); ++it )
             if ( (*it).first != "Model" ) original_region_names.insert( (*it).first.c_str() );
        
           // partitioning regions without revisiting new partitions that can inserted into region map
           for ( auto it=original_region_names.begin(); it!=original_region_names.end(); ++it )
             ansys_model.PartitionRegionIntoContiguousSubRegions( (*it).c_str() );
        }

#ifdef convertANSYS_to_CSMP_model_TESTING
    // region-by region check of partitioning
    // --------------------------------------
    ansys_model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
    // marking the perimeter nodes=10 and the interion ones 0
    const csmp::Index f_key(ansys_model.Database().StorageKey("fluid pressure"));
    for ( auto it=ansys_model.UniqueRegionsBegin(); it!=ansys_model.UniqueRegionsEnd(); it++ )
      for ( auto nit=(*it).second.PerimeterNodesBegin(); nit!=(*it).second.NodesEnd(); ++nit )
        (*nit)->Store( f_key, makeScalar(FIELD_DATA,10.) );
   
    VTK_Interface<3U>().OutputRegionByRegionToVTK( ansys_model, ansys_model_name, "fluid pressure", 0L );
#endif
    
      // 4. writing the CSMP binary
      // ------------------------------------------
      ansys_model.OutputToBinaryFile( ansys_model_name );
    
      return true;
    
  } // end convertANSYS_ToCSMP_Model



/**
       FOR TWO-DIMENSIONAL ANSYS MODELS
 */
bool convertANSYS2D_to_CSMP_model( const char* ansys_model_name, const char* var_file_name )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      Standard_IO_Handler  stdio( ansys_model_name );
      int choice(-1);
    
      if ( stdio.YesNo("convertANSYS_to_CSMP_model: Do you want to rotate the coordinate system") ) {
           cout <<"\n\tRotate from ANSYS (1) or SKUA (2) to CSMP coordinate system ? ";
           cin >> choice;
        }

      ANSYS_Model2D  ansys_model( ansys_model_name, var_file_name );
      Region<2>&     mref(ansys_model.Region("Model"));
     
      printModelDimensions( ansys_model );
      cout <<"\nread_and_configure_ANSYS_model: model volume:       "<< mref.Volume() << endl;
      cout <<"\nread_and_configure_ANSYS_model: model surface area: "<< mref.SurfaceArea() << endl;

      // 1. reading only selected blocks
      // ------------------------------------------------------
      if ( ansys_model.Database().IsDefined("run duration") ) {
          ComputationalSettings  run_settings;
          InputDataManager<2U>().ConfigureFromFile( ansys_model,
                                               ansys_model_name,
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
        InputDataManager<2U>().ConfigureFromFile( ansys_model,
                                                  ansys_model_name,
                                                  false,
                                                  true,   // 2) default prop.values
                                                  true,   // 3) group prop.values
                                                  false,  // 4) essential box-boundary conditions
                                                  false ); // 5) essential flags
    

      // 2. rotating model if so required by user
      // ------------------------------------------------------
      if ( choice == 1 ) { // the y-coordinate needs to become z and z=-y
          cout <<"\n\tTransforming ANSYS coordinate system to CSMP: turning y-coordinate into z and z into -y..."<< endl;
          for ( auto nit=mref.NodesBegin(); nit!=mref.NodesEnd(); ++nit ) {
               double old_y = (*nit)->y();
               (*nit)->y( (*nit)->z() );
               (*nit)->z( -old_y );
            }
        }
      else if ( choice == 2 ) {
          cout <<"\n\tTransforming SKUA coordinate system to CSMP: Nothing is done because this operation would violate righthand rule..."<< endl;
        }
      else cout <<"\nconvertANSYS_to_CSMP_model: rotation choice not recognized, nothing was done."<< endl;
 
          
      // 3. partitioning regions into contiguous subregions
      // --------------------------------------------------
      if ( stdio.YesNo("convertANSYS_to_CSMP_model: Do you want to partition disconnected regions into contiguous subregions") ) {
           set<string>  original_region_names;
           for ( auto it=ansys_model.UniqueRegionsBegin(); it!=ansys_model.UniqueRegionsEnd(); ++it )
             if ( (*it).first != "Model" ) original_region_names.insert( (*it).first.c_str() );
        
           // partitioning regions without revisiting new partitions that can inserted into region map
           for ( auto it=original_region_names.begin(); it!=original_region_names.end(); ++it )
             ansys_model.PartitionRegionIntoContiguousSubRegions( (*it).c_str() );
        }

#ifdef convertANSYS_to_CSMP_model_TESTING
    // region-by region check of partitioning
    // --------------------------------------
    ansys_model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
    // marking the perimeter nodes=10 and the interion ones 0
    const csmp::Index f_key(ansys_model.Database().StorageKey("fluid pressure"));
    for ( auto it=ansys_model.UniqueRegionsBegin(); it!=ansys_model.UniqueRegionsEnd(); it++ )
      for ( auto nit=(*it).second.PerimeterNodesBegin(); nit!=(*it).second.NodesEnd(); ++nit )
        (*nit)->Store( f_key, makeScalar(FIELD_DATA,10.) );
   
    VTK_Interface<2U>().OutputRegionByRegionToVTK( ansys_model, ansys_model_name, "fluid pressure", 0L );
#endif
    
      // 4. writing the CSMP binary
      // ------------------------------------------
      ansys_model.OutputToBinaryFile( ansys_model_name );
    
      return true;
    
  } // end convertANSYS3D_ToCSMP_Model








/**
    As soon as a model contains lower-dimensional elements,
    their connections to equidimensional neighbor elements are no longer uniquely defined:
    A line element, for instance, may connect to multiple other line elements.
    Also, lower-dimensional regions, the elements of which share faces 
    will be traversed by flood-fill algorithms trying to break regions into
    contiguous subregions.
    
    This method tries to remove potential manifolds and restrict neighbor
    connectivity strictly to the region provided as function argument.
    
    @author SKM 30/8/14
    
    @test ?
*/
void restrictLowerDimensionalElementConnectivityByRegion( Region<3U>& regionRef )
  {
     // looping over the region boundary, nulling neighbor element pointers
     // of the boundary faces
     for ( size_t i=regionRef.InteriorCells(); i<regionRef.Cells(); ++i )
       for ( uint32_t j=0U; j<regionRef.PerimeterFaces(i); ++j ) {
            uint32_t n_neighbor = regionRef.PerimeterFace( i, j );
            regionRef.E(i)->Unassign( regionRef.E(i)->Neighbor( n_neighbor ) );
         }

  } // end restrictLowerDimensionalElementConnectivityByRegion

void restrictLowerDimensionalElementConnectivityByRegion( Region<2U>& regionRef )
  {
     // looping over the region boundary, nulling neighbor element pointers
     // of the boundary faces
     for ( size_t i=regionRef.InteriorCells(); i<regionRef.Cells(); ++i )
       for ( uint32_t j=0U; j<regionRef.PerimeterFaces(i); ++j ) {
            uint32_t n_neighbor = regionRef.PerimeterFace( i, j );
            regionRef.E(i)->Unassign( regionRef.E(i)->Neighbor( n_neighbor ) );
         }

  } // end restrictLowerDimensionalElementConnectivityByRegion



} // end csmp

