//
//  QuadrilateratorToCSMPbinary_Example.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 21/8/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#include "QuadrilateratorToCSMPbinary_Example.h"

#include "Model.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "Quadrilaterator.h"
#include "InputDataManager.h"
#include "VSet.h"
#include "VTU_Interface.h"
#include "ModelTopology.h"

using namespace std;

namespace csmp {

void QuadrilateratorToCSMPbinary_Example::Specifications()
{
  SetTitle( "Creation of a regular quadrilateral grid and export of model to CSMP native binary file" );
  SetDifficulty( 1 );
  SetCategory( "Software Functionality" );
  AddAuthor( "SKM" );
  AddDescription( "Quadrilaterator class is used to construct the grid" );
  AddDescription( "Some grid-smoothing is illustrated" );
  AddDescription( "source in: QuadrilateratorToCSMPbinary_Example.cpp" );
  AddRequirement( "input file set: permMatrix50x25" );
  AddRequirement( "IMPES-variables.txt" );
  AddRequirement( "permMatrix50x25-configuration.txt" );
}




/**
    Builds a quadrilateral model from cross-section data and saves it as a CSMP binary file for simulation
    with the FECFVM simulator, using the variables defined in the text file "IMPES-variables.txt"
*/
void QuadrilateratorToCSMPbinary_Example::Run()
 {
    // --------------------------------------------
    // 1. Read color-coded pixel image
    // --------------------------------------------
    // permMatrix50x25
    cout <<"\nmain: Enter name of textfile (no extension) with the permeability data defining the flow geometry of the model: ";
    string file_name;
    cin >> file_name;

    Quadrilaterator    quadrilaterator; // simple FE mesher
    VSet<2U>           mesh_container;  // container to store the input mesh

    // 50 x 20 m
    cout <<"\nmain: Enter x(horizontal) and y(vertical) extent of the model: "; 
    double x_extend(10.), y_extend(5.);
    cin >> x_extend >> y_extend;
    const bool from_bitmap=true;
    // also converts integer values to permeability (see utilities/convertColorToPermeability.h)
    quadrilaterator.QuadrilateralsFromRegularGrid( mesh_container, file_name.c_str(), x_extend, y_extend, from_bitmap );


    // -----------------------------------------------------------------
    // 2. Create 2D Model of isoparametric quadrilateral finite elements
    // -----------------------------------------------------------------
    ModelTopology     topology( true /* isoparametric elements */ );
    const set<string> fem_types = {"ISOPARAMETRIC_LINEAR_QUADRILATERAL"};
    vector<size_t>    elements( mesh_container.Elements() );
    iota( begin(elements),  end(elements), 0 );
    // creating a region of all elements 'sediments'
    topology.AddDomain( "sediments",fem_types, elements );
    // constructing the rectangular model              
    const bool create_boundary_objects(false /* since there are no line elements */), box_shaped(true);
    Model<2U>  model( topology, mesh_container, "IMPES-variables.txt", true );
    recreateBoxBoundaryFlags( model );
    model.Name( "permMatrix50x25" ); 
    // give the model dimensions
    printModelDimensions( model, true );


    // ----------------------------------------------
    // 3. Computation of properties from permeability
    //    (Leverett J-analysis etc.)
    // -----------------------------------------------
    const csmp::Index k_key   = model.Database(). StorageKey("permeability");
    const csmp::Index phi_key = model.Database(). StorageKey("porosity");
    const csmp::Index pd_key  = model.Database(). StorageKey("entry pressure");
    const csmp::Index bcp_key = model.Database(). StorageKey("brooks corey parameter");
    Region<2> model_domain(model.Region("Model"));
    
    for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it ) {
          const double permeability = (*it)->Read( k_key );
          // porosity: phi = cubic root of permeability times constant factor + offset
          const double scale_factor(2.5e3);
          double       porosity = scale_factor * cbrt( permeability );
          (*it)->Store( phi_key, makeScalar(PLAIN,porosity) );
          // entry pressure: (decreases with permeability from 10 MPa to 1kPa at 1D)
          double       pd = 1.0e-9 * (1/permeability);
          (*it)->Store( pd_key, makeScalar(PLAIN,pd) );
          // Brooks-Corey parameter: (increases with permeability from 3 to 6)
          double       bcp = -log(permeability)/5. - 3.;
          (*it)->Store( bcp_key, makeScalar(PLAIN,bcp) );
      }


    // ----------------------------------------------
    // 4. Assignment of (fluid) properties from file
    // -----------------------------------------------
    // default property values
    InputDataManager<2U>().ConfigureFromFile( model,
                                              model.Name(),
                                              false,   // 1) IGNORED
                                              true,    // 2) default prop.values
                                              true,    // 3) region prop.values
                                              true,    // 4) essential box-boundary conditions
                                              true,    // 5) essential flags (inlet)
                                              false ); // 6) IGNORED - generic boundary conditions

 
    // -------------------------------------------------------------------------------
    // 5. Creating and parameterising 1-element wide high-perm inflow zone on the left
    // -------------------------------------------------------------------------------
    const bool with_inflow_region(true);
    if ( with_inflow_region ) {
          createInflowRegion( model );
          Region<2U>&  well(model.Region("inflow"));
          well.InputPropertyValue("rock type", makeScalar(ANY,0.) );
          well.InputPropertyValue("porosity", makeScalar(ANY,1.) );
          well.InputPropertyValue("permeability", makeScalar(ANY,1.0e-10) );
          well.InputPropertyValue("entry pressure", makeScalar(ANY,1.) );
          well.InputPropertyValue("residual saturation wetting phase", makeScalar(ANY,0.) );
          well.InputPropertyValue("residual saturation non-wetting phase", makeScalar(ANY,0.) );
          well.InputPropertyValue("brooks corey parameter", makeScalar(ANY,0.) );
      }
    // NB: boundary conditions were supplied to left boundary via config file

    // --------------------------------------------
    // 6. Output to VTU and CSMP binary
    // --------------------------------------------
    list<string> output_vars = { "rock type", "permeability", "porosity", "entry pressure", "brooks corey parameter" };
    VTU_Interface<2U>  vtu_output( model );
   
    vtu_output.OutputDataToVTU( string( string( model.Name() ) + "-static_model" ).c_str(), output_vars, "Model", 0 );

    // writing the CSMP binary
    model.OutputToBinaryFile( model.Name() );
 
 } // end run







// AUXILIARY FUNCTIONS

/**
    Robust detection of NO_DATA vales
*/
inline bool is_NO_DATA_Value( double value )
 {
    // checking the value range, including no-data values
    assert( !isnan(value) );
    if ( static_cast<long>(value) == -9999 || static_cast<long>(value) == -99999  ) return true;
    return false;
 }




 
 

/**
    Turns elements located on the left model boundary into a region called 'inflow'.
*/
void createInflowRegion( Model<2U>& model )
  {
     Region<2U>&  domain(model.Region("Model"));
     domain.UpdateMemberIndexes();
    
     vector<size_t> element_ids;
     element_ids.reserve(800); // there are 800 elements in the vertical
    
     for ( auto eit=domain.CellsBegin(); eit!=domain.CellsEnd(); ++eit ) {
          // for quadrilateral elements, if any of their nodes are on the left boundary the element is as well
          assert( (*eit)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL );
          bool at_left_boundary(false);
          for ( auto i{0}; i<(*eit)->Nodes(); ++i )
            if ( isLEFT( (*eit)->N(i)->AtBoundary() ) ) {
                 at_left_boundary = true;
                 break;
              }
          if ( at_left_boundary ) element_ids.push_back( (*eit)->Idx() );
       }
    
     // creating the region
     if ( !model.FormRegionFrom( "inflow", element_ids ) )
       cerr <<"\ncreateInflowRegion: region 'inflow' could not be formed.\n";
     else
       cout <<"\ncreateInflowRegion: region 'inflow' formed successfully from "<< element_ids.size() <<" elements.\n";

  } // end createInflowRegion



} // end csmp
