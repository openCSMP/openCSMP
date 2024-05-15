//
//  SplitBoundary_Example.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 12/11/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#include "SplitBoundary_Example.h"

#include "ANSYS_Model2D.h"
#include "SplitBoundary.h"
#include "VTU_Interface.h"

using namespace std;

namespace csmp {

void SplitBoundary_Example::Specifications()
  {
     SetTitle( "SplitBoundary_Example" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "SKM" );
     AddDescription( "source in: SplitBoundary_Example.cpp" );
     AddDescription( "TBA" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }

// OTHER RELATIONSHIPS

// helper
void printNeigboursOfPerimeterElements( const PropertyDatabase<2>& pbase, const Region<2>& domain );

/**
     Put CSMP code that you would like to test here and run it as part of the
     example suite.
     
     Systems Modelling and Design model that calculates the effective stress in a dam,
     using the gravitational loading and the plane stress assumption.
     
     input models:
     - slope_model1
*/
void SplitBoundary_Example::Run()
{
  // 0. converting ANSYS input file set into a CSMP model and save it to binary file
  // -------------------------------------------------------------------------------
  const bool binary_file(true);             /* true = binary, false = ascii */
  const bool use_regions_file(true);        /* true = reduce regions according to regions file, false = does not redure regions */
  // Ansys model
  const string model_name("BoxHalfs2D");
  ANSYS_Model2D ansys_model("BoxHalfs2D", "BoxHalfs2D", "THMC_shear_zone-variables.txt",
                             binary_file, use_regions_file );
                             
  // testing region insertion here
  Region<2> region1_before(ansys_model.Region("MATRIX_RIGHT"));
//  printNeigboursOfPerimeterElements( ansys_model.Database(), region1_before );
//  ansys_model.InsertSplitBoundary( "MATRIX_RIGHT", "MATRIX_LEFT" );

  // saving model into CSMP native file format
  ansys_model.OutputToBinaryFile( model_name.c_str() );

  // 1. starting the simulation with the creation of a SplitBoundary
  // -------------------------------------------------------------------------------
  // read model model from file and get started with SplitBoundary code
  Model<2> model( model_name );
  model.RegionsOut();
  model.BoundariesOut();
  // checking the regions of the model (visualising their perimeter)
  VTU_Interface<2>  vtk_out( model );
  model.InputPropertyValue( "test variable", makeScalar(ANY,0.) );
  model.Region("MATRIX_LEFT").InputPropertyValue( "test variable", makeScalar(ANY,1.), PERIMETER );
  model.Region("MATRIX_RIGHT").InputPropertyValue( "test variable", makeScalar(ANY,1.), PERIMETER );
  vtk_out.OutputDataToVTU( "region-flag", "test variable", "MATRIX_LEFT", 1 );
  vtk_out.OutputDataToVTU( "region-flag", "test variable", "MATRIX_RIGHT", 2 );

  // create SplitBoundary between the model regions
  Region<2> region1_after(model.Region("MATRIX_RIGHT"));
//  printNeigboursOfPerimeterElements( ansys_model.Database(), region1_before );
  
  model.CreateSplitBoundaryBetween( "MATRIX_RIGHT", "MATRIX_LEFT" );
  // putting a lower dimensional region inside of all split boundaries
  const int32_t material_id_for_new_elements(8);
  set<string>  newly_created_regions = model.InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id_for_new_elements );  
  assert( !newly_created_regions.empty() );
  model.RegionsOut();
  // getting a reference to the split boundary
  SplitBoundary<2>  interface( model.SplitBoundary("SPLITBOUNDARY_MATRIX_RIGHT_MATRIX_LEFT") );
  // getting a reference to the newly created region
  Region<2>  detached_surface( model.Region( (*newly_created_regions.begin()) ) );

  cout <<"\nExperimental_Example: That's it!\n";
  
} // end Run


void printNeigboursOfPerimeterElements( const PropertyDatabase<2>& pbase, const Region<2>& domain )
 {
     const csmp::Index key(pbase.StorageKey("element number"));
     map<long,const Element<2>* const>  ordered_elmts;
     
     for ( auto it=domain.PerimeterCellsBegin(); it!=domain.CellsEnd(); ++it ) {
           ordered_elmts.insert( make_pair( (*it)->Read(key), (*it) ) );
       }
       
     // printint the neighbours
     cout <<"\nprintNeigboursOfPerimeterElements: of region '"<< domain.Name() <<"'\t";
     for ( auto it=ordered_elmts.begin(); it!=ordered_elmts.end(); ++it ) {
          cout <<"\n"<< (*it).first <<": ";
          for ( auto i{0}; i<(*it).second->Faces(); ++i )
            if ( (*it).second->Neighbor(i) == nullptr )
              cout <<"NONE ";
            else
              cout << (*it).second->Neighbor(i)->Read(key) <<" ";
       }
       
 } // printNeigboursOfPerimeterElements

} // csmp
