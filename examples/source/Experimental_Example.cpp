//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "Experimental_Example.h"
#include "CompareFloats.h"
#include "ANSYS_Model3D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include <iostream>
#include "InputDataManager.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"

// include any header files that you need here...

using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
{
   SetTitle( "Experimental_Example" );
   SetDifficulty( 1 );
   SetCategory( "Software Functionality" );
   AddAuthor( "Irina Sin" );
   AddDescription( "source in: Experimental_Example.cpp" );
   AddDescription( "Empty example for the user to experiment with" );
   AddRequirement( "none" );
   AddRequirement( "no predefined model or variables file" );
}


/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.
*/
void Experimental_Example::Run()
{
  Point<3U> a(1.,2.,3.), b(3.,4.,5.);
  
  Point<3U> c = a + b;

  Point<3U> d = a + b;
  
  a.Out();
  b.Out();
  c.Out();

  return;

  // testing model construction
  const size_t   dim(3U);
  const string   model_name("prism_test");
  ANSYS_Model3D  model( model_name.c_str(), "example25.txt");

  printModelDimensions( model, true );
  
  csmp::Index pf_key = model.Database().StorageKey("fluid pressure");

  InputDataManager<dim>  model_configuration;

  model_configuration.ConfigureFromFile( model, model_name.c_str(),
                                         false, true, true, true, false );

  VTK_Interface<dim>  vtk_output;
  
  Region<dim>& frac_volume(model.Region("FRAC_VOLUMES"));
  cerr <<"\n\n\nrun: indices before modification\n";
  ScalarVariable pf;
  for ( auto it=frac_volume.NodesBegin(); it!=frac_volume.NodesEnd(); ++it ) {
       cerr << (*it)->Idx() <<" ";
       (*it)->Read( pf_key, pf );
       cerr << pf <<" ";
    }
  
  frac_volume.UpdateMemberIndexes();
  cerr <<"\n\n\nrun: indices after modification\n";
  for ( auto it=frac_volume.NodesBegin(); it!=frac_volume.NodesEnd(); ++it ) cerr << (*it)->Idx() <<" ";

  // enlisting the elements to be processed by their indices
  const Region<dim>&  model_region(model.Region("Model"));
  model_region.UpdateMemberIndexes();
  set<size_t> unprocessed_elements;
  for ( auto it=model_region.ElementsBegin(); it!=model_region.ElementsEnd(); it++ )
    unprocessed_elements.insert( (*it)->Idx() );
  
  // computing regions of elements that do not share any nodes
  set<size_t>  discovered_nodes;
  
  // creating the first region named 'region1'
  string prefix("region");
  map<string,set<size_t> > regions_by_elmts;
  set<size_t>              region;
  size_t                   regions(0U);
  
  while ( !unprocessed_elements.empty() ) {
       for ( auto i=unprocessed_elements.begin(); i!=unprocessed_elements.end(); i++ )
         {
            // checking whether any of the nodes of the element is already in the map
            size_t counter(0U);
            const size_t nodes(model_region.E(*i)->Nodes());
            for ( size_t j=0U; j<nodes; j++ )
              if ( discovered_nodes.find( model_region.E(*i)->N(j)->Idx() ) == discovered_nodes.end() ) counter++;
            // when none of the nodes was discovered previously, the element can go into the non-overlapping new set
            if ( counter == nodes ) {
                 region.insert( *i );
                 unprocessed_elements.erase( *i );
                 // marking the nodes as discovered
                 for ( size_t j=0U; j<nodes; j++ )
                   discovered_nodes.insert( model_region.E(*i)->N(j)->Idx() );
              }
         }
       // adding the region to the non-overlapping region map
       if ( !region.empty() ) {
            regions_by_elmts.insert( make_pair( string(prefix) += regions, region ) );
            regions++;
            region.clear();
         }
    }
  

  cout <<"\nExperimental_Example: That's it!\n";
  
} // end Run

} // csmp
