//
//  SKUA_FiniteElementMeshInterface_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 2/12/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "SKUA_FiniteElementMeshInterface_Test.h"
#include "SKUA_Model.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

void SKUA_FiniteElementMeshInterface_Test::run()
 {
    // 0. creating box-shaped test model with nodes on the side surfaces
    SKUA_Model  model( "V6_SOLID2",
                       "CSMP-variables.txt",
                       false,   ///< true = binary, false = ascii */
                       true,    ///< true = reduce regions according to regions file, false = does not redure regions
                       false,   ///< true = creates boundaries around model, false = does not create boundaries 
                       false,   ///< true = creates splitboundaries around model, false = does not create splitboundaries 
                       true );  ///< isoparametric is the only option if model does not only consist of simplex elements

    printModelDimensions( model, true );
    model.RegionsOut();

    // Testing
    _test( TestNeighborConnectivity( model ) );


    // checking model boundary flagging
    boxFlagsToVariable( model, "node variable", "element variable" );
    
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "porosity", "porosity", 1 );
    vtk_output.OutputDataToVTK( model, "node-flag", "node variable", 1 );
    vtk_output.OutputDataToVTK( model, "element-flag", "element variable", 1 );
    
    // 1. basic testing
    /// consecutive numbering of nodes and elements
    /// mesh consistency (Jacobians, element volume range, non-manifold vertices, triangle boxes
    /// pfverts, disambiguated neighbors of lower-dimensional elements
    /// consistent facing directions of surface elements
    /// consecutive numbering of series of line elements
    /// material IDs
    /// test assignment of box boundary flags (sides, edges, corners)
    /// property assignment
    /// use of regions file (is this still needed?)
    // 2. advanced testing
    /// generation of boundaries  
    /// generation of split boundaries                      

 } // end run





/**
   Checks whethe the assignment of element neighbors from SKUA is identical with that done by CSMP; the side boundary flags are not considered.
   
      TODO: also check neighbor connectivity of Face and InterFace objects
*/
bool SKUA_FiniteElementMeshInterface_Test::TestNeighborConnectivity( Model<3U>& model, bool verbose )
 {
    const Region<3U>& domain(model.Region("Model"));
    // recording connectivity from SKUA in an element neighbor vector
    vector<vector<Element<3U>*> > pfverts;
    pfverts.reserve( domain.Elements() );
    for ( vector<Element<3U>*>::const_iterator
          it=domain.ElementsBegin(); it!=domain.ElementsEnd(); ++it ) {
         vector<Element<3U>*> nbors( (*it)->Neighbors(), nullptr );
         for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr )
             nbors[i] = (*it)->Neighbor(i);
         pfverts.emplace_back( nbors );    
      }
    // recreating the connectivity in CSMP
    establishNeighborConnectivity( model.Region("Model").CellVector() );
    
    // comparing SKUA with CSMP connectivity
    size_t failed_comparisons(0U);
    vector<vector<Element<3U>*> >::const_iterator pfit(pfverts.begin());
    for ( vector<Element<3U>*>::const_iterator
          it=domain.ElementsBegin(); it!=domain.ElementsEnd(); ++it, ++pfit ) {
         for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != (*pfit)[i] ) {
                if ( verbose ) {
                     cerr <<"\nelement neighbor "<< (*it)->Idx() <<":"<< i <<": ";
                     cerr << (*it)->Neighbor(i)->Idx() <<" vs. "<< (*pfit)[i]->Idx() <<" ";
                  }
                failed_comparisons++;
             }
      }
    if ( verbose && failed_comparisons > 0 ) cerr << endl;
      
    if ( failed_comparisons > 0U ) return false;  
    return true;
    
 } // end TestNeighborConnectivity
 

} // end csmp
