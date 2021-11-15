//
//  SKUA_FiniteElementMeshInterface_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 2/12/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "SKUA_FiniteElementMeshInterface_Test.h"
#include "MeshManagementUtilities.h"
#include "CSMP_highLevelUtilities.h"
#include "SKUA_Model.h"
#include "Region.h"
#include "Element.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

void SKUA_FiniteElementMeshInterface_Test::run()
 {
    // 0. creating box-shaped test model with nodes on the side surfaces
    SKUA_Model  model( "V6_SOLID2_withwells_BINARY",   // V6_SOLID2_with_wells, "V6_SOLID2", ,"V1_SOLID2_ASCII", "V6_SOLID2_withwells_ASCII", "EW01_100ftBAND_10ft"
                       "CSMP-variables.txt",
                       true,   ///< true = binary, false = ascii */
                       true,    ///< true = reduce regions according to regions file, false = does not redure regions
                       false,   ///< true = creates boundaries around model, false = does not create boundaries 
                       false,   ///< true = creates splitboundaries around model, false = does not create splitboundaries 
                       true );  ///< isoparametric is the only option if model does not only consist of simplex elements

    printModelDimensions( model, true );
    model.RegionsOut();

/*
    // Testing
    _test( model.Mesh().Elements() == 21116/4 );
    _test( model.Mesh().Nodes() == 1262 );
*/   
//    _test( TestNeighborConnectivity( model ) );
    const Region<3U>&  unit_A(model.Region("UNIT_A"));
    _test( unit_A.Elements() == 1775 );
/*
    const Region<3U>&  unit_B(model.Region("UNIT_B"));
    _test( unit_B.Elements() == 1743 );

    const Region<3U>&  unit_C(model.Region("UNIT_C"));
    _test( unit_C.Elements() == 1761 );
    
    // are the property values assigned in the correct sequence?
    _test( TestNodePropertyAssignment( model ) );
*/
    // checking model boundary flagging
    boxFlagsToVariable( model, "nodal variable", "element variable" );
    
    VTK_Interface<3U>  vtk_output;
//    vtk_output.OutputDataToVTK( model, "porosity", "porosity", 1 );
    vtk_output.OutputDataToVTK( model, "node-flag", "nodal variable", 1 );
    vtk_output.OutputDataToVTK( model, "element-flag", "element variable", 1 );

    // creating 'permeability' variable if it does not already exist
    if ( !model.Database().IsDefined("permeability") ) {
         model.CreateProperty( "permeability", "m2", SCALAR, ELEMENT );
      }
    OutputRegionsToVTK( model, "permeability" );
    
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



// TODO: such a method must already exist somewhere else
void SKUA_FiniteElementMeshInterface_Test::OutputRegionsToVTK( const Model<3U>& model, const char* var_name ) const
 {
     std::cout <<"\n\nSKUA_FiniteElementMeshInterface_Test::OutputRegionsToVTK:\n";
     std::cout <<"\n\tUnique regions of model:\n";
     for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit ) {
          std::cout <<"\t\t"<< (*rit).first;
          std::cout <<" "<< (*rit).second.Elements() <<" elements,";
          std::pair<int32_t, int32_t> rdim = (*rit).second.ElementSpatialDimensions();
          if ( rdim.second == 3 )
            std::cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 2 )
            std::cout <<" surface area (m2): "<< (*rit).second.Volume() <<", perimeter length (m): "<< (*rit).second.SurfaceArea();
          std::cout <<", range of spatial dimensions: "<< rdim.first <<", highest spatial dimension "<< rdim.second << std::endl;
       }

     // writing all regions with model name - region name and variable name
     VTK_Interface<3U>().OutputRegionByRegionToVTK( model, "SKUA_CSMP_test_", var_name, 0 );

 } // end OutputRegionsToVTK



/**
   Checks whethe the assignment of element neighbors from SKUA is identical with that done by CSMP; the side boundary flags are not considered.
   
      TODO: also check neighbor connectivity of Face and InterFace objects
*/
bool SKUA_FiniteElementMeshInterface_Test::TestNeighborConnectivity( Model<3U>& model, bool verbose )
 {
    const csmp::Index eidx = model.Database().StorageKey("element number");
    const csmp::Index nidx = model.Database().StorageKey("node number");
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
      
   if ( verbose ) PrintOriginalNeighborIDs( model );
      
    // recreating the connectivity in CSMP
    establishNeighborConnectivity( model.Region("Model").CellVector() );
    if ( verbose ) PrintOriginalNeighborIDs( model );
    
    // comparing SKUA with CSMP connectivity
    size_t failed_comparisons(0U);
    vector<vector<Element<3U>*> >::const_iterator pfit(pfverts.begin());
    for ( vector<Element<3U>*>::const_iterator
          it=domain.ElementsBegin(); it!=domain.ElementsEnd(); ++it, ++pfit ) {
         for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != (*pfit)[i] ) {
                if ( verbose ) {
                     //(*it)->Out();
                     cerr <<"\nelement neighbor "<< (*it)->Read(eidx) <<":"<< i <<": ";
                     cerr << (*it)->Neighbor(i)->Idx() <<" vs. "<< (*pfit)[i]->Read(eidx) <<" ";
                  }
                failed_comparisons++;
             }
      }
    if ( verbose && failed_comparisons > 0 ) cerr << endl;
      
    if ( failed_comparisons > 0U ) return false;  
    return true;
    
 } // end TestNeighborConnectivity
 
 
 
bool SKUA_FiniteElementMeshInterface_Test::TestNodePropertyAssignment( const Model<3>& model )
 {
    // if the element number is known so that the material IDs can be assigned correctly
    if ( model.Database().IsDefined("node number") ) {
         const csmp::Index nn_key = model.Database().StorageKey("node number");
         const csmp::Index nu_key = model.Database().StorageKey("number");
         // creating a mapping between current elements in region 'Model' and the VSet from the element number
         const Region<3U>& domain = model.Region("Model");
         vector<Node<3U>*> ordered_nodes(domain.Nodes());

         for ( vector<Node<3U>*>::const_iterator nit=domain.NodesBegin(); nit!=domain.NodesEnd(); ++nit ) {
              if ( fabs((*nit)->Read(nn_key) - (*nit)->Read(nu_key)) > numeric_limits<double>::epsilon() ) { 
                   cerr <<"\n\tnode number vs. number: "<< (*nit)->Read(nn_key) <<" vs. "<< (*nit)->Read(nu_key);
                   _equal( (*nit)->Read(nn_key), (*nit)->Read(nu_key), numeric_limits<double>::epsilon() );
                }
              ordered_nodes[ static_cast<size_t>((*nit)->Read(nn_key)) ] = (*nit);
           }
           
         bool problem_found(false);  
         for ( size_t i=0U; i<ordered_nodes.size(); ++i ) {
              if ( ordered_nodes[i] == nullptr ) {
                  cerr <<"\nordered nodes map contains nullptr for node "<< i; 
                  problem_found = true;
                }
              else if ( i != static_cast<size_t>(ordered_nodes[i]->Read(nu_key)) ) {
                  cerr <<"\nmismatch: node "<< i <<": vs 'number' "<< ordered_nodes[i]->Read(nu_key);
                  problem_found = true;
                }
            }
         if ( problem_found == true ) return false;  
      }       
    return true;

 } // TestNodePropertyAssignment



void SKUA_FiniteElementMeshInterface_Test::PrintOriginalNeighborIDs( const Model<3U>& model ) const
 {
    // if the element number is known so that the material IDs can be assigned correctly
    if ( model.Database().IsDefined("element number") ) {
         const csmp::Index eid_key = model.Database().StorageKey("element number");
         // creating a mapping between current elements in region 'Model' and the VSet from the element number
         const Region<3U>& domain = model.Region("Model");
         vector<Element<3U>*> ordered_elmts(domain.Elements());
         for ( vector<Element<3U>*>::const_iterator it=domain.ElementsBegin(); it!=domain.ElementsEnd(); ++it )
           ordered_elmts[ static_cast<size_t>((*it)->Read(eid_key)) ] = (*it);
           
         // reading the VSet material record
         for ( vector<Element<3U>*>::const_iterator it=ordered_elmts.begin(); it!=ordered_elmts.end(); ++it ) {
              cerr <<"\nelmt "<< (*it)->Read( eid_key ) <<": ";
              // printing the neighbors of this element
              for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
                if ( (*it)->Neighbor(i) != nullptr )
                  cerr << (*it)->Neighbor(i)->Read( eid_key ) <<" ";
                else 
                  cerr << parseBoundary( (*it)->AtBoundary(i) ) <<" ";
           }
         cerr << endl;
      }
    
  } // end PrintOriginalNeighborIDs


} // end csmp
