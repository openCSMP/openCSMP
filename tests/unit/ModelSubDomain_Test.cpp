//
//  ModelSubDomain_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include <ranges>    // C++20 ranges
#include "ModelSubDomain_Test.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyConstraints.h"
#include "SplitBoundary.h"
#include "ANSYS_Model3D.h"
#include "CSMP_highLevelUtilities.h"
#include "meshManagementUtilities.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "vsetMakers.h"
#include "Element.h"
#include "Box.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {


/**
     Observations
     - AllElements in single domain model is non-unique, but should be unique
      
     add tests for the other important parts of the functionality.
*/
void ModelSubDomain_Test::run()
  {
     // new comprehensive tests
//     Test_Basics();
//     Test_SubDomainConstructionMethods(); // TODO: Identify perimeter for non-unique model domain with SplitBoundaries fails!
//     Test_SubDomainDiagnostics();
//     Test_GeometricOperations();
//     Test_PropertyManipulations();
       Test_PropertyTransfer();
       Test_NonMemberFunctions();
  
  
     // Test 0: methods of subdomain in live model
     // ------------------------------------------
     _test( Test_EstablishNeighborConnectivity() );
     
     Test_RebuildCellAndNodeVectors();
  
     bool test_binary_file_recovery1(true),
          verbose(true);
    
     // Test 1: subdomain storage to file and recreation in a new model
     // ---------------------------------------------------------------
     if ( test_binary_file_recovery1 )
       {
         VSet<3U>   vset;
         const bool skewed(false);
         create_Prism_Hexa_VSet( vset, skewed );
         vset.RemoveData("element number");
         Model<3U>   model1( vset );
        
         model1.CreateProperty( "box flag", "bFn", "none", SCALAR, NODE );
         model1.CreateProperty( "box flag element", "bFe", "none", SCALAR, ELEMENT );
         boxFlagsToVariable( model1, "box flag", "box flag element" );
         if ( verbose_ ) {
              VTK_Interface<3U>  vtk_output;
              vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
              vtk_output.OutputDataToVTK( model1, "element_flag", "box flag element", 0 );
           }
           
         cerr <<"\nModelSubDomain_Test::run: original model.";
         if ( verbose_ ) model1.Out();
        
         model1.OutputToBinaryFile("ModelSubDomain_Test");
        
         // testing: model2.InputFromBinaryFile("ModelSubDomain_Test");
         Model<3U>  model2( string("ModelSubDomain_Test") );
         if ( verbose_ ) {
              VTK_Interface<3U>  vtk_output;
              vtk_output.OutputDataToVTK( model2, "node_flag", "box flag", 1 );
              vtk_output.OutputDataToVTK( model2, "element_flag", "box flag element", 1 );
           }
         cerr <<"\nModelSubDomain_Test::run: model reconstructed from disk.";
         if ( verbose_ ) model2.Out();
        
        // TESTING
         _test( CompareModelSubdomains( model1.Region("Model"), model2.Region("Model"), verbose ) );
       }
        
  } // end run




void ModelSubDomain_Test::Test_RebuildCellAndNodeVectors()
 {
    // using a simple model from vsetMakers
    VSet<2> vset;
    create_TrianglePatch_VSet( vset );
    Model<2>  model( vset, "MeshManager_Test-variables.txt" );
    Region<2> model_domain = model.Region("Model");
    auto      n_interior_elmts  = distance( model_domain.CellsBegin(), model_domain.PerimeterCellsBegin() );
    auto      n_perimeter_elmts = distance( model_domain.PerimeterCellsBegin(), model_domain.CellsEnd() );
    _test( n_interior_elmts == static_cast<long>(model_domain.InteriorCells()) );
    _test( n_interior_elmts == static_cast<long>(model_domain.Cells()) - n_perimeter_elmts );
    model_domain.UpdateMemberIndexes();
    // remembering the perimeter faces
    vector<uint32_t> original_perimeter_faces;
    original_perimeter_faces.reserve( model_domain.PerimeterCells() );
    for ( size_t eid{model_domain.InteriorCells()}; eid<model_domain.Cells(); ++eid )
      for ( uint32_t i{0u}; i<model_domain.PerimeterFaces(eid); ++i )
        original_perimeter_faces.push_back( model_domain.PerimeterFace(eid,i) );

    
    // TESTING CELL VECTOR REBUILD
    
    // deleting central element (3) from model (node number stays the same)
    _test( model_domain.IsPerimeterCell( model_domain.E(3) ) == false );
    model.Mesh().Delete( next(model_domain.CellVector().begin(),3) );
    model_domain.RebuildCellAndPerimeterFaceVector();
    _test( model_domain.Cells() == model.Mesh().Elements() );
    _test( static_cast<long>(model_domain.InteriorCells()) == n_interior_elmts - 1 );
    _test( static_cast<long>(model_domain.PerimeterCells()) == n_perimeter_elmts );
    // the perimeter face vector should still be valid
    vector<uint32_t> perimeter_faces;
    perimeter_faces.reserve( model_domain.PerimeterCells() );
    for ( size_t eid{model_domain.InteriorCells()}; eid<model_domain.Cells(); ++eid )
      for ( uint32_t i{0u}; i<model_domain.PerimeterFaces(eid); ++i )
        perimeter_faces.push_back( model_domain.PerimeterFace(eid,i) );
    _test( perimeter_faces == original_perimeter_faces );

    // deleting perimeter element (5) from model (node number stays the same)
    _test( model_domain.IsPerimeterCell( model_domain.E(5) ) == true );
    model.Mesh().Delete( next(model_domain.CellVector().begin(),5) );
    model_domain.RebuildCellAndPerimeterFaceVector();
    _test( model_domain.Cells() == model.Mesh().Elements() );
    _test( static_cast<long>(model_domain.InteriorCells()) == n_interior_elmts - 1 );
    _test( static_cast<long>(model_domain.PerimeterCells()) == n_perimeter_elmts - 1 );
    // has the perimeter face vector been rebuild?
    vector<uint32_t> perimeter_faces2;
    perimeter_faces2.reserve( model_domain.PerimeterCells() );
    for ( size_t eid{model_domain.InteriorCells()}; eid<model_domain.Cells(); ++eid )
      for ( uint32_t i{0u}; i<model_domain.PerimeterFaces(eid); ++i )
        perimeter_faces2.push_back( model_domain.PerimeterFace(eid,i) );
    _test( perimeter_faces2 != original_perimeter_faces );


    // TESTING NODE VECTOR REBUILD
    
    // deleting the next perimeter element
    _test( model_domain.IsPerimeterCell( model_domain.E(4) ) == true );
    _test( model_domain.IsPerimeterCell( model_domain.E(6) ) == true );
    model.Mesh().Delete( next(model_domain.CellVector().begin(),4) );
    model.Mesh().Delete( next(model_domain.CellVector().begin(),6) );
    model_domain.RebuildCellAndPerimeterFaceVector();
    // after rebuild there should be a single orphan node (5) that needs to be deleted
    model.Mesh().Delete( next(model_domain.NodeVector().begin(),5) );
    model_domain.RebuildNodeVector();
    _test( model_domain.Nodes() == model.Mesh().Nodes() );
    // add further checks here
    
 } // end Test_RebuildCellAndNodeVectors





/** compares node locations and connectivity
    - uses point locations for node comparison

*/
template<uint32_t dim,template<uint32_t> class simplicial_complex>
bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>& domain1,
                                                  const ModelSubDomain<dim,simplicial_complex>& domain2,
                                                  bool verbose )
 {
    // domain names
    if ( domain1.Name() != domain2.Name() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: names mismatch: "<< domain1.Name() <<" vs. "<< domain2.Name() <<"\n";
         //return false;
      }
    _test( domain1.Name() == domain2.Name() );
   
    // 1. nodes
    // --------
    // number of nodes
    if ( domain1.Nodes() != domain2.Nodes() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: node numbers do not match: "<< domain1.Nodes() <<" vs "<< domain2.Nodes() <<"\n";
         //return false;
      }
    _test( domain1.Nodes() == domain2.Nodes() );
   
    // interior nodes
    if ( domain1.InteriorNodes() != domain2.InteriorNodes() ) {
         if ( verbose ) {
              cerr <<"\ncompareModelSubdomains: number of interior nodes does not match: "<< domain1.InteriorNodes() <<" vs "<< domain2.InteriorNodes() <<"\n";
              cerr <<"\tinterior nodes domain 1 vs domain 2:\n";
              for ( auto nit=domain1.NodesBegin(); nit!=domain1.PerimeterNodesBegin(); ++nit )
                cerr << (*nit)->Idx() <<" ";
              cerr << endl;
              for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); ++nit )
                cerr << (*nit)->Idx() <<" ";
              cerr << endl;
           }
         //return false;
      }
    _test( domain1.InteriorNodes() == domain2.InteriorNodes() );
   
    // node coordinates (sorted ranges of points have to be created)
    // ranges for domain 1
    set<Point<dim> >  interior_points1;
    for ( auto nit=domain1.NodesBegin(); nit!=domain1.PerimeterNodesBegin(); ++nit )
      interior_points1.insert( (*nit)->Coordinate() );
    set<Point<dim> >  perimeter_points1;
    for ( auto nit=domain1.PerimeterNodesBegin(); nit!=domain1.NodesEnd(); ++nit )
      interior_points1.insert( (*nit)->Coordinate() );
    // ranges for domain 2
    set<Point<dim> >  interior_points2;
    for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); ++nit )
      interior_points2.insert( (*nit)->Coordinate() );
    set<Point<dim> >  perimeter_points2;
    for ( auto nit=domain2.PerimeterNodesBegin(); nit!=domain2.NodesEnd(); ++nit )
      interior_points2.insert( (*nit)->Coordinate() );
    // comparisons
    if ( interior_points1 != interior_points2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: locations of interior points do not match.\n";
         //return false;
      }
    _test( interior_points1 == interior_points2 );

    if ( perimeter_points1 != perimeter_points2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: locations of perimeter points do not match.\n";
         //return false;
      }
    _test( perimeter_points1 == perimeter_points2 );
   
    interior_points1.clear();
    interior_points2.clear();
    perimeter_points1.clear();
    perimeter_points2.clear();

    // node boundary flags
    // (creating maps where the nodes are ordered by their location because
    //  they do not necessarily have the same number)
    map<Point<dim>,BOX_BOUNDARY>  boundary_flags1;
    for ( auto nit=domain1.NodesBegin(); nit!=domain1.NodesEnd(); ++nit )
      boundary_flags1.insert( make_pair( (*nit)->Coordinate(), (*nit)->AtBoundary() ) );
    //
    map<Point<dim>,BOX_BOUNDARY>  boundary_flags2;
    for ( auto nit=domain2.NodesBegin(); nit!=domain2.NodesEnd(); ++nit )
      boundary_flags1.insert( make_pair( (*nit)->Coordinate(), (*nit)->AtBoundary() ) );
    //
    // may fail due to tolerance differences in point classification _test( boundary_flags1.size() == boundary_flags2.size() );
    int node_flag_mismatches(0);
    auto nit2=boundary_flags2.begin();
    for ( auto nit1=boundary_flags1.begin(); (nit1!=boundary_flags1.end() && nit2!=boundary_flags2.end()); ++nit1, ++nit2 )
      if ( (*nit1).second != (*nit2).second ) {
            cerr <<"\n\t"<< parseBoundary((*nit1).second) <<" vs "<< parseBoundary((*nit2).second);
            node_flag_mismatches++;
        }
    if ( node_flag_mismatches > 0 ) {
          cerr <<"\ncompareModelSubdomains: the BOX_BOUNDARY flags of the domains do not match.\n";
          //return false;
      }
    _test( node_flag_mismatches == 0 );
    boundary_flags1.clear();
    boundary_flags2.clear();

   
    // 2. elements, i.e. connectivity
    // ------------------------------
    // number of elements
    if ( domain1.Cells() != domain2.Cells() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: number of elements does not match.\n";
         //return false;
      }
    _test( domain1.Cells() == domain2.Cells() );
   
    // number of interior elements
    if ( domain1.PerimeterCells() != domain2.PerimeterCells() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: mismatch in number of perimeter elements.\n";
         //return false;
      }
    _test( domain1.PerimeterCells() == domain2.PerimeterCells() );
   
    // are the interior elements the same ?
    set<size_t>  interior_elmts1, interior_elmts2;
    for ( auto it=domain1.CellsBegin(); it!=domain1.PerimeterCellsBegin(); ++it ) interior_elmts1.insert( (*it)->Idx() );
    for ( auto it=domain2.CellsBegin(); it!=domain2.PerimeterCellsBegin(); ++it ) interior_elmts2.insert( (*it)->Idx() );
    _test( interior_elmts1 == interior_elmts2 );
   
    // element connectivity (not assuming that elements are in same order)
    set<vector<size_t> > plist_entries1;
    for ( auto it=domain1.CellsBegin(); it!=domain1.CellsEnd(); ++it ) {
         vector<size_t> nodes( (*it)->Nodes() );
         for ( uint32_t i{0U}; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
         plist_entries1.insert( std::move(nodes) );
      }
    set<vector<size_t> > plist_entries2;
    for ( auto it=domain2.CellsBegin(); it!=domain2.CellsEnd(); ++it ) {
         vector<size_t> nodes( (*it)->Nodes() );
         for ( uint32_t i{0U}; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
//cerr <<"\n"<< (*it)->Idx() <<": ";
//out(nodes);
         plist_entries2.insert( std::move(nodes) );
      }
    // comparing plists
    bool plists_are_the_same(true);
    if ( plist_entries1 != plist_entries2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: (plist) mismatch in element-to-node connectivity.\n";
         plists_are_the_same = false;
         //return false;
      }
    _test( plists_are_the_same );

    plist_entries1.clear();
    plist_entries2.clear();
   
    // comparing the element neighbor connectivity
    set<vector<int64_t> > pfverts_entries1;
    for ( auto it=domain1.CellsBegin(); it!=domain1.CellsEnd(); ++it ) {
         vector<int64_t> nbors( (*it)->Neighbors(),0 );
         for ( uint32_t i{0}; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = static_cast<long>((*it)->Neighbor(i)->Idx());
             }
         pfverts_entries1.insert( std::move(nbors) );
      }
    set<vector<int64_t> > pfverts_entries2;
    for ( auto it=domain2.CellsBegin(); it!=domain2.CellsEnd(); ++it ) {
         vector<int64_t> nbors( (*it)->Neighbors(),0 );
         for ( uint32_t i{0}; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = static_cast<long>((*it)->Neighbor(i)->Idx());
             }
         pfverts_entries2.insert( std::move(nbors) );
      }
    // comparing neighbor connectivity lists
    bool pfverts_are_the_same(true);
    if ( pfverts_entries1 != pfverts_entries2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: (pfverts) mismatch in element to neighbor-element connectivity.\n";
         pfverts_are_the_same = false;
         //return false;
      }
    _test( pfverts_are_the_same );

    plist_entries1.clear();
    plist_entries2.clear();
   
    // all comparisons passed - domains have same nodes and elements.
    return true;
   
 } // end compareModelSubdomains

template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<3U,Element>&,
                                                           const ModelSubDomain<3U,Element>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<2U,Element>&,
                                                           const ModelSubDomain<2U,Element>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<1U,Element>&,
                                                           const ModelSubDomain<1U,Element>&, bool );

template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<3U,Face>&,
                                                           const ModelSubDomain<3U,Face>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<2U,Face>&,
                                                           const ModelSubDomain<2U,Face>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<1U,Face>&,
                                                           const ModelSubDomain<1U,Face>&, bool );





/**
    Using vsets from 'vset_makers' as input data, this tests weither     ModelSubDoman::EstablishNeighborConnectivity()  recreates the correct neighbor connectivity
    
        @author SKM 5/12/20
        
    @note this actually tests functionality of the MeshManager

*/
bool ModelSubDomain_Test::Test_EstablishNeighborConnectivity() 
 {
    // 0. creating the test model
    // --------------------------
    string varFileName("ModelSubDomain_Test-variables.txt");
    const bool   skewed_elements(false); // otherwise model is not a box anymore
    VSet<3U>     vset;
    create_Pyramid_Hexa_VSet( vset, skewed_elements );
    // adding 'node number' as a variable
    PropertyData node_nums( NODE, SCALAR, 3U );
    node_nums.Reserve( vset.Vertices() );
    for ( auto i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );
    // adding 'element number' as a variable
    PropertyData elmt_nums( ELEMENT, SCALAR, 3U );
    node_nums.Reserve( vset.Elements() );
    for ( size_t i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );

    Model<3U> model( vset, varFileName.c_str() );

    // 1. recreating the neighbor connectivity and comparing
    // -----------------------------------------------------
    const csmp::Index eid_key(model.Database().StorageKey("element number"));
    Region<3U>& domain = model.Region("Model");
    
    model.Mesh().UpdateConnectivity();
    //           ^^^^^^^^^^^^^^^^^^^^^
    bool no_mismatch(true);
    for ( auto it=domain.CellsBegin(); it!=domain.CellsEnd(); ++it )
      for ( uint32_t i{0}; i<(*it)->Neighbors(); ++i ) {
           if ( (*it)->Neighbor(i) != nullptr ) {
                const auto elmt_id = static_cast<size_t>((*it)->Read( eid_key ));
                const auto nbor_id = static_cast<long>((*it)->Neighbor(i)->Read( eid_key ));
                _test( nbor_id == vset.Pfvert( elmt_id, i ) );
                if ( nbor_id != vset.Pfvert( elmt_id, i ) ) {
                     cerr <<"\nelmt "<< elmt_id <<":"<< i <<": vset vs. reconstructed neighbor: ";
                     cerr << vset.Pfvert( elmt_id, i ) <<" vs. "<< nbor_id;
                     no_mismatch = false;
                  }
             }
           if ( no_mismatch == false ) cerr << endl;
        }

    return no_mismatch;
    
 } // end Test_EstablishNeighborConnectivity
     




/**
    std::string Name() const noexcept;
    void Name( const std::string& ) noexcept;
    
    void ScheduleForRebuild() noexcept;
    bool NeedsRebuild() const noexcept;

    void CreateNodePointerVector();
    void CreateNodePointerVector( std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& contacting_cells );
    void SortVectors( size_t interior_cells, size_t interior_nodes );

    int32_t DomainIndex() const noexcept;
    virtual size_t  RenumberNodes() const noexcept;
    size_t  RenumberCells() const noexcept;
    void    UpdateMemberIndexes() const noexcept;
    std::vector<size_t>  MemberCellIndexes() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator  PerimeterCellsBegin() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator  CellsEnd() const noexcept;

    bool IsUnique()
    void IsUnique( bool unique_domain )
    bool IsContiguous() const;
    std::pair<CELL_SHAPE,bool>  SingleCellShapeDomain() const;
    */
void ModelSubDomain_Test::Test_Basics()
 {
     // 2D model with everything: regions, boundaries and split boundaries
     {
       VSet<2>         vset;
       ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
       Model<2>       model( topo, vset, "CSMP-1phase-variables.txt", false /* use regions file if any */ );
       model.Name( "SPLIT22_BASIC" );
       
       if ( verbose_ ) {
           model.RegionsOut();
           model.BoundariesOut();
           model.SplitBoundariesOut();
         }
       
       Region<2>& model_domain = model.Region("Model"); // domain idx = 1 (range 1..n!)
       _test( model_domain.DomainIndex() == 1 );
       
       Region<2>& lower = model.Region("lower"); // domain idx = 2
       Region<2>& upper = model.Region("upper"); // domain idx = 3
       _test( lower.DomainIndex() == 2 );
       _test( upper.DomainIndex() == 3 );
       // is this correctly conveyed to the elements of the domain?
       // DO NOT USE DomainIndex() for non-unique regions! - for ( const auto& elmt: model_domain.CellVector() ) _test( elmt->Region_ID() == 1 );
       for ( const auto& elmt: lower.CellVector() ) _test( elmt->Region_ID() == 2 );
       for ( const auto& elmt: upper.CellVector() ) _test( elmt->Region_ID() == 3 );
      
       Boundary<2>& top = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)
       _test( top.DomainIndex() == 4 );
       
       SplitBoundary<2>& horizontal_splitboundary = model.SplitBoundary("horizontal_splitboundary"); // domain idx = 1
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2
       _test( horizontal_splitboundary.DomainIndex() == 1 );
       _test( inclined_split_boundary.DomainIndex() == 2 );
       
       // removal and re-creation of Model. Does it affect the domain indices?
       model.RemoveRegion("Model");
       model.FormModelRegion( false /* is_unique */ );
       _test( horizontal_splitboundary.DomainIndex() == 1 );
       _test( inclined_split_boundary.DomainIndex() == 2 );
       for ( const auto& elmt: lower.CellVector() ) _test( elmt->Region_ID() == 2 );
       for ( const auto& elmt: upper.CellVector() ) _test( elmt->Region_ID() == 3 );
       
      // TODO: verify that it is universally true for multiply modified models that such operations have not invalidated unique domain-indices and element region IDs
      
       // testing Name()
       _test( lower.Name() == "lower" );
       _test( top.Name() == "TOP" );
       _test( horizontal_splitboundary.Name() == "horizontal_splitboundary" );
       
       horizontal_splitboundary.Name("slit1");
       _test( horizontal_splitboundary.Name() == "slit1" );
       horizontal_splitboundary.Name("horizontal_splitboundary");
       
       // Rebuild flaggin
       _test( lower.NeedsRebuild() == false );
       lower.ScheduleForRebuild();
       _test( lower.NeedsRebuild() == true );
       
       // Node vector creation and rebuilding
       vector<Node<2>*> nptrs0 = upper.NodeVector();
       upper.RebuildNodeVector(); // sorting should be preserved
       vector<Node<2>*> nptrs1 = upper.NodeVector();
      _test( nptrs0 == nptrs1 );
       
       // loosing the original sorting
       upper.CreateNodePointerVector();
       nptrs1 = upper.NodeVector();
       sort( nptrs0.begin(), nptrs0.end() );
       sort( nptrs1.begin(), nptrs1.end() );
       _test( nptrs0 == nptrs1 );
            
       // SortVector
       size_t inner_cells = upper.InteriorCells();
       size_t inner_nodes = upper.InteriorNodes();
       upper.SortVectors( inner_cells, inner_nodes );
       _test( upper.InteriorCells() == inner_cells );
       _test( upper.InteriorNodes() == inner_nodes );
       
       // DomainIndex
       _test( lower.DomainIndex() == 2 ); // first unique region created
       _test( upper.DomainIndex() == 3 ); // second unique region created
       // boundaries
       _test( top.DomainIndex() == 4 );
       // split boundaries
       _test( horizontal_splitboundary.DomainIndex() == 1 );
       
       // RenumberNodes()
       lower.RenumberNodes();
       // pipe the nodes into a transform and then use max
       auto idx_view  = nptrs1 | views::transform([](auto* n) { return n->Idx(); });
       auto max_n_idx = ranges::max(idx_view);
       _test( nptrs1.size() == max_n_idx+1 );
       
      // RenumberCells()
      lower.RenumberCells();
      auto cidx_view  = lower.CellVector() | views::transform([](auto* n) { return n->Idx(); });
      auto max_c_idx = ranges::max(cidx_view);
      _test( lower.Cells() == max_c_idx+1 );

      // UpdateMemberIndexes (cells and nodes)
      upper.UpdateMemberIndexes();
      vector<size_t> cell_indices = upper.MemberCellIndexes();
      _test( cell_indices.size() == upper.Cells() );
      max_c_idx = *max_element( cell_indices.begin(), cell_indices.end() );
      _test( max_c_idx == upper.Cells()-1 );
      
      // Iterator - sorted vector consistency
      _test( (*lower.PerimeterCellsBegin()) == lower.E( lower.InteriorCells() ) );
      _test( (*lower.PerimeterNodesBegin()) == lower.N( lower.InteriorNodes() ) );
    }
    
    // testing CreateNodeVector( matching_cells )
    {
       ANSYS_Model3D model3D( "BoxHalfs3D", "CSMP-variables.txt" );
       Region<3>&    mref3D( model3D.Region( "HALF" ) ); // lower-dim dividing wall
       
       vector<pair<pair<Element<3>*,uint32_t>,pair<Element<3>*,uint32_t> > > contacting_cells;
       size_t shared_cell_pairs = sharedPerimeterCells( model3D.Region("MATRIX_LEFT"), model3D.Region("MATRIX_RIGHT"), contacting_cells );
                              
       _test( shared_cell_pairs == mref3D.Cells() );
       
       vector<Node<3>*> node_vec0 = mref3D.NodeVector();
       mref3D.CreateNodePointerVector( contacting_cells );
       //     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       vector<Node<3>*> node_vec1 = mref3D.NodeVector();

       sort( node_vec0.begin(), node_vec0.end() );
       sort( node_vec1.begin(), node_vec1.end() );

       _test( node_vec0 == node_vec1 );
    }

     
 } // end ModelSubDomain_Test::Test_Basics





/**
    Prints VTU file where topotype identifiers have been converted to number.
 */
template<uint32_t dim>
void ModelSubDomain_Test::BoundaryAndTopoTypeFlagsToVTU( Model<dim>& model )
 {
    model.CreateProperty( "topo type", "TT", "none");
    model.CreateProperty( "box flag", "BF", "none");
    
    // getting diagnostic values
    boundaryFlagToNumber( model, "box flag" );
    topoTypeToNumber( model, "topo type" );
 
    VTU_Interface<dim>  vtu(model);
    list<string> outvars{ "box flag", "topo type" };
    vtu.OutputDataToVTU( "BoundaryAndTopoTypeFlagsToVTU", outvars, model.Region("Model"), 0 );
    
    model.DeleteProperty( "topo type" );
    model.DeleteProperty( "box flag" );
 }

template void ModelSubDomain_Test::BoundaryAndTopoTypeFlagsToVTU( Model<3>& );
template void ModelSubDomain_Test::BoundaryAndTopoTypeFlagsToVTU( Model<2>& );
template void ModelSubDomain_Test::BoundaryAndTopoTypeFlagsToVTU( Model<1>& );




/**
    void IdentifyPerimeter();
    void UpdateTopoTypeNodeFlags();
    void BuildPerimeterFaceVector( int64_t interior_cells );
    void RebuildSubDomainAfterChangeOfCellVector();
    void RebuildSubDomainAfterChangeOfNodeVector();

    void UpdateCellMembershipApplyingConstraints( typename std::vector<CELL<dim>*>::const_iterator master_domain_start,
                                                  typename std::vector<CELL<dim>*>::const_iterator master_domain_end,
                                                  const PropertyConstraints& );
    size_t RebuildCellAndPerimeterFaceVector();
 */
void ModelSubDomain_Test::Test_SubDomainConstructionMethods()
 {
    // 2D Testing: model with everything: regions, boundaries and split boundaries
     {
       VSet<2>        vset;
       ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
       Model<2>       model( topo, vset, "CSMP-1phase-variables.txt", false /* use regions file if any */ );
       model.Name( "SPLIT22_BASIC" );
       
       if ( verbose_ ) {
           model.RegionsOut();
           model.BoundariesOut();
           model.SplitBoundariesOut();
         }
      
       // verifying that the nodes in the manifolds are not connected to each other
       for ( auto nmf=model.Mesh().NodeManifoldsBegin(); nmf!=model.Mesh().NodeManifoldsEnd(); ++nmf )
         _test( (*nmf).InterConnectedMemberNodes().second == false );
       
       // accessors
       Region<2>& model_domain = model.Region("Model"); // domain idx = 1 (range 1..n!)
       Region<2>& lower = model.Region("lower"); // domain idx = 2
       Region<2>& upper = model.Region("upper"); // domain idx = 3
      
       Boundary<2>& top = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)

       // visualisaton
       if ( verbose_ )
         {
           BoundaryAndTopoTypeFlagsToVTU( model );
           VTU_Interface<2> vtu_output( model );
           list<string>     output_props{"node number"};
           vtu_output.OutputDataToVTU( "split22_basic_", output_props, model_domain, 0 );
           SplitBoundary<2>& horizontal_splitboundary = model.SplitBoundary("horizontal_splitboundary"); // domain idx = 1
           SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2
           vtu_output.OutputDataToVTU( "split22_basic_", output_props, horizontal_splitboundary, 0 );
           vtu_output.OutputDataToVTU( "split22_basic_", output_props, inclined_split_boundary, 0 );
         }
         
       // Identification of cells and nodes that form the perimeter of the region
       model.Mesh().UpdateConnectivity(); // connectivity in VSet is broken
       
       // non-unique Region Model'
       model_domain.IdentifyPerimeter();
       //           ^^^^^^^^^^^^^^^^^
       _test( model_domain.InteriorNodes()  ==  4 );
       _test( model_domain.PerimeterNodes() == 31 );
       // checking the perimeter nodes using BOX_BOUNDARY flags
       for ( auto nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
            _test( (*nit)->AtBoundary() != NOT );
            if ( verbose_ && (*nit)->AtBoundary() == NOT ) printNodeAttibutes<2>( nit, nit+1 );
         }
       // checking interior nodes
       for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.PerimeterNodesBegin(); ++nit )
         // not on an external or an internal boundary
         _test( (*nit)->AtBoundary() == NOT && !(*nit)->IsManifold() ); // || (*nit)->AtBoundary() == INTERNAL );
       
       // testing
       const csmp::Index nkey = model.Database().StorageKey("node number");
       
       // TOPOTYPE flags: are they the same as in the input VSet?
       // -------------------------------------------------------
       // get unperturbed gflags for testing
       vector<TOPOTYPE> gflags0; gflags0.reserve( model_domain.Nodes() );
       for ( const auto& node : model_domain.NodeVector() ) gflags0.push_back( node->Attribute() );
       
       // were any errors made during model construction?
       // (NOTE: TOPOTYPE data exactly as in the VSet because VData::InitialiseNodeTopologyIdentifiers() is not called)
       if ( verbose_ ) cout <<"\n"<<"classifiers in VData vs. TOPOTYPE attribute in Model:";
       for ( size_t n{0}; n<vset.Vertices(); ++n ) {
           // finding the node with the matching node by number (that was stored in the VSet)
           const Node<2>* matching_node = nullptr;
           for ( const auto& node : model_domain.NodeVector() )
             if ( static_cast<size_t>(node->Read(nkey)) == n ) {
                  matching_node = node;
                  break;
             }
           assert( matching_node != nullptr );
           if ( verbose_ ) {
                cout <<"\n\t"<<"Node "<< n <<": "
                     << parseBoundary(static_cast<BOX_BOUNDARY>(vset.BFlag(n))) <<": "
                     << parseTopology( static_cast<TOPOTYPE>(vset.BREP_Flag(n)) )
                     <<":"<< parseTopology( matching_node->Attribute() );
             }
           _test( vset.BREP_Flag(n) == matching_node->Attribute() );
         }
       
       model_domain.UpdateTopoTypeNodeFlags();
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^
 
       // did the gflags change
       vector<TOPOTYPE> gflags1; gflags1.reserve( model_domain.Nodes() );
       for ( const auto& node : model_domain.NodeVector() ) gflags1.push_back( node->Attribute() );
       _test( gflags0 == gflags1 );

       // perimeter face vector (in presence of boundaries and split boundaries
       // ---------------------------------------------------------------------
       // equidimensional region
       // ----------------------
       // existing perimeter face vector
       for ( size_t f=lower.InteriorCells(); f<lower.Cells(); ++f )
         for ( uint32_t i{0}; i<lower.PerimeterFaces(f); ++i )
             // must either be a nullptr or a cell that is not part of the region
              _test( lower.E(f)->Neighbor( lower.PerimeterFace(f,i) ) == nullptr ||
                     lower.Contains( lower.E(f)->Neighbor( lower.PerimeterFace(f,i) ) ) == false );
                     
       // rebuild version
       lower.BuildPerimeterFaceVector( static_cast<int64_t>(lower.InteriorCells()) );
       //    ^^^^^^^^^^^^^^^^^^^^^^^^^
       for ( size_t f=lower.InteriorCells(); f<lower.Cells(); ++f )
         for ( uint32_t i{0}; i<lower.PerimeterFaces(f); ++i )
             // must either be a nullptr or a cell that is not part of the region
              _test( lower.E(f)->Neighbor( lower.PerimeterFace(f,i) ) == nullptr ||
                     lower.Contains( lower.E(f)->Neighbor( lower.PerimeterFace(f,i) ) ) == false );
                     
                     
       // rebuilding regions
       // ------------------
       // (no need to change cell membership because methods rebuild from scratch)
       
       // region upper
       // ------------
       auto copy_elmt_vec = upper.CellVector();
       auto copy_node_vec = upper.NodeVector();
       auto n_perim_cells = upper.PerimeterCells();
       auto n_perim_nodes = upper.PerimeterNodes();
                     
       upper.RebuildSubDomainAfterChangeOfCellVector();
       //    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( copy_elmt_vec == upper.CellVector() );
       _test( copy_node_vec == upper.NodeVector() );
       _test( n_perim_cells == upper.PerimeterCells() );
       _test( n_perim_nodes == upper.PerimeterNodes() );
       
       auto cells_removed = upper.RebuildCellAndPerimeterFaceVector();
       //                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( cells_removed == 0 );
       _test( copy_elmt_vec == upper.CellVector() );
       _test( copy_node_vec == upper.NodeVector() );
       _test( n_perim_cells == upper.PerimeterCells() );
       _test( n_perim_nodes == upper.PerimeterNodes() );
       
       // boundary top
       // ------------
       auto copy_face_vec = top.CellVector();
       copy_node_vec      = top.NodeVector();
       auto n_perim_faces = top.PerimeterCells();
       n_perim_nodes      = top.PerimeterNodes();
                     
       upper.RebuildSubDomainAfterChangeOfNodeVector();
       //    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( copy_face_vec == top.CellVector() );
       _test( copy_node_vec == top.NodeVector() );
       _test( n_perim_faces == top.PerimeterCells() );
       _test( n_perim_nodes == top.PerimeterNodes() );

       // changing region to only include cells with properties in given range
       // --------------------------------------------------------------------
       PropertyConstraints constraints( model.Database(), "element number", 14, 16 );
       model_domain.RenumberCells();
       
       lower.UpdateCellMembershipApplyingConstraints( model_domain.CellsBegin(), model_domain.CellsEnd(), constraints );
       //    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( lower.IsUnique() == false );
       _test( lower.Cells() == 3 );
       _test( lower.PerimeterNodes() == 9 );
    }

 } // end Test_SubDomainConstructionMethods





/**
    bool              Empty() const noexcept;
    size_t            Nodes() const noexcept;
    size_t            InteriorNodes() const noexcept;
    size_t            PerimeterNodes() const noexcept;
    size_t            IntegrationPoints() const noexcept;
    size_t            SectorIntegrationPoints() const noexcept;
    size_t            FacetIntegrationPoints() const noexcept;
    size_t            Cells() const noexcept;
    size_t            InteriorCells() const noexcept;
    size_t            PerimeterCells() const noexcept;
    bool              Contains( const CELL<dim>* const ) const noexcept;
    bool              Contains( const Node<dim>* const ) const noexcept;
    bool              IsPerimeterNode( const csmp::Node<dim>* const ) const noexcept;
    bool              IsPerimeterCell( const CELL<dim>* const ) const noexcept;
    uint32_t        PerimeterFaces( size_t cell_idx ) const;
    uint32_t        PerimeterFace( size_t cell_idx, uint32_t face ) const;
    bool              IsPerimeterNode( const size_t nidx ) const;
    bool              IsPerimeterCell( const size_t eidx ) const;
    size_t            SharedPerimeterNodes( start, end )

    csmp::Node<dim>*  N( size_t n ) const noexcept;
    CELL<dim>*        E( size_t n ) const noexcept;

*/
void ModelSubDomain_Test::Test_SubDomainDiagnostics()
 {
    // 2D Testing: model with everything: regions, boundaries and split boundaries
     {
       VSet<2>         vset;
       ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
       Model<2>       model( topo, vset, "CSMP-1phase-variables.txt", false /* use regions file if any */ );
       model.Name( "SPLIT22_BASIC" );
       
       if ( verbose_ ) {
           model.RegionsOut();
           model.BoundariesOut();
           model.SplitBoundariesOut();
         }
       
       // accessors
       Region<2>& model_domain = model.Region("Model"); // domain idx = 1 (range 1..n!)
       Region<2>& lower = model.Region("lower"); // domain idx = 2
       Region<2>& upper = model.Region("upper"); // domain idx = 3
      
       Boundary<2>& top   = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)
       Boundary<2>& right = model.Boundary("RIGHT");
       
       //SplitBoundary<2>& horizontal_splitboundary = model.SplitBoundary("horizontal_splitboundary"); // domain idx = 1
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2


       // testing using hand-coded model SPLIT22_BASIC (vsetMakers create_BoundarySplitBoundaryPatch() )
       
       // Non-unique Region 'Model'
       // -------------------------
       _test( model_domain.Empty() == false );
       _test( model_domain.Nodes() == 35 );
       _test( model_domain.InteriorNodes() == 4 );
       if ( model_domain.InteriorNodes() != 4 && verbose_ ) printNodeAttibutes<2>( model_domain.NodesBegin(), model_domain.NodesEnd() );
       _test( model_domain.PerimeterNodes() == 31 ); // fail->32 interesting question: are nodes at interior boundaries added to perimeter?
       _test( model_domain.IntegrationPoints() == (16 * 4 + 3 * 3) );       // 16 quads + 3 triangles
       model.Mesh().InitializeFiniteVolumeStencils( model.Database(), true /* assign_to_elmts */ );
       _test( model_domain.SectorIntegrationPoints() == (16 * 4 + 3 * 3) ); // elmts * nodes
       _test( model_domain.FacetIntegrationPoints() == (16 * 4 + 3 * 3) );  // elements * facets per element
       _test( model_domain.Cells() == 19 ); // elements
       _test( model_domain.InteriorCells() == 2 );
       _test( model_domain.PerimeterCells() == 17 );
       model_domain.UpdateMemberIndexes();
       _test( model_domain.Contains( model_domain.E(5) ) == true );
       _test( model_domain.Contains( model_domain.N(3) ) == true );
       _test( model_domain.IsPerimeterCell( 0ul ) == false );
       _test( model_domain.IsPerimeterCell( (*next(model_domain.PerimeterCellsBegin(),1)) ) == true );
       _test( model_domain.IsPerimeterNode( 1ul ) == false );
       _test( model_domain.IsPerimeterNode( model_domain.Nodes()-1 ) == true );
       _test( model_domain.IsPerimeterNode( (*next(model_domain.PerimeterNodesBegin(),2)) ) == true );
       // corner quad has two perimeter faces
       _test( model_domain.PerimeterFaces( model_domain.InteriorCells()) == 2 );
       // _test( model_domain.PerimeterFaces(0) == 2 ); - does not crash
       // in positions 0 & 3
       _test( model_domain.PerimeterFace( model_domain.InteriorCells(), 0 ) == 0 );
       _test( model_domain.PerimeterFace( model_domain.InteriorCells(), 1 ) == 3 );
       _test( model_domain.SharedPerimeterNodes( lower.PerimeterNodesBegin(), lower.NodesEnd() ) == 14 );

       // Unique Region 'Upper' (intersected by SB)
       // -----------------------------------------
       // tested: OK SKM 8/2/26
       _test( upper.Empty() == false );
       _test( upper.Nodes() == 21 );
       _test( upper.InteriorNodes() == 3 );
       if ( upper.InteriorNodes() != 3 && verbose_ ) printNodeAttibutes<2>( upper.NodesBegin(), upper.NodesEnd() );
       _test( upper.PerimeterNodes() == 18 ); // fail->32 interesting question: are nodes at interior boundaries added to perimeter?
       _test( upper.IntegrationPoints() == (10 * 4 + 2 * 3) );       // 10 quads + 2 triangles
       model.Mesh().InitializeFiniteVolumeStencils( model.Database(), true /* assign_to_elmts */ );
       _test( upper.SectorIntegrationPoints() == (10 * 4 + 2 * 3) ); // elmts * nodes
       _test( upper.FacetIntegrationPoints() == (10 * 4 + 2 * 3) );  // elements * facets per element
       _test( upper.Cells() == 12 ); // elements
       _test( upper.InteriorCells() == 2 );
       _test( upper.PerimeterCells() == 10 );
       upper.UpdateMemberIndexes();
       _test( upper.Contains( upper.E(2) ) == true );
       _test( upper.Contains( upper.N(3) ) == true );
       _test( upper.IsPerimeterCell( 0ul ) == false );
       _test( upper.IsPerimeterCell( (*next(upper.PerimeterCellsBegin(),1)) ) == true );
       _test( upper.IsPerimeterNode( 1ul ) == false );
       _test( upper.IsPerimeterNode( upper.Nodes()-1 ) == true );
       _test( upper.IsPerimeterNode( (*next(upper.PerimeterNodesBegin(),2)) ) == true );
       // corner quad has two perimeter faces
       _test( upper.PerimeterFaces( upper.InteriorCells()) == 2 );
       _test( upper.PerimeterFace( upper.InteriorCells(), 0 ) == 0 );
       _test( upper.PerimeterFace( upper.InteriorCells(), 1 ) == 3 );
       _test( upper.SharedPerimeterNodes( lower.PerimeterNodesBegin(), lower.NodesEnd() ) == 0 );

       // Boundary 'TOP'
       // -----------------------------------------
       // tested:
       _test( top.Empty() == false );
       _test( top.Nodes() == 4 );
       _test( top.InteriorNodes() == 2 );
       if ( top.InteriorNodes() != 2 && verbose_ ) printNodeAttibutes<2>( top.NodesBegin(), top.NodesEnd() );
       _test( top.PerimeterNodes() == 2 ); // fail->32 interesting question: are nodes at interior boundaries added to perimeter?
       _test( top.IntegrationPoints() == (3 * 2) );       // 3 line elements
       model.Mesh().InitializeFiniteVolumeStencils( model.Database(), true /* assign_to_elmts */ );
       _test( top.SectorIntegrationPoints() == (3 * 2) ); // cells * nodes
       _test( top.FacetIntegrationPoints() == 3 );  // cells * facets per element
       _test( top.Cells() == 3 ); // elements
       _test( top.InteriorCells() == 1 );
       _test( top.PerimeterCells() == 2 );
       upper.UpdateMemberIndexes();
       _test( top.Contains( top.E(2) ) == true );
       _test( top.Contains( top.N(3) ) == true );
       _test( top.IsPerimeterCell( 1ul ) == true );
       _test( top.IsPerimeterCell( (*next(top.PerimeterCellsBegin(),1)) ) == true );
       _test( top.IsPerimeterNode( 1ul ) == false );
       _test( top.IsPerimeterNode( top.Nodes()-1 ) == true );
       _test( top.IsPerimeterNode( (*next(top.PerimeterNodesBegin(),1)) ) == true );
       // corner quad has two perimeter faces
       _test( top.PerimeterFaces( top.InteriorCells()) == 1 );
       _test( top.PerimeterFaces( top.Cells()-1) == 1 );
       _test( top.PerimeterFace( top.InteriorCells(), 0 ) == 1 );
       _test( top.SharedPerimeterNodes( right.PerimeterNodesBegin(), right.NodesEnd() ) == 1 );

       // SplitBoundary 'inclined_split_boundary'
       // -----------------------------------------
       // tested:
       _test( inclined_split_boundary.Empty() == false );
       _test( inclined_split_boundary.NodeManifolds().size() == 2 );
       _test( inclined_split_boundary.InsideNodes().first.size() == 5 ); // 3 + 2 at intersection with horiz. SB
       _test( inclined_split_boundary.OutsideNodes().first.size() == 5 );
       _test( inclined_split_boundary.IntegrationPoints() == (3 * 2) );       // 3 line cells
       model.Mesh().InitializeFiniteVolumeStencils( model.Database(), true /* assign_to_elmts */ );
       _test( inclined_split_boundary.SectorIntegrationPoints() == (3 * 2) ); // cells * nodes
       _test( inclined_split_boundary.FacetIntegrationPoints() == 3 );  // cells * facets per element
       _test( inclined_split_boundary.Cells() == 3 ); // InterFace objects
       _test( inclined_split_boundary.InteriorCells() == 1 );
       _test( inclined_split_boundary.PerimeterCells() == 2 );
       upper.UpdateMemberIndexes();
       _test( inclined_split_boundary.Contains( inclined_split_boundary.E(2) ) == true );
       _test( inclined_split_boundary.IsPerimeterCell( 1ul ) == true );
       // perimeter manifolds should only exist at model exterior
       // corner quad has two perimeter faces
       _test( inclined_split_boundary.PerimeterFaces( top.InteriorCells()) == 1 );
       _test( inclined_split_boundary.PerimeterFaces( top.Cells()-1) == 1 );
       _test( inclined_split_boundary.PerimeterFace( top.InteriorCells(), 0 ) == 1 );

       if ( verbose_ ) cout <<"\nThat's it."<< endl;
     }
 } // end (diagnostics)





/**
    Tests:
     
    std::pair<int32_t,int32_t>  SpatialDimensions() const;
    Point<dim> Centroid() const;
    Point<dim> CenterOfGravity( const csmp::Index& rho_key ) const;
    void  MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const;
    void AssignNodeCoordinatesTo( const char* vector_prop );
    void AssignNodeCoordinatesTo( const char* scalar_prop, char coord );
    void AssignCellCharacteristicsTo( const char* characteristic, const char* var );
    void NodeAttributesToCSV();
*/
void ModelSubDomain_Test::Test_GeometricOperations()
 {
    // 2D Testing: model with everything: regions, boundaries and split boundaries
     {
       VSet<2>         vset;
       ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
       Model<2>       model( topo, vset, "CSMP-1phase-variables.txt", false /* use regions file if any */ );
       model.Name( "SPLIT22_BASIC" );
       
       if ( verbose_ ) {
           model.RegionsOut();
           model.BoundariesOut();
           model.SplitBoundariesOut();
         }
       
       // accessors
       Region<2>&        model_domain = model.Region("Model"); // domain idx = 1 (range 1..n!)
       Boundary<2>&      top = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2
       
       // TESTING
       // all these should be single dimension subdomains
       pair<int32_t,int32_t> model_dim = model_domain.SpatialDimensions();
       _test( model_dim.first  == 1 ); // only surfaces
       _test( model_dim.second == 2 ); // surface

       pair<int32_t,int32_t> boundary_dim = top.SpatialDimensions();
       _test( boundary_dim.first  == 1 ); // only lines
       _test( boundary_dim.second == 1 ); // line

       pair<int32_t,int32_t> split_boundary_dim = inclined_split_boundary.SpatialDimensions();
       _test( split_boundary_dim.first  == 1 ); // only lines
       _test( split_boundary_dim.second == 1 ); // line
       
       // bounding box
       Point<2> xyz_min, xyz_max;
       model_domain.MinMaxCoordinates( xyz_min, xyz_max );
       Point<2> mid_diagonal = (xyz_max-xyz_min) / 2.;
       Point<2> centroid     = model_domain.Centroid();
       _equal( mid_diagonal[0], centroid[0], 0.01 ); // 1%
       _equal( mid_diagonal[1], centroid[1], 0.01 );
       
       const csmp::Index rho_key  = model.Database().StorageKey("concentration"); // scalar node variable used as proxy
       model_domain.InputPropertyValue("concentration", makeScalar(ANY,1.) );
       Point<2> center_of_gravity = model_domain.CenterOfGravity( rho_key );
       _equal( center_of_gravity[0], centroid[0], 0.01 ); // 1%
       _equal( center_of_gravity[1], centroid[1], 0.01 );
       
       model_domain.AssignNodeCoordinatesTo("nodal velocity");
       
       const bool print_maximum{true};
       model_domain.AssignNodeCoordinatesTo("concentration", 'x' );
       _equal( printRangeOfVariable( model, "concentration",print_maximum ), xyz_max[0], 1.0e-5 );
       model_domain.AssignNodeCoordinatesTo("concentration", 'y' );
       _equal( printRangeOfVariable( model, "concentration",print_maximum ), xyz_max[1], 1.0e-5 );
       
       model_domain.AssignCellCharacteristicsTo( "area", "porosity" );
       // print minimum
       _test( printRangeOfVariable( model, "porosity", false ) >= 0.6 );
      
       // writes to text file: model.Name() + "_node_attribute.csv"
       if ( verbose_ ) model_domain.NodeAttributesToCSV(); // OK=reads correctly with Excel and Paraview
       
    } // end 2D testing
    
 } // end Test_GeometricOperations
 
 
 
 
 
/**
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );
    void ChangePropertyStatus( const char* property, VARIABLE_FLAG new_status_of_scalar, SUBDOMAIN_PART=COMPLETE );
    void ChangePropertyStatusWhere( const char* property,
                                    VARIABLE_FLAG new_status_of_scalar,
                                    double min_value_to_change,
                                    double max_value_to_change );

    void ChangePropertyStatus( const char* property,
                               const std::vector<VARIABLE_FLAG>& new_status,
                               SUBDOMAIN_PART=COMPLETE );

     void ChangePropertyStatusWhere( const char* property,
                                    const std::vector<VARIABLE_FLAG>& new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    void ChangePropertyStatus( const char* property,
                               uint32_t component,
                               VARIABLE_FLAG new_status,
                               SUBDOMAIN_PART=COMPLETE );

    void ChangePropertyStatusWhere( const char* property,
                                    uint32_t component,
                                    VARIABLE_FLAG new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    VARIABLE_FLAG  PropertyStatus( const char* variable, SUBDOMAIN_PART flag=COMPLETE , uint32_t i=0 ) const;
    void MinMaxOf( const char* property,   double& gmin, double& gmax ) const;
    void MinMaxOf( const csmp::Index&,     double& gmin, double& gmax ) const;
    void   CopyReplace( const char* from, const char* to );
*/
void ModelSubDomain_Test::Test_PropertyManipulations()
{
    // 2D Testing: model with everything: regions, boundaries and split boundaries
     {
       VSet<2>        vset;
       ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
       // remove property
       vset.RemoveData("permeability");
       Model<2>       model( topo, vset, "CSMP-variables.txt", false /* don't regions file */ );
       model.Name( "SPLIT22_BASIC" );
       
       if ( verbose_ ) {
           model.RegionsOut();
           model.BoundariesOut();
           model.SplitBoundariesOut();
         }
       
       // accessors
       Region<2>& model_domain = model.Region("Model"); // domain idx = 1 (range 1..n!)
       Region<2>& lower = model.Region("lower"); // domain idx = 2
       Region<2>& upper = model.Region("upper"); // domain idx = 3
      
       Boundary<2>& top   = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)
       Boundary<2>& bottom = model.Boundary("BOTTOM");
       
       SplitBoundary<2>& horizontal_splitboundary = model.SplitBoundary("horizontal_splitboundary"); // domain idx = 1
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2
       
       // INITIALISING VARIABLES
       
       // regions
       model_domain.InputPropertyValue("nodal variable", makeScalar(ANY,0.), COMPLETE );
       //           ^^^^^^^^^^^^^^^^^^
       upper.InputPropertyValue("element tensor", makeTensor(ANY,ANY,1.0e-15,0.,0.,1.0e-16) );
       //    ^^^^^^^^^^^^^^^^^^
       upper.InputPropertyValue("element variable", makeScalar(ANY,0.15) );
       upper.InputPropertyValue("nodal variable",   makeScalar(ANY,0.), INTERIOR );
       upper.InputPropertyValue("nodal variable",   makeScalar(ANY,0.), FIELD_DATA, INTERIOR ); // same as previous
       
       lower.InputPropertyValue("element tensor",   makeTensor(ANY,ANY,1.0e-12,0.,0.,2.0e-12) );
       lower.InputPropertyValue("element variable", makeScalar(ANY,0.25) );
       lower.InputPropertyValue("nodal variable",   makeScalar(ANY,3.), INTERIOR );
       lower.InputPropertyValue("nodal variable",   makeScalar(ANY,3.), FIELD_DATA, INTERIOR ); // same as previous
       
       top.InputPropertyValue("nodal variable",    makeScalar(DIRICH,0.), COMPLETE );
       bottom.InputPropertyValue("nodal variable", makeScalar(DIRICH,5.), COMPLETE );
       
       // boundaries (methods work only for boundaries)
       top.InputPropertyValue("face variable",    makeScalar(DIRICH,7.), COMPLETE );
       bottom.InputPropertyValue("face variable", makeScalar(DIRICH,11.), COMPLETE );

       // split boundaries (methods work only for boundaries)
       horizontal_splitboundary.InputPropertyValue("interface variable", makeScalar(DIRICH,15.), INTERIOR );
       horizontal_splitboundary.InputPropertyValue("interface variable", makeScalar(DIRICH,25.), PERIMETER );
       _test( horizontal_splitboundary.PropertyStatus( "interface variable", COMPLETE ) == DIRICH );
       horizontal_splitboundary.InputPropertyValue("interface variable", makeScalar(FIELD_DATA,25.), PERIMETER );
       const csmp::Index& pkey = model.Database().StorageKey("interface variable");
       size_t interior_matches{0}, perimeter_matches{0};
       for ( const auto& iface : horizontal_splitboundary.CellVector() ) {
            if ( iface->Status(pkey) == DIRICH ) interior_matches++;
            else if ( iface->Status(pkey) == FIELD_DATA ) perimeter_matches++;
         }
       _test( interior_matches  == 2 );
       _test( perimeter_matches == 2 );

       inclined_split_boundary.InputPropertyValue("interface variable", makeScalar(DIRICH,35.), INTERIOR );
       inclined_split_boundary.InputPropertyValue("interface variable", makeScalar(DIRICH,45.), PERIMETER );
       _test( inclined_split_boundary.PropertyStatus( "interface variable", COMPLETE ) == DIRICH );

       // InputNodePropertyValue( const char* input_node_prop, const Var&, SUBDOMAIN_PART, INTERFACE_SIDE );
       horizontal_splitboundary.InputNodePropertyValue( "nodal variable", makeScalar(ANY,9.), INTERIOR, INSIDE );
       horizontal_splitboundary.InputNodePropertyValue( "nodal variable", makeScalar(ANY,9.), PERIMETER, INSIDE );
       // TODO: to enable testing with nodal variables override base-class method in SplitBoundary
      // throw _test( horizontal_splitboundary.PropertyStatus( "nodal variable", COMPLETE ) != ANY ); // only on inside!

       // InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );
       inclined_split_boundary.InputPropertyValue("interface variable", makeScalar(DIRICH,37.), DIRICH, COMPLETE );
       inclined_split_boundary.InputPropertyValue("interface variable", makeScalar(DIRICH,31.), ANY, COMPLETE );
       _test( horizontal_splitboundary.PropertyStatus( "interface variable", COMPLETE ) != DIRICH );
       
       // ChangeNodePropertyStatus( const char* property, VARIABLE_FLAG new_status_of_scalar, SUBDOMAIN_PART, INTERFACE_SIDE );
       horizontal_splitboundary.ChangeNodePropertyStatus( "nodal variable", NEUMANN, INTERIOR, OUTSIDE );
       horizontal_splitboundary.ChangeNodePropertyStatus( "nodal variable", NEUMANN, INTERIOR, INSIDE );
       // throw _test( horizontal_splitboundary.PropertyStatus( "nodal variable", COMPLETE ) == NEUMANN );
       
       // ChangeNodePropertyStatusWhere( const char*, VARIABLE_FLAG, INTERFACE_SIDE, double, double );
       horizontal_splitboundary.ChangeNodePropertyStatusWhere( "nodal variable", ROBIN, OUTSIDE, 11., 24. );
       // throw _test( horizontal_splitboundary.PropertyStatus( "nodal variable", COMPLETE ) != NEUMANN );
       // throw _test( horizontal_splitboundary.PropertyStatus( "nodal variable", COMPLETE ) != ROBIN );

       // TESTING
       
       // assignments
       const double tolerance{1.0e-13};
       double val_min, val_max;
       model_domain.MinMaxOf( "element tensor", val_min, val_max );
       //           ^^^^^^^^
       _test( approximatelyEqual(val_min,1.0e-16,tolerance) == true );
       _test( approximatelyEqual(val_max,2.0e-12) == true );
       
       const csmp::Index prop_key = model.Database().StorageKey("element variable");
       model_domain.MinMaxOf( prop_key, val_min, val_max );
       //           ^^^^^^^^
       _test( approximatelyEqual(val_min,0.15) == true );
       _test( approximatelyEqual(val_max,0.25) == true );
       
       model_domain.MinMaxOf( "face variable", val_min, val_max );
       //           ^^^^^^^^
       _test( approximatelyEqual(val_min,0.15) == true );
       _test( approximatelyEqual(val_max,0.25) == true );
       
       // unitialised variable
       model_domain.MinMaxOf( "element array", val_min, val_max );
       //           ^^^^^^^^
       _test( isnan(val_min) == true );
       _test( isnan(val_max) == true );
 
       model_domain.CopyReplace( "element variable", "element number" );
       //           ^^^^^^^^^^^
       model_domain.MinMaxOf( "element variable", val_min, val_max );
       _test( approximatelyEqual(val_min,0.15) == true );
       _test( approximatelyEqual(val_max,0.25) == true );

       // value ranges
       model_domain.MinMaxOf( "nodal variable", val_min, val_max );
       _test( approximatelyEqual(val_min,0.) );
       _test( approximatelyEqual(val_max,9.) );
      
       
       // TESTING FLAGS AND THEIR MANIPULATION
       
       _test( model_domain.PropertyStatus( "nodal variable", COMPLETE ) == ANY ); // scalar, return ANY when there are multiple flags
       //                  ^^^^^^^^^^^^^^
       _test( top.PropertyStatus( "nodal variable", COMPLETE ) == DIRICH ); // scalar
       _test( bottom.PropertyStatus( "nodal variable", COMPLETE ) == DIRICH ); // scalar
 
       // no more ANY
       model_domain.ChangePropertyStatus( "nodal variable", PLAIN );
       //           ^^^^^^^^^^^^^^^^^^^^
       _test( model_domain.PropertyStatus( "nodal variable" ) == PLAIN );
       
       
       double min_value_to_change{0.09}, max_value_to_change{0.16};
       //                                                          new flag    min-val-affected     max-val affected
       model_domain.ChangePropertyStatusWhere( "element variable", FIELD_DATA, min_value_to_change, max_value_to_change );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^
       size_t counter{0};
       for ( const auto& cell : model_domain.CellVector() ) if ( cell->Status(prop_key) == FIELD_DATA ) counter++;
       _test( counter == 12 );

       // changing all flags from default ANY
       model_domain.InputPropertyValue("element vector", makeVector(ANY,ANY,0.,0.), COMPLETE );
       model_domain.ChangePropertyStatusWhere( "element vector", vector<VARIABLE_FLAG>{DIRICH,PLAIN}, -1e30, 1e30 );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^
       VARIABLE_FLAG component0_flag = model_domain.PropertyStatus( "element vector", COMPLETE, 0 );
       VARIABLE_FLAG component1_flag = model_domain.PropertyStatus( "element vector", COMPLETE, 1 );
       _test( component0_flag == DIRICH );
       _test( component1_flag == PLAIN );

       // for vector, tensor and flagged Array variables
       min_value_to_change = -1.0e+30;
       max_value_to_change =  1.0e+30;
       uint32_t component  = 1;
       
       model_domain.InputPropertyValue("nodal variable", makeScalar(ANY,1.), COMPLETE );
       model_domain.InputPropertyValue("nodal vector", makeVector(ANY,ANY,1.,2.), COMPLETE );
       model_domain.InputPropertyValue("element tensor", makeTensor(ANY,ANY,1.,0.,0.,1.), COMPLETE );
       //FlaggedArrayVariable flagged_array = { ANY,PLAIN,DIRICH,ROBIN,FIELD_DATA, 0.5,1.0,1.5,2.0,2.5,3.0 };
       FlaggedArrayVariable farray(5); farray = 0.;
       farray.Flag(0)=ANY; farray.Flag(1)=PLAIN; farray.Flag(2)=DIRICH; farray.Flag(3)=ROBIN; farray.Flag(4)=FIELD_DATA;
       model_domain.InputPropertyValue("nodal flagged array variable", farray, COMPLETE );
       _test( printRangeOfVariable( model, "Model", "nodal flagged array variable" ) == 0. );
       
       // scalar variable (node)
       //                                 var                                        flag
       model_domain.ChangePropertyStatus( "nodal variable", component, ROBIN, PERIMETER );
       //           ^^^^^^^^^^^^^^^^^^^^
       _test ( model_domain.PropertyStatus( "nodal variable", PERIMETER, component ) == ROBIN );
       _test ( model_domain.PropertyStatus( "nodal variable", INTERIOR, component ) != ROBIN );

       model_domain.ChangePropertyStatusWhere( "nodal variable", component, ROBIN, min_value_to_change, max_value_to_change );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^
       _test ( model_domain.PropertyStatus( "nodal variable", COMPLETE, component ) == ROBIN );

       // vector variable (node)
       upper.ChangePropertyStatus( "nodal vector", component, ROBIN, PERIMETER );
       _test ( upper.PropertyStatus( "nodal vector", PERIMETER, component ) == ROBIN );
       _test ( upper.PropertyStatus( "nodal vector", INTERIOR, component ) != ROBIN );

       upper.ChangePropertyStatusWhere( "nodal vector", component, ROBIN, min_value_to_change, max_value_to_change );
       _test ( upper.PropertyStatus( "nodal vector", COMPLETE, component ) == ROBIN );

       // tensor variable (element)
       component = 3;
       upper.ChangePropertyStatus( "element tensor", component, ROBIN, PERIMETER );
       _test ( upper.PropertyStatus( "element tensor", PERIMETER, component ) == ROBIN );
       _test ( upper.PropertyStatus( "element tensor", INTERIOR, component ) != ROBIN );

       upper.ChangePropertyStatusWhere( "element tensor", component, ROBIN, min_value_to_change, max_value_to_change );
       _test ( upper.PropertyStatus( "element tensor", COMPLETE, component ) == ROBIN );

       // flaged array variable (node)
       component = 4;
       upper.ChangePropertyStatus( "nodal flagged array variable", component, ROBIN, PERIMETER );
       _test ( upper.PropertyStatus( "nodal flagged array variable", PERIMETER, component ) == ROBIN );
       _test ( upper.PropertyStatus( "nodal flagged array variable", INTERIOR, component ) != ROBIN );

       upper.ChangePropertyStatusWhere( "nodal flagged array variable", component, ROBIN, min_value_to_change, max_value_to_change );
       _test ( upper.PropertyStatus( "nodal flagged array variable", COMPLETE, component ) == ROBIN );

    } // end 2D
   
} // end Test_PropertyManipulations





/**
    double Average( const char* property ) const;
    Point<dim> AverageUnitNormal() const;
    void InterpolateNodeToCellProperty( const char* nprop, const char* eprop );
    void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* eprop );
    void InterpolateIntegrationPointToCellProperty( const char* cprop, const char* eprop );
    void ExtrapolateCellToIntegrationPointProperty( const char* eprop, const char* cprop );
    void ExtrapolateCellToFacetIntegrationPointProperty( const char* eprop, const char* fipprop );
    void ExtrapolateCellToNodeProperty( const char* eprop, const char* nprop, bool by_distance=true );
    void ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop );
    bool   CopyGradientOfProperty_A_To_B( const char* node_prop, const char* cell_prop );
*/
void ModelSubDomain_Test::Test_PropertyTransfer()
 {
    // 2D Testing: model with everything: regions, boundaries and split boundaries
     {
       VSet<2>        vset;
       ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
       // remove property
       vset.RemoveData("permeability");
       Model<2>       model( topo, vset, "CSMP-variables.txt", false /* don't regions file */ );
       model.Name( "SPLIT22_BASIC" );

       // accessors
       Region<2>& model_domain = model.Region("Model"); // domain idx = 1 (range 1..n!)
       Region<2>& lower = model.Region("lower"); // domain idx = 2
       Region<2>& upper = model.Region("upper"); // domain idx = 3
      
       Boundary<2>& top = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)
       
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2

 
       // TESTING
       
       lower.InputPropertyValue( "nodal variable", makeScalar(ANY,3.) );
       _test( approximatelyEqual(lower.Average( "nodal variable" ), 3. ) == true );
       //                              ^^^^^^^
       upper.InputPropertyValue( "element variable", makeScalar(ANY,0.15) );
       _test( approximatelyEqual(upper.Average( "element variable" ), 0.15 ) == true );
      
       // split boundary
       Point<2> unrml = inclined_split_boundary.AverageUnitNormal();
       //                                       ^^^^^^^^^^^^^^^^^
       _test( approximatelyEqual(0.91487779690726068,unrml[0]) == true );
       _test( approximatelyEqual(0.39617101327914628,unrml[1]) == true );
       // boundary
       unrml = top.AverageUnitNormal();
       _test( approximatelyEqual(unrml[0],0.) == true );
       _test( approximatelyEqual(unrml[1],1.) == true );
       
       // region inserted into split-boundary
       pair<string,bool> dimM1 = model.InsertRegionIntoSplitBoundary( inclined_split_boundary.Name().c_str(), 3 /* material_id_for_new_elements */ );
       _test( dimM1.second == true );
       _test( model.ContainsRegion(dimM1.first) == true );
       const Region<2> dimM1_region = model.Region( dimM1.first );
       unrml = dimM1_region.AverageUnitNormal();
       _test( approximatelyEqual(0.91487779690726068,unrml[0]) == true );
       _test( approximatelyEqual(0.39617101327914628,unrml[1]) == true );
       
       // interpolation
       const char y_coordinate{'y'};
       assignNodeCoordinatesTo( model, y_coordinate, "nodal variable" );
       model_domain.InterpolateNodeToCellProperty( "nodal variable", "element variable" );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       double min_val, max_val;
       model_domain.MinMaxOf( "element variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );

       model.CreateProperty( "eipoint variable", "EIP", "none", SCALAR, ELEMENT_INTEGRATION_POINT );
       model_domain.InterpolateNodeToIntegrationPointProperty( "nodal variable", "eipoint variable" );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       model_domain.MinMaxOf( "eipoint variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );

       model_domain.InterpolateIntegrationPointToCellProperty( "eipoint variable", "element variable" );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       model_domain.MinMaxOf( "element variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );
       
        // extrapolation
        model_domain.ExtrapolateCellToIntegrationPointProperty( "element variable", "eipoint variable" );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       model_domain.MinMaxOf( "element variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );
       
       model.Mesh().InitializeFiniteVolumeStencils( model.Database(), true );
       model.CreateProperty( "fip variable", "FIP", "none", SCALAR, FACET_INTEGRATION_POINT );

       model_domain.ExtrapolateCellToFacetIntegrationPointProperty( "element variable", "fip variable" );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       model_domain.MinMaxOf( "fip variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );

       bool by_distance{ true };
       model_domain.ExtrapolateCellToNodeProperty( "element variable", "nodal variable", by_distance );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       model_domain.MinMaxOf( "nodal variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );

       model_domain.ExtrapolateCellToNodeProperty( "element variable", "nodal variable", by_distance=false );
       model_domain.MinMaxOf( "nodal variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );

       model_domain.InterpolateNodeToIntegrationPointProperty( "nodal variable", "eipoint variable" );
       model_domain.ExtrapolateIntegrationPointToNodeProperty( "eipoint variable", "nodal variable" );
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       model_domain.MinMaxOf( "nodal variable", min_val, max_val );
       _test( min_val >= 0. );
       _test( max_val <= 7. );

       // start from scratch
       assignNodeCoordinatesTo( model, y_coordinate, "nodal variable" );
       bool success = model_domain.CopyGradientOfProperty_A_To_B( "nodal variable", "element vector" );
       //                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( success == true );
        model_domain.MinMaxOf( "element vector", min_val, max_val );
       // _test( approximatelyEqual(min_val,1.) == true ); // only 0.9
       _test( approximatelyEqual(max_val,1.,1.0e-14) == true );

    } // end 2D
    
 } // end Test_PropertyTransfer




/**
      PLACEMENT modelSubdomainType( const std::string& subdomain_name ) noexcept;
      size_t  sharedNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );
      size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );
      size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>&, std::vector<Node<dim>*>& );
      size_t  sharedPerimeterCells( const ModelSubDomain<dim,CELL>& subdomain1, const ModelSubDomain<dim,CELL>& subdomain2,
                                    std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& matching_cells );
                                    
      std::set<TOPOTYPE> nodeTopologyFlags( typename std::vector<Node<dim>*>::const_iterator first,
                                            typename std::vector<Node<dim>*>::const_iterator last );
                                    
      std::set<TOPOTYPE> nodeTopologyFlags( const SplitBoundary<dim>& );

      void readDomainIndexesFromBinaryFile( std::fstream&, SubDomainInfo& ); // tested elsewhere
*/
void ModelSubDomain_Test::Test_NonMemberFunctions()
 {
 
 } // end Test_NonMemberFunctions



} // end csmp
