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
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "vsetMakers.h"
#include "Element.h"
#include "Box.h"

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
     Test_Basics();
     Test_SubDomainConstructionMethods();
  
  
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
       {
       VTU_Interface<2> vtu_output( model );
       list<string>     output_props{"node number"};
       vtu_output.OutputDataToVTU( "split22_basic_", output_props, model_domain, 0 );
       SplitBoundary<2>& horizontal_splitboundary = model.SplitBoundary("horizontal_splitboundary"); // domain idx = 1
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2
       vtu_output.OutputDataToVTU( "split22_basic_", output_props, horizontal_splitboundary, 0 );
       vtu_output.OutputDataToVTU( "split22_basic_", output_props, inclined_split_boundary, 0 );
       }
       
       // testing
       const csmp::Index nkey = model.Database().StorageKey("node number");
       
       // TOPOTYPE flags: are they the same as in the input VSet?
       // -------------------------------------------------------
       vector<TOPOTYPE> gflags0; gflags0.reserve( model_domain.Nodes() );
       if ( verbose_ ) cout <<"\n"<<"classifiers in VData vs. TOPOTYPE attribute in Model:";
       for ( const auto& node : model_domain.NodeVector() ) gflags0.push_back( node->Attribute() );
       // were any errors made during model construction?
       // (NOTE: TOPOTYPE data exactly as in the VSet because VData::InitialiseNodeTopologyIdentifiers() is not called)
       for ( size_t n{0}; n<vset.Vertices(); ++n ) {
           // finding the node with the matching node by number (that was stored in the VSet)
           size_t n_node = numeric_limits<size_t>::max();
           for ( const auto& node : model_domain.NodeVector() )
             if ( static_cast<size_t>(node->Read(nkey)) == n ) {
                  n_node = static_cast<size_t>(node->Read(nkey));
                  break;
             }
           if ( verbose_ ) cout <<"\n\t"<<"Node "<< n <<": "
                                << parseBoundary(static_cast<BOX_BOUNDARY>(vset.BFlag(n))) <<": "
                                << parseTopology( static_cast<TOPOTYPE>(vset.BREP_Flag(n)) )
                                <<":"<< parseTopology( model_domain.N(n_node)->Attribute() );
//           _test( vset.BREP_Flag(n) == model_domain.N(n_node)->Attribute() );
         }
       
       model_domain.UpdateTopoTypeNodeFlags();
       //           ^^^^^^^^^^^^^^^^^^^^^^^^^
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
       
       auto elmts = upper.RebuildCellAndPerimeterFaceVector();
       //                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( elmts == upper.Cells() );
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
       PropertyConstraints constraints( "element number", 13, 16 );
       model_domain.RenumberCells();
       
       lower.UpdateCellMembershipApplyingConstraints( model_domain.CellsBegin(), model_domain.CellsEnd(), constraints );
       //    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       _test( lower.IsUnique() == false );
       _test( lower.Cells() == 4 );
       _test( lower.PerimeterNodes() > 4 );
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
    uint32_t          PerimeterFaces( size_t cell_idx ) const;
    uint32_t          PerimeterFace( size_t cell_idx, uint32_t face ) const;
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
      
       Boundary<2>& top = model.Boundary("TOP"); // domain idx = 4 (last of 4 boundaries created)
       
       SplitBoundary<2>& horizontal_splitboundary = model.SplitBoundary("horizontal_splitboundary"); // domain idx = 1
       SplitBoundary<2>& inclined_split_boundary = model.SplitBoundary("inclined_split_boundary");   // domain idx = 2

       // testing
     }
 } // end


} // end csmp
