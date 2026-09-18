#include "ModelTopology_Test.h"
#include "ANSYS_ElementSpecifications.h"
#include "PL_Utilities.h"
#include "vsetMakers.h"

using namespace std;

namespace csmp{

void ModelTopology_Test::run()
{
  setName( "csmp::ModelTopology_Test" );
  if ( verbose_ ) cout << "\nUnit Test " << getName() << endl;
  typedef ANSYS_ElementSpecifications fem_specs;
  const bool isoparametric( true );
  const uint32_t dim( 3 );

  // .)CONSTRUCTORS
  ModelTopology topology1; // isoparametric = false
  ModelTopology topology2( false );
  ModelTopology topology3( "topology3", false );

  // .)RETURN FUNCTIONS 1
  string modelName( "topology1" );
  topology2.ModelName( modelName.c_str() );
  _test( topology2.ModelName() == modelName );

  // .)ADDREGION
  set<string>    femTypes;
  vector<size_t> elmIDS;
  const size_t   elms( 8 );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TETRA_4", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TETRA_4", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TETRA_4", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TETRA_4", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_8", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_8", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_8", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_8", isoparametric, dim ) );
  elmIDS.push_back( 0 );
  elmIDS.push_back( 1 );
  elmIDS.push_back( 2 );
  elmIDS.push_back( 3 );
  elmIDS.push_back( 4 );
  elmIDS.push_back( 5 );
  elmIDS.push_back( 6 );
  elmIDS.push_back( 7 );
  _test( topology2.AddDomain( "Region1", femTypes, elmIDS ) );
  _test( !topology2.AddDomain( "Region1", femTypes, elmIDS ) );
  _test( topology2.Cells() == elms );
  set<size_t> elmSet;
  topology2.Cells( elmSet );
  _test( elmSet.size() == elms );
  const size_t elms2( 4 );
  elmIDS.push_back( 11 );
  elmIDS.push_back( 12 );
  elmIDS.push_back( 13 );
  elmIDS.push_back( 14 );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_3", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "BAR_2", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_3_X", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "BAR_3", isoparametric, dim ) );
  _test( topology2.AddDomain( "Region2", femTypes, elmIDS ) );

  // .)RETURN FUNCTIONS 2
  _test( topology2.Cells() == elms2+elms+elms );
  _test( topology2.CellsWithinDomain( "Region2" ) == elms2+elms );
  _test( topology2.CellsWithinDomain( "Region1" ) == elms );
  auto region1begin( topology2.CellsOfDomainBegin( "Region1" ) );
//  vector<uint32_t>::const_iterator region2end( topology2.ElementsOfRegionEnd( "Region2" ) );
  _test( *region1begin == elmIDS[0] );
  set<string> femTypesReturn;
  topology2.CellTypesOfDomain( "Region2", femTypesReturn );
  _test( femTypesReturn == femTypes );
  _test( topology2.IsWithinDomain( "Region2", 12 ) );
  _test( !topology2.IsWithinDomain( "Region1", 12 ) );
  _test( topology2.Contains( "Region2" ) );
  _test( !topology2.Contains( "Region3" ) );

  // .)COPY&ASSIGNMENT CONSTRUCTORS
  ModelTopology topology4( topology2 );
  ModelTopology topology5 = topology2;
  _test( topology4.Cells() == topology2.Cells() );
  _test( topology5.Cells() == topology2.Cells() );
  _test( topology4.IsWithinDomain( "Region2", 12 ) );
  _test( !topology4.IsWithinDomain( "Region1", 12 ) );
  _test( topology4.Contains( "Region2" ) );
  _test( !topology4.Contains( "Region3" ) );
  _test( topology5.IsWithinDomain( "Region2", 12 ) );
  _test( !topology5.IsWithinDomain( "Region1", 12 ) );
  _test( topology5.Contains( "Region2" ) );
  _test( !topology5.Contains( "Region3" ) );
  _test( !topology4.AddDomain( "Region1", femTypes, elmIDS ) );
  _test( !topology5.AddDomain( "Region1", femTypes, elmIDS ) );

  // .)FEM TYPES
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_3", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_3_X", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_6", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_6_X", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TETRA_4", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TETRA_10", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "PYRA_5", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "PYRA_13", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "PYRA_14", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "PENTA_6", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "PENTA_15", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "PENTA_18", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "QUAD_4", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "QUAD_4_X", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "QUAD_8", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "QUAD_8_X", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "QUAD_9", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_8", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_20", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_27", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "BAR_2", isoparametric, dim ) );
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "BAR_3", isoparametric, dim ) );
// POLYGON does not exist:  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "POLYGON", isoparametric, dim ) );
  for( size_t i = 0; i < 23; ++i )
    elmIDS.push_back( i );
  ModelTopology topology6( topology2 );
  _test( topology6.AddDomain( "Region3", femTypes, elmIDS ) );


  // .)ELEMENT TYPES
  ModelTopology topTypes( topology6 );
  string oldType( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_27", isoparametric, dim ) );
  topTypes.ChangeCellType( oldType, fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_20", isoparametric, dim ) );
  set<string> region3Types;
  topTypes.CellTypesOfDomain( "Region3", region3Types );
  for( auto it = region3Types.begin(); it != region3Types.end(); ++it )
    _test( *it != fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName("HEXA_27", isoparametric, dim ) );

  set<string> top6types;
  topology6.CellTypesOfDomain( "Region3", top6types );

  if ( verbose_ ) cout << "\n" << getName() << ": ANSYS FEM Types(stored in ModelTopology:\n";
  for( auto it = top6types.begin(); it != top6types.end(); ++it )
    if ( verbose_ ) cout << *it << endl;
  set<string> csmpTypes;
  topology6.FiniteElementTypes( csmpTypes );
  if ( verbose_ ) cout << "\n" << getName() << ": CSMP FEM Types(converted from ModelTopology):\n";
  for( auto it = csmpTypes.begin(); it != csmpTypes.end(); ++it )
    if ( verbose_ ) cout << *it << endl;
  set<CSMP_FEM_TYPE> csmpTypesENUM;
  topology6.FiniteElementTypes( csmpTypesENUM );
  if ( verbose_ ) cout << "\n" << getName() << ": CSMP FEM Types(converted from ModelTopology):\n";
  for( auto it = csmpTypesENUM.begin(); it != csmpTypesENUM.end(); ++it )
    if ( verbose_ ) cout << parseFiniteElementType(*it) << endl;
  set<string> typesCheck;
  _test( topTypes.FiniteElementTypes( typesCheck ) == 21 );
  oldType = fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_20", isoparametric, dim );
  // adding an extra element type
  topTypes.ChangeCellType( oldType, fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_27", isoparametric, dim ) );
  _test( topTypes.FiniteElementTypes( typesCheck ) == 22 );

  // .)REGION OPS
  ModelTopology topology7( topology2 );
  topology7.ReduceToDomains( "ModelTopology_Test" );
  _test( !topology7.Contains( "Region1" ) );
  list<string> exportRegions;
  exportRegions.push_back( "Region1" );
  topology2.ExportSelectionTo( exportRegions, topology7 );
  _test( topology7.Contains( "Region1" ) );
  set<string> top7typesR1;
  set<string> top2typesR1;
  topology7.CellTypesOfDomain( "Region1", top7typesR1 );
  topology2.CellTypesOfDomain( "Region1", top2typesR1 );
  _test( top7typesR1 == top2typesR1 );


  // .) ELMT
  _test( topology7.IsWithinDomain( "Region1", 2 ) );
  _test( !topology7.IsWithinDomain( "Region1", 14 ) );
  _test( topology7.CheckCellNumbering() );
  set<string> femTypes2;
  femTypes2.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_3", isoparametric, dim ) );
  vector<size_t> elmtIDS2;
  elmtIDS2.push_back( 22 );
  elmtIDS2.push_back( 42 );
  topology7.AddDomain( "Region4", femTypes2, elmtIDS2 );
  _test( !topology7.CheckCellNumbering() );
  map<size_t,size_t> oldNewIDs;
  oldNewIDs.insert( make_pair( 22, 15 ) );
  oldNewIDs.insert( make_pair( 42, 16 ) );
  topology7.CreateNewCellNumbers( oldNewIDs );
  _test( topology7.CheckCellNumbering() );

  // .)TOPOLOGY TYPE
  _test( !topology7.BoxShapedModel() );
  _test( !topology7.RectangleShapedModel() );
  topology7.AddDomain( "LEFT", femTypes2, elmtIDS2 );
  topology7.AddDomain( "RIGHT", femTypes2, elmtIDS2 );
  topology7.AddDomain( "BOTTOM", femTypes2, elmtIDS2 );
  topology7.AddDomain( "TOP", femTypes2, elmtIDS2 );
  topology7.AddDomain( "FRONT", femTypes2, elmtIDS2 );
  topology7.AddDomain( "BACK", femTypes2, elmtIDS2 );
  _test( topology7.BoxShapedModel() );
  _test( topology7.RectangleShapedModel() );
  _test( topology7.SolidModel() );
  _test( !topology7.LineModel() );
  _test( !topology7.SurfaceModel() );
  _test( topology7.MinimumSpatialDimensionOfDomain( "Region1" ) == 3U ); // only volume elements contained
  _test( topology7.IsoparametricFiniteElements() == false );
  _test( topology7.QuadraticElementMesh() == false );
  _test( topology7.LinearElementMesh() == true );
  topology7.UseIsoparametricFiniteElementTypes();
  _test( topology7.IsoparametricFiniteElements() == true );

    // ========================================================================
    // TEST 1: Construction and basic state
    // ========================================================================
    {
        ModelTopology mt;
        _test( mt.ModelDomains() == 0 );
        _test( mt.Cells() == 0 );
        _test( !mt.IsoparametricFiniteElements() );
        _test( mt.ModelName() == "not initialized" );
    }

    // ========================================================================
    // TEST 2: AddDomain / Contains / RemoveDomain
    // ========================================================================
    {
        ModelTopology mt;
        set<string> types = { "TETRA4" };
        vector<size_t> ids = { 0, 1, 2, 3, 4 };

        _test( mt.AddDomain( "rock", types, ids ) );
        _test( mt.Contains( "rock" ) );
        _test( mt.ModelDomains() == 1 );
        _test( mt.Cells() == 5 );
        _test( mt.CellsWithinDomain( "rock" ) == 5 );

        // duplicate domain should fail
        _test( !mt.AddDomain( "rock", types, ids ) );

        mt.RemoveDomain( "rock" );
        _test( !mt.Contains( "rock" ) );
        _test( mt.ModelDomains() == 0 );
    }

    // ========================================================================
    // TEST 3: AddDomainCellId / AddDomainCellIds / AddDomainCellType
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "fractures", "TRIA3" );
        mt.AddDomainCellId( "fractures", 10 );
        mt.AddDomainCellId( "fractures", 11 );
        mt.AddDomainCellIds( "fractures", { 12, 13, 14 } );

        _test( mt.Contains( "fractures" ) );
        _test( mt.CellsWithinDomain( "fractures" ) == 5 );
        _test( mt.IsWithinDomain( "fractures", 12 ) );
        _test( !mt.IsWithinDomain( "fractures", 99 ) );
    }

    // ========================================================================
    // TEST 4: ModelName
    // ========================================================================
    {
        ModelTopology mt( "test_model", false );
        _test( mt.ModelName() == "test_model" );
        mt.ModelName( "renamed_model" );
        _test( mt.ModelName() == "renamed_model" );
    }

    // ========================================================================
    // TEST 5: SolidModel / SurfaceModel / LineModel
    // ========================================================================
    {
        ModelTopology mt_solid;
        mt_solid.AddDomainCellType( "rock", "LINEAR_TETRAHEDRON" );
        _test( mt_solid.SolidModel() );
        _test( !mt_solid.SurfaceModel() );
        _test( !mt_solid.LineModel() );

        ModelTopology mt_surface;
        mt_surface.AddDomainCellType( "surface", "LINEAR_TRIANGLE" );
        _test( !mt_surface.SolidModel() );
        _test( mt_surface.SurfaceModel() );
        _test( !mt_surface.LineModel() );

        ModelTopology mt_line;
        mt_line.AddDomainCellType( "well", "LINEAR_BAR" );
        _test( !mt_line.SolidModel() );
        _test( mt_line.SurfaceModel() );  // line is also surface (no volumes)
        _test( mt_line.LineModel() );
    }

    // ========================================================================
    // TEST 6: IsoparametricFiniteElements
    // ========================================================================
    {
        ModelTopology mt_iso( true );
        _test( mt_iso.IsoparametricFiniteElements() );

        ModelTopology mt_std( false );
        _test( !mt_std.IsoparametricFiniteElements() );
        mt_std.UseIsoparametricFiniteElementTypes();
        _test( mt_std.IsoparametricFiniteElements() );
    }

    // ========================================================================
    // TEST 7: EliminateCellType / EliminateLineCells / EliminateSurfaceCells
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "rock",     "LINEAR_TETRAHEDRON" );
        mt.AddDomainCellType( "surface",  "LINEAR_TRIANGLE" );
        mt.AddDomainCellType( "well",     "LINEAR_BAR" );

        mt.EliminateLineCells();
        _test( !mt.Contains( "well" ) );
        _test( mt.Contains( "rock" ) );
        _test( mt.Contains( "surface" ) );

        mt.EliminateSurfaceCells();
        _test( !mt.Contains( "surface" ) );
        _test( mt.Contains( "rock" ) );
    }

    // ========================================================================
    // TEST 8: EliminateVolumeCells
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "rock",    "LINEAR_TETRAHEDRON" );
        mt.AddDomainCellType( "surface", "LINEAR_TRIANGLE" );

        mt.EliminateVolumeCells();
        _test( !mt.Contains( "rock" ) );
        _test( mt.Contains( "surface" ) );
    }

    // ========================================================================
    // TEST 9: OutputAll / OutputRegions / OutputBoundaries
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "rock",              "LINEAR_TETRAHEDRON" );
        mt.AddDomainCellType( "BOUNDARY_BOTTOM",   "LINEAR_TRIANGLE" );
        mt.AddDomainCellType( "SPLITBOUNDARY_FRAC","LINEAR_TRIANGLE" );

        list<string> all;
        mt.OutputAll( all );
        _test( all.size() == 3 );

        list<string> regions;
        mt.OutputRegions( regions );
        _test( regions.size() == 1 );
        _test( regions.front() == "rock" );

        list<string> boundaries;
        mt.OutputBoundaries( boundaries );
        _test( boundaries.size() == 1 );
        _test( boundaries.front() == "BOUNDARY_BOTTOM" );

        list<string> split_boundaries;
        mt.OutputSplitBoundaries( split_boundaries );
        _test( split_boundaries.size() == 1 );
        _test( split_boundaries.front() == "SPLITBOUNDARY_FRAC" );
    }

    // ========================================================================
    // TEST 10: ReduceToDomains
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "rock",     "LINEAR_TETRAHEDRON" );
        mt.AddDomainCellType( "fracture", "LINEAR_TRIANGLE" );
        mt.AddDomainCellType( "well",     "LINEAR_BAR" );

        set<string> keep = { "rock", "fracture" };
        mt.ReduceToDomains( keep );

        _test( mt.Contains( "rock" ) );
        _test( mt.Contains( "fracture" ) );
        _test( !mt.Contains( "well" ) );
        _test( mt.ModelDomains() == 2 );
    }

    // ========================================================================
    // TEST 11: ChangeDomainName
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "old_name", "LINEAR_TETRAHEDRON" );
        mt.AddDomainCellId( "old_name", 0 );

        _test( mt.ChangeDomainName( "old_name", "new_name" ) );
        _test( mt.Contains( "new_name" ) );
        _test( !mt.Contains( "old_name" ) );
        _test( mt.CellsWithinDomain( "new_name" ) == 1 );

        // non-existent domain should return false
        _test( !mt.ChangeDomainName( "does_not_exist", "anything" ) );
    }

    // ========================================================================
    // TEST 12: CheckCellNumbering
    // ========================================================================
    {
        // Consecutive numbering: should pass
        ModelTopology mt_ok;
        mt_ok.AddDomainCellIds( "rock", { 0, 1, 2, 3, 4 } );
        _test( mt_ok.CheckCellNumbering() );

        // Non-consecutive: should fail
        ModelTopology mt_bad;
        mt_bad.AddDomainCellIds( "rock", { 0, 1, 3, 4 } );  // gap at 2
        _test( !mt_bad.CheckCellNumbering() );

        // Not starting at 0: should fail
        ModelTopology mt_bad2;
        mt_bad2.AddDomainCellIds( "rock", { 1, 2, 3, 4 } );
        _test( !mt_bad2.CheckCellNumbering() );
    }

    // ========================================================================
    // TEST 13: CreateNewCellNumbers
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellIds( "rock",     { 0, 2, 4 } );
        mt.AddDomainCellIds( "fracture", { 1, 3, 5 } );

        map<size_t,size_t> mapping;
        mt.CreateNewCellNumbers( mapping );

        // All 6 elements should be mapped to 0..5
        _test( mapping.size() == 6 );
        set<size_t> new_ids;
        for ( const auto& p : mapping ) new_ids.insert( p.second );
        _test( *new_ids.begin() == 0 );
        _test( *new_ids.rbegin() == 5 );
    }

    // ========================================================================
    // TEST 14: ExportSelectionTo
    // ========================================================================
    {
        ModelTopology mt_src;
        mt_src.AddDomainCellType( "rock",     "LINEAR_TETRAHEDRON" );
        mt_src.AddDomainCellType( "fracture", "LINEAR_TRIANGLE" );
        mt_src.AddDomainCellType( "well",     "LINEAR_BAR" );

        list<string> selection = { "rock", "fracture" };
        ModelTopology mt_dst;
        mt_src.ExportSelectionTo( selection, mt_dst );

        _test( mt_dst.Contains( "rock" ) );
        _test( mt_dst.Contains( "fracture" ) );
        _test( !mt_dst.Contains( "well" ) );
        _test( mt_dst.ModelDomains() == 2 );
    }

    // ========================================================================
    // TEST 15: BoxShapedModel / RectangleShapedModel
    // ========================================================================
    {
        ModelTopology mt_box;
        for ( const char* b : { "BOTTOM","TOP","LEFT","RIGHT","FRONT","BACK" } )
            mt_box.AddDomainCellType( b, "LINEAR_TRIANGLE" );
        _test( mt_box.BoxShapedModel() );

        ModelTopology mt_rect;
        for ( const char* b : { "BOTTOM","TOP","LEFT","RIGHT" } )
            mt_rect.AddDomainCellType( b, "LINEAR_BAR" );
        _test( mt_rect.RectangleShapedModel() );

        // Incomplete box: should fail
        ModelTopology mt_incomplete;
        for ( const char* b : { "BOTTOM","TOP","LEFT","RIGHT" } )
            mt_incomplete.AddDomainCellType( b, "LINEAR_TRIANGLE" );
        _test( !mt_incomplete.BoxShapedModel() );
    }

    // ========================================================================
    // TEST 16: Out() and InputFromTextFile() round-trip
    // ========================================================================
    {
        ModelTopology mt_out;
        mt_out.ModelName( "round_trip_test" );
        mt_out.AddDomainCellType( "rock", "LINEAR_TETRAHEDRON" );
        mt_out.AddDomainCellIds( "rock", { 0, 1, 2 } );
        mt_out.AddDomainCellType( "fracture", "LINEAR_TRIANGLE" );
        mt_out.AddDomainCellIds( "fracture", { 3, 4 } );

        const char* test_file = "test_topology.asc";
        mt_out.Out( test_file );

        ModelTopology mt_in;
        _test( mt_in.InputFromTextFile( test_file ) );
        _test( mt_in.Contains( "rock" ) );
        _test( mt_in.Contains( "fracture" ) );
        _test( mt_in.CellsWithinDomain( "rock" ) == 3 );
        _test( mt_in.CellsWithinDomain( "fracture" ) == 2 );
    }

    // ========================================================================
    // TEST 17: Erase
    // ========================================================================
    {
        ModelTopology mt;
        mt.AddDomainCellType( "rock", "LINEAR_TETRAHEDRON" );
        mt.AddDomainCellType( "fracture", "LINEAR_TRIANGLE" );
        _test( mt.ModelDomains() == 2 );

        mt.Erase();
        _test( mt.ModelDomains() == 0 );
        _test( mt.Cells() == 0 );
        _test( mt.ModelName() == "erased" );
    }

} // run



} // csmp

