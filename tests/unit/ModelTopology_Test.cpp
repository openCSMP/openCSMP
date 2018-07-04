#include "ModelTopology_Test.h"
#include "ANSYS_ElementSpecifications.h"
#include "PL_Utilities.h"

using namespace std;

namespace csmp{

void ModelTopology_Test::run()
{
  const bool verbose(false);

  setName( "csmp::ModelTopology_Test" );
  if ( verbose ) cout << "\nUnit Test " << getName() << endl;
  typedef ANSYS_ElementSpecifications fem_specs;
  const bool isoparametric( true );
  const size_t dim( 3 );

  // .)CONSTRUCTORS
  ModelTopology topology1; // isoparametric = false
  ModelTopology topology2( false );
  ModelTopology topology3( "topology3", false );

  // .)RETURN FUNCTIONS 1
  string modelName( "topology1" );
  topology2.ModelName( modelName.c_str() );
  _test( topology2.ModelName() == modelName );

  // .)ADDREGION
  set<string> femTypes;
  vector<size_t>   elmIDS;
  const size_t elms( 8 );
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
  _test( topology2.AddRegion( "Region1", femTypes, elmIDS ) );
  _test( !topology2.AddRegion( "Region1", femTypes, elmIDS ) );
  _test( topology2.Elements() == elms );
  set<size_t> elmSet;
  topology2.Elements( elmSet );
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
  _test( topology2.AddRegion( "Region2", femTypes, elmIDS ) );

  // .)RETURN FUNCTIONS 2
  _test( topology2.Elements() == elms2+elms+elms );
  _test( topology2.ElementsOfRegion( "Region2" ) == elms2+elms );
  _test( topology2.ElementsOfRegion( "Region1" ) == elms );
  vector<size_t>::const_iterator region1begin( topology2.ElementsOfRegionBegin( "Region1" ) );
//  vector<size_t>::const_iterator region2end( topology2.ElementsOfRegionEnd( "Region2" ) );
  _test( *region1begin == elmIDS[0] );
  set<string> femTypesReturn;
  topology2.ElementTypesOfRegion( "Region2", femTypesReturn );
  _test( femTypesReturn == femTypes );
  _test( topology2.IsWithinRegion( "Region2", 12 ) );
  _test( !topology2.IsWithinRegion( "Region1", 12 ) );
  _test( topology2.Contains( "Region2" ) );
  _test( !topology2.Contains( "Region3" ) );

  // .)COPY&ASSIGNMENT CONSTRUCTORS
  ModelTopology topology4( topology2 );
  ModelTopology topology5 = topology2;
  _test( topology4.Elements() == topology2.Elements() );
  _test( topology5.Elements() == topology2.Elements() );
  _test( topology4.IsWithinRegion( "Region2", 12 ) );
  _test( !topology4.IsWithinRegion( "Region1", 12 ) );
  _test( topology4.Contains( "Region2" ) );
  _test( !topology4.Contains( "Region3" ) );
  _test( topology5.IsWithinRegion( "Region2", 12 ) );
  _test( !topology5.IsWithinRegion( "Region1", 12 ) );
  _test( topology5.Contains( "Region2" ) );
  _test( !topology5.Contains( "Region3" ) );
  _test( !topology4.AddRegion( "Region1", femTypes, elmIDS ) );
  _test( !topology5.AddRegion( "Region1", femTypes, elmIDS ) );

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
  femTypes.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "POLYGON", isoparametric, dim ) );
  for( size_t i = 0; i < 23; ++i )
    elmIDS.push_back( i );
  ModelTopology topology6( topology2 );
  _test( topology6.AddRegion( "Region3", femTypes, elmIDS ) );


  // .)ELEMENT TYPES
  ModelTopology topTypes( topology6 );
  string oldType( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_27", isoparametric, dim ) );
  topTypes.ChangeElementType( oldType, fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_QUADRATIC", isoparametric, dim ) );
  set<string> region3Types;
  topTypes.ElementTypesOfRegion( "Region3", region3Types );
  for( set<string>::const_iterator it = region3Types.begin(); it != region3Types.end(); ++it )
    _test( *it != fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName("HEXA_27", isoparametric, dim ) );
  // ISSUES, check with stephan
  topology6.EliminateLineElements();
  set<string> top6types;
  topology6.ElementTypesOfRegion( "Region3", top6types );
  // ISSUES, check with stephan
  /*
  _test( top6types.find( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName("BAR_2",isoparametric,dim) ) == top6types.end() );
  _test( top6types.find( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName("BAR_3",isoparametric,dim) ) == top6types.end() );
  _test( top6types.find( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName("POLYGON",isoparametric,dim) ) == top6types.end() );
  */
  if ( verbose ) cout << "\n" << getName() << ": ANSYS FEM Types(stored in ModelTopology:\n";
  for( set<string>::const_iterator it = top6types.begin(); it != top6types.end(); ++it )
    if ( verbose ) cout << *it << endl;
  set<string> csmpTypes;
  topology6.FiniteElementTypes( csmpTypes );
  if ( verbose ) cout << "\n" << getName() << ": CSMP FEM Types(converted from ModelTopology):\n";
  for( set<string>::const_iterator it = csmpTypes.begin(); it != csmpTypes.end(); ++it )
    if ( verbose ) cout << *it << endl;
  set<int32> csmpTypesENUM;
  topology6.FiniteElementTypes( csmpTypesENUM );
  if ( verbose ) cout << "\n" << getName() << ": CSMP FEM Types(converted from ModelTopology):\n";
  for( set<int32>::const_iterator it = csmpTypesENUM.begin(); it != csmpTypesENUM.end(); ++it )
    if ( verbose ) cout << *it << endl;
  set<string> typesCheck;
  _test( 22 == topTypes.FiniteElementTypes( typesCheck ) );
  oldType = fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_QUADRATIC", isoparametric, dim );
  topTypes.ChangeElementType( oldType, fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "HEXA_27", isoparametric, dim ) );
  _test( 23 == topTypes.FiniteElementTypes( typesCheck ) );

  // .)REGION OPS
  ModelTopology topology7( topology2 );
  topology7.ReduceToRegions( "ModelTopology_Test" );
  _test( !topology7.Contains( "Region1" ) );
  list<string> exportRegions;
  exportRegions.push_back( "Region1" );
  topology2.ExportSelectionTo( exportRegions, topology7 );
  _test( topology7.Contains( "Region1" ) );
  set<string> top7typesR1;
  set<string> top2typesR1;
  topology7.ElementTypesOfRegion( "Region1", top7typesR1 );
  topology2.ElementTypesOfRegion( "Region1", top2typesR1 );
  _test( top7typesR1 == top2typesR1 );


  // .) ELMT
  _test( topology7.IsWithinRegion( "Region1", 2 ) );
  _test( !topology7.IsWithinRegion( "Region1", 14 ) );
  _test( topology7.CheckElementNumbering() );
  set<string> femTypes2;
  femTypes2.insert( fem_specs::CSMP_TypeNameFrom_ANSYS_TypeName( "TRI_3", isoparametric, dim ) );
  vector<size_t> elmtIDS2;
  elmtIDS2.push_back( 22 );
  elmtIDS2.push_back( 42 );
  topology7.AddRegion( "Region4", femTypes2, elmtIDS2 );
  _test( !topology7.CheckElementNumbering() );
  map<size_t,size_t> oldNewIDs;
  oldNewIDs.insert( make_pair( 22, 15 ) );
  oldNewIDs.insert( make_pair( 42, 16 ) );
  topology7.CreateNewElementNumbers( oldNewIDs );
  _test( topology7.CheckElementNumbering() );

  // .)TOPOLOGY TYPE
  _test( !topology7.BoxShapedModel() );
  _test( !topology7.RectangleShapedModel() );
  topology7.AddRegion( "LEFT", femTypes2, elmtIDS2 );
  topology7.AddRegion( "RIGHT", femTypes2, elmtIDS2 );
  topology7.AddRegion( "BOTTOM", femTypes2, elmtIDS2 );
  topology7.AddRegion( "TOP", femTypes2, elmtIDS2 );
  topology7.AddRegion( "FRONT", femTypes2, elmtIDS2 );
  topology7.AddRegion( "BACK", femTypes2, elmtIDS2 );
  _test( topology7.BoxShapedModel() );
  _test( topology7.RectangleShapedModel() );
  _test( topology7.SolidModel() );
  _test( !topology7.LineModel() );
  _test( !topology7.SurfaceModel() );
  _test( topology7.MinimumSpatialDimensionOfRegion( "Region1" ) == 0U );
  _test( topology7.IsoparametricElements() == false );
  _test( topology7.QuadraticElementMesh() == false );
  _test( topology7.LinearElementMesh() == true );
  topology7.TreatElementsAsIsoparametric();
  _test( topology7.IsoparametricElements() == true );


} // run

/*
PropertiesOfRegions() not tested

Container class that stores the region names and their contained elements(through IDs)
and corresponding fem types.

ModelTopology Does not check for proper input!!! (ie if vector and set are of same size...)

*/

} // csmp

