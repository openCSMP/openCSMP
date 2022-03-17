#include "VTU_Interface_Test.h"

#include "VTU_Interface.h"
#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"
#include "Region.h"
#include "Boundary.h"

using namespace std;

namespace csmp {

VTU_Interface_Test::VTU_Interface_Test()
{
}

void VTU_Interface_Test::run()
{
  enum{DIM=3U};
  vector<double> vector1; vector1.push_back(4.0); vector1.push_back(0.0); vector1.push_back(0.0);
  vector<double> vector2; vector2.push_back(0.0); vector2.push_back(0.5); vector2.push_back(0.0);
  vector<double> vector3; vector3.push_back(0.0); vector3.push_back(0.); vector3.push_back(2.0);
  vector<vector<double> > vectors;

  
  vectors.push_back( vector1 ); vectors.push_back( vector2 ); vectors.push_back( vector3 );
  _test( VTU_Interface<DIM>::OutputVectorsToVTU( "3D_arrowsTest", "effective permeability", vectors ) );

  vector<double> xyzLengths( 3, 0. );
  xyzLengths.at( 0 ) = 1.0; xyzLengths.at( 1 ) = 5.0; xyzLengths.at( 2 ) = 10.0;
  _test( VTU_Interface<DIM>::OutputPrincipalVectorsToVTU( "3D_xzyTest", "effective permeability", xyzLengths ) );
  

  vector<vector<double> > tensor;
  tensor.push_back( vector1 ); tensor.push_back( vector2 ); tensor.push_back( vector3 );
  _test( VTU_Interface<DIM>::OutputTensorToVTU( "3D_tensorTest", "permeability tensor", tensor ) );

  
  ANSYS_Model3D model( "Cube", "CSMP-variables.txt" );
  VTU_Interface<DIM> vtu( model );
  Region<DIM>& rref( model.Region( "Model" ) );
  rref.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0.5 ) );
  rref.InputPropertyValue( "element variable", makeScalar( PLAIN, 1.0E-12 ) );
  ArrayVariable av_node(5);
  ArrayVariable av_elmnt(22);
  FlaggedArrayVariable fav_node(5);
  FlaggedArrayVariable fav_elmnt(22);
  for(auto i =0;i<5;i++)
  {
      av_node(i) = i;
      fav_node(i) = 2.*i;
  }
  for(auto i =0;i<22;i++)
  {
      av_elmnt(i) = 3.*i;
      fav_elmnt(i) = 4.*i;
  }
  rref.InputPropertyValue( "nodal array variable", av_node );
  rref.InputPropertyValue( "element array variable", av_elmnt );
  rref.InputPropertyValue( "nodal flagged array variable", fav_node );
  rref.InputPropertyValue( "element flagged array variable", fav_elmnt );


  Region<DIM>& rref1( model.Region( "LAYER1" ) );
  Region<DIM>& rref2( model.Region( "LAYER2" ) );
  rref1.InputPropertyValue( "element variable", makeScalar( PLAIN, 1.0E-13 ) );
  rref1.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0.1234567890123456 ) ); // to check loss of double precision

  // regions
  _test( vtu.OmitZeroInFileName() == true );
  vtu.OmitZeroInFileName( false );
  _test( vtu.OmitZeroInFileName() == false );
  _test( vtu.OutputDataToVTU( "VTU_RegionName_Test1", "nodal variable", "Model", static_cast<int>(0) ) );
  _test( vtu.OutputDataToVTU( "VTU_RegionReference_Test1", "nodal variable", rref, static_cast<int>(0) ) );
  _test( vtu.OutputDataToVTU( "VTU_RegionName_Test2", "nodal variable", "LAYER2", static_cast<int>(0) ) );
  _test( vtu.OutputDataToVTU( "VTU_RegionReference_Test2", "nodal variable", rref2, static_cast<int>(0) ) );

  // boundary
  Boundary<DIM>& bref( model.Boundary("FRONT") );
  _test( vtu.OutputDataToVTU( "VTU_BoundaryReference_Test1", "nodal variable", bref, static_cast<int>(0) ) );

  // vector variable
  const VectorVariable<3> unitVectorX( PLAIN, PLAIN, PLAIN, 1.0, 0., 0. );
  const VectorVariable<3> vectorA( PLAIN, PLAIN, PLAIN, -1.0, 1.0, 0. );
  list<string> outputPropsVector;
  outputPropsVector.push_back("nodal vector");
  outputPropsVector.push_back("element vector");
  outputPropsVector.push_back("element variable");
  outputPropsVector.push_back("nodal variable");
  list<string> outputPropsScalar;
  outputPropsScalar.push_back("nodal variable");
  outputPropsScalar.push_back("element variable");
  outputPropsScalar.push_back("nodal array variable");
  outputPropsScalar.push_back("element array variable");
  outputPropsScalar.push_back("nodal flagged array variable");
  outputPropsScalar.push_back("element flagged array variable");
  outputPropsScalar.push_back("nodal variable");// check that duplicated names doesn't affect vtu output
  list<string> outputPropsArrays;
  outputPropsArrays.push_back("nodal array variable");
  outputPropsArrays.push_back("element array variable");
  outputPropsArrays.push_back("nodal flagged array variable");
  outputPropsArrays.push_back("element flagged array variable");

  model.InputPropertyValue( "nodal vector", unitVectorX );
  model.InputPropertyValue( "element vector", vectorA );
  vtu.OutputDataToVTU( "VTU_TestVector3D", outputPropsVector, "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "VTU_TestScalar3D", outputPropsScalar, "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "VTU_TestArray3D",  outputPropsArrays, "Model", static_cast<int>(0) );

  // planar model in 3D space (do not create boundaries, else model will fail)
  ANSYS_Model3D modelB( "BoxHalfs2D", "CSMP-variables.txt", true, true, true ); //  TODO: test does not require boundaries
  VTU_Interface<DIM> vtuB( modelB );
  vtuB.OmitZeroInFileName( true );
  modelB.InputPropertyValue( "nodal vector", vectorA );
  modelB.InputPropertyValue( "element vector", vectorA );
  vtuB.OutputDataToVTU( "VTU_TestElementVector2Din3D", outputPropsVector, "Model", static_cast<int>(0) );



  VectorVariable<2> vectorC( PLAIN, PLAIN, 1.0, 1.0 );
  ANSYS_Model2D modelC( "SplitBoundaries2D", "CSMP-variables.txt" );
  VTU_Interface<2> vtuC( modelC );
  vtuC.OmitZeroInFileName( true );
  modelC.InputPropertyValue( "nodal vector", vectorC );
  modelC.InputPropertyValue( "element vector", vectorC );
  modelC.InputPropertyValue( "element variable", makeScalar( PLAIN, 0.7 ) );
  modelC.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0.7 ) );
  vtuC.OutputDataToVTU( "VTU_Test2D", outputPropsVector, "Model", static_cast<int>(0) );


}

} // csmp
