#include "ANSYS_Model2D_Test.h"
#include "Region.h"
#include "Boundary.h"

#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "DenseMatrix.h"
#include "plf_colony.h"

using namespace std;

namespace csmp
{


void ANSYS_Model2D_Test::run()
  {
    const bool verbose(false);
    
    ANSYS_Model2D model( "BoxHalfs2D", "CSMP-variables.txt" );
    Region<2>& rref( model.Region( "Model" ) );


    DenseMatrix<DM_MIN> dm;
    const vector<Element<2>*>::const_iterator elementsEnd( rref.ElementsEnd() );
    try
      {    
        for( vector<Element<2>*>::const_iterator it = rref.ElementsBegin(); it != elementsEnd; ++it )
          (*it)->CoordinateMatrix();
      }
    catch(...)
      {
        _test(false);
      }

    if ( verbose ) cout << "\nModel Node Count: " << rref.Nodes() << endl;
    _test( rref.Nodes() == 106 );

    const Index nodalKey( model.Database().StorageKey( "nodal variable" ) );
    rref.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0. ) );
    size_t boundaryNodeCount( 0 );
    const vector<Node<2>*>::const_iterator nodesEnd( rref.NodesEnd() );
    for( vector<Node<2>*>::const_iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
      {
        const BOX_BOUNDARY boxBoundary( (*it)->AtBoundary() );
        if( boxBoundary != NOT )
          (*it)->Store( nodalKey, makeScalar( PLAIN, 1.0 ) );
        else
          continue;   
        ++boundaryNodeCount;
      }
    if ( verbose ) cout << "\nBoundary Node Count: " << boundaryNodeCount << endl;
    _test ( boundaryNodeCount == 40 );
    if ( verbose ) {
         VTU_Interface<2> vtu( model );
         vtu.OmitZeroInFileName( true );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-BoxBoundary", "nodal variable", "Model", static_cast<int>(0) );
      }
    _test( model.ContainsBoundary( "LEFT" ) );
    _test( model.ContainsBoundary( "RIGHT" ) );
    _test( model.ContainsBoundary( "BOTTOM" ) );
    _test( model.ContainsBoundary( "TOP" ) );

    _test( model.ContainsRegion( "MATRIX_LEFT" ) );
    _test( model.ContainsRegion( "MATRIX_RIGHT" ) );


    Boundary<2>& left( model.Boundary( std::string("LEFT") ) );
    Boundary<2>& right( model.Boundary( std::string("RIGHT") ) );
    Boundary<2>& bottom( model.Boundary( std::string("BOTTOM") ) );
    Boundary<2>& top( model.Boundary( std::string("TOP") ) );
    
    // testing that the corner nodes have been flagged correctly
    const Region<2>& mref(model.Region("Model"));
    size_t  counter(0U);
    for ( auto it=mref.PerimeterNodesBegin(); it!=mref.NodesEnd(); it++ ) {
         if ( (*it)->AtBoundary() == CNR1 ) counter++;
         else if ( (*it)->AtBoundary() == CNR2 ) counter++;
         else if ( (*it)->AtBoundary() == CNR3 ) counter++;
         else if ( (*it)->AtBoundary() == CNR4 ) counter++;
      }
    _test( counter == 4U );

    right.InputPropertyValue( "face variable", makeScalar( PLAIN, 2. ) );
    ArrayVariable av( "nodal array", model.Database() );
    // SKM FIX: if av is not initialized its values are NAN and a comparison with another array variable will always evaluate as false
    av = 5.;
    model.InputPropertyValue( "nodal array", av );
      
    const size_t leftNodes( left.Nodes() );
    const size_t rightNodes( right.Nodes() );
    const size_t bottomNodes( bottom.Nodes() );
    const size_t topNodes( top.Nodes() );
    const size_t leftFaces( left.Elements() );
    const size_t rightFaces( right.Elements() );
    const size_t bottomFaces( bottom.Elements() );
    const size_t topFaces( top.Elements() );
    
    if ( verbose ) cout << "\nModel Node Count: " << rref.Nodes() << endl;
    _test( rref.Nodes() == 106 );

    string bin1name("ANSYS2D_bin");
    model.OutputToBinaryFile(bin1name.c_str());
    // model.Out();

    Model<2U> modelBinIn0(bin1name);
    Index nodalArrayKey0( modelBinIn0.Database().StorageKey("nodal array") );
    Index faceVariableKey0( modelBinIn0.Database().StorageKey("face variable") );
    _test( (*modelBinIn0.Boundary("RIGHT").ElementsBegin())->Read(faceVariableKey0) == 2. );
    ArrayVariable avBin0( "nodal array", modelBinIn0.Database() );
    (*modelBinIn0.Region("Model").NodesBegin())->Read( nodalArrayKey0, avBin0 );
    _test( avBin0 == av );
    
    Model<2U> modelBinIn1(bin1name);
    Index nodalArrayKey( modelBinIn1.Database().StorageKey("nodal array") );
    Index faceVariableKey( modelBinIn1.Database().StorageKey("face variable") );
    _test( leftNodes == modelBinIn1.Boundary("LEFT").Nodes() );
    _test( rightNodes == modelBinIn1.Boundary("RIGHT").Nodes() );
    _test( topNodes == modelBinIn1.Boundary("TOP").Nodes() );
    _test( bottomNodes == modelBinIn1.Boundary("BOTTOM").Nodes() );
    _test( leftFaces == modelBinIn1.Boundary("LEFT").Elements() );
    _test( rightFaces == modelBinIn1.Boundary("RIGHT").Elements() );
    _test( topFaces == modelBinIn1.Boundary("TOP").Elements() );
    _test( bottomFaces == modelBinIn1.Boundary("BOTTOM").Elements() );

    _test( (*modelBinIn1.Boundary("RIGHT").ElementsBegin())->Read(faceVariableKey) == 2. );
    ArrayVariable avBin( "nodal array", modelBinIn1.Database() );
    (*modelBinIn1.Region("Model").NodesBegin())->Read( nodalArrayKey, avBin );
    _test( avBin == av );
  }

} // csmp
