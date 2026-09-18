#include "PointPropertyToCellMapper2D_Test.h"
#include "Region.h"
#include "Boundary.h"

#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "DenseMatrix.h"
#include "plf_colony.h"

using namespace std;

namespace csmp
{


void PointPropertyToCellMapper2D_Test::run()
  {
    const bool verbose(false);
    
    ANSYS_Model2D model( "BoxHalfs2D", "CSMP-variables.txt" );
    Region<2>& rref( model.Region( "Model" ) );


    DenseMatrix<DM_MIN> dm;
    const auto elementsEnd( rref.CellsEnd() );
    try
      {    
        for( auto it = rref.CellsBegin(); it != elementsEnd; ++it )
          (*it)->CoordinateMatrix();
      }
    catch(...)
      {
        _test(false);
      }

    if ( verbose ) cout << "\nModel Node Count: " << rref.Nodes() << endl;
    _test( rref.Nodes() == 106 );

    Index nodalKey( model.Database().StorageKey( "nodal variable" ) );
    rref.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0. ) );
    size_t boundaryNodeCount( 0 );
    const auto nodesEnd( rref.NodesEnd() );
    for( auto it = rref.NodesBegin(); it != nodesEnd; ++it )
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
    
    // set 'face variable' to a non-NaN value so that it gets preserved when output to binary
    model.InputPropertyValue( "face variable", makeScalar( PLAIN, 0. ) );
    right.InputPropertyValue( "face variable", makeScalar( PLAIN, 2. ) );
    
    ArrayVariable av( "nodal array", model.Database() );
    // SKM FIX: if av is not initialized its values are NAN and a comparison with another array variable will always evaluate as false
    av = 5.;
    model.InputPropertyValue( "nodal array", av );
      
    const size_t leftNodes( left.Nodes() );
    const size_t rightNodes( right.Nodes() );
    const size_t bottomNodes( bottom.Nodes() );
    const size_t topNodes( top.Nodes() );
    const size_t leftFaces( left.Cells() );
    const size_t rightFaces( right.Cells() );
    const size_t bottomFaces( bottom.Cells() );
    const size_t topFaces( top.Cells() );
    
    if ( verbose ) cout << "\nModel Node Count: " << rref.Nodes() << endl;
    _test( rref.Nodes() == 106 );

    string bin1name("ANSYS2D_bin");
    model.OutputToBinaryFile(bin1name.c_str());
    // model.Out();

    // Tests on the model recreated from disk
    Model<2U> modelBinIn0(bin1name);
    
    Index nodalArrayKey0 = modelBinIn0.Database().StorageKey("nodal array");
    Index faceVariableKey0 = modelBinIn0.Database().StorageKey("face variable");

    _test( (*modelBinIn0.Boundary("RIGHT").CellsBegin())->Read(faceVariableKey0) == 2. );

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
    _test( leftFaces == modelBinIn1.Boundary("LEFT").Cells() );
    _test( rightFaces == modelBinIn1.Boundary("RIGHT").Cells() );
    _test( topFaces == modelBinIn1.Boundary("TOP").Cells() );
    _test( bottomFaces == modelBinIn1.Boundary("BOTTOM").Cells() );

    _test( (*modelBinIn1.Boundary("RIGHT").CellsBegin())->Read(faceVariableKey) == 2. );
    ArrayVariable avBin( "nodal array", modelBinIn1.Database() );
    (*modelBinIn1.Region("Model").NodesBegin())->Read( nodalArrayKey, avBin );
    _test( avBin == av );
  }

} // csmp
