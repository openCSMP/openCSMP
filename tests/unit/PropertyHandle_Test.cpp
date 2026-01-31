// PropertyHandle Test
// Created by Alina Yapparova
// 21.11.2011

#include "PropertyHandle_Test.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "PropertyHandle.h"

using namespace std;

namespace csmp {

/**
        TODO: add additional test where PropertyHandle operations are restricted to specific ModelSubDomain objects
 */
PropertyHandle_Test::PropertyHandle_Test( double tolerance, bool verbose )
    : model( new ANSYS_Model3D("BoxHalfs3D", "PropertyHandle_Test-variables.txt") ),
      TOLERANCE(tolerance),
      elementVariable1( *model, "element variable 1", SCALAR, ELEMENT ),
      elementVariable2( *model, "element variable 2", SCALAR, ELEMENT ),
      nodeVariable1( *model, "node variable 1", SCALAR, NODE ),
      nodeVariable2( *model, "node variable 2", SCALAR, NODE ),
      IPVariable1( *model, "integration point variable 1", SCALAR, ELEMENT_INTEGRATION_POINT ),
      IPVariable2( *model, "integration point variable 2", SCALAR, ELEMENT_INTEGRATION_POINT ),
      elementVariable3( *model, "element variable 3", VECTOR, ELEMENT ),
      elementVariable4( *model, "element variable 4", VECTOR, ELEMENT ),
      nodeVariable3( *model, "node variable 3", VECTOR, NODE ),
      nodeVariable4( *model, "node variable 4", VECTOR, NODE ),
      IPVariable3( *model, "integration point variable 3", VECTOR, ELEMENT_INTEGRATION_POINT ),
      IPVariable4( *model, "integration point variable 4", VECTOR, ELEMENT_INTEGRATION_POINT ),
      elementVariable5( *model, "element variable 5", TENSOR, ELEMENT ),
      elementVariable6( *model, "element variable 6", TENSOR, ELEMENT ),
      nodeVariable5( *model, "node variable 5", TENSOR, NODE ),
      nodeVariable6( *model, "node variable 6", TENSOR, NODE ),
      IPVariable5( *model, "integration point variable 5", TENSOR, ELEMENT_INTEGRATION_POINT ),
      IPVariable6( *model, "integration point variable 6", TENSOR, ELEMENT_INTEGRATION_POINT ),
      verbose_(verbose)
  {
  }


/// @attention: removed major memory leak!
PropertyHandle_Test::~PropertyHandle_Test()
  {
     delete model;
  }



void PropertyHandle_Test::run()
    {
        if ( verbose_ ) {
            cout << "\n======================";
            cout << "\nTesting PropertyHandle" << endl;
            cout << "\n======================" << endl;
          }

        // scalars
        elementVariable1 = 1.5;
        elementVariable2 = 2.;
        nodeVariable1 = 3.;
        nodeVariable2 = 4.;
        IPVariable1 = 5.;
        IPVariable2 = 6.;

        // vectors
        elementVariable3 = 7.5;
        elementVariable4 = 8.;
        nodeVariable3 =  9.;
        nodeVariable4 =  10.;
        IPVariable3 = 11.;
        IPVariable4 = 12.;

        // tensors
        elementVariable5 = 13.5;
        elementVariable6 = 14.;
        nodeVariable5 =  15.;
        nodeVariable6 =  16.;
        IPVariable5 = 17.;
        IPVariable6 = 18.;

        //Testing  function operator +=
        if ( verbose_ ) {
            cout << "\nTesting operator += function for scalar variables" << endl;
            cout << "===================================================" << endl;
          }
        ScalarVariable sc;

        // NODE
        Index  nkey   = model->Database().StorageKey( "node variable 2" );
        nodeVariable2 += elementVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc(), 5.5, numeric_limits<double>::epsilon() * 5. );
        nodeVariable2 += IPVariable1;
         (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc(), 10.5, numeric_limits<double>::epsilon() * 10. ); // x10 since the number is double digit

        // ELEMENT
        Index  ekey   = model->Database().StorageKey( "element variable 2" );
        elementVariable2 += nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc(), 5., numeric_limits<double>::epsilon() );
        elementVariable2 += IPVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc(), 10., numeric_limits<double>::epsilon() );

//printRangeOfVariable( *model, "node variable 1", true );
//printRangeOfVariable( *model, "integration point variable 2", true );

        // INTEGRATION POINT
        Index  ipkey = model->Database().StorageKey( "integration point variable 2" );
        IPVariable2 += nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 9., numeric_limits<double>::epsilon() );
        IPVariable2 += elementVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 10.5, numeric_limits<double>::epsilon() );

        if ( verbose_ ) {
            cout << "\nTesting operator += function for vector variables" << endl;
            cout << "===================================================" << endl;
          }
      
        VectorVariable<3> vc;

        // NODE
        nkey   = model->Database().StorageKey( "node variable 4" );
        nodeVariable4 += elementVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0), 17.5, numeric_limits<double>::epsilon() * 10. );
        nodeVariable4  = 1.3;
        IPVariable3    = 0.7;
        nodeVariable4 += IPVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0), 2.0, numeric_limits<double>::epsilon() * 10. );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 4" );
        elementVariable4  = 0.6;
        nodeVariable3     = 0.4;
        elementVariable4 += nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), 1.0, numeric_limits<double>::epsilon() );
        IPVariable3       = 0.7;
        elementVariable4 += IPVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), 1.7, numeric_limits<double>::epsilon() );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 4" );
        IPVariable4 += nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), 12.4, numeric_limits<double>::epsilon() );
        IPVariable4 += elementVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), 19.9, numeric_limits<double>::epsilon() );

        if ( verbose_ ) {
            cout << "\nTesting operator += function for tensor variables" << endl;
            cout << "===================================================" << endl;
          }

        TensorVariable<3> ts;

        // NODE
        nkey   = model->Database().StorageKey( "node variable 6" );
        nodeVariable6 += elementVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0), 29.5, numeric_limits<double>::epsilon() * 100. );
        nodeVariable6 += IPVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0), 46.5, numeric_limits<double>::epsilon() * 100. );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 6" );
        elementVariable6 += nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0), 29., numeric_limits<double>::epsilon() );
        elementVariable6 += IPVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0), 46., numeric_limits<double>::epsilon() );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 6" );
        IPVariable6 += nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0), 33., numeric_limits<double>::epsilon() );
        IPVariable6 += elementVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0), 46.5, numeric_limits<double>::epsilon() );

        /////////////////////////////////////////////////////////////////////////

        // scalars
        elementVariable1 = 1.5;
        elementVariable2 = 2.;
        nodeVariable1 = 3.;
        nodeVariable2 = 4.;
        IPVariable1 = 5.;
        IPVariable2 = 6.;

        // vectors
        elementVariable3 = 7.5;
        elementVariable4 = 8.;
        nodeVariable3 =  9.;
        nodeVariable4 =  10.;
        IPVariable3 = 11.;
        IPVariable4 = 12.;

        // tensors
        elementVariable5 = 13.5;
        elementVariable6 = 14.;
        nodeVariable5 =  15.;
        nodeVariable6 =  16.;
        IPVariable5 = 17.;
        IPVariable6 = 18.;

        //Testing  function operator -=
        if ( verbose_ ) {
            cout << "\nTesting operator -= function for scalar variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 2" );
        nodeVariable2 -= elementVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc(), 2.5, numeric_limits<double>::epsilon() * 10. );
        nodeVariable2 -= IPVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc(), -2.5, TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 2" );
        elementVariable2 -= nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc(), -1, numeric_limits<double>::epsilon() );
        elementVariable2 -= IPVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc(), -6., numeric_limits<double>::epsilon() );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 2" );
        IPVariable2 -= nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 3., numeric_limits<double>::epsilon() );
        IPVariable2 -= elementVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 1.5, numeric_limits<double>::epsilon() );

        if ( verbose_ ) {
            cout << "\nTesting operator -= function for vector variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 4" );
        nodeVariable4 -= elementVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0), 2.5, TOLERANCE );
        nodeVariable4 -= IPVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0) , -8.5, TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 4" );
        elementVariable4 -= nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), -1., numeric_limits<double>::epsilon() );
        elementVariable4 -= IPVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), -12., numeric_limits<double>::epsilon() );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 4" );
        IPVariable4 -= nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), 3., numeric_limits<double>::epsilon() );
        IPVariable4 -= elementVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), -4.5, numeric_limits<double>::epsilon() );

        if ( verbose_ ) {
            cout << "\nTesting operator -= function for tensor variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 6" );
        nodeVariable6 -= elementVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0), 2.5, TOLERANCE );
        nodeVariable6 -= IPVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0) , -14.5, TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 6" );
        elementVariable6 -= nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0) , -1., TOLERANCE );
        elementVariable6 -= IPVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0), -18., TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 6" );
        IPVariable6 -= nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0), 3., TOLERANCE );
        IPVariable6 -= elementVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0), -10.5, TOLERANCE );

        /////////////////////////////////////////////////////////////////////////

        // scalars
        elementVariable1 = 1.5;
        elementVariable2 = 2.;
        nodeVariable1 = 3.;
        nodeVariable2 = 4.;
        IPVariable1 = 5.;
        IPVariable2 = 6.;

        // vectors
        elementVariable3 = 7.5;
        elementVariable4 = 8.;
        nodeVariable3 =  9.;
        nodeVariable4 =  10.;
        IPVariable3 = 11.;
        IPVariable4 = 12.;

        // tensors
        elementVariable5 = 13.5;
        elementVariable6 = 14.;
        nodeVariable5 =  15.;
        nodeVariable6 =  16.;
        IPVariable5 = 17.;
        IPVariable6 = 18.;

        //Testing  function operator *=
        if ( verbose_ ) {
            cout << "\nTesting operator *= function for scalar variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 2" );
        nodeVariable2 *= elementVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc() , 6. , TOLERANCE );
        nodeVariable2 *= IPVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc() , 30., TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 2" );
        elementVariable2 *= nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc(), 6., TOLERANCE );
        elementVariable2 *= IPVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc(), 30., TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 2" );
        IPVariable2 *= nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 18., TOLERANCE );
        IPVariable2 *= elementVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 27., TOLERANCE );

        if ( verbose_ ) {
            cout << "\nTesting operator *= function for vector variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 4" );
        nodeVariable4 *= elementVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0) , 75. , TOLERANCE );
        nodeVariable4 *= IPVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0) , 825. , TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 4" );
        elementVariable4 *= nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), 72., TOLERANCE );
        elementVariable4 *= IPVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), 792., TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 4" );
        IPVariable4 *= nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), 108., TOLERANCE );
        IPVariable4 *= elementVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), 810., TOLERANCE );

        if ( verbose_ ) {
            cout << "\nTesting operator *= function for tensor variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 6" );
        nodeVariable6 *= elementVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0) , 3*216. , TOLERANCE );
        nodeVariable6 *= IPVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0) , 9*3672. , TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 6" );
        elementVariable6 *= nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0) , 3*210. , TOLERANCE );
        elementVariable6 *= IPVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0) , 9*3570. , TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 6" );
        IPVariable6 *= nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0) , 3*270. , TOLERANCE );
        IPVariable6 *= elementVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0) , 9*3645. , TOLERANCE);


        /////////////////////////////////////////////////////////////////////////

        // scalars
        elementVariable1 = 1.5;
        elementVariable2 = 2.;
        nodeVariable1 = 3.;
        nodeVariable2 = 4.;
        IPVariable1 = 5.;
        IPVariable2 = 6.;

        // vectors
        elementVariable3 = 7.5;
        elementVariable4 = 8.;
        nodeVariable3 =  9.;
        nodeVariable4 =  10.;
        IPVariable3 = 11.;
        IPVariable4 = 12.;

        // tensors
        elementVariable5 = 13.5;
        elementVariable6 = 14.;
        nodeVariable5 =  15.;
        nodeVariable6 =  16.;
        IPVariable5 = 17.;
        IPVariable6 = 18.;

        //Testing  function operator /=
        if ( verbose_ ) {
            cout << "\nTesting operator /= function for scalar variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 2" );
        nodeVariable2 /= elementVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc() , 8./3. , TOLERANCE );
        nodeVariable2 /= IPVariable1;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, sc );
        _equal( sc() , 8./15. , TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 2" );
        elementVariable2 /= nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc() , 2./3. , TOLERANCE );
        elementVariable2 /= IPVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, sc );
        _equal( sc() , 2./15. , TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 2" );
        IPVariable2 /= nodeVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 2., TOLERANCE );
        IPVariable2 /= elementVariable1;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, sc );
        _equal( sc(), 4./3., TOLERANCE );

        if ( verbose_ ) {
            cout << "\nTesting operator /= function for vector variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 4" );
        nodeVariable4 /= elementVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0) , 4./3. , TOLERANCE );
        nodeVariable4 /= IPVariable3;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, vc );
        _equal( vc(0) , 4./33. , TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 4" );
        elementVariable4 /= nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0), 8./9., TOLERANCE );
        elementVariable4 /= IPVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, vc );
        _equal( vc(0) , 8./99. , TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 4" );
        IPVariable4 /= nodeVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0), 4./3., TOLERANCE );
        IPVariable4 /= elementVariable3;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, vc );
        _equal( vc(0) , 8./45., TOLERANCE );

        if ( verbose_ ) {
            cout << "\nTesting operator /= function for tensor variables" << endl;
            cout << "=======================" << endl;
          }

        // NODE
        nkey   = model->Database().StorageKey( "node variable 6" );
        nodeVariable6 /= elementVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0) , 32./27. , TOLERANCE );
        nodeVariable6 /= IPVariable5;
        (*(model->Region( "Model" ).NodesBegin()))->Read( nkey, ts );
        _equal( ts(0,0) , 32./27./17. , TOLERANCE );

        // ELEMENT
        ekey   = model->Database().StorageKey( "element variable 6" );
        elementVariable6 /= nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0) , 14./15. , TOLERANCE );
        elementVariable6 /= IPVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( ekey, ts );
        _equal( ts(0,0) , 14./15./17. , TOLERANCE );

        // INTEGRATION POINT
        ipkey   = model->Database().StorageKey( "integration point variable 6" );
        IPVariable6 /= nodeVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0) , 6./5. , TOLERANCE );
        IPVariable6 /= elementVariable5;
        (*(model->Region( "Model" ).CellsBegin()))->Read( 0U, ipkey, ts );
        _equal( ts(0,0) , 4./45. , TOLERANCE);

  }


} // end csmp


