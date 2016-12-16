//
//  IsoparametricQuadraticTetrahedron_Test.cpp
//  
//
//  Created by Stephan Matthai on 7/14/14.
//
//

#include "IsoparametricQuadraticTetrahedron_Test.h"
#include "CSMP_definitions.h"
#include "ErrorHandler.h"
#include "FiniteElement.h"
#include "Node.h"
#include "Element.h"
#include "DenseMatrix.h"

// element to test
#include "IsoparametricQuadraticTriangle.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricQuadraticTetrahedron.h"

using namespace std;

namespace csmp {

/**
    building a quadratic tetrahedron
*/
IsoparametricQuadraticTetrahedron_Test::IsoparametricQuadraticTetrahedron_Test( bool verbose )
 : /* quadratic_triangle(new IsoparametricQuadraticTriangle(dim)), */
   ltetra_(new IsoparametricLinearTetrahedron()),
   qtetra_(new IsoparametricQuadraticTetrahedron()),
   element_(new Element<3U>(qtetra_)),
   element2_(new Element<3U>(ltetra_)),
   tolerance_factor_(10.),
   verbose_(verbose)
 {
    n0.Idx( 1 );
    n1.Idx( 2 );
    n2.Idx( 3 );
    n3.Idx( 4 );
    n4.Idx( 5 );
    n5.Idx( 6 );
    n6.Idx( 7 );
    n7.Idx( 8 );
    n8.Idx( 9 );
    n9.Idx( 10 );

    element_->Idx( 0 );
    element_->Assign( 0, &n0 );
    element_->Assign( 1, &n1 );
    element_->Assign( 2, &n2 );
    element_->Assign( 3, &n3 );
    element_->Assign( 4, &n4 );
    element_->Assign( 5, &n5 );
    element_->Assign( 6, &n6 );
    element_->Assign( 7, &n7 );
    element_->Assign( 8, &n8 );
    element_->Assign( 9, &n9 );
   
    ChangeNodeCoordinatesToTestConfiguration();
   
    // setting the id of the underlying finite element as it will be used to generate VTK output
    element_->FE()->CurrentID(1U);
    element_->CoordinateMatrix();

    // linear tetrahedron
    element2_->Idx( 0 );
    element2_->Assign( 0, &n0 );
    element2_->Assign( 1, &n1 );
    element2_->Assign( 2, &n2 );
    element2_->Assign( 3, &n3 );

    element2_->FE()->CurrentID(1U);
    element2_->CoordinateMatrix();
   
    cout <<"\nconstructor: reference element volume: "<< element2_->Volume() << endl;
}




/// releasing the dynamic memory
IsoparametricQuadraticTetrahedron_Test::~IsoparametricQuadraticTetrahedron_Test()
 {
    delete element_;
    delete qtetra_;

    delete element2_;
    delete ltetra_;
 }


/**
    Derived from IsoparametricQuadraticTetrahedron model in Rhino

    IsoparametricQuadraticTetrahedron.3dm

    //0:   0. 0. 0.
    //1:  20.05217927806071 0.04595452895723895 0.
    //2:  10. 16. 0.
    //3:  10. 6. 17.
    //4:  10. 0. 0.
    //5:  15. 8. 0.
    //6:  5. 8. 0.
    //7:  5.290248577070509 3.174149146242305 8.993422581019864
    //8:  14.72708095415324 3.200084263543058 9.005676182477684
    //9:  10. 10.71037680758673 8.992359427102571

    0:  0.                  0.                      0.
    1:  20.05217927806071   0.04595452895723895     0.
    2:  10.                 16.                     0.
    3:  10.                 6.                      17.
    4:  10.026089639        0.022977264             0.
    5:  15.026089639        8.022977264             0.
    6:  5.                  8.                      0.
    7:  5.                  3.                      8.5
    8:  15.026089639        3.022977264             8.5
    9:  10.                 10.71037680758673       8.5

    Roman, 2014: Corrected the coordinates of midside nodes ( old coordinates are left until the confirmation )
*/
void IsoparametricQuadraticTetrahedron_Test::ChangeNodeCoordinatesToTestConfiguration()
 {
//    // corner nodes
//    element_->N(0)->x( 0. ),                element_->N(0)->y( 0. ),                    element_->N(0)->z( 0. );
//    element_->N(1)->x( 20.05217927806071 ), element_->N(1)->y( 0.04595452895723895 ),   element_->N(1)->z( 0. );
//    element_->N(2)->x( 10. ),               element_->N(2)->y( 16. ),                   element_->N(2)->z( 0.);
//    element_->N(3)->x( 10. ),               element_->N(3)->y( 6. ),                    element_->N(3)->z( 17. );

//    // midside nodes
//    // basis
//    element_->N(4)->x( 10. ), element_->N(4)->y( 0. ),     element_->N(4)->z( 0. );
//    element_->N(5)->x( 15. ), element_->N(5)->y( 8. ),     element_->N(5)->z( 0. );
//    element_->N(6)->x( 5. ),  element_->N(6)->y( 8. ),     element_->N(6)->z( 0. );
//    // sides
//    element_->N(7)->x( 5.290248577070509 ), element_->N(7)->y( 3.174149146242305 ), element_->N(7)->z( 8.993422581019864 );
//    element_->N(8)->x( 14.72708095415324 ), element_->N(8)->y( 3.200084263543058 ), element_->N(8)->z( 9.005676182477684 );
//    element_->N(9)->x( 10. ),               element_->N(9)->y( 10.71037680758673 ), element_->N(9)->z( 8.992359427102571 );

    /// Roman, 2014: Corrected coordinates of midside nodes

    // corner nodes
    element_->N(0)->x( 0. ),                element_->N(0)->y( 0. ),                    element_->N(0)->z( 0. );
    element_->N(1)->x( 20.05217927806071 ), element_->N(1)->y( 0.04595452895723895 ),   element_->N(1)->z( 0. );
    element_->N(2)->x( 10. ),               element_->N(2)->y( 16. ),                   element_->N(2)->z( 0.);
    element_->N(3)->x( 10. ),               element_->N(3)->y( 6. ),                    element_->N(3)->z( 17. );

    // midside nodes
    // basis
    element_->N(4)->x( 10.026089639 ),      element_->N(4)->y( 0.022977264 ),           element_->N(4)->z( 0. );
    element_->N(5)->x( 15.026089639 ),      element_->N(5)->y( 8.022977264 ),           element_->N(5)->z( 0. );
    element_->N(6)->x( 5. ),                element_->N(6)->y( 8. ),                    element_->N(6)->z( 0. );
    // sides
    element_->N(7)->x( 5. ),                element_->N(7)->y( 3. ),                    element_->N(7)->z( 8.5 );
    element_->N(8)->x( 15.026089639 ),      element_->N(8)->y( 3.022977264 ),           element_->N(8)->z( 8.5 );
    element_->N(9)->x( 10. ),               element_->N(9)->y( 11. ),                   element_->N(9)->z( 8.5 );

//    // generic way of calculating midside node coordinates

//    // midside nodes
//    // basis
//    element_->N(4)->Coordinate( ( element_->N(0)->Coordinate() + element_->N(1)->Coordinate() )/2. );
//    element_->N(5)->Coordinate( ( element_->N(1)->Coordinate() + element_->N(2)->Coordinate() )/2. );
//    element_->N(6)->Coordinate( ( element_->N(0)->Coordinate() + element_->N(2)->Coordinate() )/2. );
//    // sides
//    element_->N(7)->Coordinate( ( element_->N(0)->Coordinate() + element_->N(3)->Coordinate() )/2. );
//    element_->N(8)->Coordinate( ( element_->N(1)->Coordinate() + element_->N(3)->Coordinate() )/2. );
//    element_->N(9)->Coordinate( ( element_->N(2)->Coordinate() + element_->N(3)->Coordinate() )/2. );

}


void IsoparametricQuadraticTetrahedron_Test::ChangeNodeCoordinatesToParametric()
 {
   // corner nodes
   element_->N(0)->x(0.); element_->N(0)->y(0.); element_->N(0)->z(0.);
   element_->N(1)->x(1.); element_->N(1)->y(0.); element_->N(1)->z(0.);
   element_->N(2)->x(0.); element_->N(2)->y(1.); element_->N(2)->z(0.);
   element_->N(3)->x(0.); element_->N(3)->y(0.); element_->N(3)->z(1.);
   // midside nodes
   // basis
   element_->N(4)->x(0.5); element_->N(4)->y(0.); element_->N(4)->z(0.);
   element_->N(5)->x(0.5); element_->N(5)->y(0.5); element_->N(5)->z(0.);
   element_->N(6)->x(0.); element_->N(6)->y(0.5); element_->N(6)->z(0.);
   // side edges
   element_->N(7)->x(0.); element_->N(7)->y(0.); element_->N(7)->z(0.5);
   element_->N(8)->x(0.5); element_->N(8)->y(0.); element_->N(8)->z(0.5);
   element_->N(9)->x(0.); element_->N(9)->y(0.5); element_->N(9)->z(0.5);
 }




/**
    Running the tests:
    
    1. Correct representation of element in physical space (node numbering and VTK output)
    
    2. Correct location of integration points
    
    3. Interpolation functions
       - must be 1 at their node and 0 everywhere else
       - their sum should be 1 at any point inside of element
       
    4. Interpolation function derivatives
       - correct location (consistency with interpolation functions)
       - correct orientation in parametric space
       - correct orientation in physical space
       
    5. Jacobian transformation
       - all scaling parameters should be one for an element in physical space with same coordinates as in parametric space
       - lower to higher dim Jacobians
       - efficiency of Jacobian calculations
       
    6. Volume correctness in parametric and physical space
    
    7. Interpolation of integration point values to nodes
*/
void IsoparametricQuadraticTetrahedron_Test::run()
 {
   cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: testing: "<< typeid(IsoparametricQuadraticTetrahedron()).name() <<"\n";
   
   // 0. Output of element for visual check
   // -------------------------------------
   if ( verbose_ ) {
        // coordinates
        DenseMatrix<DM_MIN>  XY(element_->Nodes(),3U);
        element_->CoordinateMatrix( XY );
        cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: Element coordinate matrix:"<< endl;
        XY.Out();
        // nodes
        vector<size_t>  NN(element_->Nodes());
        for ( size_t i=0U; i<NN.size(); ++i ) NN[i] = element_->N(i)->Idx();
        //cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: Node number vector (element):"<< endl;
        //out(NN);
        // dummy data (nodes 1-10 translating integers to numbers)
        DenseMatrix<DM_MIN>  DATA(1,element_->Nodes());
        DATA.Zero();
        DATA(0,0) = 1.0;
        DATA(0,1) = 2.0;
        DATA(0,2) = 3.0;
        DATA(0,3) = 4.0;
        DATA(0,4) = 5.0;
        DATA(0,5) = 6.0;
        DATA(0,6) = 7.0;
        DATA(0,7) = 8.0;
        DATA(0,8) = 9.0;
        DATA(0,9) = 10.0;
        element_->FE()->OutputNodeDataToVTK( "etests", "dummy", DATA );
     
        // placement of integration points
        OutputIntegrationPointsToVTK( "integration_point_data", *element_ );
     }

  // 1. Testing the interpolation functions
  // ----------------------------------------------------
  if ( verbose_ ) cout <<"\n\nIsoparametricQuadraticTetrahedron_Test::run: testing interpolation functions:\n";
  // one at the corresponding node and zero everywhere else?
  TestInterpolationFunctionValues( *element_ );
  // do ipols sum up to one at integration points and barycentre
  TestSumOfInterpolationFunctionValuesEqualTo1( *element_ );

  if ( verbose_ ) {
       //setName( "interactive user-defined interpolation function test" );
       vector<double64> IP(element_->Nodes()), xyz(3U);

        cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: Interpolation functions at point xyz:"<< endl;
        // cout <<"Enter point: ";
        // cin >> xyz[0] >> xyz[1] >> xyz[2];
        xyz[0] = 0.1; xyz[1] = 0.1; xyz[2] = 0.1;
        element_->N_AtGlobalPoint( IP, xyz );
        out(IP);
        cout <<"\tsum interpolation functions: "<< (IP[0]+IP[1]+IP[2]+IP[3]+IP[4]+IP[5]+IP[6]+IP[7]+IP[8]+IP[9]) << endl;
        
        element_->N_AtIntegrationPoint( 2, IP );
        out(IP);
        cout <<"\tsum intpol.f. at integration point: "<< (IP[0]+IP[1]+IP[2]+IP[3]+IP[4]+IP[5]+IP[6]+IP[7]+IP[8]+IP[9]) << endl;

        // computing position of element barycenter
        csmp::Point<3U> bcenter(0.,0.,0.);
        for ( size_t i=0; i<element_->Nodes(); i++ ) bcenter += element_->N(i)->Coordinate();
        bcenter /= static_cast<double64>(element_->Nodes());
        cout <<"\nelement barycenter: "; bcenter.Out();
        // versus 4-node approximation
        bcenter = 0.;
        for ( size_t i=0; i<4; i++ ) bcenter += element_->N(i)->Coordinate();
        bcenter /= 4.;
        cout <<"\nelement barycenter (as based on corner nodes): "; bcenter.Out();
        xyz = bcenter.Coordinates();
        element_->N_AtGlobalPoint( IP, xyz );
        cout <<"\nShape function values calculated at barycentre (global coordinates): "<< endl;
        out(IP);
        cout <<"\tsum Ni: "<< (IP[0]+IP[1]+IP[2]+IP[3]+IP[4]+IP[5]+IP[6]+IP[7]+IP[8]+IP[9]) << endl;

        cout <<"\nShape function values at barycentre (in parametric space): "<< endl;
        element_->N_AtBaryCenter( IP );
        out(IP);
        cout <<"\tsum Ni: "<< (IP[0]+IP[1]+IP[2]+IP[3]+IP[4]+IP[5]+IP[6]+IP[7]+IP[8]+IP[9]) << endl;
     }


  // 2. Testing interpolation function derivatives
  // ---------------------------------------------
  ChangeNodeCoordinatesToParametric();
  element_->Idx( 2 ); // to prompt update of coordinate matrix
  
  //setName( "visual interpolation function-derivative test (parametric space)" );
  // visual test
  DenseMatrix<DM_MIN>  DN(3U,element_->Nodes());
  element_->dN( DN );
  DenseMatrix<DM_MIN> DATA2;
  DATA2 = DN;
  cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: Writing local interpolation function derivatives to file 'etestd2'.";
  element_->FE()->OutputNodeDataToVTK( "etestd", "derivatives", DATA2 );

  element_->dN_AtNode( DN, 7 );
  cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: Element1: Calculated interpolation function derivatives, at node 7 (n=0..n-1):";
  DN.Out();

  ChangeNodeCoordinatesToTestConfiguration();
  element_->Idx( 3 ); // to prompt update of coordinate matrix
  
  //setName( "visual interpolation function-derivative test (physical space)" );
  element_->dN( DN );
  DATA2 = DN;
  cout <<"\nIsoparametricQuadraticTetrahedron_Test::run: Writing global interpolation function derivatives to file 'etestd3'.";
  element_->FE()->OutputNodeDataToVTK( "etestd", "derivatives", DATA2 );


  // 3. Element volume (first in parametric, then in physical space) and barycenter location
  // ---------------------------------------------------------------------------------------
  //setName( "element volume in physical space test" );

  // Volume: +/- 1.0e-7, but calculated from corner nodes only

  //const double64 correct_volume(907.730082);

  const double64 a1 = 20.615528128; // edge 0-3
  const double64 a2 = 20.627577796; // edge 1-3
  const double64 a3 = 19.723082923; // edge 2-3
  const double64 a4 = 20.052231936; // edge 0-1
  const double64 a5 = 18.856772659; // edge 1-2
  const double64 a6 = 18.867962264; // edge 0-2

  //  // generic way of calculating edge lengths

  //  const double64 a1 = (element_->N(0)->Coordinate()-element_->N(3)->Coordinate()).Length(); // edge 0-3
  //  const double64 a2 = (element_->N(1)->Coordinate()-element_->N(3)->Coordinate()).Length(); // edge 1-3
  //  const double64 a3 = (element_->N(2)->Coordinate()-element_->N(3)->Coordinate()).Length(); // edge 2-3
  //  const double64 a4 = (element_->N(0)->Coordinate()-element_->N(1)->Coordinate()).Length(); // edge 0-1
  //  const double64 a5 = (element_->N(1)->Coordinate()-element_->N(2)->Coordinate()).Length(); // edge 1-2
  //  const double64 a6 = (element_->N(0)->Coordinate()-element_->N(2)->Coordinate()).Length(); // edge 0-2

  const double64 s1 = (a1+a2+a4)/2.;
  const double64 s2 = (a2+a3+a5)/2.;
  const double64 s3 = (a3+a6+a1)/2.;
  const double64 s4 = (a4+a5+a6)/2.;

  double64 correct_surface_area = 676.82886432673; // should be around this value
  double64 correct_volume       = 907.73008226913; // should be around this value

  // generic way of calculating volume and surface area

  correct_surface_area = std::sqrt(s1*(s1-a1)*(s1-a2)*(s1-a4))+
                         std::sqrt(s2*(s2-a2)*(s2-a3)*(s2-a5))+
                         std::sqrt(s3*(s3-a3)*(s3-a6)*(s3-a1))+
                         std::sqrt(s4*(s4-a4)*(s4-a5)*(s4-a6));

  correct_volume = std::sqrt( 1./144.*(a1*a1*a5*a5*(a2*a2+a3*a3+a4*a4+a6*a6-a1*a1-a5*a5) +
                                       a2*a2*a6*a6*(a1*a1+a3*a3+a4*a4+a5*a5-a2*a2-a6*a6) +
                                       a3*a3*a4*a4*(a1*a1+a2*a2+a5*a5+a6*a6-a3*a3-a4*a4) -
                                       a1*a1*a2*a2*a4*a4 -
                                       a2*a2*a3*a3*a5*a5 -
                                       a1*a1*a3*a3*a6*a6 -
                                       a4*a4*a5*a5*a6*a6
                                      )
                              );
  // taking all node locations into account
  _equal( element_->Volume(), correct_volume, correct_volume/1000. );
  
  // location of the barycenter in physical space
  //const double64 bc[] = { 10.007, 5.51306, 4.39915 };
  const double64 bc[] = { 10.01304482, 5.511488632, 4.25 };
  Point<3U> bcenter = element_->BaryCenter();
  //setName( "location of barycenter test (x-coordinate)" );
  _equal( bcenter[0], bc[0], 1.0e-5 );
  //setName( "location of barycenter test (y-coordinate)" );
  _equal( bcenter[1], bc[1], 1.0e-5 );
  //setName( "location of barycenter test (z-coordinate)" );
  _equal( bcenter[2], bc[2], 1.0e-5 );

  ChangeNodeCoordinatesToParametric();
  element_->Idx( 4 ); // to prompt update of coordinate matrix
  //setName( "element volume in parametric space test" );
  _equal( element_->Volume(), 1/6., 1.0e-7 );


  // 4. Jacobian transformation (matrix should be unity as we are in parametric space)
  // ---------------------------------------------------------------------------------
  Point<3U> rst_center(0.25,0.25,0.25);
  element_->FE()->JacobianAtIntegrationPoint( 0 ); // near the origin
  const double64 detJ = element_->FE()->JacobianInverse();
  cout <<"\ninverted Jacobian matrix (should be unit-diagonal):";
  element_->FE()->JINV.Out();
  //setName( "testing inverted Jacobian matrix" );
  _equal( element_->FE()->JINV(0,0), 1., 1.0e-7 );
  _equal( element_->FE()->JINV(1,1), 1., 1.0e-7 );
  _equal( element_->FE()->JINV(2,2), 1., 1.0e-7 );
  //setName( "testing Jacobian scaling factor for element in parametric space" );
  _equal( detJ, 1., 1.0e-7 );


  // 5. interpolation and extrapolation of integration point variables to nodes
  // ---------------------------------------------------------------------------------
  ChangeNodeCoordinatesToTestConfiguration();
  element_->Idx( 5 ); // to prompt update of coordinate matrix

  // interpolation from the nodes to the element integration points
  CheckInterpolation( *element_ );

  vector<double64> IPVF(4), NVF(10);
  fill( IPVF.begin(), IPVF.end(), 2. );
        
  // 5.1 linear extrapolation of integration point values to the nodes
  //setName( "integration-point to node extrapolation test 1 (const values=2)" );
  element_->ExtrapolateIntegrationPointVariableToNodes( 1, IPVF, NVF );
  cout <<"\nExtrapolated constant integration point values: ";
  out( NVF );
  for ( int i=0; i<NVF.size(); i++ ) _equal( NVF[i], 2., 1.0e-7 );
  
  // 5.2 variable exhibiting a linear variation (y-coordinate)
  // (going back to the local coordinate system to get precise integration point locations)
  ChangeNodeCoordinatesToParametric();
  element_->Idx( 6 );
  IPVF[0] = (element_->IntegrationPoint(0))[1];
  IPVF[1] = (element_->IntegrationPoint(1))[1];
  IPVF[2] = (element_->IntegrationPoint(2))[1];
  IPVF[3] = (element_->IntegrationPoint(3))[1];
  //setName( "integration-point to node extrapolation test 2 (linear variation)" );
  element_->ExtrapolateIntegrationPointVariableToNodes( 1, IPVF, NVF );
  cout <<"\nExtrapolated gradient integration point values: ";
  out( NVF );
  for ( size_t i=0U; i<element_->Nodes(); i++ )
    _equal( NVF[i], element_->N(i)->y(), 1.0e-7 );


  // 6. consistency between face numbering / normals and interpolation
  // -------------------------------------------------------------------------------------------------------
  ChangeNodeCoordinatesToTestConfiguration();
  element_->Idx( 7 ); // to prompt update of coordinate matrix
  element_->CoordinateMatrix();
  DATA2.Resize(1, element_->Nodes());
  for ( size_t i=0; i<element_->Nodes(); i++ ) DATA2(0,i) = NVF[i];
  element_->FE()->OutputNodeDataToVTK( "encoords", "dummy", DATA2 );
  // checking the normals of the faces for their correct orientation
  OutputFaceNormalsToVTK( "enormals", *element_ );

  // checking the flux balance for element for constant velocity parameter
  //setName( "face-flux consistency test for constant element velocity" );
  CheckElementFaceConsistency( *element_ );

  // report results
  //setName("IsoparametricQuadraticTetrahedron_Test");
  //report();

} // end run







/** Testing the interpolation functions:

    Are the 1 at the corresponding node and zero everywhere else?
*/
void IsoparametricQuadraticTetrahedron_Test::TestInterpolationFunctionValues( const Element<3U>& e )
 {
    //setName( "are ipol-function values 1 at corresponding nodes and zero everywhere else?" );
   
    vector<double64> IPOL(element_->Nodes()), xyz(3U);

    for ( size_t i=0; i<e.Nodes(); i++ ) {
         Point<3U> pt = e.N(i)->Coordinate();
         xyz = pt.Coordinates();
         e.N_AtGlobalPoint( IPOL, xyz );
         // testing
         // 1 at point
         _equal( IPOL[i], 1., tolerance_factor_ * numeric_limits<double64>::epsilon() );
         // zero everywhere else
         for ( size_t j=0; j<IPOL.size(); j++ )
           if ( j != i )
             _equal( IPOL[j], 0., tolerance_factor_ * numeric_limits<double64>::epsilon() );
         // sum = 1 (is given
      }

 } // end




/**
    Tests that interpolation function values sum up to one at:
    - integration points
    - barycentre
    
    @attention IPOL values at the nodes are tested by TestInterpolationFunctionValues
*/
void IsoparametricQuadraticTetrahedron_Test::TestSumOfInterpolationFunctionValuesEqualTo1( const Element<3U>& e )
 {
    //setName( "do ipol-function sum to 1 at ips and barycentre ?" );

    vector<double64> IPOL(element_->Nodes());

    // at integration points
    for ( size_t i=0; i<e.IntegrationPoints(); i++ ) {
         e.N_AtIntegrationPoint( i, IPOL );
         double64 sum(0.);
         for ( size_t j=0; j<IPOL.size(); j++ ) sum += IPOL[j];
        _equal( sum, 1., tolerance_factor_ * numeric_limits<double64>::epsilon() );
      }
   
    // at element barycenter
    e.N_AtBaryCenter( IPOL );
    double64 sum(0.);
    for ( size_t j=0; j<IPOL.size(); j++ ) sum += IPOL[j];
  _equal( sum, 1., tolerance_factor_ * numeric_limits<double64>::epsilon() );
   
 } // end
 







/**
     Assigns a nodal field of variable values that then gets interpolated 
     across the element. The values of this field variable are equivalent to the
     position ot the integration points in physical space.
     This relationship is used for the testing.
     
     @test OK SKM 3/9/14
*/
void IsoparametricQuadraticTetrahedron_Test::CheckInterpolation( const Element<3U>& e )
 {
    //setName( "interpolation function node-value- to integration-point interpolation test" );
   
    // 0. generating values at integration points from coordinates
    // ------------------------------------------------------------
    // using the negative sum of the global coordinates as an interpolant:
    // At all nodes, val = -x + -y + -z.
    vector<double64>  node_vals(e.Nodes());
    for ( size_t i=0U; i<e.Nodes(); i++ )
      node_vals[i] = -e.N(i)->x() + -e.N(i)->y() + -e.N(i)->z();

    // interpolating nodal values to integration points and checking these against theoretic values
    vector<double64> IPOL;
    for ( size_t i=0U; i<element_->IntegrationPoints(); i++ ) {
         // interpolation
         e.N_AtIntegrationPoint( i, IPOL );
         double64 ivalue(0.);
         for ( size_t j=0; j<e.Nodes(); j++ )
           ivalue += IPOL[j] * node_vals[j];
         // checking ivalue (IntegrationPoint returns the location of the integration point in physical space)
         Point<3U> xyz(element_->IntegrationPoint(i));
         double64 correct_value = -xyz[0] + -xyz[1] + -xyz[2];
        _equal( ivalue, correct_value, tolerance_factor_ * numeric_limits<double64>::epsilon() );
      }
   
 } // end CheckInterpolation






/**
    Checks that the face normals are pointing correctly and computes face flux balance 
    over all faces of the element. This should be zero.
    
    @todo check how accurate this result can be expected to be?
*/
void IsoparametricQuadraticTetrahedron_Test::CheckElementFaceConsistency( const Element<3U>& e )
 {
    vector<double64>   nrml;
    VectorVariable<3>  flux(PLAIN,PLAIN,PLAIN,0.3,0.5,1.0);
    vector<double64>   IPOL(6);
    double64           flux_balance(0.);
 
    cout <<"\nElement "<< e.Idx() <<", volume: "<< e.Volume();
 
    // 0. generating face normals, scaling and storing them
    // -----------------------------------------------------------
    for ( size_t i=0; i<e.Faces(); i++ ) {
         // computing and storing normal to face
         e.FE()->UnitNormalAtFaceBarycenter( i, nrml );
         VectorVariable<3U>  fn(nrml);
         // projecting flux onto normal
         // dot product fn . vc
         double64 projflux  = fn & flux;
         cout <<"\nFace "<< i <<", projected flux:         "<< projflux;
         
         // integrating flux over the face
         double64 fflux = projflux * e.FaceArea( i );
      
         // summing flux balance
         flux_balance += fflux;  
         cout <<"\nFace "<< i <<", integrated normal flux: "<< fflux;         
      }
    cout <<"\n";

   _equal( flux_balance, 0., tolerance_factor_ * numeric_limits<double64>::epsilon() );
    cout <<"\nFlux balance of element "<< e.Idx() <<": "<< flux_balance << endl;
      
 } // end CheckElementFaceConsistency





/**
    The integration points are coloured by number and are then output to VTK file.
    Point 1=0, point2=1...n-1.
*/
void IsoparametricQuadraticTetrahedron_Test::OutputIntegrationPointsToVTK( const char* file, const Element<3U>& e ) const
 {
     // 0. generating values, first integration point=0, second=1...
     // ------------------------------------------------------------
    vector<pair<Point<3U>,double64> >  ipoint_data(element_->IntegrationPoints()); // locations in physical space
    // convention: integration points get values equivalent to their number 0..n-1
    double64 display_value(0.);
    for ( size_t i=0U; i<element_->IntegrationPoints(); i++ ) {
         Point<3U> xyz = element_->IntegrationPoint(i);
         ipoint_data[i] = make_pair( xyz, display_value );
         display_value += 1.;
      }
   
     // 1. writing the file header
     // --------------------------
     string  outfile(file);
     outfile +=".vtk";
     ofstream  ofs(outfile.c_str());
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"integration points colored by number: '";
     ofs << "point_number" <<"'"<< endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3. writing point coordinates in order
     // --------------------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << ipoint_data.size() <<" float"<< endl;
     for ( vector<pair<Point<3U>,double64> >::const_iterator
           it=ipoint_data.begin(); it!=ipoint_data.end(); it++ ) {
          ofs << (*it).first[0] <<" "<< (*it).first[1] <<" "<< (*it).first[2] << endl;
       }
     ofs << endl;  
       
     // 4. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     // (each record has a cell number and the corresponding point specifier)
     ofs <<"CELLS "<< ipoint_data.size() <<" "<< (ipoint_data.size() * 2) << endl;
     for ( size_t counter(0); counter<ipoint_data.size(); counter++ ) {
          ofs << 1 <<" "<< counter << endl;
       }
     ofs << endl;
     
     // 5. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< ipoint_data.size() << endl;
     // cell type VTK_POINT (=1)
     for ( size_t counter(0); counter<ipoint_data.size(); counter++ ) ofs << 1 << endl;
     ofs << endl;

     // 6. writing POINT_DATA point-type data values
     // ----------------------------------------------------
     ofs <<"POINT_DATA "<< ipoint_data.size() << endl;
     ofs <<"SCALARS "<< "point_number" <<" float"<< endl;
     ofs <<"LOOKUP_TABLE default" << endl;
     for ( vector<pair<Point<3U>,double64> >::const_iterator
           it=ipoint_data.begin(); it!=ipoint_data.end(); it++ ) {
          ofs << (*it).second <<" "<< endl;
       }
     ofs << endl;
     if ( verbose_ ) cout <<"\nVTK file '"<< outfile <<"' written successfully."<< endl;
      
 } // end OutputIntegrationPointsToVTK






/**
    Loops over the faces of the element, collecting the (outward-pointing) face normals
    for display in VTK.
    Each normal is scaled by the sqrt of the element volume and placed 
    with its origin on the barycenter of the corresponding face.
*/
void IsoparametricQuadraticTetrahedron_Test::OutputFaceNormalsToVTK( const char* file, const Element<3U>& e ) const
 {
    vector<double64>    nrml;
    size_t              counter(0);
    map<size_t,pair<Point<3>,Point<3> > >  normals;

    // 0. generating face normals, scaling and storing them
    // -----------------------------------------------------------
    vector<size_t> fnids;
    for ( size_t i=0; i<e.Faces(); i++ ) {
         // computing and storing normal to face
         e.FE()->UnitNormalAtFaceBarycenter( i, nrml );
         Point<3U> unormal(nrml);
         // scaling the normals (by empirical factor)
         unormal *= (sqrt(e.Volume()) / 10.);
         // finding root points for the normals = barycenters of faces
         Point<3U> bctr = e.FaceBaryCenter( i );
         // storing starting pointsd and normals in map
         normals[ counter++ ] = make_pair( bctr, unormal );
      }  
   
     // 1. writing the file header
     // --------------------------
     string  outfile(file);
     outfile +=".vtk";
     ofstream  ofs(outfile.c_str());
     ofs <<"# vtk DataFile Version 4.2"<< endl;
     ofs <<"normals to faces of element: '";
     ofs << "normals" <<"'"<< endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3. writing point coordinates in order
     // --------------------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << e.Faces() <<" double"<< endl;
     for ( map<size_t,pair<Point<3>,Point<3> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) {
          ofs << (*it).second.first[0] <<" "<<  (*it).second.first[1] <<" "<<  (*it).second.first[2] << endl;
       }
     ofs << endl;  
       
     // 4. writing CELLS (cell-size and member face (points))
     // -----------------------------------------------------
     // (each normal=Edge has two points plus a specifier)
     ofs <<"CELLS "<< e.Faces() <<" "<< (e.Faces() * 2) << endl;
     counter=0;
     for ( map<size_t,pair<Point<3>,Point<3> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) {
          ofs <<"1 "<< counter++ << endl;
       }
     ofs << endl;
     
     // 5. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< normals.size() << endl;
     // cell type VTK_VERTEX (=1)
     for ( map<size_t,pair<Point<3>,Point<3> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) ofs << 1 << endl;
     ofs << endl;

     // 6. writing POINT_DATA point-type data values
     // ----------------------------------------------------
     ofs <<"POINT_DATA "<< normals.size() << endl;
     ofs <<"VECTORS "<< "normals" <<" double"<< endl;
     for ( map<size_t,pair<Point<3>,Point<3> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) {
          ofs << (*it).second.second[0] <<" "<<  (*it).second.second[1] <<" "<<  (*it).second.second[2] << endl;
       }
     ofs << endl;
     if ( verbose_ ) cout <<"\nVTK file '"<< outfile <<"' written successfully."<< endl;
      
 } // end OutputFaceNormalsToVTK






} // end csmp
