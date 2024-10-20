//
//  IsoparametricLinearPyramid_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 15/10/2024.
//

#include "IsoparametricLinearPyramid_Test.h"
#include "IsoparametricLinearPyramid.h"
#include "Element.h"
#include "Node.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

void IsoparametricLinearPyramid_Test::run()
 {
    InterpolationFunctionAtIntegrationPointsTest();
    InterpolationFunctionDerivative_Test();
    InterpolationFunctionDerivative_Test();
    Volume_Test();
 }
 
 

bool IsoparametricLinearPyramid_Test::SumOfInterpolationFunctionTest() // at integration points
 {
    return true;
 }
 
 
Element<3U> IsoparametricLinearPyramid_Test::CreateLocalCoordinatePyramidElement( int n_integration_points )
 {
    // creating a pyramid with the same coordinates as in parametric space
    static IsoparametricLinearPyramid pyra( n_integration_points );
    Element<3> elmt( &pyra );
    elmt.Idx( 0 );
    Node<3> n1, n2, n3, n4, n5;
    
    DenseMatrix<DM_MIN> NXYZ(5,3);
    NXYZ(0,0) =-1.0;NXYZ(0,1) =-1.0;NXYZ(0,2) =  0.0;
    NXYZ(1,0) = 1.0;NXYZ(1,1) =-1.0;NXYZ(1,2) =  0.0;
    NXYZ(2,0) = 1.0;NXYZ(2,1) = 1.0;NXYZ(2,2) =  0.0;
    NXYZ(3,0) =-1.0;NXYZ(3,1) = 1.0;NXYZ(3,2) =  0.0;
    NXYZ(4,0) = 0.0;NXYZ(4,1) = 0.0;NXYZ(4,2) =  1.0;

    n1.Coordinate( Point<3>( NXYZ(0,0), NXYZ(0,1), NXYZ(0,2) ) ); n1.Idx(0);
    n2.Coordinate( Point<3>( NXYZ(1,0), NXYZ(1,1), NXYZ(1,2) ) ); n2.Idx(1);
    n3.Coordinate( Point<3>( NXYZ(2,0), NXYZ(2,1), NXYZ(2,2) ) ); n3.Idx(2);
    n4.Coordinate( Point<3>( NXYZ(3,0), NXYZ(3,1), NXYZ(3,2) ) ); n4.Idx(3);
    n5.Coordinate( Point<3>( NXYZ(4,0), NXYZ(4,1), NXYZ(4,2) ) ); n5.Idx(4);
    
    elmt.Assign( 0, &n1 );
    elmt.Assign( 1, &n2 );
    elmt.Assign( 2, &n3 );
    elmt.Assign( 3, &n4 );
    elmt.Assign( 4, &n5 );
    
    elmt.CoordinateMatrix();
    
    if ( verbose_ ) {
         cout <<"\n\n"<<"IsoparametricLinearPyramid_Test::CreateLocalCoordinatePyramidElement: ";
         cout <<"\n\t"<<"printing pyramid as in reference coordinates"<< endl;
         elmt.Out();
         DenseMatrix<DM_MIN> DATA;
         elmt.dN(DATA);
         pyra.OutputNodeDataToVTK( "Pyramid_in_ref_coords", "nodal_derivatives", DATA );
      }

    return elmt;
 }

 
/**
        As constructed in Rhino file:  SkewedPyramidFlatBase.3dm
        
        volumer = 833.3333333333333
*/
Element<3U> IsoparametricLinearPyramid_Test::CreateDistortedPyramidElement( int n_integration_points )
 {
    // creating a pyramid with the same coordinates as in parametric space
    static IsoparametricLinearPyramid pyra( n_integration_points );
    Element<3> elmt( &pyra );
    elmt.Idx( 1 );
    Node<3> n1, n2, n3, n4, n5;
    
    DenseMatrix<DM_MIN> NXYZ(5,3);
    NXYZ(0,0) =  0.0; NXYZ(0,1) =  0.0; NXYZ(0,2) =  0.0;
    NXYZ(1,0) = 13.0; NXYZ(1,1) =  2.0; NXYZ(1,2) =  0.0;
    NXYZ(2,0) = 14.0; NXYZ(2,1) = 11.0; NXYZ(2,2) =  0.0;
    NXYZ(3,0) =  1.0; NXYZ(3,1) =  9.0; NXYZ(3,2) =  0.0;
    NXYZ(4,0) = 24.0; NXYZ(4,1) =  2.0; NXYZ(4,2) = 22.0;

    n1.Coordinate( Point<3>( NXYZ(0,0), NXYZ(0,1), NXYZ(0,2) ) ); n1.Idx(0);
    n2.Coordinate( Point<3>( NXYZ(1,0), NXYZ(1,1), NXYZ(1,2) ) ); n2.Idx(1);
    n3.Coordinate( Point<3>( NXYZ(2,0), NXYZ(2,1), NXYZ(2,2) ) ); n3.Idx(2);
    n4.Coordinate( Point<3>( NXYZ(3,0), NXYZ(3,1), NXYZ(3,2) ) ); n4.Idx(3);
    n5.Coordinate( Point<3>( NXYZ(4,0), NXYZ(4,1), NXYZ(4,2) ) ); n5.Idx(4);

    elmt.Assign( 0, &n1 );
    elmt.Assign( 1, &n2 );
    elmt.Assign( 2, &n3 );
    elmt.Assign( 3, &n4 );
    elmt.Assign( 4, &n5 );

    elmt.CoordinateMatrix();

    if ( verbose_ ) {
         cout <<"\n\n"<<"IsoparametricLinearPyramid_Test::CreateLocalCoordinatePyramidElement: ";
         cout <<"\n\t"<<"printing distored pyramid with planar base"<< endl;
         elmt.Out();
         DenseMatrix<DM_MIN> DATA;
         elmt.dN(DATA);
         pyra.OutputNodeDataToVTK( "Pyramid_distorted1", "nodal_derivatives", DATA );
      }

    return elmt;
 }
 
 
 
 
bool IsoparametricLinearPyramid_Test::InterpolationFunctionAtIntegrationPointsTest() // at integration points
 {
    Element<3> elmt = CreateLocalCoordinatePyramidElement(5);
    
    // N at quadrature points 5 ips, 5 N-function values
//    vector<vector<double>> qpoint_vals = { {0.103647450843758,0.103647450843758,0.103647450843758,0.103647450843758,0.585410196624968},
//                                           {0.484014,0.161338,0.053779,0.161338,0.138197},
//                                           {0.161338,0.484014,0.161338,0.053779,0.138197},
//                                           {0.053779,0.161338,0.484014,0.161338,0.138197},
//                                           {0.161338,0.053779,0.161338,0.484014,0.138197} };
    vector<vector<double>> qpoint_vals = { {0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000},
                                           {0.455341801262, 0.122008467929, 0.033333333333, 0.122008467929, 0.267307929547},
                                           {0.122008467929, 0.455341801262, 0.122008467929, 0.033333333333, 0.267307929547},
                                           {0.033333333333, 0.122008467929, 0.455341801262, 0.122008467929, 0.267307929547},
                                           {0.122008467929, 0.033333333333, 0.122008467929, 0.455341801262, 0.267307929547} };
     vector<double> nrst(5,0.);
     for ( int i{0}; i<elmt.Nodes(); ++i ) {
          elmt.N_AtIntegrationPoint( i, nrst );
          // testing the value of each interpolation function at each integration point
          _test( approximatelyEqual( nrst[0], qpoint_vals[i][0]) );
          _test( approximatelyEqual( nrst[1], qpoint_vals[i][1]) );
          _test( approximatelyEqual( nrst[2], qpoint_vals[i][2]) );
          _test( approximatelyEqual( nrst[3], qpoint_vals[i][3]) );
          _test( approximatelyEqual( nrst[4], qpoint_vals[i][4]) );
          // testing that the functions sum-up to 1
          double sum{0.};
          for ( const auto& nf : nrst ) sum += nf;
          _test( approximatelyEqual( sum,  1. ) );
       }

    return true;
 }


 
bool IsoparametricLinearPyramid_Test::InterpolationFunctionDerivative_Test() // at the nodes
 {
     // starting in local coordinates
     Element<3U> elmt = CreateLocalCoordinatePyramidElement(5);
     
     // derivatives in double precision at the 5 quadrature points (one for each row
     DenseMatrix<DM_MIN>  dN1 = { {-0.103647453303, -0.103647453303, -0.103647453303},
                                  {-0.161337922884, -0.161337922884, -0.484013761081},
                                  {-0.053779530972, -0.161337922884, -0.161337922884},
                                  {-0.161337922884, -0.053779530972, -0.161337922884},
                                  {-0.484013761081, -0.161337922884, -0.161337922884} };

     DenseMatrix<DM_MIN>  dN2 = { {0.103647453303, -0.103647453303, -0.103647453303},
                                  {0.161337922884, -0.053779530972, -0.161337922884},
                                  {0.053779530972, -0.484013761081, -0.161337922884},
                                  {0.161337922884, -0.161337922884, -0.161337922884},
                                  {0.161337922884, -0.053779530972, -0.161337922884} };

     DenseMatrix<DM_MIN>  dN3 = { {0.103647453303, 0.103647453303, -0.103647453303},
                                  {0.053779530972, 0.053779530972, -0.161337922884},
                                  {0.484013761081, 0.161337922884, -0.161337922884},
                                  {0.161337922884, 0.484013761081, -0.161337922884},
                                  {0.053779530972, 0.161337922884, -0.161337922884} };

     DenseMatrix<DM_MIN>  dN4 = { {-0.103647453303, 0.103647453303, -0.103647453303},
                                  {-0.053779530972, 0.161337922884, -0.161337922884},
                                  {-0.161337922884, 0.053779530972, -0.161337922884},
                                  {-0.053779530972, 0.161337922884, -0.161337922884},
                                  {-0.161337922884, 0.484013761081, -0.161337922884} };

     DenseMatrix<DM_MIN>  dN5 = { {0.000000000000, 0.000000000000, 1.000000000000},
                                  {0.000000000000, 0.000000000000, 1.000000000000},
                                  {0.000000000000, 0.000000000000, 1.000000000000},
                                  {0.000000000000, 0.000000000000, 1.000000000000},
                                  {0.000000000000, 0.000000000000, 1.000000000000} };
     // testing
     const uint32_t dimensions{3u};
     DenseMatrix<DM_MIN> DN;
     // integration point 1
     elmt.dN_AtIntegrationPoint( DN, 0, dimensions );
     for ( uint32_t i{0u}; i<elmt.Nodes(); i++ )
       for ( uint32_t j{0u}; j<dimensions; j++ )
         _test( approximatelyEqual( DN(i,j), dN1(i,j) ) );

     // integration point 2
     elmt.dN_AtIntegrationPoint( DN, 1, dimensions );
     for ( uint32_t i{0u}; i<elmt.Nodes(); i++ )
       for ( uint32_t j{0u}; j<dimensions; j++ )
         _test( approximatelyEqual( DN(i,j), dN1(i,j) ) );

     // integration point 3
     elmt.dN_AtIntegrationPoint( DN, 2, dimensions );
     for ( uint32_t i{0u}; i<elmt.Nodes(); i++ )
       for ( uint32_t j{0u}; j<dimensions; j++ )
         _test( approximatelyEqual( DN(i,j), dN1(i,j) ) );

     // integration point 4
     elmt.dN_AtIntegrationPoint( DN, 3, dimensions );
     for ( uint32_t i{0u}; i<elmt.Nodes(); i++ )
       for ( uint32_t j{0u}; j<dimensions; j++ )
         _test( approximatelyEqual( DN(i,j), dN1(i,j) ) );

     // integration point 5
     elmt.dN_AtIntegrationPoint( DN, 4, dimensions );
     for ( uint32_t i{0u}; i<elmt.Nodes(); i++ )
       for ( uint32_t j{0u}; j<dimensions; j++ )
         _test( approximatelyEqual( DN(i,j), dN1(i,j) ) );

     return true;
 }
 
 
 
 
bool IsoparametricLinearPyramid_Test::Volume_Test()
 {
    { // reference element
      Element<3U> elmt = CreateLocalCoordinatePyramidElement(8);

      const double vol_ref_elmt{ 4./3. };
      const double computed_pyra_volume{ elmt.Volume() };

      _test( approximatelyEqual( computed_pyra_volume, vol_ref_elmt ) );
    }
    { // distorted element
      Element<3U> elmt = CreateDistortedPyramidElement(8);

      const double vol_distorted_pyra_planar_base{ 843.3333333333333 };
      const double computed_pyra_volume{ elmt.Volume() };
      _test( approximatelyEqual( computed_pyra_volume, vol_distorted_pyra_planar_base ) );
      
      // checking integration points and weights
      vector<double> xyz;
      cout <<"\n"<<"current integration points\n";
      for ( uint32_t i{0u}; i<elmt.IntegrationPoints(); i++ ) {
           elmt.FE()->IntegrationPoint( i, xyz );
           cout <<"\n\t"<<"ip: "<< i <<"("<< xyz[0] <<","<< xyz[1] <<","<< xyz[2] <<"), ";
           cout <<"weight: "<< elmt.FE()->WeightAtIntegrationPoint( i );
        }
    }    
    return true;
 }


} // end csmp
