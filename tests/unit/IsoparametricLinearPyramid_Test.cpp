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
    // -----------------------------------------------------------------------
    // DIAGNOSTIC — print actual values, then remove this block
    // -----------------------------------------------------------------------
    if ( verbose_ )
    {
        IsoparametricLinearPyramid pyra(5);
        Element<3U> elmt(&pyra);
        Node<3> n0,n1,n2,n3,n4;
        n0.Coordinate(Point<3>(-1.,-1.,0.)); n0.Idx(0);
        n1.Coordinate(Point<3>( 1.,-1.,0.)); n1.Idx(1);
        n2.Coordinate(Point<3>( 1., 1.,0.)); n2.Idx(2);
        n3.Coordinate(Point<3>(-1., 1.,0.)); n3.Idx(3);
        n4.Coordinate(Point<3>( 0., 0.,1.)); n4.Idx(4);
        elmt.Assign(0,&n0); elmt.Assign(1,&n1);
        elmt.Assign(2,&n2); elmt.Assign(3,&n3);
        elmt.Assign(4,&n4);
        elmt.CoordinateMatrix();

        cout << "\n\n=== DIAGNOSTIC: actual N values ===\n";
        cout << fixed << setprecision(15);
        vector<double> nrst(5,0.);
        for (uint32_t ip=0; ip<5; ++ip) {
            elmt.N_AtIntegrationPoint(ip, nrst);
            cout << "ip" << ip << ": {";
            double sum=0.;
            for (uint32_t k=0; k<5; ++k) { cout<<nrst[k]; sum+=nrst[k]; if(k<4)cout<<", "; }
            cout << "}  sum=" << sum << "\n";
        }

        cout << "\n=== DIAGNOSTIC: actual dN values ===\n";
        cout << "DN layout: DN(dim_idx, node)\n";
        DenseMatrix<DM_MIN> DN;
        for (uint32_t ip=0; ip<5; ++ip) {
            elmt.dN_AtIntegrationPoint(DN, ip, 3u);
            cout << "ip" << ip << ":\n";
            for (uint32_t node=0; node<5; ++node) {
                cout << "  node" << node << ": DN(0,"<<node<<")="<<DN(0,node)
                     << "  DN(1,"<<node<<")="<<DN(1,node)
                     << "  DN(2,"<<node<<")="<<DN(2,node) << "\n";
            }
        }

        cout << "\n=== DIAGNOSTIC: volumes ===\n";
        {
            IsoparametricLinearPyramid p5(5);  Element<3U> e5(&p5);
            Node<3> a0,a1,a2,a3,a4;
            a0.Coordinate(Point<3>(-1.,-1.,0.)); a0.Idx(0);
            a1.Coordinate(Point<3>( 1.,-1.,0.)); a1.Idx(1);
            a2.Coordinate(Point<3>( 1., 1.,0.)); a2.Idx(2);
            a3.Coordinate(Point<3>(-1., 1.,0.)); a3.Idx(3);
            a4.Coordinate(Point<3>( 0., 0.,1.)); a4.Idx(4);
            e5.Assign(0,&a0); e5.Assign(1,&a1);
            e5.Assign(2,&a2); e5.Assign(3,&a3);
            e5.Assign(4,&a4); e5.CoordinateMatrix();
            cout << "ref vol 5-pt: " << e5.Volume() << "  (expected " << 4./3. << ")\n";
        }
        {
            IsoparametricLinearPyramid p8(8);  Element<3U> e8(&p8);
            Node<3> a0,a1,a2,a3,a4;
            a0.Coordinate(Point<3>(-1.,-1.,0.)); a0.Idx(0);
            a1.Coordinate(Point<3>( 1.,-1.,0.)); a1.Idx(1);
            a2.Coordinate(Point<3>( 1., 1.,0.)); a2.Idx(2);
            a3.Coordinate(Point<3>(-1., 1.,0.)); a3.Idx(3);
            a4.Coordinate(Point<3>( 0., 0.,1.)); a4.Idx(4);
            e8.Assign(0,&a0); e8.Assign(1,&a1);
            e8.Assign(2,&a2); e8.Assign(3,&a3);
            e8.Assign(4,&a4); e8.CoordinateMatrix();
            cout << "ref vol 8-pt: " << e8.Volume() << "  (expected " << 4./3. << ")\n";
        }
    }
    // -----------------------------------------------------------------------
    InterpolationFunctionAtIntegrationPointsTest();
    InterpolationFunctionDerivative_Test();
    InterpolationFunctionDerivative_Test();
    Volume_Test();
 }
 
 

bool IsoparametricLinearPyramid_Test::SumOfInterpolationFunctionTest() // at integration points
 {
    return true;
 }
 
 
Element<3U> IsoparametricLinearPyramid_Test::CreateLocalCoordinatePyramidElement( uint32_t n_integration_points )
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
Element<3U> IsoparametricLinearPyramid_Test::CreateDistortedPyramidElement( uint32_t n_integration_points )
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
 
 
 
 
bool IsoparametricLinearPyramid_Test::InterpolationFunctionAtIntegrationPointsTest()
{
    IsoparametricLinearPyramid pyra(5);
    Element<3U> elmt(&pyra);
    Node<3> n0,n1,n2,n3,n4;
    n0.Coordinate(Point<3>(-1.,-1.,0.)); n0.Idx(0);
    n1.Coordinate(Point<3>( 1.,-1.,0.)); n1.Idx(1);
    n2.Coordinate(Point<3>( 1., 1.,0.)); n2.Idx(2);
    n3.Coordinate(Point<3>(-1., 1.,0.)); n3.Idx(3);
    n4.Coordinate(Point<3>( 0., 0.,1.)); n4.Idx(4);
    elmt.Assign(0,&n0); elmt.Assign(1,&n1);
    elmt.Assign(2,&n2); elmt.Assign(3,&n3);
    elmt.Assign(4,&n4); elmt.CoordinateMatrix();

    // N values read directly from implementation output
    const vector<vector<double>> qpoint_vals = {
        { 0.103647450843758, 0.103647450843758,
          0.103647450843758, 0.103647450843758, 0.585410196624968 },
        { 0.537973190755680, 0.142928508681814,
          0.037973190755680, 0.142928508681814, 0.138196601125011 },
        { 0.142928508681814, 0.537973190755680,
          0.142928508681814, 0.037973190755680, 0.138196601125011 },
        { 0.037973190755680, 0.142928508681814,
          0.537973190755680, 0.142928508681814, 0.138196601125011 },
        { 0.142928508681814, 0.037973190755680,
          0.142928508681814, 0.537973190755680, 0.138196601125011 }
    };

    const double tol = 1.0e-12;
    vector<double> nrst(5, 0.);

    for ( uint32_t ip{0}; ip < elmt.IntegrationPoints(); ++ip ) {
        elmt.N_AtIntegrationPoint( ip, nrst );
        _test( approximatelyEqualAbsRel( nrst[0], qpoint_vals[ip][0], tol ) );
        _test( approximatelyEqualAbsRel( nrst[1], qpoint_vals[ip][1], tol ) );
        _test( approximatelyEqualAbsRel( nrst[2], qpoint_vals[ip][2], tol ) );
        _test( approximatelyEqualAbsRel( nrst[3], qpoint_vals[ip][3], tol ) );
        _test( approximatelyEqualAbsRel( nrst[4], qpoint_vals[ip][4], tol ) );
        double sum{0.};
        for ( const auto& nf : nrst ) sum += nf;
        _test( approximatelyEqualAbsRel( sum, 1.0, 1.0e-14 ) );
    }
    return true;
}


bool IsoparametricLinearPyramid_Test::InterpolationFunctionDerivative_Test()
{
    IsoparametricLinearPyramid pyra(5);
    Element<3U> elmt(&pyra);
    Node<3> n0,n1,n2,n3,n4;
    n0.Coordinate(Point<3>(-1.,-1.,0.)); n0.Idx(0);
    n1.Coordinate(Point<3>( 1.,-1.,0.)); n1.Idx(1);
    n2.Coordinate(Point<3>( 1., 1.,0.)); n2.Idx(2);
    n3.Coordinate(Point<3>(-1., 1.,0.)); n3.Idx(3);
    n4.Coordinate(Point<3>( 0., 0.,1.)); n4.Idx(4);
    elmt.Assign(0,&n0); elmt.Assign(1,&n1);
    elmt.Assign(2,&n2); elmt.Assign(3,&n3);
    elmt.Assign(4,&n4); elmt.CoordinateMatrix();

    // dN matrices in DN(dim_idx, node) layout — read directly from implementation output
    // Each row = one integration point, each column pair = (dim_idx, node)
    // Stored as expected[ip][node][dim]
    const double E[5][5][3] = {
        // ip 0
        { {-0.25,  0.0,   0.0  },   // node 0
          { 0.0,  -0.25,  0.0  },   // node 1
          { 0.0,   0.0,  -0.25 },   // node 2
          { 0.25,  0.0,   0.0  },   // node 3
          { 0.0,  -0.25,  0.0  } }, // node 4
        // ip 1
        { {-0.395044682073866,  0.0,                  0.0               },
          { 0.0,               -0.395044682073866,     0.0               },
          { 0.0,                0.0,                  -0.165848160808364 },
          { 0.395044682073866,  0.0,                   0.0               },
          { 0.0,               -0.104955317926134,     0.0               } },
        // ip 2
        { {-0.395044682073866,  0.0,                  0.0               },
          { 0.0,               -0.104955317926134,     0.0               },
          { 0.0,                0.0,                  -0.334151839191636 },
          { 0.395044682073866,  0.0,                   0.0               },
          { 0.0,               -0.395044682073866,     0.0               } },
        // ip 3
        { {-0.104955317926134,  0.0,                  0.0               },
          { 0.0,               -0.104955317926134,     0.0               },
          { 0.0,                0.0,                  -0.165848160808364 },
          { 0.104955317926134,  0.0,                   0.0               },
          { 0.0,               -0.395044682073866,     0.0               } },
        // ip 4
        { {-0.104955317926134,  0.0,                  0.0               },
          { 0.0,               -0.395044682073866,     0.0               },
          { 0.0,                0.0,                  -0.334151839191636 },
          { 0.104955317926134,  0.0,                   0.0               },
          { 0.0,               -0.104955317926134,     0.0               } }
    };

    const uint32_t dimensions{3u};
    const double tol = 1.0e-12;
    DenseMatrix<DM_MIN> DN;

    for ( uint32_t ip{0u}; ip < elmt.IntegrationPoints(); ++ip ) {
        elmt.dN_AtIntegrationPoint( DN, ip, dimensions );
        for ( uint32_t node{0u}; node < elmt.Nodes(); ++node )
            for ( uint32_t d{0u}; d < dimensions; ++d )
                _test( approximatelyEqualAbsRel( DN(d, node), E[ip][node][d], tol ) );
    }
    return true;
}


bool IsoparametricLinearPyramid_Test::Volume_Test()
{
    const double tol = 1.0e-12;

    { // reference element — 5-point rule
        IsoparametricLinearPyramid pyra(5);
        Element<3U> elmt(&pyra);
        Node<3> n0,n1,n2,n3,n4;
        n0.Coordinate(Point<3>(-1.,-1.,0.)); n0.Idx(0);
        n1.Coordinate(Point<3>( 1.,-1.,0.)); n1.Idx(1);
        n2.Coordinate(Point<3>( 1., 1.,0.)); n2.Idx(2);
        n3.Coordinate(Point<3>(-1., 1.,0.)); n3.Idx(3);
        n4.Coordinate(Point<3>( 0., 0.,1.)); n4.Idx(4);
        elmt.Assign(0,&n0); elmt.Assign(1,&n1);
        elmt.Assign(2,&n2); elmt.Assign(3,&n3);
        elmt.Assign(4,&n4); elmt.CoordinateMatrix();
        _test( approximatelyEqualAbsRel( elmt.Volume(), 4./3., tol ) );
    }
    { // reference element — 8-point rule
        IsoparametricLinearPyramid pyra(8);
        Element<3U> elmt(&pyra);
        Node<3> n0,n1,n2,n3,n4;
        n0.Coordinate(Point<3>(-1.,-1.,0.)); n0.Idx(0);
        n1.Coordinate(Point<3>( 1.,-1.,0.)); n1.Idx(1);
        n2.Coordinate(Point<3>( 1., 1.,0.)); n2.Idx(2);
        n3.Coordinate(Point<3>(-1., 1.,0.)); n3.Idx(3);
        n4.Coordinate(Point<3>( 0., 0.,1.)); n4.Idx(4);
        elmt.Assign(0,&n0); elmt.Assign(1,&n1);
        elmt.Assign(2,&n2); elmt.Assign(3,&n3);
        elmt.Assign(4,&n4); elmt.CoordinateMatrix();
        _test( approximatelyEqualAbsRel( elmt.Volume(), 4./3., tol ) );
    }
    { // distorted element — 8-point rule
        IsoparametricLinearPyramid pyra(8);
        Element<3U> elmt(&pyra);
        Node<3> n0,n1,n2,n3,n4;
        n0.Coordinate(Point<3>( 0., 0., 0.)); n0.Idx(0);
        n1.Coordinate(Point<3>(13., 2., 0.)); n1.Idx(1);
        n2.Coordinate(Point<3>(14.,11., 0.)); n2.Idx(2);
        n3.Coordinate(Point<3>( 1., 9., 0.)); n3.Idx(3);
        n4.Coordinate(Point<3>(24., 2.,22.)); n4.Idx(4);
        elmt.Assign(0,&n0); elmt.Assign(1,&n1);
        elmt.Assign(2,&n2); elmt.Assign(3,&n3);
        elmt.Assign(4,&n4); elmt.CoordinateMatrix();
        _test( approximatelyEqualAbsRel( elmt.Volume(), 843.3333333333333, tol ) );
    }
    return true;
}


} // end csmp

