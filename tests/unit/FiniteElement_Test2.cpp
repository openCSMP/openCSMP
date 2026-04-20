#include "FiniteElement_Test2.h"
#include "Point.h"
#include "Node.h"
#include "Element.h"
#include "IsoparametricLinearPrism.h"
#include "compareFloats.h"
#include "meshManagementUtilities.h"

using namespace std;

namespace csmp {

/// tests fem operations against provided data in text file
void FiniteElement_Test2::run()
{
    constexpr uint32_t dim{3};
    
    // testing that the normals to isoparametric linear prim are outward pointing
    // --------------------------------------------------------------------------
    // creating prism in reference coordinates
    uint32_t integrationPoints{6};
    IsoparametricLinearPrism prism( integrationPoints );
    LocalVariables lvs;
    
    // reference element
    /*
     NXYZ(0,0) = 0.0; NXYZ(0,1) = 0.0; NXYZ(0,2) = -1.0;
     NXYZ(1,0) = 1.0; NXYZ(1,1) = 0.0; NXYZ(1,2) = -1.0;
     NXYZ(2,0) = 0.0; NXYZ(2,1) = 1.0; NXYZ(2,2) = -1.0;
     NXYZ(3,0) = 0.0; NXYZ(3,1) = 0.0; NXYZ(3,2) =  1.0;
     NXYZ(4,0) = 1.0; NXYZ(4,1) = 0.0; NXYZ(4,2) =  1.0;
     NXYZ(5,0) = 0.0; NXYZ(5,1) = 1.0; NXYZ(5,2) =  1.0;
    */
    Node<3> n0( 0, Point<dim>(0.,0.,-1.), lvs ), n1( 1, Point<dim>(1.,0.,-1.), lvs ), n2( 2, Point<dim>(0.,1.,-1.), lvs ),
            n3( 3, Point<dim>(0.,0.,1.), lvs ),  n4( 4, Point<dim>(1.,0.,1.), lvs ),  n5( 5, Point<dim>(0.,1.,1.), lvs );
    
    Element<dim> e0( &prism );
    e0.Idx(999);
    e0.Assign( 0, &n0 );
    e0.Assign( 1, &n1 );
    e0.Assign( 2, &n2 );
    e0.Assign( 3, &n3 );
    e0.Assign( 4, &n4 );
    e0.Assign( 5, &n5 );
    e0.CoordinateMatrix();
    DenseMatrix<DM_MIN> DATA(e0.Nodes(),1);
    for ( uint32_t i=0; i<e0.Nodes(); ++i ) DATA(i,0) = static_cast<double>(i);
    prism.OutputNodeDataToVTK( "isoprism", "nodevar", DATA );
    
    // testing normals to top and bottom faces
    vector<double> unrml;
    e0.FE()->UnitNormalToFace( 0, unrml );
    Point<3> bottom_nrml(0.,0.,-1), face0_nrml(unrml[0],unrml[1],unrml[2]);
    _test( approximatelyEqual( dotProduct(bottom_nrml,face0_nrml), 1.) == true );

    e0.FE()->UnitNormalToFace( 4, unrml );
    Point<3> top_nrml(0.,0.,1), face4_nrml(unrml[0],unrml[1],unrml[2]);
    _test( approximatelyEqual( dotProduct(top_nrml,face4_nrml), 1.) == true );

    // testing all faces
    _test( areUnitNormalsToFacesAreOutwardPointing(&e0) == true );
    
    cout <<"\n"<<"FiniteElement_Test2::run: complete."<< endl;

} // run



} // csmp
