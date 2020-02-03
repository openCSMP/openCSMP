#include "IsoparametricLinearTetrahedron.h"
//#include <climits>
#include "Exception.h"
#include "TriangularFacet.h"

using namespace std;

namespace csmp {

IsoparametricLinearTetrahedron::IsoparametricLinearTetrahedron	(
                                                    size_t integrationPoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_TETRAHEDRON, true, true, 1U ),
      NXYZ(4,3),
      IP(integrationPoints,3),
      projectionCalledNTimes(0),
      accDistance(0.0),
      totIterations(0),
      nonConvergenceOfProjections(0)
 {
    //AAM, 11.10.02
    dim = 3;
    itp = 1;
    npf = 3;
    npe = 4;
    fpe = 4;
    spe = 6;
    epe = 4;
    nne = 8;
    cne = 0;
    // Integration points
    gpe = integrationPoints;

    M.Resize(npe,npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_LINEAR_TETRAHEDRON);

      // base class matrices
    XY.Resize(npe,dim);
    JAC.Resize(dim,dim);
    JINV.Resize(dim,dim);

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    // initializing local node coordinates
    // local node coordinates are defined as r==ksi, s==nu, t==mu.
    // similar to the base of the tetrahedron
    NXYZ(0,0) = 0.0, NXYZ(0,1) = 0.0, NXYZ(0,2) =  0.0;
    NXYZ(1,0) = 1.0, NXYZ(1,1) = 0.0, NXYZ(1,2) =  0.0;
    NXYZ(2,0) = 0.0, NXYZ(2,1) = 1.0, NXYZ(2,2) =  0.0;
    NXYZ(3,0) = 0.0, NXYZ(3,1) = 0.0, NXYZ(3,2) =  1.0;

    W.resize(gpe);

    if(integrationPoints==1)
    {
        IP(0,0)=1./4.; IP(0,1)=1./4.; IP(0,2)=1./4.;
        W[0]=1./6.;
    }
    else if (integrationPoints==4)
    {
        IP(0,0)=0.585410196624969; IP(0,1)=0.138196601125011; IP(0,2)=0.138196601125011;
        W[0]=1./24.;
        IP(1,0)=0.138196601125011; IP(1,1)=0.585410196624969; IP(1,2)=0.138196601125011;
        W[1]=1./24.;
        IP(2,0)=0.138196601125011; IP(2,1)=0.138196601125011; IP(2,2)=0.585410196624969;
        W[2]=1./24.;
        IP(3,0)=0.138196601125011; IP(3,1)=0.138196601125011; IP(3,2)=0.138196601125011;
        W[3]=1./24.;
     }
    else {
        cerr<<"\nIsoparametricLinearTetrahedron::IsoparametricLinearTetrahedron: Number of integration points should be 1 or 4"<<endl;
        throw std::range_error
        ("***ERROR: IsoparametricLinearTetrahedron::IsoparametricLinearTetrahedron: N of integration points should be 1 or 4");
     }

}



IsoparametricLinearTetrahedron::~IsoparametricLinearTetrahedron()
 {
//   cout<<" IsoparametricLinearTetrahedron::~IsoparametricLinearTetrahedron: during lifetime of the element: "<<endl;
//   cout<<" Physical->Parametric projections called "<<projectionCalledNTimes<<" times, non-convergent="<<nonConvergenceOfProjections<<endl;
//   cout<<" Total number of N-R Physical->Parametric projection iterations ="<<totIterations<<endl;
//   cout<<" Accumulated projection error:                                  ="<<accDistance<<endl;
 }




/**

Returns the value of the element interpolation functions at the point 'rst'
in local coordinates.

@section arguments Input Arguments

The floating point coordinates 'r', 's' and 't' of the point at which the
interpolation shall be carried out.

@param N The fourth method argument is the vector into which the values of the
n=nodes interpolation functions at the point 'rst' will be returned.

@section implementation Implementation

See in source code header file.

@section application Application

Method is used to compute property values at the integration points of
the element.

*/
void IsoparametricLinearTetrahedron::Nrst(
                    double64 r,
                    double64 s,
                    double64 t,
                    std::vector<double64>& N ) const
{
    N.resize(npe);
    N[0] = 1. - r - s - t; // L1
    N[1] = r; // L2
    N[2] = s; // L3
    N[3] = t; // L4
 }

void IsoparametricLinearTetrahedron::Nrst(
                    double64 r,
                    double64 s,
                    double64 t,
                    double64* N ) const
{
    N[0] = 1. - r - s - t; // L1
    N[1] = r; // L2
    N[2] = s; // L3
    N[3] = t; // L4
 }

/**

This method and the complementary method dNs() compute the shape function
derivates with respect to the local coordinate axis 'r', and 's'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNR The shape function derivatives are returned into the second method
argument which is a floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.

*/
void IsoparametricLinearTetrahedron::dNr (
                double64,
                double64,
                double64,
                std::vector<double64>& DNR ) const
{
   DNR.resize(npe);
   DNR[0] = -1.;
   DNR[1] =  1.;
   DNR[2] =  0.;
   DNR[3] =  0.;
 }


void IsoparametricLinearTetrahedron::dNs(
                double64,
                double64,
                double64,
                std::vector<double64>& DNS ) const
{
   DNS.resize(npe);
   DNS[0] = -1.0;
   DNS[1] =  0.0;
   DNS[2] =  1.0;
   DNS[3] =  0.0;
}



void IsoparametricLinearTetrahedron::dNt(
                double64,
                double64,
                double64,
                std::vector<double64>& DNT ) const
{
   DNT.resize(npe);
   DNT[0] = -1.;
   DNT[1] =  0.;
   DNT[2] =  0.;
   DNT[3] =  1.;
}




 void
 IsoparametricLinearTetrahedron::JacobianAtIntegrationPoint( size_t gauss_point )
 {
      assert( gauss_point < gpe );

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }




void
IsoparametricLinearTetrahedron::N_AtBaryCenter( std::vector<double64>& N )
 {
    N.resize(npe);
    const double64 OneFourth=0.25;
    Nrst(OneFourth,OneFourth,OneFourth,N);
 }






/** Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local corner node ID numbers
for the element.
*/
void
IsoparametricLinearTetrahedron::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
}



/** Returns the counter-clockwise local node numbering for the element.

@param ids Returns the counter-clockwise local node numbering for the element.

*/
void
IsoparametricLinearTetrahedron::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
}



void
IsoparametricLinearTetrahedron::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    snids.resize(2);
    if ( segm_id == 0 ) {
         snids[0] = 0;
         snids[1] = 1;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 1;
         snids[1] = 2;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 2;
         snids[1] = 0;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 0;
         snids[1] = 3;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 1;
         snids[1] = 3;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 2;
         snids[1] = 3;
      }
    else
    std::cout <<"\nIsoparametricLinearTetrahedron::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment



/// local node ids in counter-clockwise order from the outside looking into the faces
void IsoparametricLinearTetrahedron::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    fnids.resize(3);
    if  ( face_id == 0  )
      {
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 3;
      }
    else if ( face_id == 1 )
      {
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 2;
     }
    else if ( face_id == 2 )
      {
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 3;
     }
    else if ( face_id == 3 )
      {
         fnids[0] = 0;
         fnids[1] = 2;
         fnids[2] = 1;
      }
    else
    std::cout <<"\nIsoparametricLinearTetrahedron::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }



CSMP_FEM_TYPE  IsoparametricLinearTetrahedron::ElementTypeOfFace( size_t )  const
 {
    return ISOPARAMETRIC_LINEAR_TRIANGLE;
 }

double64 IsoparametricLinearTetrahedron::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }




/** Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN4 The interpolation-function derivative matrix is returned into the
second method argument.
*/
void
IsoparametricLinearTetrahedron::dN( DenseMatrix<DM_MIN>& DN4 )
{
    DN4.Resize(dim,npe);
    M.Resize(dim,1);

     // Jacobian transformation to global coordinate system
     for ( size_t i=0; i<npe; i++ )
       {
          // here the global coordinates come in
          dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
          dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
          dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

          // compute Jacobian matrix, its determinant and inversex
          Jacobian( DNR, DNS, DNT );
          JacobianInverse();
          ////////// Debug Printout////////////////////////////////
          //  cout<<" At vertice #: "<<i<<endl;
          //	cout<<" Jacobian Matrix: "<<endl;
          //	JAC.Out();
          //	cout<<" Jacobian Inverse Matrix: "<<endl;
          //JINV.Out();
          /////////////////////////////////////////////////////////

          M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];

          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN4(0,i) = JINV(0,0);
          DN4(1,i) = JINV(1,0);
          DN4(2,i) = JINV(2,0);
          JINV.Resize(dim,dim);
       }
} // end dN




/** Computes derivatives of shape functions at a point XY given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.

*/
double64
IsoparametricLinearTetrahedron::dN( DenseMatrix<DM_MIN>& DN2,
                                    const vector<double64>& xyz )
 {
    vector<double64> rst(dim);

    PhysicalToParametric(rst, xyz);

    // here find derivatives
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    DN2.Resize(dim,dim);
    DN2  = JINV;

    //Instead of using previos version: dN(DN);
    DenseMatrix<DM_MIN>DN(dim,npe);
    for(int i=0;i<npe;i++){
        DN(0,i)=DNR[i];
        DN(1,i)=DNS[i];
        DN(2,i)=DNT[i];
    }

    DN2*=DN;

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
//    cout<<" IsoparametricLinearTetrahedron::dN  For given xyz=("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<"), rst=("<<
//    rst[0]<<","<<rst[1]<<","<<rst[2]<<")"<<endl;
//    cout<<" IsoparametricLinearTetrahedron::dN  Matrix DN2: "<<endl;
//    DN2.Out();
//    getchar();
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

   return detJ;
 }





/** Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
*/
void
IsoparametricLinearTetrahedron::N( vector<double64>& N,const vector<double64>& xyz )
{

      DNR.resize(dim);

      PhysicalToParametric(DNR, xyz);

      Nrst(DNR[0],DNR[1],DNR[2], N );

}



/**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J at the desired node
point.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the
point i (i=0...i=gauss points-1).

@section arguments Input Arguments

The first argument is a reference to the parent Element object,
the second argument is the result matrix, the third argument indicates
the gauss point at which the derivative matrix is computed and the
fourth argument specifies the degrees of freedom per node, for which
the B matrix shall be transformed.

@return The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.

*/
double64  IsoparametricLinearTetrahedron::dN_AtNode( DenseMatrix<DM_MIN>& B,
                                                      size_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 6 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];

    B = JINV * B;

    return detJ;
 }




/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricLinearTetrahedron::IntegrationPoint( size_t ip,
                                                        vector<double64>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(3U); xyz[0]=xyz[1]=xyz[2]=0.;
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for( size_t i=0U; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint





/** Projection function from rst->xyz
*/
void IsoparametricLinearTetrahedron::ParametricToPhysical( vector<double64>& rst,
                                                           vector<double64>& xyz )
{
    Nrst(rst[0],rst[1],rst[2], DNR );

    xyz.resize(dim);
    xyz[0]=xyz[1]=xyz[2]=0.;

    for(size_t i=0; i<npe; i++) {
         xyz[0]+=XY(i,0)*DNR[i];
         xyz[1]+=XY(i,1)*DNR[i];
         xyz[2]+=XY(i,2)*DNR[i];
      }
}



/** Projection function from xyz->rst
*/
void IsoparametricLinearTetrahedron::PhysicalToParametric(
                                       vector<double64>& rSt,
                                       const vector<double64>& xyz )
{

    const double64 x12( XY(0,0) - XY(1,0) );
    const double64 x13( XY(0,0) - XY(2,0) );
    const double64 x14( XY(0,0) - XY(3,0) );
    const double64 x23( XY(1,0) - XY(2,0) );
    const double64 x24( XY(1,0) - XY(3,0) );
    const double64 x34( XY(2,0) - XY(3,0) );
    const double64 x21( -x12 );
    const double64 x31( -x13 );
    const double64 x32( -x23 );
    //const double64 x42( -x24 );
    const double64 x43( -x34 );

    const double64 y12( XY(0,1) - XY(1,1) );
    const double64 y13( XY(0,1) - XY(2,1) );
    const double64 y14( XY(0,1) - XY(3,1) );
    const double64 y23( XY(1,1) - XY(2,1) );
    const double64 y24( XY(1,1) - XY(3,1) );
    const double64 y34( XY(2,1) - XY(3,1) );
    const double64 y21( -y12 );
    const double64 y31( -y13 );
    //const double64 y32( -y23 );
    //const double64 y42( -y24 );
    const double64 y43( -y34 );

    const double64 z12( XY(0,2) - XY(1,2) );
    const double64 z13( XY(0,2) - XY(2,2) );
    const double64 z14( XY(0,2) - XY(3,2) );
    const double64 z23( XY(1,2) - XY(2,2) );
    const double64 z24( XY(1,2) - XY(3,2) );
    const double64 z34( XY(2,2) - XY(3,2) );
    const double64 z21( -z12 );
    const double64 z31( -z13 );
    //const double64 z32( -z23 );
    //const double64 z42( -z24 );
    const double64 z43( -z34 );

    //const double64 a1( y42*z32 - y32*z42 );
    const double64 a2( y31*z43 - y34*z13 );
    const double64 a3( y24*z14 - y14*z24 );
    const double64 a4( y13*z21 - y12*z31 );

    //const double64 b1( x32*z42 - x42*z32 );
    const double64 b2( x43*z31 - x13*z34 );
    const double64 b3( x14*z24 - x24*z14 );
    const double64 b4( x21*z13 - x31*z12 );

    //const double64 c1( x42*y32 - x32*y42 );
    const double64 c2( x31*y43 - x34*y13 );
    const double64 c3( x24*y14 - x14*y24 );
    const double64 c4( x13*y21 - x12*y31 );

    const double64 V00( x21*( y23*z34 - y34*z23) + x32*(y34*z12 - y12*z34) + x43*(y12*z23 - y23*z12) );

    //const double64 V01( XY(1,0) * ( XY(2,1)*XY(3,2) - XY(3,1)*XY(2,2) ) + XY(2,0) * ( XY(3,1)*XY(1,2) - XY(1,1)*XY(3,2) ) + XY(3,0) * ( XY(1,1)*XY(2,2) - XY(2,1)*XY(1,2) ) );
    const double64 V02( XY(0,0) * ( XY(3,1)*XY(2,2) - XY(2,1)*XY(3,2) ) + XY(2,0) * ( XY(0,1)*XY(3,2) - XY(3,1)*XY(0,2) ) + XY(3,0) * ( XY(2,1)*XY(0,2) - XY(0,1)*XY(2,2) ) );
    const double64 V03( XY(0,0) * ( XY(1,1)*XY(3,2) - XY(3,1)*XY(1,2) ) + XY(1,0) * ( XY(3,1)*XY(0,2) - XY(0,1)*XY(3,2) ) + XY(3,0) * ( XY(0,1)*XY(1,2) - XY(1,1)*XY(0,2) ) );
    const double64 V04( XY(0,0) * ( XY(2,1)*XY(1,2) - XY(1,1)*XY(2,2) ) + XY(1,0) * ( XY(0,1)*XY(2,2) - XY(2,1)*XY(0,2) ) + XY(2,0) * ( XY(1,1)*XY(0,2) - XY(0,1)*XY(1,2) ) );

    //N[0] = ( V01 + a1*xyz[0] + b1*xyz[1] + c1*xyz[2])/V00;
    //N[1] = ( V02 + a2*xyz[0] + b2*xyz[1] + c2*xyz[2])/V00;
    //N[2] = ( V03 + a3*xyz[0] + b3*xyz[1] + c3*xyz[2])/V00;
    //N[3] = ( V04 + a4*xyz[0] + b4*xyz[1] + c4*xyz[2])/V00;

    rSt[0] = ( V02 + a2*xyz[0] + b2*xyz[1] + c2*xyz[2])/V00;
    rSt[1] = ( V03 + a3*xyz[0] + b3*xyz[1] + c3*xyz[2])/V00;
    rSt[2] = ( V04 + a4*xyz[0] + b4*xyz[1] + c4*xyz[2])/V00;

}

// tested: o.k.
inline size_t IsoparametricLinearTetrahedron::n( size_t i, size_t a ) const
{
     if ( i+a >= 4 ) return i+a-4;
     return i+a;
}


/**

Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64  IsoparametricLinearTetrahedron::AspectRatio()
{
   NRST.resize(spe);

   EdgeLengths( NRST );

   double64 seg_max(NRST[0]), seg_min(NRST[0]);

   // find largest segment
   for ( size_t i=1; i<spe; i++ ) {
        if ( NRST[i] > seg_max ) seg_max = NRST[i];
        if ( NRST[i] < seg_min ) seg_min = NRST[i];
     }

   return seg_max / seg_min;
}




/** Segment length is calculated directly as the length of linear segments.

@param len The lengths of the 12 segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as simple Euclidian distance between corner
vertices of the edge.
*/
void IsoparametricLinearTetrahedron::EdgeLengths( std::vector<double64>& len )
{
    double64 sum;
    len.resize(spe);

    // segment 1
    sum     = (XY(1,0)-XY(0,0)) * (XY(1,0)-XY(0,0));
    sum    += (XY(1,1)-XY(0,1)) * (XY(1,1)-XY(0,1));
    sum    += (XY(1,2)-XY(0,2)) * (XY(1,2)-XY(0,2));
    len[0]  = sqrt(sum);
    // segment 2
    sum     = (XY(2,0)-XY(1,0)) * (XY(2,0)-XY(1,0));
    sum    += (XY(2,1)-XY(1,1)) * (XY(2,1)-XY(1,1));
    sum    += (XY(2,2)-XY(1,2)) * (XY(2,2)-XY(1,2));
    len[1]  = sqrt(sum);
    // segment 3
    sum     = (XY(0,0)-XY(2,0)) * (XY(0,0)-XY(2,0));
    sum    += (XY(0,1)-XY(2,1)) * (XY(0,1)-XY(2,1));
    sum    += (XY(0,2)-XY(2,2)) * (XY(0,2)-XY(2,2));
    len[2]  = sqrt(sum);
    // segment 4
    sum     = (XY(3,0)-XY(0,0)) * (XY(3,0)-XY(0,0));
    sum    += (XY(3,1)-XY(0,1)) * (XY(3,1)-XY(0,1));
    sum    += (XY(3,2)-XY(0,2)) * (XY(3,2)-XY(0,2));
    len[3]  = sqrt(sum);
    // segment 5
    sum     = (XY(3,0)-XY(1,0)) * (XY(3,0)-XY(1,0));
    sum    += (XY(3,1)-XY(1,1)) * (XY(3,1)-XY(1,1));
    sum    += (XY(3,2)-XY(1,2)) * (XY(3,2)-XY(1,2));
    len[4]  = sqrt(sum);
    // segment 6
    sum     = (XY(3,0)-XY(2,0)) * (XY(3,0)-XY(2,0));
    sum    += (XY(3,1)-XY(2,1)) * (XY(3,1)-XY(2,1));
    sum    += (XY(3,2)-XY(2,2)) * (XY(3,2)-XY(2,2));
    len[5]  = sqrt(sum);
 } // end EdgeLengths



/**

Computes Integral N dA = sum_1...n Wi Ji Ni. This formulation also
gives a meaningful area if the element boundaries are curved. For a
description of numerical integration  of isoparametric quadratic
triangular elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The area (m2) of the finite element.
*/
double64 IsoparametricLinearTetrahedron::Volume()
{
    double64   area; // determinant
    size_t  i;

    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( area=0.0, i=0; i<gpe; i++ )
      {
         dNr( IP(i,0), IP(i,1), IP(i,2), DNR );
         dNs( IP(i,0), IP(i,1), IP(i,2), DNS );
         dNt( IP(i,0), IP(i,1), IP(i,2), DNT );

         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix
         Jacobian( DNR, DNS, DNT );
         area += JacobianInverse() * W[i];
      }

   //	cout<<" Volume of Tetra"<<area<<endl;
   return area;
}




/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

@param N The interpolation function values are returned into the third argument.

*/
void IsoparametricLinearTetrahedron::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )
 {
    assert( ip < gpe );
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), N );
 } // end N_AtIntegrationPoint







 /**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the gauss
point i (i=0...i=gauss points-1).

@section arguments Input Arguments

The first argument is a reference to the parent Element object,
the second argument is the result matrix, the third argument indicates
the gauss point at which the derivative matrix is computed and the
fourth argument specifies the degrees of freedom per node, for which
the B matrix shall be transformed.

@return The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.

*/
double64
IsoparametricLinearTetrahedron::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, size_t gauss_point )
 {
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    // cout<<" IsoparametricLinearTetrahedron::dN_AtIntegrationPoint: No integration points used ..."<<endl;
    assert( gauss_point < gpe );
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 8 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];

    B = JINV * B;

    return detJ;
 }


double64
IsoparametricLinearTetrahedron::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    double64 OneFourth=0.25;

    dNr( OneFourth,OneFourth,OneFourth, DNR );
    dNs( OneFourth,OneFourth,OneFourth, DNS );
    dNt( OneFourth,OneFourth,OneFourth, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];

    B = JINV * B;

    return detJ;
 }


 /**

This method attempts to calculate the inner radius of the curved-sided
quad by talking 1/6 of the perimeter and dividing this measure by
the area of the quad.

The procedure is empirically based giving a good match if the quadss
are relatively even-sided and have straight edges.

@section arguments Input Arguments

The parent element is queried for its node coordinates.
*/
double64
IsoparametricLinearTetrahedron::InnerRadius()
{
   NRST.resize(spe);
   double64  sum(0.);

   EdgeLengths( NRST );
   for ( size_t i=0; i<spe; i++ ) sum += NRST[i];

   if(AspectRatio()>4.)
   cout<<" IsoparametricLinearTetrahedron::InnerRadius: ***WARNING: function not applicable for CURRENT HAR element"<<endl;

   return Volume() / (sum/2.);
}




/** Returns 0 local node ids of the nodes located at the midsides of
the element.

Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

@attention method is not implemented.

*/
void IsoparametricLinearTetrahedron::MidSideNodes(std::vector<size_t>&) const
 {
    throw csmp::Exception( WARNING, "IsoparametricLinearTetrahedron::MidSideNodes:",
                                    "MidSideNodes not present",
                                     "Probably unintended use of function." );
 }


void
IsoparametricLinearTetrahedron::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                            vector<size_t>& fnids )
 {
    fnids.resize(bnodes.size());
    cout <<"\nIsoparametricLinearTetrahedron::ConsecutiveNodesAtBoundary: not implemented."<< endl;
     if ( bnodes.size() != 4 && bnodes.size() != 3 )
       throw csmp::Exception( ERROR, "IsoparametricLinearTetrahedron::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );


 } // end ConsecutiveNodesAtBoundary


 /**

As the linear character of the element do not expect high accuracy, one
integration point is enough and used in most of the applications (FE, V1
by Zienkiwicz).

Because in the local coordinate system, the pyramid is always straight
sided, and the derivatives of the interpolation functions are bilinear,
variable values computed at integration points can simply be extrapolated
to the nodes.  This involves the steps:

1. Compute the interpolation function coefficients for the pyramid which
   is defined by the one integration point.

2. For each node point compute the values of the interpolation functions
   and use these to extrapolate the values of the variables at the nodes.
*/
void
IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                            const vector<double64>& IVAR,
                                                                            vector<double64>& NVAR	)
const
{
   assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

   if ( gpe == 1U ) {
        // Define nodal values as bi-linear variation of the integration points values
        // See Zienkewitch, pp. 351, for example
        for ( size_t i=0; i<npe; i++ )
            for ( size_t k=0; k<nvars; k++ )	NVAR[i*nvars + k] = IVAR[k];

        return;
     }

   static double64  a[4], b[4], c[4], d[4], intpol[4], volume6;
   static bool      first_call(true);
   size_t i;
   int32     j;

   vector<double64>  sum(nvars);

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   if ( IVAR.size() != (gpe*nvars) )
     throw csmp::Exception( FATAL_ERROR, "IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
                                  "Input vector must have 'nvars' x 4 entries");

   if ( NVAR.size() != (npe*nvars) )
     throw csmp::Exception( FATAL_ERROR, "IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
                                  "Output vector must have 'nvars' x nodes entries");

   // 1. Compute the local interpolation function coefficients for the tetrahedron which
   //    is defined by the four integration points.
   if ( first_call ) {
       // checking starting conditions
       if ( gpe != 4 )
       throw csmp::Exception( FATAL_ERROR, "IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
         "This method expects four integration points on which extrapolation functions will be based on" );

        // test function coefficients
        for ( i=0, j=1; i<4; i++ )
          {
             // a(i)
             a[i]  = -IP(n(i,1),0) * (IP(n(i,3),1)*IP(n(i,2),2)-IP(n(i,3),2)*IP(n(i,2),1));
             a[i] -=  IP(n(i,2),0) * (IP(n(i,1),1)*IP(n(i,3),2)-IP(n(i,1),2)*IP(n(i,3),1));
             a[i] -=  IP(n(i,3),0) * (IP(n(i,2),1)*IP(n(i,1),2)-IP(n(i,2),2)*IP(n(i,1),1));
             // b(i)
             b[i]  = IP(n(i,3),1)*IP(n(i,2),2) - IP(n(i,3),2)*IP(n(i,2),1);
             b[i] += IP(n(i,1),1)*IP(n(i,3),2) - IP(n(i,1),2)*IP(n(i,3),1);
             b[i] += IP(n(i,2),1)*IP(n(i,1),2) - IP(n(i,2),2)*IP(n(i,1),1);
             // c(i)
             c[i]  = IP(n(i,3),2)*IP(n(i,2),0) - IP(n(i,3),0)*IP(n(i,2),2);
             c[i] += IP(n(i,1),2)*IP(n(i,3),0) - IP(n(i,1),0)*IP(n(i,3),2);
             c[i] += IP(n(i,2),2)*IP(n(i,1),0) - IP(n(i,2),0)*IP(n(i,1),2);
             // d(i)
             d[i]  = IP(n(i,3),0)*IP(n(i,2),1) - IP(n(i,3),1)*IP(n(i,2),0);
             d[i] += IP(n(i,1),0)*IP(n(i,3),1) - IP(n(i,1),1)*IP(n(i,3),0);
             d[i] += IP(n(i,2),0)*IP(n(i,1),1) - IP(n(i,2),1)*IP(n(i,1),0);
             //
             a[i] *=  static_cast<double64>(j);
             b[i] *=  static_cast<double64>(j);
             c[i] *=  static_cast<double64>(j);
             d[i] *=  static_cast<double64>(j);
             j    *= -1;
          }
        // computing local element volume x 6
        volume6 = a[0] + a[1] + a[2] + a[3];

        first_call = false;
     }

   // 2. For each node point compute the values of the linear extrapolation functions
   //    and use these to extrapolate the values of the variables at the nodes.
   for ( i=0; i<npe; i++ )
     {
        // compute interpolation function values at node i
        intpol[0] = (a[0] + b[0] * NXYZ(i,0) + c[0] * NXYZ(i,1) + d[0] * NXYZ(i,2)) / volume6;
        intpol[1] = (a[1] + b[1] * NXYZ(i,0) + c[1] * NXYZ(i,1) + d[1] * NXYZ(i,2)) / volume6;
        intpol[2] = (a[2] + b[2] * NXYZ(i,0) + c[2] * NXYZ(i,1) + d[2] * NXYZ(i,2)) / volume6;
        intpol[3] = (a[3] + b[3] * NXYZ(i,0) + c[3] * NXYZ(i,1) + d[3] * NXYZ(i,2)) / volume6;

        // carry out extrapolation
        fill( sum.begin(), sum.end(), static_cast<double64>(0.0) );
        for ( size_t j=0; j<gpe; j++ )
          for ( size_t k=0; k<nvars; k++ ) sum[k] += intpol[j] * IVAR[j*nvars + k];

        // store result in output vector
        for ( size_t k=0; k<nvars; k++ ) NVAR[i*nvars + k] = sum[k];
     }

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)




// Adriana 2005
void IsoparametricLinearTetrahedron::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords = NXYZ;
}


// Adriana 2005
void IsoparametricLinearTetrahedron::JacobianAt( const std::vector<double64>& rst )
{
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR, DNS, DNT );
}



/**
     uses the TriangularFacet to compute the normals to its faces.

     @author SKM 15/2/2016
     
     @test OK 
*/
void  IsoparametricLinearTetrahedron::UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const
 {
     assert( face < Faces() );
     unrml.resize(3);
   
     // if face lies opposite to node 0
     if ( face == 0 ) { // OK
          Point<3> nrml = normalOfTriangle( Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                            Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                            Point<3>(XY(3,0),XY(3,1),XY(3,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 1 ) { // OK - counter-clockwise nodes (ouside looking in): 0-3-2
          Point<3> nrml = normalOfTriangle( Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                            Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                            Point<3>(XY(0,0),XY(0,1),XY(0,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 2 ) { // OK
          Point<3> nrml = normalOfTriangle( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                            Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                            Point<3>(XY(3,0),XY(3,1),XY(3,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 3 ) { // OK
          Point<3> nrml = normalOfTriangle( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                            Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                            Point<3>(XY(1,0),XY(1,1),XY(1,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }
   
 } // end UnitNormalToFace





/**
     the matrix DATA (dimensions x nodes) can either contain scalars (1-row)
     vectors (dim rows). The method may be extended to contain tensor data.
*/
void IsoparametricLinearTetrahedron::OutputNodeDataToVTK( const char* file_name,
                                                          const char* var_name,
                                                          DenseMatrix<DM_MIN>& DATA ) const
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     sprintf( elmt, "%lu", CurrentID() ); // ID may not be initialised
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cout <<"\nIsoparametricLinearTetrahedron::OutputNodeDataToVTK ";
           cout <<"Output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" double"<< endl;
     for ( size_t i=0; i<npe; i++ )
       {
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< 5 << endl;

     ofs << 4 <<" 0 1 2 3" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 10 << endl; // VTK_TETRA
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     // Unfortunately the data can only be output as nodal variables
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );

     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" double"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x9
           for ( size_t i=0; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" double"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 10
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
               for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricLinearTetrahedron::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

 } // end OutputNodeDataToVTK





} // end namespace csmp
