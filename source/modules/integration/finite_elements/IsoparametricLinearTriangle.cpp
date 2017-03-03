#include "IsoparametricLinearTriangle.h"
//#include <set>
//#include <fstream>
//#include <cstring>
#include <climits>
#include "Exception.h"
#include "MJL_Edge.h"

using namespace std;

namespace csmp {

/**

Default constructor initializes quadrilateral to use local coordinates
and a four-point integration scheme where the Gauss points are located
inside quad at convergence points . The weights of the integration points
are set to 1 such that the integrated contributions add up to 1.  .
*/
IsoparametricLinearTriangle::IsoparametricLinearTriangle( size_t dimensions,
                                                          size_t ipoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_TRIANGLE, true, true, 1U ),
   use2Dto3Djacobi( !(dimensions==2) ),
   DN(dimensions,3),
   NXY(3,dimensions),
   JMAT(dimensions,dimensions),
   RS(2),
   LXY(3)
{    
    dim = dimensions;
    itp = 1;
    npf = 2;
    npe = 3;
    fpe = 3;
    spe = 3;
    epe = 3;
    nne = 6;
    cne = 6;
    gpe = ipoints;

    // base class matrices
    XY.Resize(npe,dim);
    JAC.Resize(dim,dim);
    JINV.Resize(dim,dim);
    JMAT.Resize(dim,dim);

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    W.resize( gpe );
    rr.resize( gpe );
    ss.resize( gpe );

    // full 4-point integration basis
    if( ipoints == 1 )
    {
        rr[0] = 1./3.;
        ss[0] = 1./3.;
        W[0]  = .5;
    }
    else if( ipoints == 3 )
    {
        rr[0] = 0.5,    rr[1] = 0.5,   rr[2] = 0.0;
        ss[0] = 0.0,    ss[1] = 0.5,   ss[2] = 0.5;
        W[0]  = 1./6.,  W[1]  = 1./6., W[2]  = 1./6.;
    }
    else if( ipoints == 4 )
    {
        rr[0] = 1./3.,  rr[1] = 0.2,     rr[2] = 0.6,     rr[3] = 0.2;
        ss[0] = 1./3.,  ss[1] = 0.2,     ss[2] = 0.2,     ss[3] = 0.6;
        W[0]  = -27./96., W[1] = 25./96.,  W[2] = 25./96., W[3]  = 25./96.;
    }
    else
    {
        cout<<" IsoparametricLinearTriangle::IsoparametricLinearTriangle Number of integration points not supported: "<<ipoints<<endl;
        throw std::range_error
        ("***ERROR: IsoparametricLinearTriangle::IsoparametricLinearTriangle N of integration points should be 1,3 or 4");
    }


    // local node coordinates are defined as r==ksi, s==nu, t==mu.
    NXY(0,0) = 0.0;NXY(0,1) = 0.0;
    NXY(1,0) = 1.0;NXY(1,1) = 0.0;
    NXY(2,0) = 0.0;NXY(2,1) = 1.0;

    UsesLocalCoordinates(true);
    Isoparametric(true);
    SurfaceElement();
    ElementType(ISOPARAMETRIC_LINEAR_TRIANGLE);
}

IsoparametricLinearTriangle::~IsoparametricLinearTriangle()
 {
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
void
IsoparametricLinearTriangle::Nrs(
                    double64 r,
                    double64 s,
                    std::vector<double64>& N ) const
{
   N.resize(npe);
   N[0] = 1-r-s;
   N[1] = r;
   N[2] = s;
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
void
IsoparametricLinearTriangle::dNr (
                double64,
                double64,
                std::vector<double64>& DNR ) const
{
   DNR.resize(npe);
   DNR[0] = -1.;
   DNR[1] =  1.;
   DNR[2] =  0.;
}


void IsoparametricLinearTriangle::dNs(
                double64 ,
                double64 ,
                std::vector<double64>& DNS ) const
{
   DNS.resize(npe);
   DNS[0] = -1.;
   DNS[1] =  0.;
   DNS[2] =  1.;
}


/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

@param N The interpolation function values are returned into the third argument.

*/
void
IsoparametricLinearTriangle::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )
 {
    assert( ip < gpe );

    // local interpolation function values
    Nrs( rr[ip], ss[ip], N );

 } // end N_AtIntegrationPoint




/// segments are numbered like faces
void
IsoparametricLinearTriangle::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
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
    else
    std::cout <<"IsoparametricLinearTriangle::NodesOfSegment: Segment ID out of range: "<< segm_id << std::endl;
 }



/**

Gives the local node numbers (0...5) of the nodes which make up the
requested face of the element (possible numbers 0, 1, 2).
For this element, the faces are numbered such that face 0 lies opposite of
node 0, face 1 node 1 etc.

@section arguments Input Arguments

The number of the desired face (0, 1, or 2).

@param fnids The local node numbers are returned into the integer vector 'fnids'.

@section application Application

When boundary conditions shall be applied it is necessary to determine
the properties associated with the nodes of that element face.
*/
void
IsoparametricLinearTriangle::NodesOfFace( size_t face_id,
                                          std::vector<size_t>& fnids ) const
 {
    fnids.resize(2);

    if( face_id == 0 )
      {
         fnids[0] = 1;
         fnids[1] = 2;
      }
    else if ( face_id == 1 )
      {
         fnids[0] = 2;
         fnids[1] = 0;
      }
    else if ( face_id == 2 )
      {
         fnids[0] = 0;
         fnids[1] = 1;
      }
    else
    std::cerr <<"\nIsoparametricLinearTriangle::NodesOfFace: Erratic input face ID: "<< face_id << std::endl;
 }

/** Returns local node ids of the nodes which sit on the corner of
the quadratic triagular element.

@param ids Returns a Meschpp integer vector with the 3 local corner node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).
*/
void
IsoparametricLinearTriangle::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }




/** Returns local node ids of the nodes which sit on the midsides of
the quadratic triagular element.

Returns an integer vector with the 3 local midside-node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).

@todo (3) Check what should be assigned
*/
void  IsoparametricLinearTriangle::MidSideNodes( std::vector<size_t>& ) const
 {
     std::cout <<"\nIsoparametricLinearTriangle::MidSideNodes Not present..."<<std::endl;
 }




void  IsoparametricLinearTriangle::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }


CSMP_FEM_TYPE IsoparametricLinearTriangle::ElementTypeOfFace( size_t )  const
 {
    return ISOPARAMETRIC_LINEAR_BAR;
 }


double64 IsoparametricLinearTriangle::WeightAtIntegrationPoint( size_t i )
const { return W[i]; }


void  IsoparametricLinearTriangle::N_AtBaryCenter( std::vector<double64>& N )
 {
    N.resize(npe);
    Nrs(1./3.,1./3.,N);
 }



void IsoparametricLinearTriangle::Dimensions( size_t dimensions )
 {
   dim = dimensions;
   if ( dim == 3 ) use2Dto3Djacobi = true;
   else            use2Dto3Djacobi = false;
   // base class matrices
   XY.Resize(npe,dim);
   JAC.Resize(dim,dim);
   JINV.Resize(dim,dim);
   JMAT.Resize(dim,dim);
 }



/** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64
IsoparametricLinearTriangle::AspectRatio()
{
   vector<double64> vec(spe);

   EdgeLengths( vec );

   // order segment
   set<double64> segms;

   for ( size_t i=0; i<spe; i++ ) segms.insert( vec[i] );

   double64 segm1 = (*segms.begin()),
             segm2 = (*segms.rbegin());

   return segm2 / segm1;
}




/**

This method attempts to calculate the inner radius of the straight-sided
triangle by talking 1/2 of the perimeter and dividing this measure by
the area of the tri.

The procedure is empirically based giving a good match if triangles
are relatively even-sided and have straight edges.

@section arguments Input Arguments

The parent element is queried for its node coordinates.
*/
double64  IsoparametricLinearTriangle::InnerRadius()
{
   // Alternative radius
/*
   vector<double64> rsCenter(parametricDimehsions), rsBaseR(parametricDimehsions);
   rsCenter[0]=1-0.5*sqrt(2.0); rsCenter[1]=1-0.5*sqrt(2.0);
   rsBaseR [0]=1-0.5*sqrt(2.0); rsBaseR[1]=0.0;

   vector<double64> xyzCenter(dim),xyzBaseR(dim);

   ParametricToPhysical(rsCenter,xyzCenter);
   ParametricToPhysical(rsBaseR,xyzBaseR);

   double64 radius1=0.;

   if(dim==2){

   radius1=sqrt( (xyzCenter[0]-xyzBaseR[0])*(xyzCenter[0]-xyzBaseR[0]) +
                 (xyzCenter[1]-xyzBaseR[1])*(xyzCenter[1]-xyzBaseR[1]) );
   }
   else
    if(dim==3){
    radius1=sqrt( (xyzCenter[0]-xyzBaseR[0])*(xyzCenter[0]-xyzBaseR[0]) +
                  (xyzCenter[1]-xyzBaseR[1])*(xyzCenter[1]-xyzBaseR[1]) +
                  (xyzCenter[2]-xyzBaseR[2])*(xyzCenter[2]-xyzBaseR[2])
                );
    }

    //cout<<" IsoparametricLinearTriangle::InnerRadius r1="<<vol<<", r2="<<radius1<<endl;
    //getchar();

   if( AspectRatio() > 10. )
     throw csmp::Exception( CSMP_WARNING, "IsoparametricLinearTriangle::InnerRadius",
                    "Approximate radius - not valid for given large aspect-ratio element");

   return radius1;
  */

    NRST.resize(spe);
    double64  sum(0.);

    EdgeLengths( NRST );
    for ( size_t i=0; i<spe; i++ ) sum += NRST[i];

    if(AspectRatio()>4.)
    cout<<" IsoparametricLinearTetrahedron::InnerRadius: ***WARNING: function not applicable for CURRENT HAR element"<<endl;

    return Volume() / (sum/2.);


}


/**

Segment length is calculated as the sum of the length of two linear
segments making up each face of the triangle.

@param len The lengths of the segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as the sum of the straight line
segments which make up the element face.

@section application Application

For example, when stresses are to be applied at the element side, the
(area=length at unit thickness) must be taken into account.
*/
void
IsoparametricLinearTriangle::EdgeLengths( std::vector<double64>& len )
{
    double64 sum=0.0;
    len.resize(spe);

    if(dim==3)
    {
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
    }
    else
    if(dim==2)
    {
    // segment 1
    sum     = (XY(1,0)-XY(0,0)) * (XY(1,0)-XY(0,0));
    sum    += (XY(1,1)-XY(0,1)) * (XY(1,1)-XY(0,1));
    len[0]  = sqrt(sum);
    // segment 2
    sum     = (XY(2,0)-XY(1,0)) * (XY(2,0)-XY(1,0));
    sum    += (XY(2,1)-XY(1,1)) * (XY(2,1)-XY(1,1));
    len[1]  = sqrt(sum);
    // segment 3
    sum     = (XY(0,0)-XY(2,0)) * (XY(0,0)-XY(2,0));
    sum    += (XY(0,1)-XY(2,1)) * (XY(0,1)-XY(2,1));
    len[2]  = sqrt(sum);
    }

 } // end EdgeLengths



/**

Computes Integral N dA = 1/2 sum_1...n Wi Ji Ni. This formulation also
gives a meaningful area if the element boundaries are curved. For a
description of numerical integration  of isoparametric quadratic
triangular elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The area (m2) of the finite element.
*/
double64  IsoparametricLinearTriangle::Volume()
{
    double64  area;
    size_t  i;

    // numerical integration:
    // looping over the Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( area=0.0, i=0; i<gpe; i++ )
      {
         // getting interpolation function derivatives
         RS[0] = rr[i];
         RS[1] = ss[i];
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix
         if ( use2Dto3Djacobi ) area += Jacobi( RS ) * W[i];
         else {
              dNr( rr[i], ss[i], DNR );
              dNs( rr[i], ss[i], DNS );
              // getting global intpol. function derivative matrix and determinant of
              // byproduct Jacobian matrix
              Jacobian( DNR, DNS );
              area += JacobianInverse() * W[i];
           }
      }
  return area;
}

/**

Compute derivatives of interpolation functions at corresponding nodes with
respect to the global coordinate system. If a transformation from 2D local
to 3D global space is applied, the function Jacobi is called.

@section arguments Input Arguments

The element is used to obtain the global interpolation of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN3 The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricLinearTriangle::dN( DenseMatrix<DM_MIN>& DN3 )
  {
     DN3.Resize(dim,npe);
     double64 detJ(std::numeric_limits<double>::quiet_NaN());

     // Jacobian transformation to global coordinate system
     if ( use2Dto3Djacobi )
       {
          vector<double64> EFG(dim);
          double64 det_inverse;

          for ( size_t i=0; i<npe; i++ ) {
                RS[0] = NXY(i,0);
                RS[1] = NXY(i,1);
                // generate the 2 x 3 shape function derivative matrix
                detJ = Jacobi( RS, EFG, JMAT );
                det_inverse = 1. / ( detJ * detJ );
                DN3(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN3(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN3(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN3(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN3(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN3(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             }
           return;
        }

     // if the triangle 2D
     DenseMatrix<DM_MIN> TEMP(dim,1);
     for ( size_t i=0; i<npe; i++ )
       {
           dNr( NXY(i,0), NXY(i,1), DNR );
           dNs( NXY(i,0), NXY(i,1), DNS );
           Jacobian( DNR, DNS );
           JacobianInverse();

           TEMP(0,0)=DNR[i];TEMP(1,0)=DNS[i];

           JINV *= TEMP;

           DN3(0,i) = JINV(0,0);
           DN3(1,i) = JINV(1,0);

           JINV.Resize(dim,dim);
       }
  } // end dN






/// Projection function from rs->xy(z)
void
IsoparametricLinearTriangle::ParametricToPhysical(std::vector<double64> &rst,
                                                  std::vector<double64>& xyz)
{
    vector<double64> N(npe);
    Nrs(rst[0],rst[1], N );

    for(size_t i=0; i<dim; i++)
        xyz[i]=0.0;

    if(dim==2)
    {
        for(size_t i=0; i<npe; i++)
        {
            xyz[0]+=XY(i,0)*N[i];
            xyz[1]+=XY(i,1)*N[i];
        }
    }
    else if(dim==3)
    {
        for(size_t i=0; i<npe; i++)
        {
            xyz[0]+=XY(i,0)*N[i];
            xyz[1]+=XY(i,1)*N[i];
            xyz[2]+=XY(i,2)*N[i];
        }
    }
}

/** Projection function from xyz->rst
*/
void
IsoparametricLinearTriangle::PhysicalToParametric(
                                                std::vector<double64>& rSt,
                                                const std::vector<double64>& xyz
                                                 )
{

    if( dim == 2 )
    {

        // no sign is taken, thus method works with cw and ccw node numbering
        double64 ae2 = 1.0 /( ( XY(1,0)*XY(2,1) + XY(0,0)*XY(1,1) +
                                XY(0,1)*XY(2,0) - XY(2,1)*XY(0,0) -
                                XY(2,0)*XY(1,1) - XY(1,0)*XY(0,1) ) );

        // calculation the exponents of the interpolation functions
        //const double64 a0 = XY(1,0) * XY(2,1) - XY(2,0) * XY(1,1);
        const double64 a1 = XY(2,0) * XY(0,1) - XY(0,0) * XY(2,1);
        const double64 a2 = XY(0,0) * XY(1,1) - XY(1,0) * XY(0,1);

        //const double64 b0 = XY(1,1) - XY(2,1);
        const double64 b1 = XY(2,1) - XY(0,1);
        const double64 b2 = XY(0,1) - XY(1,1);

        //const double64 c0 = XY(2,0) - XY(1,0);
        const double64 c1 = XY(0,0) - XY(2,0);
        const double64 c2 = XY(1,0) - XY(0,0);

        //N[0] = ae2 * (a0 + b0 * xyz[0] + c0 * xyz[1]);
        //N[1] = ae2 * (a1 + b1 * xyz[0] + c1 * xyz[1]);
        //N[2] = ae2 * (a2 + b2 * xyz[0] + c2 * xyz[1]);

        rSt.resize(2U);
        rSt[0] = ae2 * (a1 + b1 * xyz[0] + c1 * xyz[1]);
        rSt[1] = ae2 * (a2 + b2 * xyz[0] + c2 * xyz[1]);


    }else{

        //Edges of triangle 123
        const double64  X12 = XY(1,0) - XY(0,0);
        const double64  X23 = XY(2,0) - XY(1,0);
        const double64  X31 = XY(0,0) - XY(2,0);
        const double64  Y12 = XY(1,1) - XY(0,1);
        const double64  Y23 = XY(2,1) - XY(1,1);
        const double64  Y31 = XY(0,1) - XY(2,1);
        const double64  Z12 = XY(1,2) - XY(0,2);
        const double64  Z23 = XY(2,2) - XY(1,2);
        const double64  Z31 = XY(0,2) - XY(2,2);

        //Vectors from point P to the points of triangle 123
        const double64  XP1 = XY(0,0) - xyz[0];
        //const double64  XP2 = XY(1,0) - xyz[0];
        const double64  XP3 = XY(2,0) - xyz[0];
        const double64  YP1 = XY(0,1) - xyz[1];
        //const double64  YP2 = XY(1,1) - xyz[1];
        const double64  YP3 = XY(2,1) - xyz[1];
        const double64  ZP1 = XY(0,2) - xyz[2];
        //const double64  ZP2 = XY(1,2) - xyz[2];
        const double64  ZP3 = XY(2,2) - xyz[2];

        // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
        const double64  XNRM = -Y12*Z31 + Z12*Y31;
        const double64  YNRM = -Z12*X31 + X12*Z31;
        const double64  ZNRM = -X12*Y31 + Y12*X31;

        // XN12,YN12,ZN12 is normal to the edge 12 in triangle plane
        const double64  XN12 = Y12*ZNRM - Z12*YNRM;
        const double64  YN12 = Z12*XNRM - X12*ZNRM;
        const double64  ZN12 = X12*YNRM - Y12*XNRM;

        // XN23,YN23,ZN23 is normal to the edge 23 in triangle plane
        //const double64  XN23 = Y23*ZNRM - Z23*YNRM;
        //const double64  YN23 = Z23*XNRM - X23*ZNRM;
        //const double64  ZN23 = X23*YNRM - Y23*XNRM;

        // XN31,YN31,ZN31 is normal to the edge 31 in triangle plane
        const double64  XN31 = Y31*ZNRM - Z31*YNRM;
        const double64  YN31 = Z31*XNRM - X31*ZNRM;
        const double64  ZN31 = X31*YNRM - Y31*XNRM;

        //Scalar product of the edge and correponding normal to the neighbor edge = 2* Square of triangle
        //const double64  A123 = XN23*X12 + YN23*Y12 + ZN23*Z12;
        const double64  A231 = XN31*X23 + YN31*Y23 + ZN31*Z23;
        const double64  A312 = XN12*X31 + YN12*Y31 + ZN12*Z31;

        //Scalar product of the vector from point P and correponding normal to the edge = 2* Square of triangle P(edge): (P23), (P31), (P12)
        //const double64  AP23  = XN23*XP2 + YN23*YP2 + ZN23*ZP2;
        const double64  AP31  = XN31*XP3 + YN31*YP3 + ZN31*ZP3;
        const double64  AP12  = XN12*XP1 + YN12*YP1 + ZN12*ZP1;

        // Finds the values of three basis functions  at point 'xyz' given
        // points of triangle in 3d, X(1-3),Y(1-3),Z(1-3).
        //N.resize(dim);
        //N[0] = AP23 / A123;
        //N[1] = AP31 / A231;
        //N[2] = AP12 / A312;

        rSt.resize(2U);
        rSt[0] = AP31 / A231;
        rSt[1] = AP12 / A312;

        /*
        // alternative approach
        //Edges of triangle 123
        const double64  X12 = XY(1,0) - XY(0,0);
        //const double64  X23 = XY(2,0) - XY(1,0);
        const double64  X31 = XY(0,0) - XY(2,0);
        const double64  Y12 = XY(1,1) - XY(0,1);
        //const double64  Y23 = XY(2,1) - XY(1,1);
        const double64  Y31 = XY(0,1) - XY(2,1);
        const double64  Z12 = XY(1,2) - XY(0,2);
        //const double64  Z23 = XY(2,2) - XY(1,2);
        const double64  Z31 = XY(0,2) - XY(2,2);

        //Vectors from point P to the points of triangle 123
        const double64  XP1 = XY(0,0) - xyz[0];
        const double64  XP2 = XY(1,0) - xyz[0];
        const double64  XP3 = XY(2,0) - xyz[0];
        const double64  YP1 = XY(0,1) - xyz[1];
        const double64  YP2 = XY(1,1) - xyz[1];
        const double64  YP3 = XY(2,1) - xyz[1];
        const double64  ZP1 = XY(0,2) - xyz[2];
        const double64  ZP2 = XY(1,2) - xyz[2];
        const double64  ZP3 = XY(2,2) - xyz[2];

        // 2 * ( Area of Triangle 123 )
        const double64  XN312 = -Y12*Z31 + Z12*Y31;
        const double64  YN312 = -Z12*X31 + X12*Z31;
        const double64  ZN312 = -X12*Y31 + Y12*X31;
        const double64  A123  = sqrt(XN312*XN312 + YN312*YN312 + ZN312*ZN312);

        // 2 * ( Area of Triangle P31 )
        const double64  XN3P1 = -YP1*ZP3 + ZP1*YP3;
        const double64  YN3P1 = -ZP1*XP3 + XP1*ZP3;
        const double64  ZN3P1 = -XP1*YP3 + YP1*XP3;
        const double64  AP31  = sqrt(XN3P1*XN3P1 + YN3P1*YN3P1 + ZN3P1*ZN3P1);

        // 2 * ( Area of Triangle P12 )
        const double64  XN1P2 = -YP1*ZP2 + ZP1*YP2;
        const double64  YN1P2 = -ZP1*XP2 + XP1*ZP2;
        const double64  ZN1P2 = -XP1*YP2 + YP1*XP2;
        const double64  AP12  = sqrt(XN1P2*XN1P2 + YN1P2*YN1P2 + ZN1P2*ZN1P2);

        // 2 * ( Area of Triangle P23 )
        //const double64  XN2P3 = -YP2*ZP3 + ZP2*YP3;
        //const double64  YN2P3 = -ZP2*XP3 + XP2*ZP3;
        //const double64  ZN2P3 = -XP2*YP3 + YP2*XP3;
        //const double64  AP23  = sqrt(XN2P3*XN2P3 + YN2P3*YN2P3 + ZN2P3*ZN2P3);

        // Finds the values of three basis functions  at point 'xyz' given
        //N.resize(dim);
        //N[0] = AP23 / A123;
        //N[1] = AP31 / A123;
        //N[2] = AP12 / A123;

        rSt.resize(2U);
        rSt[0] = AP31 / A123;
        rSt[1] = AP12 / A123;
        */
    }
}



/**

Computes derivatives of shape functions at a point XY given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.
*/
double64
IsoparametricLinearTriangle::dN( DenseMatrix<DM_MIN>& DN2,
                                 const vector<double64>& xyz  )
{
  vector<double64> rst(dim);
  double64 detJ(std::numeric_limits<double>::quiet_NaN());

  PhysicalToParametric(rst, xyz);

    // here find derivatives
  dNr( rst[0], rst[1], DNR );
  dNs( rst[0], rst[1], DNS );

  if ( use2Dto3Djacobi ){
     vector<double64> EFG(dim);
     double64 det_inverse;

       for ( size_t i=0; i<npe; i++ ) {
              detJ = Jacobi( rst, EFG, JMAT );
              det_inverse = 1.0 / ( detJ * detJ );
              DN2(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
              DN2(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
              DN2(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
              DN2(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
              DN2(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
              DN2(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
      }
     return detJ;
  }

  Jacobian( DNR, DNS );
  detJ = JacobianInverse();

  DN2.Resize(dim,dim);
  DN2  = JINV;

  DenseMatrix<DM_MIN>DN(dim,npe);
  for(int i=0;i<npe;i++)
  {
      DN(0,i)=DNR[i];
      DN(1,i)=DNS[i];
  }

  DN2*=DN;

  /////////////////////////////// Debug printout ///////////////////////////////////////////////
  //cout<<" IsoparametricLinearTriangle3D::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
  //rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
  //cout<<" IsoparametricLinearTriangle3D::dN  Matrix DN2: "<<endl;
  //DN2.Out();
  /////////////////////////////// Debug printout ///////////////////////////////////////////////
  return detJ;
}

/**

When a surface element is located in a 3D space, the Jacobian matrix reduces
to a 2x3 matrix. This function allows to calculate the jacobi (essentially
the determinant of the 2x3 Jacobian matrix), which can be used to extrapolate
the values of the interpolation functions and interpolation function derivatives
that were calculated on a 2D reference element in local coordinates to the
surface element that is located in the 3D space in global coordinates

@section arguments Input Arguments

A reference to the parent element and a vector containing the local coordinates
at which the Jacobian shall be evaluated are the input function arguments. The
method returns the jacobi. In the function argument, it also returns a vector
containing the L2 norms (squared) of the 2x3 Jacobian matrix (returned in the
next argument), and the interpolation function derivatives with respect to
the local coordinates r and s.

A second overloaded function exists that returns only the jacobi and uses the
reference to the parent element and the vector containing the local coordinates
as input arguments.
*/
double64
IsoparametricLinearTriangle::Jacobi(    const vector<double64>& rs,
                                        vector<double64>& EFG,
                                        DenseMatrix<DM_MIN>& J )
 {
    // 1. Compute 3x2 Jacobian Matrix
    dNr( rs[0], rs[1], DNR );
    dNs( rs[0], rs[1], DNS );

    // compute 2x3 Jacobian matrix (expects 6x3 x-y-z coordinate matrix)
    J.Resize(2,dim);
    J.Zero();
    for ( size_t i=0; i<XY.Cols(); i++ )
      for ( size_t j=0; j<XY.Rows(); j++ ) {
           J(0,i) += DNR[j] * XY(j,i);
           J(1,i) += DNS[j] * XY(j,i);
        }
    // 2. Compute Jacobi J' := "determinant" of the 3x2 Jacobian
    // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)

    // compute E, F, and g
    EFG.resize(3);
    fill( EFG.begin(), EFG.end(), 0.0 );

    for ( size_t i=0; i<dim; i++ ) {
        EFG[0] += J(0,i) * J(0,i);
        EFG[1] += J(0,i) * J(1,i);
        EFG[2] += J(1,i) * J(1,i);
      }

   // compute J'
   return sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );

 }  // end Jacobi


double64 IsoparametricLinearTriangle::Jacobi( const vector<double64>& rs )
{
    // 1. Compute 2x3 Jacobian Matrix
    dNr( rs[0], rs[1], DNR );
    dNs( rs[0], rs[1], DNS );

    // compute 2x3 Jacobian matrix
    JMAT.Resize(2,dim);
    JMAT.Zero();

    for ( size_t i=0; i<XY.Cols(); i++ )
      for ( size_t j=0; j<XY.Rows(); j++ ) {
           JMAT(0,i) += DNR[j] * XY(j,i);
           JMAT(1,i) += DNS[j] * XY(j,i);
        }
    // 2. Compute Jacobi J' := "determinant" of the 4x2 Jacobian
    // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)

    // compute E, F, and g
    double64 E(0.0), F(0.0), G(0.0);

    for ( size_t i=0; i<dim; i++ ) {
         E += JMAT(0,i) * JMAT(0,i);
         F += JMAT(0,i) * JMAT(1,i);
         G += JMAT(1,i) * JMAT(1,i);
      }

   // compute J'
   return  sqrt( E * G - F * F );

}  // end Jacobi


/**

Computes the values of the interpolation functions at the global coordinate 'xy'.
Since a Jacobian matrix can only be obtained for points which lie inside
the local triangle, such a transformation cannot be performed here.
Because of this reason global area-coordinate functions are employed
to obtain the interpolation function values at the global point. These
functions are borrowed from the NaturalQuadratic triangle.

Note, that the use of global interpolation functions implies that the method
only gives highly accurate results if the element boundaries are
not curved.

@section arguments Input Arguments

The parent Element, the vector which will hold the interpolation function
values, and the coordinate vector 'xy' in the global coordinate system.
Note that if the global coordinates given by 'xy' do not lie within the
global extent of the quadratic triangular element, the interpolation function values
will no longer add to one and the interpolation will then be erroneous.

@param N The second method argument, the vector FN, will hold the interpolation function
values as computed at the point 'xy'.

@section implementation Implementation

The inverse of the Jacobian matrix is found for the global point 'xy'.
Then 'xy' transposed is pre-multiplied with the Jacobian to find the
local coordinate pair 'rs' corresponding to 'xy'. Using these local
coordinates the interpolation function values are found.

@section application Application

To interpolate a property value withing the quadratic triangular
element.

 */
void
IsoparametricLinearTriangle::N( vector<double64>& N, const vector<double64>& xyz )
{
   vector<double64> rs(parametricDimehsions);

   PhysicalToParametric(rs, xyz);

   Nrs(rs[0],rs[1], N);

} // end


/**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the gauss
point i (i=0...i=gauss points-1). If a transformation from 2D local
to 3D global space is applied, the function Jacobi is called.

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
IsoparametricLinearTriangle::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B,
                                                    size_t gauss_point )
 {
    double64 det;

    assert( gauss_point < gpe );

    if ( use2Dto3Djacobi ){
         vector<double64> EFG(dim);
         RS[0] = rr[gauss_point];
         RS[1] = ss[gauss_point];

         det = Jacobi( RS, EFG, JMAT );
         double64 det_inverse = 1.0 / ( det * det );

         B.Resize(dim,npe);

         for ( size_t i=0; i<npe; i++ )
           {
             B(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
             B(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             B(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
             B(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             B(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
             B(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
           }
      return det;
    }

    // compute local test-function derivative matrix at gauss point
    // get local interpolation function derivatives at Gauss point
    dNr( rr[gauss_point], ss[gauss_point], DNR );
    dNs( rr[gauss_point], ss[gauss_point], DNS );

    // compute Jacobian matrix, its determinant and inverse
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2];

    B = JINV * B;

  return det;
 }

/**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J at the desired node
point.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the
point i (i=0...i=gauss points-1). If a transformation from 2D local
to 3D global space is applied, the function Jacobi is called.

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
IsoparametricLinearTriangle::dN_AtNode(
                                        DenseMatrix<DM_MIN>& B,
                                        size_t nd
                                      )
{
  double64  det;

  assert( nd < npe );
    if ( use2Dto3Djacobi ){
     vector<double64> EFG(dim);
     RS[0] = NXY(nd,0);
     RS[1] = NXY(nd,1);

       det = Jacobi( RS, EFG, JMAT );
       double64 det_inverse = 1.0 / ( det * det );

       B.Resize(dim,npe);

        for ( size_t i=0; i<npe; i++ ){
          B(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
          B(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
          B(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
          B(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
          B(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
          B(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
        }

        return det;
   } //if


    // normal 2D case
    // --------------------------------------
    dNr( NXY(nd,0), NXY(nd,1), DNR );
    dNs( NXY(nd,0), NXY(nd,1), DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2];

    B = JINV * B;

  return det;
}


double64 IsoparametricLinearTriangle::dN_AtBarycenter(
                                                        DenseMatrix<
                                                        DM_MIN>& B
                                                      )
{
  double64   det;
  const double64 OneThird=1./3.0;

  if ( use2Dto3Djacobi ){
   vector<double64> EFG(dim);
   RS[0] = OneThird;
   RS[1] = OneThird;

   det = Jacobi( RS, EFG, JMAT );
   double64 det_inverse = 1.0 / ( det * det );

   B.Resize(dim,npe);

     for ( size_t i=0; i<npe; i++ ){
       B(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
       B(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
       B(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
       B(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
       B(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
       B(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
     }
   return det;
  }

    // 2D case
    // ----------------------------------------------------
    dNr( OneThird, OneThird, DNR );
    dNs( OneThird, OneThird, DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2];

    B = JINV * B;

  return det;
}


/**

Method returns a vector with a size of 3, containing the consecutively
ordered local node numbers of the side of the element which lies at the
indicated model boundary. If the element is not on the model boundary an
error is reported and the vector is initialized to unspecified.

@section arguments Input Arguments

The parent Element, the target boundary of the the Model of elements,
and a Meschach vector which will hold the the local indices of the identified
nodes (0...5).

@param fnids The resulting local node id's are returned into the third method argument.

@section implementation Implementation

While only the element knows which nodes are located at the mode boundary,
the FiniteElement knows in which order these appear.

@section application Application

To assign Neumann boundary conditions with a PDE operator for surface
integrals.
*/
void
IsoparametricLinearTriangle::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                         vector<size_t>& fnids )
{
   fnids.resize(3);

  // in 3D, all nodes may be at model boundary
  if ( use2Dto3Djacobi && bnodes.size() != 2 )
     throw csmp::Exception( CSMP_FATAL_ERROR, "IsoparametricLinearTriangle::ConsecutiveNodesAtBoundary",
                            "Cannot resolve node sequence for element boundary",
                            "Probably because element lies at two boundaries simultaneously" );

  if ( bnodes.size() > 2 )
     throw csmp::Exception( CSMP_ERROR, "IsoparametricLinearTriangle::ConsecutiveNodesAtBoundary",
                              "Cannot resolve node sequence for element boundary",
                              "Probably because element lies at two boundaries simultaneously" );
/*
     // the midside nodes determine the boundaries so they are sought after
     if ( bnodes[0] == 3 || bnodes[1] == 3 || bnodes[2] == 3 ) {
           fnids[0] = 0, fnids[1] = 3, fnids[2] = 1;
           return;
       }
     if ( bnodes[0] == 4 || bnodes[1] == 4 || bnodes[2] == 4 ) {
           fnids[0] = 1, fnids[1] = 4, fnids[2] = 2;
           return;
       }
     if ( bnodes[0] == 5 || bnodes[1] == 5 || bnodes[2] == 5 ) {
           fnids[0] = 2, fnids[1] = 5, fnids[2] = 0;
           return;
       }

*/
} // end ConsecutiveNodesAtBoundary



void
IsoparametricLinearTriangle::OutputNodeDataToVTK( const char* file_name,
                                                  const char* var_name,
                                                  DenseMatrix<DM_MIN>& DATA ) const
 {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     sprintf( elmt, "%lu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cout <<"\nIsoparametricLinearTriangle3D::OutputNodeDataToVTK ";
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
     ofs <<"POINTS " << npe <<" float"<< endl;
     for ( size_t i=0; i<npe; i++ ){
       if(dim==2){
        for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
           ofs<<0.0;
           ofs << endl;
        }
        else
        if(dim==3){
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
        }
     }
     ofs << endl;

    // 3. writing CELLS (cell-size and member nodes (point)) -corect
    // Cells
    // -----------------------------------------------------
    ofs <<"CELLS "<< 1 <<" "<< 4 << endl;
    ofs << 3 <<" 0 1 2" << endl;
    ofs << endl;

    // 4. writing CELL_TYPES - for the
    // ---------------------
    ofs <<"CELL_TYPES "<< 1 << endl;
    ofs << 5 << endl; // VTK_TRI
    ofs << endl;

    // 5. writing POINT_DATA point-type data values
    // --------------------------------------------
    // Unfortunately the data can only be output as nodal variables
    ofs <<"POINT_DATA "<< npe << endl;
    ofs.setf( ios::scientific );

    if ( DATA.Rows() == 1 ){
      ofs <<"SCALARS "<< var_name <<" float"<< endl;
      ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
      // matrix DATA is 1x9
         for ( size_t i=0; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
      ofs << endl;
    }
     else {
       ofs <<"VECTORS "<< var_name <<" float"<< endl;
       // variables have always 3 components since view screen is 3D
       // matrix DATA is vec-dim x 10
       if(dim==2) {
         for ( size_t i=0; i<DATA.Cols(); i++ ) {
           for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
              ofs<<0.0<<" ";
              ofs << endl;
         }
       }
       else {
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
            for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
            ofs << endl;
          }
       }
     }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricLinearTriangle3D::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;
 } // end OutputNodeDataToVTK


/**

Because, in the local coordinate system, the quadrilateral is always straight
sided, and the derivatives of the interpolation functions are bilinear,
variable values computed at integration points can simply be extrapolated
to the nodes.  This involves the steps:

1. Compute the interpolation function coefficients for the quad which
   is defined by the four integration points.

2. For each node point compute the values of the interpolation functions
   and use these to extrapolate the values of the variables at the nodes.
*/
void
IsoparametricLinearTriangle::ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                        const vector<double64>& IVAR,
                                                              vector<double64>& NVAR )
const
{
   static double64 a[3], b[3], c[3], intpol[3], ae2;
   static bool      first_call(true);

   vector<double64>  sum(nvars);

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   //assert( gpe == 3 );
   assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

   // 1. Compute the interpolation function coefficients for the triangle which
   //    is defined by the three integration points.
   if ( first_call )
   {
        a[0] = rr[1] * ss[2] - rr[2] * ss[1];
        a[1] = rr[2] * ss[0] - rr[0] * ss[2];
        a[2] = rr[0] * ss[1] - rr[1] * ss[0];

        b[0] = ss[1] - ss[2];
        b[1] = ss[2] - ss[0];
        b[2] = ss[0] - ss[1];

        c[0] = rr[2] - rr[1];
        c[1] = rr[0] - rr[2];
        c[2] = rr[1] - rr[0];

        // compute the area*2 of the triangle spanned by the integration points
        ae2 = 1.0 / ( ( rr[1]*ss[2] + rr[0]*ss[1] +  ss[0]*rr[2] - ss[2]*rr[0] -
                                                     rr[2]*ss[1] - rr[1]*ss[0] ) );
        first_call = false;
   }

   // 2. For each node point compute the values of the interpolation functions
   //    and use these to extrapolate the values of the variables at the nodes.
   for ( size_t i=0; i<npe; i++ )
     {
        // compute interpolation function values at node i
        for ( size_t j=0; j<gpe; j++ )
          intpol[j] = ae2 * (a[j] + b[j] * NXY(i,0) + c[j] * NXY(i,1));

        // carry out extrapolation
        for ( size_t k=0; k<nvars; k++ ) sum[k] = 0.0;
        for ( size_t j=0; j<gpe; j++ )
          for ( size_t k=0; k<nvars; k++ ) sum[k] += intpol[j] * IVAR[j*nvars + k];

        // store result in output vector
        for ( size_t k=0; k<nvars; k++ ) NVAR[i*nvars + k] = sum[k];
     }

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)



/// overwrites the base class method in order to get 2D to 3D mapping functionality
double64  IsoparametricLinearTriangle::JacobianInverse()
{
    if ( use2Dto3Djacobi )
      {
        // Compute Jacobi J' := "determinant" of the 3x2 Jacobian
        // ------------------------------------------------------
        // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)
          vector<double64> EFG(dim);
          fill( EFG.begin(), EFG.end(), 0.0 );

             for ( size_t i=0; i<dim; i++ ) {
                 EFG[0] += JAC(0,i) * JAC(0,i);
                 EFG[1] += JAC(0,i) * JAC(1,i);
                 EFG[2] += JAC(1,i) * JAC(1,i);
               }

         double64 detJ=sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );
         double64 det1 = 1.0 / detJ;


         //Standard inverse Jacobian definition with 3rd stroka all 1.0
         //JAC(2,0)=JAC(2,1)=JAC(2,2)=1.0;

         JINV(0,0) = ( JAC(1,1) * 1.0 - JAC(1,2) * 1.0 ) *  det1;
         JINV(1,0) = ( JAC(1,0) * 1.0 - JAC(1,2) * 1.0 ) * -det1;
         JINV(2,0) = ( JAC(1,0) * 1.0 - JAC(1,1) * 1.0 ) *  det1;
         JINV(0,1) = ( JAC(0,1) * 1.0 - JAC(0,2) * 1.0 ) * -det1;
         JINV(1,1) = ( JAC(0,0) * 1.0 - JAC(0,2) * 1.0 ) *  det1;
         JINV(2,1) = ( JAC(0,0) * 1.0 - JAC(0,1) * 1.0 ) * -det1;
         JINV(0,2) = ( JAC(0,1) * JAC(1,2) - JAC(0,2) * JAC(1,1) ) *  det1;
         JINV(1,2) = ( JAC(0,0) * JAC(1,2) - JAC(0,2) * JAC(1,0) ) * -det1;
         JINV(2,2) = ( JAC(0,0) * JAC(1,1) - JAC(0,1) * JAC(1,0) ) *  det1;


          // compute J'
          return detJ;
       }

    // the 2D case
    // compute determinant
    double64 detJ = JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);

    // inversion of J
    double64 dum = JAC(0,0) / detJ;
    JINV(0,0)  =  JAC(1,1) / detJ;
    JINV(0,1)  = -JAC(0,1) / detJ;
    JINV(1,0)  = -JAC(1,0) / detJ;
    JINV(1,1)  =  dum;

    if ( detJ <= 0 ) {
         std::cout <<"\nIsoparametricLinearTriangle::JacobianInverse: Erroneous determinant of Jacobian matrix: ";
         std::cout << detJ << std::endl;
         throw std::range_error("IsoparametricLinearTriangle::JacobianInverse");
    }

 return detJ;
}


/**
     normal points to front when looking at element with nodes 
     numbered counter clockwise.
     
     @attention in 2D, a 3D unit normal is returned, with the z-component
     equal to 1.
     
     @test OK - SKM 3/3/2016
*/
void IsoparametricLinearTriangle::UnitNormal( std::vector<double64>& vc ) const
 {
    vc.resize(3);

    // normal points into non-existant third dimension
    if ( dim == 2U ) {
         vc[0] = static_cast<double64>(0.);
         vc[1] = static_cast<double64>(0.);
         vc[2] = static_cast<double64>(1.);
         return;
      }

    double64 X12 = XY(1,0) - XY(0,0), // X
             X31 = XY(0,0) - XY(2,0),
             Y12 = XY(1,1) - XY(0,1), // Y
             Y31 = XY(0,1) - XY(2,1),
             Z12 = XY(1,2) - XY(0,2), // Z
             Z31 = XY(0,2) - XY(2,2);

    // normal to triangle (*- to flip to outside)
    vc[0]  = -Y12*Z31 + Z12*Y31;
    vc[1]  = -Z12*X31 + X12*Z31;
    vc[2]  = -X12*Y31 + Y12*X31;
   
    // normalization to unit length
    double64 length = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
    vc[0] /= length;
    vc[1] /= length;
    vc[2] /= length;
 }




/**
    outward pointing normals to the faces that lie on the opposite nodes
    
    @author SKM 18/2/2016
    
    @test OK - for 3D version
*/
void  IsoparametricLinearTriangle::UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const
 {
     assert( face < Faces() );
     // if this is a planar element in a 2D model
     if ( Dim() == 2 ) {
         unrml.resize(2);
         // nodes 1 and 2
         if ( face == 0 ) {
              mjl::Edge  normal( mjl::Point(XY(1,0),XY(1,1)), mjl::Point(XY(2,0),XY(2,1)) );
              // rotating edge clockwise to find outward pointing normal to face
              normal.Rot();
              normal.NormalizeTo( 1. );
              unrml[0] = normal.Destination()[0];
              unrml[1] = normal.Destination()[1];
              return;
           }
         // nodes 2 and 0
         if ( face == 1 ) {
              mjl::Edge  normal( mjl::Point(XY(2,0),XY(2,1)), mjl::Point(XY(0,0),XY(0,1)) );
              normal.Rot();
              normal.NormalizeTo( 1. );
              unrml[0] = normal.Destination()[0];
              unrml[1] = normal.Destination()[1];
              return;
           }
         // nodes 0 and 1
         if ( face == 2 ) {
              mjl::Edge  normal( mjl::Point(XY(0,0),XY(0,1)), mjl::Point(XY(1,0),XY(1,1)) );
              normal.Rot();
              normal.NormalizeTo( 1. );
              unrml[0] = normal.Destination()[0];
              unrml[1] = normal.Destination()[1];
              return;
           }
          return;
       }
   
     // if the triangle is suspended into 3D space
     unrml.resize(3);
     Point<3> face1(XY(2,0)-XY(1,0),XY(2,1)-XY(1,1),XY(2,2)-XY(1,2));
     Point<3> face2(XY(0,0)-XY(2,0),XY(0,1)-XY(2,1),XY(0,2)-XY(2,2));
     Point<3> enrml( crossProduct( face1, face2 ) );

     // the face normals are found as cross-products between element normal and edges
     // nodes 1 and 2
     if ( face == 0 ) {
          Point<3> nrml( crossProduct( enrml, face1 ) );
          nrml.NormalizeLengthTo(/* 1 */);
          // (-) to flip the normal to the outside
          unrml[0] = -nrml[0];
          unrml[1] = -nrml[1];
          unrml[2] = -nrml[2];
          return;
       }
     // nodes 2 and 0
     if ( face == 1 ) {
          Point<3> nrml( crossProduct( enrml, face2 ) );
          nrml.NormalizeLengthTo(/* 1 */);
          unrml[0] = -nrml[0];
          unrml[1] = -nrml[1];
          unrml[2] = -nrml[2];
          return;
       }
     // nodes 0 and 1
     if ( face == 2 ) {
          Point<3> face3(XY(1,0)-XY(0,0),XY(1,1)-XY(0,1),XY(1,2)-XY(0,2));
          Point<3> nrml( crossProduct( enrml, face3 ) );
          nrml.NormalizeLengthTo(/* 1 */);
          unrml[0] = -nrml[0];
          unrml[1] = -nrml[1];
          unrml[2] = -nrml[2];
       }
   
 } // end UnitNormalToFace




/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be updated
void  IsoparametricLinearTriangle::IntegrationPoint( size_t ip, vector<double64>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(NXY.Cols());
    xyz[0]=xyz[1]=0.;

    // local interpolation function values
    Nrs( rr[ip], ss[ip], NRST );

    // 2D
    if ( NXY.Cols() == 2U ) {
          for( size_t i=0U; i<npe; i++ ) {
               xyz[0] += XY(i,0) * NRST[i];
               xyz[1] += XY(i,1) * NRST[i];
            }
         return;
      }

    // 3D case
    xyz[2]=0.;
    for( size_t i=0U; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint




void
IsoparametricLinearTriangle::JacobianAtIntegrationPoint( size_t gauss_point )
 {
    if ( !use2Dto3Djacobi ) {
         dNr( rr[gauss_point], ss[gauss_point], DNR );
         dNs( rr[gauss_point], ss[gauss_point], DNS );
         Jacobian( DNR, DNS );
         return;
    }

   vector<double64> EFG(3);
   RS[0] = rr[gauss_point];
   RS[1] = ss[gauss_point];

   Jacobi( RS, EFG, JAC );

 }



void
IsoparametricLinearTriangle::IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS)
 {
    IPPHYS.Resize(gpe,dim);
    vector<double64> outxyz(dim);
    const size_t numberOfParametricDims=2;
    vector<double64> rst(numberOfParametricDims);

    for(size_t i=0;i<gpe;i++){
        rst[0]=rr[i];rst[1]=ss[i];
        ParametricToPhysical( rst, outxyz);
        for(size_t j=0;j<dim;j++) IPPHYS(i,j)=outxyz[j];
    }
}

void
IsoparametricLinearTriangle::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords(0,0) = NXY(0,0);
    matCoords(0,1) = NXY(0,1);
    matCoords(1,0) = NXY(1,0);
    matCoords(1,1) = NXY(1,1);
    matCoords(2,0) = NXY(2,0);
    matCoords(2,1) = NXY(2,1);
}


double64  IsoparametricLinearTriangle::JacobianDeterminant()
{
  if (dim == 3)
    {
      // compute E, F, and g
      std::vector<double64> EFG(3);
      fill( EFG.begin(), EFG.end(), 0.0 );

      for ( size_t i=0; i<dim; i++ )
       {
         EFG[0] += JAC(0,i) * JAC(0,i);
         EFG[1] += JAC(0,i) * JAC(1,i);
         EFG[2] += JAC(1,i) * JAC(1,i);
       }

      // compute J'
      return sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );
     }
   //dim == 2
   return JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);
}


void IsoparametricLinearTriangle::JacobianAt( const std::vector<double64>& rst )
{
    if ( !use2Dto3Djacobi )
    {
         dNr( rst[0], rst[1], DNR );
         dNs( rst[0], rst[1], DNS );
         Jacobian( DNR, DNS );
         return;
    }

  vector<double64> EFG;
  Jacobi( rst, EFG, JAC );
}

} // end namespace csmp
