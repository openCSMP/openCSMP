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
IsoparametricLinearTriangle::IsoparametricLinearTriangle( uint32_t dimensions,
                                                          uint32_t ipoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_TRIANGLE, true, true, 1U ),
    use2Dto3Djacobi_( !(dimensions==2) ),
    NXY_(3,dimensions),
    JMAT_(dimensions,dimensions),
    RS_(2)
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

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    W_.resize( gpe );
    rr_.resize( gpe );
    ss_.resize( gpe );

    // full 4-point integration basis
    if( ipoints == 1 )
    {
        rr_[0] = 1./3.;
        ss_[0] = 1./3.;
        W_[0]  = .5;
    }
    else if( ipoints == 3 )
    {
        rr_[0] = 0.5;    rr_[1] = 0.5;   rr_[2] = 0.0;
        ss_[0] = 0.0;    ss_[1] = 0.5;   ss_[2] = 0.5;
        W_[0]  = 1./6.;  W_[1]  = 1./6.; W_[2]  = 1./6.;
    }
    else if( ipoints == 4 )
    {
        rr_[0] = 1./3.;    rr_[1] = 0.2;    rr_[2] = 0.6;    rr_[3] = 0.2;
        ss_[0] = 1./3.;    ss_[1] = 0.2;    ss_[2] = 0.2;    ss_[3] = 0.6;
        W_[0]  = -27./96.; W_[1] = 25./96.; W_[2] = 25./96.; W_[3]  = 25./96.;
    }
    else
    {
        cerr<<" IsoparametricLinearTriangle::IsoparametricLinearTriangle Number of integration points not supported: "<<ipoints<<endl;
        throw std::range_error
        ("***ERROR: IsoparametricLinearTriangle::IsoparametricLinearTriangle N of integration points should be 1,3 or 4");
    }

    // local node coordinates are defined as r==ksi, s==nu, t==mu.
    NXY_(0,0) = 0.0; NXY_(0,1) = 0.0;
    NXY_(1,0) = 1.0; NXY_(1,1) = 0.0;
    NXY_(2,0) = 0.0; NXY_(2,1) = 1.0;

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
                    double r,
                    double s,
                    std::vector<double>& N ) const
{
   N.resize(npe);
   N[0] = 1-r-s;
   N[1] = r;
   N[2] = s;
}

void
IsoparametricLinearTriangle::Nrs(
                    double r,
                    double s,
                    double* N ) const
{
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
                double,
                double,
                std::vector<double>& DNR ) const
{
   DNR.resize(npe);
   DNR[0] = -1.;
   DNR[1] =  1.;
   DNR[2] =  0.;
}


void IsoparametricLinearTriangle::dNs(
                double ,
                double ,
                std::vector<double>& DNS ) const
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
IsoparametricLinearTriangle::N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N )
 {
    assert( ip < gpe );

    // local interpolation function values
    Nrs( rr_[ip], ss_[ip], N );

 } // end N_AtIntegrationPoint




/// segments are numbered like faces
void
IsoparametricLinearTriangle::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
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
    std::cerr <<"IsoparametricLinearTriangle::NodesOfSegment: Segment ID out of range: "<< segm_id << std::endl;
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
/*
void IsoparametricLinearTriangle::NodesOfFace( uint32_t face_id, vector<uint32_t>& fnids ) const
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
*/



vector<uint32_t>  IsoparametricLinearTriangle::NodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nIsoparametricLinearTriangle::NodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }


vector<uint32_t>  IsoparametricLinearTriangle::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nIsoparametricLinearTriangle::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }



/// returns the local  numbers of the nodes at the other end of the sgment that the argument node is on
std::vector<uint32_t>  IsoparametricLinearTriangle::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nIsoparametricLinearTriangle::NodesConnectedTo: node "<< node_id <<" does not exist.";
    return vector<uint32_t>{};
  }



/** Returns local node ids of the nodes which sit on the corner of
the quadratic triagular element.

@param ids Returns a Meschpp integer vector with the 3 local corner node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).
*/
void
IsoparametricLinearTriangle::CornerNodes( std::vector<uint32_t>& ids ) const
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
void  IsoparametricLinearTriangle::MidSideNodes( std::vector<uint32_t>& ) const
 {
     std::cerr <<"\nIsoparametricLinearTriangle::MidSideNodes Not present..."<<std::endl;
 }




void  IsoparametricLinearTriangle::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }


CSMP_FEM_TYPE IsoparametricLinearTriangle::ElementTypeOfFace( uint32_t )  const
 {
    return ISOPARAMETRIC_LINEAR_BAR;
 }


double IsoparametricLinearTriangle::WeightAtIntegrationPoint( uint32_t i )
const { return W_[i]; }


void  IsoparametricLinearTriangle::N_AtBaryCenter( std::vector<double>& N )
 {
    N.resize(npe);
    Nrs(1./3.,1./3.,N);
 }



void IsoparametricLinearTriangle::Dimensions( uint32_t dimensions )
 {
   dim = dimensions;
   if ( dim == 3 ) use2Dto3Djacobi_ = true;
   else            use2Dto3Djacobi_ = false;
   // base class matrices
   XY.Resize(npe,dim);
   JAC.Resize(dim,dim);
   JINV.Resize(dim,dim);
   JMAT_.Resize(dim,dim);
 }



/** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double
IsoparametricLinearTriangle::AspectRatio()
{
   vector<double> vec(spe);

   EdgeLengths( vec );

   // order segment
   set<double> segms;

   for ( uint32_t i{0U}; i<spe; i++ ) segms.insert( vec[i] );

   double segm1 = (*segms.begin()),
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
double  IsoparametricLinearTriangle::InnerRadius()
{
   // Alternative radius
/*
   vector<double> rsCenter(parametricDimensions_), rsBaseR(parametricDimensions_);
   rsCenter[0]=1-0.5*sqrt(2.0); rsCenter[1]=1-0.5*sqrt(2.0);
   rsBaseR [0]=1-0.5*sqrt(2.0); rsBaseR[1]=0.0;

   vector<double> xyzCenter(dim),xyzBaseR(dim);

   ParametricToPhysical(rsCenter,xyzCenter);
   ParametricToPhysical(rsBaseR,xyzBaseR);

   double radius1=0.;

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
     throw csmp::Exception( WARNING, "IsoparametricLinearTriangle::InnerRadius",
                    "Approximate radius - not valid for given large aspect-ratio element");

   return radius1;
  */

    NRST.resize(spe);
    double  sum(0.);

    EdgeLengths( NRST );
    for ( uint32_t i{0U}; i<spe; i++ ) sum += NRST[i];

    if(AspectRatio()>4.)
     cerr<<"\nIsoparametricLinearTriangle::InnerRadius: WARNING: function not applicable for this high element aspect ratio.\n"<<endl;

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
IsoparametricLinearTriangle::EdgeLengths( std::vector<double>& len )
{
    double sum=0.0;
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
double  IsoparametricLinearTriangle::Volume()
{
    double  area{0.};

    // numerical integration:
    // looping over the Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( uint32_t i{0U}; i<gpe; i++ )
      {
         // getting interpolation function derivatives
         RS_[0] = rr_[i];
         RS_[1] = ss_[i];
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix
         if ( use2Dto3Djacobi_ ) area += Jacobi( RS_ ) * W_[i];
         else {
              dNr( rr_[i], ss_[i], DNR );
              dNs( rr_[i], ss_[i], DNS );
              // getting global intpol. function derivative matrix and determinant of
              // byproduct Jacobian matrix
              Jacobian( DNR, DNS );
              area += JacobianInverse() * W_[i];
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

@param DN  The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricLinearTriangle::dN( DenseMatrix<DM_MIN>& DN )
  {
     DN.Resize(dim,npe);
     double detJ(std::numeric_limits<double>::quiet_NaN());

     // Jacobian transformation to global coordinate system
     if ( use2Dto3Djacobi_ )
       {
          vector<double> EFG(dim);
          for ( uint32_t i{0U}; i<npe; i++ ) {
                RS_[0] = NXY_(i,0);
                RS_[1] = NXY_(i,1);
                // generate the 2 x 3 shape function derivative matrix
                detJ = Jacobi( RS_, EFG, JMAT_ );
                double det_inverse = 1. / ( detJ * detJ );
                DN(0,i)  = det_inverse * JMAT_(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN(0,i) += det_inverse * JMAT_(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN(1,i)  = det_inverse * JMAT_(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN(1,i) += det_inverse * JMAT_(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN(2,i)  = det_inverse * JMAT_(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN(2,i) += det_inverse * JMAT_(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             }
           return;
        }

     // if the triangle 2D
     DenseMatrix<DM_MIN> TEMP(dim,1);
     for ( uint32_t i{0U}; i<npe; i++ )
       {
           dNr( NXY_(i,0), NXY_(i,1), DNR );
           dNs( NXY_(i,0), NXY_(i,1), DNS );
           Jacobian( DNR, DNS );
           JacobianInverse();

           TEMP(0,0)=DNR[i];
           TEMP(1,0)=DNS[i];

           JINV *= TEMP;

           DN(0,i) = JINV(0,0);
           DN(1,i) = JINV(1,0);

           JINV.Resize(dim,dim);
       }
  } // end dN






/// Projection function from rs->xy(z)
void
IsoparametricLinearTriangle::ParametricToPhysical(std::vector<double> &rst,
                                                  std::vector<double>& xyz)
{
    vector<double> N(npe);
    Nrs(rst[0],rst[1], N );

    for( uint32_t i{0U}; i<dim; i++)
        xyz[i]=0.0;

    if(dim==2)
    {
        for( uint32_t i{0U}; i<npe; i++)
        {
            xyz[0]+=XY(i,0)*N[i];
            xyz[1]+=XY(i,1)*N[i];
        }
    }
    else if(dim==3)
    {
        for( uint32_t i{0U}; i<npe; i++)
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
                                                std::vector<double>& rSt,
                                                const std::vector<double>& xyz
                                                 )
{

    if( dim == 2 )
    {

        // no sign is taken, thus method works with cw and ccw node numbering
        double ae2 = 1.0 /( ( XY(1,0)*XY(2,1) + XY(0,0)*XY(1,1) +
                                XY(0,1)*XY(2,0) - XY(2,1)*XY(0,0) -
                                XY(2,0)*XY(1,1) - XY(1,0)*XY(0,1) ) );

        // calculation the exponents of the interpolation functions
        //const double a0 = XY(1,0) * XY(2,1) - XY(2,0) * XY(1,1);
        const double a1 = XY(2,0) * XY(0,1) - XY(0,0) * XY(2,1);
        const double a2 = XY(0,0) * XY(1,1) - XY(1,0) * XY(0,1);

        //const double b0 = XY(1,1) - XY(2,1);
        const double b1 = XY(2,1) - XY(0,1);
        const double b2 = XY(0,1) - XY(1,1);

        //const double c0 = XY(2,0) - XY(1,0);
        const double c1 = XY(0,0) - XY(2,0);
        const double c2 = XY(1,0) - XY(0,0);

        //N[0] = ae2 * (a0 + b0 * xyz[0] + c0 * xyz[1]);
        //N[1] = ae2 * (a1 + b1 * xyz[0] + c1 * xyz[1]);
        //N[2] = ae2 * (a2 + b2 * xyz[0] + c2 * xyz[1]);

        rSt.resize(2U);
        rSt[0] = ae2 * (a1 + b1 * xyz[0] + c1 * xyz[1]);
        rSt[1] = ae2 * (a2 + b2 * xyz[0] + c2 * xyz[1]);


    }else{

        //Edges of triangle 123
        const double  X12 = XY(1,0) - XY(0,0);
        const double  X23 = XY(2,0) - XY(1,0);
        const double  X31 = XY(0,0) - XY(2,0);
        const double  Y12 = XY(1,1) - XY(0,1);
        const double  Y23 = XY(2,1) - XY(1,1);
        const double  Y31 = XY(0,1) - XY(2,1);
        const double  Z12 = XY(1,2) - XY(0,2);
        const double  Z23 = XY(2,2) - XY(1,2);
        const double  Z31 = XY(0,2) - XY(2,2);

        //Vectors from point P to the points of triangle 123
        const double  XP1 = XY(0,0) - xyz[0];
        //const double  XP2 = XY(1,0) - xyz[0];
        const double  XP3 = XY(2,0) - xyz[0];
        const double  YP1 = XY(0,1) - xyz[1];
        //const double  YP2 = XY(1,1) - xyz[1];
        const double  YP3 = XY(2,1) - xyz[1];
        const double  ZP1 = XY(0,2) - xyz[2];
        //const double  ZP2 = XY(1,2) - xyz[2];
        const double  ZP3 = XY(2,2) - xyz[2];

        // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
        const double  XNRM = -Y12*Z31 + Z12*Y31;
        const double  YNRM = -Z12*X31 + X12*Z31;
        const double  ZNRM = -X12*Y31 + Y12*X31;

        // XN12,YN12,ZN12 is normal to the edge 12 in triangle plane
        const double  XN12 = Y12*ZNRM - Z12*YNRM;
        const double  YN12 = Z12*XNRM - X12*ZNRM;
        const double  ZN12 = X12*YNRM - Y12*XNRM;

        // XN23,YN23,ZN23 is normal to the edge 23 in triangle plane
        //const double  XN23 = Y23*ZNRM - Z23*YNRM;
        //const double  YN23 = Z23*XNRM - X23*ZNRM;
        //const double  ZN23 = X23*YNRM - Y23*XNRM;

        // XN31,YN31,ZN31 is normal to the edge 31 in triangle plane
        const double  XN31 = Y31*ZNRM - Z31*YNRM;
        const double  YN31 = Z31*XNRM - X31*ZNRM;
        const double  ZN31 = X31*YNRM - Y31*XNRM;

        //Scalar product of the edge and correponding normal to the neighbor edge = 2* Square of triangle
        //const double  A123 = XN23*X12 + YN23*Y12 + ZN23*Z12;
        const double  A231 = XN31*X23 + YN31*Y23 + ZN31*Z23;
        const double  A312 = XN12*X31 + YN12*Y31 + ZN12*Z31;

        //Scalar product of the vector from point P and correponding normal to the edge = 2* Square of triangle P(edge): (P23), (P31), (P12)
        //const double  AP23  = XN23*XP2 + YN23*YP2 + ZN23*ZP2;
        const double  AP31  = XN31*XP3 + YN31*YP3 + ZN31*ZP3;
        const double  AP12  = XN12*XP1 + YN12*YP1 + ZN12*ZP1;

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
        const double  X12 = XY(1,0) - XY(0,0);
        //const double  X23 = XY(2,0) - XY(1,0);
        const double  X31 = XY(0,0) - XY(2,0);
        const double  Y12 = XY(1,1) - XY(0,1);
        //const double  Y23 = XY(2,1) - XY(1,1);
        const double  Y31 = XY(0,1) - XY(2,1);
        const double  Z12 = XY(1,2) - XY(0,2);
        //const double  Z23 = XY(2,2) - XY(1,2);
        const double  Z31 = XY(0,2) - XY(2,2);

        //Vectors from point P to the points of triangle 123
        const double  XP1 = XY(0,0) - xyz[0];
        const double  XP2 = XY(1,0) - xyz[0];
        const double  XP3 = XY(2,0) - xyz[0];
        const double  YP1 = XY(0,1) - xyz[1];
        const double  YP2 = XY(1,1) - xyz[1];
        const double  YP3 = XY(2,1) - xyz[1];
        const double  ZP1 = XY(0,2) - xyz[2];
        const double  ZP2 = XY(1,2) - xyz[2];
        const double  ZP3 = XY(2,2) - xyz[2];

        // 2 * ( Area of Triangle 123 )
        const double  XN312 = -Y12*Z31 + Z12*Y31;
        const double  YN312 = -Z12*X31 + X12*Z31;
        const double  ZN312 = -X12*Y31 + Y12*X31;
        const double  A123  = sqrt(XN312*XN312 + YN312*YN312 + ZN312*ZN312);

        // 2 * ( Area of Triangle P31 )
        const double  XN3P1 = -YP1*ZP3 + ZP1*YP3;
        const double  YN3P1 = -ZP1*XP3 + XP1*ZP3;
        const double  ZN3P1 = -XP1*YP3 + YP1*XP3;
        const double  AP31  = sqrt(XN3P1*XN3P1 + YN3P1*YN3P1 + ZN3P1*ZN3P1);

        // 2 * ( Area of Triangle P12 )
        const double  XN1P2 = -YP1*ZP2 + ZP1*YP2;
        const double  YN1P2 = -ZP1*XP2 + XP1*ZP2;
        const double  ZN1P2 = -XP1*YP2 + YP1*XP2;
        const double  AP12  = sqrt(XN1P2*XN1P2 + YN1P2*YN1P2 + ZN1P2*ZN1P2);

        // 2 * ( Area of Triangle P23 )
        //const double  XN2P3 = -YP2*ZP3 + ZP2*YP3;
        //const double  YN2P3 = -ZP2*XP3 + XP2*ZP3;
        //const double  ZN2P3 = -XP2*YP3 + YP2*XP3;
        //const double  AP23  = sqrt(XN2P3*XN2P3 + YN2P3*YN2P3 + ZN2P3*ZN2P3);

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
double
IsoparametricLinearTriangle::dN( DenseMatrix<DM_MIN>& DN,
                                 const vector<double>& xyz  )
{
  vector<double> rst(dim);
  double detJ(std::numeric_limits<double>::quiet_NaN());

  PhysicalToParametric( rst, xyz );

    // here find derivatives
  dNr( rst[0], rst[1], DNR );
  dNs( rst[0], rst[1], DNS );

  if ( use2Dto3Djacobi_ ){
     vector<double> EFG(dim);
     for ( uint32_t i{0U}; i<npe; i++ ) {
            detJ = Jacobi( rst, EFG, JMAT_ );
            double det_inverse = 1.0 / ( detJ * detJ );
            DN(0,i)  = det_inverse * JMAT_(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
            DN(0,i) += det_inverse * JMAT_(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
            DN(1,i)  = det_inverse * JMAT_(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
            DN(1,i) += det_inverse * JMAT_(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
            DN(2,i)  = det_inverse * JMAT_(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
            DN(2,i) += det_inverse * JMAT_(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
      }
     return detJ;
  }

  Jacobian( DNR, DNS );
  detJ = JacobianInverse();

  DN.Resize(dim,dim);
  DN = JINV;

  DenseMatrix<DM_MIN> DN_TEMP(dim,npe);
  for( uint32_t i{0U};i<npe;i++ )
    {
      DN_TEMP(0,i)=DNR[i];
      DN_TEMP(1,i)=DNS[i];
    }
  DN *= DN_TEMP;

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
double
IsoparametricLinearTriangle::Jacobi( const vector<double>& rs,
                                     vector<double>& EFG,
                                     DenseMatrix<DM_MIN>& J )
 {
    // 1. Compute 3x2 Jacobian Matrix
    dNr( rs[0], rs[1], DNR );
    dNs( rs[0], rs[1], DNS );

    // compute 2x3 Jacobian matrix (expects 6x3 x-y-z coordinate matrix)
    J.Resize(2,dim);
    J.Zero();
    for ( uint32_t i{0U}; i<XY.Cols(); i++ )
      for ( uint32_t j{0U}; j<XY.Rows(); j++ ) {
           J(0,i) += DNR[j] * XY(j,i);
           J(1,i) += DNS[j] * XY(j,i);
        }
    // 2. Compute Jacobi J' := "determinant" of the 3x2 Jacobian
    // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)

    // compute E, F, and g
    EFG.resize(3);
    fill( EFG.begin(), EFG.end(), 0.0 );

    for ( uint32_t i{0U}; i<dim; i++ ) {
        EFG[0] += J(0,i) * J(0,i);
        EFG[1] += J(0,i) * J(1,i);
        EFG[2] += J(1,i) * J(1,i);
      }

   // compute J'
   return sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );

 }  // end Jacobi


double IsoparametricLinearTriangle::Jacobi( const vector<double>& rs )
{
    // 1. Compute 2x3 Jacobian Matrix
    dNr( rs[0], rs[1], DNR );
    dNs( rs[0], rs[1], DNS );

    // compute 2x3 Jacobian matrix
    JMAT_.Resize(2,dim);
    JMAT_.Zero();

    for ( uint32_t i{0U}; i<XY.Cols(); i++ )
      for ( uint32_t j{0U}; j<XY.Rows(); j++ ) {
           JMAT_(0,i) += DNR[j] * XY(j,i);
           JMAT_(1,i) += DNS[j] * XY(j,i);
        }
        
    // 2. Compute Jacobi J' := "determinant" of the 4x2 Jacobian
    // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)

    // compute E, F, and g
    double E(0.0), F(0.0), G(0.0);

    for ( uint32_t i{0U}; i<dim; i++ ) {
         E += JMAT_(0,i) * JMAT_(0,i);
         F += JMAT_(0,i) * JMAT_(1,i);
         G += JMAT_(1,i) * JMAT_(1,i);
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
IsoparametricLinearTriangle::N( vector<double>& N, const vector<double>& xyz )
{
   vector<double> rs(parametricDimensions_);

   PhysicalToParametric( rs, xyz );

   Nrs( rs[0], rs[1], N );

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
double
IsoparametricLinearTriangle::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B,
                                                    uint32_t gauss_point )
 {
    assert( gauss_point < gpe );

    if ( use2Dto3Djacobi_ ) {
         vector<double> EFG(dim);
         RS_[0] = rr_[gauss_point];
         RS_[1] = ss_[gauss_point];

         const double det = Jacobi( RS_, EFG, JMAT_ );
         const double det_inverse = 1.0 / ( det * det );

         B.Resize(dim,npe);
         for ( uint32_t i{0U}; i<npe; i++ )
           {
             B(0,i)  = det_inverse * JMAT_(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
             B(0,i) += det_inverse * JMAT_(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             B(1,i)  = det_inverse * JMAT_(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
             B(1,i) += det_inverse * JMAT_(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             B(2,i)  = det_inverse * JMAT_(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
             B(2,i) += det_inverse * JMAT_(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
           }
        return det;
      }

    // compute local test-function derivative matrix at gauss point
    // get local interpolation function derivatives at Gauss point
    dNr( rr_[gauss_point], ss_[gauss_point], DNR );
    dNs( rr_[gauss_point], ss_[gauss_point], DNS );

    // compute Jacobian matrix, its determinant and inverse
    Jacobian( DNR, DNS );
    const double det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2];

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
double IsoparametricLinearTriangle::dN_AtNode( DenseMatrix<DM_MIN>& B, uint32_t nd )
 {
    assert( nd < npe );
    if ( use2Dto3Djacobi_ ) {
         vector<double> EFG(dim);
         RS_[0] = NXY_(nd,0);
         RS_[1] = NXY_(nd,1);

         const double det = Jacobi( RS_, EFG, JMAT_ );
         const double det_inverse = 1. / ( det * det );

         B.Resize(dim,npe);

          for ( uint32_t i{0U}; i<npe; i++ ){
            B(0,i)  = det_inverse * JMAT_(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
            B(0,i) += det_inverse * JMAT_(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
            B(1,i)  = det_inverse * JMAT_(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
            B(1,i) += det_inverse * JMAT_(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
            B(2,i)  = det_inverse * JMAT_(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
            B(2,i) += det_inverse * JMAT_(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
          }

          return det;

     } // if


    // normal 2D case
    // --------------------------------------
    dNr( NXY_(nd,0), NXY_(nd,1), DNR );
    dNs( NXY_(nd,0), NXY_(nd,1), DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    double det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2];

    B = JINV * B;

  return det;
}


double IsoparametricLinearTriangle::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
{
  const double OneThird=1./3.0;

  if ( use2Dto3Djacobi_ ) {
       vector<double> EFG(dim);
       RS_[0] = OneThird;
       RS_[1] = OneThird;

       double det = Jacobi( RS_, EFG, JMAT_ );
       double det_inverse = 1.0 / ( det * det );

       B.Resize(dim,npe);

       for ( uint32_t i{0U}; i<npe; i++ ) {
           B(0,i)  = det_inverse * JMAT_(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
           B(0,i) += det_inverse * JMAT_(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
           B(1,i)  = det_inverse * JMAT_(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
           B(1,i) += det_inverse * JMAT_(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
           B(2,i)  = det_inverse * JMAT_(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
           B(2,i) += det_inverse * JMAT_(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
         }
       return det;
    }

    // 2D case
    // ----------------------------------------------------
    dNr( OneThird, OneThird, DNR );
    dNs( OneThird, OneThird, DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    const double det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2];

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
/*
void IsoparametricLinearTriangle::ConsecutiveNodesAtBoundary( const vector<uint32_t>& bnodes,
                                                              vector<uint32_t>& fnids )
{
   fnids.resize(3);

  // in 3D, all nodes may be at model boundary
  if ( use2Dto3Djacobi_ && bnodes.size() != 2 )
     throw csmp::Exception( FATAL_ERROR, "IsoparametricLinearTriangle::ConsecutiveNodesAtBoundary",
                            "Cannot resolve node sequence for element boundary",
                            "Probably because element lies at two boundaries simultaneously" );

  if ( bnodes.size() > 2 )
     throw csmp::Exception( ERROR, "IsoparametricLinearTriangle::ConsecutiveNodesAtBoundary",
                              "Cannot resolve node sequence for element boundary",
                              "Probably because element lies at two boundaries simultaneously" );

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

} // end ConsecutiveNodesAtBoundary
*/



void IsoparametricLinearTriangle::OutputNodeDataToVTK( const char* file_name,
                                                       const char* var_name,
                                                       DenseMatrix<DM_MIN>& DATA ) const
 {
     assert( file_name != nullptr );
     string  outfile(file_name), elmt(to_string(CurrentID()));
     outfile += elmt;
     outfile +=".vtk";

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cerr <<"\nIsoparametricLinearTriangle3D::OutputNodeDataToVTK ";
           cerr <<"Output file could not be opened."<< endl;
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
     for ( uint32_t i{0U}; i<npe; i++ ){
       if(dim==2){
        for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
           ofs<<0.0;
           ofs << endl;
        }
        else
        if(dim==3){
          for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
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
      // matrix DATA is 1 x nodes
         for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
      ofs << endl;
    }
     else {
       ofs <<"VECTORS "<< var_name <<" float"<< endl;
       // variables have always 3 components since view screen is 3D
       // matrix DATA is vec-dim x 10
       if(dim==2) {
         for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
           for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
              ofs<<0.0<<" ";
              ofs << endl;
         }
       }
       else {
          for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
            for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
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
IsoparametricLinearTriangle::ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                        const vector<double>& IVAR,
                                                              vector<double>& NVAR )
const
{
   static double a[3], b[3], c[3], intpol[3], ae2;
   static bool      first_call(true);

   vector<double>  sum(nvars);

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   //assert( gpe == 3 );
   assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

   // 1. Compute the interpolation function coefficients for the triangle which
   //    is defined by the three integration points.
   if ( first_call )
   {
        a[0] = rr_[1] * ss_[2] - rr_[2] * ss_[1];
        a[1] = rr_[2] * ss_[0] - rr_[0] * ss_[2];
        a[2] = rr_[0] * ss_[1] - rr_[1] * ss_[0];

        b[0] = ss_[1] - ss_[2];
        b[1] = ss_[2] - ss_[0];
        b[2] = ss_[0] - ss_[1];

        c[0] = rr_[2] - rr_[1];
        c[1] = rr_[0] - rr_[2];
        c[2] = rr_[1] - rr_[0];

        // compute the area*2 of the triangle spanned by the integration points
        ae2 = 1. / ( ( rr_[1]*ss_[2] + rr_[0]*ss_[1] + ss_[0]*rr_[2] - ss_[2]*rr_[0] - rr_[2]*ss_[1] - rr_[1]*ss_[0] ) );
        first_call = false;
   }

   // 2. For each node point compute the values of the interpolation functions
   //    and use these to extrapolate the values of the variables at the nodes.
   for ( uint32_t i{0U}; i<npe; i++ )
     {
        // compute interpolation function values at node i
        for ( uint32_t j{0U}; j<gpe; j++ )
          intpol[j] = ae2 * (a[j] + b[j] * NXY_(i,0) + c[j] * NXY_(i,1));

        // carry out extrapolation
        for ( uint32_t k{0U}; k<nvars; k++ ) sum[k] = 0.;
        for ( uint32_t j{0U}; j<gpe; j++ )
          for ( uint32_t k{0U}; k<nvars; k++ ) sum[k] += intpol[j] * IVAR[j*nvars + k];

        // store result in output vector
        for ( uint32_t k{0U}; k<nvars; k++ ) NVAR[i*nvars + k] = sum[k];
     }

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)



/// overwrites the base class method in order to get 2D to 3D mapping functionality
double  IsoparametricLinearTriangle::JacobianInverse()
{
    if ( use2Dto3Djacobi_ )
      {
        // Compute Jacobi J' := "determinant" of the 3x2 Jacobian
        // ------------------------------------------------------
        // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)
          vector<double> EFG(dim);
          fill( EFG.begin(), EFG.end(), 0.0 );

             for ( uint32_t i{0U}; i<dim; i++ ) {
                 EFG[0] += JAC(0,i) * JAC(0,i);
                 EFG[1] += JAC(0,i) * JAC(1,i);
                 EFG[2] += JAC(1,i) * JAC(1,i);
               }

          // compute J'
         const double detJ=sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );
         const double det1 = 1.0 / detJ;

         //Standard inverse Jacobian definition with 3rd Stroka et al. 1.0
         //JAC(2,0)=JAC(2,1)=JAC(2,2)=1.0;
         JINV.Resize( 3, 3 );
         JINV(0,0) = ( JAC(1,1) * 1.0 - JAC(1,2) * 1.0 ) *  det1;
         JINV(1,0) = ( JAC(1,0) * 1.0 - JAC(1,2) * 1.0 ) * -det1;
         JINV(2,0) = ( JAC(1,0) * 1.0 - JAC(1,1) * 1.0 ) *  det1;
         JINV(0,1) = ( JAC(0,1) * 1.0 - JAC(0,2) * 1.0 ) * -det1;
         JINV(1,1) = ( JAC(0,0) * 1.0 - JAC(0,2) * 1.0 ) *  det1;
         JINV(2,1) = ( JAC(0,0) * 1.0 - JAC(0,1) * 1.0 ) * -det1;
         JINV(0,2) = ( JAC(0,1) * JAC(1,2) - JAC(0,2) * JAC(1,1) ) *  det1;
         JINV(1,2) = ( JAC(0,0) * JAC(1,2) - JAC(0,2) * JAC(1,0) ) * -det1;
         JINV(2,2) = ( JAC(0,0) * JAC(1,1) - JAC(0,1) * JAC(1,0) ) *  det1;
         return detJ;
       }

    // the 2D case
    // compute determinant
    const double detJ = JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);

    // inversion of J
    double dum = JAC(0,0) / detJ;
    JINV.Resize( 2, 2 );
    JINV(0,0)  =  JAC(1,1) / detJ;
    JINV(0,1)  = -JAC(0,1) / detJ;
    JINV(1,0)  = -JAC(1,0) / detJ;
    JINV(1,1)  =  dum;

    if ( detJ <= 0 ) {
         cerr <<"\n\nIsoparametricLinearTriangle::JacobianInverse: Erroneous determinant of Jacobian matrix: ";
         cerr << detJ << std::endl;
         cerr <<"(are the nodes perhaps numbered clockwise?), node coordinate matrix:";
         this->XY.Out();
         string file_name{ parseFiniteElementType( ElementType() ) };
         // to print a scalar, the data matrix only needs 1 row
         DenseMatrix<DM_MIN> DATA( 1, Nodes() );
         // values increase linearly from first to last node (so that node numbering direction can be seen)
         for ( uint32_t j{0u}; j<Nodes(); ++j ) DATA(0,j) = j;
         OutputNodeDataToVTK( file_name.c_str(), "error_code", DATA );
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
vector<double>  IsoparametricLinearTriangle::UnitNormal() const
 {
    // creates a fictious third dimension perpendicular to triangle plane
    if ( dim == 2U ) return vector<double>{ 0., 0., 1. };

    const double X12 = XY(1,0) - XY(0,0), // X
                 X31 = XY(0,0) - XY(2,0),
                 Y12 = XY(1,1) - XY(0,1), // Y
                 Y31 = XY(0,1) - XY(2,1),
                 Z12 = XY(1,2) - XY(0,2), // Z
                 Z31 = XY(0,2) - XY(2,2);

    // normal to triangle (*- to flip to outside)
    vector<double> vc{ -Y12*Z31 + Z12*Y31,
                       -Z12*X31 + X12*Z31,
                       -X12*Y31 + Y12*X31 };
   
    // normalization to unit length
    double length = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
    vc[0] /= length;
    vc[1] /= length;
    vc[2] /= length;
    
    return vc;
 }




/**
    outward pointing normals to the faces that lie on the opposite nodes
    
    @author SKM 18/2/2016
    
    @test OK - for 3D version
*/
void  IsoparametricLinearTriangle::UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const
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
void  IsoparametricLinearTriangle::IntegrationPoint( uint32_t ip, vector<double>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(NXY_.Cols());
    xyz[0]=xyz[1]=0.;

    // local interpolation function values
    Nrs( rr_[ip], ss_[ip], NRST );

    // 2D
    if ( NXY_.Cols() == 2U ) {
          for( auto i{0U}; i<npe; i++ ) {
               xyz[0] += XY(i,0) * NRST[i];
               xyz[1] += XY(i,1) * NRST[i];
            }
         return;
      }

    // 3D case
    xyz[2]=0.;
    for( auto i{0U}; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint




void
IsoparametricLinearTriangle::JacobianAtIntegrationPoint( uint32_t gauss_point )
 {
    if ( !use2Dto3Djacobi_ ) {
         dNr( rr_[gauss_point], ss_[gauss_point], DNR );
         dNs( rr_[gauss_point], ss_[gauss_point], DNS );
         Jacobian( DNR, DNS );
         return;
    }

   vector<double> EFG(3);
   RS_[0] = rr_[gauss_point];
   RS_[1] = ss_[gauss_point];

   Jacobi( RS_, EFG, JAC );

 }



void
IsoparametricLinearTriangle::IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS)
 {
    IPPHYS.Resize(gpe,dim);
    vector<double> outxyz(dim);
    vector<double> rst(parametricDimensions_);

    for( uint32_t i{0U};i<gpe;i++){
        rst[0]=rr_[i];
        rst[1]=ss_[i];
        ParametricToPhysical( rst, outxyz);
        for( uint32_t j{0U};j<dim;j++) IPPHYS(i,j)=outxyz[j];
    }
}


void
IsoparametricLinearTriangle::ReferenceCoordinates( DenseMatrix<DM_MIN> & matCoords ) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);
    matCoords(0,0) = NXY_(0,0);
    matCoords(0,1) = NXY_(0,1);
    matCoords(1,0) = NXY_(1,0);
    matCoords(1,1) = NXY_(1,1);
    matCoords(2,0) = NXY_(2,0);
    matCoords(2,1) = NXY_(2,1);
}


double  IsoparametricLinearTriangle::JacobianDeterminant()
{
  if (dim == 3)
    {
      // compute E, F, and g
      std::vector<double> EFG(3);
      fill( EFG.begin(), EFG.end(), 0.0 );

      for ( uint32_t i{0U}; i<dim; i++ )
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



void IsoparametricLinearTriangle::JacobianAt( const std::vector<double>& rst )
{
    if ( !use2Dto3Djacobi_ )
    {
         dNr( rst[0], rst[1], DNR );
         dNs( rst[0], rst[1], DNS );
         Jacobian( DNR, DNS );
         return;
    }

  vector<double> EFG;
  Jacobi( rst, EFG, JAC );
}

} // end namespace csmp
