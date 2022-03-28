#include "IsoparametricQuadraticTriangle.h"
#include "ErrorHandler.h"
#include <set>
#include <fstream>
#include <cstring>
#include <climits>
#include "Exception.h"

using namespace std;

namespace csmp {

/**

Default constructor initializes triangle to use local coordinates
and a three-point integration scheme where the Gauss points are located
at the midside nodes. The weights of the integration points are set to
1/3 such that the integrated contributions add up to 1.  .
 */
IsoparametricQuadraticTriangle::IsoparametricQuadraticTriangle( uint32_t dimensions )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_TRIANGLE, true, true, 2U ),
    use2Dto3Djacobi( !(dimensions==2) ),
   DN(dimensions,6),
   BEE(dimensions,6),
   NXY(6,dimensions),
   JMAT(dimensions,dimensions),
   RS(2),
   EFG(dimensions),
   LXY(3)
 {

   dim = dimensions;
   itp = 2;
   npf = 3;
   npe = 6;
   fpe = 3;
   spe = 3;
   epe = 3;
   nne = 6;
   cne = 6;
   gpe = 3;

   // only used for the N-coefficients for interpolation functions in global
   // coordinates
   M.Resize(3,3);

   // base class matrices
   XY.Resize(npe,dim);
   JAC.Resize(dim,dim);
   JINV.Resize(dim,dim);

   BEE.Resize(dim,npe);
   JMAT.Resize(dim,dim);

   // base class vectors
   NRST.resize(npe);
   DNR.resize(npe);
   DNS.resize(npe);


   W.resize( gpe );
   rr.resize( gpe );
   ss.resize( gpe );

   // see also discussion in Smith & Griffiths 98, p. 71, Fig. 3.9, but they have a clockwise
   // node numbering
   // (fourth integration point is only used if the Gauss Radau rule is switched on)
   rr[0] = 0.5;   rr[1] = 0.5;   rr[2] = 0.0;
   ss[0] = 0.0;   ss[1] = 0.5;   ss[2] = 0.5;
   W[0]  = 1./6.; W[1]  = 1./6.; W[2]  = 1./6.;


   // initializing local node coordinates
   // nodal x-coordinates
   NXY(0,0) = 0.0;
   NXY(1,0) = 1.0;
   NXY(2,0) = 0.0;
   NXY(3,0) = 0.5;
   NXY(4,0) = 0.5;
   NXY(5,0) = 0.0;
   // nodal y-coordinates
   NXY(0,1) = 0.0;
   NXY(1,1) = 0.0;
   NXY(2,1) = 1.0;
   NXY(3,1) = 0.0;
   NXY(4,1) = 0.5;
   NXY(5,1) = 0.5;

   // intializing fixed size matrix with the interpolation function derivatives
   // at each node, see for instance Akin, p. 113
   // verified with Mathematica
   DN(0,0) = -3.0; DN(1,0) = -3.0; // w.r.t. x
   DN(0,1) =  3.0; DN(1,1) =  0.0;
   DN(0,2) =  0.0; DN(1,2) =  3.0;
   DN(0,3) =  0.0; DN(1,3) = -2.0;
   DN(0,4) =  2.0; DN(1,4) =  2.0;
   DN(0,5) = -2.0; DN(1,5) =  0.0;

   UsesLocalCoordinates(true);
   Isoparametric(true);
   SurfaceElement();
   ElementType(ISOPARAMETRIC_QUADRATIC_TRIANGLE);
}



IsoparametricQuadraticTriangle::~IsoparametricQuadraticTriangle()
 {
 }


void IsoparametricQuadraticTriangle::Dimensions( uint32_t dimensions )
 {
   dim = dimensions;
   if ( dim == 3 ) use2Dto3Djacobi = true;
   else            use2Dto3Djacobi = false;
   // base class matrices
   XY.Resize(npe,dim);
   JAC.Resize(dim,dim);
   JINV.Resize(dim,dim);
   BEE.Resize(dim,npe);
   JMAT.Resize(dim,dim);
 }



/**

While the default Gauss points are located at the midside nodes,
this function allows to switch them to the distance of 1/6 near
the corner nodes.
*/
void IsoparametricQuadraticTriangle::GaussPointsNearCorners()
 {
    gpe = 3;
    W.resize( gpe );
    rr.resize( gpe );
    ss.resize( gpe );
    // integration points a la Jean Braun
    rr[0] = 1./6.; rr[1] = 2./3.; rr[2] = 1./6.;
    ss[0] = 1./6.; ss[1] = 1./6.; ss[2] = 2./3.;
    W[0]  = 1./6.; W[1]  = 1./6.; W[2]  = 1./6.;
 }



/**

Sets the integration scheme to an asymmetric four point Gauss-Radau rule
with variable weights depending on the location of the Gauss points. This
integration scheme is described for the unit triangle in Akin 1982,
Fig. 6.6a, p. 103.
 */
void IsoparametricQuadraticTriangle::GaussRadau4PointIntegration()
 {
    gpe = 4;
    W.resize( gpe );
    rr.resize( gpe );
    ss.resize( gpe );
   // Gauss Radau 4 point quadrature rule for unit triangle Fig. 6.6a, p. 103, Akin 1982.
   rr[0] = 0.2800199155; rr[1] = 0.6663902460; rr[2] = 0.0750311102; rr[3] = 0.1785587283;
   ss[0] = 0.0750311102; ss[1] = 0.1785587283; ss[2] = 0.2800199155; ss[3] = 0.6663902460;
   W[0]  = 0.0909793091;  W[1] = 0.1590206909;  W[2] = 0.0909793091;  W[3] = 0.1590206909;
 }




/**

Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment. This rule differs from that presented
for the straight-sided triangle which is is discussed in the
LinearTriangle subclass.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double  IsoparametricQuadraticTriangle::AspectRatio()
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

This method attempts to calculate the inner radius of the curved-sided
triangle by talking 1/2 of the perimeter and dividing this measure by
the area of the triangle.

The procedure is empirically based giving a good match if the triangles
are relatively even-sided and have straight edges.

@section arguments Input Arguments

The parent element is queried for its node coordinates.
*/
double  IsoparametricQuadraticTriangle::InnerRadius()
{
   vector<double> segms(npe);
   double         sum(0.0), vol;

   EdgeLengths( segms );
   for ( uint32_t i{0U}; i<segms.size(); i++ ) sum += segms[i];
   sum /= 2.0;
   vol  = Volume();
   vol /= sum;
   return vol;
}



/**

Segment length is approximated as the sum of the length of two linear
segments making up each face of the quadratic triangle. The first segment
connects corner node 1 with midpoint node 3 and midpoint node 3 with
corner node 2 and so forth.

@param len The lengths of the segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as the sum of the straight line
segments which make up the element face.

@section application Application

For example, when stresses are to be applied at the element side, the
(area=length at unit thickness) must be taken into account.

*/
void  IsoparametricQuadraticTriangle::EdgeLengths( std::vector<double>& len )
{
    double sum;
    len.resize(spe);

    // segment 1
    sum     = (XY(3,0)-XY(0,0)) * (XY(3,0)-XY(0,0));
    sum    += (XY(3,1)-XY(0,1)) * (XY(3,1)-XY(0,1));
    len[0]  = sqrt(sum);
    sum     = (XY(1,0)-XY(3,0)) * (XY(1,0)-XY(3,0));
    sum    += (XY(1,1)-XY(3,1)) * (XY(1,1)-XY(3,1));
    len[0] += sqrt(sum);

    // segment 2
    sum     = (XY(2,0)-XY(4,0)) * (XY(2,0)-XY(4,0));
    sum    += (XY(2,1)-XY(4,1)) * (XY(2,1)-XY(4,1));
    len[1]  = sqrt(sum);
    sum     = (XY(4,0)-XY(1,0)) * (XY(4,0)-XY(1,0));
    sum    += (XY(4,1)-XY(1,1)) * (XY(4,1)-XY(1,1));
    len[1] += sqrt(sum);

    // segment 3
    sum     = (XY(0,0)-XY(5,0)) * (XY(0,0)-XY(5,0));
    sum    += (XY(0,1)-XY(5,1)) * (XY(0,1)-XY(5,1));
    len[2]  = sqrt(sum);
    sum     = (XY(5,0)-XY(2,0)) * (XY(5,0)-XY(2,0));
    sum    += (XY(5,1)-XY(2,1)) * (XY(5,1)-XY(2,1));
    len[2] += sqrt(sum);

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
double  IsoparametricQuadraticTriangle::Volume()
{
    double  area{0.};

    // numerical integration:
    // looping over the 3 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( auto i{0U}; i<gpe; i++ )
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

  void IsoparametricQuadraticTriangle::Nrs( double r, double s, vector<double>& N ) const



Description:

Returns the value of the element interpolation functions at the point 'rs'
in local coordinates.

@section arguments Input Arguments

The floating point coordinates 'r' and 's' of the point at which the
interpolation shall be carried out.

@param nrst The third method argument is the vector into which the values of the
n=nodes interpolation functions at the point 'rs' will be returned.

@section implementation Implementation

See in source code header file.

@section application Application

Method is used to compute property values at the integration points of
the element.
*/
void IsoparametricQuadraticTriangle::Nrs( double r, double s, std::vector<double>& nrst ) const
{
   nrst.resize(npe);
   nrst[3] = 4. * r * (1. - r - s);
   nrst[4] = 4. * r * s;
   nrst[5] = 4. * s * (1. - r - s);
   nrst[0] = 1. - r - s - nrst[3]/2. - nrst[5]/2.;
   nrst[1] = r - nrst[3]/2. - nrst[4]/2.;
   nrst[2] = s - nrst[4]/2. - nrst[5]/2.;
}

void IsoparametricQuadraticTriangle::Nrs( double r, double s, double* nrst ) const
{
   nrst[3] = 4. * r * (1. - r - s);
   nrst[4] = 4. * r * s;
   nrst[5] = 4. * s * (1. - r - s);
   nrst[0] = 1. - r - s - nrst[3]/2. - nrst[5]/2.;
   nrst[1] = r - nrst[3]/2. - nrst[4]/2.;
   nrst[2] = s - nrst[4]/2. - nrst[5]/2.;
}


/**

  void IsoparametricQuadraticTriangle::N_AtIntegrationPoint( uint32_t ip, vector<double>& N )



Description:

Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

@param N The interpolation function values are returned into the third argument.

*/
void IsoparametricQuadraticTriangle::N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N )
 {
    assert( ip < gpe );

    N.resize(npe);
    // local interpolation function values
    Nrs( rr[ip], ss[ip], N );

 } // end N_AtIntegrationPoint




/**

This method and the complementary method dNs() compute the shape function
derivates with respect to the local coordinate axis 'r', and 's'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNR The shape function derivatives are returned into the second method
argument which is a Meschpp floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.

*/
void IsoparametricQuadraticTriangle::dNr( double s, double t, std::vector<double>& DNR ) const
{
   DNR.resize(npe);
   DNR[0] = -3. + 4.*s + 4.*t;
   DNR[1] = -1. + 4.*s;
   DNR[2] =  0.;
   DNR[3] =  4. - 8.*s - 4.*t;
   DNR[4] =  4.*t;
   DNR[5] = -4.*t;
}


void IsoparametricQuadraticTriangle::dNs(  double s, double t, std::vector<double>& DNS ) const
{
   DNS.resize(npe);
   DNS[0] = -3. + 4.*s + 4.*t;
   DNS[1] =  0.;
   DNS[2] = -1. + 4.*t;
   DNS[3] = -4.*s;
   DNS[4] =  4.*s;
   DNS[5] =  4. - 4.*s - 8.*t;
}








/**

void IsoparametricQuadraticTriangle::NodesOfFace( uint32_t face_id,
                                     vector<uint32_t>& fnids )



Description:

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
void IsoparametricQuadraticTriangle::NodesOfFace( uint32_t face_id,
                                                  std::vector<uint32_t>& fnids ) const
 {
    fnids.resize(3);

    if ( face_id == 0 )
      {
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 4;
      }
    else if ( face_id == 1 )
      {
         fnids[0] = 2;
         fnids[1] = 0;
         fnids[2] = 5;
      }
    else if ( face_id == 2 )
      {
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 3;
      }
    else
    std::cerr <<"\nIsoparametricQuadraticTriangle::NodesOfFace: Erratic input face ID: "<< face_id << std::endl;
 }



/// returns the local  numbers of the nodes at the other end of the sgment that the argument node is on
std::vector<uint32_t>  IsoparametricQuadraticTriangle::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nIsoparametricQuadraticTriangle::NodesConnectedTo: node "<< node_id <<" does not exist.";
    return vector<uint32_t>{};
  }



vector<uint32_t>  IsoparametricQuadraticTriangle::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nIsoparametricQuadraticTriangle::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }



void IsoparametricQuadraticTriangle::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
 {
    snids.resize(3);
    if ( segm_id == 0 ) {
         snids[0] = 1;
         snids[1] = 2;
         snids[2] = 4;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 0;
         snids[1] = 2;
         snids[2] = 5;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 0;
         snids[1] = 1;
         snids[2] = 3;
      }
    else
    std::cerr <<"\nIsoparametricQuadraticTriangle::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment





/**

Returns local node ids of the nodes which sit on the corner of
the quadratic triagular element.

@param ids Returns an integer vector with the 3 local corner node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).

*/
void  IsoparametricQuadraticTriangle::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(3);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }




/**

Returns local node ids of the nodes which sit on the midsides of
the quadratic triagular element.

@param ids Returns an integer vector with the 3 local midside-node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).

*/
void  IsoparametricQuadraticTriangle::MidSideNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npf);
    ids[0] = 3;
    ids[1] = 4;
    ids[2] = 5;
 }



void  IsoparametricQuadraticTriangle::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 3;
    ids[2] = 1;
    ids[3] = 4;
    ids[4] = 2;
    ids[5] = 5;
 }

double IsoparametricQuadraticTriangle::WeightAtIntegrationPoint( uint32_t i ) const { return W[i]; }



CSMP_FEM_TYPE IsoparametricQuadraticTriangle::ElementTypeOfFace( uint32_t ) const
 {
    return ISOPARAMETRIC_QUADRATIC_BAR;
 }





void  IsoparametricQuadraticTriangle::JacobianAtIntegrationPoint( uint32_t gauss_point )
 {
    if ( !use2Dto3Djacobi ) {
         dNr( rr[gauss_point], ss[gauss_point], DNR );
         dNs( rr[gauss_point], ss[gauss_point], DNS );
         Jacobian( DNR, DNS );
         return;
      }

    RS[0] = rr[gauss_point];
    RS[1] = ss[gauss_point];

    Jacobi( RS, EFG, JAC );
 }


/**
     Assumes that RS, EFG and the local JAC matrix were initialised before.
 */
double IsoparametricQuadraticTriangle::JacobianDeterminant()
  {
        if ( !use2Dto3Djacobi )
          return FiniteElement::JacobianDeterminant();
    
        return Jacobi( RS, EFG, JAC );
    }
    
    
    

void  IsoparametricQuadraticTriangle::N_AtBaryCenter( std::vector<double>& N )
 {
    N.resize(npe);
    N[0] = -1./9.;
    N[1] =  N[0];
    N[2] =  N[0];
    N[3] =  4./9.;
    N[4] =  N[3];
    N[5] =  N[3];
 }



/**

Compute derivatives of interpolation functions at corresponding nodes with
respect to the global coordinate system. If a transformation from 2D local
to 3D global space is applied, the function Jacobi is called.

@section arguments Input Arguments

The element is used to obtain the global interpolation of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN6 The interpolation-function derivative matrix is returned into the
second method argument.

 */
void IsoparametricQuadraticTriangle::dN( DenseMatrix<DM_MIN>& DN6 )
  {
     DN6.Resize(dim,npe);

     double detJ;

     // Jacobian transformation to global coordinate system
     if ( use2Dto3Djacobi )
       {
          for ( uint32_t i{0U}; i<npe; i++ ) {
                RS[0] = NXY(i,0);
                RS[1] = NXY(i,1);

                // generate the 3 x 6 shape function derivative matrix
                detJ = Jacobi( RS, EFG, JMAT );
                double det_inverse = 1. / ( detJ * detJ );
                DN6(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN6(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN6(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN6(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN6(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN6(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             }
           return;
        }

     // if the quadratic triangle is 2D
     DenseMatrix<DM_MIN> TEMP(dim,1);
     DN6 = DN;

     for ( uint32_t i{0U}; i<npe; i++ )
       {
           dNr( NXY(i,0), NXY(i,1), DNR );
           dNs( NXY(i,0), NXY(i,1), DNS );
           Jacobian( DNR, DNS );
           JacobianInverse();

           TEMP(0,0) = DN6(0,i);
           TEMP(1,0) = DN6(1,i);

           // 2x2 * 2x1 = 2x1 gives the global DN entries
           JMAT.Resize(dim,dim);
           JMAT = JINV;
           JMAT *= TEMP;
           DN6(0,i) = JMAT(0,0);
           DN6(1,i) = JMAT(1,0);
       }

  } // end dN





/** Computes derivatives of interpolation functions at a point XY given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global interpolation of the triangle and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.
*/
double IsoparametricQuadraticTriangle::dN_At( DenseMatrix<DM_MIN>& mDN, const vector<double>& xyz  )
 {
    if ( dim == 3 )
      throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticTriangle::dN",
                      "Function can only be used if dim = 2");

    mDN.Resize(dim,npe);
    // update the interpolation function coefficients
    //
    //      | ai aj ak |
    //      | bi bj bk |  = coeffs
    //      | ci cj ck |
    //

    DNR.resize(dim);

    PhysicalToParametric(DNR,xyz);

    LXY[0]=DNR[2];
    LXY[1]=DNR[0];
    LXY[2]=DNR[1];

    // compute 1 / (2 x area) assuming that the triangle is straight sided
    double area2 = 1.0 / ( XY(1,0)*XY(2,1) + XY(0,0)*XY(1,1) +
                             XY(0,1)*XY(2,0) - XY(2,1)*XY(0,0) -
                             XY(2,0)*XY(1,1) - XY(1,0)*XY(0,1) );

    // compute interpolation function derivatives wrt x as given in Cook page 158
    mDN(0,0) = (4.0 * LXY[0] - 1.0 ) * M(1,0);
    mDN(0,1) = (4.0 * LXY[1] - 1.0 ) * M(1,1);
    mDN(0,2) = (4.0 * LXY[2] - 1.0 ) * M(1,2);
    mDN(0,3) = 4.0 * ( LXY[1] * M(1,0) + LXY[0] * M(1,1) );
    mDN(0,4) = 4.0 * ( LXY[2] * M(1,1) + LXY[1] * M(1,2) );
    mDN(0,5) = 4.0 * ( LXY[0] * M(1,2) + LXY[2] * M(1,0) );

    // compute interpolation function derivatives wrt y as given in Cook page 158
    mDN(1,0) = (4.0 * LXY[0] - 1.0 ) * M(2,0);
    mDN(1,1) = (4.0 * LXY[1] - 1.0 ) * M(2,1);
    mDN(1,2) = (4.0 * LXY[2] - 1.0 ) * M(2,2);
    mDN(1,3) = 4.0 * ( LXY[1] * M(2,0) + LXY[0] * M(2,1) );
    mDN(1,4) = 4.0 * ( LXY[2] * M(2,1) + LXY[1] * M(2,2) );
    mDN(1,5) = 4.0 * ( LXY[0] * M(2,2) + LXY[2] * M(2,0) );

    return 1.0 / ( 2.0 * area2 );

  } // end dN











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
double IsoparametricQuadraticTriangle::Jacobi( const vector<double>& rs, vector<double>& vEFG, DenseMatrix<DM_MIN>& J )
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
    vEFG.resize(3);
    fill( vEFG.begin(), vEFG.end(), 0. );

    for ( uint32_t i{0U}; i<dim; i++ ) {
        vEFG[0] += J(0,i) * J(0,i);
        vEFG[1] += J(0,i) * J(1,i);
        vEFG[2] += J(1,i) * J(1,i);
      }

   // compute J'
   return sqrt( vEFG[0] * vEFG[2] - vEFG[1] * vEFG[1] );

 }  // end Jacobi





double IsoparametricQuadraticTriangle::Jacobi( const vector<double>& rs )
 {
    // 1. Compute 2x3 Jacobian Matrix
    dNr( rs[0], rs[1], DNR );
    dNs( rs[0], rs[1], DNS );

    // compute 2x3 Jacobian matrix
    JMAT.Resize(2,dim);
    JMAT.Zero();

    for ( uint32_t i{0U}; i<XY.Cols(); i++ )
      for ( uint32_t j{0U}; j<XY.Rows(); j++ ) {
           JMAT(0,i) += DNR[j] * XY(j,i);
           JMAT(1,i) += DNS[j] * XY(j,i);
        }

    // 2. Compute Jacobi J' := "determinant" of the 3x2 Jacobian
    // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)

    // compute E, F, and g
    double E(0.0), F(0.0), G(0.0);

    for ( uint32_t i{0U}; i<dim; i++ ) {
         E += JMAT(0,i) * JMAT(0,i);
         F += JMAT(0,i) * JMAT(1,i);
         G += JMAT(1,i) * JMAT(1,i);
      }

   // compute J'
   return  sqrt( E * G - F * F );

 }  // end Jacobi


/// Projection function from rs->xy(z)
void
IsoparametricQuadraticTriangle::ParametricToPhysical( std::vector<double> &rst,
                                                      std::vector<double>& xyz)
{
  vector<double> N(npe);
  Nrs(rst[0],rst[1], N );

  for( uint32_t i{0U}; i<dim; i++) xyz[i]=0.0;

    if(dim==2){
      for( uint32_t i{0U}; i<npe; i++){
        xyz[0]+=XY(i,0)*N[i];
        xyz[1]+=XY(i,1)*N[i];
       }
    }else if(dim==3){
        for( uint32_t i{0U}; i<npe; i++){
          xyz[0]+=XY(i,0)*N[i];
          xyz[1]+=XY(i,1)*N[i];
          xyz[2]+=XY(i,2)*N[i];
        }
      }
}

/** Projection function from xyz->rst


In first part over regular grid looks for closest point to the given one;
Then the point is used as initial guess for a Newton-Pahson iterration with
constant Jacobian matrix
*/
void IsoparametricQuadraticTriangle::PhysicalToParametric( std::vector<double>& rSt,
                                                           const std::vector<double>& xyz )
{
    vector<double> outxyz(dim);
    vector<double> rstHatK(2U);

    std::vector<double> distanceFromGivenPointLinf(3U,0.0);
    double distanceFromGivenPointL2;

    // Find largest and smallest segments in order to define precision
    vector<double> vec(spe);
    EdgeLengths( vec );
    double seg_max(vec[0]), seg_min(vec[0]);
    for ( uint32_t i=1; i<spe; i++ )
    {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
    }

    const double geometricTolerance = 0.005*seg_min;

    // First guess as BaryCenter
    rstHatK[0] = 1.0/3.0;
    rstHatK[1] = 1.0/3.0;

    ParametricToPhysical( rstHatK, outxyz);

    if( dim == 3 )
    {
        distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
        distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
        distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

        distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                         distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                         distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );
    }else{

        distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
        distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );

        distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                         distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] );
    }

    if( distanceFromGivenPointL2 > geometricTolerance )
    //if( (distanceFromGivenPointLinf[0] > geometricTolerance) ||
    //    (distanceFromGivenPointLinf[1] > geometricTolerance) ||
    //    (distanceFromGivenPointLinf[2] > geometricTolerance)  )
    {

        vector<double> rstHatK_PlusOne(2U);
        double minDistanceFromGivenPoint;

        const double constantMu               = 1.0;
        const uint32_t numberOfFirstIterrations   = 5;
        const uint32_t maxNumberOfIterrations     = 20;
        const uint32_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double incrementR = 2.0/(numberOfFirstIterrations-1);
        const double incrementS = 2.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = -1.0;
        for( uint32_t i{0U};i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = -1.0;
            for( uint32_t j{0U};j<numberOfFirstIterrations;j++)
            {
                ParametricToPhysical( rstHatK_PlusOne, outxyz);

                if( dim == 3 )
                {
                    distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
                    distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
                    distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

                    distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                                     distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                                     distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );
                }else{

                    distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
                    distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );

                    distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                                     distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] );
                }

                if( minDistanceFromGivenPoint > distanceFromGivenPointL2 )
                {
                    for( uint32_t l=0; l<2U; l++)
                        rstHatK[l] = rstHatK_PlusOne[l];

                    minDistanceFromGivenPoint = distanceFromGivenPointL2;
                }
                rstHatK_PlusOne[1] += incrementS;
            }
            rstHatK_PlusOne[0] += incrementR;
        }

        ParametricToPhysical( rstHatK, outxyz);

        if( dim == 3 )
        {
            distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
            distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
            distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

            distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                             distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                             distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );
        }else{

            distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
            distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );

            distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                             distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] );
        }

        uint32_t iteration = 1;

        /// Newton-Raphson iterations
        while ( ( distanceFromGivenPointL2 > geometricTolerance ) && ( iteration < maxNumberOfIterrations ) )
        //while ( ( ( distanceFromGivenPointLinf[0] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[1] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[2] > geometricTolerance )    )
        //      &&  ( iteration < maxNumberOfIterrations )                )
        {
            // Out of range check
            if( (rstHatK[0]<0.) || (rstHatK[0]>1.) ||
                (rstHatK[1]<0.) || (rstHatK[1]>1.)  )
            {
                rstHatK[2] = 0.0;
                rstHatK[2] = 2.0;

                break;
            }

            // Inverse Jacobian calculations
            if( iteration < numberOfIterationsWhenJacobiIsNotConstant )
            {
                dNr( rstHatK[0], rstHatK[1], DNR );
                dNs( rstHatK[0], rstHatK[1], DNS );
                Jacobian( DNR, DNS );
                // Check whether Jacobian is positive ( might be not true for the point outside the element )
                //const double detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            if ( dim == 3 )
            {
                rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
                rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
                rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

                for( uint32_t i{0U}; i<2; i++)
                    rstHatK[i] = rstHatK_PlusOne[i];

                ParametricToPhysical( rstHatK, outxyz);
                distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
                distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
                distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );
                distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                                 distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                                 distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );
            }else{

                rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) );
                rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) );

                for( uint32_t i{0U}; i<2; i++)
                    rstHatK[i] = rstHatK_PlusOne[i];

                ParametricToPhysical( rstHatK, outxyz);
                distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
                distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
                distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                                 distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] );

            }

            iteration++;

        }

        JAC.Zero();

        if( iteration == maxNumberOfIterrations )
        {
            cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            if( dim ==3 )
            {
                cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Real Point:\t"
                    <<"x = "<<xyz[0]<<" ;\t"
                    <<"y = "<<xyz[1]<<" ;\t"
                    <<"z = "<<xyz[2]<<"\n";
                cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Found Point:\t"
                   <<"x = "<<outxyz[0]<<" ;\t"
                   <<"y = "<<outxyz[1]<<" ;\t"
                   <<"z = "<<outxyz[2]<<"\n";
                cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Distance:\t"
                   <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
                   <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
                   <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            }else{
                cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Real Point:\t"
                    <<"x = "<<xyz[0]<<" ;\t"
                    <<"y = "<<xyz[1]<<"\n";
                cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Found Point:\t"
                    <<"x = "<<outxyz[0]<<" ;\t"
                    <<"y = "<<outxyz[1]<<"\n";
                cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Distance:\t"
                    <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
                    <<"y = "<<distanceFromGivenPointLinf[1]<<endl;
            }

            cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<"\n";

            std::vector<double> N(npe,0.0);
            Nrs(rstHatK[0], rstHatK[1], N );
            cout<<" IsoparametricQuadraticTriangle::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<" ;\t"
                <<"N[5] = "<<N[5]<<"\n";

            csmp::Exception( WARNING, "IsoparametricQuadraticTriangle::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for( uint32_t i{0U}; i<2; i++)
        rSt[i] = rstHatK[i];
}



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

@param N The method argument, the vector FN, will hold the interpolation function
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
void IsoparametricQuadraticTriangle::N( vector<double>& N, const vector<double>& xyz )
 {
    vector<double> rs(parametricDimensions);

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

@param B The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.

*/
double IsoparametricQuadraticTriangle::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, uint32_t gauss_point )
 {
    double det;

    assert( gauss_point < gpe );

    if ( use2Dto3Djacobi )
      {
         vector<double> vEFG(dim);
         RS[0] = rr[gauss_point];
         RS[1] = ss[gauss_point];

         det = Jacobi( RS, vEFG, JMAT );
         double det_inverse = 1.0 / ( det * det );

         B.Resize(dim,npe);

         for ( uint32_t i{0U}; i<npe; i++ )
           {
             B(0,i)  = det_inverse * JMAT(0,0) * ( vEFG[2] * DNR[i] - vEFG[1] * DNS[i] );
             B(0,i) += det_inverse * JMAT(1,0) * ( vEFG[0] * DNS[i] - vEFG[1] * DNR[i] );
             B(1,i)  = det_inverse * JMAT(0,1) * ( vEFG[2] * DNR[i] - vEFG[1] * DNS[i] );
             B(1,i) += det_inverse * JMAT(1,1) * ( vEFG[0] * DNS[i] - vEFG[1] * DNR[i] );
             B(2,i)  = det_inverse * JMAT(0,2) * ( vEFG[2] * DNR[i] - vEFG[1] * DNS[i] );
             B(2,i) += det_inverse * JMAT(1,2) * ( vEFG[0] * DNS[i] - vEFG[1] * DNR[i] );
           }

         return det;
      }

    // compute local test-function derivative matrix at gauss point
    // get local interpolation function derivatives at Gauss point
    dNr( rr[gauss_point], ss[gauss_point], DNR );
    dNs( rr[gauss_point], ss[gauss_point], DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2];
    B(0,3) = DNR[3]; B(1,3) = DNS[3];
    B(0,4) = DNR[4]; B(1,4) = DNS[4];
    B(0,5) = DNR[5]; B(1,5) = DNS[5];

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
double IsoparametricQuadraticTriangle::dN_AtNode( DenseMatrix<DM_MIN>& B, uint32_t nd )
 {
    double  det;

    assert( nd < npe );

    if ( use2Dto3Djacobi )
      {
        RS[0] = NXY(nd,0);
        RS[1] = NXY(nd,1);

        det = Jacobi( RS, EFG, JMAT );
        double det_inverse = 1.0 / ( det * det );

        B.Resize(dim,npe);

        for ( uint32_t i{0U}; i<npe; i++ )
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


    // normal 2D case
    // --------------------------------------
    dNr( NXY(nd,0), NXY(nd,1), DNR );
    dNs( NXY(nd,0), NXY(nd,1), DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2];
    B(0,3) = DNR[3]; B(1,3) = DNS[3];
    B(0,4) = DNR[4]; B(1,4) = DNS[4];
    B(0,5) = DNR[5]; B(1,5) = DNS[5];

    B = JINV * B;

    return det;
 }




double IsoparametricQuadraticTriangle::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    double  one_third(1.0/3.0), det;

    if ( use2Dto3Djacobi )
      {
        RS[0] = one_third;
        RS[1] = one_third;

        det = Jacobi( RS, EFG, JMAT );
        double det_inverse = 1.0 / ( det * det );

        B.Resize(dim,npe);

        for ( uint32_t i{0U}; i<npe; i++ )
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

    // 2D case
    // ----------------------------------------------------
    dNr( one_third, one_third, DNR );
    dNs( one_third, one_third, DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2];
    B(0,3) = DNR[3]; B(1,3) = DNS[3];
    B(0,4) = DNR[4]; B(1,4) = DNS[4];
    B(0,5) = DNR[5]; B(1,5) = DNS[5];

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
void  IsoparametricQuadraticTriangle::ConsecutiveNodesAtBoundary( const vector<uint32_t>& bnodes,
                                                                  vector<uint32_t>& fnids )
 {
     fnids.resize(3);

     // in 3D, all nodes may be at model boundary
     if ( use2Dto3Djacobi && bnodes.size() != 3 )
       throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticTriangle::ConsecutiveNodesAtBoundary",
                              "Cannot resolve node sequence for element boundary",
                              "Probably because element lies at two boundaries simultaneously" );

     if ( bnodes.size() > 3 )
       throw csmp::Exception( ERROR, "IsoparametricQuadraticTriangle::ConsecutiveNodesAtBoundary",
                              "Cannot resolve node sequence for element boundary",
                              "Probably because element lies at two boundaries simultaneously" );

     // the midside nodes determine the boundaries so they are sought after
     if ( bnodes[0] == 3 || bnodes[1] == 3 || bnodes[2] == 3 ) {
          fnids[0] = 0; fnids[1] = 3; fnids[2] = 1;
          return;
       }
     if ( bnodes[0] == 4 || bnodes[1] == 4 || bnodes[2] == 4 ) {
          fnids[0] = 1; fnids[1] = 4; fnids[2] = 2;
          return;
       }
     if ( bnodes[0] == 5 || bnodes[1] == 5 || bnodes[2] == 5 ) {
          fnids[0] = 2; fnids[1] = 5; fnids[2] = 0;
          return;
       }

 } // end ConsecutiveNodesAtBoundary










/** Outputs nodal values of a triangle into user-specified output file.


@section arguments Input Arguments

The parent element, the output file name to which the extension ".vtk"
will be appended automatically, the name of the variable,
a data matrix which will contain n-colums=nodes and n-rows = dimensions
of data (n-rows=1 = scalar data, n-rows=3 = vector data etc.).

@section application Application

Single elements are output to VTK as polygons in order to test the
quality of interpolation and the computed properties directly.

@section messages Messages

The method will indicate if there is a problem in opening the output
file.

@attention SKM fixed for quadratic triangle VTK cell type (=22)

 */
void IsoparametricQuadraticTriangle::OutputNodeDataToVTK( const char* file_name,
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
           cerr <<"\nQuadraticTriangle::OutputNodeDataToVTK ";
           cerr <<"Output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 4.2"<< endl;
     ofs <<"Finite-element dataset (CSMP++): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" double"<< endl;
     for ( uint32_t i{0U}; i<6; i++ ) {
         if ( use2Dto3Djacobi ) {
             for ( uint32_t j{0U}; j<3; j++ ) ofs << COORD(i,j) <<" ";
             ofs << endl;
           }
         else {
             for ( uint32_t j{0U}; j<2; j++ ) ofs << COORD(i,j) <<" ";
             ofs << 0.0 << endl;
           }
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< npe+1 << endl;
     // number of points per cell, point1, ... point n
//     ofs << 6 <<" 0 3 1 4 2 5";
     ofs << 6 <<" 0 1 2 3 4 5";
     ofs << endl;
     ofs << endl;

     // 4. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 22 << endl; // VTK_QUADRATIC_TRIANGLE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     ofs.setf( ios::scientific );

     ofs <<"POINT_DATA "<< npe << endl;
// if cell data shall be output, there is just a single value per element:  ofs <<"CELL_DATA "<< 1 << endl;

     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" double"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x6
           for ( uint32_t i{0U}; i<6; i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" double"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 6
          for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
               for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               if ( use2Dto3Djacobi )
                 ofs << endl;
               else
                 ofs << 0.0 << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nQuadraticTriangle::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK




/**

Because, in the local coordinate system, the triangle is always straight
sided, and the derivatives of the interpolation functions are bilinear,
variable values computed at integration points can simply be extrapolated
to the nodes using an approach as for the linear triangular element. This
involves the steps:

0. Decide which case is dealt with in terms of the integration points
   which are used.

1. Compute the interpolation function coefficients for the triangle which
   is defined by the three integration points.

2. For each node point compute the values of the interpolation functions
   and use these to extrapolate the values of the variables at the nodes.
 */
void IsoparametricQuadraticTriangle::ExtrapolateIntegrationPointVariableToNodes(
                                                   uint32_t nvars,
                                                   const vector<double>& IVAR,
                                                   vector<double>& NVAR )
const
{
   static double  a[3], b[3], c[3], intpol[3], ae2;
   static bool       first_call(true);

   vector<double>  sum(nvars);

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   assert( gpe == 3 );
   assert( IVAR.size() >= (gpe*nvars) );
   NVAR.resize( npe * nvars );

   // 1. Compute the interpolation function coefficients for the triangle which
   //    is defined by the three integration points.
   if ( first_call ) {
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
   for ( uint32_t i{0U}; i<npe; i++ )
     {
        // compute interpolation function values at node i
        for ( uint32_t j{0U}; j<gpe; j++ )
          intpol[j] = ae2 * (a[j] + b[j] * NXY(i,0) + c[j] * NXY(i,1));

        // carry out extrapolation
        for ( uint32_t k{0U}; k<nvars; k++ ) sum[k] = 0.0;
        for ( uint32_t j{0U}; j<gpe; j++ )
          for ( uint32_t k{0U}; k<nvars; k++ ) sum[k] += intpol[j] * IVAR[j*nvars + k];

        // store result in output vector
        for ( uint32_t k{0U}; k<nvars; k++ ) NVAR[i*nvars + k] = sum[k];
     }

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)








/**
    Specific subclass method that implements 2D-to-3D Jacobian mapping functionality
    that is needed to represent the quadratic triangle if it is placed in the 
    three-dimensional model.
    
    @attention if a negative Jacobian is detected, method provides additional diagnostics to 
    locate the problem in the mesh.
*/
double  IsoparametricQuadraticTriangle::JacobianInverse()
 {
    if ( use2Dto3Djacobi )
      {
         // Compute Jacobi J' := "determinant" of the 3x2 Jacobian
         // ------------------------------------------------------
         // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)
         static vector<double> EFG(dim);
           fill( EFG.begin(), EFG.end(), 0.0 );

           for ( uint32_t i{0U}; i<dim; i++ ) {
               EFG[0] += JAC(0,i) * JAC(0,i);
               EFG[1] += JAC(0,i) * JAC(1,i);
               EFG[2] += JAC(1,i) * JAC(1,i);
             }

         // compute J'
         return sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );
      }

    // the 2D case
    // compute determinant
    double detJ = JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);

    // inversion of J
    double dum = JAC(0,0) / detJ;
    JINV(0,0)  =  JAC(1,1) / detJ;
    JINV(0,1)  = -JAC(0,1) / detJ;
    JINV(1,0)  = -JAC(1,0) / detJ;
    JINV(1,1)  =  dum;

    if ( detJ <= 0. ) {
       ErrorHandler&  csmp_error( ErrorHandler::Instance() );
       cerr <<"\nIsoparametricQuadraticTriangle::JacobianInverse: element "<< CurrentID() <<": erroneous determinant of 2D Jacobian matrix: ";
        cerr << detJ << endl;
        cerr <<"\ncaused by element of type: "<< parseFiniteElementType(ElementType()) << endl;
        for ( uint32_t i{0U}; i<Nodes(); i++ )
          {
            cerr<<" Node( "<<i<<" ): "<<endl;
            for ( uint32_t j{0U}; j<XY.Cols(); j++ )
                cerr << XY(i,j) <<" ";
            cerr<<endl;
          }
        csmp_error.notice( WARNING, "soparametricQuadraticTriangle::JacobianInverse:",
                          "the value of the Jacobian is negative; check node-numbering.");

        return fabs(detJ);
     }

    return detJ;
   
 } // end JacobianInverse




/** Computes the normal to the triangle surface assuming that the triangle
is perfectly planar.


@section implementation Implementation

Method assumes that the triangle is planar. For the case that the
triangle is warped, normals must be calculated for each integration
point.
*/
void IsoparametricQuadraticTriangle::UnitNormal( std::vector<double>& vc ) const
 {
    // normal only exists in 3D
    vc.resize(3U);

    if ( dim == 2U ) {
         vc[0] = vc[1] = static_cast<double>(0.0);
         vc[2] = static_cast<double>(1.0);
         return;
      }

    double X12 = XY(1,0) - XY(0,0), // X
              X31 = XY(0,0) - XY(2,0),
              Y12 = XY(1,1) - XY(0,1), // Y
              Y31 = XY(0,1) - XY(2,1),
              Z12 = XY(1,2) - XY(0,2), // Z
              Z31 = XY(0,2) - XY(2,2);

  // normal to triangle (but not unit normal!)
    vc[0]  = -Y12*Z31 + Z12*Y31;
    vc[1]  = -Z12*X31 + X12*Z31;
    vc[2]  = -X12*Y31 + Y12*X31;
    // normalization to unit length
    double length = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
    vc[0] /= length;
    vc[1] /= length;
    vc[2] /= length;
 }




/** The integration point location is transformed and output into the
second argument in global coordinates.

@warning To work correctly the matrix XYZ must be uptodate.
*/
void  IsoparametricQuadraticTriangle::IntegrationPoint( uint32_t ip,
                                                        vector<double>& xyz ) const
 {
    assert( &xyz != &NRST );
    assert( ip < gpe );
    xyz.resize(dim);

    // local interpolation function values
    Nrs( rr[ip], ss[ip], NRST );

    // 2D (unrolled loop)
    if ( dim == 2U ) {
       xyz[0] = XY(0,0) * NRST[0];
       xyz[1] = XY(0,1) * NRST[0];
       xyz[0] += XY(1,0) * NRST[1];
       xyz[1] += XY(1,1) * NRST[1];
       xyz[0] += XY(2,0) * NRST[2];
       xyz[1] += XY(2,1) * NRST[2];
       xyz[0] += XY(3,0) * NRST[3];
       xyz[1] += XY(3,1) * NRST[3];
       xyz[0] += XY(4,0) * NRST[4];
       xyz[1] += XY(4,1) * NRST[4];
       xyz[0] += XY(5,0) * NRST[5];
       xyz[1] += XY(5,1) * NRST[5];
         return;
      }

    // 3D case
    xyz[0] = XY(0,0) * NRST[0];
    xyz[1] = XY(0,1) * NRST[0];
    xyz[2] = XY(0,2) * NRST[0];
    xyz[0] += XY(1,0) * NRST[1];
    xyz[1] += XY(1,1) * NRST[1];
    xyz[2] += XY(1,2) * NRST[1];
    xyz[0] += XY(2,0) * NRST[2];
    xyz[1] += XY(2,1) * NRST[2];
    xyz[2] += XY(2,2) * NRST[2];
    xyz[0] += XY(3,0) * NRST[3];
    xyz[1] += XY(3,1) * NRST[3];
    xyz[2] += XY(3,2) * NRST[3];
    xyz[0] += XY(4,0) * NRST[4];
    xyz[1] += XY(4,1) * NRST[4];
    xyz[2] += XY(4,2) * NRST[4];
    xyz[0] += XY(5,0) * NRST[5];
    xyz[1] += XY(5,1) * NRST[5];
    xyz[2] += XY(5,2) * NRST[5];

 } // end IntegrationPoint




// AP 2006
void IsoparametricQuadraticTriangle::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords(0,0) = NXY(0,0);
    matCoords(0,1) = NXY(0,1);
    matCoords(1,0) = NXY(1,0);
    matCoords(1,1) = NXY(1,1);
    matCoords(2,0) = NXY(2,0);
    matCoords(2,1) = NXY(2,1);
    matCoords(3,0) = NXY(3,0);
    matCoords(3,1) = NXY(3,1);
    matCoords(4,0) = NXY(4,0);
    matCoords(4,1) = NXY(4,1);
    matCoords(5,0) = NXY(5,0);
    matCoords(5,1) = NXY(5,1);
}


} // end namespace csmp
