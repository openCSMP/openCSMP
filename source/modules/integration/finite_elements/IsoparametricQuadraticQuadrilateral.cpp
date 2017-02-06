#include "IsoparametricQuadraticQuadrilateral.h"
#include <set>
#include <fstream>
#include <cstring>
#include "Exception.h"
#include <climits>

using namespace std;

namespace csmp {

/**

Default constructor initializes quadrilateral to use local coordinates
and a four-point integration scheme where the Gauss points are located
inside quad at convergence points. The weights of the integration points
are set to 1 such that the integrated contributions add up to 1.  .
*/

IsoparametricQuadraticQuadrilateral::IsoparametricQuadraticQuadrilateral( size_t dimensions )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9, true, true, 2U ),
   use2Dto3Djacobi( !(dimensions==2) ),
   DN(dimensions,9),
   BEE(dimensions,9),
   NXY(9,dimensions),
   JMAT(dimensions,dimensions),
   EFG(dimensions),
   RS(2)
 {
   dim = dimensions;
   itp = 2;
   npf = 3;
   npe = 9;
   fpe = 4;
   spe = 4;
   epe = 4;
   nne = 4;
   cne = 4;
   gpe = 4;
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

   // full 4-point integration basis
   double64 sqRootOneThird=sqrt(1.0/3.0);
   rr[0] = sqRootOneThird,   rr[1] = -sqRootOneThird,   rr[2] = -sqRootOneThird,   rr[3] = sqRootOneThird;
   ss[0] = sqRootOneThird,   ss[1] =  sqRootOneThird,   ss[2] = -sqRootOneThird,   ss[3] =-sqRootOneThird;

   W[0]  = W[1]  = W[2]  =  W[3]  = 1.0;

   // local node coordinates are defined as r==ksi, s==nu, t==mu.
   NXY(0,0) =-1.0;NXY(0,1) =-1.0;
   NXY(1,0) = 1.0;NXY(1,1) =-1.0;
   NXY(2,0) = 1.0;NXY(2,1) = 1.0;
   NXY(3,0) =-1.0;NXY(3,1) = 1.0;

   NXY(4,0) = 0.0;NXY(4,1) =-1.0;
   NXY(5,0) = 1.0;NXY(5,1) = 0.0;
   NXY(6,0) = 0.0;NXY(6,1) = 1.0;
   NXY(7,0) =-1.0;NXY(7,1) = 0.0;

   NXY(8,0) = 0.0;NXY(8,1) = 0.0;

   UsesLocalCoordinates(true);
   Isoparametric(true);
   SurfaceElement();
   ElementType(ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9);
}



IsoparametricQuadraticQuadrilateral::~IsoparametricQuadraticQuadrilateral()
 {
 }



/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

The interpolation function values are returned into the third argument.

*/
void
IsoparametricQuadraticQuadrilateral::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )
 {
    assert( ip < gpe );

    // local interpolation function values
    Nrs( rr[ip], ss[ip], N );

 } // end N_AtIntegrationPoint


/**

Gives the local node numbers (0...5) of the nodes which make up the
requested face of the element (possible numbers 0, 1, 2).
For this element, the faces are numbered such that face 0 lies opposite of
node 0, face 1 node 1 etc.

@section arguments Input Arguments

The number of the desired face (0, 1, or 2).

The local node numbers are returned into the integer vector 'fnids'.

@section application Application

When boundary conditions shall be applied it is necessary to determine
the properties associated with the nodes of that element face.
*/
void
IsoparametricQuadraticQuadrilateral::NodesOfFace( size_t face_id,
                                    std::vector<size_t>& fnids ) const
 {
    fnids.resize(3);

    if ( face_id == 0 )
      {
         fnids[0] = 0;
         fnids[1] = 4;
         fnids[2] = 1;

      }
    else if ( face_id == 1 )
      {
         fnids[0] = 1;
         fnids[1] = 5;
         fnids[2] = 2;
      }
    else if ( face_id == 2 )
      {
         fnids[0] = 2;
         fnids[1] = 6;
         fnids[2] = 3;
      }
    else if ( face_id == 3 )
      {
         fnids[0] = 3;
         fnids[1] = 7;
         fnids[2] = 0;
      }
    else
    std::cout <<"\nIsoparametricQuadraticQuadrilateral::NodesOfFace: Erratic input face ID: "<< face_id << std::endl;
 }




void
IsoparametricQuadraticQuadrilateral::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    snids.resize(3);
    if ( segm_id == 0 ) {
         snids[0] = 0;
         snids[1] = 1;
         snids[2] = 4;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 1;
         snids[1] = 2;
         snids[2] = 5;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 2;
         snids[1] = 3;
         snids[2] = 6;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 3;
         snids[1] = 0;
         snids[2] = 7;
      }
    else
    std::cout <<"\nIsoparametricQuadraticQuadrilateral::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment




/** Returns local node ids of the nodes which sit on the corner of
the quadratic triagular element.

Returns a Meschpp integer vector with the 3 local corner node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).

*/
void
IsoparametricQuadraticQuadrilateral::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(4);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
 }



/** Returns local node ids of the nodes which sit on the midsides of
the quadratic triagular element.

This is an integer vector with the 3 local midside-node ID numbers
for the element (These integers are germane to all quadratic triangular
elements in the mesh).

*/
void  IsoparametricQuadraticQuadrilateral::MidSideNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(4);
    ids[0] = 4;
    ids[1] = 5;
    ids[2] = 6;
    ids[3] = 7;
 }

void  IsoparametricQuadraticQuadrilateral::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 4;
    ids[2] = 1;
    ids[3] = 5;
    ids[4] = 2;
    ids[5] = 6;
    ids[6] = 3;
    ids[7] = 7;
    ids[8] = 8;
 }

CSMP_FEM_TYPE IsoparametricQuadraticQuadrilateral::ElementTypeOfFace( size_t ) const
 {
    return ISOPARAMETRIC_QUADRATIC_BAR;
 }


double64 IsoparametricQuadraticQuadrilateral::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }


void  IsoparametricQuadraticQuadrilateral::JacobianAtIntegrationPoint( size_t gauss_point )
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

void  IsoparametricQuadraticQuadrilateral::N_AtBaryCenter( std::vector<double64>& N )
 {
    N.resize(npe);
      Nrs(0.0,0.0,N);
 }


void IsoparametricQuadraticQuadrilateral::Dimensions( size_t dimensions )
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



/** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64
IsoparametricQuadraticQuadrilateral::AspectRatio()
{
   vector<double64> vec(spe);

   EdgeLengths( vec );

   // order segment
   set<double64> segms;

   for ( size_t i=0U; i<spe; i++ ) segms.insert( vec[i] );

   double64 segm1 = (*segms.begin()),
            segm2 = (*segms.rbegin());

   return segm2 / segm1;
}


/**

This method attempts to calculate the inner radius of the curved-sided
triangle by talking 1/2 of the perimeter and dividing this measure by
the area of the quad.

The procedure is empirically based giving a good match if quads
are relatively even-sided and have straight edges.

@section arguments Input Arguments

The parent element is queried for its node coordinates.
*/
double64  IsoparametricQuadraticQuadrilateral::InnerRadius()
{

   vector<double64> segms(npe);
   double64         sum(0.0), vol;

   EdgeLengths( segms );
   for ( size_t i=0; i<segms.size(); i++ ) sum += segms[i];
   sum /= 2.0;
   vol  = Volume();
   vol /= sum;

   cout<<" IsoparametricQuadraticQuadrilateral::InnerRadius():  ***WARNING: Approximate radius - not valid for HAR elements"<<endl;

   return vol;
}


/**

Segment length is approximated as the sum of the length of two linear
segments making up each face of the quadratic triangle. The first segment
connects corner node 1 with midpoint node 3 and midpoint node 3 with
corner node 2 and so forth.

The lengths of the segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as the sum of the straight line
segments which make up the element face.

@section application Application

For example, when stresses are to be applied at the element side, the
(area=length at unit thickness) must be taken into account.
*/
void
IsoparametricQuadraticQuadrilateral::EdgeLengths( std::vector<double64>& len )
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
    sum     = (XY(3,0)-XY(2,0)) * (XY(3,0)-XY(2,0));
    sum    += (XY(3,1)-XY(2,1)) * (XY(3,1)-XY(2,1));
    sum    += (XY(3,2)-XY(2,2)) * (XY(3,2)-XY(2,2));
    len[2]  = sqrt(sum);

    // segment 4
    sum     = (XY(0,0)-XY(3,0)) * (XY(0,0)-XY(3,0));
    sum    += (XY(0,1)-XY(3,1)) * (XY(0,1)-XY(3,1));
    sum    += (XY(0,2)-XY(3,2)) * (XY(0,2)-XY(3,2));
    len[3]  = sqrt(sum);
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
    sum     = (XY(3,0)-XY(2,0)) * (XY(3,0)-XY(2,0));
    sum    += (XY(3,1)-XY(2,1)) * (XY(3,1)-XY(2,1));
    len[2]  = sqrt(sum);

    // segment 4
    sum     = (XY(0,0)-XY(3,0)) * (XY(0,0)-XY(3,0));
    sum    += (XY(0,1)-XY(3,1)) * (XY(0,1)-XY(3,1));
    len[3]  = sqrt(sum);
    }
 } // end EdgeLengths




/**

Computes Integral N dA = 1/2 sum_1...n Wi Ji Ni. This formulation also
gives a meaningful area if the element boundaries are curved. For a
description of numerical integration  of isoparametric quadratic
elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The area (m2) of the finite element.
*/
double64  IsoparametricQuadraticQuadrilateral::Volume()
{
    double64  area;
    size_t  i;

    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
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




/** Returns the value of the element interpolation functions at the point 'rst'
in local coordinates.

@section arguments Input Arguments

The floating point coordinates 'r', 's' and 't' of the point at which the
interpolation shall be carried out.

The fourth method argument is the vector into which the values of the
n=nodes interpolation functions at the point 'rst' will be returned.

@section implementation Implementation

See in source code header file.

@section application Application

Method is used to compute property values at the integration points of
the element.
*/
void
IsoparametricQuadraticQuadrilateral::Nrs(
                    double64 r,
                    double64 s,
                    std::vector<double64>& N ) const
{
   N.resize(npe);
   const double64 rPlus=1.0+r;
   const double64 sPlus=1.0+s;
   const double64 rMinus=1.0-r;
   const double64 sMinus=1.0-s;

   N[0] = 0.25*r*s*rMinus*sMinus;
   N[1] =-0.25*r*s*rPlus *sMinus;
   N[2] = 0.25*r*s*rPlus *sPlus;
   N[3] =-0.25*r*s*rMinus*sPlus;

   N[4] =-0.5*s*sMinus*(1.-r*r);
   N[5] = 0.5*r*rPlus *(1.-s*s);
   N[6] = 0.5*s*sPlus *(1.-r*r);
   N[7] =-0.5*r*rMinus*(1.-s*s);
   N[8] = (1.-r*r)*(1.-s*s);

}



/**

This method and the complementary method dNs() compute the shape function
derivates with respect to the local coordinate axis 'r', and 's'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

The shape function derivatives are returned into the second method
argument which is a floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.
*/
void
IsoparametricQuadraticQuadrilateral::dNr (
                                            double64 r,
                                            double64 s,
                                            std::vector<double64>& DNR ) const
{
   DNR.resize(npe);
   DNR[0] = 0.25*s*(1.0 - r)*(1.0 - s) - 0.25*r*s*(1.0 - s);
   DNR[1] = -0.25*s*(1.0 + r)*(1.0 - s) - 0.25*r*s*(1.0 - s);
   DNR[2] = 0.25*s*(1.0 + r)*(1.0 + s) + 0.25*r*s*(1.0 + s);
   DNR[3] = -0.25*s*(1.0 - r)*(1.0 + s) + 0.25*r*s*(1.0 + s);
   DNR[4] = 1.0*r*s*(1.0 - s);
   DNR[5] = 0.5*(1.0 + r)*(1. - s*s) + 0.5*r*(1. - s*s);
   DNR[6] = -1.0*r*s*(1.0 + s);
   DNR[7] = -0.5*(1.0 - r)*(1. - s*s) + 0.5*r*(1. - s*s);
   DNR[8] = -2.*r*(1. - s*s);
}

void IsoparametricQuadraticQuadrilateral::dNs(
                                                        double64 r,
                                                        double64 s,
                                                        std::vector<double64>& DNS ) const
{
   DNS.resize(npe);
   DNS[0] = .25*r*(1-r)*(1-s)-.25*r*(1-r)*s;
   DNS[1] =-.25*r*(1+r)*(1-s)+.25*r*(1+r)*s;
   DNS[2] = .25*r*(1+r)*(1+s)+.25*r*(1+r)*s;
   DNS[3] =-.25*r*(1-r)*(1+s)-.25*r*(1-r)*s;

   DNS[4] = -.5*(1-s)*(1-r*r)+.5*s*(1-r*r);
   DNS[5] = -1.0*r*(1+r)*s;
   DNS[6] =  0.5*(1+s)*(1-r*r)+.5*s*(1-r*r);
   DNS[7] =  1.0*r*(1-r)*s;

   DNS[8] = -2.0*s*(1-r*r);

}





/**

Compute derivatives of interpolation functions at corresponding nodes with
respect to the global coordinate system. If a transformation from 2D local
to 3D global space is applied, the function Jacobi is called.

@section arguments Input Arguments

The element is used to obtain the global interpolation of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricQuadraticQuadrilateral::dN( DenseMatrix<DM_MIN>& DN9 )
  {
     DN9.Resize(dim,npe);
     double64 detJ(std::numeric_limits<double>::quiet_NaN());

     // Jacobian transformation to global coordinate system
     if ( use2Dto3Djacobi )
       {
          for ( size_t i=0; i<npe; i++ ) {
                RS[0] = NXY(i,0);
                RS[1] = NXY(i,1);

                // generate the 3 x 4 shape function derivative matrix
                detJ = Jacobi( RS, EFG, JMAT );
                double64 det_inverse = 1.0 / ( detJ * detJ );
                DN9(0,i)  = det_inverse * JMAT(0,0) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN9(0,i) += det_inverse * JMAT(1,0) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN9(1,i)  = det_inverse * JMAT(0,1) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN9(1,i) += det_inverse * JMAT(1,1) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
                DN9(2,i)  = det_inverse * JMAT(0,2) * ( EFG[2] * DNR[i] - EFG[1] * DNS[i] );
                DN9(2,i) += det_inverse * JMAT(1,2) * ( EFG[0] * DNS[i] - EFG[1] * DNR[i] );
             }
           return;
        }

     // if the quadralateral is 2D
     DenseMatrix<DM_MIN> TEMP(dim,1);
     for ( size_t i=0; i<npe; i++ )
       {
           dNr( NXY(i,0), NXY(i,1), DNR );
           dNs( NXY(i,0), NXY(i,1), DNS );
           Jacobian( DNR, DNS );
           JacobianInverse();

           TEMP(0,0)=DNR[i];TEMP(1,0)=DNS[i];

           JINV *= TEMP;

           DN9(0,i) = JINV(0,0);
           DN9(1,i) = JINV(1,0);

            JINV.Resize(dim,dim);
       }
  } // end dN





/** Projection function from rs->xy(z)
*/
void
IsoparametricQuadraticQuadrilateral::ParametricToPhysical(std::vector<double64> &rst, std::vector<double64>& xyz)
{
Nrs(rst[0],rst[1], NRST );

for(size_t i=0; i<dim; i++) xyz[i]=0.0;

    if(dim==2)
    {
        for(size_t i=0; i<npe; i++)
        {
        xyz[0]+=XY(i,0)*NRST[i];
        xyz[1]+=XY(i,1)*NRST[i];
        }
    }
    else if(dim==3)
    {
    for(size_t i=0; i<npe; i++)
        {
        xyz[0]+=XY(i,0)*NRST[i];
        xyz[1]+=XY(i,1)*NRST[i];
        xyz[2]+=XY(i,2)*NRST[i];
        }
    }
}

/** Projection function from xyz->rst

In first part over regular grid looks for closest point to the given one;
Then the point is used as initial guess for a Newton-Pahson iterration with
constant Jacobian matrix
*/
void
IsoparametricQuadraticQuadrilateral::PhysicalToParametric(
                                       std::vector<double64>& rSt,
                                       const std::vector<double64>& xyz
                                        )
{
    vector<double64> outxyz(dim);
    vector<double64> rstHatK(2U);

    std::vector<double64> distanceFromGivenPointLinf(3U,0.0);
    double64 distanceFromGivenPointL2;

    // Find largest and smallest segments in order to define precision
    vector<double64> vec(spe);
    EdgeLengths( vec );
    double64 seg_max(vec[0]), seg_min(vec[0]);
    for ( size_t i=1; i<spe; i++ )
    {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
    }

    const double64 geometricTolerance = 0.005*seg_min;

    // First guess as BaryCenter
    rstHatK[0] = 0.0;
    rstHatK[1] = 0.0;

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

        vector<double64> rstHatK_PlusOne(2U);
        double64 minDistanceFromGivenPoint;

        const double64 constantMu               = 1.0;
        const size_t numberOfFirstIterrations   = 5;
        const size_t maxNumberOfIterrations     = 20;
        const size_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double64 incrementR = 2.0/(numberOfFirstIterrations-1);
        const double64 incrementS = 2.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = -1.0;
        for(size_t i=0;i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = -1.0;
            for(size_t j=0;j<numberOfFirstIterrations;j++)
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
                    for(size_t i=0; i<2U; i++)
                        rstHatK[i] = rstHatK_PlusOne[i];

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

        size_t iteration = 1;

        /// Newton-Raphson iterations
        while ( ( distanceFromGivenPointL2 > geometricTolerance ) && ( iteration < maxNumberOfIterrations ) )
        //while ( ( ( distanceFromGivenPointLinf[0] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[1] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[2] > geometricTolerance )    )
        //      &&  ( iteration < maxNumberOfIterrations )                )
        {
            // Out of range check
            if( (rstHatK[0]<-1.) || (rstHatK[0]>1.) ||
                (rstHatK[1]<-1.) || (rstHatK[1]>1.)  )
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
                //const double64 detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            if ( dim == 3 )
            {
                rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
                rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
                rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

                for(size_t i=0; i<2; i++)
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

                for(size_t i=0; i<2; i++)
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
            cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            if( dim ==3 )
            {
                cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Real Point:\t"
                    <<"x = "<<xyz[0]<<" ;\t"
                    <<"y = "<<xyz[1]<<" ;\t"
                    <<"z = "<<xyz[2]<<"\n";
                cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Found Point:\t"
                   <<"x = "<<outxyz[0]<<" ;\t"
                   <<"y = "<<outxyz[1]<<" ;\t"
                   <<"z = "<<outxyz[2]<<"\n";
                cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Distance:\t"
                   <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
                   <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
                   <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            }else{
                cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Real Point:\t"
                    <<"x = "<<xyz[0]<<" ;\t"
                    <<"y = "<<xyz[1]<<"\n";
                cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Found Point:\t"
                    <<"x = "<<outxyz[0]<<" ;\t"
                    <<"y = "<<outxyz[1]<<"\n";
                cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Distance:\t"
                    <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
                    <<"y = "<<distanceFromGivenPointLinf[1]<<endl;
            }

            cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<"\n";

            std::vector<double64> N(npe,0.0);
            Nrs(rstHatK[0], rstHatK[1], N );
            cout<<" IsoparametricQuadraticQuadrilateral::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<" ;\t"
                <<"N[5] = "<<N[5]<<" ;\t"
                <<"N[6] = "<<N[6]<<" ;\t"
                <<"N[7] = "<<N[7]<<" ;\t"
                <<"N[8] = "<<N[8]<<"\n";

            csmp::Exception( WARNING, "IsoparametricQuadraticQuadrilateral::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for(size_t i=0; i<2; i++)
        rSt[i] = rstHatK[i];

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
IsoparametricQuadraticQuadrilateral::dN_At( DenseMatrix<DM_MIN>& DN2,
                                                            const vector<double64>& xyz  )
  {

    vector<double64> rst(dim);
    // As for the time being par->phys mapping is fully defined for 2D only, excude 3D case
    if ( dim == 3 )
         {
        csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticQuadrilateral::PhysicalToParametric",
                          "Function can only be used if dim = 2");
            return -1.0;
        }

    PhysicalToParametric(rst, xyz);

    // here find derivatives
    dNr( rst[0], rst[1], DNR );
    dNs( rst[0], rst[1], DNS );

    Jacobian( DNR, DNS );
    double64 detJ = JacobianInverse();

    DN2.Resize(dim,dim);
    DN2  = JINV;

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    cout<<" IsoparametricQuadraticQuadrilateral3D::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
    rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
    cout<<" IsoparametricQuadraticQuadrilateral3D::dN  Matrix DN2: "<<endl;
    DN2.Out(cout);
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
IsoparametricQuadraticQuadrilateral::Jacobi(	const vector<double64>& rs,
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


double64 IsoparametricQuadraticQuadrilateral::Jacobi( const vector<double64>& rs )
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

The second method argument, the vector FN, will hold the interpolation function
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
IsoparametricQuadraticQuadrilateral::N( vector<double64>& N, const vector<double64>& xyz )
 {
    vector<double64> rst(parametricDimensions);

    PhysicalToParametric(rst, xyz);

    Nrs(rst[0],rst[1], N);

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
IsoparametricQuadraticQuadrilateral::dN_AtIntegrationPoint(
                                            DenseMatrix<
                                            DM_MIN>& B,
                                            size_t gauss_point )
 {
    double64 det;

    assert( gauss_point < gpe );

    if ( use2Dto3Djacobi )
      {
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

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3];

    B(0,4) = DNR[4], B(1,4) = DNS[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5];
    B(0,6) = DNR[6], B(1,6) = DNS[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7];
    B(0,8) = DNR[8], B(1,8) = DNS[8];

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
IsoparametricQuadraticQuadrilateral::dN_AtNode(
                                DenseMatrix<
                                DM_MIN>& B,
                                size_t nd
                                )
 {
    double64  det;

    assert( nd < npe );

    if ( use2Dto3Djacobi )
      {
        static vector<double64> EFG(dim);
        RS[0] = NXY(nd,0);
        RS[1] = NXY(nd,1);

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
    B(0,3) = DNR[3], B(1,3) = DNS[3];

    B(0,4) = DNR[4], B(1,4) = DNS[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5];
    B(0,6) = DNR[6], B(1,6) = DNS[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7];
    B(0,8) = DNR[8], B(1,8) = DNS[8];

    B = JINV * B;

    return det;
 }

double64
IsoparametricQuadraticQuadrilateral::dN_AtBarycenter(
                                    DenseMatrix<
                                    DM_MIN>& B
                                    )
 {
    double64   det;

    if ( use2Dto3Djacobi )
      {
        RS[0] = 0.0;
        RS[1] = 0.0;

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

    // 2D case
    // ----------------------------------------------------
    dNr( 0.0, 0.0, DNR );
    dNs( 0.0, 0.0, DNS );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS );
    det = JacobianInverse();

    // Compute global DN by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3];

    B(0,4) = DNR[4], B(1,4) = DNS[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5];
    B(0,6) = DNR[6], B(1,6) = DNS[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7];
    B(0,8) = DNR[8], B(1,8) = DNS[8];

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

The resulting local node id's are returned into the third method argument.

@section implementation Implementation

While only the element knows which nodes are located at the mode boundary,
the FiniteElement knows in which order these appear.

@section application Application

To assign Neumann boundary conditions with a PDE operator for surface
integrals.

@todo (3) Check index of bnodes
*/
void
IsoparametricQuadraticQuadrilateral::ConsecutiveNodesAtBoundary(
                                            const vector<size_t>& bnodes,
                                            vector<size_t>& fnids )
 {
     fnids.resize(3);

     // in 3D, all nodes may be at model boundary
     if ( use2Dto3Djacobi && bnodes.size() != 2 )
       throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticQuadrilateral::ConsecutiveNodesAtBoundary",
                              "Cannot resolve node sequence for element boundary",
                              "Probably because element lies at two boundaries simultaneously" );

     if ( bnodes.size() > 2 )
       throw csmp::Exception( ERROR, "IsoparametricQuadraticQuadrilateral::ConsecutiveNodesAtBoundary",
                              "Cannot resolve node sequence for element boundary",
                              "Probably because element lies at two boundaries simultaneously" );

     // the midside nodes determine the boundaries so they are sought after
     if ( bnodes[0] == 2 || bnodes[1] == 2 || bnodes[2] == 2 ) {
           fnids[0] = 0, fnids[1] = 4, fnids[2] = 1;
           return;
       }
     if ( bnodes[0] == 3 || bnodes[1] == 3 || bnodes[2] == 3 ) {
           fnids[0] = 1, fnids[1] = 5, fnids[2] = 2;
           return;
       }
     if ( bnodes[0] == 4 || bnodes[1] == 4 || bnodes[2] == 4 ) {
           fnids[0] = 2, fnids[1] = 6, fnids[2] = 3;
           return;
       }

      if ( bnodes[0] == 5 || bnodes[1] == 5 || bnodes[2] == 5 ) {
          fnids[0] = 3, fnids[1] = 7, fnids[2] = 0;
          return;
      }

 } // end ConsecutiveNodesAtBoundary


void
IsoparametricQuadraticQuadrilateral::OutputNodeDataToVTK( const char* file_name,
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
           cout <<"\nIsoparametricQuadraticQuadrilateral3D::OutputNodeDataToVTK ";
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
     for ( size_t i=0; i<npe; i++ )
       {

        if(dim==2)
        {
            for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
            ofs<<0.0;
            ofs << endl;
        }
        else
        if(dim==3)
        {
            for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
            ofs << endl;
        }
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point)) -corect
     // Cells
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< 10 << endl;
     ofs << 9 <<" 0 4 1 5 2 6 3 7 0" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     // Unfortunately the data can only be output as nodal variables
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );

     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" float"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x9
           for ( size_t i=0; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 10
          if(dim==2)
          {
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
               for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs<<0.0<<" ";
               ofs << endl;
            }
          }
          else
          {
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
               for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
          }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricQuadraticQuadrilateral3D::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

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
IsoparametricQuadraticQuadrilateral::ExtrapolateIntegrationPointVariableToNodes(
                                                                              size_t nvars,
                                                        const vector<double64>& IVAR,
                                                  vector<double64>& NVAR	)
const
{
 assert( gpe == 4);
 assert( IVAR.size() >= (gpe*nvars) );
   NVAR.resize( npe * nvars );

  // Define nodal values as bi-linear variation of the integration points values
  // See Zienkevitch, pp. 351, for example

  const size_t corner_nodes(4);

  DenseMatrix<DM_MIN> MATRIX_A(gpe,gpe), TEMP_IP(gpe,1), TEMP_N(corner_nodes,1);
  // coefficients of the matrix A
  double64 a=1+sqrt(3.0)/2.0, b=-0.5, c=1-sqrt(3.0)/2.0;

     MATRIX_A(0,0)=MATRIX_A(1,1)=MATRIX_A(2,2)=MATRIX_A(3,3) =a;
     MATRIX_A(0,1)=MATRIX_A(0,3)=MATRIX_A(1,0)=MATRIX_A(1,2)
     =MATRIX_A(2,1)=MATRIX_A(2,3)=MATRIX_A(3,0)=MATRIX_A(3,2)=b;
     MATRIX_A(0,2)=MATRIX_A(1,3)=MATRIX_A(2,0)=MATRIX_A(3,1) =c;

   for ( size_t i=0; i<nvars; i++ )
    {
         for ( size_t k=0; k<gpe; k++ ) TEMP_IP(k,0)=IVAR[k*nvars +i];
         TEMP_N = MATRIX_A * TEMP_IP;
         // form corner nodes part of vector
         for ( size_t j=0; j<corner_nodes; j++ ) NVAR[j*nvars+i]=TEMP_N(j,0);

         // average midside nodes
         NVAR[4*nvars + i]=0.5*(TEMP_N(0,0)+TEMP_N(1,0));
         NVAR[5*nvars + i]=0.5*(TEMP_N(1,0)+TEMP_N(2,0));
         NVAR[6*nvars + i]=0.5*(TEMP_N(2,0)+TEMP_N(3,0));
         NVAR[7*nvars + i]=0.5*(TEMP_N(3,0)+TEMP_N(1,0));
             NVAR[8*nvars + i]=0.25*(NVAR[7*nvars + i]+NVAR[6*nvars + i]+NVAR[5*nvars + i]+NVAR[4*nvars + i]);
    }

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)







/// overwrites the base class method in order to get 2D to 3D mapping functionality
double64
IsoparametricQuadraticQuadrilateral::JacobianInverse()
 {
    if ( use2Dto3Djacobi )
      {
         // Compute Jacobi J' := "determinant" of the 3x2 Jacobian
         // ------------------------------------------------------
         // using equation J' = ( E * g - F^2 )^0.5 see FEM-development-in-CSP.doc equation (15)
         fill( EFG.begin(), EFG.end(), 0.0 );

         for ( size_t i=0; i<dim; i++ ) {
             EFG[0] += JAC(0,i) * JAC(0,i);
             EFG[1] += JAC(0,i) * JAC(1,i);
             EFG[2] += JAC(1,i) * JAC(1,i);
           }

         // compute J'
         return sqrt( EFG[0] * EFG[2] - EFG[1] * EFG[1] );
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
         std::cout <<"\nIsoparametricQuadraticQuadrilateral::JacobianInverse: Erroneous determinant of Jacobian matrix: ";
         std::cout << detJ << std::endl;
         throw std::range_error("IsoparametricQuadraticQuadrilateral::JacobianInverse");
      }

    return detJ;
 }


/**

Computes the normal to the triangle surface assuming that the quad
is perfectly planar, based on 3 points definition of plane.

@section implementation Implementation

Method assumes that the quad is planar. For the case that the
quad is warped, normals must be calculated for each integration
point.
*/
void
IsoparametricQuadraticQuadrilateral::UnitNormal( std::vector<double64>& vc ) const
 {
    // normal only exists in 3D
    vc.resize(3);

    if ( dim == 2 ) {
         vc[0] = vc[1] = static_cast<double64>(0.0);
         vc[2] = static_cast<double64>(1.0);
         return;
      }

    double64 X12 = XY(1,0) - XY(0,0), // X
              X31 = XY(0,0) - XY(2,0),
              Y12 = XY(1,1) - XY(0,1), // Y
              Y31 = XY(0,1) - XY(2,1),
              Z12 = XY(1,2) - XY(0,2), // Z
              Z31 = XY(0,2) - XY(2,2);

  // normal to quad
    vc[0]  = -Y12*Z31 + Z12*Y31,
    vc[1]  = -Z12*X31 + X12*Z31,
    vc[2]  = -X12*Y31 + Y12*X31;
    // normalization to unit length
    double64 length = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
    vc[0] /= length;
    vc[1] /= length;
    vc[2] /= length;
 }





/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricQuadraticQuadrilateral::IntegrationPoint( size_t ip,
                                                             vector<double64>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(NXY.Cols());
    xyz[0]=xyz[1]=0;

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
IsoparametricQuadraticQuadrilateral::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
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
  matCoords(6,0) = NXY(6,0);
  matCoords(6,1) = NXY(6,1);
  matCoords(7,0) = NXY(7,0);
  matCoords(7,1) = NXY(7,1);
}
}// namespace csp

