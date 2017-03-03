#include "IsoparametricQuadraticLineElement.h"
#include "Exception.h"
#include "MJL_Edge.h"

using namespace std;

namespace csmp {


IsoparametricQuadraticLineElement::IsoparametricQuadraticLineElement( size_t dimensions )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_BAR, true, true, 2U ),
    current_detJ(std::numeric_limits<double64>::quiet_NaN())
 {
    dim = dimensions;
    itp = 2;
    npf = 3;
    npe = 3;
    fpe = 2;
    spe = 1;
    epe = 2;
    nne = 2;
    cne = 0;
    gpe = 2;
    XY.Resize(npe, dim);
    M.Resize(npe,npe);
    JAC.Resize(dim,dim);
    JINV.Resize(dim,dim);

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);

    UsesLocalCoordinates( true );
    Isoparametric( true );
    LineElement();
    ElementType(ISOPARAMETRIC_QUADRATIC_BAR);

    IP.resize( gpe );
    W.resize( gpe );

    // Cook et al. p. 171, table 6.4-1
    IP[0]   = -1. / sqrt(3.), IP[1]     = 1. / sqrt(3.);
    W[0]    =  1.,            W[1] = 1.; // since length = 2.0 (-1,1)

    // local node coordinates (node order 1,2,3)
    NX[0] = -1.;
    NX[1] =  1.;
    NX[2] =  0.;
 }


IsoparametricQuadraticLineElement::~IsoparametricQuadraticLineElement() 
 {
 }




/// Cheung et al. p. 26, Cook et al., p. 164
void IsoparametricQuadraticLineElement::Nr( double64 r, std::vector<double64>& N ) const
{
   N.resize(npe);
   double64 r2 = r * r;
   N[0] = 0.5 * (-r + r2);
   N[1] = 0.5 * ( r + r2);
   N[2] = 1. - r2;
}


/// Cheung et al. p. 26
void IsoparametricQuadraticLineElement::dNr( double64 r, std::vector<double64>& dnr ) const
{
   double64 rr = 2. * r;
   dnr[0] = 0.5 * (-1. + rr);
   dnr[1] = 0.5 * ( 1. + rr);
   dnr[2] = -rr;
}


double64 IsoparametricQuadraticLineElement::JacobianFor( const std::vector<double64>& DNR, size_t coord ) const
 {
    assert( coord < dim );
    return DNR[0] * XY(0,coord) + DNR[1] * XY(1,coord) + DNR[2] * XY(2,coord);
 }


// these next 3 auxiliary methods are all used by the Volume() function
double64 IsoparametricQuadraticLineElement::Jacobian1D( const std::vector<double64>& DNR ) const
 {
    return DNR[0] * XY(0,0) + DNR[1] * XY(1,0) + DNR[2] * XY(2,0);
 }


double64 IsoparametricQuadraticLineElement::Jacobian2D( const std::vector<double64>& DNR ) const
 {
    double64 jac_x = DNR[0] * XY(0,0) + DNR[1] * XY(1,0) + DNR[2] * XY(2,0);
    double64 jac_y = DNR[0] * XY(0,1) + DNR[1] * XY(1,1) + DNR[2] * XY(2,1);
    
    // compute a length-like expression for the transformation
    return std::sqrt( jac_x*jac_x + jac_y*jac_y );
    
 } // end 


double64 IsoparametricQuadraticLineElement::Jacobian3D( const std::vector<double64>& DNR ) const
 {
    double64 jac_x = DNR[0] * XY(0,0) + DNR[1] * XY(1,0) + DNR[2] * XY(2,0);
    double64 jac_y = DNR[0] * XY(0,1) + DNR[1] * XY(1,1) + DNR[2] * XY(2,1);
    double64 jac_z = DNR[0] * XY(0,2) + DNR[1] * XY(1,2) + DNR[2] * XY(2,2);
    
    return std::sqrt( jac_x*jac_x + jac_y*jac_y + jac_z*jac_z );
    
 } // end 


// Next 2 methods belong together
/// return determinant J for values of previous function
double64 IsoparametricQuadraticLineElement::JacobianInverse() { return current_detJ; }


void IsoparametricQuadraticLineElement::JacobianAtIntegrationPoint( size_t ipoint )
 {
    if ( ipoint > 1 ) {
         std::cerr <<"\nIsoparametricQuadraticLineElement::JacobianAtIntegrationPoint: ";
         std::cerr <<"Element has only 2 integration points.\n";
         current_detJ = std::numeric_limits<double64>::quiet_NaN();
         return;
      }
      
    dNr( IP[ipoint], DNR );
    if      ( dim == 1 ) current_detJ = Jacobian1D( DNR );
    else if ( dim == 2 ) current_detJ = Jacobian2D( DNR );
    else if ( dim == 3 ) current_detJ = Jacobian3D( DNR );
 }
 
 


double64 IsoparametricQuadraticLineElement::WeightAtIntegrationPoint( size_t i ) const
 {
    if      ( i == 0 ) return W[0];
    else if ( i == 1 ) return W[1];
    
    return std::numeric_limits<double64>::quiet_NaN();
 }



void IsoparametricQuadraticLineElement::N_AtIntegrationPoint( size_t gauss_point, std::vector<double64>& N )
 {
    assert( gauss_point < gpe );

    // local interpolation function values 
    Nr( IP[gauss_point], N );
 }

void IsoparametricQuadraticLineElement::N_AtBaryCenter( std::vector<double64>& N )
  {
    // barycenter is located at 0. as the element extends from -1 to 1 in local coordinate space
    Nr( 0., N );
  }


void  IsoparametricQuadraticLineElement::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 2;
    ids[2] = 1;
 }


void  IsoparametricQuadraticLineElement::CornerNodes( std::vector<size_t>& ids ) const   
 {
    ids.resize(2);
    ids[0] = 0;
    ids[1] = 1;
 }



void  IsoparametricQuadraticLineElement::MidSideNodes( std::vector<size_t>& ids ) const   
 {
    ids.resize(1);
    ids[0] = 2;
 }


/// segments are numbered like faces
void IsoparametricQuadraticLineElement::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    if ( segm_id > 1U )
      std::cout <<"\nIsoparametricQuadraticLineElement::NodesOfSegment: There is only one segment present."<< std::endl;
    snids.resize(3);
    snids[0] = 0;
    snids[1] = 1;  
    snids
      [2] = 2; // midside node comes last
 }


/**
    Method returns into its argument vector the local node number of either of its 2 faces located at its nodes.
    Convention: face 0 has only one node which is the first node of the element and face 1 contains the second node.
*/
void  IsoparametricQuadraticLineElement::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    if ( face_id > 1U )
      std::cout <<"\nIsoparametricQuadraticLineElement::NodesOfFace: There are only 2 faces present."<< std::endl;
    fnids.resize(1U);
    fnids[0] = face_id;
 }


/// assuming straight segments between the nodes
void  IsoparametricQuadraticLineElement::EdgeLengths( vector<double64>& vec ) 
 {
    vec.resize(2);
    
    if ( dim == 1U ) {
         vec[0] = fabs(XY(2,0)-XY(0,0)); 
         vec[1] = fabs(XY(1,0)-XY(2,0)); 
         return;
      }
    
    if ( dim == 2U ) {
	     // segment node 1,3
	     double64 lenx = XY(2,0) - XY(0,0);
	     double64 leny = XY(2,1) - XY(0,1);
	     vec[0] = hypot( lenx, leny );
	     // segment node 3,2
	     lenx = XY(1,0) - XY(2,0);
	     leny = XY(1,1) - XY(2,1);
	     vec[1] = hypot( lenx, leny );
         return;
      }

    if ( dim == 3U ) {
	     // segment node 1,3
	     double64 lenx = XY(2,0) - XY(0,0);
	     double64 leny = XY(2,1) - XY(0,1);
	     double64 lenz = XY(2,2) - XY(0,2);
	     vec[0] = sqrt( lenx*lenx + leny*leny + lenz*lenz );
	     // segment node 3,2
	     lenx = XY(1,0) - XY(2,0);
	     leny = XY(1,1) - XY(2,1);
	     lenz = XY(1,2) - XY(2,2);
	     vec[1] = sqrt( lenx*lenx + leny*leny + lenz*lenz );
         return;
      }

 } // end EdgeLengths





/** Compute derivatives of interpolation functions at corresponding nodes with
respect to the global coordinate system.  

It is assumed that the element is a straight line in space.

@section arguments Input Arguments 

The element is used to obtain the global interpolation of the tetrahedron and the 
interpolation function derivatives at each node are returned into the 
matrix M of dimensions rows = spatial dimensions x columns = nodes.  

@param DN The interpolation-function derivative matrix is returned into the
second method argument. 
 
 */
void IsoparametricQuadraticLineElement::dN( DenseMatrix<DM_MIN>& DN ) 
  { 
     DN.Resize(dim,npe);
     // using M matrix for temporary storage
     dN_AtNode( M, 0 );
     // creating entries for first column of DN
     for ( size_t i=0; i<dim; i++ ) DN(i,0) = M(i,0);
     dN_AtNode( M, 1 );
     // creating entries for first column of DN
     for ( size_t i=0; i<dim; i++ ) DN(i,1) = M(i,1);
     dN_AtNode( M, 2 );
     // creating entries for first column of DN
     for ( size_t i=0; i<dim; i++ ) DN(i,2) = M(i,2);
             
  } // end dN
  




/**
 
@section arguments Input Arguments 

The finite element which supplies the nodal coordinates from which
the area is computed. 

@return The volume (m3) of the finite element.  
 */
double64 IsoparametricQuadraticLineElement::Volume()
{
   double64  len(0.);

   // 2-dimensional models
   if ( dim == 2 )
     for ( size_t i=0; i<gpe; i++ ) {
          dNr( IP[i], DNR );
          len += W[i] * Jacobian2D( DNR );
       }
     
   // 3-dimensional models
   else if ( dim == 3 )
     for ( size_t i=0; i<gpe; i++ ) {
          dNr( IP[i], DNR );
          len += W[i] * Jacobian3D( DNR );
       }
       
   // in this case the length of the element is equivalent to
   // the sum of the Jacobian determinants at the integration points
   else if ( dim == 1 ) 
     for ( size_t i=0; i<gpe; i++ ) {
          dNr( IP[i], DNR );
          len += W[i] * Jacobian1D( DNR );
       }
      
    return len;
}

 
 
 



/**

JacobianFor( DNR, x) functions return the determinant of the Jacobian matrices
for the spatial dimension x

The inverse of a matrix A-1 is equal to 1/detJ(A) * A^T (transposed).
It is used to compute the interpolation function derivatives in
global coordinates as

DN_global = DN^T J-1^T  ops J-1 DN
*/
double64 IsoparametricQuadraticLineElement::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& DN,
                                                       size_t gauss_point )
 {
     DN.Resize(dim,npe);
     
     // see Cook et al., p. 166 and discussion Cheung et al. p. 26
     if ( dim == 1 ) {
          // node 1
          dNr( IP[gauss_point], DNR );
          double64 dxdr = JacobianFor( DNR, 0 );
          double64 detJinv = 1. / dxdr;
          // as in Cook et al. p. 166, eqn 6.2-7
          DN(0,0) = DNR[0] * detJinv; 
          // node 2
          DN(0,1) = DNR[1] * detJinv; 
          // node 3
          DN(0,2) = DNR[2] * detJinv; 
          return 1. / detJinv; // returns the determinant of J
       }
     
     if ( dim == 2 ) {
          // the dxi_i need to be multiplied with global coordinates
          dNr( IP[gauss_point], DNR );
          // getting dx/dr
          double64 dxdr = JacobianFor( DNR, 0 );
          double64 dydr = JacobianFor( DNR, 1 );
          // 1 / detJ = detJinv
          double64 detJinv = 1. / hypot( dxdr, dydr );
          double64 cosa = detJinv * dxdr;
          double64 sina = detJinv * dydr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv; 
          DN(0,1) = DNR[1] * cosa * detJinv;
          DN(0,2) = DNR[2] * cosa * detJinv; 
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv; 
          DN(1,1) = DNR[1] * sina * detJinv; 
          DN(1,2) = DNR[2] * sina * detJinv; 
          return 1. / detJinv; // detJ
       }
       
     if ( dim == 3 ) {
          dNr( IP[gauss_point], DNR );
          double64 dxdr = JacobianFor( DNR, 0 );
          double64 dydr = JacobianFor( DNR, 1 );
          double64 dzdr = JacobianFor( DNR, 2 );
          double64 detJinv = 1. / sqrt( dxdr*dxdr + dydr*dydr + dzdr*dzdr );
          double64 cosa = detJinv * dxdr;
          double64 sina = detJinv * dydr;
          double64 sinb = detJinv * dzdr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv; 
          DN(0,1) = DNR[1] * cosa * detJinv;
          DN(0,2) = DNR[2] * cosa * detJinv; 
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv; 
          DN(1,1) = DNR[1] * sina * detJinv; 
          DN(1,2) = DNR[2] * sina * detJinv; 
          // row 3 = z-derivatives
          DN(2,0) = DNR[0] * sinb * detJinv; 
          DN(2,1) = DNR[1] * sinb * detJinv; 
          DN(2,2) = DNR[2] * sinb * detJinv; 
          return 1. / detJinv;
       }
             
    throw logic_error("IsoparametricQuadraticLineElement::dN_AtIntegrationPoint: Dimension not recognised");
    
 } // end dN_AtIntegrationPoint







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
double64 IsoparametricQuadraticLineElement::dN_AtNode( DenseMatrix<DM_MIN>& DN, size_t nd )
 {
     DN.Resize(dim,npe);
     
     // see Cook et al., p. 166 and discussion Cheung et al. p. 26
     if ( dim == 1 ) {
          // node 1
          dNr( NX[nd], DNR );
          double64 detJinv = 1. / JacobianFor( DNR, 0 );
          DN(0,0) = DNR[0] * detJinv; 
          // node 2
          DN(0,1) = DNR[1] * detJinv; 
          // node 3
          DN(0,2) = DNR[2] * detJinv; 
          return 1. / detJinv;
       }
       
     if ( dim == 2 ) {
          dNr( NX[nd], DNR );
          // getting dx/dr
          double64 dxdr = JacobianFor( DNR, 0 );
          double64 dydr = JacobianFor( DNR, 1 );
          // 1 / detJ = detJinv
          double64 detJinv = 1. / hypot( dxdr, dydr );
          double64 cosa = detJinv * dxdr;
          double64 sina = detJinv * dydr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv; 
          DN(0,1) = DNR[1] * cosa * detJinv;
          DN(0,2) = DNR[2] * cosa * detJinv; 
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv; 
          DN(1,1) = DNR[1] * sina * detJinv; 
          DN(1,2) = DNR[2] * sina * detJinv; 
          return 1. / detJinv; // detJ
       }
       
     if ( dim == 3 ) {
          dNr( NX[nd], DNR );
          double64 dxdr = JacobianFor( DNR, 0 );
          double64 dydr = JacobianFor( DNR, 1 );
          double64 dzdr = JacobianFor( DNR, 2 );
          double64 detJinv = 1. / sqrt( dxdr*dxdr + dydr*dydr + dzdr*dzdr );
          double64 cosa = detJinv * dxdr;
          double64 sina = detJinv * dydr;
          double64 sinb = detJinv * dzdr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv; 
          DN(0,1) = DNR[1] * cosa * detJinv;
          DN(0,2) = DNR[2] * cosa * detJinv; 
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv; 
          DN(1,1) = DNR[1] * sina * detJinv; 
          DN(1,2) = DNR[2] * sina * detJinv; 
          // row 3 = z-derivatives
          DN(2,0) = DNR[0] * sinb * detJinv; 
          DN(2,1) = DNR[1] * sinb * detJinv; 
          DN(2,2) = DNR[2] * sinb * detJinv; 
          return 1. / detJinv;
       }
       
    throw logic_error("IsoparametricQuadraticLineElement::dN_AtNode: Dimension not recognised");
       
 } // end dN_AtNode
 
 
 



double64 IsoparametricQuadraticLineElement::dN_AtBarycenter( DenseMatrix<DM_MIN>& DN )
 {
    // the middle node is equivalent to the barycenter
    return dN_AtNode( DN, 2 );
 }



/// here the normal is calculated using the derivative of the shape function
/// at the middle node
void  IsoparametricQuadraticLineElement::UnitNormal( vector<double64>& vc ) const
 {
    if ( dim == 1 ) {
         vc[0] = 1.;
         return;
      }
    
    // dx = (N1'x1 + N2'x2 + N3'x3) dr
    if ( dim == 2 ) {
         // finding tangent at mid-point node
         dNr( NX[2], DNR );
         vc[0] = JacobianFor( DNR, 0 ); // dx
         vc[1] = JacobianFor( DNR, 1 ); // dy
         mjl::Edge  normal( mjl::Point(0.,0.), mjl::Point(vc[0],vc[1]) ); 
         // rotating tangent edge counter-clockwise to find normal to face
         normal.Rot();
         normal.NormalizeTo( 1. );
         vc[0] = normal.Destination()[0];
         vc[1] = normal.Destination()[1];
         return;
      }

    if ( dim == 3 ) {
         // using slope at mid-point node
         dNr( NX[2], DNR );
         vc[0] = JacobianFor( DNR, 0 ); // dx
         vc[1] = JacobianFor( DNR, 1 ); // dy
         vc[2] = JacobianFor( DNR, 2 ); // dz
         double64 sum = vc[0] + vc[1] + vc[2];
         // normalizing the normal
         vc[0] /= sum;
         vc[1] /= sum;
         vc[2] /= sum;
         cout <<"\nIsoparametricQuadraticLineElement::UnitNormal: In 3D a reference direction is needed to find normal.\n";
         return;
      }

 } // end UnitNormal



/**

This assumes that the interpolated data are supplied as
a symmetric material element property matrix (this analytic expression
exists only if the values do not vary over the nodes)
EPROP can also contain normal and tangential forces such that

           ft fn
  EPROP = -fn ft

analytic solution, eqn 2.9, Cheung et al. p. 26 */
void  IsoparametricQuadraticLineElement::IntegralN( DenseMatrix<DM_MIN>& EPROP )
 {
    // transpose EPROP
    EPROP.Transposed( M );
    
    // using EPROP to hold coefficients for multiplication
    if ( dim == 2 ) {
		// no Jacobian required since this is an analytical solution
		EPROP.Resize(dim,npe);
		EPROP(0,0) = 4.*XY(2,0)-XY(0,0)-XY(1,0); 
		EPROP(0,1) = 4.*(XY(1,0)-XY(0,0)); 
		EPROP(1,0) = 4.*XY(2,1)-XY(0,1)-XY(1,1); 
		EPROP(1,1) = 4.*(XY(1,1)-XY(0,1)); 
        EPROP(0,2) = XY(0,0)-4.*XY(2,0)+3.*XY(1,0);
        EPROP(1,2) = XY(0,1)-4.*XY(2,1)+3.*XY(1,1);
        
		M *= EPROP;

		EPROP.Resize(dim,npe);
		for ( size_t i=0; i<dim; i++ )
		  for ( size_t j=0; j<npe; j++ ) EPROP(i,j) = M(i,j) * 1. / 6.;
  
        return;
     }
      
 } // end IntegralN
 

 
void  IsoparametricQuadraticLineElement::ExtrapolateIntegrationPointVariableToNodes( size_t nvars, 
                                                                  const std::vector<double64>& IVAR, 
                                                                  std::vector<double64>& NVAR ) const
{
   // 0. Decide which case is dealt with in terms of the integration points 
   //    which are used (rr and ss contain the integr.p. locations)
   assert( IVAR.size() >= (gpe*nvars) );
   NVAR.resize( npe * nvars );

   // 1. Compute the slope of the linear interpolation function defined by the 2 integration points.
   double64 dx  = IP[1] - IP[0];
   double64 xn1 = fabs(NX[0]) - fabs(IP[0]);
   double64 xn2 = NX[1] - IP[1];
   double64 xn3 = NX[2] - IP[0];
   
   // 2. For each node point compute the values of the interpolation functions
   //    and use these to extrapolate the values of the variables at the nodes.
        // compute interpolation function values at node i
   // node 1
   for ( size_t k=0; k<nvars; k++ ) {
        // compute the slope of the interpolation function
        double64 m = (IVAR[1*nvars + k] - IVAR[0*nvars + k]) / dx;
        // extrapolate value from IP1 to node 1
        NVAR[0*nvars + k] = IVAR[0*nvars + k] - m * xn1;
        // node 2
        NVAR[1*nvars + k] = IVAR[1*nvars + k] + m * xn2;
        // node 3
        NVAR[2*nvars + k] = IVAR[0*nvars + k] + m * xn3;
     }

 } // end ExtrapolateIntegrationPointVariableToNodes 






/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricQuadraticLineElement::IntegrationPoint( size_t i, vector<double64>& xyz ) const
 {
    assert( i < gpe );
    xyz.resize(dim); 
    xyz[0]=0.;  

    // local interpolation function values 
    Nr( IP[i], NRST );

    // 1D
    if ( dim == 1U ) {
         for( size_t i=0U; i<npe; i++ )
           xyz[0] += XY(i,0) * NRST[i];
	     return;
      }
    
    // 2D
    if ( dim == 2U ) {
         xyz[1]=0;
         for( size_t i=0U; i<npe; i++ ) {
	          xyz[0] += XY(i,0) * NRST[i];
	          xyz[1] += XY(i,1) * NRST[i];
	       }
	     return;
      }
   
    // 3D case   
    xyz[1]=xyz[2]=0.; 
    for( size_t i=0U; i<npe; i++ ) {
	      xyz[0] += XY(i,0) * NRST[i];
	      xyz[1] += XY(i,1) * NRST[i];
	      xyz[2] += XY(i,2) * NRST[i];
	  }

 } // end IntegrationPoint





void  IsoparametricQuadraticLineElement::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                        vector<size_t>& fnids )
 {
    fnids.resize(bnodes.size());

     // no ordering has to be established in these cases
     if ( bnodes.size() == npe || bnodes.size() == 1 ) {
          fnids = bnodes;
          return;
       }
     
     if ( bnodes.size() == 2 ) {
          if ( bnodes[0] == 0 && bnodes[1] == 1 )
            throw csmp::Exception( CSMP_ERROR, "IsoparametricQuadraticLineElement::ConsecutiveNodesAtBoundary", 
                           "While corner nodes lie on boundary, midside node does not ?" );
          if ( bnodes[0] == 0 && bnodes[1] == 2 )
            throw csmp::Exception( CSMP_ERROR, "IsoparametricQuadraticLineElement::ConsecutiveNodesAtBoundary", 
                           "First 2 nodes lie on boundary, last node does not ?" );
          if ( bnodes[0] == 1 && bnodes[1] == 2 )
            throw csmp::Exception( CSMP_ERROR, "IsoparametricQuadraticLineElement::ConsecutiveNodesAtBoundary", 
                           "Last 2 nodes lie on boundary, first node does not ?" );
          return; 
       }
            
 } // end ConsecutiveNodesAtBoundary









void IsoparametricQuadraticLineElement::OutputNodeDataToVTK( const char* file_name, 
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
           cout <<"\nIsoparametricQuadraticLineElement::OutputNodeDataToVTK "; 
           cout <<"Output file could not be opened."<< endl;
           return;
       }  
       
     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): IsoparametricQuadraticLineElement: Variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" float"<< endl;
     for ( size_t i=0; i<npe; i++ ) {
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" "; 
          if ( dim == 2 ) ofs << 0.0 <<" ";
          ofs << endl;
       }
     ofs << endl;  
       
     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< npe+1 << endl;
     // number of points per cell, point1, ... point n
     ofs << npe <<" 0 2 1";
     ofs << endl;
     ofs << endl;
     
     // 4. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );
     
     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" float"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1 x n_nodes
           for ( size_t i=0; i<npe; i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x n_nodes
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
               for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               if ( dim == 3 )
                 ofs << endl;
               else
                 ofs << 0.0 << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricQuadraticLineElement::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK

void IsoparametricQuadraticLineElement::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
	matCoords.Resize(npe, dim);
	matCoords.Fill(0.);
 	
 	matCoords(0,0) = NX[0];
 	matCoords(1,0) = NX[1];
 	matCoords(2,0) = NX[2];
}

} // end namespace csp



