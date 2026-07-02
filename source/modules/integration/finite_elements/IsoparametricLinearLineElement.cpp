#include "IsoparametricLinearLineElement.h"

using namespace std;

namespace csmp {


IsoparametricLinearLineElement::IsoparametricLinearLineElement( uint32_t dimensions, uint32_t ips )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_BAR, true, true, 1U ),
    current_detJ(std::numeric_limits<double>::quiet_NaN())
 {
    dim = dimensions;
    itp = 1;
    npf = 2;
    npe = 2;
    fpe = 2; // 2-faces located at the two nodes matching the location of the neighbors
    spe = 1;
    epe = 2;
    nne = 2;
    cne = 0;
    gpe = ips;

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
    ElementType(ISOPARAMETRIC_LINEAR_BAR);

    IP.resize( gpe );
    W.resize( gpe );

    if(ips==1)
    {
        IP[0]   =  0.0;
        W[0]    =  2.0;
    }
    else if(ips==2)
     {
        IP[0]   = -1. / sqrt(3.);
        IP[1]   = 1. / sqrt(3.);
        W[0]    =  1.;            
        W[1]    = 1.; // since length = 2.0 (-1,1)
     }
    else
        cout<<" IsoparametricLinearLineElement::IsoparametricLinearLineElement "<<ips<<" integration points not supported " <<endl;

    // local node coordinates (node order 1,2)
    NX[0] = -1;
    NX[1] =  1.;
}


CSMP_FEM_TYPE  IsoparametricLinearLineElement::ElementTypeOfFace( uint32_t ) const
 {
    return ZERO_DIMENSIONAL_FACE;
 }

// tested: OK1
void IsoparametricLinearLineElement::Nr( double r, std::vector<double>& N ) const
{
   N.resize(npe);
   N[0] = 0.5 * (1.-r );
   N[1] = 0.5 * (1.+r );
}

void IsoparametricLinearLineElement::Nr( double r, double* N ) const
{
   N[0] = 0.5 * (1.-r );
   N[1] = 0.5 * (1.+r );
}



// tested: OK1
void IsoparametricLinearLineElement::dNr( double, std::vector<double>& dnr ) const
{
   dnr[0] = -0.5;
   dnr[1] =  0.5;
}



// tested: OK1
double IsoparametricLinearLineElement::JacobianFor( const std::vector<double>& DNR, uint32_t coord ) const
 {
    assert( coord < dim );
    return DNR[0] * XY(0,coord) + DNR[1] * XY(1,coord);
 }



// these next 3 auxiliary methods are all used by the Volume() function
// tested:
double IsoparametricLinearLineElement::Jacobian1D( const std::vector<double>& DNR ) const
 {
    return DNR[0] * XY(0,0) + DNR[1] * XY(1,0);
 }



double IsoparametricLinearLineElement::Jacobian2D( const std::vector<double>& DNR ) const
 {
    double jac_x = DNR[0] * XY(0,0) + DNR[1] * XY(1,0);
    double jac_y = DNR[0] * XY(0,1) + DNR[1] * XY(1,1);

    // compute a length-like expression for the transformation
    return std::sqrt( jac_x*jac_x + jac_y*jac_y );

 } // end



double IsoparametricLinearLineElement::Jacobian3D( const std::vector<double>& DNR ) const
 {
    double jac_x = DNR[0] * XY(0,0) + DNR[1] * XY(1,0);
    double jac_y = DNR[0] * XY(0,1) + DNR[1] * XY(1,1);
    double jac_z = DNR[0] * XY(0,2) + DNR[1] * XY(1,2);

    return std::sqrt( jac_x*jac_x + jac_y*jac_y + jac_z*jac_z );

 } // end


// Next 2 methods belong together:
// return determinant J for values of previous function
double IsoparametricLinearLineElement::JacobianInverse() { return current_detJ; }

void IsoparametricLinearLineElement::JacobianAtIntegrationPoint( uint32_t ipoint )
 {
    if ( ipoint > 1U ) {
         std::cout <<"\nIsoparametricLinearLineElement::JacobianAtIntegrationPoint: ";
         std::cout <<"Element has only 2 integration points.\n";
         current_detJ = std::numeric_limits<double>::quiet_NaN();
         return;
      }

    dNr( IP[ipoint], DNR );
    if      ( dim == 3 ) current_detJ = Jacobian3D( DNR );
    else if ( dim == 2 ) current_detJ = Jacobian2D( DNR );
    else if ( dim == 1 ) current_detJ = Jacobian1D( DNR );
 }




double IsoparametricLinearLineElement::WeightAtIntegrationPoint( uint32_t i ) const
 {

    if(gpe==1) {
      if      ( i == 0 ) return W[0];
    }
    else {
     if      ( i == 0 ) return W[0];
     else if ( i == 1 ) return W[1];
    }

    return std::numeric_limits<double>::signaling_NaN();
 }


// tested:
void IsoparametricLinearLineElement::N_AtIntegrationPoint( uint32_t gauss_point, std::vector<double>& N )
 {
    assert( gauss_point < gpe );

    // local interpolation function values
    Nr( IP[gauss_point], N );
 }



void IsoparametricLinearLineElement::N_AtBaryCenter( std::vector<double>& N )
 {
    // barycenter is located at 0. as the element extends from -1 to 1 in local coordinate space
    Nr( 0., N );
 }



void  IsoparametricLinearLineElement::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
 }


void  IsoparametricLinearLineElement::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
 }

uint32_t    IsoparametricLinearLineElement::CornerNodes() const { return 2U; }

uint32_t    IsoparametricLinearLineElement::MidSideNodes() const { return 0U; }



void  IsoparametricLinearLineElement::MidSideNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(1);
    std::cout <<"\nIsoparametricLinearLineElement::MidSideNodes No midside nodes present"<<std::endl;
 }

// segments are numbered like faces
void IsoparametricLinearLineElement::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
 {
    if ( segm_id > 1U )
      std::cout <<"\nIsoparametricLinearLineElement::NodesOfSegment: There is only one segment present."<< std::endl;
    snids.resize(npe);
    snids[0] = 0;
    snids[1] = 1;
 }




/**
    Method returns into its argument vector the local node number of either of its 2 faces located at its nodes.
    Conventions: 
       - the faces have only a single node
       - node 0 corresponds to first face 0 opposite to it
       - node 1 corresponds to second face 1 opposite to node 1
*/
/*
void  IsoparametricLinearLineElement::NodesOfFace( uint32_t face_id,
                                                   vector<uint32_t>& fnids ) const
 {
    assert( face_id <= 1U );
    fnids.resize(1U);
    fnids[0] = (face_id == 1U) ? 0U : 1U;
 }
*/


vector<uint32_t>  IsoparametricLinearLineElement::NodesOfFace( uint32_t face_id ) const {
    return (face_id == 1U) ? vector<uint32_t>{0U} : vector<uint32_t>{1U};
 }


vector<uint32_t>  IsoparametricLinearLineElement::CornerNodesOfFace( uint32_t face_id ) const
 {
    assert( face_id <= 1 );
    if ( face_id == 1U ) return vector<uint32_t>{0U};
    return vector<uint32_t>{1U};
 }



vector<uint32_t>  IsoparametricLinearLineElement::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        // local corner node numbers are returned in ascending order
        case 0: return vector<uint32_t>{1};
        case 1: return vector<uint32_t>{0};
        default:
          cerr <<"\nIsoparametricLinearLineElement::NodesConnectedTo: node "<< node_id <<" does not exist.";
      }
    return vector<uint32_t>{};
  }


// assuming straight segments between the nodes
// tested: OK3
void  IsoparametricLinearLineElement::EdgeLengths( vector<double>& vec )
 {
    vec.resize(1);

    if ( dim == 1U ) {
         vec[0] = fabs(XY(1,0)-XY(0,0));
         return;
      }

    if ( dim == 2U ) {
         // segment node 1,2
         const double lenx = XY(1,0) - XY(0,0);
         const double leny = XY(1,1) - XY(0,1);
         vec[0] = hypot( lenx, leny );
         return;
      }

    if ( dim == 3U ) {
         // segment node 1,2
         const double lenx = XY(1,0) - XY(0,0);
         const double leny = XY(1,1) - XY(0,1);
         const double lenz = XY(1,2) - XY(0,2);
         vec[0] = sqrt( lenx*lenx + leny*leny + lenz*lenz);
         return;
      }

 } // end EdgeLengths


/**

  double IsoparametricLinearLineElement::Volume()



Description:


@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The volume (m3) of the finite element.

@test OK 

*/
double IsoparametricLinearLineElement::Volume()
{
   double  len(0.);

   // 2-dimensional models
   if ( dim == 2U )
     for ( uint32_t i{0U}; i<gpe; i++ ) {
          dNr( IP[i], DNR );
          len += W[i] * Jacobian2D( DNR );
       }

   // 3-dimensional models
   else if ( dim == 3U )
     for ( uint32_t i{0U}; i<gpe; i++ ) {
          dNr( IP[i], DNR );
          len += W[i] * Jacobian3D( DNR );
       }

   // in this case the length of the element is equivalent to
   // the sum of the Jacobian determinants at the integration points
   else if ( dim == 1U )
     for ( uint32_t i{0U}; i<gpe; i++ ) {
          dNr( IP[i], DNR );
          len += W[i] * Jacobian1D( DNR );
       }

    return len;
}


// assuming unit (missing) dimensions, the volume is equivalent to the length/1
double IsoparametricLinearLineElement::AspectRatio()
 {
    return Volume();
 }



double IsoparametricLinearLineElement::InnerRadius()
{

    return Volume()/2.0;

}
/**

  void LinearLineElement::N( const Element& e, vector<double>& FN, const vector<double>& xy )



Description:

Computes the values of the interpolation functions at the global coordinate 'xyz'.
Since a Jacobian matrix can only be obtained for points which lie inside
the local tetrahedron, such a transformation cannot be performed here.
Because of this reason global area-coordinate functions are employed
to obtain the interpolation function values at the global point.

Note, that the use of global interpolation functions implies that the method
only gives highly accurate results if the element boundaries are
not curved.

@section arguments Input Arguments

The parent Element, the vector which will hold the interpolation function
values, and the coordinate vector 'xyz' in the global coordinate system.
Note that if the global coordinates given by 'xyz' do not lie within the
global extent of the quadratic triangular element, the interpolation function values
will no longer add to one and the interpolation will then be erroneous.

The second method argument, the vector FN, will be filled with the interpolation function
values as computed at the point 'xyz'.

@section implementation Implementation

The inverse of the Jacobian matrix is found for the global point 'xyz'.
Then 'xyz' transposed is pre-multiplied with the Jacobian to find the
local coordinate pair 'rs' corresponding to 'xyz'. Using these local
coordinates the interpolation function values are found.

@section application Application

To interpolate a property value withing the quadratic triangular
element.


tested:   */
void IsoparametricLinearLineElement::N(vector<double>& FN, const vector<double>& xy)
{
    // The element has two nodes (npe = 2)
    FN.resize(2); 

    double length_sq = 0.0;
    double dist_to_node1_sq = 0.0;

    // Calculate squared distances to avoid multiple sqrt calls if possible,
    // but for basis functions, we need the actual linear ratio.
    for (uint32_t i = 0; i < dim; ++i) {
        double d_edge = XY(1, i) - XY(0, i); // Vector from node 0 to node 1
        double d_p1   = XY(1, i) - xy[i];    // Vector from input point to node 1
        
        length_sq        += d_edge * d_edge;
        dist_to_node1_sq += d_p1 * d_p1;
    }

    double len = std::sqrt(length_sq);
    double l1  = std::sqrt(dist_to_node1_sq);

    // Safeguard against zero-length elements
    if (len < 1e-14) {
        FN[0] = 0.5;
        FN[1] = 0.5;
        return;
    }

    // Basis functions for a linear line element:
    // N0 is 1 at node 0 and 0 at node 1 (l1/len)
    // N1 is 0 at node 0 and 1 at node 1 ((len-l1)/len)
    FN[0] = l1 / len;
    FN[1] = 1.0 - FN[0]; 
}




/**

  void IsoparametricLinearLineElement::dN( DenseMatrix<DM_MIN>& M )



Description:

Compute derivatives of interpolation functions at corresponding nodes with
respect to the global coordinate system.

It is assumed that the element is a straight line in space.

@section arguments Input Arguments

The element is used to obtain the global interpolation of the bar and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

The interpolation-function derivative matrix is returned into the
second method argument.

@test: OK
*/
void IsoparametricLinearLineElement::dN( DenseMatrix<DM_MIN>& DN )
  {
     DN.Resize(dim,npe);
     // using M matrix for temporary storage
     dN_AtNode( M, 0 );
     // creating entries for first column of DN
     for ( uint32_t i{0U}; i<dim; i++ ) DN(i,0) = M(i,0);
     dN_AtNode( M, 1 );
     // creating entries for first column of DN
     for ( uint32_t i{0U}; i<dim; i++ ) DN(i,1) = M(i,1);

  } // end dN






// NOTE:  JacobianFor( DNR, x) functions return the determinant of the Jacobian matrices
//        for the spatial dimension x
//        The inverse of a matrix A-1 is equal to 1/detJ(A) * A^T (transposed).
//        It is used to compute the interpolation function derivatives in
//        global coordinates as
//
//        DN_global = DN^T J-1^T  ops J-1 DN
//
// tested: OK2
double IsoparametricLinearLineElement::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& DN,
                                                                 uint32_t gauss_point )
 {
     DN.Resize(dim,npe);

     // see Cook et al., p. 166 and discussion Cheung et al. p. 26
     if ( dim == 1U ) {
          // node 1
          dNr( IP[gauss_point], DNR );
          double dxdr = JacobianFor( DNR, 0 );
          double detJinv = 1. / dxdr;
          // as in Cook et al. p. 166, eqn 6.2-7
          DN(0,0) = DNR[0] * detJinv;
          // node 2
          DN(0,1) = DNR[1] * detJinv;
          return 1. / detJinv; // returns the determinant of J
       }

     if ( dim == 2U ) {
          // the dxi_i need to be multiplied with global coordinates
          dNr( IP[gauss_point], DNR );
          // getting dx/dr
          double dxdr = JacobianFor( DNR, 0 );
          double dydr = JacobianFor( DNR, 1 );
          // 1 / detJ = detJinv
          double detJinv = 1. / hypot( dxdr, dydr );
          double cosa = detJinv * dxdr;
          double sina = detJinv * dydr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv;
          DN(0,1) = DNR[1] * cosa * detJinv;
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv;
          DN(1,1) = DNR[1] * sina * detJinv;
          return 1. / detJinv; // detJ
       }

     if ( dim == 3U ) {
          dNr( IP[gauss_point], DNR );
          double dxdr = JacobianFor( DNR, 0 );
          double dydr = JacobianFor( DNR, 1 );
          double dzdr = JacobianFor( DNR, 2 );
          double detJinv = 1. / sqrt( dxdr*dxdr + dydr*dydr + dzdr*dzdr );
          double cosa = detJinv * dxdr;
          double sina = detJinv * dydr;
          double sinb = detJinv * dzdr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv;
          DN(0,1) = DNR[1] * cosa * detJinv;
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv;
          DN(1,1) = DNR[1] * sina * detJinv;
          // row 3 = z-derivatives
          DN(2,0) = DNR[0] * sinb * detJinv;
          DN(2,1) = DNR[1] * sinb * detJinv;
          return 1. / detJinv;
       }

    throw logic_error("IsoparametricLinearLineElement::dN_AtIntegrationPoint: Dimension not recognised");

 } // end dN_AtIntegrationPoint



/**

  double IsoparametricLinearLineElement::dN_AtNode( DenseMatrix<DM_MIN>& B, uint32_t nd );                                       int dof )



Description:

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

The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.

@test ok

*/
double IsoparametricLinearLineElement::dN_AtNode( DenseMatrix<DM_MIN>& DN, uint32_t nd )
 {
     DN.Resize(dim,npe);

     // see Cook et al., p. 166 and discussion Cheung et al. p. 26
     if ( dim == 1 ) {
          // node 1
          dNr( NX[nd], DNR );
          double detJinv = 1. / JacobianFor( DNR, 0 );
          DN(0,0) = DNR[0] * detJinv;
          // node 2
          DN(0,1) = DNR[1] * detJinv;
          return 1. / detJinv;
       }
     if ( dim == 2 ) {
          dNr( NX[nd], DNR );
          // getting dx/dr
          double dxdr = JacobianFor( DNR, 0 );
          double dydr = JacobianFor( DNR, 1 );
          // 1 / detJ = detJinv
          double detJinv = 1. / hypot( dxdr, dydr );
          double cosa = detJinv * dxdr;
          double sina = detJinv * dydr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv;
          DN(0,1) = DNR[1] * cosa * detJinv;
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv;
          DN(1,1) = DNR[1] * sina * detJinv;
          return 1. / detJinv; // detJ
       }

     if ( dim == 3 ) {
          dNr( NX[nd], DNR );
          double dxdr = JacobianFor( DNR, 0 );
          double dydr = JacobianFor( DNR, 1 );
          double dzdr = JacobianFor( DNR, 2 );
          double detJinv = 1. / sqrt( dxdr*dxdr + dydr*dydr + dzdr*dzdr );
          double cosa = detJinv * dxdr;
          double sina = detJinv * dydr;
          double sinb = detJinv * dzdr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv;
          DN(0,1) = DNR[1] * cosa * detJinv;
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv;
          DN(1,1) = DNR[1] * sina * detJinv;
          // row 3 = z-derivatives
          DN(2,0) = DNR[0] * sinb * detJinv;
          DN(2,1) = DNR[1] * sinb * detJinv;
          return 1. / detJinv;
       }

    throw logic_error("IsoparametricLinearLineElement::dN_AtNode: Dimension not recognised");

 } // end dN_AtNode





// tested OK1
double IsoparametricLinearLineElement::dN_AtBarycenter( DenseMatrix<DM_MIN>& DN )
 {
     DN.Resize(dim,npe);

     const double rAtBaryCenter(0.0);

     // see Cook et al., p. 166 and discussion Cheung et al. p. 26
     if ( dim == 1U ) {
          // node 1
          dNr( rAtBaryCenter, DNR );
          double dxdr = JacobianFor( DNR, 0 );
          double detJinv = 1. / dxdr;
          // as in Cook et al. p. 166, eqn 6.2-7
          DN(0,0) = DNR[0] * detJinv;
          // node 2
          DN(0,1) = DNR[1] * detJinv;
          return 1. / detJinv; // returns the determinant of J
       }

     if ( dim == 2U ) {
          // the dxi_i need to be multiplied with global coordinates
          dNr( rAtBaryCenter, DNR );
          // getting dx/dr
          double dxdr = JacobianFor( DNR, 0 );
          double dydr = JacobianFor( DNR, 1 );
          // 1 / detJ = detJinv
          double detJinv = 1. / hypot( dxdr, dydr );
          double cosa = detJinv * dxdr;
          double sina = detJinv * dydr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv;
          DN(0,1) = DNR[1] * cosa * detJinv;
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv;
          DN(1,1) = DNR[1] * sina * detJinv;
          return 1. / detJinv; // detJ
       }

     if ( dim == 3U ) {
          dNr( rAtBaryCenter, DNR );
          double dxdr = JacobianFor( DNR, 0 );
          double dydr = JacobianFor( DNR, 1 );
          double dzdr = JacobianFor( DNR, 2 );
          double detJinv = 1. / sqrt( dxdr*dxdr + dydr*dydr + dzdr*dzdr );
          double cosa = detJinv * dxdr;
          double sina = detJinv * dydr;
          double sinb = detJinv * dzdr;
          // row 1 = x-derivatives
          DN(0,0) = DNR[0] * cosa * detJinv;
          DN(0,1) = DNR[1] * cosa * detJinv;
          // row 2 = y-derivatives
          DN(1,0) = DNR[0] * sina * detJinv;
          DN(1,1) = DNR[1] * sina * detJinv;
          // row 3 = z-derivatives
          DN(2,0) = DNR[0] * sinb * detJinv;
          DN(2,1) = DNR[1] * sinb * detJinv;
          return 1. / detJinv;
       }

    throw logic_error("IsoparametricLinearLineElement::dN_AtBarycenter: Dimension not recognised");
 }



/**
    here the normal is calculated using the derivative of the shape function at the middle node
    the normal is the 90o counter-clockwise rotated origin to destination vector.
*/
vector<double> IsoparametricLinearLineElement::UnitNormal() const
{
    if (dim == 1U) return { 1.0 };

    const double rAtBaryCenter = 0.0;
    dNr(rAtBaryCenter, DNR);

    double dx = JacobianFor(DNR, 0);
    double dy = JacobianFor(DNR, 1);
    
    if (dim == 2U) {
        // Updated to CW Rotation: (dx, dy) -> (dy, -dx)
        // This ensures consistency with the Quadratic element and the InterFace logic.
        double nx = dy;
        double ny = -dx;

        double len = std::sqrt(nx * nx + ny * ny);
        if (len > 1e-14) {
            return { nx / len, ny / len };
        }
        return { 0.0, 0.0 }; 
    }

    // 3D Case: Unit Tangent
    double dz = JacobianFor(DNR, 2);
    double len = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    if (len > 1e-14) {
        return { dx / len, dy / len, dz / len };
    }
    return { 0.0, 0.0, 0.0 };
}


/**
    Line element has 2 faces located at the nodes.
    The normals point outward in the direction of the element.
    
    @note convention: face 1 is located at the first node.
*/
void IsoparametricLinearLineElement::UnitNormalToFace(uint32_t face, std::vector<double>& unrml) const
{
    assert(face < 2); // A line has exactly 2 faces (endpoints)
    
    const uint32_t d = Dim();
    unrml.assign(d, 0.0);

    // Identify indices for "this" node and the "other" node
    // If face is 0, we use (Node 0 - Node 1)
    // If face is 1, we use (Node 1 - Node 0)
    uint32_t this_node  = (face == 0) ? 0 : 1;
    uint32_t other_node = (face == 0) ? 1 : 0;

    double length_sq = 0.0;
    for (uint32_t i = 0; i < d; ++i) {
        unrml[i] = XY(this_node, i) - XY(other_node, i);
        length_sq += unrml[i] * unrml[i];
    }

    double length = std::sqrt(length_sq);
    if (length > 1e-14) {
        for (uint32_t i = 0; i < d; ++i) {
            unrml[i] /= length;
        }
    } else {
        // Fallback for degenerate element
        if (d > 0) {
            unrml[0] = (face == 0) ? -1.0 : 1.0;
            // Optionally set others to 0.0, though .assign(d, 0.0) already did that.
        }
    }
}


/**
    This assumes that the interpolated data are supplied as
    a symmetric material element property matrix (this analytic expression
    exists only if the values do not vary over the nodes)
    EPROP can also contain normal and tangential forces such that

             ft fn
    EPROP = -fn ft

   analytic solution, eqn 2.9, Cheung et al. p. 26

   @test  OK1
*/
void  IsoparametricLinearLineElement::IntegralN( DenseMatrix<DM_MIN>& EPROP )
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

        M *= EPROP;

        EPROP.Resize(dim,npe);
        for ( uint32_t i{0U}; i<dim; i++ )
          for ( uint32_t j{0U}; j<npe; j++ ) EPROP(i,j) = M(i,j) * 1. / 6.;

        return;
     }
     else
       cerr<<"\nIsoparametricLinearLineElement::IntegralN: dim="<<dim<<" is not supported."<<endl;

 } // end IntegralN





 // not tested yet !!!
 // OK1
void  IsoparametricLinearLineElement::ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                  const std::vector<double>& IVAR,
                                                                  std::vector<double>& NVAR ) const
{
   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

   if(gpe==2) {
     // 1. Compute the slope of the linear interpolation function defined by the 2 integration points.
     double dx  = IP[1] - IP[0];
     double xn1 = fabs(NX[0]) - fabs(IP[0]);
     double xn2 = NX[1] - IP[1];

     // 2. For each node point compute the values of the interpolation functions
     //    and use these to extrapolate the values of the variables at the nodes.
        // compute interpolation function values at node i
     // node 1
     for ( uint32_t k=0; k<nvars; k++ ) {
          // compute the slope of the interpolation function
          double m = (IVAR[1*nvars + k] - IVAR[0*nvars + k]) / dx;
          // extrapolate value from IP1 to node 1
          NVAR[0*nvars + k] = IVAR[0*nvars + k] - m * xn1;
          // node 2
          NVAR[1*nvars + k] = IVAR[1*nvars + k] + m * xn2;
     }
   }
   else
   cerr<<"\tIsoparametricLinearLineElement::ExtrapolateIntegrationPointVariableToNodes: operation not supported for IPS="<< gpe << endl;

} // end ExtrapolateIntegrationPointVariableToNodes





// integration point location transformed into global coordinates
// NB: the matrix XYZ must be uptodate
void  IsoparametricLinearLineElement::IntegrationPoint( uint32_t i,
                                                        vector<double>& xyz ) const
 {
    assert( i < gpe );
    xyz.resize(dim);
    xyz[0]=0.;

    // local interpolation function values
    Nr( IP[i], NRST );

    // 1D
    if ( dim == 1U ) {
         for( i=0; i<npe; i++ )
           xyz[0] += XY(i,0) * NRST[i];
         return;
      }

    // 2D
    if ( dim == 2U ) {
         xyz[1]=0.;
         for( i=0; i<npe; i++ ) {
              xyz[0] += XY(i,0) * NRST[i];
              xyz[1] += XY(i,1) * NRST[i];
           }
         return;
      }

    // 3D case
    xyz[1]=xyz[2]=0.;
    for( i=0; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint





// tested OK1
/*
void  IsoparametricLinearLineElement::ConsecutiveNodesAtBoundary( const vector<uint32_t>& bnodes,
                                                                  vector<uint32_t>& fnids )
 {
     fnids.resize(bnodes.size());

     // no ordering has to be established in these cases
     if ( bnodes.size() == npe || bnodes.size() == 1 ) {
          fnids = bnodes;
          return;
       }

 } // end ConsecutiveNodesAtBoundary
*/



// tested 0k3
void IsoparametricLinearLineElement::OutputNodeDataToVTK( const char* file_name,
                                                const char* var_name,
                                                DenseMatrix<DM_MIN>& DATA ) const
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     snprintf( elmt, sizeof(elmt), "%zu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cerr <<"\nIsoparametricLinearLineElement::OutputNodeDataToVTK ";
           cerr <<"Output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): IsoparametricLinearLineElement: Variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" float"<< endl;
     for ( uint32_t i{0U}; i<npe; i++ ) {
          for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          if ( dim == 2 ) ofs << 0.0 <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< npe+1 << endl;
     // number of points per cell, point1, ... point n
     ofs << npe <<" 0 1";
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
           for ( auto i{0U}; i<npe; i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x n_nodes
          for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
               for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               if ( dim == 3 )
                 ofs << endl;
               else
                 ofs << 0.0 << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricLinearLineElement::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK


void
IsoparametricLinearLineElement::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords(0,0) = NX[0];
    matCoords(1,0) = NX[1];
}

void
IsoparametricLinearLineElement::JacobianAt( const std::vector<double>& rst )
{
    dNr( rst[0], DNR );
    switch(dim)
    {
        case 3:
             current_detJ = Jacobian3D( DNR );
          return;
        case 2:
             current_detJ = Jacobian2D( DNR );
          return;
        case 1:
             current_detJ = Jacobian1D( DNR );
    }
}

double  IsoparametricLinearLineElement::JacobianDeterminant()
{
    return current_detJ;
}

} // end namespace csmp



