// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "LinearLineElement.h"
#include "Point.h"

using namespace std;

namespace csmp {


LinearLineElement::LinearLineElement( uint32_t dimensions )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( LINEAR_BAR, false, false, 1U )
 {
   dim = dimensions;
   itp = 1;
   npf = 2;
   npe = 2;
   fpe = 2;
   spe = 1;
   epe = 2;
   nne = 2;
   cne = 2;
   gpe = 0;
   XY.Resize(npe, dim);
   M.Resize(npe,npe);

   UsesLocalCoordinates(false);
   Isoparametric( false );
   LineElement();
   ElementType(LINEAR_BAR);
 }


/*
void  LinearLineElement::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
 }
*/


void  LinearLineElement::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
 }


/// segments are numbered like faces
void LinearLineElement::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
 {
    if ( segm_id > 1 )
      std::cerr <<"\nLinearLineElement::NodesOfSegment: There is only one segment present."<< std::endl;
    snids.resize(2);
    snids[0] = 0;
    snids[1] = 1;
 }


/// the Face of a line element is located opposite to the node with the same number
/*
void  LinearLineElement::NodesOfFace( uint32_t face_id, std::vector<uint32_t>& fnids ) const
 {
    assert( face_id <= 1 );
    fnids.resize(1U);
    fnids[0] = ( face_id == 0U ) ? 1U : 0U;
 }
*/


vector<uint32_t>  LinearLineElement::NodesOfFace( uint32_t face_id ) const
 {
    assert( face_id <= 1 );
    if ( face_id == 1U ) return vector<uint32_t>{0U};
    return vector<uint32_t>{1U};
 }


vector<uint32_t>  LinearLineElement::CornerNodesOfFace( uint32_t face_id ) const
 {
    assert( face_id <= 1 );
    if ( face_id == 1U ) return vector<uint32_t>{0U};
    return vector<uint32_t>{1U};
 }


vector<uint32_t>  LinearLineElement::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        // local corner node numbers are returned in ascending order
        case 0: return vector<uint32_t>{1};
        case 1: return vector<uint32_t>{0};
        default:
          cerr <<"\nLinearLineElement::NodesConnectedTo: node "<< node_id <<" does not exist.";
      }
    return vector<uint32_t>{};
  }




/**

Compute derivatives of interpolation functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global interpolation of the tetrahedron and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN The interpolation-function derivative matrix is returned into the
second method argument.
*/
void LinearLineElement::dN( DenseMatrix<DM_MIN>& DN )
{
    // Initialize with size 1 even though it's not used, 
    // just to prevent any internal vector access crashes.
    static const vector<double> dummy_coords(1, 0.0);
    dN_At( DN, dummy_coords );
}



/**

Computes derivatives of interpolation functions at a point XYZ given in
global coordinates within the element (also given in global coordinates).
Since this element is linear, the derivatives are constant in the element
and the point location therefore has no influence on the result.

@section arguments Input Arguments

The element is used to obtain the global interpolation of the tetrahedron and the
interpolation function derivatives at the point XYZ are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument and the length of the Line element is returned
as a double variable.
 */
double LinearLineElement::dN_At( DenseMatrix<DM_MIN>& DN, const vector<double>& )
{
    // npe = 2 (nodes per element), dim = 2 or 3
    DN.Resize(dim, npe); 
    
    // Create vector from Org to Dest
    // Using Point<3> as a container even for 2D is often safer in generic code,
    // but here we follow your dim logic.
    double dx = XY(1,0) - XY(0,0);
    double dy = XY(1,1) - XY(0,1);
    double dz = (dim == 3) ? (XY(1,2) - XY(0,2)) : 0.0;

    double lenSq = dx*dx + dy*dy + dz*dz;
    double len = std::sqrt(lenSq);

    // Shape function derivatives (Constant over the element)
    // Node 0 (N1): Derivative is -delta / L^2
    // Node 1 (N2): Derivative is +delta / L^2
    
    DN(0,0) = -dx / lenSq;  DN(0,1) = dx / lenSq; // dN/dx
    DN(1,0) = -dy / lenSq;  DN(1,1) = dy / lenSq; // dN/dy
    
    if (dim == 3) {
        DN(2,0) = -dz / lenSq;  DN(2,1) = dz / lenSq; // dN/dz
    }

    return len;
}




/**

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The volume (m3) of the finite element.
*/
double LinearLineElement::Volume()
{
    // The "Volume" of a line element is its length (L)
    // potentially multiplied by a cross-sectional area if defined.
    
    double dx = XY(1,0) - XY(0,0);
    double dy = XY(1,1) - XY(0,1);
    double dz = (dim == 3) ? (XY(1,2) - XY(0,2)) : 0.0;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}


/**

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

@param FN The second method argument, the vector FN, will hold the interpolation function
values as computed at the point 'xyz'.

@section implementation Implementation

The inverse of the Jacobian matrix is found for the global point 'xyz'.
Then 'xyz' transposed is pre-multiplied with the Jacobian to find the
local coordinate pair 'rs' corresponding to 'xyz'. Using these local
coordinates the interpolation function values are found.

@section application Application

To interpolate a property value withing the quadratic triangular
element.
*/
void LinearLineElement::N( vector<double>& FN, const vector<double>& xy )
 {
     double  len, l1, l2;

     // 2-dimensional edge
     if ( dim == 2 ) {
          Point<2> org(XY(0,0),XY(0,1)), dest(XY(1,0),XY(1,1)), p(xy[0],xy[1]);
          len = (dest - org).Length();
          l1  = (dest - p).Length();
          l2  = len - l1;
       }
     // 3-dimensional edge
     else
       {
          Point<3> org(XY(0,0),XY(0,1),XY(0,2)), dest(XY(1,0),XY(1,1),XY(1,2)), p(xy[0],xy[1],xy[2]);
          len = (dest - org).Length();
          l1  = (dest - p).Length();
          l2  = len - l1;
       }

     FN.resize(npe);
     //before:
     //FN[0] = l2 / len;
     //FN[1] = l1 / len;
     //after:
     FN[0] = l1 / len;
     FN[1] = l2 / len;

 } // end











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
double LinearLineElement::dN_AtNode( DenseMatrix<DM_MIN>& DN, uint32_t )
 {
    dN( DN );
/*
    if ( dof == 2 ) {
         DN.resize(dim+1,npe*2);
         DN(0,3)=0.0,     DN(0,2)=DN(0,1), DN(0,1)=0.0,     DN(0,0)=DN(0,0);
         DN(1,3)=DN(1,1), DN(1,2)=0.0,     DN(1,1)=DN(1,0), DN(1,0)=0.0;
         DN(2,3)=DN(0,1), DN(2,2)=DN(1,1), DN(2,1)=DN(0,0), DN(2,0)=DN(1,0);
      }
*/
    if ( dim == 3U ) {
         cout <<"LinearLineElement::dN_AtNode: Method not implemented for 3D case. "<< endl;
         terminate();
      }

    return 0.;
 }



double LinearLineElement::dN_AtBarycenter( DenseMatrix<DM_MIN>& DN )
 {
    dN( DN );
/*
    if ( dof == 2 ) {
         DN.resize(dim+1,npe*2);
         DN(0,3)=0.0,     DN(0,2)=DN(0,1), DN(0,1)=0.0,     DN(0,0)=DN(0,0);
         DN(1,3)=DN(1,1), DN(1,2)=0.0,     DN(1,1)=DN(1,0), DN(1,0)=0.0;
         DN(2,3)=DN(0,1), DN(2,2)=DN(1,1), DN(2,1)=DN(0,0), DN(2,0)=DN(1,0);
      }
*/
    if ( dim == 3 ) {
         cerr <<"LinearLineElement::dN_AtNode: Method not implemented for 3D case. "<< endl;
         terminate();
      }

    return 0.0;
 }

void LinearLineElement::N_AtBaryCenter(std::vector<double>& N)
{
	vector<double> p;
	p.resize(dim);
	for ( uint32_t i = 0; i < dim; ++i) p[i] = 0.5*(XY(0, i) + XY(1, i));
	LinearLineElement::N(N, p);
}



/**
    In 3D, the input vector must contain the rotation axis around which
    the edge is turned to create the unit normal.
*/
vector<double> LinearLineElement::UnitNormal() const
{
    // 2D Case
    if (dim == 2) {
        // 1. Calculate direction vector (Tangent)
        double dx = XY(1,0) - XY(0,0);
        double dy = XY(1,1) - XY(0,1);

        // 2. Rotate 90 degrees CLOCKWISE (dy, -dx)
        // To match the InterFace test results and Quadratic consistency
        double nx = dy;
        double ny = -dx;

        // 3. Normalize
        double length = std::sqrt(nx * nx + ny * ny);
        if (length > 1e-14) {
            return { nx / length, ny / length };
        }
        return { 0.0, 0.0 };
    }

    // 3D Case
    // Standardizing to return the Unit Tangent for 1D elements in 3D space,
    // as a unique normal is undefined without a reference orientation.
    double dx = XY(1,0) - XY(0,0);
    double dy = XY(1,1) - XY(0,1);
    double dz = XY(1,2) - XY(0,2);

    double length = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (length > 1e-14) {
        return { dx / length, dy / length, dz / length };
    }
    
    // If it's truly degenerate, we'll throw or return a default
    throw std::runtime_error("LinearLineElement::UnitNormal: Element has zero length.");
}



void  LinearLineElement::IntegralN( DenseMatrix<DM_MIN>& M )
 {
    M.Zero();
    M(0,0) = M(1,1) = Volume();
 }


void  LinearLineElement::IntegralNN( DenseMatrix<DM_MIN>& M )
 {
    M(0,1) = M(1,0) = Volume();
    M(0,0) = M(1,1) = 2.0 * M(0,1);
 }




/**
    Line element has 2 faces located at the nodes.
    The normals point outward in the direction of the element.
    
    @note convention: face 1 is located at the first node.
*/
void LinearLineElement::UnitNormalToFace(uint32_t face, std::vector<double>& unrml) const
{
    assert(face < 2); // Line elements have exactly 2 faces (endpoints)

    const uint32_t d = Dim();
    unrml.assign(d, 0.0);

    // If face is 0, vector is (Node 0 - Node 1)
    // If face is 1, vector is (Node 1 - Node 0)
    const uint32_t this_node  = (face == 0) ? 0 : 1;
    const uint32_t other_node = (face == 0) ? 1 : 0;

    double length_sq = 0.0;
    for (uint32_t i = 0; i < d; ++i) {
        unrml[i] = XY(this_node, i) - XY(other_node, i);
        length_sq += unrml[i] * unrml[i];
    }

    const double length = std::sqrt(length_sq);
    
    if (length > 1e-14) {
        for (uint32_t i = 0; i < d; ++i) {
            unrml[i] /= length;
        }
    } else {
        // Fallback for degenerate (zero-length) element
        if (d > 0) unrml[0] = (face == 0) ? -1.0 : 1.0;
    }
}



//
// AAM 17.01.2003
// Tested

void
LinearLineElement::OutputNodeDataToVTK( const char* file_name,
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
           cerr <<"\nLinearLineElement::OutputNodeDataToVTK ";
           cerr <<"Output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): LinearLineElement: Variable: "<< var_name << endl;
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
           for ( uint32_t i{0U}; i<npe; i++ ) ofs << DATA(0,i) <<" ";
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
     cout <<"\nLinearLineElement::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK



} // end namespace csmp



