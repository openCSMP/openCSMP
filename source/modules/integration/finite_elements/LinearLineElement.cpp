#include "LinearLineElement.h"

using namespace std;

namespace csmp {


LinearLineElement::LinearLineElement( size_t dimensions )
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


LinearLineElement::~LinearLineElement()
 {
 }



void  LinearLineElement::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
 }



void  LinearLineElement::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
 }


/// segments are numbered like faces
void LinearLineElement::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    if ( segm_id > 1 )
      std::cout <<"\nLinearLineElement::NodesOfSegment: There is only one segment present."<< std::endl;
    snids.resize(2);
    snids[0] = 0;
    snids[1] = 1;
 }


/// SKM fixed 10/6/2014
void  LinearLineElement::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    if ( face_id > 1U )
      std::cerr <<"\nLinearLineElement::NodesOfFace: There are only 2 faces present, corresponding to the nodes."<< std::endl;
    fnids.resize(1U);
    fnids[0] = face_id;
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
     double64  len;
     DN.Resize(dim,npe);

     // 2-dimensional edge
     if ( dim == 2 ) {
          org2.Set(XY(0,0),XY(0,1));
          dest2.Set(XY(1,0),XY(1,1));
          edge2.Set( org2, dest2 );
          len = edge2.Length();
          // x-derivatives
          DN(0,0) = -(dest2[0]-org2[0]) / len, DN(0,1) = (dest2[0]-org2[0]) / len;
          // y-derivatives
          DN(1,0) = -(dest2[1]-org2[1]) / len, DN(1,1) = (dest2[1]-org2[1]) / len;
       }
     // 3-dimensional edge
     else
       {
          org3.Set(XY(0,0),XY(0,1),XY(0,2));
          dest3.Set(XY(1,0),XY(1,1),XY(1,2));
          edge3.Set( org3, dest3 );
          len = edge3.Length();
          // x-derivatives
          DN(0,0) = -(dest3[0]-org3[0]) / len, DN(0,1) = (dest3[0]-org3[0]) / len;
          // y-derivatives
          DN(1,0) = -(dest3[1]-org3[1]) / len, DN(1,1) = (dest3[1]-org3[1]) / len;
          // z-derivatives
          DN(2,0) = -(dest3[2]-org3[2]) / len, DN(2,1) = (dest3[2]-org3[2]) / len;
       }

  } // end dN




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
as a double64 variable.
 */
double64 LinearLineElement::dN_At( DenseMatrix<DM_MIN>& DN, const vector<double64>& )
 {
     double64  len(std::numeric_limits<double>::quiet_NaN());
     DN.Resize(dim,npe);

     // 2-dimensional edge
     if ( dim == 2 ) {
          org2.Set(XY(0,0),XY(0,1));
          dest2.Set(XY(1,0),XY(1,1));
          edge2.Set( org2, dest2 );
          len = edge2.Length();
          // x-derivatives
          DN(0,0) = -(dest2[0]-org2[0]) / len, DN(0,1) = (dest2[0]-org2[0]) / len;
          // y-derivatives
          DN(1,0) = -(dest2[1]-org2[1]) / len, DN(1,1) = (dest2[1]-org2[1]) / len;
       }
     // 3-dimensional edge
     else
       {
          org3.Set(XY(0,0),XY(0,1),XY(0,2));
          dest3.Set(XY(1,0),XY(1,1),XY(1,2));
          edge3.Set( org3, dest3 );
          len = edge3.Length();
          // x-derivatives
          DN(0,0) = -(dest3[0]-org3[0]) / len, DN(0,1) = (dest3[0]-org3[0]) / len;
          // y-derivatives
          DN(1,0) = -(dest3[1]-org3[1]) / len, DN(1,1) = (dest3[1]-org3[1]) / len;
          // z-derivatives
          DN(1,0) = -(dest3[2]-org3[2]) / len, DN(1,1) = (dest3[2]-org3[2]) / len;
       }

     return len;

 } // end






/**

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The volume (m3) of the finite element.
*/
double64 LinearLineElement::Volume()
{
   // 2-dimensional models
   if ( dim == 2 ) {
         org2.Set(XY(0,0),XY(0,1));
         dest2.Set(XY(1,0),XY(1,1));
         return (dest2 - org2).Length();
     }

   // 3-dimensional models
   org3.Set(XY(0,0),XY(0,1),XY(0,2));
   dest3.Set(XY(1,0),XY(1,1),XY(1,2));
   return (dest3 - org3).Length();
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
void LinearLineElement::N( vector<double64>& FN, const vector<double64>& xy )
 {
     double64  len, l1, l2;

     // 2-dimensional edge
     if ( dim == 2 ) {
          org2.Set(XY(0,0),XY(0,1));
          dest2.Set(XY(1,0),XY(1,1));
          p2.Set(xy[0],xy[1]);
          len = (dest2 - org2).Length();
          l1  = (dest2 - p2).Length();
          l2  = len - l1;
       }
     // 3-dimensional edge
     else
       {
          org3.Set(XY(0,0),XY(0,1),XY(0,2));
          dest3.Set(XY(1,0),XY(1,1),XY(1,2));
          p3.Set(xy[0],xy[1],xy[2]);
          len = (dest3 - org3).Length();
          l1  = (dest3 - p3).Length();
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
double64 LinearLineElement::dN_AtNode( DenseMatrix<DM_MIN>& DN, size_t )
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



double64 LinearLineElement::dN_AtBarycenter( DenseMatrix<DM_MIN>& DN )
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
         cout <<"LinearLineElement::dN_AtNode: Method not implemented for 3D case. "<< endl;
         terminate();
      }

    return 0.0;
 }

void LinearLineElement::N_AtBaryCenter(std::vector<double64>& N)
{
	vector<double64> p;
	p.resize(dim);
	for (size_t i = 0; i < dim; ++i) p[i] = 0.5*(XY(0, i) + XY(1, i));
	LinearLineElement::N(N, p);
}



/**
    In 3D, the input vector must contain the rotation axis around which
    the edge is turned to create the unit normal.
*/
void  LinearLineElement::UnitNormal( vector<double64>& vc ) const
 {
    vc.resize(dim);
    static bool first_call(true);

    if ( dim == 2 ) {
        edge2.Set( mjl::Point(XY(0,0),XY(0,1)), mjl::Point(XY(1,0),XY(1,1)) );
        edge2.Rot();
        edge2.NormalizeTo(1.0);
        vc[0] = edge2.Destination()[0];
        vc[1] = edge2.Destination()[1];
      }
    else
      {
        if ( first_call )  cout <<"\nLinearLineElement::UnitNormal: Using input vector as reference axis."<< endl;
        edge3.Set( mjl::Point3D(XY(0,0),XY(0,1),XY(0,2)),
                   mjl::Point3D(XY(1,0),XY(1,1),XY(1,2)) );
        edge3.Rot( mjl::Point3D( vc[0], vc[1], vc[2] ) );
        edge3.Normalize();
        vc[0] = edge3.dest_.X();
        vc[1] = edge3.dest_.Y();
        vc[2] = edge3.dest_.Z();
        first_call = false;
      }
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
void  LinearLineElement::UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const
 {
     assert( face < Faces() );
   
     if ( Dim() == 1 ) {
         unrml.resize(1);
         if ( face == 0 ) {
              unrml[0] = -1.;
              return;
           }
         if ( face == 1 ) {
              unrml[0] = 1.;
              return;
           }
         return;
       }
   
     // the normal lies in the plane of the model
     // point in the direction of the element
     if ( Dim() == 2 ) {
         unrml.resize(2);
         if ( face == 0 ) {
              unrml[0] = XY(0,0) - XY(1,0);
              unrml[1] = XY(0,1) - XY(1,1);
              // normalise to length
              const double64 length = hypot(unrml[0], unrml[1]);
              unrml[0] /= length;
              unrml[1] /= length;
              return;
           }
         if ( face == 1 ) {
              unrml[0] = XY(1,0) - XY(0,0);
              unrml[1] = XY(1,1) - XY(0,1);
              const double64 length = hypot(unrml[0], unrml[1]);
              unrml[0] /= length;
              unrml[1] /= length;
              return;
           }
         return;
       }

     // 3D
     if ( Dim() == 3 ) {
         unrml.resize(3);
         if ( face == 0 ) {
              unrml[0] = XY(0,0) - XY(1,0);
              unrml[1] = XY(0,1) - XY(1,1);
              unrml[2] = XY(0,2) - XY(1,2);
              // normalise to length
              const double64 length = sqrt(unrml[0]*unrml[0] + unrml[1]*unrml[1] + unrml[2]*unrml[2]);
              unrml[0] /= length;
              unrml[1] /= length;
              unrml[2] /= length;
              return;
           }
         if ( face == 1 ) {
              unrml[0] = XY(1,0) - XY(0,0);
              unrml[1] = XY(1,1) - XY(0,1);
              unrml[2] = XY(1,2) - XY(0,2);
              const double64 length = sqrt(unrml[0]*unrml[0] + unrml[1]*unrml[1] + unrml[2]*unrml[2]);
              unrml[0] /= length;
              unrml[1] /= length;
              unrml[2] /= length;
              return;
           }
         return;
       }
   
 } // end UnitNormalToFace




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
     sprintf( elmt, "%lu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cout <<"\nLinearLineElement::OutputNodeDataToVTK ";
           cout <<"Output file could not be opened."<< endl;
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
     cout <<"\nLinearLineElement::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK



} // end namespace csmp



