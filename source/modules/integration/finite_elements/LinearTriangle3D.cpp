#include "LinearTriangle3D.h"
#include "Point.h"
#include <fstream>

using namespace std;

namespace csmp {

LinearTriangle3D::LinearTriangle3D()
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( LINEAR_TRIANGLE3D, false, false, 1U )
 {
   dim = 3;
   itp = 1;
   npf = 2;
   npe = 3;
   fpe = 3;
   spe = 3;
   epe = 3;
   nne = 6;
   cne = 0;
   gpe = 0;
   XY.Resize(npe,dim);
   M.Resize(npe,npe);
   
   UsesLocalCoordinates(false);
   Isoparametric( false );
   SurfaceElement();
   ElementType(LINEAR_TRIANGLE3D);
 }



LinearTriangle3D::~LinearTriangle3D() 
 {
 }




/// segments are numbered like faces
void 
LinearTriangle3D::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
 {
    snids.resize(2);
    if ( segm_id == 0 ) {
         snids[0] = 1;
         snids[1] = 2;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 2;
         snids[1] = 0;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 0;
         snids[1] = 1;
      }
 }



/// For this element, the faces are numbered such that face 0 lies opposite of
/// node 0, face 1 node 1 etc.
/*
void LinearTriangle3D::NodesOfFace( uint32_t face_id, std::vector<uint32_t>& fnids ) const
 {
    fnids.resize(2);

    if      ( face_id == 0 )
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
    std::cerr <<"\nLinearTriangle3D::NodesOfFace: Erratic face ID: "<< face_id << std::endl;
 }
*/




vector<uint32_t>  LinearTriangle3D::NodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nLinearTriangle3D::NodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }



vector<uint32_t>  LinearTriangle3D::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nLinearTriangle3D::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }





/// returns the local  numbers of the nodes at the other end of the sgment that the argument node is on
std::vector<uint32_t>  LinearTriangle3D::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nLinearTriangle3D::NodesConnectedTo: node "<< node_id <<" does not exist.";
    return vector<uint32_t>{};
  }





/// The linear triangle is numbered counter-clockwise by default.
/*
void LinearTriangle3D::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }
*/


/// clearly all nodes are corner nodes
void LinearTriangle3D::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }




CSMP_FEM_TYPE LinearTriangle3D::ElementTypeOfFace( uint32_t ) const
 {
    return LINEAR_BAR;
 }


/**
 
The aspect ratio of the triangle is computed as the ratio of the lengths
of the longest and shortest boundary segments of the triangle.  

@section arguments Input Arguments 

The node coordinates are retrieved from the parent element.  

@return The method returns the computed aspect ratio as a double value.  

@section application Application

The aspect ratio is a useful measure for determining whether an element
is too skewed in shape for the interpolation functions to give realistic
values of the interpolated property. It is therefore used to assess the
quality of a finiteelement mesh.  
 */
double  LinearTriangle3D::AspectRatio()
{
   vector<double>  vec(spe);

   EdgeLengths( vec );

   // order segment
   set<double> segms;

   for ( uint32_t i{0U}; i<spe; i++ ) segms.insert( vec[i] );

   double segm1 = (*segms.begin()), 
             segm2 = (*segms.rbegin());
   
   return segm2 / segm1;
}




/**
 
Method computes the radius of a circle which just fits into the inside
of the triangular element.  

@section arguments Input Arguments 

The method retrieves the node coordinates from the parent element.  

@return The method returns the radius of the inscribed circle of the triangle.  

@section implementation Implementation

The method computes the inner circle using the lengths of the segments 
bounding the triangle.  

@section application Application

The inner circle is a useful measure for defining the coarseness of grids
onto which finite element data are mapped and vice versa.  
*/
double  LinearTriangle3D::InnerRadius()
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
 
The triangles planar area in 3D is computed as the length of the 
cross-product of the vectors which define the edges 12 and 13 of the 
triangle.  

@section arguments Input Arguments 

The parent 'Element' object is queried for the node coordinates.  

@return The method returns the computed area.  
 */
double LinearTriangle3D::Volume()
{
  double X12 = XY(1,0) - XY(0,0), // X
           X31 = XY(0,0) - XY(2,0),
           Y12 = XY(1,1) - XY(0,1), // Y
           Y31 = XY(0,1) - XY(2,1),
           Z12 = XY(1,2) - XY(0,2), // Z
           Z31 = XY(0,2) - XY(2,2);

  // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
  double XNRM = -Y12*Z31 + Z12*Y31,
           YNRM = -Z12*X31 + X12*Z31,
           ZNRM = -X12*Y31 + Y12*X31;

  return 0.5 * sqrt( XNRM*XNRM + YNRM*YNRM + ZNRM*ZNRM );
}





void LinearTriangle3D::UnitNormal( std::vector<double>& vc ) const
 {
    // triangle is defined by two vectors a and b
    double a1 = XY(1,0) - XY(0,0), 
			  a2 = XY(1,1) - XY(0,1),
			  a3 = XY(1,2) - XY(0,2),
			  b1 = XY(2,0) - XY(0,0), 
			  b2 = XY(2,1) - XY(0,1),
			  b3 = XY(2,2) - XY(0,2);

  // normal to triangle (but not unit normal!)
    vc.resize(dim);
    vc[0]  = a2*b3 - a3*b2;
    vc[1]  = a3*b1 - a1*b3;
    vc[2]  = a1*b2 - a2*b1;
    // normalization to unit length and orienting normal CCW
    double length = -sqrt( vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2] );
    vc[0] /= length;
    vc[1] /= length;
    vc[2] /= length;

//    cout <<"\nLinearTriangle3D::UnitNormal: is: "<< endl;
//    out( vc );
 }







/** Computes the basis functions for this linear triangle element in 3D space
at the point 'xyz' within the plane of the triangle.  

@section arguments Input Arguments 

The method uses the data of the 'Element' to retrieve nodal coordinates,
and the point location supplied as third method argument.  

@param N The values of the 3 linear interpolation functions at the point 'xyz' are
returned into the second method argument which is a Meschpp double 
vector.  

@section implementation Implementation

The interpolation functions are calculated using an approach by Adrian
Umpleby (Imperial College) which is based on normals to planes through
faces perpendicular to the triangles plane.  
 */
void LinearTriangle3D::N( vector<double>& N, const vector<double>& xyz ) 
{
  const double X12 = XY(1,0) - XY(0,0), // X
         X23 = XY(2,0) - XY(1,0),
         X31 = XY(0,0) - XY(2,0),
         Y12 = XY(1,1) - XY(0,1), // Y
         Y23 = XY(2,1) - XY(1,1),
         Y31 = XY(0,1) - XY(2,1),
         Z12 = XY(1,2) - XY(0,2), // Z
         Z23 = XY(2,2) - XY(1,2),
         Z31 = XY(0,2) - XY(2,2);

  const double XP1 = XY(0,0) - xyz[0],
         XP2 = XY(1,0) - xyz[0],
         XP3 = XY(2,0) - xyz[0],
         YP1 = XY(0,1) - xyz[1],
         YP2 = XY(1,1) - xyz[1],
         YP3 = XY(2,1) - xyz[1],
         ZP1 = XY(0,2) - xyz[2],
         ZP2 = XY(1,2) - xyz[2],
         ZP3 = XY(2,2) - xyz[2];

  // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
  const double XNRM = -Y12*Z31 + Z12*Y31,
         YNRM = -Z12*X31 + X12*Z31,
         ZNRM = -X12*Y31 + Y12*X31,

         // XN12,YN12,ZN12 is normal to plane which goes through points 1 & 2 (not unit normal!)
         XN12 = Y12*ZNRM - Z12*YNRM,
         YN12 = Z12*XNRM - X12*ZNRM,
         ZN12 = X12*YNRM - Y12*XNRM,

         XN23 = Y23*ZNRM - Z23*YNRM,
         YN23 = Z23*XNRM - X23*ZNRM,
         ZN23 = X23*YNRM - Y23*XNRM,

         XN31 = Y31*ZNRM - Z31*YNRM,
         YN31 = Z31*XNRM - X31*ZNRM,
         ZN31 = X31*YNRM - Y31*XNRM,

         A123 = XN23*X12 + YN23*Y12 + ZN23*Z12,
         A231 = XN31*X23 + YN31*Y23 + ZN31*Z23,
         A312 = XN12*X31 + YN12*Y31 + ZN12*Z31,

         B23  = XN23*XP2 + YN23*YP2 + ZN23*ZP2,
         B31  = XN31*XP3 + YN31*YP3 + ZN31*ZP3,
         B12  = XN12*XP1 + YN12*YP1 + ZN12*ZP1;

   // Finds the values of three basis functions  at point 'xyz' given
   // points of triangle in 3d, X(1-3),Y(1-3),Z(1-3).
   N.resize(npe);
   N[0] = B23 / A123;
   N[1] = B31 / A231;
   N[2] = B12 / A312;
     
} // end N




void LinearTriangle3D::N_AtBaryCenter( std::vector<double>& IPOL )
 {
    vector<double>  xyz(3);
    xyz[0] = (XY(0,0) + XY(1,0) + XY(2,0)) / 3.;
    xyz[1] = (XY(0,1) + XY(1,1) + XY(2,1)) / 3.;
    xyz[2] = (XY(0,2) + XY(1,2) + XY(2,2)) / 3.;
    
    N( IPOL, xyz ); 
    
 } // N_AtBaryCenter



/**
 
Reports the derivatives of the 3 linear element interpolation functions.
These are stored into a matrix with m=rows=dimensions and n=cols=nodes
of the triangle.  

@section arguments Input Arguments 

The method uses the Element to get the nodal coordinates.  

@param B The partial derivatives of the element interpolation functions are
stored in the supplied matrix, where each node corresponds to a column
of x, y, and z partial derivatives.  

@section implementation Implementation

Method uses an unconventional approach designed by Adrian Umpleby at 
Imperial College. The derivatives are found using normals, n, to planes 
which are perpendicular to the triangle and go through the 2 nodes 
opposite of the node for which the derivative is computed. Thus the 
normal faces the node, lies in the plane of the triangle and is equivalent
in position to the height of the triangle. The formula for node 2 is:
 
dN1/dx = nx / n . -(p1 - p2)
dN1/dy = ny / n . -(p1 - p2)
dN1/dz = nz / n . -(p1 - p2)
 
Where . refers to the dot (scalar) product of the normal to plane through
nodes 2 and 3, and (p1-p2) is the vector giving by the difference of the
node points 12.  
*/
void LinearTriangle3D::dN( DenseMatrix<DM_MIN>& B )
{
  double X12 = XY(1,0) - XY(0,0), // X
           X23 = XY(2,0) - XY(1,0),
           X31 = XY(0,0) - XY(2,0),
           Y12 = XY(1,1) - XY(0,1), // Y
           Y23 = XY(2,1) - XY(1,1),
           Y31 = XY(0,1) - XY(2,1),
           Z12 = XY(1,2) - XY(0,2), // Z
           Z23 = XY(2,2) - XY(1,2),
           Z31 = XY(0,2) - XY(2,2);

  // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
  double XNRM = -Y12*Z31 + Z12*Y31,
           YNRM = -Z12*X31 + X12*Z31,
           ZNRM = -X12*Y31 + Y12*X31,

           // XN12 is normal to plane which goes through points 1 & 2 (not unit normal!)
           XN12 = Y12*ZNRM - Z12*YNRM,
           YN12 = Z12*XNRM - X12*ZNRM,
           ZN12 = X12*YNRM - Y12*XNRM,

           XN23 = Y23*ZNRM - Z23*YNRM,
           YN23 = Z23*XNRM - X23*ZNRM,
           ZN23 = X23*YNRM - Y23*XNRM,

           XN31 = Y31*ZNRM - Z31*YNRM,
           YN31 = Z31*XNRM - X31*ZNRM,
           ZN31 = X31*YNRM - Y31*XNRM;

   B.Resize(dim,npe);
   // Node 1: dN1/dx   dN1/dy   dN1/dz (plane 23)
   B(0,0) = XN23 / (XN23*-X12 + YN23*-Y12 + ZN23*-Z12);
   B(1,0) = YN23 / (XN23*-X12 + YN23*-Y12 + ZN23*-Z12);
   B(2,0) = ZN23 / (XN23*-X12 + YN23*-Y12 + ZN23*-Z12);

   // Node 2: dN2/dx   dN2/dy   dN2/dz (plane 31)
   B(0,1) = XN31 / (XN31*X12 + YN31*Y12 + ZN31*Z12);
   B(1,1) = YN31 / (XN31*X12 + YN31*Y12 + ZN31*Z12);
   B(2,1) = ZN31 / (XN31*X12 + YN31*Y12 + ZN31*Z12);

   // Node 3: dN3/dx   dN3/dy   dN3/dz (plane 12)
   B(0,2) = XN12 / (XN12*-X31 + YN12*-Y31 + ZN12*-Z31);
   B(1,2) = YN12 / (XN12*-X31 + YN12*-Y31 + ZN12*-Z31);
   B(2,2) = ZN12 / (XN12*-X31 + YN12*-Y31 + ZN12*-Z31);
   
} // end dN



/**
    Exactly equivalent to dN because the interpolation function derivatives are 
    constant inside the linear triangular element.
    
    Method returns the element volume.
*/
double  LinearTriangle3D::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    double X12 = XY(1,0) - XY(0,0), // X
             X23 = XY(2,0) - XY(1,0),
             X31 = XY(0,0) - XY(2,0),
             Y12 = XY(1,1) - XY(0,1), // Y
             Y23 = XY(2,1) - XY(1,1),
             Y31 = XY(0,1) - XY(2,1),
             Z12 = XY(1,2) - XY(0,2), // Z
             Z23 = XY(2,2) - XY(1,2),
             Z31 = XY(0,2) - XY(2,2);

    // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
    double XNRM = -Y12*Z31 + Z12*Y31,
             YNRM = -Z12*X31 + X12*Z31,
             ZNRM = -X12*Y31 + Y12*X31,

             // XN12 is normal to plane which goes through points 1 & 2 (not unit normal!)
             XN12 = Y12*ZNRM - Z12*YNRM,
             YN12 = Z12*XNRM - X12*ZNRM,
             ZN12 = X12*YNRM - Y12*XNRM,

             XN23 = Y23*ZNRM - Z23*YNRM,
             YN23 = Z23*XNRM - X23*ZNRM,
             ZN23 = X23*YNRM - Y23*XNRM,

             XN31 = Y31*ZNRM - Z31*YNRM,
             YN31 = Z31*XNRM - X31*ZNRM,
             ZN31 = X31*YNRM - Y31*XNRM;

     B.Resize(dim,npe);
     // Node 1: dN1/dx   dN1/dy   dN1/dz (plane 23)
     B(0,0) = XN23 / (XN23*-X12 + YN23*-Y12 + ZN23*-Z12);
     B(1,0) = YN23 / (XN23*-X12 + YN23*-Y12 + ZN23*-Z12);
     B(2,0) = ZN23 / (XN23*-X12 + YN23*-Y12 + ZN23*-Z12);

     // Node 2: dN2/dx   dN2/dy   dN2/dz (plane 31)
     B(0,1) = XN31 / (XN31*X12 + YN31*Y12 + ZN31*Z12);
     B(1,1) = YN31 / (XN31*X12 + YN31*Y12 + ZN31*Z12);
     B(2,1) = ZN31 / (XN31*X12 + YN31*Y12 + ZN31*Z12);

     // Node 3: dN3/dx   dN3/dy   dN3/dz (plane 12)
     B(0,2) = XN12 / (XN12*-X31 + YN12*-Y31 + ZN12*-Z31);
     B(1,2) = YN12 / (XN12*-X31 + YN12*-Y31 + ZN12*-Z31);
     B(2,2) = ZN12 / (XN12*-X31 + YN12*-Y31 + ZN12*-Z31);

     // volume calculation
     // XNRM,YNRM,ZNRM is normal to triangle (not unit normal!)
     return 0.5 * sqrt( XNRM*XNRM + YNRM*YNRM + ZNRM*ZNRM );
 }


/** 
    Note that the derivative is constant so that the point location is ignored.
    method returns the element volume.
*/
double  LinearTriangle3D::dN_At( DenseMatrix<DM_MIN>& B, const std::vector<double>& )
 {
    return dN_AtBarycenter( B );
 }
 










/**
 
Returns the lengths of the triangles segments into the supplied Meschpp
vector. The segments are defined as follows: Segment 1 is from node1 to
node2, segment 2 from node2 to 3, and segment 3 from node3 to 1.  

@param len vector into which the lengths will be stored.  

*/
void  LinearTriangle3D::EdgeLengths( vector<double>& len )
{
    double sum;
    len.resize(npe);

    // segment 1
    sum    = (XY(1,0)-XY(0,0)) * (XY(1,0)-XY(0,0));
    sum   += (XY(1,1)-XY(0,1)) * (XY(1,1)-XY(0,1));
    sum   += (XY(1,2)-XY(0,2)) * (XY(1,2)-XY(0,2));
    len[0] = sqrt(sum);

    // segment 2
    sum    = (XY(2,0)-XY(1,0)) * (XY(2,0)-XY(1,0));
    sum   += (XY(2,1)-XY(1,1)) * (XY(2,1)-XY(1,1));
    sum   += (XY(2,2)-XY(1,2)) * (XY(2,2)-XY(1,2));
    len[1] = sqrt(sum);

    // segment 3
    sum    = (XY(2,0)-XY(0,0)) * (XY(2,0)-XY(0,0));
    sum   += (XY(2,1)-XY(0,1)) * (XY(2,1)-XY(0,1));
    sum   += (XY(2,2)-XY(0,2)) * (XY(2,2)-XY(0,2));
    len[2] = sqrt(sum); 

 } // end EdgeLengths





/**
 
Computes the 'mass' matrix for the element for the case where the factor
which is used to multiply the matrix with is constant on the element and
can therefore be put in front of the integral.  

@section arguments Input Arguments 

The coordinates of the Element are retrieved in order to calculate
its area.  

@param CE The mass matrix (nodes x nodes) is returned into the second argument
of the method. This is a Meschach++ matrix which will be automatically
resized if it has not got the correct dimensions.  

@section implementation Implementation

The mass matrix = Integral Ni Nj dx dy dz is explained in many FE text 
books. For instance in 'Groundwater Modeling by the Finite Element Method"
AGU Monograph 13 by Jonathan Istok. The resulting matrix is output
in its consistent form.  

@section application Application

For example, to compute the element 'capacitance' matrix for transient
fluid flow.  

*/
void LinearTriangle3D::IntegralNN( DenseMatrix<DM_MIN>& CE )
{
   double area = Volume();
   
   CE.Resize(npe,npe);
   
   // diagonal elements: A/12 * 2
   CE(0,0) = CE(1,1) = CE(2,2) = area / 6.;
   
   // off-diagonal elements: A/12 * 1
   CE(0,1) = CE(0,2) = CE(1,0) = CE(1,2) = CE(2,0) = CE(2,1) = area / 12.;
}





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
*/
void LinearTriangle3D::OutputNodeDataToVTK( const char* file_name, 
                                            const char* var_name, 
                                            DenseMatrix<DM_MIN>& DATA ) const
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     snprintf( elmt, sizeof(elmt), "%lu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );
       
     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cerr <<"\nLinearTriangle3D::OutputNodeDataToVTK "; 
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
     for ( uint32_t i{0U}; i<npe; i++ )
       {
          for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;  
       
     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< npe+1 << endl;
     // number of points per cell, point1, ... point n
     ofs << npe <<" 0 1 2";
     ofs << endl;
     ofs << endl;
     
     // 4. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 5 << endl; // VTK_TRIANGLE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );
     
     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" float"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x6
           for ( uint32_t i{0U}; i<npe; i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          assert( DATA.Rows() == dim );
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 6
          for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
               for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nLinearTriangle3D::OutputNodeDataToVTK: file '"<< outfile;
     cout <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK




/**
 
Outputs triangle to VTK file for visualization. The scalar variable
called 'sum_N' sums the values of the interpolation functions at the 
nodes, and the vector variable 'derivative_N' reports the shape function 
derivatives at each node as vectors.

@section arguments Input Arguments 

The parent 'Element' and the name of the file to which the triangle
shall be output.  

@section application Application

Use this together with a tcl script newElementViewer.tcl if this is
available.  

 */
void LinearTriangle3D::OutputToVTK( const char* file_name )
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     snprintf( elmt, sizeof(elmt), "%lu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );
       
     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cerr <<"\nLinearTriangle3D::OutputToVTK "; 
           cerr <<"Output file could not be opened."<< endl;
           return;
       }  
       
     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): Element properties: "<< endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" float"<< endl;
     for ( uint32_t i{0U}; i<npe; i++ )
       {
          for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;  
       
     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< npe+1 << endl;
     // number of points per cell, point1, ... point n
     ofs << npe <<" 0 1 2";
     ofs << endl;
     ofs << endl;
     
     // 4. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 5 << endl; // VTK_TRIANGLE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );
     
     // sum of shape functions
     vector<double> IPOL(npe), xyz(dim);
     double  sum;
     
     ofs <<"SCALARS "<< "sum_N" <<" float"<< endl;
     ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
     for ( uint32_t i{0U}; i<npe; i++ ) {
          xyz[0] = XY(i,0); xyz[1] = XY(i,1); xyz[2] = XY(i,2); sum = 0.0;
          N( IPOL, xyz );
          for ( uint32_t j{0U}; j<npe; j++ ) sum += IPOL[j];
          ofs << sum <<" ";
       }
     ofs << endl;
     
     // shape function derivatives
     ofs <<"VECTORS "<< "derivative_N" <<" float"<< endl;

     DenseMatrix<DM_MIN> DN(dim,npe);
     dN( DN );
     
     for ( uint32_t i{0U}; i<DN.Cols(); i++ ) {
          for ( uint32_t j{0U}; j<DN.Rows(); j++ ) ofs << DN(j,i) <<"  ";
          ofs << endl;
       }
     ofs << endl;

     ofs.close();
     cout <<"\nLinearTriangle3D::OutputNodeDataToVTK: file '"<< outfile;
     cout <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK


/**
    outward pointing normals to the faces that lie on the opposite nodes
    
    @author SKM 18/2/2016
    
    @test OK - for 3D version
*/
void  LinearTriangle3D::UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const
 {
     assert( face < Faces() );
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





} // end namespace csmp

