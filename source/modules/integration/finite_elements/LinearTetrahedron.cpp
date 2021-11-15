#include "LinearTetrahedron.h"
#include "TriangularFacet.h"

using namespace std;

namespace csmp {

LinearTetrahedron::LinearTetrahedron()
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( LINEAR_TETRAHEDRON, false, false, 1U )
 {
   dim = 3;
   itp = 1;
   npf = 3;
   npe = 4;
   fpe = 4;
   spe = 6;
   epe = 4;
   nne = 24;
   cne = 0;
   gpe = 0;
   XY.Resize(npe,dim);
   M.Resize(npe,npe);
   UsesLocalCoordinates(false);
   Isoparametric(false);
   VolumeElement();
   ElementType(LINEAR_TETRAHEDRON);
 }


LinearTetrahedron::~LinearTetrahedron()
 {
 }



/** Outputs the shape function derivatives at the point 'xyz' which must lie
within the tetrahedron.
 */
void LinearTetrahedron::N( std::vector<double>& N, const std::vector<double>& xyz )
{

   const double x12( XY(0,0) - XY(1,0) );
   const double x13( XY(0,0) - XY(2,0) );
   const double x14( XY(0,0) - XY(3,0) );
   const double x23( XY(1,0) - XY(2,0) );
   const double x24( XY(1,0) - XY(3,0) );
   const double x34( XY(2,0) - XY(3,0) );
   const double x21( -x12 );
   const double x31( -x13 );
   const double x32( -x23 );
   const double x42( -x24 );
   const double x43( -x34 );

   const double y12( XY(0,1) - XY(1,1) );
   const double y13( XY(0,1) - XY(2,1) );
   const double y14( XY(0,1) - XY(3,1) );
   const double y23( XY(1,1) - XY(2,1) );
   const double y24( XY(1,1) - XY(3,1) );
   const double y34( XY(2,1) - XY(3,1) );
   const double y21( -y12 );
   const double y31( -y13 );
   const double y32( -y23 );
   const double y42( -y24 );
   const double y43( -y34 );

   const double z12( XY(0,2) - XY(1,2) );
   const double z13( XY(0,2) - XY(2,2) );
   const double z14( XY(0,2) - XY(3,2) );
   const double z23( XY(1,2) - XY(2,2) );
   const double z24( XY(1,2) - XY(3,2) );
   const double z34( XY(2,2) - XY(3,2) );
   const double z21( -z12 );
   const double z31( -z13 );
   const double z32( -z23 );
   const double z42( -z24 );
   const double z43( -z34 );

   const double a1( y42*z32 - y32*z42 );
   const double a2( y31*z43 - y34*z13 );
   const double a3( y24*z14 - y14*z24 );
   const double a4( y13*z21 - y12*z31 );

   const double b1( x32*z42 - x42*z32 );
   const double b2( x43*z31 - x13*z34 );
   const double b3( x14*z24 - x24*z14 );
   const double b4( x21*z13 - x31*z12 );

   const double c1( x42*y32 - x32*y42 );
   const double c2( x31*y43 - x34*y13 );
   const double c3( x24*y14 - x14*y24 );
   const double c4( x13*y21 - x12*y31 );

   const double V00( x21*( y23*z34 - y34*z23) + x32*(y34*z12 - y12*z34) + x43*(y12*z23 - y23*z12) );

   const double V01( XY(1,0) * ( XY(2,1)*XY(3,2) - XY(3,1)*XY(2,2) ) + XY(2,0) * ( XY(3,1)*XY(1,2) - XY(1,1)*XY(3,2) ) + XY(3,0) * ( XY(1,1)*XY(2,2) - XY(2,1)*XY(1,2) ) );
   const double V02( XY(0,0) * ( XY(3,1)*XY(2,2) - XY(2,1)*XY(3,2) ) + XY(2,0) * ( XY(0,1)*XY(3,2) - XY(3,1)*XY(0,2) ) + XY(3,0) * ( XY(2,1)*XY(0,2) - XY(0,1)*XY(2,2) ) );
   const double V03( XY(0,0) * ( XY(1,1)*XY(3,2) - XY(3,1)*XY(1,2) ) + XY(1,0) * ( XY(3,1)*XY(0,2) - XY(0,1)*XY(3,2) ) + XY(3,0) * ( XY(0,1)*XY(1,2) - XY(1,1)*XY(0,2) ) );
   const double V04( XY(0,0) * ( XY(2,1)*XY(1,2) - XY(1,1)*XY(2,2) ) + XY(1,0) * ( XY(0,1)*XY(2,2) - XY(2,1)*XY(0,2) ) + XY(2,0) * ( XY(1,1)*XY(0,2) - XY(0,1)*XY(1,2) ) );

   N[0] = ( V01 + a1*xyz[0] + b1*xyz[1] + c1*xyz[2])/V00;
   N[1] = ( V02 + a2*xyz[0] + b2*xyz[1] + c2*xyz[2])/V00;
   N[2] = ( V03 + a3*xyz[0] + b3*xyz[1] + c3*xyz[2])/V00;
   N[3] = ( V04 + a4*xyz[0] + b4*xyz[1] + c4*xyz[2])/V00;

}







/** Outputs shape function derivative matrix of the form:


    dN1dx dN2dx dN3dx dN4dx
B = dN1dy dN2dy dN3dy dN4dy
    dN1dz dN2dz dN3dz dN4dz


@section arguments Input Arguments

The parent Element and the matrix to hold the shape function derivatives.

@param B The shape function derivative matrix with the dimensions: 3 x 4.

*/
void LinearTetrahedron::dN( DenseMatrix<DM_MIN>& B )
{
   UpdateFor();
   B.Resize(dim,npe);

   // volume
   double vol6 = M(0,0)+M(0,1)+M(0,2)+M(0,3);

   for ( size_t i=0; i<npe; i++ )
    {
       B(0,i) = M(1,i) / vol6;
       B(1,i) = M(2,i) / vol6;
       B(2,i) = M(3,i) / vol6;
    }

//   cout <<"\nLinearTetrahedron::dN: volume: "<< vol6/6.0 << endl;
//   nicePrint( B );
}





/** Outputs the volume of the current tetrahedron.

@section arguments Input Arguments

A reference to the parent element.

@return The volume of the tetrahedron.

*/
double  LinearTetrahedron::Volume()
{
   UpdateFor();

   // volume
//   return std::fabs( (M(0,0)+M(0,1)+M(0,2)+M(0,3)) / 6.0 );
   return (M(0,0)+M(0,1)+M(0,2)+M(0,3)) / 6.;
}



void LinearTetrahedron::IntegralNN( DenseMatrix<DM_MIN>& M )
 {
     double vol = Volume();

     // consistent formulation: see users guide chapter 7
     // off-diagonal elements
     M.Resize(npe,npe);
     M(0,1) = M(0,2) = M(0,3) = M(1,2) = M(1,3) = M(2,3) = 0.05 * vol; // 1/20
     M(1,0) = M(2,0) = M(2,1) = M(3,0) = M(3,1) = M(3,2) = 0.05 * vol; // 1/20

     // diagonal (all elements must be multiplied by volume)
     M(0,0) = M(1,1) = M(2,2) = M(3,3) = 0.1 * vol; // 1/10
 }




/** Returns the counter-clockwise local node numbering for the element.

@param ids Returns the counter-clockwise local node numbering for the element.

*/
void LinearTetrahedron::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
 }


void LinearTetrahedron::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
 }



void
LinearTetrahedron::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
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
    std::cerr <<"\nLinearTetrahedron::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment




/**
    The faces are numbered such that face 0 lies opposite of
    node 0, face 1 node 1 etc., see CSMP FEM_conventions.pdf.
    
    @note the nodes are ordered counter-clockwise from the outside looking in.
*/
void LinearTetrahedron::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    fnids.resize(3);

    if ( face_id == 0 )
      {
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 3;
         return;
      }
    if ( face_id == 1 )
      {
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 2;
         return;
      }
    if ( face_id == 2 )
      {
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 3;
         return;
      }
    if ( face_id == 3 )
      {
         fnids[0] = 0;
         fnids[1] = 2;
         fnids[2] = 1;
         return;
      }
    else
    std::cerr <<"\nLinearTetrahedron::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }



vector<size_t>  LinearTetrahedron::CornerNodesOfFace( size_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<size_t>{1,2,3};
        case 1: return vector<size_t>{0,3,2};
        case 2: return vector<size_t>{0,1,3};
        case 3: return vector<size_t>{0,2,1};
      }
    cerr <<"\nLinearTetrahedron::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<size_t>{};
 }



CSMP_FEM_TYPE LinearTetrahedron::ElementTypeOfFace( size_t ) const
 {
    return LINEAR_TRIANGLE3D;
 }



/// returns rolling index that rolls over when i+a>=4
size_t LinearTetrahedron::n( size_t i, size_t a )
  {
     if ( i+a >= 4U ) return i + a - 4U;
     return i+a;
  }


/** Updates the Coordinate matrix XY and the testfunction coefficient
matrix M if the element ID has changed.
*/
void LinearTetrahedron::UpdateFor()
{
   static size_t tetrahedron_id(ULONG_MAX);

   if ( tetrahedron_id != CurrentID() )
     {
        // test function coefficients
        int32_t j(1);
        for ( size_t i=0; i<npe; i++ )
          {
             // a(i)
             M(0,i)  = -XY(n(i,1),0) * (XY(n(i,3),1)*XY(n(i,2),2)-XY(n(i,3),2)*XY(n(i,2),1));
             M(0,i) -=  XY(n(i,2),0) * (XY(n(i,1),1)*XY(n(i,3),2)-XY(n(i,1),2)*XY(n(i,3),1));
             M(0,i) -=  XY(n(i,3),0) * (XY(n(i,2),1)*XY(n(i,1),2)-XY(n(i,2),2)*XY(n(i,1),1));
             // b(i)
             M(1,i)  = XY(n(i,3),1)*XY(n(i,2),2) - XY(n(i,3),2)*XY(n(i,2),1);
             M(1,i) += XY(n(i,1),1)*XY(n(i,3),2) - XY(n(i,1),2)*XY(n(i,3),1);
             M(1,i) += XY(n(i,2),1)*XY(n(i,1),2) - XY(n(i,2),2)*XY(n(i,1),1);
             // c(i)
             M(2,i)  = XY(n(i,3),2)*XY(n(i,2),0) - XY(n(i,3),0)*XY(n(i,2),2);
             M(2,i) += XY(n(i,1),2)*XY(n(i,3),0) - XY(n(i,1),0)*XY(n(i,3),2);
             M(2,i) += XY(n(i,2),2)*XY(n(i,1),0) - XY(n(i,2),0)*XY(n(i,1),2);
             // d(i)
             M(3,i)  = XY(n(i,3),0)*XY(n(i,2),1) - XY(n(i,3),1)*XY(n(i,2),0);
             M(3,i) += XY(n(i,1),0)*XY(n(i,3),1) - XY(n(i,1),1)*XY(n(i,3),0);
             M(3,i) += XY(n(i,2),0)*XY(n(i,1),1) - XY(n(i,2),1)*XY(n(i,1),0);
             //
             M(0,i) *=  j;
             M(1,i) *=  j;
             M(2,i) *=  j;
             M(3,i) *=  j;
             j      *= -1;
          }
        tetrahedron_id = CurrentID();
     }

} // end UpdateFor


/**
     uses the TriangularFacet to compute the normals to its faces.

     @author SKM 15/2/2016
     
     @test OK 
*/
void  LinearTetrahedron::UnitNormalToFace( size_t face, std::vector<double>& unrml ) const
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



} // end namespace csmp



