#include "Triangulator.h"
#include "Box.h"
#include "VSet.h"
#include "Matrix.h"
#include "Exception.h"
#include "LinearTriangle.h"
#include "TextInterface.h"
#include "convertColorToPermeability.h"

using namespace std;

namespace csmp {

#define DOWN_DIAG    0 	   // default square splitting |\|   
#define UP_DIAG      1     // square splitting |/|              

Triangulator::Triangulator()
 {
    m_mtrx = n_mtrx = 0;
 } 

Triangulator::~Triangulator()
 {
 }


/**
 
The permeability is output into el_perm for both triangles.
Method expects m and n as matrix ij indexes ranging from 1...m and 1...n.
The method returns whether a pixed should be split into upward diagonal
or downward diagonal triangles.  
*/
short Triangulator::TestOutline( unsigned int m, unsigned int n, 
                                 const Matrix& perm, 
                                 double& el_perm )
 {
     short c = DOWN_DIAG; // default 
     // clock-wise, inner and outer square numbering 
     double k1, k2, k3, k4, k5, k6, k7;
     
     // Reading the required permeabilities/hydraul. conductivities 
     // security for array indices avoiding to go beyond boundaries 
     //
     //    k6  k5 
     //    k7  k1  k4 
     //        k2  k3 
     //
     k1 = k2 = k3 = k4 = numeric_limits<double>::infinity();
     if ( m  > 0 && n  > 0 )     k1  = perm( (m-1), (n-1) );
     if ( m<m_mtrx && n  > 0 )   k2  = perm( (m),   (n-1) ); 
     if ( n<n_mtrx && m<m_mtrx ) k3  = perm( (m),   (n) );
     if ( n<n_mtrx && m  > 0 )   k4  = perm( (m-1), (n) );
     if ( m  > 1 && n  > 0 )     k5  = perm( (m-2), (n-1) );
     if ( n  > 1 && m  > 1 )     k6  = perm( (m-2), (n-2) );
     if ( n  > 1 && m  > 0 )     k7  = perm( (m-1), (n-2) );

     // 1. Default case for: if four k's are the same
     // The permeability of the two elements is always set to that of the first 
     // pixel that is considered
     el_perm = k1;
     
     // 2. The pixel square is split into up-diagonal triangles only
     //    if there is an updiagonal permeability boundary 
     if ( k2 == k4 ) if ( k1 != k2 || k3 != k2 ) c = UP_DIAG;
     if ( m>1 && n>1 && k5 == k7 ) if ( k1 != k5 || k5 != k6 ) c = UP_DIAG;
       
    return c;  

 } // end TestOutline



/**
 
Vertices are numbered 0..n-1 in accordance with the C/C++ array indexing
conventions.  
*/
unsigned long Triangulator::MapVertex( unsigned int m, unsigned int n )
 {
    unsigned long vert;
    
    // for the first and second line of data  
    if ( m == 1 ) {
         if ( n == 1 ) vert = 1; 
         else          vert = n * 2;
      }
    else if ( m == 2 ) {
       if ( n <= 2 ) vert = n + 1;  // first two elements in line 
       if ( n > 2 )  vert = 2 * n - 1;
     }
    // for the third line to the rest 
    else if ( m >= 3 ) vert = (m-1) * (n_mtrx+1) + n;

    // safety 
    if ( m >m_mtrx+1 || m < 1 ) 
      throw csmp::Exception( FATAL_ERROR, "Triangulator::MapVertex", 
                                   "n-array index out of range");

    if ( n > n_mtrx+1 || n < 1 ) 
      throw csmp::Exception( FATAL_ERROR, "Triangulator::MapVertex", 
                                   "n-array index out of range");
    return vert - 1U;
     
 } // end MapVertex



/**

A pixel bitmap image is used to generate triangular elements to which a 
permeability is assigned. The permeability is chosen to be that of the 
upper left pixel of four pixels that are converted into 4 nodes of two 
triangular elements generated from the square. 
 
By default, squares are broken as |\|. However, if the permeability of the 
upward diagonal pixels is the same, and the permeability of the down-diagonal 
pixels (left to right) differs, the square is split into two
upward-diagonal triangles. 

In the course of the mesh generation, nodes which lie at the
model boundary are flagged as such. In doing this, a coordinate system
is adhered to in which Y points downward while X points to the right. 
This schemen is in accordance with the pixel data that is used as input 
as the first row of the pixel matrix represents the top of the image. 

@param grid A Matrix with 256-integer values representing the color
codes for the permeability of the triangular elements which will be 
generated is the first method argument.
@param vset a VSet, a data container in which the finite-element data are stored as px, py, pz,
plist etc. arrays. 
 */
void Triangulator::TrianglesFromRegularGrid( const Matrix& grid, VSet<2U>& vset )
 {
   unsigned int i, j;
   short        split, t;
   int64_t        fed1, fed2, fed3;
   const short  bleft   = LEFT_OUTSIDE,
                bright  = RIGHT_OUTSIDE,
                btop    = TOP_OUTSIDE,
                bbottom = BOTTOM_OUTSIDE;
   
   // calculate overall parameters
   m_mtrx     = grid.Rows();
   n_mtrx     = grid.Cols();
   n_elements = m_mtrx * n_mtrx * 2;
   n_vertices = (m_mtrx+1) * (n_mtrx+1);
                                   
   vset.Resize( LinearTriangle().Nodes(),
                LinearTriangle().Neighbors(),
                LinearTriangle().ElementType(), 
                n_vertices, n_elements );
                                   
   deque<vector<int64_t> >     plist( n_elements, vector<int64_t>(3) );
   deque<vector<int64_t> >     pfvert( n_elements, vector<int64_t>(3) );
   // FEM_Data<ScalarVariable >  edata( ELEMENT, n_elements );
   PropertyData               edata( ELEMENT, SCALAR, 2U );
   edata.Resize( n_elements );
   
   // n counts the elements
   size_t  n(0U);
   
   for( i=1; i<=m_mtrx; i++ )
     for( j=1; j<=n_mtrx; j++ )
       {
         // 1. establishing how perm boundaries can be smoothed 
         split = TestOutline( i, j, grid, elperm );
         
         // 2. two cases for up- and down-diagonal elements 
         // FIRST CASE (DEFAULT): DOWN-DIAGONAL ELEMENTS 
         if ( split == DOWN_DIAG ) 
           { 		
     	     // 2.1.1 odd elements (0,1,3,5,n): 
     	     // -------------------------------
     	       // assigning the permeability  
             edata.Value(n) = elperm;
             // assigning vertex numbers
             plist[n][0] = MapVertex( i, j );   
             plist[n][1] = MapVertex( (i+1), j ); 
             plist[n][2] = MapVertex( (i+1), (j+1) ); 
 
             // 2.1.3 the vertex coordinates are taken 
             vset.Py( plist[n][0], i );      	
             vset.Px( plist[n][0], j );		    
             vset.Py( plist[n][1], i + 1 );
             vset.Px( plist[n][1], j );
             vset.Py( plist[n][2], i + 1 );
             vset.Px( plist[n][2], j + 1 );
         
             // 2.1.5 assigning boundary flag (node) values 
             if ( j == 1 )      vset.AddBFlag( plist[n][0], bleft );
             if ( j == 1 )      vset.AddBFlag( plist[n][1], bleft );
             if ( i == m_mtrx ) vset.AddBFlag( plist[n][1], bbottom );
             if ( i == m_mtrx ) vset.AddBFlag( plist[n][2], bbottom );
             if ( j == n_mtrx ) vset.AddBFlag( plist[n][2], bright );
 
             // 2.1.6 getting face-edge values, setting 
             // face edges are the elements that sit on the opposite site of a certain
             // node of the triangle 
             fed1 = ((n+1) + 2 * n_mtrx + 1);
             if ( i<m_mtrx and (t=TestOutline( (i+1), j, grid, elperm2 )) != DOWN_DIAG )
               fed1 = ((n+1) + 2 * n_mtrx);
             fed2 = ((n+1) + 1);
             if( j == 1 ) fed3 = LEFT_OUTSIDE;                      
             else fed3 = (n);
             if( fed1 > n_elements ) fed1 = BOTTOM_OUTSIDE;
             if ( fed1 > 0 ) fed1--;
             if ( fed2 > 0 ) fed2--;
             if ( fed3 > 0 ) fed3--;
             pfvert[n][0] = fed1;
             pfvert[n][1] = fed2;
             pfvert[n][2] = fed3;
             n++;  // one more element
  
             // 2.2.1 even elements 
             // ----------------------         
    	       // assigning the permeability
             edata.Value(n) = elperm;
	           plist[n][0] = MapVertex( (i+1), (j+1) ); // assign vertex numbers
	           plist[n][1] = MapVertex( i, (j+1) );
	           plist[n][2] = MapVertex( i, j ); 

             // 2.2.3 the vertex coordinates are taken 
             vset.Py( plist[n][0], i + 1 ); // node 1   
             vset.Px( plist[n][0], j + 1 );
             vset.Py( plist[n][1], i );     // node 2
             vset.Px( plist[n][1], j + 1 );
             vset.Py( plist[n][2], i );     // node 3
             vset.Px( plist[n][2], j );
          
             // 2.2.5 assigning boundary flags 
             if ( i == 1 )      vset.AddBFlag( plist[n][2], btop );
             if ( i == 1 )      vset.AddBFlag( plist[n][1], btop );
             if ( j == 1 )      vset.AddBFlag( plist[n][2], bleft );
             if ( i == m_mtrx ) vset.AddBFlag( plist[n][0], bbottom );
             if ( j == n_mtrx ) vset.AddBFlag( plist[n][0], bright );
             if ( j == n_mtrx ) vset.AddBFlag( plist[n][1], bright );

             // 2.2.6 getting face-edge values, setting 
             if ( i == 1 ) fed1 = TOP_OUTSIDE;  // O.K.
             else  {
                  if ( (t=TestOutline( (i-1), j, grid, elperm2 )) != DOWN_DIAG )
                    fed1 = static_cast<int32_t>((n+1) - 2 * n_mtrx);       
                  else fed1 = static_cast<int32_t>((n+1) - 2 * n_mtrx - 1); 
               } 
             fed2 = n; // O.K.
             if( j >= n_mtrx ) fed3 = RIGHT_OUTSIDE; // O.K.
             else fed3 = ((n+1) + 1);  // O.K.
             if ( fed1 > 0 ) fed1--;
             if ( fed2 > 0 ) fed2--;
             if ( fed3 > 0 ) fed3--;
             pfvert[n][0] = fed1;
             pfvert[n][1] = fed2;
             pfvert[n][2] = fed3;
             n++;  // one more element

             // end even 
         
          } // end down diagonal element case 
         
         
         // 3. SECOND CASE: UP-DIAGONAL ELEMENTS 
         //  -----
         //  | / |
         //  -----
         if ( split == UP_DIAG  ) 
          { 		  
             // 3.1 odd elements 
             // -------------------
    	       // 3.1.1 assigning the permeability
             edata.Value(n) = elperm;
             plist[n][0] = MapVertex( i, j );  // assign vertex numbers
             plist[n][1] = MapVertex( (i+1), j );
             plist[n][2] = MapVertex( i, (j+1) );

             // 3.1.3 the vertex coordinates are taken 
             vset.Py( plist[n][0], i );      	
             vset.Px( plist[n][0], j );		
             vset.Py( plist[n][1], i + 1 );
             vset.Px( plist[n][1], j );
             vset.Py( plist[n][2], i );
             vset.Px( plist[n][2], j + 1 );
         
             // 3.1.5 assigning boundary flags
             if ( i == 1 )      vset.AddBFlag( plist[n][0], btop );
             if ( i == 1 )      vset.AddBFlag( plist[n][2], btop );
             if ( j == 1 )      vset.AddBFlag( plist[n][0], bleft );
             if ( j == 1 )      vset.AddBFlag( plist[n][1], bleft );
             if ( i == m_mtrx ) vset.AddBFlag( plist[n][1], bbottom );
             if ( j == n_mtrx ) vset.AddBFlag( plist[n][2], bright );

             // 3.1.6 getting face-edge values, setting 
             fed1 = static_cast<int32_t>((n+1) + 1);  // O.K.
             if ( i == 1 ) fed2 = TOP_OUTSIDE;      // O.K.
             else {
                  if ( (t=TestOutline( (i-1), j, grid, elperm2 )) != DOWN_DIAG )
                    fed2 = static_cast<int32_t>((n+1) - 2 * n_mtrx + 1);   
		              else fed2 = static_cast<int32_t>((n+1) - 2 * n_mtrx);  
		           }
             if( j == 1 ) fed3 = LEFT_OUTSIDE; // O.K.
             else fed3 = static_cast<int32_t>(n); // O.K.
             if ( fed1 > 0 ) fed1--;
             if ( fed2 > 0 ) fed2--;
             if ( fed3 > 0 ) fed3--;
             pfvert[n][0] = fed1;
             pfvert[n][1] = fed2;
             pfvert[n][2] = fed3;
             n++;  // one more element
  
             // 3.2.1 even elements 
             // ----------------------
    	       // assigning the permeability
             edata.Value(n) = elperm;
             plist[n][0] = MapVertex( (i+1), j );  // assign vertex numbers
             plist[n][1] = MapVertex( (i+1), (j+1) );
             plist[n][2] = MapVertex( i, (j+1) );
        
             // 3.2.3 the vertex coordinates are taken 
             vset.Py( plist[n][0], i + 1 );    
             vset.Px( plist[n][0], j );
             vset.Py( plist[n][1], i + 1 );
             vset.Px( plist[n][1], j + 1 );
             vset.Py( plist[n][2], i );
             vset.Px( plist[n][2], j + 1 );

             // 3.2.5 reading the boundary values as we go along 
             if ( i == 1 )      vset.AddBFlag( plist[n][2], btop );
             if ( j == 1 )      vset.AddBFlag( plist[n][0], bleft );
             if ( i == m_mtrx ) vset.AddBFlag( plist[n][0], bbottom );
             if ( i == m_mtrx ) vset.AddBFlag( plist[n][1], bbottom );
             if ( j == n_mtrx ) vset.AddBFlag( plist[n][1], bright );
             if ( j == n_mtrx ) vset.AddBFlag( plist[n][2], bright );

             // 3.2.6 getting face-edge values, setting 
             if ( j >= n_mtrx ) fed1 = RIGHT_OUTSIDE;  
             else fed1 = ((n+1) + 1);
             fed2 = n;
             if ( i >= m_mtrx ) fed3 = BOTTOM_OUTSIDE; 
             else  {
                  if ( i<m_mtrx and (t=TestOutline( (i+1), j, grid, elperm2 )) != DOWN_DIAG )
                    fed3 = ((n+1) + 2 * n_mtrx - 1);
		              else fed3 = ((n+1) + 2 * n_mtrx);     
    		       }
             if ( fed1 > 0 ) fed1--;
             if ( fed2 > 0 ) fed2--;
             if ( fed3 > 0 ) fed3--;
             pfvert[n][0] = fed1;
             pfvert[n][1] = fed2;
             pfvert[n][2] = fed3;
             n++;  // one more element

             // end even 
         
	     } // end up-diagonal case 
	   
     } // end i-loop 
   
   // flagging the 4 corner elements, overwriting the side flags
   vset.AddBFlag( m_mtrx * (n_mtrx+1), CNR_MIN );   // CNR1 
   vset.AddBFlag( n_vertices-1U,   CNR_MIN_MAXX );  // CNR2 
   vset.AddBFlag( 2 * n_mtrx + 1U, CNR_MAX_MAXX );  // CNR3 
   vset.AddBFlag( 0U,              CNR_MAX_MINXZ ); // CNR4 
   
   // setting up the vset
   vset.AddPlist( plist.begin(), plist.end() );
   vset.AddPfverts( pfvert.begin(), pfvert.end() );
   vset.AddData( "permeability", edata );
   
   // Flipping the Y-axis
   double ymax = vset.Py(0);
   for ( i=1; i<vset.Vertices(); i++ ) if ( vset.Py(i) > ymax ) ymax = vset.Py(i);
   for ( i=0; i<vset.Vertices(); i++ ) {
        vset.Py( i, ymax - vset.Py(i) );
        vset.Pz( i, 0. );
     }

  } // end finite_element_mesh 




/**

Reads a 264-colour coded text-input file (matrix format) and uses Triangulator to generate a simple
FE mesh with uniform triangles that is returned as a VSet.

*/
VSet<2U> readTextPixelData()
{
    cout <<"\nreadTextPixelData: Please specify 256-color-coded geometry input file (text): ";
    string model_name;
    cin >> model_name;

    // 0.1 Reading a pixelated permeability image ASCII file into the
    //    new Matrix 'pixelcolors'
    size_t m, n;
    TextInterface().SizeofPixelTextImage256( model_name.c_str(), m, n );
    cout <<"\nreadTextPixelData: The size (in pixels) of the input image is: "<< n <<"h x "<< m;
    cout <<"v"<< endl;
    Matrix  pixelcolors(m,n);
    TextInterface().ReadPixelTextImage256( model_name.c_str(), pixelcolors );

    // 0.2 Converting the 256-color values into permeabilities
    convertColorToPermeability( 1., pixelcolors );

    // 0.3 Building a 2d mesh of triangular elements, incorporating the
    //    permeability data and the boundary conditions. The mesh is
    //    stored in the VSet object 'vset'
    VSet<2U>  vset;
    Triangulator().TrianglesFromRegularGrid( pixelcolors, vset );

    // 0.4 Scaling the geometrical input object that will become the Region
    //    The origin of the object is assumed to be zero.
    cout <<"\nreadTextPixelData: Please enter the horizontal and vertical dimensions of the model (m): ";
    double extent1, extent2;
    cin >> extent1 >> extent2;
    double zero(0.);
    vset.ScaleCoordinateToRange( 'x', zero, extent1 );
    vset.CoordinateRange( 'x', zero, extent1 );
    cout <<"\nAssigned X range: "<< zero <<" to "<< extent1 << " meter." << endl;
    vset.ScaleCoordinateToRange( 'y', zero, extent2 );
    vset.CoordinateRange( 'y', zero, extent2 );
    cout <<"\nAssigned Y range: "<< zero <<" to "<< extent2 << " meter." << endl;

    return vset;
} //


} // end namespace csmp


