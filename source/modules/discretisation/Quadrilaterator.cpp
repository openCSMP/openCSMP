#include "Quadrilaterator.h"
#include "Box.h"
#include "PropertyDatabase.h"
#include "convertColorToPermeability.h"
#include "TextInterface.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "VSet.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/**
 
Constructor. By setting the boolean to true (default is false), the bitmap 
entires corresponding to a nodal permeability value are harmonically averaged 
for each quadrilateral finite element. 
 */
Quadrilaterator::Quadrilaterator( bool harmonic_permeability_averaging )
 : harmonic(harmonic_permeability_averaging)
 {
   cols = rows = n_vertices = n_elements = 0;
 } 



Quadrilaterator::~Quadrilaterator()
 {
 }


/**
 
Calculates the harmonic permeability average of the four corner nodes of the 
quadrilateral. 
*/
double Quadrilaterator::HarmonicPermeabilityAverage( unsigned int m, unsigned int n, 
                                                       const Matrix& perm ) const
 {
     double kinv  = 1. / perm( m, n  );
     kinv += 1.0/perm(m+1,n  );
     kinv += 1.0/perm(m+1,n+1);
     kinv += 1.0/perm(m,  n+1);
     
     double k = 4. / kinv;
       
     return k;  

 }



/** Calculates the global node ID for the quadrilateral (0..n-1).
*/
size_t  Quadrilaterator::MapVertex( unsigned int m, unsigned int n ) const
 {
    size_t vert(UNSPECIFIED);
    
    // for the first and second line of data  
    if ( m == 1 ) {
         if ( n == 1 ) vert = 1;
         else vert = n * 2;
      }
    if ( m == 2 ) {
        if ( n <= 2 ) vert = n + 1;  // first two elements in line
        if ( n > 2 )  vert = 2 * n - 1;
     }
    // for the third line to the rest 

    if ( m >= 3 ) vert = m * cols - cols + n;

    // safety 
    if ( m > rows || m < 1 ) 
      throw csmp::Exception( FATAL_ERROR, "Quadrilaterator::MapVertex", 
                                   "n-array index out of range");

    if ( n > cols || n < 1 ) 
      throw csmp::Exception( FATAL_ERROR, "Quadrilaterator::MapVertex", 
                                   "n-array index out of range");
    return vert - 1U;
     
 } // end MapVertex



void Quadrilaterator::ReadPixelMatrix( bool from_bitmap )
 {
    size_t            m, n;
    TextInterface        IO;
    char                 name[200];
    
    cout << "\nQuadrilaterator< dim>::ReadPixelMatrix: ";
    cout <<" Enter the name of the 256 color bitmap input file: " << endl;
    cin >> name;
    
    IO.SizeofPixelTextImage256( name, m, n );
    matrix_.Resize(m,n);
    IO.ReadPixelTextImage256( name, matrix_ );
    if ( from_bitmap ) convertColorToPermeability( 1., matrix_ );
 }



void Quadrilaterator::ReadPixelMatrix( const char* name, bool from_bitmap )
 {
    size_t            m, n;
    TextInterface        IO;

    IO.SizeofPixelTextImage256( name, m, n );
    matrix_.Resize(m,n);
    IO.ReadPixelTextImage256( name, matrix_ );
    if ( from_bitmap ) convertColorToPermeability( 1., matrix_ );
 }




void Quadrilaterator::ScaleModelRange( VSet<2U>& vset_scaled ) const
 {
    double zero(0.0),  extent1,  extent2;
    // Scaling the geometrical input object that will become the Model
    // The origin of the object is assumed to be zero.
    cout <<"\nQuadrilaterator<"<< 2U <<">::ScaleModelRange: ";
    cout <<"\nPlease enter the horizontal (" << cols << " cols) and vertical (" << rows;
    cout << " rows) dimensions of the model (m): " << endl;
    cin >> extent1 >> extent2;
    vset_scaled.ScaleCoordinateToRange( 'x', zero, extent1 );
    vset_scaled.CoordinateRange( 'x', zero, extent1 );
    cout <<"\nAssigned X range: "<< zero <<" to "<< extent1 << " meter." << endl;
    vset_scaled.ScaleCoordinateToRange( 'y', zero, extent2 );
    vset_scaled.CoordinateRange( 'y', zero, extent2 );
    cout <<"\nAssigned Y range: "<< zero <<" to "<< extent2 << " meter." << endl;

 }
 
 
 

void Quadrilaterator::ScaleModelRange( VSet<2U>& vset_scaled, double x_dim, double y_dim ) const
 {
    double zero(0.0);
    // Scaling the geometrical input object that will become the Model
    // The origin of the object is assumed to be zero.
    cout <<"\nQuadrilaterator<"<< 2U <<">::ScaleModelRange: ";
    vset_scaled.ScaleCoordinateToRange( 'x', zero, x_dim );
    vset_scaled.CoordinateRange( 'x', zero, x_dim );
    cout <<"\nScaling the model horizontally (" << cols << " cols) from " << zero << " to " << x_dim << " meters ";
    vset_scaled.ScaleCoordinateToRange( 'y', zero, y_dim );
    vset_scaled.CoordinateRange( 'y', zero, y_dim );
    cout <<"\nScaling the model vertically (" << rows << " rows) from " << zero << " to " << y_dim << " meters ";
 }




/**
 
A pixel bitmap image is used to generate 2D straightsided isoparametric 
quadrilateral elements to which a permeability is assigned. 
The permeability is chosen to be that of the 
lower left pixel of four pixels that are converted into 4 nodes of the 
quadrilateral elements. Alternatively, by setting the boolean of the constructor
to true, the harmonic permeability average is calculated for each element from
the four nodes.  
 
The method automatically assigns the boundary flags for nodes and elements and
prompts the user to scale the model to the physical
dimensions. The finite element mesh is stored in a VSet object that can
be used to construct the Model.  

The nodes are ordered in a counter-clockwise fashion with the first node being 
the lower left corner of the quadrilateral.

@param vset The reference to the VSet that will be used to construct the Model
object later on.
*/
void Quadrilaterator::QuadrilateralsFromRegularGrid( VSet<2U>& vset, bool from_bitmap )
 {
   // read the color-coded permeability textfile
   ReadPixelMatrix( from_bitmap );
   
   cout << "\nQuadrilaterator< dim>::QuadrilateralsFromRegularGrid: \
           Generating Vset for finite element mesh... " << endl;
   
   // calculate overall parameters
   rows     = matrix_.Rows();
   cols     = matrix_.Cols();
   n_elements = (rows-1) * (cols-1);
   n_vertices = rows * cols;
   
   // initialize vset                               
   vset.Resize( IsoparametricLinearQuadrilateral().Nodes(),
                IsoparametricLinearQuadrilateral().Neighbors(),
                IsoparametricLinearQuadrilateral().ElementType(), 
                n_vertices, n_elements );
                 
   // generate vset                                
   GenerateVSet( vset );

   // scale the vset to the physical property range
   ScaleModelRange( vset );
   vset.EstablishZeroBasedNumbering();

  } // end finite_element_mesh 





void Quadrilaterator::QuadrilateralsFromRegularGrid( VSet<2U>& vset, const char* file_name, 
                                                     double x_extend, double y_extend, bool from_bitmap )
 {
   // read the colorcoded permeability textfile
   ReadPixelMatrix( file_name, from_bitmap );
   
   cout << "\nQuadrilaterator< dim>::QuadrilateralsFromRegularGrid: Generating Vset for finite element mesh... " << endl;
   
   // calculate overall parameters
   rows     = matrix_.Rows();
   cols     = matrix_.Cols();
   n_elements = (rows-1) * (cols-1);
   n_vertices = rows * cols;
                                   
   // initialize vset                               
   vset.Resize( IsoparametricLinearQuadrilateral().Nodes(),
                IsoparametricLinearQuadrilateral().Neighbors(),
                IsoparametricLinearQuadrilateral().ElementType(),
                n_vertices, n_elements );
                 
   // generate vset                                
   GenerateVSet( vset );

   //vset.Out();
   
   // scale the vset to the physical property range
   ScaleModelRange( vset, x_extend, y_extend );
   
   vset.EstablishZeroBasedNumbering();

  } // end finite_element_mesh 


/**
    @attention x_nodes actually refers to the number of rows in the matrix, i.e. nodes along the Y (vertical) axes,
    y_nodes are the number of nodes in the horizontal=x direction.
*/
void Quadrilaterator::QuadrilateralsFromRegularGrid( VSet<2U>& vset, double x_extend, double y_extend, size_t x_nodes, size_t y_nodes )
  {
    cout << "\nQuadrilaterator< dim>::QuadrilateralsFromRegularGrid: Generating Vset for finite element mesh... " << endl;

    matrix_.Resize(x_nodes,y_nodes); 

    // calculate overall parameters
    rows     = matrix_.Rows();
    cols     = matrix_.Cols();
    n_elements = (rows-1) * (cols-1);
    n_vertices = rows * cols;

    // initialize vset                               
    vset.Resize( IsoparametricLinearQuadrilateral().Nodes(),
                 IsoparametricLinearQuadrilateral().Neighbors(),
                 IsoparametricLinearQuadrilateral().ElementType(),
                 n_vertices, n_elements );

    // generate vset                                
    GenerateVSet( vset );

    // scale the vset to the physical property range
    ScaleModelRange( vset, x_extend, y_extend );
    cout <<"\nQuadrilaterator::QuadrilateralsFromRegularGrid: adjusted coordinates of the nodes of "<< vset.Elements() <<" elements.\n";

    vset.EstablishZeroBasedNumbering();
  }





void Quadrilaterator::GenerateVSet( VSet<2U>& vset ) const
 {
   deque<vector<int64_t> >     plist( n_elements, vector<int64_t>(4) ); 
   deque<vector<int64_t> >     pfvert( n_elements, vector<int64_t>(4) );
   //FEM_Data<ScalarVariable >  edata( ELEMENT, n_elements );
   PropertyData               edata( ELEMENT, SCALAR, 2U );
   
   unsigned int    i, j, n; 
   int64_t           fed1, fed2, fed3, fed4;
  
   // n counts the elements
   n = 0;
   
   edata.Resize( n_elements );
   
   for( i=1; i<rows; i++ )
     for( j=1; j<cols; j++ ) {
 
      	 // 1.0 assigning the permeability 
      	 // either taking the permeability of the local node 0
      	 // or using a harmonic averaging of all 4 nodes 
         if ( !harmonic ) edata.Value(n) = matrix_( i, j-1 );
         else             edata.Value(n) = HarmonicPermeabilityAverage( i-1, j-1, matrix_ );
         
         // 2.0 assign the node ids
         // local node numbering (ccw)
         // in accordance with IsoparametricLinearQuadrilateral
         // 3 -- 2
         // |    |
         // |    |
         // 0 -- 1
         plist[n][0] = MapVertex( (i+1),  j ); 
         plist[n][1] = MapVertex( (i+1), (j+1) ); 
         plist[n][2] = MapVertex(  i,    (j+1) ); 
         plist[n][3] = MapVertex(  i,     j );  

 
         // 3.0 get the node coordinates 
         vset.Py( plist[n][0], i + 1 );
         vset.Px( plist[n][0], j );
         vset.Py( plist[n][1], i + 1 );
         vset.Px( plist[n][1], j + 1 );
         vset.Py( plist[n][2], i );    
         vset.Px( plist[n][2], j + 1 );
         vset.Py( plist[n][3], i );      	
         vset.Px( plist[n][3], j );		    

        // 5. assigning boundary flags (node) 
        if ( i == 1 ) {
            if ( j == 1 )        vset.BFlag( plist[n][3], CNR_MAX_MINXZ );  // CNR4
            else                 vset.BFlag( plist[n][3], TOP_OUTSIDE );
            if ( j == (cols-1) ) vset.BFlag( plist[n][2], CNR_MAX_MAXX );   // CNR3
            else                 vset.BFlag( plist[n][2], TOP_OUTSIDE );
          }
        if ( i == (rows-1) ) {
            if ( j == 1 )        vset.BFlag( plist[n][0], CNR_MIN );        // CNR1
            else                 vset.BFlag( plist[n][0], BOTTOM_OUTSIDE );
            if ( j == (cols-1) ) vset.BFlag( plist[n][1], CNR_MIN_MAXX );   // CNR2
            else                 vset.BFlag( plist[n][1], BOTTOM_OUTSIDE );
          }
        if ( j == 1 ) {
            if ( i == 1 )        vset.BFlag( plist[n][3], CNR_MAX_MINXZ );  // CNR4
            else                 vset.BFlag( plist[n][3], LEFT_OUTSIDE );
            if ( i == (rows-1) ) vset.BFlag( plist[n][0], CNR_MIN );        // CNR1
            else                 vset.BFlag( plist[n][0], LEFT_OUTSIDE );
          }
        if ( j == (cols-1) ) {
            if ( i == 1 )        vset.BFlag( plist[n][2], CNR_MAX_MAXX );   // CNR3
            else                 vset.BFlag( plist[n][2], RIGHT_OUTSIDE );
            if ( i == (rows-1) ) vset.BFlag( plist[n][1], CNR_MIN_MAXX );   // CNR2
            else                 vset.BFlag( plist[n][1], RIGHT_OUTSIDE );
          }

 
        // 6.getting the ids of the elements neighboring the current element
        // fed1 = element neighboring local nodes 0 and 1 (bottom)
        // fed2 = element neighboring local nodes 1 and 2 (right)
        // fed3 = element neighboring local nodes 2 and 3 (top)
        // fed4 = element neighboring local nodes 3 and 4 (left)        
        if ( i == 1 )        fed3 = TOP_OUTSIDE;
        else                 fed3 = static_cast<int64_t>((n+1)-(rows-1) - 1);
        if ( i == (rows-1) ) fed1 = BOTTOM_OUTSIDE;
        else                 fed1 = static_cast<int64_t>((n+1)+(rows-1) - 1);
        if ( j == 1 )        fed4 = LEFT_OUTSIDE;
        else                 fed4 = static_cast<int64_t>((n+1)-2);
        if ( j == (cols-1) ) fed2 = RIGHT_OUTSIDE;
        else                 fed2 = static_cast<int64_t>((n+1));
        pfvert[n][0] = fed1;
        pfvert[n][1] = fed2;
        pfvert[n][2] = fed3;
        pfvert[n][3] = fed4;
        
        // 7. increment the elemend id
        n++;
	   
     } // end j-loop 
   
   // setting up the vset
   vset.AddPlist( plist.begin(), plist.end() );
   vset.AddPfverts( pfvert.begin(), pfvert.end() );
   vset.AddData( "permeability", edata );
   
   // Flipping the Y-axis
   double ymax = vset.Py(0U);
   for ( i=1; i<vset.Vertices(); i++ ) 
       if ( vset.Py(i) > ymax ) ymax = vset.Py(i);
   ymax += 1.0; // scale from 1 to rows      
   for ( i=0; i<vset.Vertices(); i++ ) {
        vset.Py( i, ymax - vset.Py(i) );
        vset.Pz( i, 0. );
     }
 }

} // end namespace csmp


