#include "FemFromGridVisitor.h"
#include "Exception.h"
#include "MJL_Point.h"
#include "MJL_Edge.h"
#include "PropertyDatabase.h"
#include "FiniteDifferenceGrid.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"

using namespace std;

namespace csmp {

template<size_t dim>
FemFromGridVisitor<dim>::FemFromGridVisitor( const PropertyDatabase<dim>& p, 
                                             const FiniteDifferenceGrid& g, 
                                             const char* var,
                                             size_t elements )
    : pref(p), grid(g), 
      egrids(elements), 
      XY(3,2),
      NN(3,3),
      key(pref.StorageKey(var))
  { 
     if ( dim == 3U ) 
        throw csmp::Exception( FATAL_ERROR, "FemFromGridVisitor(constructor 3D)",
                       "Thus far, 2D-models only, can be written to grids."); 

     if ( key.type != SCALAR )
       {
          cout <<"\nFemFromGridVisitor: only scalar properties can be ";
          cout <<"mapped from regular grid so far. Now exciting..." << endl;
          if ( key.type != SCALAR )
            throw csmp::Exception( ERROR, "FemFromGridVisitor<dim>(constructor)",
                            var, "basic property must be a scalar.");
       }
     if ( key.place == INTER_FACE )
       {
          cout <<"\nFemFromGridVisitor: only nodal and element properties can be ";
          cout <<"mapped from regular grid so far. Now exciting..." << endl;
          if ( key.place == INTER_FACE )
            throw csmp::Exception( ERROR, "FemFromGridVisitor<dim>(constructor)",
                            var, "basic property must not be placed on Face.");
       }
  }




      
template<size_t dim>
FemFromGridVisitor<dim>::~FemFromGridVisitor() 
 {  
 }

template<size_t dim>
void FemFromGridVisitor<dim>::MinMaxCoordinates( double64& min_x, double64& max_x, 
                                                    double64& min_y, double64& max_y )
 {
    min_x = max_x = XY(0,0);
    min_y = max_y = XY(0,1);
    for ( size_t i=1; i<XY.Rows(); i++ )
      {
         if ( XY(i,0) < min_x ) min_x = XY(i,0);
         if ( XY(i,0) > max_x ) max_x = XY(i,0);
         if ( XY(i,1) < min_y ) min_y = XY(i,1);
         if ( XY(i,1) > max_y ) max_y = XY(i,1);
      }
 }



template<size_t dim>
bool  FemFromGridVisitor<dim>::IsInsideTriangle( double64 x, double64 y, bool update )
  {
     mjl::Point p1(XY(0,0),XY(0,1)), 
               p2(XY(1,0),XY(1,1)), 
               p3(XY(2,0),XY(2,1)), 
               p(x,y), mp(XY(0,0),XY(0,1));
               
     mjl::Edge a[3];
 
     if ( update )
       {
          mp += p2;
          mp += p3;
          mp /= 3.0;
          a[0].Set( p1, p2 );
          a[1].Set( p2, p3 );
          a[2].Set( p3, p1 );
          
          // flipping segments if triangles are numbered counter-clockwise
          if ( a[0].Classify(mp) == mjl::RIGHT )
            for ( int32 r=0; r<3; r++ ) a[r].Flip();
       }
     // TEST: if the midpoint does not lie to the right of each edge
     // the edges are flipped to change the sense of rotation
     // of the triangle
     for ( size_t q=0; q<3; q++ )
       if ( a[q].Classify(p) == mjl::RIGHT ) return false;
       
     return true;
  }


template<size_t dim>
bool  FemFromGridVisitor<dim>::IsInsideQuadrilateral( double64 x, double64 y )
  {
     // this method currently works only for regular rectangles
     // test if the quadrilateral is straightsided and regular
     // if MJL_Points and MJL_Edges are used to identify if points
     // lie within the quadrilateral, as in the function IsInsideTriangle()
     // sometimes FD points that lie on the edge of the quadrilateral 
     // are flagged as 'outside' and the JPEG output is crappy
     if ( XY(0,1) != XY(1,1) || XY(2,1) != XY(3,1) ||  // y coordinate
          XY(0,0) != XY(3,0) || XY(1,0) != XY(2,0) )   // x coordinate
       { 
         cout << "\nFemToGridVisitor<dim>::IsInsideQuadrilateral: ";
         cout << "\nApparently quadrilateral is not straightsided and rectangular";
         cout << "\nMethod cannot identify if point lies within quadrilateral";
         cout << "\nNothing is done..." << endl;
         return false;
       }
          
     if ( x >= XY(0,0) && x <= XY(1,0) && y >= XY(0,1) && y <= XY(3,1) ) return true;
     return false;
  }


template<size_t dim>
void FemFromGridVisitor<dim>::InitializeElementGrid( size_t idx, CSMP_FEM_TYPE fe_type )
 {
    double64     min_x, max_x, min_y, max_y;
    int32  i, j, i_min, j_min, i_max, j_max;

    MinMaxCoordinates( min_x, max_x, min_y, max_y );      
    grid.ClosestGridPointTo( min_x, min_y, i_min, j_min ); 
    grid.ClosestGridPointTo( max_x, max_y, i_max, j_max ); 

    // accumulating grid points into ElementGrid
    for ( i=i_min; i<=i_max; i++ )
      for ( j=j_min; j<=j_max; j++ ) 
        {
           // testing whether point lies within triangle
           if ( fe_type == ISOPARAMETRIC_LINEAR_TRIANGLE ||
                fe_type == LINEAR_TRIANGLE ||
                fe_type == ISOPARAMETRIC_QUADRATIC_TRIANGLE )
               if ( IsInsideTriangle(grid.X(j),grid.Y(i)) ) 
                 egrids[ idx ].AddPoint( i, j, grid(i,j) );
           if ( fe_type == ISOPARAMETRIC_LINEAR_QUADRILATERAL || fe_type == LINEAR_RECTANGLE ||
                fe_type == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) 
               if ( IsInsideQuadrilateral(grid.X(j),grid.Y(i)) ) 
                 egrids[ idx ].AddPoint( i, j, grid(i,j) );
        }
    if ( egrids[ idx ].Empty() )
      {
         cout <<"\nFemFromGridVisitor::InitializeElementGrid: "; 
         cout <<" No grid point values will be assigned to element: "<< idx << endl;
      }
 }


// Visit(Element* n) 
// ------------------------
template<size_t dim>
void FemFromGridVisitor<dim>::Visit( Element<dim>* n )   
  { 
     // getting to element
     size_t e_id = n->Idx();
     n->NodeCoordinateMatrix( XY );

     // reading values from grid
     InitializeElementGrid( e_id, n->FE_Type() );
     
     // getting grid-point indices and values for element 
     // egrids[ e_id ].Out();
     val() = egrids[ e_id ].GridAverage();
     val.Flag() = PLAIN;

     if ( n->Status( key ) == PLAIN )
       n->Store( key, val );

} // end Visit


template class FemFromGridVisitor<2U>;

} // end namespace csmp
