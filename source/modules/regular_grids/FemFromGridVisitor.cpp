#include "FemFromGridVisitor.h"
#include "Exception.h"
#include "PropertyDatabase.h"
#include "FiniteDifferenceGrid.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
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




      
template<uint32_t dim>
FemFromGridVisitor<dim>::~FemFromGridVisitor() 
 {  
 }

template<uint32_t dim>
void FemFromGridVisitor<dim>::MinMaxCoordinates( double& min_x, double& max_x, 
                                                 double& min_y, double& max_y )
 {
    min_x = max_x = XY(0,0);
    min_y = max_y = XY(0,1);
    for ( uint32_t i=1u; i<XY.Rows(); i++ )
      {
         if ( XY(i,0) < min_x ) min_x = XY(i,0);
         if ( XY(i,0) > max_x ) max_x = XY(i,0);
         if ( XY(i,1) < min_y ) min_y = XY(i,1);
         if ( XY(i,1) > max_y ) max_y = XY(i,1);
      }
 }



template<uint32_t dim>
bool FemFromGridVisitor<dim>::IsInsideTriangle(double x, double y, bool /*update*/ )
{
    // 1. Setup points
    Point<2> p1(XY(0,0), XY(0,1)), 
             p2(XY(1,0), XY(1,1)), 
             p3(XY(2,0), XY(2,1)), 
             p(x,y);

    // 2. Cross product helper function (2D version)
    // Returns positive if P is to the LEFT of vector AB, negative if to the RIGHT
    auto cross_product = [](const Point<2>& a, const Point<2>& b, const Point<2>& p) {
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
    };

    // 3. Determine orientation (CCW vs CW)
    // The proprietary code used a midpoint test to handle orientation.
    // We can simply calculate the signed area (cross product) of the triangle.
    double area = cross_product(p1, p2, p3);
    
    // 4. Perform the "Point in Triangle" test
    double d1 = cross_product(p1, p2, p);
    double d2 = cross_product(p2, p3, p);
    double d3 = cross_product(p3, p1, p);

    if (area > 0) {
        // Counter-Clockwise triangle: Point must be to the left of all edges
        return (d1 >= 0 && d2 >= 0 && d3 >= 0);
    } else {
        // Clockwise triangle: Point must be to the right of all edges
        return (d1 <= 0 && d2 <= 0 && d3 <= 0);
    }
}


template<uint32_t dim>
bool FemFromGridVisitor<dim>::IsInsideQuadrilateral(double x, double y)
{
    // 1. Setup points for the 4 corners
    Point<2> p[4] = {
        Point<2>(XY(0,0), XY(0,1)),
        Point<2>(XY(1,0), XY(1,1)),
        Point<2>(XY(2,0), XY(2,1)),
        Point<2>(XY(3,0), XY(3,1))
    };
    Point<2> p_test(x, y);

    // 2. Cross product helper: returns positive if P is to the LEFT of AB
    auto cross_product = [](const Point<2>& a, const Point<2>& b, const Point<2>& p) {
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
    };

    // 3. Determine winding order using the first three points
    double area = cross_product(p[0], p[1], p[2]);
    
    // 4. Check against all 4 edges
    // We use a small epsilon or simply 0.0 with '=' to include edges
    if (area > 0) {
        // Counter-Clockwise: Point must be Left/On all edges
        for (int i = 0; i < 4; ++i) {
            if (cross_product(p[i], p[(i + 1) % 4], p_test) < -1e-12) return false;
        }
    } else {
        // Clockwise: Point must be Right/On all edges
        for (int i = 0; i < 4; ++i) {
            if (cross_product(p[i], p[(i + 1) % 4], p_test) > 1e-12) return false;
        }
    }

    return true;
}


template<uint32_t dim>
void FemFromGridVisitor<dim>::InitializeElementGrid( size_t idx, CSMP_FEM_TYPE fe_type )
 {
    double     min_x, max_x, min_y, max_y;
    int32_t  i, j, i_min, j_min, i_max, j_max;

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
template<uint32_t dim>
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
