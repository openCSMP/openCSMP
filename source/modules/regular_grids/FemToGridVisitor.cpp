// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FemToGridVisitor.h"
#include "Exception.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "FiniteDifferenceGrid.h"
#include "Element.h"

using namespace std;

namespace csmp {


template<uint32_t dim>
FemToGridVisitor<dim>::FemToGridVisitor( const PropertyDatabase<dim>&  p, 
                                         FiniteDifferenceGrid& g, 
                                         const char*           var,
                                         size_t                elements )
    : pref(p), grid(g), 
      visited(elements,false), 
      egrids(elements), 
      NF(3),
      xy(2),
      P(3),
      XY(3,2), 
      NN(3,3),
      prop_key(p.StorageKey(var)),
      overwrite(true)
  { 
     if ( dim == 3U )
        throw csmp::Exception( FATAL_ERROR, "FemToGridVisitor(constructor 3D)",
                              "Thus far, 2D-models only, can be written to grids."); 

     if ( prop_key.type != SCALAR )
       {
          cout <<"\nFemToGridVisitor: only scalar properties can be ";
          cout <<"mapped on regular grid so far. Now exciting..." << endl;
       }
     if ( prop_key.place != ELEMENT && prop_key.place != NODE )
       {
          cout <<"\nFemToGridVisitor: only nodal and element properties can be ";
          cout <<"mapped on regular grid so far. Now exciting..." << endl;
       }
  }



      
template<uint32_t dim>
FemToGridVisitor<dim>::~FemToGridVisitor() 
 {  
 }



template<uint32_t dim>
FemToGridVisitor<dim>&  FemToGridVisitor<dim>::operator=( const FemToGridVisitor<dim>& vis )
 {
    if ( &vis != this )
      {
         cout <<"\nFemToGridVisitor::operator=: Assignment made did not change any references"<< endl;
         egrids    = vis.egrids; 
         visited   = vis.visited;
         XY        = vis.XY;
         P         = vis.P;
         NN        = vis.NN;
         overwrite = vis.overwrite;
         assert( prop_key == vis.prop_key );
      }
    return *this;
 }






template<uint32_t dim>
bool FemToGridVisitor<dim>::OverWrite() const { return overwrite; }



template<uint32_t dim>
void FemToGridVisitor<dim>::OverWrite( bool ow ) { overwrite=ow; }


/* 
template<uint32_t dim>
void FemToGridVisitor<dim>::ChangeOutputProperty( const char* var )
 {
     prop_key = pref.StorageKey( var );

     if ( prop_key.type != SCALAR )
       {
          cout <<"\nFemToGridVisitor::ChangeOutputProperty: only scalar properties can be ";
          cout <<"mapped on regular grid so far. Now exciting..." << endl;
       }
     if ( prop_key.place == EDGE || prop_key.place == INTER_FACE )
       {
          cout <<"\nFemToGridVisitor::ChangeOutputProperty:";
          cout <<" only nodal and element properties can be ";
          cout <<"mapped on regular grid so far. Now exciting..." << endl;
       }
 }
*/


template<uint32_t dim>
void FemToGridVisitor<dim>::MinMaxCoordinates( double& min_x, double& max_x, 
                                               double& min_y, double& max_y )
 {
    min_x = max_x = XY(0,0);
    min_y = max_y = XY(0,1);
    for ( uint32_t i{1u}; i<XY.Rows(); i++ )
      {
         if ( XY(i,0) < min_x ) min_x = XY(i,0);
         if ( XY(i,0) > max_x ) max_x = XY(i,0);
         if ( XY(i,1) < min_y ) min_y = XY(i,1);
         if ( XY(i,1) > max_y ) max_y = XY(i,1);
      }
 }


template<uint32_t dim>
bool FemToGridVisitor<dim>::IsInsideTriangle(double x, double y, bool /*update*/ )
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
bool FemToGridVisitor<dim>::IsInsideQuadrilateral(double x, double y)
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
void FemToGridVisitor<dim>::InitializeElementGrid( size_t idx, CSMP_FEM_TYPE fe_type  )
 {
    double  min_x, max_x, min_y, max_y;
    int32_t      i, j, i_min, j_min, i_max, j_max;

    MinMaxCoordinates( min_x, max_x, min_y, max_y ); 
    //                                     row    column     
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
         cout <<"\nFemToGridVisitor::InitializeElementGrid: "; 
         cout <<" No grid points were found in element: "<< idx << endl;
      }
    else visited[ idx ] = true;
 }



/// if the visitor is passed to a region, all its elements are processed
template<uint32_t dim>
void FemToGridVisitor<dim>::Visit( Region<dim>* r )
 {
    r->RenumberCells();
    for ( auto it=r->CellsBegin(); it!=r->CellsEnd(); ++it )
      Visit( *it );
 }

/**

Using the element interpolation functions, the property is interpolated
onto the FiniteDifferenceGrid object.  

@param n A pointer to the Element in which the values shall be interpolated.

@section application Application 

To output model results to JPEG images or other gridded data.  
*/
template<uint32_t dim>
void FemToGridVisitor<dim>::Visit( Element<dim>* n )   
  { 
     double val;
     size_t   e_id, i;
     
     // getting to neighbor elements
     e_id = n->Idx();
     if ( prop_key.place == NODE )
       n->NodePropertyVector( prop_key, P );
     else
       n->Read( prop_key, sc );
     
     // initializing egrid of element if this has not been done before
     if ( !visited[ e_id ] ) 
       {
          n->NodeCoordinateMatrix( XY );
          InitializeElementGrid( e_id, n->FE_Type() );
       }
     
     // getting coordinates 
     eit    = egrids[ e_id ].Begin();
     it_end = egrids[ e_id ].End();
     
     while ( eit != it_end )
       {
          if ( prop_key.place == NODE )
            {
               // getting the point in question
               // x,y-coordinates
               xy[0] = grid.X( (*eit).first.second );
               xy[1] = grid.Y( (*eit).first.first );

               // computing test-function values at xy, using the current finite element 
               (*n).N_AtGlobalPoint( NF, xy );

               // gettin 'val' by summing up testfunction values at the point
               for ( val=0., i=0U; i<n->Nodes(); i++ ) val += NF[i] * P[i]();
            }
          else // if ELEMENT property
          val = sc();

          // storing the (element)-interpolated value on the grid
          if ( overwrite )
            grid( (*eit).first.first, (*eit).first.second ) = val;
          else
          grid( (*eit).first.first, 
                (*eit).first.second ) = grid( (*eit).first.first, (*eit).first.second ) + val;
          eit++; 
       }
} // end Visit
    

  
    
        
    
template<uint32_t dim>
void FemToGridVisitor<dim>::OutputProperty( const char* prop )
 {
     Index tmp_key = pref.StorageKey( prop );

     if ( prop_key.type != SCALAR ) {
          cout <<"\nFemToGridVisitor::OutputProperty: only scalar properties can be ";
          cout <<"mapped on regular grid so far. No assignment was made" << endl;
          return;
       }
     prop_key = tmp_key;
 }
 



template<uint32_t dim>
const char* FemToGridVisitor<dim>::OutputProperty() const
 {
    return pref.Name( prop_key );
 }
    
    
    
template class FemToGridVisitor<2U>;
    
    
} // end namespace csmp
