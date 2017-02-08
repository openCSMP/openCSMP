#include "FemToGridVisitor.h"
#include "Exception.h"
#include "MJL_Point.h"
#include "MJL_Edge.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "FiniteDifferenceGrid.h"
#include "Element.h"

using namespace std;

namespace csmp {


template<size_t dim>
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
        throw csmp::Exception( CSMP_FATAL_ERROR, "FemToGridVisitor(constructor 3D)",
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



      
template<size_t dim>
FemToGridVisitor<dim>::~FemToGridVisitor() 
 {  
 }



template<size_t dim>
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






template<size_t dim>
bool FemToGridVisitor<dim>::OverWrite() const { return overwrite; }



template<size_t dim>
void FemToGridVisitor<dim>::OverWrite( bool ow ) { overwrite=ow; }


/* 
template<size_t dim>
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


template<size_t dim>
void FemToGridVisitor<dim>::MinMaxCoordinates( double64& min_x, double64& max_x, 
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
bool  FemToGridVisitor<dim>::IsInsideTriangle( double64 x, double64 y, bool update )
  {
     mjl::Point p1(XY(0,0),XY(0,1)), 
               p2(XY(1,0),XY(1,1)), 
               p3(XY(2,0),XY(2,1)), 
               p(x,y), mp(XY(0,0),XY(0,1));
               
     static mjl::Edge  a[3];
 
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
     for ( int32 q=0; q<3; q++ )
       if ( a[q].Classify(p) == mjl::RIGHT ) return false;

     return true;
  }

template<size_t dim>
bool  FemToGridVisitor<dim>::IsInsideQuadrilateral( double64 x, double64 y )
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
void FemToGridVisitor<dim>::InitializeElementGrid( size_t idx, CSMP_FEM_TYPE fe_type  )
 {
    double64  min_x, max_x, min_y, max_y;
    int32      i, j, i_min, j_min, i_max, j_max;

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
           if ( fe_type == ISOPARAMETRIC_LINEAR_QUADRILATERAL || 
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
template<size_t dim>
void FemToGridVisitor<dim>::Visit( Region<dim>* r )
 {
    r->RenumberElements();
    for ( auto it=r->ElementsBegin(); it!=r->ElementsEnd(); ++it )
      Visit( *it );
 }

/**

Using the element interpolation functions, the property is interpolated
onto the FiniteDifferenceGrid object.  

@param n A pointer to the Element in which the values shall be interpolated.

@section application Application 

To output model results to JPEG images or other gridded data.  
*/
template<size_t dim>
void FemToGridVisitor<dim>::Visit( Element<dim>* n )   
  { 
     double64 val;
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
          n->CoordinateMatrix( XY );
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
    

  
    
        
    
template<size_t dim>
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
 



template<size_t dim>
const char* FemToGridVisitor<dim>::OutputProperty() const
 {
    return pref.Name( prop_key );
 }
    
    
    
template class FemToGridVisitor<2U>;
    
    
} // end namespace csmp
