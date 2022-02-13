#include "TransportVisitor2D.h"
#include "ErrorHandler.h"
#include "InterFace.h"
#include "Point.h"
#include "Region.h"
#include "Model.h"
#include "MJL_Edge.h"
#include "Element.h"
#include "FemToGridVisitor.h"
#include "FemFromGridVisitor.h"
#include "CSMP_highLevelUtilities.h"
#include "Exception.h"

using namespace std;

namespace csmp {


TransportVisitor2D::TransportVisitor2D( Model<2>& sg, 
                                        const char* advected_prop, 
                                        const char* advecting_prop )
    : pref(sg.Database()),
      egrids(sg.Region("Model").Elements()), 
      visited(sg.Region("Model").Elements()),
      XY(3,2), NN(3,3), P(3), xy(2),
      time_increment(1),      
      v_key(pref.StorageKey(advecting_prop)),
      prop_key(pref.StorageKey(advected_prop)),
      max_increments(500), 
      interpolate_only_within_grid(false)
  { 
     visited.SetAll(false);
     strcpy( transp_prop, advecting_prop );
     strcpy( adv_prop,    advected_prop );
    
     if ( v_key.place != NODE )
       throw csmp::Exception( FATAL_ERROR, "TransportVisitor2D::(constructor)", 
          "The transport variable must be placed on the Node." ); 		       

     if ( v_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "TransportVisitor2D::(constructor)", 
          "The transport variable must be a vector variable." );
           		       
     if ( prop_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "TransportVisitor2D::(constructor)", 
          "The transported variable must be a scalar variable." );
          
     if ( sg.Mesh().HybridElementMesh() ||
          (!sg.FE_Manager().ContainsElementType( LINEAR_TRIANGLE ) &&
           !sg.FE_Manager().ContainsElementType( ISOPARAMETRIC_QUADRATIC_TRIANGLE ) &&
           !sg.FE_Manager().ContainsElementType( ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) ) )
       throw csmp::Exception( FATAL_ERROR, "TransportVisitor2D::(constructor)", 
          "This Visitor only works for triangular element meshes." );

     // initializing the FiniteDifferenceGrid
     sg.AssignElementCharacteristicsTo("inner radius", "inner radius");
     double          rmin, rmax;
     sg.MinMaxOf("inner radius", rmin, rmax );
     resolution = rmin;

     csmp::Point<2>  xyz_min, xyz_max;
     sg.MinMaxCoordinates( xyz_min, xyz_max );
     //                            x-max                   y-max
     grid.Initialize(  xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], resolution, resolution );
  }
         
      
TransportVisitor2D::~TransportVisitor2D() 
 {  
 }



size_t TransportVisitor2D::MaximumIncrements() const { return max_increments; }

void TransportVisitor2D::MaximumIncrements( size_t maxi ) { max_increments = maxi; }


void TransportVisitor2D::InterpolateOnlyWithinGrid( bool do_so )
 {
    interpolate_only_within_grid = do_so;
 }
 
 


void TransportVisitor2D::Visit( Model<2U>* n )   
     { 
        char name[NAME_STRING];
        // Model has no id's therefore give element id's
        // get into the target element
        cout <<"\nTransportVisitor2D::VisitModel: please enter name of ";
        cout <<" region (Region) from which you would like to spark streamlines\n";
        cin >> name;
        n->Region(name).Accept( *this );
     }   





void TransportVisitor2D::MinMaxCoordinates( double& min_x, double& max_x, 
                                            double& min_y, double& max_y )
 {
    min_x = max_x = XY(0,0);
    min_y = max_y = XY(0,1);
    for ( size_t i=1U; i<XY.Rows(); i++ )
      {
         if ( XY(i,0) < min_x ) min_x = XY(i,0);
         if ( XY(i,0) > max_x ) max_x = XY(i,0);
         if ( XY(i,1) < min_y ) min_y = XY(i,1);
         if ( XY(i,1) > max_y ) max_y = XY(i,1);
      }
 }



bool  TransportVisitor2D::IsInsideTriangle( double x, double y, bool update )
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
            for ( int32_t r=0; r<3; r++ ) a[r].Flip();
       }
     // TEST: if the midpoint does not lie to the right of each edge
     // the edges are flipped to change the sense of rotation
     // of the triangle
     for ( int32_t q=0; q<3; q++ )
       if ( a[q].Classify(p) == mjl::RIGHT ) return false;
       
     return true;
  }




void TransportVisitor2D::InitializeElementGrid( size_t idx )
 {
    double  min_x, max_x, min_y, max_y;
    int32_t      i_min, j_min, i_max, j_max;

    MinMaxCoordinates( min_x, max_x, min_y, max_y );      
    grid.ClosestGridPointTo( min_x, min_y, i_min, j_min ); 
    grid.ClosestGridPointTo( max_x, max_y, i_max, j_max ); 

    // accumulating grid points into ElementGrid
    for ( int32_t i=i_min; i<=i_max; i++ )
      for ( int32_t j=j_min; j<=j_max; j++ ) 
        {
           // testing whether point lies within triangle
           if ( IsInsideTriangle(grid.X(j),grid.Y(i)) )
             egrids[ idx ].AddPoint( i, j, grid(i,j) );
        }
    if ( egrids[ idx ].Empty() ) {
         cout <<"\nTransportVisitor2D::InitializeElementGrid: "; 
         cout <<" No grid point values will be assigned to element: "<< idx << endl;
      }
    else visited.SetBit( idx, true );
 }






void TransportVisitor2D::Visit( Element<2U>* n )   
  { 
     // getting the node properties
     n->NodePropertyVector( v_key, P );
     
     // initializing egrid of element if this has not been done before
     if ( !visited.GetBit( n->Idx() ) ) {
          n->NodeCoordinateMatrix( XY );
          InitializeElementGrid( n->Idx() );
       }
     
     // getting coordinates 
     map<pair<int32_t,int32_t>,double>::iterator eit    = egrids[ n->Idx() ].Begin();
     map<pair<int32_t,int32_t>,double>::iterator it_end = egrids[ n->Idx() ].End();
     
     while ( eit != it_end )
       {
          // getting the point in question
          // x,y-coordinates
          xy[0] = grid.X( (*eit).first.second );
          xy[1] = grid.Y( (*eit).first.first );
          
          // computing test-function values
          (*n).N_AtGlobalPoint( n->FE()->NRST, xy );

          // getting x and y velocity components at point
          // by summing up the testfunction values at the point
          double  dvx(0.), dvy(0.);
          for ( auto i{0}; i<n->Nodes(); i++ ) {
               dvx += n->FE()->NRST[i] * P[i](0);
               dvy += n->FE()->NRST[i] * P[i](1);
            }
          // getting place from where fluid comes (upstream) and interpolating 
          // the corresponding concentration value
          double px = xy[0] - dvx * time_increment;
          double py = xy[1] - dvy * time_increment;
          
          if ( interpolate_only_within_grid ) 
            (*eit).second = grid.InterpolateWithin( px, py );
          else 
            (*eit).second = grid( px, py, false );
          eit++; 
       }
       
} // end VisitElement
    
    
    
bool TransportVisitor2D::AdvectUntil( Model<2U>& sg, double final_time )
 {
    int32_t  n(0);
    double  time(0.), vmin, vmax, p0min, p0max, p1min, p1max;

    Region<2>&  sgroup(sg.Region("Model"));
    sgroup.RenumberNodes();

    // getting range of transported property before advection
    sg.MinMaxOf( adv_prop, p0min, p0max );
    cout <<"\nTransportVisitor2D::AdvectUntil:"<< endl;
    cout <<"\n    Value range of advected property: ";
    cout << p0min <<" - "<< p0max << endl;

    // getting velocity range
    sg.MinMaxOf( transp_prop, vmin, vmax );
    cout <<"\n    Transport velocity range [m s-1]: "<< vmin <<" - "<< vmax << endl;
    time_increment = (resolution*2.) / vmax;       
    cout <<"\n    Advection increment: "<< time_increment <<" s ";
    cout <<"("<< (time_increment/final_time) * 100. <<"\n% of total transport time interval)."<< endl;
           
    // advecting on grid; limiting advection increment
    WritePropertyToGrid( sg, adv_prop );
  
    cout <<"\n    Computing advection steps..."<< endl;
    while ( time < final_time )
      {
         if ( (final_time-time) < time_increment ) time_increment = final_time - time;
         cout <<"      advection step: "<< ++n << endl;
         sg.Accept( *this );
         StoreResultsInGrid();
         grid.SaveToJPG( "test-grid", n );
         time += time_increment;
         if ( n >= static_cast<int32_t>(max_increments) ) break;
      }

    InputPropertyFromGrid( sg, adv_prop );

    // getting range of transported property after advection
    sg.MinMaxOf( adv_prop, p1min, p1max );
    
    // Checking whether advection did not screw up advected property range
    if ( p1min < p0min ) {
        cerr <<"\noriginal minimum: "<< p0min <<" vs. new minimum: "<< p1min << endl;
        throw csmp::Exception( INFO,"TransportVisitor2D::AdvectUntil","advection changed property minimum");
        return false;
      }
    if ( p1max > p0max ) {
        cerr <<"\noriginal maximum: "<< p0max <<" vs. new maximum: "<< p1max << endl;
        throw csmp::Exception( INFO,"TransportVisitor2D::AdvectUntil","advection changed property maximum");
        return false;
      }
    return true;

 } // end AdvectUntil




void TransportVisitor2D::StoreResultsInGrid()
 {
    for ( vector<ElementGrid>::const_iterator
          it=egrids.begin(); it!=egrids.end(); it++ )
      (*it).TransferDataToGrid( grid );
    
 } // end  




/**
 
Assign a property stored on a FiniteDifferenceGrid to nodes, constraint
points or elements of a finite-element mesh inside the Model object.
To find node and constraint points values, bilinear interpolation on 
the grid is used. If element properties are desired, a FemFromGridVisitor 
class object will determine the grid points which lie inside each element 
and assigns their arithmetic mean property value as element property.
The property which was input from the FiniteDifferenceGrid class object
will have the flag 'PLAIN'. 

Ideally, the grid from which the property is input should have exactly the
same dimensions as the Model finite-element mesh. The method will 
however work also, if the grid is smaller than the finite-element mesh.
In this case only some of the finite-elements will be assigned values 
from the grid. 

@section arguments Input Arguments 

The name of the input property and the FiniteDifferenceGrid object from
which the property shall be input. 

@section implementation Implementation

InputPropertyFromGrid() works only for two-dimensional CSMP models. 

Only scalar properties can be assigned using from FiniteDifferenceGrids.
These are mapped to CSMP model physical variables using the 
FemFromGridVisitor class object. 

If one wants to assign vector or tensor properties from grid. Several grids
must be supplied to constrain the different dimensions of each variable.
In this case you can create temporary scalar variables for each dimension
and assign these to the different components of a vector or tensor 
variable through a newly defined Interrelation subclass. 

@section application Application

InputPropertyFromGrid() was designed to meet two objectives. Firstly, CSP model 
properties need to be transferred to and from regular grids for visualization
and advection calculation purposes. Secondly, input properties for models
may be obtained by analysis of photographs coupled with color-coding of 
the thresholded and manipulated image data (see for instance NIH Image 
documentation to get an idea of what one can do). There may be a need to
map such data onto a mesh. 

@section messages Messages 

If any of the grid dimensions are larger than the dimensions of the 
finite-element mesh, the simulation will be terminated by the C assert 
function.  

Thus far, only scalar node and element properties can be mapped from
a grid. Therefore, if the target property is not a scalar or is a constraint 
point variable, an error will be raised and the method wil return without 
a mapping. 
*/
void  TransportVisitor2D::InputPropertyFromGrid( Model<2U>& sg, 
                                                 const char* prop ) const
 {
    
    csmp::Index                prop_key1 = sg.Database().StorageKey( prop );
    static bool                checked  = false;
    Region<2>&       super_group(sg.Region("Model"));
    ScalarVariable  res;
    
    // checking whether grid dimension match that of Model mesh
    if ( !checked ) {
         csmp::Point<2U>  xyz_min, xyz_max;
         sg.MinMaxCoordinates( xyz_min, xyz_max );
         if ( grid.MinX() < xyz_min[0] ||
              grid.MinY() < xyz_min[1] ||
              grid.MaxX() > xyz_max[0] ||
              grid.MaxY() > xyz_max[1] ) {
              cerr <<"\nTransportVisitor2D::WritePropertyToGrid: ";
              cerr <<"Grid extent is larger than that of Model.\nGrid therefore cannot be ";
              cerr <<"mapped to Model." << endl;
              throw domain_error("TransportVisitor2D::WritePropertyToGrid");
           }
         checked = true;
      }

    if ( prop_key1.type != SCALAR ) {
         throw csmp::Exception( ERROR, "TransportVisitor2D::InputPropertyFromGrid",
                                    "only scalar properties can be read from a single grid.");
         return;
      }
    if ( prop_key1.place == NODE ) {
        for ( auto nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
             res.Flag() = (*nit)->Status( prop_key1 );
             res = grid( (*nit)->x(), (*nit)->y(), false );
             (*nit)->Store( prop_key1, res ); 
          }
        return;
      }

    // visitor gets arithmetic means from grid if element properties need to be calculated
    FemFromGridVisitor<2U>  reader( sg.Database(), grid, prop, super_group.Elements() );

    if ( prop_key1.place == ELEMENT ) sg.Accept( reader );
    else
      throw csmp::Exception( ERROR, "TransportVisitor2D::InputPropertyFromGrid",
                        "only node and element properties can be read from grid.");

 } // end InputPropertyFromGrid
 
 
 
 
 
 
/**
 
Writes a distributed scalar physical variable to a FiniteDifferenceGrid
class object. If the variable is a node variable, the interpolation
functions of the current finite elements are used to interpolate the 
property onto the grid. Element properties are constant on each element 
and are mapped as such onto the FiniteDifferenceGrid. WritePropertyToGrid() 
depends on a FiniteDifferenceGrid which is fine enough to allow at least 
one grid point to be assigned to each element. Also the grid must be 
smaller or of the same size as the finite-element mesh.  

@section arguments Input Arguments 

The names of the input property and the FiniteDifferenceGrid object to
which the property shall be output. 

@section implementation Implementation

WritePropertyToGrid() works only for two-dimensional CSP models. 

The method uses a FemToGridVisitor class object to map the physical 
variable to the grid. 

@section application Application

If one depends on a regular-gridded dataset for the visualization of a
variable (as is the case if one uses Spyglass Transform or other similar 
visualization tools), this method allows to output the physical variable 
to the grid in an efficient way. The grid, in turn, can be output as a 
regular-gridded dataset, HDF file or directly as a JPEG image file. 

Another application is the mapping of variables onto a grid on which the
latter can then be advected using a finite difference method or a 
TransportVisitor class object.  

@section messages Messages 

The method reports if the grid is too small or if an erratic attempt
is made to write a 3D model to the two-dimensional grid.  
*/
void  TransportVisitor2D::WritePropertyToGrid( Model<2U>& sg, 
                                               const char* prop )
 {
    static bool                checked(false);
    const Region<2>&  super_group(sg.Region("Model"));
    
    // checking whether grid dimension match that of Model mesh
    if ( !checked ) {
         csmp::Point<2U>  xyz_min, xyz_max;
         sg.MinMaxCoordinates( xyz_min, xyz_max );
         if ( grid.MinX() < xyz_min[0] ||
              grid.MinY() < xyz_min[1] ||
              grid.MaxX() > xyz_max[0] ||
              grid.MaxY() > xyz_max[1] ) {
              cerr <<"\nTransportVisitor2D::WritePropertyToGrid: ";
              cerr <<"Grid extent is smaller  than that of Model.\nGrid therefore cannot be used ";
              cerr <<"for the mapping." << endl;
              throw domain_error("TransportVisitor2D::WritePropertyToGrid");
           }
         checked = true;
      }
      
     FemToGridVisitor<2U>  writer( sg.Database(), grid, prop, super_group.Elements() );
     writer.OverWrite( true );
     writer.OutputProperty( prop );
     sg.Accept( writer );
    
 } // end WritePropertyToGrid

    
    
} // csmp    
    
    
    
    
    
    
    
