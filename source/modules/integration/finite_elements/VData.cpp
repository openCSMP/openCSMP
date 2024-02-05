#include "VData.h"
#include "Box.h"
#include "FiniteElement.h"
#include "binaryReadWrite.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "TextFileIO.h"
#include "TriangularFacet.h"
#include "CSMP_mathUtilities.h"
#include "CSMP_highLevelUtilities.h"
#include "CSMP_ElementSpecifications.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

/// default constructor: not hybrid, no nodes, nor elements
VData::VData()
 : hybrid_mesh_(false),
   first_interface_(0U),
   first_face_(0U)
 {
 }


/**
   constructor for hybrid element meshes
   
   @param npes 1 stores the number of nodes of each element;
   thus, so its size is equal to the total number
   of elements.

   @param epes 2 stores the number of faces=neighbors per element
   so its size is equal to the total number of elements.
   
   @note it is assumed that there are no faces nor interfaces

*/
VData::VData( const deque<uint32_t>& npes,
              const deque<uint32_t>& epes,
              size_t nodes )
  : px(nodes),
    py(nodes),
    pz(nodes),
    bflags(nodes),
    gflags_(nodes),
    pelmt(epes.size()),
    plist(epes.size()),
    pfverts(epes.size()),
    hybrid_mesh_(true),
    first_face_(npes.size()),
    first_interface_(npes.size())
 {
    assert( !npes.empty() );
    const size_t n_cells{ epes.size() };
    for ( size_t i{0U}; i<n_cells; i++ )
      plist[i]   = vector<int64_t>(npes[i]);
         
    assert( !epes.empty() );
    for ( size_t i{0U}; i<n_cells; i++ )
      pfverts[i] = vector<int64_t>(epes[i]);
 }




/** 
    constructor for meshes with a single type of element

   @note it is assumed that there are no faces nor interfaces
*/
VData::VData( uint32_t nodes_per_element, uint32_t nbors_per_element, size_t nodes, size_t elmts )
  : px(nodes), py(nodes), pz(nodes), bflags(nodes), gflags_(nodes), pelmt(elmts),
    hybrid_mesh_(false),
    first_face_(elmts),
    first_interface_(elmts)
{
    for ( size_t i{0U}; i<elmts; i++ )
      plist.push_back( vector<int64_t>(nodes_per_element) );
    for ( size_t k{0U}; k<elmts; k++ )
      pfverts.push_back( vector<int64_t>(nbors_per_element) );
}




VData::~VData()
{
}


VData::VData( const VData& vd )
 : px(vd.px), py(vd.py), pz(vd.pz), 
   plist(vd.plist), pfverts(vd.pfverts),
   pelmt(vd.pelmt),
   bflags(vd.bflags),
   gflags_(vd.gflags_),
   hybrid_mesh_(vd.hybrid_mesh_),
   first_face_( vd.first_face_ ),
   first_interface_(vd.first_interface_ ),
   pmanifolds_( vd.pmanifolds_ )
{
}





VData& VData::operator=( const VData& a )
{
  if ( &a == this ) return *this;
  
  px               = a.px;
  py               = a.py;
  pz               = a.pz;
  pelmt            = a.pelmt;
  plist            = a.plist;
  pfverts          = a.pfverts;
  bflags           = a.bflags;
  gflags_          = a.gflags_;
  hybrid_mesh_     = a.hybrid_mesh_;
  first_face_      = a.first_face_;
  first_interface_ = a.first_interface_;
  pmanifolds_      = a.pmanifolds_;

  return *this;
}







void  VData::Px( size_t i, double val ) 
{ assert( i<px.size() ); px[i] = val; }

void  VData::Py( size_t i, double val ) 
{ assert( i<py.size() ); py[i] = val; }

void  VData::Pz( size_t i, double val ) 
{ assert( i<pz.size() ); pz[i] = val; }

double VData::Px( size_t i ) const 
{ assert( i<px.size() ); return px[i]; }

double VData::Py( size_t i ) const 
{ assert( i<py.size() ); return py[i]; }

double VData::Pz( size_t i ) const 
{ assert( i<pz.size() ); return pz[i]; }



void  VData::P( uint32_t coordinate_axis, size_t i, double val )
  {
    if( coordinate_axis == 0U ) {
        assert( i<px.size() );
        px[i] = val;
        return;
    }
    
    if( coordinate_axis == 1U ) {
        assert( i<py.size() );
        py[i] = val;
        return;
    }
    
    if( coordinate_axis == 2U ) {
        assert( i<pz.size() );
        pz[i] = val;
        return;
    }
    throw overflow_error( "VData::P: coordinate axis is out of range");
  }



double VData::P( uint32_t coordinate_axis, size_t i ) const
{
    if( coordinate_axis == 0U )
    {
        assert( i<px.size() );
        return px[i];
    }
    if( coordinate_axis == 1U )
    {
        assert( i<py.size() );
        return py[i];
    }
    if( coordinate_axis == 2U )
    {
        assert( i<pz.size() );
        return pz[i];
    }

  throw overflow_error( "VData::P: coordinate axis is out of range");
}


size_t VData::Vertices() const { return px.size(); }
    
size_t VData::Elements() const { return plist.size() - (plist.size() - first_face_); }

size_t VData::Faces() const { return first_interface_ - first_face_; }

size_t VData::Interfaces() const { return plist.size() - first_interface_; }

size_t VData::NodeManifolds() const { return pmanifolds_.size(); }

/// the plist contains all: elements, faces and interfaces
size_t VData::TotalNumberOfCells() const { return plist.size(); }

size_t VData::ElementTypes() const {
//     return pelmt.size();
     return set<int8_t>( pelmt.begin(), pelmt.end() ).size();
  }

size_t VData::ElementNeighbors() const { return pfverts.size(); }

size_t VData::BFlags() const {
   set<int8_t> range_of_bflags( bflags.begin(), bflags.end() );
   return range_of_bflags.size();
}

bool VData::HybridElementTypeMesh() const { return hybrid_mesh_; }

void VData::HybridElementTypeMesh( bool hybrid_mesh ) { hybrid_mesh_ = hybrid_mesh; }

void VData::AddElementTypes( vector<int8_t>::const_iterator first,
                             vector<int8_t>::const_iterator last )
 { pelmt.assign( first, last ); }

  void VData::AddElementTypes( deque<int8_t>::const_iterator first,
                               deque<int8_t>::const_iterator last )
  { pelmt.assign( first, last ); }
  


void VData::AddNodeManifolds( VData::manifoldContainer::const_iterator first,
                              VData::manifoldContainer::const_iterator last )
 {
    assert( Interfaces() > 0 );
    assert( distance(first,last) < Vertices() );
    pmanifolds_.assign( first, last );
 }
                          



    /// reports whether the model contains only isoparametric element types
bool  VData::IsoparametricElementMesh() const
 {
     for ( auto it : pelmt )
       if ( !CSMP_ElementSpecifications::IsIsoparametric(it) )
         return false;
     return true;
 }



/// using element types, coordinate range, and boundary flags, asesses whether this is a 1D, 2D , or three dimensional model
int VData::SpatialDimension() const
 {
   // 1. looking at the element types
   if ( !HybridElementTypeMesh() )
     return CSMP_ElementSpecifications::MinimumSpatialDimension( pelmt[0] );
     
   uint32_t spatial_dim{1};
   for ( auto i : pelmt ) {
        spatial_dim = max( spatial_dim, CSMP_ElementSpecifications::MinimumSpatialDimension(i) );
        if (  spatial_dim == 3 ) break;
     }
     
   // 2. taking the Z-coordinate range as an additional criterion
   const pair<double,double>  Z_range = Z_Range();
   if ( approximatelyEqual(Z_range.first,Z_range.second) && spatial_dim == 2 ) return 2;
   
   if ( spatial_dim == 3 && approximatelyEqual(Z_range.first,Z_range.second) )
     throw csmp::Exception( ERROR, "VData::SpatialDimension",
                           "while mesh contains volume elements, Z-coordinate range is zero");
   return spatial_dim;
   
 } // end SpatialDimension




//  aelement ID's 0...n-1              
vector<int64_t>::iterator  VData::PlistBegin( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw overflow_error( "VData::PlistBegin: input Element ID out of range");
#endif
     return plist[eidx].begin();
  } 
  
vector<int64_t>::iterator  VData::PlistEnd( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw overflow_error( "VData::PlistEnd: input Element ID out of range");
#endif
     return plist[eidx].end();
  }  
     
vector<int64_t>::iterator  VData::PfvertsBegin( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw overflow_error( "VData::PfvertsBegin: input Element ID out of range");
#endif
     return pfverts[eidx].begin();
  } 
  
vector<int64_t>::iterator  VData::PfvertsEnd( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw overflow_error( "VData::PfvertsEnd: input Element ID out of range");
#endif
     return pfverts[eidx].end();
  }

// constant versions
vector<int64_t>::const_iterator  VData::PlistBegin( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw overflow_error( "VData::PlistBegin: input Element ID out of range");
#endif
     return plist[eidx].begin();
  } 
  
vector<int64_t>::const_iterator  VData::PlistEnd( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw overflow_error( "VData::PlistEnd: input Element ID out of range");
#endif
     return plist[eidx].end();
  }  
     
vector<int64_t>::const_iterator  VData::PfvertsBegin( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > pfverts.size() )
       throw overflow_error( "VData::PfvertsBegin: input Element ID out of range");
#endif
     return pfverts[eidx].begin();
  } 
  
vector<int64_t>::const_iterator  VData::PfvertsEnd( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > pfverts.size() )
       throw overflow_error( "VData::PfvertsEnd: input Element ID out of range");
#endif
     return pfverts[eidx].end();
  }


void   VData::ElementType( size_t eidx, int8_t etype )
 {
    assert( eidx < pelmt.size() );
    pelmt[ eidx ] = etype;
 }

int8_t  VData::ElementType( size_t eidx ) const
 {
    assert( eidx < pelmt.size() );
    return pelmt[ eidx ];
 }




// Accessors (iterators)

deque<vector<int64_t> >::iterator VData::PlistBegin()
 { return plist.begin(); }

deque<vector<int64_t> >::iterator VData::PlistEnd()
 { return plist.end(); }

deque<vector<int64_t> >::iterator VData::PfvertsBegin()
 { return pfverts.begin(); }

deque<vector<int64_t> >::iterator VData::PfvertsEnd()
 { return pfverts.end(); }

vector<int8_t>::iterator VData::BFlagsBegin()
 { return bflags.begin(); }

vector<int8_t>::iterator VData::BFlagsEnd()
 { return bflags.end(); }


// uses actual element IDs (1...n)
uint32_t  VData::PlistSize( size_t eidx ) const
 { return static_cast<uint32_t>(plist[eidx].size()); }


// const iterators
vector<int8_t>::const_iterator  VData::PelmtBegin() const
 { return pelmt.begin(); }
    
vector<int8_t>::const_iterator  VData::PelmtEnd() const
 { return pelmt.end(); }

deque<vector<int64_t> >::const_iterator VData::PlistBegin() const
 { return plist.begin(); }

deque<vector<int64_t> >::const_iterator VData::PlistEnd() const
 { return plist.end(); }

uint32_t  VData::PfvertsSize( size_t eidx ) const
 { return static_cast<uint32_t>(pfverts[eidx].size()); }

deque<vector<int64_t> >::const_iterator VData::PfvertsBegin() const
 { return pfverts.begin(); }

deque<vector<int64_t> >::const_iterator VData::PfvertsEnd() const
 { return pfverts.end(); }

vector<int8_t>::const_iterator VData::BFlagsBegin() const
 { return bflags.begin(); }

vector<int8_t>::const_iterator VData::BFlagsEnd() const
 { return bflags.end(); }



vector<int8_t>::iterator VData::BREP_FlagsBegin()
 { return gflags_.begin(); }

vector<int8_t>::iterator VData::BREP_FlagsEnd()
 { return gflags_.end(); }

vector<int8_t>::const_iterator VData::BREP_FlagsBegin() const
 { return gflags_.begin(); }

vector<int8_t>::const_iterator VData::BREP_FlagsEnd() const
 { return gflags_.end(); }



void VData::BFlag( size_t node_id, int8_t bflag )
 {
    assert( node_id < bflags.size() );
    BOX_BOUNDARY flag = static_cast<BOX_BOUNDARY>(bflag);
    // if the boundary flag integer value is outside of the range of defined values
    if ( flag < MULTIPLE_BOUNDARIES - 1 ) {
         cerr <<"\nVData::BFlag: boundary flag "<< parseBoundary( flag );
         cerr <<" is uninterpretable; no assignment was made.\n";
         return;
      }
    bflags[ node_id ] = bflag;
 }
 
 
 
     /// returns the boundary flag of the node
int8_t VData::BFlag( size_t node_id ) const
 {
    assert( node_id < bflags.size() );
    return bflags[node_id];
 }

 
 
 void VData::BREP_Flag( size_t node_id, int8_t bflag )
 {
    assert( node_id < gflags_.size() );
    TOPOTYPE flag = static_cast<TOPOTYPE>(bflag);
    // if the boundary flag integer value is outside of the range of defined values
    if ( flag < 0 ) {
         cerr <<"\nVData::BREP_Flag: BREP topology / geometry entity flag "<< parseTopology( flag );
         cerr <<" is uninterpretable; no assignment was made.\n";
         return;
      }
    gflags_[ node_id ] = bflag;
 }
 
 
 
     /// returns the boundary flag of the node
int8_t VData::BREP_Flag( size_t node_id ) const
 {
    assert( node_id < gflags_.size() );
    return gflags_[node_id];
 }

 
 
/**
      Returns boundary identifier if any.
*/
int8_t VData::BoundaryFlag( size_t vertex ) const
 {
    if ( bflags.empty() ) {
         cerr <<"\nVData::ABoundaryFlag: cannot determine boundary flag because VData contains no boundary identifiers.\n";
         return NOT;
      }
    else if ( vertex >= bflags.size() ) {
         cerr <<"\nVData::ABoundaryFlag: input 'vertex' is out of range.\n";
         return NOT;
      }
    return bflags[vertex];
 }

 
 
 
 
// FACE AND INTERFACE-RELATED ITERATORS

// specific element, face and interface iterators
/// iterator to CSMP finite element type of first face stored in mesh
vector<int8_t>::const_iterator  VData::PelmtFacesBegin() const {
    return next( pelmt.begin(), first_face_ );
 }
 
/// iterator to CSMP finite element type of first interface stored in mesh
vector<int8_t>::const_iterator  VData::PelmtInterfacesBegin() const {
    return next( pelmt.begin(), first_interface_ );
 }

deque<vector<int64_t> >::const_iterator  VData::PlistElmtsBegin() const {
    return plist.begin();
 }

deque<vector<int64_t> >::const_iterator  VData::PlistElmtsEnd() const {
    return next( plist.begin(), first_face_ );
 }

deque<vector<int64_t> >::const_iterator  VData::PlistFacesBegin() const {
    return next( plist.begin(), first_face_ );
 }

deque<vector<int64_t> >::const_iterator  VData::PlistFacesEnd() const {
    return next( plist.begin(), first_interface_ );
 }

deque<vector<int64_t> >::const_iterator  VData::PlistInterFacesBegin() const {
    return next( plist.begin(), first_interface_ );
 }

deque<vector<int64_t> >::const_iterator  VData::PlistInterFacesEnd() const {
    return plist.end();
 }

/// neighbor iterator for first face in plist (use PlistInterFacesBegin() to find last one)
deque<vector<int64_t> >::const_iterator  VData::PfvertsFacesBegin() const {
    return next( pfverts.begin(), first_face_ );
 }
 
/// neighbor iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
deque<vector<int64_t> >::const_iterator  VData::PfvertsInterfaceBegin() const {
    return next( pfverts.begin(), first_interface_ );
 }
 
 
 

/// read / write access to the stored node manifolds
VData::manifoldContainer::iterator  VData::PmanifoldsBegin()
 {
    return pmanifolds_.begin();
 }
 
VData::manifoldContainer::iterator  VData::PmanifoldsEnd()
 {
    return pmanifolds_.end();
 }

VData::manifoldContainer::const_iterator  VData::PmanifoldsBegin() const
 {
    return pmanifolds_.begin();
 }
 
VData::manifoldContainer::const_iterator  VData::PmanifoldsEnd() const
 {
    return pmanifolds_.end();
 }

 
 
 

void VData::Plist( size_t eidx, uint32_t nidx, size_t val )
  { 
     assert( eidx < plist.size() ); 
     assert( nidx < plist[eidx].size() ); 
     plist[eidx][nidx] = val; 
  }

int64_t  VData::Plist( size_t eidx, uint32_t nidx ) const
  { 
     assert( eidx < plist.size() );
    if (nidx >= plist[eidx].size()) {
      cerr << "plist[" << eidx << ".size() = " <<plist[eidx].size() << '\n';
    }
     assert( nidx < plist[eidx].size() ); 
     return plist[eidx][nidx]; 
  }

/// nidx is neighbour index
void VData::Pfvert( size_t eidx, uint32_t nidx, int64_t  val )
  { 
     assert( eidx < pfverts.size() ); 
     assert( nidx < pfverts[eidx].size() ); 
     pfverts[eidx][nidx] = val; 
  }

/// nidx is neighbour index
int64_t  VData::Pfvert( size_t eidx, uint32_t nidx ) const
  { 
     assert( eidx < pfverts.size() ); 
     assert( nidx < pfverts[eidx].size() ); 
     return pfverts[eidx][nidx]; 
  }




/**
    Sets the element type of a mono-type element mesh
*/
void VData::SingleElementType( int8_t etype )
 {
    if ( pelmt.size() > 1U )
      cout <<"\nVData::SingleElementType: changing VSet from hybrid element to single-element container."<< endl;
    pelmt.clear();
    pelmt.reserve(1U);
    pelmt.push_back( etype );
    hybrid_mesh_ = false;
    // since faces are lower dimensional, at least 2 etypes would be required
    first_face_      = static_cast<int64_t>(plist.size());
    first_interface_ = static_cast<int64_t>(plist.size());
 }
    
    
    

/**
    Re-initialises the element-type vector that is associated
    with the VSet.
    
    @note the element type vector contains either:
    
    1. a single element (mono-element type mesh)
 
    2. element type entries for each element of the mesh
       including potential faces and interfaces
*/
void VData::ElementTypes( const vector<int8_t>& elmt_types )
 {
    assert( !elmt_types.empty() );
    pelmt.clear();
    pelmt.reserve(elmt_types.size());
    pelmt.assign( elmt_types.begin(), elmt_types.end() );
    if ( elmt_types.size() > 1 ) hybrid_mesh_ = true;
 }



/**
 
Resize() deletes all the contents of the curren VData and reinitializes
the connectivity storage for the single element type specified as 
argument. 

@section arguments Input Arguments

A reference to the FiniteElement type for which connectivity storage shall
be created and the quantities of nodes, constraint points, and elements
in the mesh that shall be stored. 

@section application Application 

Resize() allows a VSet to be re-used in a computation. 

@attention method clears the VSet: any previously stored information is lost.

 */
void VData::Resize( size_t nodes_per_element, 
                    size_t nbors_per_element, 
                    int8_t etype,
                    size_t nodes, size_t elmts )
 {
    px.resize(nodes);  vector<double>( px ).swap( px );
    py.resize(nodes);  vector<double>( py ).swap( py );
    pz.resize(nodes);  vector<double>( pz ).swap( pz );
    ResizeBFlags();
    ResizeBREP_Flags();
    
    pelmt.clear();
    plist.clear();
    pfverts.clear();

    for ( size_t i{0U}; i<elmts; ++i ) {
         plist.push_back( vector<int64_t>(nodes_per_element) );
      }
    for ( size_t k{0U}; k<elmts; ++k )
      pfverts.push_back( vector<int64_t>(nbors_per_element) );
      
    hybrid_mesh_ = false;
    
    SingleElementType( etype );

    // single element meshes have no faces nor interfaces
    first_face_      = static_cast<int64_t>(plist.size());
    first_interface_ = static_cast<int64_t>(plist.size());
 }



/**
 
Resize() deletes all the contents of the curren VData and reinitializes
the connectivity storage for a mixed element mesh with coordination 
information as stored in the argument vectors.  

@section arguments Input Arguments

Resize() uses its deque arguments to determine how many nodes and
neighbor elements will be connected to each element.

The total quantities of nodes, elements, faces and interfaces
that shall be stored in the mesh are input as last three arguments. 

@attention the nodes-per-element and neighbors-per-element deques
must already contain the correct number of neighbors for the 
faces (=nbors+2 higher-dim elements) and interfaces (nbors*2 + 2)!

@section application Application 

Resize() allows a VSet to be re-used in a computation.

@attention regarding faces and interfaces, these are assumed 
to be included into the supplied deques.
 
*/
  void VData::Resize( const deque<int8_t>& etypes,
                      const deque<uint32_t>& npes,
                      const deque<uint32_t>& epes,
                      size_t nodes,
                      size_t faces,
                      size_t interfaces )
  {
    assert(npes.size() == epes.size());
    const size_t nrCells = npes.size();
    assert( faces + interfaces < nrCells );
    
    px.resize(nodes);  vector<double>( px ).swap( px );
    py.resize(nodes);  vector<double>( py ).swap( py );
    pz.resize(nodes);  vector<double>( pz ).swap( pz );
    ResizeBFlags();
    ResizeBREP_Flags();
    
    plist.clear();
    pfverts.clear();
    pelmt.clear();
    
    pelmt.assign( etypes.begin(), etypes.end() );
    
    for ( auto i{0U}; i<nrCells; i++ )
      plist.emplace_back( vector<int64_t>(npes[i],UINT_MAX) );
    
    for ( auto i{0U}; i<nrCells; i++ )
      pfverts.emplace_back( vector<int64_t>(epes[i],IRREGULAR) );
    
    if ( etypes.size() > 1U ) hybrid_mesh_ = true;
    else hybrid_mesh_ = false;
    
    // faces
    first_face_ = static_cast<int64_t>(nrCells - faces - interfaces);
    
    // interfaces
    first_interface_ = static_cast<int64_t>(nrCells - interfaces);
  }






/**
 
Changes the size of the storage arrays for the x,y and z coordinates of
the nodes to the indicated size. 

@param nodes The new number of node coordinates stored in the Vdata.

@section implementation Implementation 

Calls Resize() on its container class. Thus, all previously stored 
coordinate numbers are lost. 

@section application Application 

To change the order of the elements in parent VSets. For instance, to 
turn a linear triangle mesh in a quadratic triangle mesh. 

@section messages Messages 

The method tests whether the new node numbers are compatible with the
node indentification numbers in the 'plist' array. If not, a warning 
is issued. 
*/
void VData::ResizeNodes( size_t nodes )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    size_t  node_max(0U);
    
    // getting highest node ID in the plist array
    if ( !plist.empty() ) {
       for (  deque<vector<int64_t> >::const_iterator
              it=plist.begin(); it!=plist.end(); it++ )
         for ( vector<int64_t>::const_iterator
               i=(*it).begin(); i!=(*it).end(); i++ )
           if ( (*i) > node_max ) node_max = (*i);
        
       if ( nodes < node_max )
         csmp_error.Note( WARNING, "VData::ResizeNodes", "'plist' contains node numbers larger than desired size" );
      }
      
    if ( px.size() > 0U )
      csmp_error.Note( WARNING, "VData::ResizeNodes", "erasing node coordinates in 'vdata'" );
    px.resize( nodes );  vector<double>( px ).swap( px );
    py.resize( nodes );  vector<double>( py ).swap( py );
    pz.resize( nodes );  vector<double>( pz ).swap( pz );
    ResizeBFlags();
    ResizeBREP_Flags();
        
 } // end ResizeNodes




/**
    Updates length of vector to that of the node vector, setting potential new flags to zero.
    preserving the original values if the size of the VSet increases.
*/
void VData::ResizeBFlags()
 {
    const size_t orig_size{ bflags.size() };
    bflags.resize( px.size() );
    if ( bflags.size() > orig_size )
    fill( next(bflags.begin(),orig_size), bflags.end(), NOT );
    bflags.shrink_to_fit();
 }


/**
     If the record grows in size, the original values are preserved and the new ones are initialised to MESH_VERTEX.
*/
void VData::ResizeBREP_Flags()
 {
    const size_t orig_size{ gflags_.size() };
    gflags_.resize( px.size() );
    if ( gflags_.size() > orig_size )
    fill( next(gflags_.begin(),orig_size), gflags_.end(), MESH_VERTEX );
    gflags_.shrink_to_fit();
 }


/**
     Resizes 'plist' so that it can contain a different number of
     elements.
     
     @note assumes that the number of faces and interfaces is not 
     affected by the change in the number of elements.
*/
void VData::ResizePlist( size_t elements )
 {
    const size_t old_n_elements(plist.size());
  
    vector<int64_t>  empty_vec;
    const size_t old_size( plist.size() );
    assert( elements > old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i{0U}; i<new_elmts; i++ )
      plist.push_back( empty_vec );
   
    // if there were no faces or interfaces initially
    if ( first_face_ == old_n_elements or first_interface_ == old_n_elements ) {
         if ( first_face_ == old_n_elements )
           first_face_  = plist.size();
         if ( first_interface_ == old_n_elements )
           first_interface_ = plist.size();
         return;
      }

    // if there were faces or interfaces respective offsets are adjusted
    const size_t change_in_n_elements(plist.size() - old_n_elements);
   
    first_face_      += change_in_n_elements;
    first_interface_ += change_in_n_elements;
 }






/**
    increasing the storage for element types to the new size without invalidating existing types unless the storage is shrunk
*/
void VData::ResizeElementTypes( size_t elements )
 {
    if ( !HybridElementTypeMesh() ) return;
    pelmt.resize( elements );
    
 } // end ResizeElementTypes
 
 
 


/**
     Resizes 'plist' so that it can contain a different number of
     elements and neighbors per element.
     
     @attention use only for meshes with a single element type.
     
     @note assumes that the number of faces and interfaces is not 
     affected by the change in the number of elements.
*/
void VData::ResizePlist( size_t elements, uint32_t nperelmt )
{
    assert( Faces() == 0 );
    assert( Interfaces() == 0 );
    vector<int64_t>  empty_vec(nperelmt,0);
    const size_t old_size( plist.size() );
    assert( elements > old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i{0U}; i<new_elmts; i++ )
        plist.push_back( empty_vec );

    // mono-element meshes can have no faces nor interfaces
    first_face_      = elements;
    first_interface_ = elements;
}



/**
    For mono element-type mesh change number of nodes per element.
*/
void VData::ResizeElementNodes( size_t eid, uint32_t nperelmt )
{
    assert( !hybrid_mesh_ );
    plist[eid].resize( nperelmt, 0 );
}



/**
    Resizes the plist without touching the number of faces or interfaces.
*/
void VData::ResizePlist( const deque<uint32_t>& mixed_ele_plist )
 {
    const size_t  old_n_elements(plist.size());

    // plist: vector<vector<uint32_t> >
    if ( !plist.empty() ) {
         cout <<"\nVData::ResizePlist: Warning: Erasing existing 'plist'"<< endl;
         plist.erase( plist.begin(), plist.end() );
      }

    hybrid_mesh_ = false;

    vector<int64_t>  empty_vec;
    size_t          n_last(mixed_ele_plist[0]);

    for ( auto i{0U}; i<mixed_ele_plist.size(); ++i )
      {
         if ( !hybrid_mesh_  and  n_last != mixed_ele_plist[i] ) hybrid_mesh_ = true;
         // extending member vectors in place (avoiding copying)
         plist.push_back( empty_vec );
         plist[i].reserve( (n_last=mixed_ele_plist[i]) );
         for ( size_t j{0U}; j<n_last; ++j ) plist[i].push_back( 0 );
      }

    // if there were no faces or interfaces initially
    if ( first_face_ == old_n_elements or first_interface_ == old_n_elements ) {
         if ( first_face_ == old_n_elements )
           first_face_      = static_cast<int64_t>(plist.size());
         if ( first_interface_ == old_n_elements )
           first_interface_ = static_cast<int64_t>(plist.size());
         return;
      }

    // the offsets to first face  and interface are adjusted
    const size_t change_in_n_elements(plist.size() - old_n_elements);
    first_face_      += change_in_n_elements;
    first_interface_ += change_in_n_elements;

 } // end ResizePlist




/**
    increases the length of the 'pfverts' record
*/
void VData::ResizePfverts( size_t elements )
{
    if ( first_face_  != plist.size() || first_interface_ != plist.size() )
      throw logic_error( "VData::ResizePfverts(size_t): method does not handle VSets with Face and InterFace objects yet");

    vector<int64_t>  empty_vec;
    const size_t old_size( pfverts.size() );
    assert( elements >= old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i{0U}; i<new_elmts; i++ )
        pfverts.push_back( empty_vec );
}


/**
    increases the number of neighbor elements in a single element-type mesh
    
        @attention this works only of the size of the pfverts record is increased
 */
void VData::ResizePfverts( size_t elements, uint32_t nperelmt )
{
    if ( first_face_  != plist.size() || first_interface_ != plist.size() )
      throw logic_error( "VData::ResizePfverts(size_t,size_t): method does not handle VSets with Face and InterFace objects yet");

    vector<int64_t>  empty_vec(nperelmt,0);
    const size_t old_size( pfverts.size() );
    assert( elements >= old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i{0U}; i<new_elmts; i++ )
        pfverts.push_back( empty_vec );
}



/**
     change the size of the storage for the number of neighbors of target element
*/
void VData::ResizeElementNeighbors( size_t eid, uint32_t nperelmt )
{
    pfverts[eid].resize( nperelmt, 0 );
}




/**
   just the size is changed but no data are transferred
*/
void VData::ResizePfverts( const deque<uint32_t>& mixed_ele_pfverts )
 {
    vector<int64_t>  empty_vec;
    const int64_t    flag_other(0);

    // pfverts: vector<vector<uint32_t> >
    if ( !pfverts.empty() ) {
         cout <<"\nVData::ResizePfverts: Warning: Erasing existing 'pfverts'"<< endl;
         pfverts.erase( pfverts.begin(), pfverts.end() );
      }

    for ( size_t i{0U}; i<mixed_ele_pfverts.size(); ++i )
      {
         // extending member vectors in place (avoiding copying)
         pfverts.push_back( empty_vec );
         pfverts[i].reserve( mixed_ele_pfverts[i] );
         for ( size_t j{0U}; j<mixed_ele_pfverts[i]; j++ ) pfverts[i].push_back( flag_other );
      }

 } // end ResizePfverts


/**
    Empties the map that contains the which nodes lie at the model boundary
    and the BOX_BOUNDARY flags that these nodes had.
*/
void VData::RemoveBflags()
 {
    bflags.erase( bflags.begin(), bflags.end() );
 }



/**
 
Checkthe data arrays inside the VData for whether they contain any
data that are out of range (node numbers which are greater than the size
of the containing array or zero). The datasets which are checked are: 

px, py, pz (Arrays)
plist, pfverts (maps of vectors)      
bflags (map)
bvals (map)

Non-existant boundary identifier flags (negative integers) in the 'pfverts'
array are replaced by the REGION_BOUNDARY (-28) flag.

@return If no violations are detected the method returns true, else false.

@section application Application 

VSets and the contained VData represent the interface to external 
programs. Thus, before erratic data are used to build Model objects
it is useful to check the input data.  

@section messages Messages 

Check() tries to be as specific as possible in reporting where and what 
kind of error occurred. Warnings are issued when datasets that are not
essential for the building of a mesh are corrupted. 
*/
bool VData::CheckFix()
 {
   // give it the benefit of the doubt
   bool ok(true);

   cout <<"\nVData::Check: checking 'vdata'..."<< endl;

   // plist
   // -----
   for ( const auto& pit : plist )
     for ( const auto& it : pit )
       if ( it > Vertices() )
         {
            cerr <<"\nVData::Check: Wrong node number (index) in plist: ";
            cerr << it <<" versus n-nodes: "<< Vertices() << endl;
            ok = false;
         }

   // pfverts
   // -------
   for ( auto& fpit : pfverts )
     for ( auto& fit : fpit ) {
       if ( fit > static_cast<int64_t>(Elements()) )
         {
            cerr <<"\nVData::Check: Non-existant high element number (index) in pfverts: ";
            cerr << fit <<" versus n-elements: "<< Elements() << endl;
            ok = false;
         }
       if ( fit < REGION_BOUNDARY )
         {
            if ( ok ) {
                 cerr <<"\n\n\n\nVData::Check: Too low boundary flag number (index) in pfverts: ";
                 cerr << fit <<" versus n-elements: "<< Elements() << endl << endl;
              }
            fit = REGION_BOUNDARY;
            ok = false;
         }
       }  

   // boundary flags
   // --------------
   if ( bflags.empty() ) 
     cout <<"\nVData::Check: No boundary nodes could be identified..."<< endl;
   else {
	    size_t  counter(0U);
	    for ( auto bf=bflags.begin(); bf!=bflags.end(); bf++ )
	      if ( (*bf) < REGION_BOUNDARY ) {
	            (*bf) = REGION_BOUNDARY;
	            counter++;
	         }
        if ( counter > 0 ) 
	      cout <<"\nVData::CheckFix: Replaced "<< counter <<" 'bflags' by REGION_BOUNDARY flag.\n";
     }

   return ok;

 } // end Check




/**
     If this is not already the case, this method will
     convert the vdata such that nodes, elements, and neighbors
     are numbered from 0..n-1.
*/
void VData::EstablishZeroBasedNumbering()
 {
    // test whether nodes are not already numbered 0..1
    for ( auto& it : plist )
      for ( auto& n : it )
        if ( n == 0U ) {
             cerr << "\nVData::EstablishZeroBasedNumbering: Numbering already is 0..n-1 based.\n";
             return;
          }
    // to convert: 
    // plist
    for ( auto& it : plist )
      for ( auto& n : it )
        n -= 1U;

    // pfverts
    for ( auto& it : pfverts )
      for ( auto& n : it )
        // only the neighbor element ids, not the boundary flags must be decremented
        if ( n > 0 ) n -= 1;

    // bflags & bconds - nothing needs to be done
    
 } // end EstablishZeroBasedNumbering







/** 
    Writes VData into a binary output file.
*/ 
void VData::OutBinary( fstream& fp ) const
 {
    static_assert( sizeof(size_t) == sizeof(streamsize), "VData::OutBinary: on this system 'streamsize' is not equal to size_t" );
    size_t n0(0U), n1(1U), records;
    double* ptr(0);
    
    // 1. writing whether we are dealing with a mixed mesh
    // ---------------------------------------------------
    if ( hybrid_mesh_ ) fp.write( reinterpret_cast<const char*>(&n1), sizeof(size_t) );
    else                fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
    
    // 2. writing all the p,c arrays or length identifiers = 0
    // -------------------------------------------------------
     {
       BinaryFileSectionWrite sect(fp, "VSETCORD");
       // px
      if ( (records=px.size()) > 0 && (ptr=const_cast<double*>( &(*px.begin()) )) != nullptr )
        {
           fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
           fp.write( reinterpret_cast<const char*>(ptr), sizeof(double) * records );
        }
      else fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
      // py
      if ( (records=py.size()) > 0 && (ptr=const_cast<double*>( &(*py.begin()) )) != nullptr )
        {
           fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
           fp.write( reinterpret_cast<const char*>(ptr), sizeof(double) * records );
        }
      else fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
      // pz
      if ( (records=pz.size()) > 0 && (ptr=const_cast<double*>( &(*pz.begin()) )) != nullptr )
        {
           fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
           fp.write( reinterpret_cast<const char*>(ptr), sizeof(double) * records );
        }
      else fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
     }

    // 3. writing pelmt, plist, pfverts, bflags, gflags
    // ------------------------------------------------
     {
       BinaryFileSectionWrite sect(fp, "VSETPELT");
       binaryFileWrite( fp, pelmt );
     }
     {
       BinaryFileSectionWrite sect(fp, "VSETPLST");
       binaryFileWrite( fp, plist );
     }
     {
       BinaryFileSectionWrite sect(fp, "VSETPFVT");
       binaryFileWrite( fp, pfverts );
     }
     {
       BinaryFileSectionWrite sect(fp, "VSETBFLG");
       binaryFileWrite( fp, bflags );
     }
     {
       BinaryFileSectionWrite sect(fp, "VSETGFLG");
       binaryFileWrite( fp, gflags_ );
     }
   
     // writing the manifold classifiers first
     if ( (records=pmanifolds_.size()) > 0 )
        {
           fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
           
           // splitting 'pmanifold' into a vector of ManifoldTypes and a deque of node-number vectors
           // ---------------------------------------------------------------------------------------
           // (vector< pair< vector<int64_t>, ManifoldType > > pmanifolds_)
           deque<vector<size_t> > manifold_nodes;
           vector<ManifoldType>    manifold_types; manifold_types.reserve( records );
           for ( const auto& nmf : pmanifolds_ ) {
                manifold_nodes.push_back( nmf.first );
                manifold_types.push_back( nmf.second );
             }
           // 4. node manifolds, if any
           // --------------------------
           { // node numbers
             BinaryFileSectionWrite sect(fp, "VSETNMND");
             binaryFileWrite( fp, manifold_nodes );
           }
           { // node flags
             BinaryFileSectionWrite sect(fp, "VSETNMMT");
             binaryFileWrite( fp, manifold_types );
           }
       }
      else // indicating that there are no records
        fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );

    // 5. offsets for faces and interfaces
    // -----------------------------------
    fp.write( reinterpret_cast<const char*>(&first_face_), sizeof(int64_t ) );
    fp.write( reinterpret_cast<const char*>(&first_interface_), sizeof(int64_t ) );
    
    cout <<"\nVData::OutBinary: Mesh has been successfully written to file."<< endl;

 } // end OutBinary
 
 
 
 
 
/** 
     Initialises VData from binary input file.
*/
void VData::InBinary( fstream& fp )
 {
    static_assert( sizeof(size_t) == sizeof(streamsize), "VData::InBinary: on this system 'streamsize' is not equal to size_t" );
    size_t mixed(0U), records(0U);
    
    // 1. reading whether we are dealing with a mixed mesh
    // ---------------------------------------------------
    fp.read( reinterpret_cast<char*>(&mixed), sizeof(size_t) );
    if ( mixed ) hybrid_mesh_ = true;
    else         hybrid_mesh_ = false;
    
    // 2. reading all the p,c arrays or length identifiers = 0
    // -------------------------------------------------------
     {
       BinaryFileSectionRead sect(fp, "VSETCORD");

      // px
      fp.read( reinterpret_cast<char*>(&records), sizeof(size_t));
      if ( records > 0U ) {
           px.resize( records );
           vector<double>( px ).swap( px );
           fp.read( reinterpret_cast<char*>(&(*px.begin())), sizeof(double) * records );
        }

      // py
      fp.read( (char*) &records, sizeof(size_t));
      if ( records > 0U ) {
           py.resize( records );
           vector<double>( py ).swap( py );
           fp.read( reinterpret_cast<char*>(&(*py.begin())), sizeof(double) * records );
        }

      // pz
      fp.read( (char*) &records, sizeof(size_t));
      if ( records > 0U ) {
           pz.resize( records );
           vector<double>( pz ).swap( pz );
           fp.read( reinterpret_cast<char*>(&(*pz.begin())), sizeof(double) * records );
        }
     }

    // 3. reading pelmt, plist, pfverts, bflags
    // ----------------------------------------
     {
        BinaryFileSectionRead sect(fp, "VSETPELT");
        binaryFileRead( fp, pelmt );
     }
     {
        BinaryFileSectionRead sect(fp, "VSETPLST");
        binaryFileRead( fp, plist );
     }
     {
        BinaryFileSectionRead sect(fp, "VSETPFVT");
        binaryFileRead( fp, pfverts );
     }
     {
        BinaryFileSectionRead sect(fp, "VSETBFLG");
        binaryFileRead( fp, bflags );
     }
     {
       BinaryFileSectionRead sect(fp, "VSETGFLG");
       binaryFileRead( fp, gflags_ );
     }
   
     // 4. node manifolds, if any
     // --------------------------
     fp.read( reinterpret_cast<char*>(&records), sizeof(size_t) );
     if ( records > 0 )
       {
          deque<vector<size_t> > manifold_nodes;
          vector<ManifoldType>   manifold_types;
          {
            BinaryFileSectionRead sect(fp, "VSETNMND");
            binaryFileRead( fp, manifold_nodes );
          }
          {
            BinaryFileSectionRead sect(fp, "VSETNMMT");
            binaryFileRead( fp, manifold_types );
          }
          // initialising the manifold container in VData
          if ( !pmanifolds_.empty() ) pmanifolds_.clear();
          pmanifolds_.reserve( manifold_nodes.size() );
          for ( size_t i{0U}; i<manifold_nodes.size(); i++ )
            pmanifolds_.push_back( make_pair( manifold_nodes[i], manifold_types[i] ) );
       }

    // 5. offsets for faces and interfaces
    // -----------------------------------
    fp.read( reinterpret_cast<char*>(&first_face_), sizeof(int64_t ) );
    fp.read( reinterpret_cast<char*>(&first_interface_), sizeof(int64_t ) );
    
    cout <<"\nVData::InBinary: Mesh has been successfully read from file."<< endl;
        
 } // end InBinary

 
 
 

/**
    @note this just prints the vdata not the entire VSet
*/
void VData::OutASCII( const char* file ) const
  {
     string  file_name(file); file_name += "-vdata.txt";
     ofstream  ofs( file_name.c_str(), ios::out );
     if ( !ofs.is_open() ) {
           string  error("VData::OutASCII");
           error += "unable to open ASCII file with name : ";
           error += file_name;
           throw invalid_argument(error.c_str());
       }
     ofs <<"\nVData::Out: "<< endl;

     // mixed mesh
     // ----------
     if ( hybrid_mesh_ ) ofs <<"\tmesh contains different type of elements..."<< endl;
     else ofs <<"\tmesh contains only one type of element..."<< endl;
    
     // faces or interfaces
     // -------------------
     if ( first_face_ < plist.size() and first_face_ < first_interface_  )
       ofs <<"\tmesh contains "<< first_interface_ - first_face_ <<" descriptors of Face objects."<< endl;
     if ( first_interface_ < plist.size() )
       ofs <<"\tmesh contains "<< plist.size() - first_interface_ <<" descriptors of InterFace objects."<< endl;

     // px, py, pz
     // ----------
     assert( px.size() == py.size() );
     assert( py.size() == pz.size() );
     ofs <<"\n'px, py, pz' coordinates of "<< px.size() <<" nodes:"<< endl;
     for ( auto i{0U}; i<px.size(); i++ )
       ofs << i <<": \t"<< px[i] <<"\t"<< py[i] <<"\t"<< pz[i] << endl;

     // pelmt
     // -----
     ofs <<"\n'pelmt' finite element types:";
     size_t i(0U);
     for ( vector<int8_t>::const_iterator
           eit=pelmt.begin(); eit!=pelmt.end(); eit++, i++ )
       ofs <<"\n"<< i <<" = "<< static_cast<int>(*eit) <<" = CSMP type: "<< parseFiniteElementType( static_cast<CSMP_FEM_TYPE>(*eit) );
     ofs << endl;

     // plist
     // -----
     ofs <<"\n'plist' entries of "<< plist.size() <<" elements:"<< endl;
     i = 0U;
     for ( deque<vector<int64_t> >::const_iterator
           pt=plist.begin(); pt!=plist.end(); pt++, i++ )
       {
          ofs << i <<": \t";
          for ( size_t j{0U}; j<(*pt).size(); j++ ) ofs << (*pt)[j] <<"\t";
          ofs << endl;
       }

     // pfverts
     // -------
     ofs <<"\n'pfverts':"<< endl;
     i = 0U;
     for ( deque<vector<int64_t> >::const_iterator
           ft=pfverts.begin(); ft!=pfverts.end(); ft++, i++ )
       {
          ofs << i <<": \t";
          for ( size_t j{0U}; j<(*ft).size(); j++ ) ofs << (*ft)[j] <<"\t";
          ofs << endl;
       }

     // bflags
     // ------
     if ( !bflags.empty() ) ofs <<"\nBoundary flags 'bflags':"<< endl;
     size_t n_node(0U);
     for ( auto bf=bflags.begin(); bf!=bflags.end(); bf++ )
       ofs << n_node++ <<": \t"<< static_cast<int>(*bf) << endl;

     // gflags
     // ------
     if ( !gflags_.empty() ) ofs <<"\nBREP geometry flags 'gflags_':"<< endl;
     n_node = 0U;
     for ( auto bf=gflags_.begin(); bf!=gflags_.end(); bf++ )
       ofs << n_node++ <<": \t"<< static_cast<int>(*bf) << endl;

     // node manifolds
     // --------------
     if ( !pmanifolds_.empty() ) {
         ofs <<"\nNode manifolds (topologically collocated nodes) in model 'pmanifolds_':"<< endl;
         n_node = 0U;
         for ( size_t n{0U}; n<pmanifolds_.size(); ++n )
           {
              ofs << n_node++ <<": ";
              for ( auto j{0U}; j<pmanifolds_[i].first.size(); j++ )
                 cerr  << pmanifolds_[i].first[j] <<" ";
               cerr <<", manifold type: "<< parse( pmanifolds_[i].second );
           }
         cerr << endl << endl;
       }
     ofs.close();
     cout <<"\nVData::OutASCII: ascii file '"<< file_name <<"' written successfully."<< endl;

  } // end Out








    /// writes initialiser lists for the current VData in C++17 format
void VData::OutCPP17( ofstream& ofs ) const
 {
    //--------------------------ELEMENT TYPES
    // vector<int8_t> pelmt;
    const size_t n_elmts{pelmt.size()};
    ofs <<"\ndeque<int8_t>  etypes{";
    for ( auto i{0U}; i<n_elmts; ++i ) {
         ofs << static_cast<int>(pelmt[i]);
         if ( i < n_elmts-1 ) ofs <<",";
      }
    ofs <<"};\n";
 
    //--------------------------NUMBERS OF NODES & NEIGHBORS
    string delim = "";
    ofs <<"\ndeque<uint32_t> npes{";
    for ( auto e : plist ) {
        ofs << delim << e.size();
        delim =",";
      }
    ofs <<"};\n";
    delim = "";
    ofs <<"\ndeque<uint32_t> epes{";
    for ( auto e : pfverts ) {
        ofs << delim << e.size();
        delim =",";
      }
    ofs <<"};\n";
    
    ofs <<"\nvset.Resize( etypes, npes, epes, ";
    ofs << Vertices() <<", "<< Faces() <<", "<< Interfaces() <<" );\n";
 
  	//--------------------------NODE COORDINATES
    // px, py, pz
    size_t counter{0};
    const size_t n_nodes{Vertices()};
    ofs <<"\ndeque<double> px{";
    for ( auto i : px ) {
         ofs << i;
         if ( counter < n_nodes-1 ) ofs <<",";
         counter++;
      }
    ofs <<"};\n";
    counter = 0;
    ofs <<"\ndeque<double> py{";
    for ( auto i : py ) {
         ofs << i;
         if ( counter < n_nodes-1 ) ofs <<",";
         counter++;
      }
    ofs <<"};\n";
    counter = 0;
    ofs <<"\ndeque<double> pz{";
    for ( auto i : pz ) {
         ofs << i;
         if ( counter < n_nodes-1 ) ofs <<",";
         counter++;
      }
    ofs <<"};\n";
  	ofs <<"\n\nvset.AddXYZ( px, py, pz );";

    //------------------------------BOUNDARY FLAGS
    ofs <<"\n\nvector<int8_t>  bflags{";
    for ( size_t i{0U}; i<Vertices(); ++i ) {
         ofs << static_cast<int>(BFlag(i));
         if ( i <Vertices()-1 ) ofs <<",";
      }
    ofs <<"};\n";
    ofs <<"\n\nvset.AddBFlags( bflags.begin(), bflags.end());\n";
      
     //------------------------------GEOMETRY FLAGS
    ofs <<"\n\nvector<int8_t>  gflags{";
    for ( size_t i{0U}; i<Vertices(); ++i ) {
         ofs << static_cast<int>(BREP_Flag(i));
         if ( i <Vertices()-1 ) ofs <<",";
      }
    ofs <<"};\n";
    ofs <<"\n\nvset.AddBREP_Flags( gflags_.begin(), gflags_.end());\n";

    //------------------------------PLIST
    ofs <<"\n//'plist' nodes that make up the elements";
    ofs <<"\ndeque<vector<int64_t>> plist( "<< pelmt.size() <<" );";
    // Plist( size_t eidx, size_t node, size_t val );
    for ( size_t eidx{0}; eidx<plist.size(); ++eidx ) {
        ofs <<"\n\tplist["<< eidx <<"] = { ";
        for ( auto node{0U}; node<plist[eidx].size(); ++node ) {
             ofs << plist[eidx][node];
             if ( node < plist[eidx].size()-1 ) ofs <<", ";
          }
        ofs <<" };";
      }
    ofs <<"\n\nvset.AddPlist( plist.begin(), plist.end());\n";

    //------------------------------PFVERTS (neighbor information)
    ofs <<"\n//'pfverts' neighbors of the faces of each element";
    ofs <<"\ndeque<vector<int64_t>> pfverts( "<< pelmt.size() <<" );";
    // Pfverts( size_t eidx, size_t neighbor, size_t val );
    for ( size_t eidx{0}; eidx<pfverts.size(); ++eidx ) {
        ofs <<"\n\tpfverts["<< eidx <<"] = { ";
        for ( auto nbor{0U}; nbor<pfverts[eidx].size(); ++nbor ) {
             ofs << pfverts[eidx][nbor];
             if ( nbor < pfverts[eidx].size()-1 ) ofs <<", ";
          }
        ofs <<" };";
      }
    ofs <<"\n\nvset.AddPfverts( pfverts.begin(), pfverts.end());";

    //--------------------------------node manifolds
    // TODO: TEST this output
    if ( !pmanifolds_.empty() ) {
        ofs <<"\n//''pmanifolds_' node manifolds (topologically collocated nodes):"<< endl;
        ofs <<"\ndeque<pair<vector<int64_t>,ManifoldType> > pmanifolds( "<< pmanifolds_.size() <<" );";
        for ( size_t n{0U}; n<pmanifolds_.size(); ++n )
          {
             ofs <<"\n\tpmanifolds_["<< n <<"] = { {";
             for ( auto j{0U}; j<pmanifolds_[n].first.size(); j++ ) {
                  ofs  << pmanifolds_[n].first[j];
                  if ( j < pmanifolds_[n].first.size()-1 ) ofs <<", ";
               }
             ofs <<"}, "<< parse( pmanifolds_[n].second );
          }
        ofs <<" };";
      }

 } // end OutCPP17







/**
    Prints contents of VSet to console;
*/
void VData::Out() const
  {
     cout<<"\nVData::Out: "<< endl;

     // mixed mesh
     // ----------
     if ( hybrid_mesh_ ) cout <<"mesh contains different type of elements..."<< endl;
     else cout <<"mesh contains only one type of element..."<< endl;

     // faces or interfaces
     // -------------------
     if ( Elements() > 0  )
       cout <<"\tmesh contains "<< Elements() <<" descriptors of Element objects."<< endl;
     if ( first_face_ < plist.size() and first_face_ < first_interface_  )
       cout <<"\tmesh contains "<< first_interface_ - first_face_ <<" descriptors of Face objects."<< endl;
     if ( first_interface_ < plist.size() )
       cout <<"\tmesh contains "<< plist.size() - first_interface_ <<" descriptors of InterFace objects."<< endl;

     // px, py, pz
     // ----------
     assert( px.size() == py.size() );
     assert( py.size() == pz.size() );
     cout<<"\n'px, py, pz' coordinates of "<< px.size() <<" nodes:"<< endl;
     for ( auto i{0U}; i<px.size(); i++ )
       cout << i <<": \t"<< px[i] <<"\t"<< py[i] <<"\t"<< pz[i] << endl;

     // pelmt
     // -----
     cout <<"\n'pelmt' finite element types:";
     size_t i(0U);
     for ( vector<int8_t>::const_iterator
           eit=pelmt.begin(); eit!=pelmt.end(); eit++, i++ )
       cout <<"\n"<< i <<": "<< parseFiniteElementType( static_cast<CSMP_FEM_TYPE>(*eit) );
     cout << endl;

     // plist
     // -----
     cout <<"\n'plist' entries of "<< plist.size() <<" elements:"<< endl;
     i = 0U;
     for ( deque<vector<int64_t> >::const_iterator
           pt=plist.begin(); pt!=plist.end(); pt++, i++ )
       {
          cout << i <<": \t";
          for ( size_t j{0U}; j<(*pt).size(); j++ ) cout << (*pt)[j] <<"\t";
          cout << endl;
       }

     // pfverts
     // -------
     cout <<"\n'pfverts':"<< endl;
     i = 0U;
     for ( deque<vector<int64_t> >::const_iterator
           ft=pfverts.begin(); ft!=pfverts.end(); ft++, i++ )
       {
          cout << i <<": \t";
          for ( size_t j{0U}; j<(*ft).size(); j++ ) {
               if ( (*ft)[j] >= 0 )
                 cout << (*ft)[j] <<"\t ";
               else {
                    BOX_BOUNDARY bflag = static_cast<BOX_BOUNDARY>((*ft)[j]);
                    cout << parseBoundary( bflag ) <<"\t ";
                 }
            }
          cout << endl;
       }

     // bflags
     // ------
     if ( !bflags.empty() ) cout <<"\nBoundary flags 'bflags':"<< endl;
     size_t n_node(0U);
     for ( auto bf=bflags.begin(); bf!=bflags.end(); bf++ )
       cout << n_node++ <<": \t"<< parseBoundary( static_cast<BOX_BOUNDARY>(*bf) ) << endl;
       
     // gflags
     // ------
     if ( !bflags.empty() ) cout <<"\nGeometry classification flags 'gflags' (TOPOTYPE):"<< endl;
     n_node = 0U;
     for ( auto gf=gflags_.begin(); gf!=gflags_.end(); gf++ )
       cout << n_node++ <<": \t"<< parseTopology( static_cast<TOPOTYPE>(*gf) ) << endl;
       
     // node manifolds, if any
     // ----------------------
     // vector< pair< vector<int64_t>, ManifoldType > > pmanifolds_
     if ( !pmanifolds_.empty() ) {
          cout <<"\nNode manifolds 'pmanifolds' connecting node-matched mesh patches at SplitBoundaries:"<< endl;
          size_t manifold{0U};
          for ( const auto& nit : pmanifolds_ ) {
               cout << manifold++ <<": \t"<< parse( nit.second ) <<", connecting nodes: ";
               for ( const auto& n : nit.first )
                 cout << n <<" ";
               cout << endl;
            }
          cout << endl;
       }

  } // end Out





/**
    permits initialisation of VSet from file stream.
    
    @note no faces or interfaces are considered
*/
void VData::InText( ifstream& ifs )
 {
     assert( ifs.is_open() );
     
     // erasing any previous records
     Erase();
     
     // 1. read node locations: 1 headlines gives number of nodes and dimensions of model
     // ---------------------------------------------------------------------------------
     char text_line[NAME_STRING];
     const char* const delims =" ,\t,:,\n,\r";
 
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     size_t nnodes = atoi(strtok( text_line, delims ));
     uint32_t dim    = atoi(strtok(NULL,delims));
     assert( nnodes > 0 );
     assert( dim >= 1 && dim <=3 );
     px.resize( nnodes );
     vector<double>( px ).swap( px );
     if ( dim  > 1 ) {
          py.resize( nnodes );
          vector<double>( py ).swap( py );
       }
     if ( dim == 3 ) {
          pz.resize( nnodes );
          vector<double>( pz ).swap( pz );
       }

     for ( size_t i{0U}; i<nnodes; i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          px[i] = atof(strtok( text_line, delims ));
          if ( dim  > 1U ) py[i] = atof(strtok(NULL,delims));
          if ( dim == 3U ) pz[i] = atof(strtok(NULL,delims));
       }
       
     // 2. pelmt and plist, mixed mesh indicated by third number  
     //    (number indicates how many element types, thus !1= mixed mesh)
     // -----------------------------------------------------------------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     size_t nelements = atoi(strtok( text_line, delims ));
     size_t netypes   = atoi(strtok(NULL,delims));
     assert( nelements > 0 );
     assert( netypes < 20 );
     plist.resize(nelements);
     deque<vector<int64_t> >( plist ).swap( plist );
     if ( netypes == 1 ) hybrid_mesh_ = false;
     else                hybrid_mesh_ = true;
     if ( hybrid_mesh_ ) pelmt.reserve(nelements);
     
     for ( size_t i{0U}; i<nelements; i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          int64_t  etype = atoi(strtok( text_line, delims ));
          assert( etype >= 0 && etype < 50 );
          pelmt.push_back( static_cast<int8_t>(etype) );
          int64_t  npe  = atoi(strtok(NULL,delims));
          assert( npe > 2 );
          plist[i].resize( static_cast<uint32_t>(npe) );
          vector<int64_t>( plist[i] ).swap( plist[i] );
          for ( size_t j{0U}; j<static_cast<uint32_t>(npe); j++ )
            plist[i][j] = atoi(strtok(NULL,delims));
       }
     
     
     // 3. pfverts
     // ----------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     int64_t  npe = atoi(strtok( text_line, delims ));
     assert( nelements == npe );
     pfverts.resize( static_cast<int64_t>(nelements) );
     deque<vector<int64_t> >( pfverts ).swap( pfverts );
     for ( auto i{0U}; i<static_cast<uint32_t>(nelements); i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          npe = atoi(strtok( text_line, delims )); // neigbor elements per element (=faces)
          pfverts[i].resize( static_cast<uint32_t>(npe) );
          vector<int64_t>( pfverts[i] ).swap( pfverts[i] );
          for ( size_t j{0U}; j<static_cast<uint32_t>(npe); j++ )
            pfverts[i][j] = atoi(strtok(NULL,delims));
       }
     
     
     // 4. boundary flags
     // -----------------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     size_t nbnodes = atoi(strtok( text_line, delims )); // neigbor elements per element (=faces)
     assert( static_cast<int>(nbnodes) <= nnodes );
     bflags.resize(nnodes,0);
     for ( auto i{0U}; i<nbnodes; i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          size_t nid   = atoi(strtok( text_line, delims ));
          int64_t  bflag = atoi(strtok(NULL,delims));
          assert( bflag <= 0 );
          bflags[nid] = static_cast<int8_t>(bflag);
       }
   
    // 5. where the first face or interface - if any start in the records
    // ------------------------------------------------------------------
    do ifs.getline( text_line, 256 );
    while ( (isCommentLine(text_line) && !ifs.eof()) );
    int64_t  face = atoi(strtok( text_line, delims ));
    assert( face <= nelements );
    first_face_ = face;
    face  = atoi(strtok( text_line, delims ));
    assert( face <= nelements );
    first_interface_ = face;
   
 } // end InText





 
 
/**
    Clear() for VSet
*/
void VData::Erase()
 {
    hybrid_mesh_ = false;
    px.resize(0); vector<double>().swap( px );
    py.resize(0); vector<double>().swap( py );
    pz.resize(0); vector<double>().swap( pz );
    pelmt.clear();
    plist.clear(); 
    pfverts.clear();    
    bflags.clear();
    first_face_ = 0U;
    first_interface_ = 0U;
 }
  
 
 
/**
    Returns dimensions edge lengths of bounding box.
*/
void VData::CoordinateRange( char coordinate_axis, double& cmin, double& cmax ) const
 {
    if ( coordinate_axis == 'x' || coordinate_axis == 'X' ) {
         cmax = cmin = px[0];
         for ( size_t i=1; i<px.size(); i++ )
           { 
              if ( cmin > px[i] ) cmin = px[i];
              if ( cmax < px[i] ) cmax = px[i];
           }
      }
    else if ( coordinate_axis == 'y' || coordinate_axis == 'Y' ) {
         cmax = cmin = py[0];
         for ( size_t i=1; i<py.size(); i++ )
           { 
              if ( cmin > py[i] ) cmin = py[i];
              if ( cmax < py[i] ) cmax = py[i];
           }
      }
    else if ( coordinate_axis == 'z' || coordinate_axis == 'Z' ) {
         cmax = cmin = pz[0];
         for ( size_t i=1; i<pz.size(); i++ )
           { 
              if ( cmin > pz[i] ) cmin = pz[i];
              if ( cmax < pz[i] ) cmax = pz[i];
           }
      }
    else {
       cout <<"\nVData::CoordinateRange: Coordinate axis specifier not recognized: ";
       cout << coordinate_axis << endl;
       cout <<"\n\tNot returning any values."<< endl;
     }
 } 
 
 


/**
    Scales the chosen coordinate range to the supplied minimum - maximum range.
*/
void VData::ScaleCoordinateToRange( char coordinate_axis, double cmin, double cmax )
 {
    double oldmin(-1.0e30), oldmax(1.0e30);
    CoordinateRange( coordinate_axis, oldmin, oldmax );
    
    // test whether we are already O.K.
    if ( cmax == oldmax && cmin == oldmin ) return;
  
    double  old_range = oldmax - oldmin;
    double  new_range = cmax   - cmin;  

    if ( coordinate_axis == 'x' || coordinate_axis == 'X' )
      for ( size_t i{0U}; i<px.size(); i++ )
        px[i] = cmin + ((px[i] - oldmin)/old_range) * new_range;

    else if ( coordinate_axis == 'y' || coordinate_axis == 'Y' )
      for ( size_t i{0U}; i<py.size(); i++ )
        py[i] = cmin + ((py[i] - oldmin)/old_range) * new_range;

    else if ( coordinate_axis == 'z' || coordinate_axis == 'Z' )
      for ( size_t i{0U}; i<pz.size(); i++ )
        pz[i] = cmin + ((pz[i] - oldmin)/old_range) * new_range;
 } 
 
 

/**
    comparitor for VSets.
    
    @attention node coordinates are ignored because such a comparison would depend on precision.
*/
bool  VData::operator==( const VData& vd ) const
 {
    bool return_value(true);
    if ( !(hybrid_mesh_ == vd.hybrid_mesh_) ) return false;
    if ( first_face_      != vd.first_face_ ) return false;
    if ( first_interface_ != vd.first_interface_ ) return false;
    
//    if ( !(px == vd.px) ) { cerr<<"\nVData::operator== failed 'px' comparison."; return_value = false; }
//    if ( !(py == vd.py) ) { cerr<<"\nVData::operator== failed 'py' comparison."; return_value = false; }
//    if ( !(pz == vd.pz) ) { cerr<<"\nVData::operator== failed 'pz' comparison."; return_value = false; }

    if ( !(pelmt == vd.pelmt) ) {
         cerr<<"\nVData::operator== failed 'pelmt' (finite-element type) comparison.";
         if ( pelmt.size() != vd.pelmt.size() ) cerr <<"\nelement type records have different sizes.\n";
         for ( size_t i{0U}; i<pelmt.size(); ++i )
           if ( pelmt[i] != vd.pelmt[i] )
             cerr <<"\n\t\t"<< parseFiniteElementType(pelmt[i]) <<" "<< i <<": "<< static_cast<int>(pelmt[i]) <<" vs "<< static_cast<int>(vd.pelmt[i]);
         cerr << endl << endl;
         return_value = false;
      }
    
    if ( !(plist == vd.plist) ) {
         cerr<<"\nVData::operator== failed 'plist' (nodes of element) comparison.";
         for ( size_t i{0U}; i<plist.size(); ++i )
           if ( plist[i] != vd.plist[i] ) {
                string csmp_class{"Element"};
                if ( i >= Elements() && i<Elements() + Faces() ) csmp_class = "Face";
                if ( i >= Elements() + Faces() ) csmp_class = "InterFace";
                if ( plist[i].size() != vd.plist[i].size() ) {
                     cerr <<"\n\t\t"<< csmp_class <<" "<< i <<" nodes-per-element records have different sizes";
                     cerr <<" ("<< plist[i].size() <<" vs "<< vd.plist[i].size() <<").";
                  }
                else {
                     cerr <<"\n"<< csmp_class <<" "<< parseFiniteElementType(pelmt[i]) <<" "<< i <<": 'plist' member comparison: ";
                     for ( size_t j{0}; j<plist[i].size(); ++j )
                       cerr <<"\n\t\t"<< j <<": "<< plist[i][j] <<" vs "<< vd.plist[i][j];
                  }
             }
         cerr << endl << endl;
         return_value = false;
      }
    
    if ( !(pfverts == vd.pfverts) ) {
         cerr<<"\nVData::operator== failed 'pfverts' (element neighbor) comparison.\n";
         for ( size_t i{0U}; i<pfverts.size(); ++i )
           if ( pfverts[i] != vd.pfverts[i] ) {
                string csmp_class{"Element"};
                if ( i >= Elements() && i<Elements() + Faces() ) csmp_class = "Face";
                if ( i >= Elements() + Faces() ) csmp_class = "InterFace";
                if ( pfverts[i].size() != vd.pfverts[i].size() ) {
                     cerr <<"\n\t\t"<< csmp_class <<" "<< i <<" neighbor records have different sizes";
                     cerr <<" ("<< pfverts[i].size() <<" vs "<< vd.pfverts[i].size() <<").";
                  }
                else {
                     cerr <<"\n"<< csmp_class <<" "<< parseFiniteElementType(pelmt[i]) <<" "<< i <<": 'pfvert' member comparison: ";
                     for ( size_t j{0}; j<pfverts[i].size(); ++j )
                       cerr <<"\n\t\t"<< j <<": "<< pfverts[i][j] <<" vs "<< vd.pfverts[i][j];
                  }
             }
         cerr << endl << endl;
         return_value = false;
      }
      
    if ( !(bflags == vd.bflags) ) {
         cerr<<"\nVData::operator== failed 'bflags' (box-boundary flag) comparison.";
         if ( bflags.size() != vd.bflags.size() ) cerr <<"\n\tboundary flag records have different sizes.\n";
         for ( size_t i{0U}; i<bflags.size(); ++i )
           if ( bflags[i] != vd.bflags[i] )
             cerr <<"\n\t\tnode flag "<< i <<": "<< parseBoundary(intToBOX_BOUNDARY(bflags[i])) <<" vs "<< parseBoundary(intToBOX_BOUNDARY(vd.bflags[i]));
         cerr << endl << endl;
         return_value = false;
      }

    if ( !(gflags_ == vd.gflags_) ) {
         cerr<<"\nVData::operator== failed 'gflags' (node-parent BREP entity identifier) comparison.";
         if ( gflags_.size() != vd.gflags_.size() ) cerr <<"\n\tTOPOTYPE flag records have different sizes.\n";
         for ( size_t i{0U}; i<gflags_.size(); ++i )
           if ( gflags_[i] != vd.gflags_[i] )
             cerr <<"\n\t\tlefthandside node geometry flag "<< i <<": "<< parseTopology(static_cast<TOPOTYPE>(gflags_[i])) <<" vs "<< parseTopology(static_cast<TOPOTYPE>(vd.gflags_[i]));
         cerr << endl << endl;
         return_value = false;
      }

    if ( !(pmanifolds_ == vd.pmanifolds_) ) {
         cerr<<"\nVData::operator== failed 'pmanifolds_' (node-manifold) comparison.";
         if ( pmanifolds_.size() != vd.pmanifolds_.size() ) cerr <<"\n\tnode-manifold records have different sizes.\n";
         for ( size_t i{0U}; i<pmanifolds_.size(); ++i )
           if ( pmanifolds_[i] != vd.pmanifolds_[i] ) {
               cerr <<"\n\t\tnode records in manifold "<< i <<" are not the same: (";
               for ( auto j{0U}; j<pmanifolds_[i].first.size(); j++ )
                 cerr  << pmanifolds_[i].first[j] <<" ";
               cerr <<") vs (";
               for ( auto j{0U}; j<vd.pmanifolds_[i].first.size(); j++ )
                 cerr  << vd.pmanifolds_[i].first[j] <<" ";
               cerr <<").";
             }
         cerr << endl << endl;
         return_value = false;
      }

    return return_value;
    
 } // end operator==




/**
 
Using the Plist a unique set of nodes is created.
If this set contains a larger number of nodes than is stored in the 
px/py/pz vectors, a FATAL_ERROR is reported and the method returs true.
 
A fatal error is reported also if the node numbering as identified from
the Plist is discontiguous.  

If the numbering of nodes is contiguous but the number of nodes in the 
px/py/pz vectors is larger than in the node set created from the Plist,
the extra nodes are deleted. This step is optional.
   

@section application Application 

Test the Vdata for internal consistency.  
*/
bool VData::DetectAndEliminateOrphanNodes( bool eliminate_orphan_nodes )
 {
    set<size_t>  node_set;
 
    // creating the unique nodeset
    for ( deque<vector<int64_t> >::const_iterator
          it=plist.begin(); it!=plist.end(); it++ )
      for ( vector<int64_t>::const_iterator
            vt=(*it).begin(); vt!=(*it).end(); vt++ )
        node_set.insert( *vt );
        
    // if nothing can be done because there are not enough node coordinates
    if ( node_set.size() > px.size() )
      throw csmp::Exception( FATAL_ERROR, "VData::CheckForOrphanNodes",
                     "plist contains more nodes node coordinate arrays");
 
    // detecting discontinuities in the node numbering
    if ( (*max_element( node_set.begin(), node_set.end() )) >= node_set.size() )
      throw csmp::Exception( FATAL_ERROR, "VData::CheckForOrphanNodes",
                     "there is a discontinuity in the node-numbers stored in 'plist'");
    
    // if there are extra nodes but everything else is OK, these can be removed 
    // swap trick is used to trim extra capacity from vectors                
    if ( node_set.size() < px.size() ) {
          // if the nodes are not stored in the plist they are orphan and can be removed
          if ( eliminate_orphan_nodes ) {
               px.resize( node_set.size() );
               vector<double>( px ).swap( px );
               if ( !py.empty() ) {
                    py.resize( node_set.size() );
                    vector<double>( py ).swap( py );
                 }
               if ( !pz.empty() ) {
                    pz.resize( node_set.size() );
                    vector<double>( pz ).swap( pz );
                 }
               //throw csmp::Exception( INFO, "VData::CheckForOrphanNodes",
                              //"there were orphan nodes but they have been removed as requested");
                 cout << "\nVData::CheckForOrphanNodes -- there were orphan nodes but they have been removed as requested\n";
             }
          return true;
      }

    return false;
    
 } // end CheckForOrphanNodes




/**
 
The original VSet is condensed to the elements stored as keys (first arg)
in the supplied map. Also, only the nodes connected to these elements
are retained.  

As a second step the elements and nodes are renumbered contiguously
and corresponding amendments are made to the pfverts, plist and 
node arrays.    

@section arguments Input Arguments

The supplied map must contain the new element IDs (0..n) accessible via
the old element IDs (x..y) depending on which regions were chopped 
away or selected from the entire model. In the original design, the map 
is output from the ModelTopology class.  

@section application Application 

Used to make region selections from the supplied input meshes.  
*/
void VData::ReduceTo( const map<size_t,size_t>& o_n_elmt_ids, ///< the element  idx mapping from old to new
                      map<size_t,size_t>& o_n_node_ids )      ///< the old->new node idx mapping but only for the nodes that are kept
 {
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    // checks
    bool with_connectivity( !pfverts.empty() );
    if ( plist.empty() ) {
         cerr <<"\nVData::ReduceTo: The 'plist' ";
         cerr <<"vector is empty. Nothing was done."<< endl;
         return;
      }
    if ( o_n_elmt_ids.size() > plist.size() ) {
         cerr <<"\nVData::ReduceTo: Warning: The VData's 'plist' contains data for less elements ";
         cerr <<"than the map with the translation of old to new element indices.";
         cerr <<" This may indicate an inconsistency; nothing was done."<< endl;
         return;
      }
    // if there are faces or interfaces, the mesh cannot be handled by model
    if ( plist.size() > first_face_ or first_face_ != first_interface_ ) {
         cerr <<"\nVData::ReduceTo: Warning: The VData contains faces or interfaces. ";
         cerr <<" This method cannot handle such VData; nothing was done."<< endl;
         return;
      }
    const size_t old_nelements = plist.size();
   
   
    // 0. resizing the element type container 'pelmt'
    // ----------------------------------------------
    if ( HybridElementTypeMesh() ) {
         vector<int32_t> new_element_types( o_n_elmt_ids.size() );
         for ( const auto& it : o_n_elmt_ids ) {
              // range-checked access of the 2 containers using at()
              new_element_types.at( it.second ) = pelmt.at( it.first );
           }
	       pelmt.assign( new_element_types.begin(), new_element_types.end() );
      }
    // else nothing needs to be done
   
    // checking wether the mesh is still mixed element type and resizing 
    // the element type deque if it is not
    set<int8_t>  n_etypes( pelmt.begin(), pelmt.end() );
    if ( n_etypes.size() == 1U ) {
         pelmt.resize(1U);
         vector<int8_t>( pelmt ).swap( pelmt );
      }


    // 1. erasing abandonned plist & pfverts entries
    // ---------------------------------------------
    // copying those entries which are to be retained
    // and creating a node correspondance map with the new node IDs
    vector<vector<int64_t> >  new_plist( o_n_elmt_ids.size() );
    vector<vector<int64_t> >  new_pfverts( with_connectivity ? o_n_elmt_ids.size() : 0 );
    size_t                   n_node(0U);
    //  old   & new node numbers
    o_n_node_ids.clear();
    assert( plist.size() == pfverts.size() || !with_connectivity );
    for ( const auto& it : o_n_elmt_ids )
      {
         const size_t  eidx( it.first );
         const size_t  neidx( it.second );
         if ( eidx > plist.size()-1U ) {
              cerr <<"\nVData::ReduceTo: attempt to read element "<< eidx;
              cerr <<" of original 'plist' with size: "<< plist.size();
           }
         if ( neidx > new_plist.size()-1U ) {
              cerr <<"\nVData::ReduceTo: attempt to read element "<< neidx;
              cerr <<" of renumbered 'plist' with size: "<< new_plist.size();
           }
         new_plist[neidx]   = plist[eidx];
         if( with_connectivity ) new_pfverts[neidx] = pfverts[eidx];
         // recording which nodes are retained and creating an access map for them
         // in which the old node IDs are used as keys for the new ones
         for ( vector<int64_t>::const_iterator
               nit=plist[eidx].begin(); nit!=plist[eidx].end(); nit++ ) {
              pair<map<size_t,size_t>::iterator,bool>
              iit = o_n_node_ids.insert( make_pair( *nit, n_node ) );
              // if a new node was added the node number is incremented
              if ( iit.second == true ) n_node++;
           }
      }
    // verifying that plist and pfverts do not contain empty elements
    for ( const auto& it : new_plist )
      if ( it.empty() )
        throw csmp::Exception( FATAL_ERROR, "VData::ReduceTo",
                              "The reduced 'plist' deque contains empty entries. Unable to continue");
    if( with_connectivity )
        for ( const auto& it : new_pfverts )
          if ( it.empty() )
            throw csmp::Exception( FATAL_ERROR, "VData::ReduceTo",
                                  "The reduced 'pverts' deque contains empty entries. Unable to continue");

    // reassigning the plist and pfverts now
    plist.assign( new_plist.begin(), new_plist.end() );
    new_plist.clear();
    
    if( with_connectivity )
      {
          pfverts.assign( new_pfverts.begin(), new_pfverts.end() );
          new_pfverts.clear();

          // 2. 'pfverts': assigning new contiguous element IDs to 'pfverts' map
          // -------------------------------------------------------------------
          bool first_incidence(true);
          for ( deque<vector<int64_t> >::iterator it=pfverts.begin(); it!=pfverts.end(); it++ )
            for ( vector<int64_t>::iterator pit=(*it).begin(); pit!=(*it).end(); pit++ )
              // only if there was a neighboring element before its ID is updated
              if ( (*pit) > 0 ) {
                   map<size_t,size_t>::const_iterator
                     eit=o_n_elmt_ids.find( *pit );
                   if ( eit == o_n_elmt_ids.end() ) {
                         if ( first_incidence ) {
                               if ( csmp_error.Verbose() ) {
                                   cerr <<"\nVData::ReduceTo 'pfvert' neighbor element ID could not be updated for element ";
                                   cerr << *pit <<" and possible others.\n Treating them as 'REGION_BOUNDARY'"<< endl;
                                 }
                               first_incidence = false;
                           }
                         *pit = REGION_BOUNDARY;
                     }
                   else
                   *pit = static_cast<int64_t>((*eit).second);
                }
      }

    // 3. 'plist': updating node IDs if necessary
    // ------------------------------------------
    if ( px.size() != o_n_node_ids.size() ) {
         for ( deque<vector<int64_t> >::iterator it=plist.begin(); it!=plist.end(); it++ )
           for ( vector<int64_t>::iterator pit=(*it).begin(); pit!=(*it).end(); pit++ )
             // only if there was a neighboring element before its ID is updated
              {
                  map<size_t,size_t>::const_iterator nit=o_n_node_ids.find( *pit );
                  if ( nit == o_n_node_ids.end() )
                     throw csmp::Exception( ERROR, "VData::ReduceTo",
                                                   "'plist' node ID could not be updated");
                  else *pit = (*nit).second;
               }

         // 4. Eliminating orphan nodes
         // ---------------------------
         if( csmp_error.Verbose() )
             cout <<"\nVData::ReduceTo: eliminating orphan nodes..."<< endl;
         // copying the retained nodes over into new vectors
         // px
         vector<double>  new_px( o_n_node_ids.size() );
         //        old_ID    new_ID
         for ( map<size_t,size_t>::const_iterator
               nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
           new_px[ (*nit).second ] = px[ (*nit).first ];
         // updating px
         px = new_px;
         new_px.clear();
         // py
         if ( !py.empty() ) {
              vector<double>  new_py( o_n_node_ids.size() );
              //        old_ID    new_ID
              for ( map<size_t,size_t>::const_iterator
                    nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
                new_py[ (*nit).second ] = py[ (*nit).first ];
              // updating px
              py = new_py;
              new_py.clear();
           }
         // pz
         if ( !pz.empty() ) {
              vector<double>  new_pz( o_n_node_ids.size() );
              //        old_ID    new_ID
              for ( map<size_t,size_t>::const_iterator
                    nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
                new_pz[ (*nit).second ] = pz[ (*nit).first ];
              // updating px
              pz = new_pz;
              new_pz.clear();
           }

         // 5. updating boundary flags
         // --------------------------
         vector<int8_t> new_bflags(px.size(),0);
         //        old_ID    new_ID
         for ( map<size_t,size_t>::const_iterator
               nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
           new_bflags[ (*nit).second ] = bflags[ (*nit).first ];
         bflags = new_bflags;
         new_bflags.clear();

         // 6. updating geometry flags
         // --------------------------
         vector<int8_t> new_gflags(px.size(),0);
         //        old_ID    new_ID
         for ( map<size_t,size_t>::const_iterator
               nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
           new_gflags[ (*nit).second ] = gflags_[ (*nit).first ];
         gflags_ = new_gflags;
         new_gflags.clear();
      }

    // 6. updating the 'mixed_mesh' boolean variable
    // ---------------------------------------------
    // (this check will not work if there are two element types with the same
    //  number of nodes)
    size_t  nodes_per_element = plist[0].size();

    hybrid_mesh_=false;
    for ( const auto& it : plist )
      if ( it.size() != nodes_per_element ) {
           hybrid_mesh_ = true;
           break;
        }

    // 7. updating face and interface information
    // ------------------------------------------
    // it is assumed that only elements as opposed to faces or interfaces were eliminated
    const size_t n_elements_eliminated(old_nelements - plist.size());
    first_face_      -= n_elements_eliminated;
    first_interface_ -= n_elements_eliminated;

    // 8, reporting
    // ------------
    if ( plist.size() < old_nelements ) {
        if ( csmp_error.Verbose() ) {
             cout <<"\nVData::ReduceTo: Successfully reduced number of elements by "<< n_elements_eliminated <<" elements from ";
             cout << old_nelements <<" to "<< plist.size() << endl;
          }
    }

 } // end ReduceTo
 
 
 
 
/**
    Considering all known element types, returns their interpolation order.
    
    @attention assumption is made that all elements in the mesh use the same
    order element interpolation functions.
*/
uint32_t  VData::OrderOfFiniteElementInterpolationFunctions() const
 {
     for ( vector<int8_t>::const_iterator
           it=PelmtBegin(); it!=PelmtEnd(); it++ ) {
		  if ( (*it) == LINEAR_CUBOID ) return 1;
		  if ( (*it) == LINEAR_RECTANGLE ) return 1;
          if ( (*it) == LINEAR_BAR ) return 1;
          if ( (*it) == QUADRATIC_BAR ) return 2; 
          if ( (*it) == CUBIC_BAR ) return 3;  
          if ( (*it) == LINEAR_TRIANGLE ) return 1;
          if ( (*it) == LINEAR_TRIANGLE3D ) return 1;
          if ( (*it) == BARYCENTRIC_LINEAR_TRIANGLE ) return 1;
          if ( (*it) == QUADRATIC_TRIANGLE ) return 2;
          if ( (*it) == BARYCENTRIC_QUADRATIC_TRIANGLE ) return 2;
          if ( (*it) == CUBIC_TRIANGLE ) return 3;
          if ( (*it) == LINEAR_TETRAHEDRON ) return 1;
          if ( (*it) == QUADRATIC_TETRAHEDRON ) return 2;
          if ( (*it) == BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return 2;
          if ( (*it) == CUBIC_TETRAHEDRON ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_BAR ) return 1;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_BAR ) return 2;
          if ( (*it) == ISOPARAMETRIC_CUBIC_BAR ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_TRIANGLE ) return 1;
          if ( (*it) == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ) return 1;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return 2; 
          if ( (*it) == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return 2;
          if ( (*it) == ISOPARAMETRIC_CUBIC_TRIANGLE ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return 1;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return 2;
          if ( (*it) == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return 2;
          if ( (*it) == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_PYRAMID ) return 1; 
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return 2;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return 2; 
          if ( (*it) == ISOPARAMETRIC_CUBIC_PYRAMID ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_PRISM ) return 1;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return 2;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return 2;
          if ( (*it) == ISOPARAMETRIC_CUBIC_PRISM ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return 1;
          if ( (*it) == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return 1; 
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return 2;
          if ( (*it) == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return 2; 
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return 2;
          if ( (*it) == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return 3;
          if ( (*it) == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) return 1;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return 2;
          if ( (*it) == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return 2;
          if ( (*it) == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return 3;
       }
       
    cerr <<"\nVData::OrderOfFiniteElementInterpolationFunctions: No indications of order found; returning 0!"<< endl;
    return 0;
       
 } // end 



/**
    Runs a series of tests to establish whether there is a plausible PFverts array
    without checking the actual inter-element connectivity.
    Certainty is built via negative discrimination.
*/
bool VData::WithNeighbourConnectivity() const
 {
     // empty
     if ( pfverts.empty() ) return false;
    return true;
 }



/// range of the X coordinate (increasing to the east)
pair<double,double>  VData::X_Range() const
 {
   if ( px.empty() ) return make_pair(0.,0.);
   vector<double>::const_iterator xt_max = max_element(px.begin(), px.end());
   vector<double>::const_iterator xt_min = min_element(px.begin(), px.end());
   return make_pair( (*xt_min), (*xt_max) );
 }

/// range of the Y coordinate (increasing upward and often indicating elevation above sealevel)
pair<double,double>  VData::Y_Range() const
 {
    if ( py.empty() ) return make_pair(0.,0.);
    vector<double>::const_iterator yt_max = max_element(py.begin(), py.end());
    vector<double>::const_iterator yt_min = min_element(py.begin(), py.end());
    return make_pair( (*yt_min), (*yt_max) ); 
    //bool y_zero = ( fabs( (*yt_max)-(*yt_min) ) <= numeric_limits<double>::epsilon() ) ? true : false; 
 }

/// range of the Z coordinate (from north to south)
pair<double,double>  VData::Z_Range() const
 {
   if ( pz.empty() ) return make_pair(0.,0.);
   vector<double>::const_iterator zt_max = max_element(pz.begin(), pz.end());
   vector<double>::const_iterator zt_min = min_element(pz.begin(), pz.end());
   return make_pair( (*zt_min), (*zt_max) );
   //bool z_zero = ( fabs( (*zt_max)-(*zt_min) ) <= numeric_limits<double>::epsilon() ) ? true : false; 
 }




/** 
   Flips clockwise-numbered elements, into counter-clockwise right-hand rule compliant orientation; lower dimensional elements are made consistent; returns how many were flipped
   
      @attention method expects elements to be of CSMP_FEM_TYPE
      
      @attention works only for linear and quadratic elements
      
      @attention method does not touch the line elements
*/
size_t VData::RenumberElementsCounterClockwise2D()
 {
    size_t orientation_changes(0U);
    CSMP_FEM_TYPE   etype = (HybridElementTypeMesh()==false) ? parseFiniteElementTypeEnum( ElementType(0) ) : UNKNOWN;
    bool      is_triangle = (etype != UNKNOWN) ? isTriangularElement( etype ) : false;
    bool is_quadrilateral = (etype != UNKNOWN) ? isQuadrilateralElement( etype ) : false;
    
    // 1. processing the elements that have the same dimension as the 2D mesh (with at least 3 nodes)
    // ----------------------------------------------------------------------------------------------
    const vector<int8_t>::const_iterator end = PelmtEnd();
    vector<int8_t>::const_iterator       eit = PelmtBegin();
    deque<vector<int64_t> >::iterator     pls = PlistBegin();
    deque<vector<int64_t> >::iterator     pfv = PfvertsBegin();
    
    while( eit != end ) {
         // if this is a triangle or quadrilateral
         if ( HybridElementTypeMesh() ) {
              etype = parseFiniteElementTypeEnum(*eit);
              is_triangle      = isTriangularElement( etype );
              is_quadrilateral = isQuadrilateralElement( etype ); 
           }
         if ( is_triangle or is_quadrilateral ) {
              // computing (signed) cross product=normal of vector from origin node to next (clockwise) and origin backwards (ccw)
              const double v1x = px[ (*pls)[2] ] - px[ (*pls)[1] ];
              const double v1y = py[ (*pls)[2] ] - py[ (*pls)[1] ];
              const double v2x = px[ (*pls)[0] ] - px[ (*pls)[1] ];
              const double v2y = py[ (*pls)[0] ] - py[ (*pls)[1] ];
              const double crossProduct = v1x * v2y - v1y * v2x;
              // if this product is not positive, the element node numbering must be flipped
              if ( crossProduct < 0. ) {
                   // linear elements (less than or equal to 4 nodes)
                   const size_t nodes = (*pls).size();
                   // reverse plist and pfverts of linear elements
                   if ( nodes <= 4 ) {
                        reverse( (*pls).begin(), (*pls).end() );  
                     }
                   // else quadratic triangles or quads thus far 
                   else {
                        if ( is_triangle ) {
                             // reverse corner nodes, then midside nodes, bubble node stays in barycentric element 
                             assert( nodes <= 7 );
                             reverse( (*pls).begin(), next((*pls).begin(),3) );
                             reverse( next((*pls).begin(),4), next((*pls).begin(),6) );
                          }
                        else { // quadratic quadrilateral with and without bubble node
                             assert( nodes <= 9 );
                             reverse( (*pls).begin(), next((*pls).begin(),8) );  
                          }
                     }
                   // in all cases, the neighbor numbering has to be reversed as well
                   if ( !(*pfv).empty() ) reverse( (*pfv).begin(), (*pfv).end() );
                   orientation_changes++;  
                }
           }
         eit++;
         pls++;
         pfv++;
      }
      
    if ( orientation_changes > 0 ) {
         cerr <<"\nVData::RenumberElementsCounterClockwise2D: flipped node numbering of "<< orientation_changes;
         cerr <<" surface elements from clockwise to counter-clockwise.\n";
      }
    
    return orientation_changes;
 
 } // end RenumberElementsCounterClockwise2D
 
 
 
 
 
 
 /**
   Line elements are oriented in a consistent way for internal 1D regions within a 2D domain.
   
   @brief To ensure that the normals pointing out of a lower (1D) dimensional elements located at the model boundary,
   and that they all point in the same direction for internal boundaries, the nodes of such elements must be locally ordered consistently;
   meaning these line elements are aligned with the node-numbering of the faces of 2D surface elements located at the
   model boundary whose nodes are numbered counter clockwise.
   
   The normals of the line elements are found by 90-degrees clockwise rotation of the tangent vector of the line elements,
   which points from node 0 to node 1 (in a Linear or Quadratic Line Element).
   
   @attention this method expects that there already is connectivity information for line elements
   
   @attention if the line elements do not know their neighbors already, this method will not work!
   
   @attention assumes node numbering in 'plist' follows the convention corner nodes first, then Interior nodes.
   
   @attention method assumes that all manifolds (connections of more than 2 line elements at a node have been removed
   =disambiguated before. Thus it relies on a valid neighbor connectivity of the line elements, which must stored in pfverts.
   
   @attention the neighbors of line elements are situated opposite to its corner nodes (like for other simplex elements)
   
   @note algorithmic procedure: 1) chains of connected line elements are found and stored in a map, using their "root" nodes as keys.
   2) these are then traversed, flipping nodes whenever successive elements are not aligned. 3) at the model boundaries,  the chains of
   line elements are flipped if their direction is not matched with that of the faces of the surface elements at the model boundary.
   4) line element chains inside of the model are kept in whichever, yet consistent orientation.
   
  @author SKM
  @date 3/22020
  @test 2/10/21
  
  */
 void VData::CreateConsistentLineElementOrientations2D()
   {
      if ( bflags.empty() ) {
            throw csmp::Exception( ERROR, "VData::CreateConsistentLineElementOrientations2D",
                                  "Method relies on correct flagging of boundary nodes. This is not present in VData.");
        }
      if ( pfverts.empty() ) {
            cerr <<"\nVData::CreateConsistentLineElementOrientations2D: 'pfverts' is empty; nothing could be done.\n";
            return;
         }
      assert( plist.size() == pfverts.size() ); 
       
      // assuming that the mesh consists of surface and line elements
      assert( HybridElementTypeMesh()==true ); 
       
      // 1. Finding line elements and edges of surface elements located at the model boundary
      // ------------------------------------------------------------------------------------
      set<int64_t>  line_elmts, boundary_line_elmts;
      // recording the corner-node ids of the surface elements for later searching
      // face-nd-id key, boundary face node ids in correct sequence
      deque<pair<int64_t ,int64_t> > surf_elmt_face_nd_ids;

      const vector<int8_t>::const_iterator  end = PelmtEnd();
      vector<int8_t>::const_iterator        eit = PelmtBegin();
      int64_t                               elmt_idx{0U};

       while( eit != end ) {
            CSMP_FEM_TYPE etype = parseFiniteElementTypeEnum( (*eit) );
            // only recording line elements that are located at the beginning or end of a polyline = line element chain
            if ( isLineElement( etype ) ) {
                  // line elements with a single neighbor at the beginning or end of polyline
                  if ( ((pfverts[elmt_idx][0] < 0 && pfverts[elmt_idx][1] >= 0) || (pfverts[elmt_idx][0] >= 0 && pfverts[elmt_idx][1] < 0)) )
                    // but only, if that neighbor is located on the inside of the model
                    if ( bflags[ plist[elmt_idx][0] ] == NOT || bflags[ plist[elmt_idx][1] ] == NOT ) {
                         line_elmts.insert( elmt_idx );
                      }
                  // line elements on the model boundary
                 if ( bflags[ plist[elmt_idx][0] ] < 0 && bflags[ plist[elmt_idx][1] ] < 0 ) {
                       boundary_line_elmts.insert( elmt_idx );
                    }
              }
            // if this is a surface element  
            else {
                 // looping over the faces (edges), finding those at the boundary
                 // and storing their nodes in correct order for later comparisons
                 if ( isTriangularElement( etype ) ) {
                      const size_t faces{3};
                      for ( auto i{0U}; i < faces; ++i )
                        if ( pfverts[elmt_idx][i] < 0 )
                          {
                             if ( i == 0 ) { // face 0
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][1], plist[elmt_idx][2] ) );
                               }
                             else if ( i == 1 ) { // face 1
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][2], plist[elmt_idx][0] ) );
                               }
                             else { // face 2
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][0], plist[elmt_idx][1] ) );
                               }
                          }
                   }
                 // quadrilaterals
                 else if ( isQuadrilateralElement( etype ) ) {
                      const size_t faces{4};
                      for ( auto i{0U}; i < faces; ++i )
                        if ( pfverts[elmt_idx][i] < 0 )
                          {
                             if ( i == 0 ) { // face 0
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][0], plist[elmt_idx][1] ) );
                               }
                             else if ( i == 1 ) { // face 1
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][1], plist[elmt_idx][2] ) );
                               }
                             else if ( i == 2 ) { // face 2
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][2], plist[elmt_idx][3] ) );
                               }
                             else { // face 3
                                  surf_elmt_face_nd_ids.push_back( make_pair( plist[elmt_idx][3], plist[elmt_idx][0] ) );
                               }
                         }
                   } 
                 else {
                      cerr <<"\n\telement "<< elmt_idx <<": type "<< parseFiniteElementType(*eit) <<" cannot be processed.\n";
                      throw csmp::Exception( ERROR, "VData::CreateConsistentLineElementOrientations2D", "element type not recognised");
                   }
              }
            elmt_idx++;
            eit++;
         }
         
      // 2. following potential line element chains from beginning to end, re-orientating elements as necessary
      //    and recording individual polylines
      // ------------------------------------------------------------------------------------------------------
      // (after the potential re-orientation of surface elements by ModelTopology) line element orientations might be inconsistent with surface ones)
      // root-elmt & numbers of interconnected line elements in discovered chain
      map<int64_t ,deque<int64_t> > polylines;
      set<int64_t>                  processed_elmts;
      
      // looping over the line elements that are missing one neighbor, starting at the beginning or end of a chain
      for ( auto it=line_elmts.begin(); it!= line_elmts.end(); ++it )
        {
           // getting index of element and skipping elements that have already been processed
           elmt_idx = (*it);
           if ( processed_elmts.find(elmt_idx) != processed_elmts.end() ) continue;

           // isolated line elements (no neighbors) are skipped as well
           if ( pfverts[elmt_idx][0] < 0 && pfverts[elmt_idx][1] < 0 ) {
                polylines.insert( make_pair( elmt_idx, deque<int64_t>{elmt_idx} ) );
                processed_elmts.insert( elmt_idx );
                continue;
             }

         // storing the element as the first in a new line element sequence
          auto chain_it=polylines.insert( make_pair( elmt_idx, deque<int64_t>{elmt_idx} ) );
          // making sure that the element was indeed inserted (else it is a duplicate)
          assert( chain_it.second == true );
          processed_elmts.insert( elmt_idx );
          
          // traversing the line element chain until running out of neighbors
          // ----------------------------------------------------------------
          // NB: like for other simplex elements, the line element neighbor is located opposite to the node with the same number.
          //     nbor1 | [nd0]-line element0-[nd1] | nbor0
          
          // -----------------------
          // forward chain traversal
          // -----------------------
          // (the neighbor element first node will is shared with the zeroth node of the following line element)
          if ( pfverts[elmt_idx][0] >= 0L && pfverts[elmt_idx][1] < 0L ) {
              bool end_of_polyline(false);
              while ( end_of_polyline == false )
                {
                   // if this next element's first neighbor is element 'elmt_idx', everything is fine and no flip is required
                   bool flip = ( pfverts[ pfverts[elmt_idx][0] ][1] == elmt_idx ) ? false : true;
                   // else we flip this element's nodes and neighbors
                   if ( flip == true ) {
                        swap( plist[ pfverts[elmt_idx][0] ][0],  plist[ pfverts[elmt_idx][0] ][1] );
                        swap( pfverts[ pfverts[elmt_idx][0] ][0], pfverts[ pfverts[elmt_idx][0] ][1] );
                     }
                   // we move to the next element
                   elmt_idx = pfverts[elmt_idx][0];
                   // then we store this element in the polyline
                   (*chain_it.first).second.push_back( elmt_idx );
                   // and record it as processed
                   auto inserted = processed_elmts.insert( elmt_idx );
                   // checking that we are not visiting a previously visited element again
                   assert( inserted.second == true );
                   // exit condition (if the far node of line element is on the BOX_BOUNDARY)
                   if ( pfverts[elmt_idx][0] < 0L ) // TODO: do we need this extra condition?  || bflags[ plist[elmt_idx][0] ] < 0 )
                     end_of_polyline = true;
                }
               continue;
            }
          // -------------------------------
          // backward chain traversal
          // -------------------------------
          if ( pfverts[elmt_idx][0] < 0L && pfverts[elmt_idx][1] >= 0L ) {
              bool beginning_of_polyline(false);
              while ( beginning_of_polyline == false )
                {
                   // if the next neighbor's first neighbor element is element 'elmt_idx', everything is fine and no flip is required,
                   bool flip = ( pfverts[ pfverts[elmt_idx][1] ][0] == elmt_idx ) ? false : true;
                   // and, if necessary, we flip its nodes and neighbors.
                   if ( flip == true ) {
                        swap( plist[ pfverts[elmt_idx][1] ][0],  plist[ pfverts[elmt_idx][1] ][1] );
                        swap( pfverts[ pfverts[elmt_idx][1] ][0], pfverts[ pfverts[elmt_idx][1] ][1] );
                     }
                   // we move to this next element
                   elmt_idx = pfverts[elmt_idx][1];
                   // then we store this element in the polyline.
                   (*chain_it.first).second.push_back( elmt_idx );
                   auto inserted = processed_elmts.insert( elmt_idx );
                   assert( inserted.second == true );
                   // exit condition (if the far node of line element is on the BOX_BOUNDARY)
                   if ( pfverts[elmt_idx][1] < 0L ) // TODO: do we need this extra condition?  || bflags[ plist[elmt_idx][0] ] < 0 )
                     beginning_of_polyline = true;
                }
            }
       } // for line_elements
       
      // NOTE: since every chain has a beginning and an end, it would normally be stored twice, but:
      //       - of the interior chains, only the ones starting with the lower element number are kept
      //       - for the ones surrounding the model the one consistent with the counter-clockwise numbering of the higher-dim. ele is stored
      //       - the only line elements left now, are those forming part of loops, these must be consistent in their orientation with
      //         the faces of the elements that they surround.
      if ( !polylines.empty() ) {
          cout <<"\n\nVData::CreateConsistentLineElementOrientations2D: processed "<< polylines.size() <<" polylines with the line elements:"<< endl;
          for ( auto it : polylines ) {
               cout <<"\t"<< it.first <<": ";
               for ( auto idx : it.second ) cout << idx <<" ";
               cout << endl;
            }
          cout << endl;
        }
      else {
           cerr <<"\nVData::CreateConsistentLineElementOrientations: WARNING: No changes were made. Unable to process line element chains. ";
           cerr <<"\n\n\t\t"<<"line elements processed: "<< processed_elmts.size() << endl;
        }
      
      // 3. Reordering chains that are located on the model boundary to make them consistent with the counter-clockwise element numbering
      // --------------------------------------------------------------------------------------------------------------------------------
      // using: //  face-nd-ids,      elmt, face   to verify that elements are indeed oriented correctly
      //        map<set<uint32_t>,pair<size_t,size_t> > surf_elmt_face_nd_ids;
      //
      // getting the surface element deque ready for binary_search
      sort( surf_elmt_face_nd_ids.begin(), surf_elmt_face_nd_ids.end() );
      
      // looping over the line elements that are missing one neighbor, i.e., are at the beginning of a chain
      for ( auto it=boundary_line_elmts.begin(); it!= boundary_line_elmts.end(); ++it )
        {
           // searching for the corresponding face of a higher dimensional element
           // --------------------------------------------------------------------
           // if a surface element face with same node numbering is found the line element is already correctly oriented
           if ( binary_search( surf_elmt_face_nd_ids.begin(), surf_elmt_face_nd_ids.end(), make_pair( plist[*it][0], plist[*it][1]) ) )
             continue;
           // if the face has the opposite orientation, the line elements in the polyline have to be flipped
           if ( binary_search( surf_elmt_face_nd_ids.begin(), surf_elmt_face_nd_ids.end(), make_pair( plist[*it][1], plist[*it][0]) ) ) {
                swap( plist[*it][0], plist[*it][1] );     // swapping nodes
                swap( pfverts[*it][0], pfverts[*it][1] ); // swapping neighbors
             }
           // this is a line element with no surface element next to it?
           else {
                cerr <<"\nVData::CreateConsistentLineElementOrientations: detected detached line element "<< *it <<" at border with the nodes:\n\t\t";
                cerr << plist[*it][0] <<"("<< parseBoundary(intToBOX_BOUNDARY(bflags[plist[*it][0]])) <<"), ";
                cerr << plist[*it][1] <<"("<< parseBoundary(intToBOX_BOUNDARY(bflags[plist[*it][1]])) <<"), ";
                cerr <<" element has no higher-dimensional neighbor; its orientation was left untouched.\n";
             }
           
        } // end boundary_line_elmts
              
  } // end CreateConsistentLineElementOrientations



 // checking the node flagging for the case that the line elements are Faces or InterFaces
 /*
 int bflag{MULTIPLE}, no_nbor{0};
 if      ( pfverts[elmt_idx][0] < 0 ) { bflag = bflags[ plist[elmt_idx][1] ]; no_nbor = 0; }
 else if ( pfverts[elmt_idx][1] < 0 ) { bflag = bflags[ plist[elmt_idx][0] ]; no_nbor = 1; }
 // the element's neighbor-free side must be at an internal or external boundary (bflag<0), else
 if ( bflag >= 0 ) {
      cerr <<"\nVData::CreateConsistentLineElementOrientations2D: ";
      cerr <<"'pfvert["<< elmt_idx <<"]["<< no_nbor <<"]' entry for missing neighbour element "<< no_nbor;
      cerr <<": "<< parseBoundary( static_cast<BOX_BOUNDARY>(pfverts[elmt_idx][no_nbor]) ) <<" ("<< pfverts[elmt_idx][no_nbor] <<")";
      cerr <<" should be equivalent to boundary flag of the adjacent node: ";
      cerr << parseBoundary( static_cast<BOX_BOUNDARY>(bflag) ) <<" ("<< bflag <<")"<< endl;
   }
 assert( bflag >= MULTIPLE );
*/
 



 
 /// finds angle between 2 line elements determined by their number in the plist                                            
 double VData::AngleBetweenLineElements2D( size_t elmt1, size_t elmt2 )
  {
     assert( elmt1 < pelmt.size() );
     assert( elmt2 < pelmt.size() );
     assert( isLineElement( parseFiniteElementTypeEnum( pelmt[elmt1] ) ) );
     assert( isLineElement( parseFiniteElementTypeEnum( pelmt[elmt2] ) ) );
     
     // representing line elements as vectors
     const double v1x = px[ plist[elmt1][1U] ] - px[ plist[elmt1][0U] ];  
     const double v1y = py[ plist[elmt1][1U] ] - py[ plist[elmt1][0U] ];  
     const double v2x = px[ plist[elmt2][1U] ] - px[ plist[elmt2][0U] ];  
     const double v2y = py[ plist[elmt2][1U] ] - py[ plist[elmt2][0U] ];  
     
     // cos theta = dot-product over cross-product (length1 * length2)
     double cos_theta = (v1x*v2x + v1y*v2y) / (sqrt(v1x*v1x+v1y*v1y) * sqrt(v2x*v2x+v2y*v2y));     

     return radiansToDegrees( acos(cos_theta) );
     
  } // end AngleBetweenLineElements2D




 /// finds angle between 2 line elements determined by their number in the plist
 double VData::AngleBetweenLineElements3D( size_t elmt1, size_t elmt2 )
  {
     assert( elmt1 < pelmt.size() );
     assert( elmt2 < pelmt.size() );
     assert( isLineElement( parseFiniteElementTypeEnum( pelmt[elmt1] ) ) );
     assert( isLineElement( parseFiniteElementTypeEnum( pelmt[elmt2] ) ) );
     
     // representing line elements as vectors
     const double v1x = px[ plist[elmt1][1U] ] - px[ plist[elmt1][0U] ];
     const double v1y = py[ plist[elmt1][1U] ] - py[ plist[elmt1][0U] ];
     const double v1z = pz[ plist[elmt1][1U] ] - pz[ plist[elmt1][0U] ];
     const double v2x = px[ plist[elmt2][1U] ] - px[ plist[elmt2][0U] ];
     const double v2y = py[ plist[elmt2][1U] ] - py[ plist[elmt2][0U] ];
     const double v2z = pz[ plist[elmt2][1U] ] - pz[ plist[elmt2][0U] ];
     
     // cos theta = dot-product over cross-product (length1 * length2)
     double cos_theta = (v1x*v2x + v1y*v2y + v1z*v2z) / (sqrt(v1x*v1x+v1y*v1y+v1z*v1z) * sqrt(v2x*v2x+v2y*v2y+v2z*v2z));

     return radiansToDegrees( acos(cos_theta) );
     
  } // end AngleBetweenLineElements2D

 
 
 
 
 /**
      Computes the unit normals to the supplied surface elements and then returns the angle between them in degrees.
      
      @attention the surface elements are assumed to be planar; only the first 3 nodes are considered.
 */
double VData::AngleBetweenSurfaceElements3D( size_t elmt1, size_t elmt2 )
  {
     assert( elmt1 < pelmt.size() );
     assert( elmt2 < pelmt.size() );
     assert( CSMP_ElementSpecifications::SurfaceElement( parseFiniteElementTypeEnum( pelmt[elmt1] ) ) );
     assert( CSMP_ElementSpecifications::SurfaceElement( parseFiniteElementTypeEnum( pelmt[elmt2] ) ) );
     // are there at least 3 nodes
     assert( plist[elmt1].size() >= 3 );
     assert( plist[elmt2].size() >= 3 );

     // getting the node points of the triangle to construct the normal for
     const Point<3U> p1_0( px[plist[elmt1][0]], py[plist[elmt1][0]], pz[plist[elmt1][0]] );
     const Point<3U> p1_1( px[plist[elmt1][1]], py[plist[elmt1][1]], pz[plist[elmt1][1]] );
     const Point<3U> p1_2( px[plist[elmt1][2]], py[plist[elmt1][2]], pz[plist[elmt1][2]] );
     const Point<3U> normal_e1 = normalOfTriangle( p1_0, p1_1, p1_2 );

     const Point<3U> p2_0( px[plist[elmt2][0]], py[plist[elmt2][0]], pz[plist[elmt2][0]] );
     const Point<3U> p2_1( px[plist[elmt2][1]], py[plist[elmt2][1]], pz[plist[elmt2][1]] );
     const Point<3U> p2_2( px[plist[elmt2][2]], py[plist[elmt2][2]], pz[plist[elmt2][2]] );
     const Point<3U> normal_e2 = normalOfTriangle( p2_0, p2_1, p2_2 );

     // cos theta = dot-product over cross-product (length1 * length2)
     // already in degrees
     return angleBetweenEdges( make_pair( Point<3U>{0.,0.,0.}, normal_e1 ),
                               make_pair( Point<3U>{0.,0.,0.}, normal_e2 ) );

  } // end angleBetweenLineSegments

 

 
 

/**
   (Re)bBuild 'pfverts' container for the case it is empty or may contain unreliable information.
   computes 'pfverts' connectivity between equidimensional elements, faces and interfaces and replaces existing connectivity with it.
   
   The connectivity of line elements is established as well. Where there are line element manifolds, those of the line elements that are sharing a node
   are connected which are most closely aligned.
   
   All line elements around the perimeter of a model are connected with one another.
   
   Method calls  CreateConsistentLineElementOrientations2D()  in order to make edge orientations consistent.
   
   @attention this method assumes that all elements are numbered in counter-clockwise direction
   
   @attention this method does not change the BOX_BOUNDARY flags of the model
   
  @author SKM
  @date 2/2/2020
  @test 4/10/2021
  
  @todo SKM 4/10/21: the line-element functionality is duplicated elsewhere in CSMP and must also be fixed because it is probably not working.
  
*/
void  VData::EstablishElementConnectivity2D()
 {
   ErrorHandler& csmp_error( ErrorHandler::Instance() );
   
   if ( !pfverts.empty() ) {
         cout <<"\nVData::EstablishElementConnectivity2D: 'pfverts' not empty; deleting original content before reconstruction.\n";
         pfverts.clear();
         assert( !plist.empty() );
      }
    
   CSMP_FEM_TYPE   etype = (HybridElementTypeMesh()==false) ? parseFiniteElementTypeEnum( ElementType(0) ) : UNKNOWN;
   bool      is_triangle = (etype != UNKNOWN) ? isTriangularElement( etype ) : false;
   bool is_quadrilateral = (etype != UNKNOWN) ? isQuadrilateralElement( etype ) : false;
   bool triangle_with_all_nodes_on_boundary(false);
   
   if ( !plist.empty() )
     {
        // building face keys for surface and line elements types, collecting them into maps,
        // including the elements that they belong to with corresponding face IDs.
        //  face-nd-ids    elmt-id, face-of-element
        map<set<int64_t>,map<int64_t ,size_t> >  surf_elmt_nbors;
        map<int64_t ,set<int64_t> >              line_elmt_that_share_node;

        const deque<vector<int64_t> >::iterator end = PlistEnd();
        deque<vector<int64_t> >::iterator       eit = PlistBegin();
        int64_t  elmt_idx(0U);

        // 1. establishing connectivity
        // ----------------------------
        pfverts.resize( plist.size() );
        while( eit != end ) {
             // if this is a triangle or quadrilateral
              if ( HybridElementTypeMesh() ) {
                   etype = parseFiniteElementTypeEnum( pelmt[elmt_idx] );
                   assert ( etype != CUBIC_TRIANGLE );
                   assert ( etype != ISOPARAMETRIC_CUBIC_TRIANGLE );
                   assert ( etype != CUBIC_BAR );
                   assert ( etype != ISOPARAMETRIC_CUBIC_BAR );
                   assert ( etype != ISOPARAMETRIC_CUBIC_QUADRILATERAL );
                   is_triangle      = isTriangularElement( etype );
                   is_quadrilateral = isQuadrilateralElement( etype ); 
                }
             // triangular elements have their faces opposite of the nodes 
             // and the corner nodes are sufficient to identify them
             if ( is_triangle ) {
                  pfverts[elmt_idx].resize(3U,IRREGULAR);
                  const set<int64_t> face1({ (*eit)[1], (*eit)[2] });
                  const set<int64_t> face2({ (*eit)[2], (*eit)[0] });
                  const set<int64_t> face3({ (*eit)[0], (*eit)[1] });
                  // inserting the faces
                  // face 1
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face1_it =
                    surf_elmt_nbors.insert( make_pair( face1, map<int64_t ,size_t>({{elmt_idx,{0}}}) ) );
                  // if the face already is in the map, the element ID is added
                  if ( face1_it.second == false ) (*face1_it.first).second.insert( make_pair( elmt_idx, 0 ) );
                  // face 2
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face2_it =
                    surf_elmt_nbors.insert( make_pair( face2, map<int64_t ,size_t>({{elmt_idx,{1}}}) ) );
                  if ( face2_it.second == false ) (*face2_it.first).second.insert( make_pair( elmt_idx, 1 ) );
                  // face 3
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face3_it =
                    surf_elmt_nbors.insert( make_pair( face3, map<int64_t ,size_t>({{elmt_idx,{2}}}) ) );
                  if ( face3_it.second == false ) (*face3_it.first).second.insert( make_pair( elmt_idx, 2 ) );
               }
             else if ( is_quadrilateral ) {
                  pfverts[elmt_idx].resize(4U,IRREGULAR);
                  // see CSMP fem specifications for these conventions 
                  const set<int64_t> face1{ (*eit)[0], (*eit)[1] };
                  const set<int64_t> face2{ (*eit)[1], (*eit)[2] };
                  const set<int64_t> face3{ (*eit)[2], (*eit)[3] };
                  const set<int64_t> face4{ (*eit)[3], (*eit)[0] };
                  // inserting the faces
                  // face 1
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face1_it =
                    surf_elmt_nbors.insert( make_pair( face1, map<int64_t ,size_t>({{elmt_idx,{0}}}) ) );
                  // if the face already is in the map, the element ID is added
                  if ( face1_it.second == false ) (*face1_it.first).second.insert( make_pair( elmt_idx, 0 ) );
                  // face 2
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face2_it =
                    surf_elmt_nbors.insert( make_pair( face2, map<int64_t ,size_t>({{elmt_idx,{1}}}) ) );
                  if ( face2_it.second == false ) (*face2_it.first).second.insert( make_pair( elmt_idx, 1 ) );
                  // face 3
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face3_it =
                    surf_elmt_nbors.insert( make_pair( face3, map<int64_t ,size_t>({{elmt_idx,{2}}}) ) );
                  if ( face3_it.second == false ) (*face3_it.first).second.insert( make_pair( elmt_idx, 2 ) );               
                  // face 4
                  pair<map<set<int64_t>,map<int64_t ,size_t> >::iterator,bool> face4_it =
                    surf_elmt_nbors.insert( make_pair( face4, map<int64_t ,size_t>({{elmt_idx,{3}}}) ) );
                  if ( face4_it.second == false ) (*face4_it.first).second.insert( make_pair( elmt_idx, 3 ) );               
               }
             // any line elements (making a map of which ones share a node; this map also identifies manifolds)
             else {
                  assert( isLineElement( etype ) );
                  pfverts[elmt_idx].resize(2U,IRREGULAR);
                  pair<map<int64_t ,set<int64_t> >::iterator,bool> // node 1 (on the side of neighbor 0)
                    it0 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][1], set<int64_t>{elmt_idx} ) );
                  // if there is already an entry for the node, the elmt id is added to the set
                  if ( !it0.second )
                    (*it0.first).second.insert( elmt_idx );
                  
                  pair<map<int64_t ,set<int64_t> >::iterator,bool> // node 0 (opposite neighbor 1)
                    it1 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][0], set<int64_t>{elmt_idx} ) );
                  // if there is already an entry for the node, the elmt id is added to the set
                  if ( !it1.second )
                    (*it1.first).second.insert( elmt_idx );
               }
             elmt_idx++;
             eit++;
          }
 
        // 2. reconnecting surface type elements
        // -------------------------------------
        // where there is no neighbor an attempt is made to assign a face to a particular boundary,
        // based on its alignment with boundary normal
        size_t node0(NULL_IDX), node1(NULL_IDX);
        for ( map<set<int64_t>,map<int64_t ,size_t> >::const_iterator
              it=surf_elmt_nbors.begin(); it!=surf_elmt_nbors.end(); ++it ) {
             // there should be no more than 2 entries per face
             assert( (*it).second.size() <= 2U );
             // if there is only a single entry, the face is at the model boundary
             // ------------------------------------------------------------------ 
             if ( (*it).second.size() == 1U ) {
                  elmt_idx    = (*(*it).second.begin()).first;
                  size_t face = (*(*it).second.begin()).second;
                  // nodes of face
                  etype = (HybridElementTypeMesh()==false) ? parseFiniteElementTypeEnum( ElementType(0) ) : 
                                                             parseFiniteElementTypeEnum( ElementType(elmt_idx) );
                  is_triangle      = isTriangularElement( etype );
                  is_quadrilateral = isQuadrilateralElement( etype );
                  assert( is_triangle || is_quadrilateral );
                  if ( is_triangle ) {
                       if      ( face==0 ) { node0=1; node1=2; }
                       else if ( face==1 ) { node0=2; node1=0; }
                       else if ( face==2 ) { node0=0; node1=1; }
                    }
                  else if ( is_quadrilateral ) {
                       if      ( face==0 ) { node0=0; node1=1; }
                       else if ( face==1 ) { node0=1; node1=2; }
                       else if ( face==2 ) { node0=2; node1=3; }
                       else if ( face==3 ) { node0=3; node1=0; }
                    }
                  // creating face unit normal
                  const double dx = px[ plist[elmt_idx][node1] ] - px[ plist[elmt_idx][node0] ]; // dx=x2-x1 
                  const double dy = py[ plist[elmt_idx][node1] ] - py[ plist[elmt_idx][node0] ]; // dy=y2-y1, 
                  const double length = sqrt( dx*dx + dy*dy );
                  // normals are (-dy, dx)=clockwise and (dy, -dx)=counter-clockwise (USED HERE)
                  const double unrml[2U] = { dy/length, -dx/length };
                  // determining BOX_BOUNDARY
                  // TODO: does not work yet
                  BOX_BOUNDARY boundary(IRREGULAR);
                  bool boundary_found{false};
                  // bottom (within +/-20o of side-boundary normal (dot-product >=0.94)
                  double dotproduct = /* unrml[0U] * 0. + */ unrml[1U] * -1.;
                  if ( dotproduct >= 0.94 ) {
                       boundary       = BOTTOM;
                       boundary_found = true;
                    }
                  // right
                  else dotproduct = unrml[0U] * 1. /* + unrml[1U] * 0. */;
                  if ( dotproduct >= 0.94 && !boundary_found ) {
                       boundary       = RIGHT;
                       boundary_found = true;
                    }
                  // top
                  else dotproduct = /* unrml[0U] * 0. + */ unrml[1U] * 1.;
                  if ( dotproduct >= 0.94 && !boundary_found ) {
                       boundary       = TOP;
                       boundary_found = true;
                    }
                  // left
                  else dotproduct = unrml[0U] * -1. /* + unrml[1U] * 0. */;
                  if ( dotproduct >= 0.94 && !boundary_found )
                    boundary = LEFT;
                  // all other cases remain IRREGULAR
                  
                  // assigning boundary to the face without neighbor
                  pfverts[elmt_idx][face] = boundary;
               }
             // if there are 2 elements that share the face
             // -------------------------------------------
             else {
                  // the neighbors are recorded in the 'pfverts' map
                  const size_t elmt1      = (*(*it).second.begin()).first;
                  const size_t face_elmt1 = (*(*it).second.begin()).second;
                  const size_t elmt2      = (*next((*it).second.begin(),1)).first;
                  const size_t face_elmt2 = (*next((*it).second.begin(),1)).second;
                  pfverts[elmt1][face_elmt1] = elmt2;
                  pfverts[elmt2][face_elmt2] = elmt1;
               } 
          }
        
        // detecting elements with more than one face on boundary (these need to be fixed)
        elmt_idx = 0U;
        for ( auto pft=pfverts.begin(); pft!=pfverts.end(); ++pft, ++elmt_idx ) {
              auto etype2D = (HybridElementTypeMesh()==false) ? pelmt[0] : pelmt[elmt_idx];
              if ( !isLineElement( parseFiniteElementTypeEnum( etype2D ) ) )
                {
                   size_t boundaries_per_element(0U);
                   for ( vector<int64_t>::const_iterator pt=(*pft).begin(); pt!=(*pft).end(); ++pt )
                     // the face is on the boundary
                     if ( (*pt) < 0 ) boundaries_per_element++;
                   if ( boundaries_per_element > 1U &&
                       isTriangularElement( parseFiniteElementTypeEnum( etype ) ) )
                     {
                       cerr <<"\n\n\telement "<< elmt_idx <<" ("<< parseFiniteElementType( pelmt[elmt_idx] ) <<") ";
                       cerr <<" has "<< boundaries_per_element <<" faces on model boundary.\n";
                       csmp_error.Note( WARNING, "VData::EstablishElementConnectivity2D:",
                                         "triangular element with  2 faces on boundary ");
                       triangle_with_all_nodes_on_boundary = true;
                     }
                }
            }
        
        // 3. reconnecting line elements
        // -----------------------------
        // map<size_t,set<size_t> >  line_elmt_that_contain_node;
        if ( HybridElementTypeMesh() )
          for ( const auto& it : line_elmt_that_share_node )
            {
                const size_t n_connections(it.second.size()-1);
                
                // 1. isolated line elements terminating either at an inside node (NOT) or at the BOX_BOUNDARY
                // --------------------------------------------------------------------------------------------------
                if ( n_connections == 0U ) {
                     const size_t elmt = (*it.second.begin());
                     // if there is no neighbor element opposite to the zeroeth node
                     if ( it.first == plist[elmt][1] ) {
                          // identifying the boundary that the missing neighbor is located at
                          // (note: although we are not at a BOX_BOUNDARY we need a negative flag, else CSMP infers elmt #0=NOT)
                          pfverts[elmt][0] = (bflags[ it.first ]==0) ? INTERNAL : bflags[ it.first ];
                       }
                     else if ( it.first == plist[elmt][0] ) {
                          // the missing neighbor may be located at a model boundary
                          pfverts[elmt][1] = (bflags[ plist[elmt][0] ]==0) ? INTERNAL : bflags[ it.first ];
                       }
                     else throw csmp::Exception( ERROR, "VData::EstablishElementConnectivity2D", "orphan line element node");
                  }
                // 2. two line elements sharing one node
                // --------------------------------------------------------------------------------------------------
                else if ( n_connections == 1U ) {
                     const size_t elmt1 = (*it.second.begin());
                     const size_t elmt2 = (*it.second.rbegin());
                     // processing the neighbors
                     // line element 1
                     if ( it.first      == plist[elmt1][1] ) pfverts[elmt1][0] = elmt2;
                     else if ( it.first == plist[elmt1][0] ) pfverts[elmt1][1] = elmt2;
                     // line element 2
                     if ( it.first      == plist[elmt2][1] ) pfverts[elmt2][0] = elmt1;
                     else if ( it.first == plist[elmt2][0] ) pfverts[elmt2][1] = elmt1;
                  }
                // 3. line element manifolds (multiple line elements)
                // --------------------------------------------------------------------------------------------------
                // off the possible neighbors, the aligned elements are picked
                else { // n_connections > 1 )
                     // finding the pair of most closely aligned line elements starting at node
                     // establish element combinations
                     vector<int64_t>         joint_line_elmts( it.second.begin(), it.second.end() ); // actual element ids
                     const int64_t           n_elmts_to_combine(2U);
                     deque<vector<int64_t> > combinations;
                     if ( createUniqueCombinations( joint_line_elmts, n_elmts_to_combine, combinations ) == 0 )
                       csmp_error.Note( ERROR, "EstablishElementConnectivity2D", "no combinations between elements available");
                     // finding inter-element angle for all combinations
                     //             angle, combination number
                     vector<pair<double,size_t> > inter_element_angles;
                     inter_element_angles.reserve( combinations.size() );
                     size_t n_combi{0};
                     for ( auto& cit : combinations ) {
                          const double angle = AngleBetweenLineElements2D( cit[0], cit[1] );
                          // ignoring edge direction
                          const double acute_angle = ( angle > 90. ) ? 180. - angle : angle;
                          inter_element_angles.push_back( make_pair( acute_angle, n_combi++ ) );
                       }
                     // sorting the angles to find the edges that are closest to a straight continuation
                     // (= smallest angles for aligned, edges and closest to 180o for ones greater that 90o)
                     sort( inter_element_angles.begin(), inter_element_angles.end(),
                           [](auto& a, auto& b) -> bool { return a.first < b.first; } );
                      // for any 2 edges unique connections are made until there are no more elements to connect
                      set<size_t> assigned_elements;
                      for ( auto& aet : inter_element_angles ) {
                           // connecting the pair of line elements
                           // ------------------------------------
                           const size_t elmt1 = combinations[aet.second][0];
                           const size_t elmt2 = combinations[aet.second][1];
                           // only if both elements in the combination have not been assigned already
                           if ( assigned_elements.find(elmt1) == assigned_elements.end() &&
                                assigned_elements.find(elmt2) == assigned_elements.end() )
                             {
                                // finding the correct side of edge1
                                if      ( it.first == plist[elmt1][1] ) pfverts[elmt1][0] = elmt2;
                                else if ( it.first == plist[elmt1][0] ) pfverts[elmt1][1] = elmt2;
                                assigned_elements.insert( elmt1 );
                                // and edge2
                                if      ( it.first == plist[elmt2][1] ) pfverts[elmt2][0] = elmt1;
                                else if ( it.first == plist[elmt2][0] ) pfverts[elmt2][1] = elmt1;
                                assigned_elements.insert( elmt2 );
                             }
                        }
                      // assigning boundary flag to left-over neighbor elements at manifolds
                      if ( assigned_elements.size() < joint_line_elmts.size() ) {
                           // making sure that there only is a single unassigned element
                           assert( joint_line_elmts.size() - 1 == assigned_elements.size() );
                           // finding the yet-to-be-assigned element
                           size_t unassigned_elmt = numeric_limits<size_t>::max();
                           for ( auto& leit : joint_line_elmts )
                             if ( assigned_elements.find(leit) == assigned_elements.end() ) {
                                  unassigned_elmt = leit;
                                  break;
                               }
                           assert ( unassigned_elmt != numeric_limits<size_t>::max() );
                           // finding the correct side of the line element and assigning the vertex bflag to irt
                           if   ( it.first == plist[unassigned_elmt][1] ) {
                                 pfverts[unassigned_elmt][0] = ( bflags[it.first] < 0 ) ? bflags[it.first] : INTERNAL;
                              }
                           if ( it.first == plist[unassigned_elmt][0] ) {
                                 pfverts[unassigned_elmt][1] = ( bflags[it.first] < 0 ) ? bflags[it.first] : INTERNAL;
                              }
                        }
                  }
                              
            } // processing the line element neighbors

     } // if plist empty

   assert( !pfverts.empty() );

   // 4. Elements with all nodes on the boundary can degrade solver convergence for certain boundary conditions
   //    This method eliminates such elements by node swapping
   if ( triangle_with_all_nodes_on_boundary ) SwitchCornerTriangles2D();

   // 5. reorienting line-element chains (done in other method)
   if ( HybridElementTypeMesh() ) CreateConsistentLineElementOrientations2D();

 } // end EstablishElementConnectivity2D









/**
    Reconnects triangular elements with 3 nodes on the model boundary by switching nodes with their only neighbor; @note needs valid 'pfverts.'
    
    @return the number of replaced triangles;
*/
size_t VData::SwitchCornerTriangles2D()
 {
    if ( HybridElementTypeMesh() ) {
         cout <<"\n\n"<<"VData::SwitchCornerTriangles2D: WARNING: while mesh contains corner-spanning triangles ";
         cout <<" this method cannot be applied because it cannot simultaneously fix the adjacent line elements yet that are contained in the VSet.\n";
         cout <<" Please improve method. or fix corner elements in meshing tool."<< endl;
         return 0U;
      }
      
    size_t switched_triangles{0};
    
    if ( pfverts.empty() || pfverts.size() != plist.size() ) {
         cerr <<"\nVData::SwitchCornerTriangles2D: method needs valid 'pfverts' (neighbor connectivity) for its operation. ";
         cerr <<" no modifications made.\n";
         return 0U;
      }
      
    const size_t n_elements{pfverts.size()};
    for ( size_t elmt_idx{0}; elmt_idx < n_elements; ++elmt_idx )
      {
         if ( isTriangularElement( parseFiniteElementTypeEnum( pelmt[elmt_idx] ) ) ) {
              // count valid neighbors
              size_t nbors = count_if( pfverts[elmt_idx].begin(), pfverts[elmt_idx].end(),
                                       []( int64_t  nbor )->bool { return (nbor >= 0) ? true : false; } );
              assert( nbors != 0 );
              if ( nbors == 1 ) {
                   // finding the only valid neigbor and its shared face
                   pair<size_t,size_t> face_nds;
                   size_t              cnr_nd(UINT_MAX);
                   int64_t               nb_idx(UINT_MAX);
                   for ( auto i{0U}; i<3; ++i ) {
                        if ( pfverts[elmt_idx][i] >= 0 ) {
                             nb_idx = pfverts[elmt_idx][i];
                             switch ( i ) {
                                 case 0: face_nds = make_pair( plist[elmt_idx][1],plist[elmt_idx][2] ); cnr_nd=0; break;
                                 case 1: face_nds = make_pair( plist[elmt_idx][2],plist[elmt_idx][0] ); cnr_nd=1; break;
                                 case 2: face_nds = make_pair( plist[elmt_idx][0],plist[elmt_idx][1] ); cnr_nd=2;
                               }
                             break;
                          }
                     }
                   assert( isTriangularElement( parseFiniteElementTypeEnum( pelmt[nb_idx] ) ) );
                   //           face nodes, other node, but order in pair flipped so that opposite face can be matched
                   map<pair<size_t,size_t>,size_t> fnids{ {{plist[nb_idx][2],plist[nb_idx][1]},0},
                                                          {{plist[nb_idx][0],plist[nb_idx][2]},1},
                                                          {{plist[nb_idx][1],plist[nb_idx][0]},2} };
                   // searching reversed nodes in map
                   auto face_it = fnids.find( face_nds );
                   if ( face_it != fnids.end() ) {
                        // reassigning nodes to new elements (first nodes are taken as those on the far sides (not shared ones))
                        array<int64_t, 3> elmt1_nds{ plist[elmt_idx][0], plist[elmt_idx][1], plist[elmt_idx][2] };
                        array<int64_t, 3> elmt2_nds{ plist[nb_idx][0],   plist[nb_idx][1],   plist[nb_idx][2] };
                        // rotating these vectors such that the unshared nodes become the corner nodes
                        rotate( begin(elmt1_nds), begin(elmt1_nds) + cnr_nd, end(elmt1_nds) );
                        rotate( begin(elmt2_nds), begin(elmt2_nds) + (*face_it).second, end(elmt2_nds) );
                        // new element 1
                        plist[elmt_idx][0] = elmt1_nds[0];
                        plist[elmt_idx][1] = elmt1_nds[1];
                        plist[elmt_idx][2] = elmt2_nds[0];
                        // new element 2 (nb_idx)
                        plist[nb_idx][0]   = elmt2_nds[0];
                        plist[nb_idx][1]   = elmt2_nds[1];
                        plist[nb_idx][2]   = elmt1_nds[0];
                        // reassigning neighbors
                        // old neighbors
                        array<int64_t, 3> elmt1_nbors{ pfverts[elmt_idx][0], pfverts[elmt_idx][1], pfverts[elmt_idx][2] };
                        array<int64_t, 3> elmt2_nbors{ pfverts[nb_idx][0],   pfverts[nb_idx][1],   pfverts[nb_idx][2] };
                        // rotating these vectors such that the unshared nodes become the corner nodes
                        rotate( elmt1_nbors.begin(), elmt1_nbors.begin() + cnr_nd, elmt1_nbors.end() );
                        rotate( elmt2_nbors.begin(), elmt2_nbors.begin() + (*face_it).second, elmt2_nbors.end() );
                        // neighbors of new element 1
                        pfverts[elmt_idx][0] = elmt2_nbors[1];
                        pfverts[elmt_idx][1] = nb_idx;
                        pfverts[elmt_idx][2] = elmt1_nbors[2];
                        // new element 2
                        pfverts[nb_idx][0]   = elmt1_nbors[1];
                        pfverts[nb_idx][1]   = elmt_idx;
                        pfverts[nb_idx][2]   = elmt2_nbors[2];
                        // TODO: check: these new elements should only have a single boundary face left
                        switched_triangles += 2U;
                     }
                }
           }
      }
      
    if ( switched_triangles > 0 )
      cout <<"\n\nVData::SwitchCornerTriangles2D: reconnected: "<< switched_triangles <<" triangular elements."<< endl;
    return switched_triangles;
 
 } // end SwitchCornerTriangles2D









/**
    Rebuilds 3D  'pfverts' from scratch
    
    @attention method relies on correct boundary flags
    
    @todo needs to take into account potentials Faces and Interfaces.
*/
void VData::EstablishElementConnectivity3D()
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( plist.empty() ) {
         csmp_error.Note( WARNING, "VData::EstablishElementConnectivity3D:", "supplied cell vector is empty; nothing was done." );
         return;
      }
    if ( Faces() > 0 ) {
         csmp_error.Note( WARNING, "VData::EstablishElementConnectivity3D:", "Face objects not handled yet" );
         return;
      }
    if ( Interfaces() > 0 ) {
         csmp_error.Note( WARNING, "VData::EstablishElementConnectivity3D:", "InterFace objects not handled yet" );
         return;
      }
 
    // 1. making separate search vectors of face keys for volume, surface and line elements
    // ------------------------------------------------------------------------------------
    cout << "\nEstablishElementConnectivity3D: Establishing CSMP FE neighbor connectivity...\n";
    cout << "  Building a multimap of the faces of the cells...\n";
    // face-node key  element idx, face number
    map<set<size_t>,map<size_t,uint32_t> >  volume_neighbor_keys, surface_neighbor_keys;
    //  node key, elements connected at that node
    map<size_t,set<size_t> >  line_elmt_that_share_node;

    const size_t n_elements(plist.size());
    pfverts.resize( plist.size() );

    for ( size_t elmt_idx{0}; elmt_idx < n_elements; ++elmt_idx )
      {
         // getting the element type (unfortunately this is known only at runtime)
         const auto etype = (HybridElementTypeMesh()==true) ? static_cast<CSMP_FEM_TYPE>(pelmt[elmt_idx]) : static_cast<CSMP_FEM_TYPE>(pelmt[0]);
         assert( parseFiniteElementTypeEnum( pelmt[elmt_idx] ) != UNKNOWN );
         
          const auto faces(CSMP_ElementSpecifications::FacesPerElementOfType(etype));
          pfverts[elmt_idx].resize(faces,IRREGULAR);
          for ( auto face=0U; face<faces; ++face )
            {
               // creating face key of node pointers from indices of face nodes
               set<size_t> key;
               const auto nodes(CSMP_ElementSpecifications::NodesPerFaceForElementOfType( etype, face ) );
               for ( auto j{0U}; j<nodes; ++j ) {
                    // inserting the global  node numbers into the key
                    const size_t face_node = plist[elmt_idx][ CSMP_ElementSpecifications::FaceNodeForElementOfType( etype, face, j ) ];
                    key.insert( face_node );
                 }
               // inserting newly generated keys into map
               if ( CSMP_ElementSpecifications::VolumeElement( etype ) ) {
                    pair<map<set<size_t>,map<size_t,uint32_t> >::iterator,bool>
                      vit = volume_neighbor_keys.insert( make_pair( key, map<size_t,uint32_t>{{elmt_idx,face}} ) );
                    if ( !vit.second ) (*vit.first).second.insert( make_pair(elmt_idx,face) );
                 }
               else if ( CSMP_ElementSpecifications::SurfaceElement( etype ) ) {
                    pair<map<set<size_t>,map<size_t,uint32_t> >::iterator,bool>
                      vit = surface_neighbor_keys.insert( make_pair( key, map<size_t,uint32_t>{{elmt_idx,face}} ) );
                    if ( !vit.second ) (*vit.first).second.insert( make_pair(elmt_idx,face) );
                 }
               else {
                    // for all line elements
                    // node 0
                    pair<map<size_t,set<size_t> >::iterator,bool>
                      it0 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][0], set<size_t>{elmt_idx} ) );
                    // if there is already an entry for the node, the elmt id is added to the set
                    if ( !it0.second ) (*it0.first).second.insert( elmt_idx );
                    // node 1
                    pair<map<size_t,set<size_t> >::iterator,bool>
                      it1 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][1], set<size_t>{elmt_idx} ) );
                    // if there is already an entry for the node, the elmt id is added to the set
                    if ( !it1.second ) (*it1.first).second.insert( elmt_idx );
                 }
            }
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    cout << "  (Re)building neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_elmt_that_share_node.empty() )
      {
        cout << "\n\t\tline elements...\n";
        // map<size_t,set<uint32_t> >  line_elmt_that_contain_node;
        for ( const auto& it : line_elmt_that_share_node )
          {
              const size_t n_connections(it.second.size()-1);
              
              // 1. isolated line elements terminating either at an inside node (INTERNAL) or at the BOX_BOUNDARY
              // --------------------------------------------------------------------------------------------------
              if ( n_connections == 0U ) {
                   const size_t elmt = (*it.second.begin());
                   // if there is no neighbor element opposite the first node
                   if ( it.first == plist[elmt][0] ) {
                        // identifying the boundary that the missing neighbor is located at
                        pfverts[elmt][0] = (bflags[ it.first ]==0) ? INTERNAL : bflags[ it.first ];
                     }
                   else if ( it.first == plist[elmt][1] ) {
                        // the missing neighbor is located at a boundary
                        pfverts[elmt][1] = (bflags[ plist[elmt][1] ]==0) ? INTERNAL : bflags[ plist[elmt][1] ];
                     }
                   else throw csmp::Exception( ERROR, "VData::EstablishElementConnectivity3D", "orphan line element node");
                }
              // 2. two line elements sharing one node
              // --------------------------------------------------------------------------------------------------
              else if ( n_connections == 1U ) {
                   const size_t elmt1 = (*it.second.begin());
                   const size_t elmt2 = (*it.second.rbegin());
                   // processing the neighbors
                   // line element 1
                   if ( it.first      == plist[elmt1][0] ) pfverts[elmt1][0] = elmt2;
                   else if ( it.first == plist[elmt1][1] ) pfverts[elmt1][1] = elmt2;
                   // line element 2
                   if ( it.first      == plist[elmt2][0] ) pfverts[elmt2][0] = elmt1;
                   else if ( it.first == plist[elmt2][1] ) pfverts[elmt2][1] = elmt1;
                }
              // 3. line element manifolds (multiple line elements)
              // --------------------------------------------------------------------------------------------------
              // off the possible neighbors, the aligned elements are picked
              else { // n_connections > 1 )
                   // finding the pair of most closely aligned line elements starting at node
                   // establish element combinations
                   vector<int64_t>         joint_line_elmts( it.second.begin(), it.second.end() ); // actual element ids
                   const int64_t           n_elmts_to_combine(2U);
                   deque<vector<int64_t> > combinations;
                   if ( createUniqueCombinations( joint_line_elmts, n_elmts_to_combine, combinations ) == 0 )
                     csmp_error.Note( ERROR, "EstablishElementConnectivity3D", "no combinations between elements available");
                   // finding inter-element angle for all combinations
                   //             angle, combination number
                   vector<pair<double,size_t> > inter_element_angles;
                   inter_element_angles.reserve( combinations.size() );
                   size_t n_combi{0};
                   for ( auto cit : combinations ) {
                        const double angle = AngleBetweenLineElements3D( cit[0], cit[1] );
                        // ignoring edge direction
                        const double acute_angle = ( angle > 90. ) ? 180. - angle : angle;
                        inter_element_angles.push_back( make_pair( acute_angle, n_combi++ ) );
                     }
                   // sorting the angles to find the edges that are closest to a straight continuation
                   // (= smallest angles for aligned, edges and closest to 180o for ones greater that 90o)
                   sort( inter_element_angles.begin(), inter_element_angles.end(),
                         [](auto& a, auto& b) -> bool { return a.first < b.first; } );
                    // for any 2 edges unique connections are made until there are no more elements to connect
                    set<size_t> assigned_elements;
                    for ( auto& aet : inter_element_angles ) {
                         // connecting the pair of line elements
                         // ------------------------------------
                         const size_t elmt1 = combinations[aet.second][0];
                         const size_t elmt2 = combinations[aet.second][1];
                         // only if both elements in the combination have not been assigned already
                         if ( assigned_elements.find(elmt1) == assigned_elements.end() &&
                              assigned_elements.find(elmt2) == assigned_elements.end() )
                           {
                              // finding the correct side of edge1
                              if      ( it.first == plist[elmt1][0] ) pfverts[elmt1][0] = elmt2;
                              else if ( it.first == plist[elmt1][1] ) pfverts[elmt1][1] = elmt2;
                              assigned_elements.insert( elmt1 );
                              // and edge2
                              if      ( it.first == plist[elmt2][0] ) pfverts[elmt2][0] = elmt1;
                              else if ( it.first == plist[elmt2][1] ) pfverts[elmt2][1] = elmt1;
                              assigned_elements.insert( elmt2 );
                           }
                      }
                    // assigning boundary flag to left-over neighbor elements at manifolds
                    if ( assigned_elements.size() < joint_line_elmts.size() ) {
                         // making sure that there only is a single unassigned element
                         assert( joint_line_elmts.size() - 1 == assigned_elements.size() );
                         // finding the yet-to-be-assigned element
                         size_t unassigned_elmt{UINT_MAX};
                         for ( auto& eit : joint_line_elmts )
                           if ( assigned_elements.find(eit) == assigned_elements.end() ) {
                                unassigned_elmt = eit;
                                break;
                             }
                         assert ( unassigned_elmt != UINT_MAX );
                         // finding the correct side of the line element and assigning the vertex bflag to irt
                         if      ( it.first == plist[unassigned_elmt][0] ) pfverts[unassigned_elmt][0] = bflags[it.first];
                         else if ( it.first == plist[unassigned_elmt][1] ) pfverts[unassigned_elmt][1] = bflags[it.first];
                      }
                }
                            
          } // processed the neighbors of the line elements that share a node
    
       }
       
    // 2.2 surface elements
    // --------------------
    if ( !surface_neighbor_keys.empty() )
      {
        cout << "\n\t\tsurface elements...\n";
        for ( auto it=surface_neighbor_keys.begin(); it!=surface_neighbor_keys.end(); ++it )
          {
             assert( !(*it).second.empty() );
             // 2.2.1 if there is only a single entry, the face is at the model boundary
             // ------------------------------------------------------------------------
             if ( (*it).second.size() == 1U ) {
                  const size_t elmt_idx = (*(*it).second.begin()).first;
                  const uint32_t face   = (*(*it).second.begin()).second;
                  // finding out on which boundary the edge of the face is
                  // relying on appropriate box-boundary flagging
                  // (using only the two first nodes of the face)
                  const CSMP_FEM_TYPE etype = static_cast<CSMP_FEM_TYPE>( pelmt[elmt_idx] );
                  assert( parseFiniteElementTypeEnum( pelmt[elmt_idx] ) != UNKNOWN );
                  const size_t n0 = plist[elmt_idx][CSMP_ElementSpecifications::FaceNodeForElementOfType( etype, face, 0 ) ];
                  const size_t n1 = plist[elmt_idx][CSMP_ElementSpecifications::FaceNodeForElementOfType( etype, face, 1 ) ];
                  pfverts[elmt_idx][face] = whichBoundary( static_cast<BOX_BOUNDARY>( bflags[n0] ),
                                                           static_cast<BOX_BOUNDARY>( bflags[n1] ) );
                }
             // 2.2.2 if there is only one matching neighbor it gets recorded
             // ------------------------------------------------------------------------
             else if ( (*it).second.size() == 2U ) {
                  // the neighbors are recorded in the 'pfverts' map
                  const size_t elmt1      = (*(*it).second.begin()).first;
                  const size_t face_elmt1 = (*(*it).second.begin()).second;
                  const size_t elmt2      = (*next((*it).second.begin(),1)).first;
                  const size_t face_elmt2 = (*next((*it).second.begin(),1)).second;
                  pfverts[elmt1][face_elmt1] = elmt2;
                  pfverts[elmt2][face_elmt2] = elmt1;
                }

             // 2.2.3 if there is a surface element manifold that needs to be disambiguated
             // ----------------------------------------------------------------------------
             else {
                   // finding the pair of most elements with the most closely aligned normals
                   // establish element combinations
                   vector<int64_t>  joint_surf_elmts;
                   joint_surf_elmts.reserve( (*it).second.size() );
                   for ( const auto& i : (*it).second ) joint_surf_elmts.push_back( i.first ); // actual element ids
                   const int64_t           n_elmts_to_combine(2U);
                   deque<vector<int64_t> > combinations;
                   if ( createUniqueCombinations( joint_surf_elmts, n_elmts_to_combine, combinations ) == 0 )
                     csmp_error.Note( ERROR, "EstablishElementConnectivity3D", "no combinations between elements available");
                   // finding inter-element angle for all combinations
                   //             angle, combination number
                   vector<pair<double,size_t> > inter_element_normal_angles;
                   inter_element_normal_angles.reserve( combinations.size() );
                   size_t n_combi{0};
                   for ( auto& cit : combinations ) {
                        const double angle = AngleBetweenSurfaceElements3D( cit[0], cit[1] );
                        // ignoring edge direction
                        const double acute_angle = ( angle > 90. ) ? 180. - angle : angle;
                        inter_element_normal_angles.push_back( make_pair( acute_angle, n_combi++ ) );
                     }
                   // sorting the angles to find the edges that are closest to a straight continuation
                   // (= smallest angles for aligned, edges and closest to 180o for ones greater that 90o)
                   sort( inter_element_normal_angles.begin(), inter_element_normal_angles.end(),
                         [](auto& a, auto& b) -> bool { return a.first < b.first; } );
                    // for any 2 edges unique connections are made until there are no more elements to connect
                    set<size_t> assigned_elements;
                    for ( auto aet : inter_element_normal_angles )
                      {
                         // connecting the pair of surface elements
                         // ---------------------------------------
                         const size_t elmt1 = combinations[aet.second][0];
                         const size_t elmt2 = combinations[aet.second][1];
                         // only if both elements in the combination have not been assigned already
                         if ( assigned_elements.find(elmt1) == assigned_elements.end() &&
                              assigned_elements.find(elmt2) == assigned_elements.end() )
                           {
                              // retrieve the faces of the elements that will be interconnected
                              const size_t nbor_face_e1 = (*(*it).second.find(elmt1)).second;
                              const size_t nbor_face_e2 = (*(*it).second.find(elmt2)).second;
                              pfverts[elmt1][nbor_face_e1] = elmt2;
                              pfverts[elmt2][nbor_face_e2] = elmt1;
                              assigned_elements.insert( elmt1 );
                              assigned_elements.insert( elmt2 );
                           }
                      }
                    // assigning boundary flag to left-over neighbor elements at manifolds
                    if ( assigned_elements.size() < joint_surf_elmts.size() ) {
                         // making sure that there only is a single unassigned element
                         assert( joint_surf_elmts.size() - 1 == assigned_elements.size() );
                         // finding the yet-to-be-assigned element
                         int64_t  unassigned_elmt{-1};
                         for ( auto& i : joint_surf_elmts )
                           if ( assigned_elements.find(i) == assigned_elements.end() ) {
                                unassigned_elmt = i;
                                break;
                             }
                         assert ( unassigned_elmt >= 0 );
                         // finding the correct side of the surface element and assigning a bflag to it
                         const size_t boundary_face = (*(*it).second.find(unassigned_elmt)).second;
                         pfverts[unassigned_elmt][boundary_face] = IRREGULAR;
                      }
               }
          }
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( !volume_neighbor_keys.empty() )
      {
        cout << "\t\tvolume elements...\n";
        for ( auto it=volume_neighbor_keys.begin(); it!=volume_neighbor_keys.end(); ++it )
          {
              // 2.3.0 - there must be at least 1 face entry
              assert( (*it).second.size() >= 0 );
              
              // 2.3.1 interior faces (there is a pair of valid adjacent element faces so that neighbor assignments can be made)
              if ( (*it).second.size() == 2 )
                {
                  // the neighbors are recorded in the 'pfverts' map
                  const size_t elmt1      = (*(*it).second.begin()).first;
                  const size_t face_elmt1 = (*(*it).second.begin()).second;
                  const size_t elmt2      = (*next((*it).second.begin(),1)).first;
                  const size_t face_elmt2 = (*next((*it).second.begin(),1)).second;
                  pfverts[elmt1][face_elmt1] = elmt2;
                  pfverts[elmt2][face_elmt2] = elmt1;
                }
              // 2.3.2 boundary faces
              else if ( (*it).second.size() == 1 ) {
                  const size_t elmt          = (*(*it).second.begin()).first;
                  const size_t boundary_face = (*(*it).second.begin()).second;
                  // TODO: one could narrow down which boundary this is, but it will not be used later
                  pfverts[elmt][boundary_face] = IRREGULAR;
                }
          }

      } // volume elements

    cout << "\n"<<"EstablishElementConnectivity3D: Established neighbor connectivity of "<< pfverts.size() <<" cells successfully.\n"<< endl;

 } // end EstablishElementConnectivity3D





/**
   Replaces corner-spanning tetrahedra (all nodes at the model boundary) and their interior neighbors
   with  3 tedrahedra each to establish the necessary degrees of freedom to assign boundary conditions.
   
   @return the number of new tetrahedra created.
   
   @attention method requires an intact neighbor connectivity. Run EstablishElementConnectivity3D before.
*/
size_t VData::RemeshCornerSpanningTetrahedra()
 {

throw csmp::Exception( ERROR, "VData::RemeshCornerSpanningTetrahedra",
                      "SKM: gave up trying to get this to work, given that there are so many special cases; fix problem when creating the mesh");
 
     // Remeshing tetrahedra that span a corner of the model
     // ----------------------------------------------------------
     // replacing the corner and inner elements with 3 new tetrahedra
     const size_t n_elements{pelmt.size()};
     
     // if this a tetrahedral only mesh no checks have to be performed
     if ( !HybridElementTypeMesh() && isTetrahedral( parseFiniteElementTypeEnum(ElementType(0)) ) )
       {
          for ( size_t i{0U}; i<n_elements; ++i ) {
               long n_neighbors{0}, neighbor(UNSPECIFIED);
               for ( size_t j{0}; j<PfvertsSize(i); ++j )
                 if ( Pfvert(i,j) >= 0 ) {
                      neighbor=Pfvert(i,j);
                      n_neighbors++;
                   }
               // if the tetrahedron has only one neighbor, it and its nbor need to be replaced
               if ( n_neighbors == 1 )
                 splitCornerTetrahedron( *this, i, neighbor );
            }
       }
     else { // for hybrid element type meshes the element type needs to be checked
          for ( size_t i{0U}; i<n_elements; ++i ) {
               long n_neighbors{0}, neighbor(UNSPECIFIED);
               for ( size_t j{0}; j<PfvertsSize(i); ++j )
                 if ( Pfvert(i,j) >= 0 ) {
                      neighbor=Pfvert(i,j);
                      n_neighbors++;
                   }
               // if the element has only one neighbor and is a tetrahedron, it and its nbor need to be replaced
               if ( n_neighbors == 1 &&
                    isTetrahedral( parseFiniteElementTypeEnum(ElementType(i)) ) ) {
                    splitCornerTetrahedron( *this, i, neighbor );
                 }
            }
       }
       
    return Elements() - n_elements;

 } // end RemeshCornerSpanningTetrahedra






/**
      Finds the neighbors of each node and returns them into the argument vector..
      
      This method is equivalent to creating a sparsity pattern for matrix accumulation.
      
      @test OK SKM 8/12/21
*/
void VData::EstablishNodeNeighborConnectivity( vector<set<size_t>>& pnode ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( plist.empty() ) {
         csmp_error.Note( WARNING, "VData::EstablishNodeNeighborConnectivity:", "plist is empty; nothing was done." );
         return;
      }

    if ( !pnode.empty() ) pnode.clear();
    pnode.resize( px.size() );
    
    bool higher_order_elements{false};
 
    const size_t n_elements{ plist.size() };
    // looping over the elements to get the corner nodes of their segments
    for ( size_t elmt{0}; elmt < n_elements; ++elmt ) {
         // getting the element type
         const auto CSMP_FE_type = (HybridElementTypeMesh()) ? static_cast<CSMP_FEM_TYPE>(pelmt[elmt]) : static_cast<CSMP_FEM_TYPE>(pelmt[0]);
         if ( !higher_order_elements && CSMP_ElementSpecifications::InterpolationOrder(CSMP_FE_type) > 1 ) {
              csmp_error.Note( ERROR, "VData::EstablishNodeNeighborConnectivity", "connectivity of midside nodes not tested yet; check!" );
              higher_order_elements = true;
           }
         // for each segment
         const auto n_segments{ CSMP_ElementSpecifications::SegmentsPerElementOfType( CSMP_FE_type ) };
         for ( auto segm_id{0}; segm_id < n_segments; ++segm_id ) {
              pair<size_t,size_t>
                segm_nodes = CSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType( CSMP_FE_type, segm_id );
              // replacing local with global node ids
              segm_nodes.first  = plist[ elmt ][ segm_nodes.first ];
              segm_nodes.second = plist[ elmt ][ segm_nodes.second ];
              // storing the node-to-node connections avoiding duplicates
              pnode[ segm_nodes.first ].insert( segm_nodes.second );
              pnode[ segm_nodes.second ].insert( segm_nodes.first );
           }
         // if this is an interface with nodes on the inside and outside
         if ( elmt > first_interface_ ) {
             // do the whole thing again, for the second lot of node entries in the plist
             const auto n_nodes = CSMP_ElementSpecifications::NodesPerElementOfType(CSMP_FE_type);
             for ( auto segm_id{0}; segm_id < n_segments; ++segm_id ) {
                  pair<size_t,size_t>
                    segm_nodes = CSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType( CSMP_FE_type, segm_id );
                  // replacing local with global node ids
                  segm_nodes.first  = plist[ elmt ][ segm_nodes.first  + n_nodes ];
                  segm_nodes.second = plist[ elmt ][ segm_nodes.second + n_nodes ];
                  // storing the node-to-node connections avoiding duplicates
                  pnode[ segm_nodes.first ].insert( segm_nodes.second );
                  pnode[ segm_nodes.second ].insert( segm_nodes.first );
               }
           }
      }
      
    // finite elements with midside nodes and bubble functions
    if ( higher_order_elements ) {
      for ( size_t elmt{0}; elmt < n_elements; ++elmt )
        if ( pnode[elmt].empty() )
          {
             // getting the element type
             const auto CSMP_FE_type = (HybridElementTypeMesh()) ? pelmt[elmt] : pelmt[0];
             if ( CSMP_ElementSpecifications::InterpolationOrder(CSMP_FE_type) == 2 ) {
                 // dealing with quadratic elements that have midside nodes
                 // relying on the numbering convention that midside nodes follow the corner nodes in the same order
                 // and that there is one midside node per segment
                 const auto n_nodes = CSMP_ElementSpecifications::NodesPerElementOfType(CSMP_FE_type);
                 // for each segment
                 const auto n_segments{ CSMP_ElementSpecifications::SegmentsPerElementOfType( CSMP_FE_type ) };
                 const auto first_midside_node = n_nodes - n_segments;
                 // assigning the segment corner nodes as neighbors of the midside node
                 for ( auto segm_id{0}; segm_id < n_segments; ++segm_id ) {
                      pair<size_t,size_t>
                        segm_nodes = CSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType( CSMP_FE_type, segm_id );
                      // replacing local with global node ids
                      const size_t midside_node = plist[ elmt ][ first_midside_node + segm_id ];
                      segm_nodes.first          = plist[ elmt ][ segm_nodes.first ];
                      segm_nodes.second         = plist[ elmt ][ segm_nodes.second ];
                      // storing the node-to-node connections avoiding duplicates
                      pnode[ midside_node ].insert( segm_nodes.first );
                      pnode[ midside_node ].insert( segm_nodes.second );
                   }
               }
             else {
                  csmp_error.Note( ERROR, "VData::EstablishNodeNeighborConnectivity", "connectivity of midside nodes for O>2 meshes not done yet" );
                  break;
               }
           }
           
      } // end higher-order elements

#ifdef MESH_MANAGER_DEBUG
    // printing the node-neighbor vector for testing
    cout <<"\n\nVData::EstablishNodeNeighborConnectivity: connectivity created for "<< pnode.size() <<" nodes:";
    size_t node{0};
    for ( auto nit : pnode ) {
         cout <<"\n\t" << node <<": ";
         for ( auto i : nit ) cout << i <<" ";
         cout <<" ("<< parseBoundary( intToBOX_BOUNDARY( bflags[node] ) ) <<")";
         node++;
      }
#endif

 } // end EstablishNodeNeighborConnectivity






/* DEPRECATED

    vertex manifolds: pairs of nodes and their INSIDE,OUTSIDE, MIDDLE classifers
    typedef deque<set<pair<size_t,int8_t> > > vertexManifoldIndices;
    
    Checks for collocated vertices into transfer data structure.

size_t VData::ExtractNodeManifolds( vertexManifoldIndices& indexes ) const
 {
    // 0.  if there are no interfaces, method returns false
    if (  first_interface_ == first_face_ ) {
         cerr <<"\nVData::ExtractNodeManifolds: does not contain any manifolds.\n";
         return 0;
      }

    // 1. creating some temporary containers
    // 1.1 all manifold vertices encountered with their side classifiers
    map<size_t,INTERFACE_SIDE>  manifold_vertex_classifications;
    // 1.2 all collocated manifold vertices are stored in sets
    set<pair<size_t,size_t> >   manifold_vertices;
      
   // 2. creating manifold data from the interface node indices stored in plist
   for ( auto it=PlistInterFacesBegin(); it!=PlistInterFacesEnd(); ++it )
     {
        // the first half of the interface vertices represent inside nodes in CCW order outside looking in
        // the second half of the indices represent the outside nodes in matching albeit reverse order
        assert( !(*it).empty() );
        const auto plist_entries((*it).size());
        assert( (plist_entries & 1) == 0  ); // even number
        const auto iface_nodes(plist_entries/2);
        
        // creating a new entry in the manifold map or getting an iterator to an existing one
        // map<size_t,set<pair<size_t,int8_t> > > vertexManifoldIndices
        for ( auto i{0U}; i<iface_nodes; ++i ) {
              // recording which vertices are collocated
              pair<size_t,size_t> vertex_pair = { static_cast<size_t>((*it)[i]), static_cast<size_t>((*it)[plist_entries-i-1]) };
              if ( vertex_pair.first > vertex_pair.second )
                swap( vertex_pair.first, vertex_pair.second );
              manifold_vertices.insert( vertex_pair );
              // collecting the manifold vertices together with their classifiers
              manifold_vertex_classifications.insert( make_pair( (*it)[i], INSIDE ) );
              manifold_vertex_classifications.insert( make_pair( (*it)[plist_entries-i-1], OUTSIDE ) );
          }
      }
      
   // 3. matching the InterFace vertices to potential intervening elements
   size_t interface{Elements()+Faces()};
   for ( auto it=next(PfvertsBegin(),first_interface_); it!=PfvertsEnd(); ++it ) {
         assert( (*it).size() > 6U ); // 6=minimum nbors, +2 higher-dim elmts +2 face_ids
         // the last entry in each pfvert record is the index of the intervening element or bflag
         const auto elmt_idx = (*it).back();
         if ( elmt_idx >= 0 ) {// if there is an intervening element
             auto nodes_of_elmt{ plist[interface].size() / 2 };
             // the nodes of the intervening element match (in position and order) those of the INSIDE of the Interface
             for ( auto i{0U}; i<nodes_of_elmt; ++i ) {
                  manifold_vertex_classifications.insert( make_pair( plist[elmt_idx][i], MIDDLE ) );
                  pair<size_t,size_t> vertex_pair = { plist[interface][i], plist[elmt_idx][i] };
                  if ( vertex_pair.first > vertex_pair.second ) swap( vertex_pair.first, vertex_pair.second );
                  manifold_vertices.insert( vertex_pair );
                }
           }
         interface++;
     }
     
   // 4. consolidating the manifold vertex pairs into singles, pairs, and multiples needed for later classification
   map<size_t,set<size_t> > manifold_node_clusters;
   for ( const auto& it : manifold_vertices ) {
        // adding all permutations of indices into map while elimating duplicates at the same time
        auto insert_it = manifold_node_clusters.insert( make_pair( it.first, set<size_t>{ it.first, it.second } ) );
        if ( !insert_it.second )
          (*insert_it.first).second.insert( it.second );
     }
   // up to here we only get a maximum of three entries before the original pairs do not know of each other
   // this is fixed in the following second pass
   for ( const auto& it : manifold_vertices )
     // if the second number of the pair is not a cluster key
     if ( manifold_node_clusters.count( it.second ) == 0 )
       {  // the first number is searched for in the cluster sets (it must be in one of them)
          for ( auto& cit : manifold_node_clusters )
            if ( cit.second.find( it.first ) != cit.second.end() ) {
                 // if it is found, the missing number 2 is inserted
                 cit.second.insert( it.second );
                 break;
              }
       }
       
   // 5. populating the output datastructure: map<size_t,set<pair<size_t,INTERFACE_SIDE> > >
   if ( !indexes.empty() ) indexes.clear();
   for ( auto& mit : manifold_node_clusters ) {
        set<size_t>  vertices;
        // for all the stored vertices
        for ( auto& sit : mit.second ) {
              auto attribute_it = manifold_vertex_classifications.find( sit );
              // the vertex must be present in map
              assert( attribute_it != manifold_vertex_classifications.end() );
              //vertices_with_attributes.insert( make_pair( sit, (*attribute_it).second ) );
          }
        indexes.insert( make_pair( mit.first, vertices ) );
     }
   
   return manifold_node_clusters.size();
   
 } // end ExtractNodeManifolds

*/





 /**
     Initialises the TOPOTYPE flags, indicating the role that the nodes play in defining the Boundary Representation (BREP) of the model,
     i.e., whether they represent essential crossing or corner points,
     define lines or surfaces on the inside or the outer model boundaries,
     or whether they are just mesh vertices that can be changed during remeshing.
     
     Depending on whether the VSet contains only 'Element' objects and 'Faces' or whether it already contains split mesh patches
     joint by 'InterFace' objects and therefore NodeManifolds, the geometry flags are initialised differently.
     
     All previous flags are overwritten.
     
     @attention this must done before line and surface elements dividing the volume mesh into different regions are removed but AFTER the BOX_BOUNDARY flags have been assigned.
     
     @attention successful applicaton of the method requires and intact element connectivity with disambiguated manifolds of lower dimensional elements.
     
     TODO: not tested yet, use Anne-Laure's VSets to perform testing
*/
void VData::InitialiseNodeTopologyIdentifiers()
 {
    const uint32_t dim = SpatialDimension();
    
    // 1. resetting all TOPOTYPE flags of the nodes in the mesh to a default chosen as MESH_VERTEX
    // -------------------------------------------------------------------------------------------
    for ( auto nit=BREP_FlagsBegin(); nit!=BREP_FlagsEnd(); ++nit )
      (*nit) = MESH_VERTEX; // default
      
    // 2. processing nodes on the basis of BOX_BOUNDARY information
    // ------------------------------------------------------------
    for ( size_t i{0U}; i<Vertices(); ++i ) {
         const BOX_BOUNDARY flag{ BFlag(i) };
         // nodes on the exterior of the model
         if ( flag != NOT ) {
              if ( dim == 1U ) {
                   // opposite boundaries depending whether 1D model is horizontal or vertical
                   if ( (!belongsToSide(LEFT,flag) || !belongsToSide(RIGHT,flag)) ||
                        (!belongsToSide(BOTTOM,flag) || !belongsToSide(TOP,flag)) )
                     BREP_Flag( i, INTERSECTION_POINT );
                   else if ( isCorner( flag ) )
                     BREP_Flag( i, EXTERIOR_POINT );
                }
              else if ( dim == 2U ) {
                   if ( isCorner( flag ) ) BREP_Flag( i, EXTERIOR_POINT );
                   else if ( isSide( flag ) )
                     BREP_Flag( i, EXTERIOR_LINE );
                }
              else if ( dim == 3U ) {
                   if ( isCorner( flag ) ) BREP_Flag( i, EXTERIOR_POINT );
                   else if ( isEdge( flag ) )
                     BREP_Flag( i, EXTERIOR_LINE );
                   else if ( isSide( flag ) )
                     BREP_Flag( i, EXTERIOR_SURFACE );
                }
              // internal boundaries (defaults; perimeter points/lines/surfaces are considered later)
              if ( flag == INTERNAL ) {
                   if      ( dim == 1U ) BREP_Flag( i, PERIMETER_POINT );
                   else if ( dim == 2U ) BREP_Flag( i, INTERIOR_LINE );
                   else if ( dim == 3U ) BREP_Flag( i, INTERIOR_SURFACE );
                }
              // else MESH_VERTEX flagging remains
           }
      }

    // 3. If InterFace objects are already present the stored node manifold information in the model is used
    // -----------------------------------------------------------------------------------------------------
    if ( Interfaces() > 0 ) {
         assert( NodeManifolds() > 0 );
         // 3.1 looping over the node manifolds to identify crossing points
         for ( auto nit=PmanifoldsBegin(); nit!=PmanifoldsEnd(); ++nit  ) {
              // if the manifolds are located at the model boundary, they have already been dealt with
              bool at_outer_boundary{false};
              for ( const auto& n : (*nit).first )
                if ( BFlag(n) != NOT && BFlag(n) != INTERNAL ) {
                     at_outer_boundary = true;
                     break;
                  }
              // processing the manifolds on the inside of the model
              if ( !at_outer_boundary ) {
                   // crossings
                   if ( (*nit).second == ManifoldType::SPLIT_BOUNDARY_CROSSING ||
                        (*nit).second == ManifoldType::MULTI_SB_CROSSING ||
                        (*nit).second == ManifoldType::SPLIT_BOUNDARY_TERMINATION ) // T-intersection
                     for ( const auto& n : (*nit).first )
                       BREP_Flag( n, INTERSECTION_POINT );
                   // intersection lines (interior lines) already handled
                }
           }
         // 3.2 dealing with perimeter nodes on SplitBoundaries
         for ( auto it=PlistInterFacesBegin(); it!=PlistInterFacesEnd(); ++it ) {
               // flagging those nodes which are shared between the inside and the outside of the interfaces
               const auto n_inside_nodes{ (*it).size() / 2 };
               for ( auto i{0U}; i<n_inside_nodes; ++i )
                 // node that outside nodes are in reverse order!
                 if ( (*it)[i] == (*it)[(*it).size()-1-i] ) {
                      if     ( dim == 2U ) BREP_Flag( (*it)[i], PERIMETER_POINT );
                      else if( dim == 3U ) BREP_Flag( (*it)[i], PERIMETER_LINE );
                   }
           }
         return;
      }
      

    // 4. If there are only Elements or Elements + Faces, meaning that the nodes are not duplicated,
    //    the lower-dimensional finite elements that the node forms part of are recorded in maps
    // ------------------------------------------------------------------------------------------
    CSMP_FEM_TYPE   elmt_type = ( HybridElementTypeMesh() ) ? UNKNOWN : static_cast<CSMP_FEM_TYPE>(ElementType(0U));
    //  node       elmt or face
    map<size_t,set<size_t> > node_parent_line_elmts;
    map<size_t,set<size_t> > node_parent_surf_elmts;
    
    for ( size_t i{0U}; i<TotalNumberOfCells(); i++ )
      {
         // getting the element type
         if ( HybridElementTypeMesh() ) elmt_type = static_cast<CSMP_FEM_TYPE>(ElementType(i));
         // line elements
         if ( isLineElement( elmt_type ) ) {
              // finding end-points and points somewhere on a polyline
              // -----------------------------------------------------
              // getting the indices of the corner nodes
              const uint32_t n0 = static_cast<uint32_t>(Plist( i, 0 ));
              const uint32_t n1 = static_cast<uint32_t>(Plist( i, 1 ));
              // recording parent elements of nodes on INTERNAL boundaries
              auto nit0 = node_parent_line_elmts.insert( make_pair( n0, set<size_t>{i} ) );
              if ( nit0.second == false ) (*nit0.first).second.insert( i );
              auto nit1 = node_parent_line_elmts.insert( make_pair( n1, set<size_t>{i} ) );
              if ( nit1.second == false ) (*nit1.first).second.insert( i );
           }
         // surface elements
         else if ( dim == 3U && isSurfaceElement( elmt_type ) ) {
             // creating the map
             for ( auto j{0U}; j<PlistSize(i); j++ ) {
                  const int64_t node{ Plist( i, j ) };
                  auto nit = node_parent_surf_elmts.insert( make_pair( node, set<size_t>{i} ) );
                  if ( nit.second == false ) (*nit.first).second.insert( (i) );
               }
          }
      }
    
    
    // 4. classification of BREP type based on the line elements discovered
    // --------------------------------------------------------------------
   //     finding crossing points of two or multiple polylines
   // (assuming that these are flagged INTERNAL
   for ( const auto& nit : node_parent_line_elmts )
     {
        // processing the nodes of the lower dimensional elements to discern intersection and perimeter points
        for ( auto& ne : nit.second ) {
              // getting the node Idx of the corner nodes of the line element
              const size_t n0 = static_cast<size_t>(Plist( ne, 0 ));
              const size_t n1 = static_cast<size_t>(Plist( ne, 1 ));
              // if the line element has no neighbor opposite the key node, the node is located on the perimeter point of an internal boundary
              if ( nit.first == n0 && BFlag(n0) == INTERNAL ) {
                   // if there is no neighbor line element
                   if ( Pfvert( ne, 1 ) < 0 )
                     BREP_Flag( nit.first, PERIMETER_POINT );
                }
              if ( nit.first == n1 && BFlag(n1) == INTERNAL ) {
                   if ( Pfvert( ne, 0 ) < 0 )
                     BREP_Flag( nit.first, PERIMETER_POINT );
                }
              // intersections of internal line elements with the model boundary
              // need to be flagged EXTERIOR_POINT
              if ( nit.first == n0 && BFlag(n0) != INTERNAL && BFlag(n0) != NOT ) {
                   // if there is no neighbor line element
                   if ( Pfvert( ne, 1 ) < 0 )
                     BREP_Flag( nit.first, EXTERIOR_POINT );
                }
              if ( nit.first == n1 && BFlag(n1) != INTERNAL && BFlag(n0) != NOT ) {
                   if ( Pfvert( ne, 0 ) < 0 )
                     BREP_Flag( nit.first, EXTERIOR_POINT );
                }
          }

         // T-intersections are treated the same as crossing lines
         size_t n_connected_lines{ nit.second.size() };
         //    T-intersection             lines crossing
         if (  n_connected_lines == 3U || n_connected_lines >= 4U ) {
              if ( BFlag(nit.first) == INTERNAL )
                BREP_Flag( nit.first, INTERSECTION_POINT );
           }
         // MULTIPLE_INTERSECTIONS else if ( nit.second.size() >= 5U )
     }
       
       
   // 5. processing surface-related topology in three-dimensional models:
   //    - perimeter lines of internal surfaces
   //    - intersection lines of internal surfaces
   //    - touching points of internal surfaces
   //    - intersection points of surface and curves made of line elements
   //      (not for lines that terminate at surfaces because their end-points are already flagged)
   // --------------------------------------------------------------------------------------------
   // TODO: test with a predefined 3D model (SKUA?)
   if ( dim == 3U ) {
       for ( const auto& nit : node_parent_surf_elmts )
         // internal surfaces
         if ( BFlag(nit.first) == INTERNAL || BFlag(nit.first) == NOT )
           {
              // recording the lower-dimensional parent elements of the node
              set<CSMP_FEM_TYPE> etypes;
              for ( const auto& eit : nit.second ) {
                   const CSMP_FEM_TYPE etype = static_cast<CSMP_FEM_TYPE>(ElementType(eit));
                   if ( !isVolumeElement(etype) )
                     etypes.insert( etype );
                }
              // if there are only surface elements
              if ( etypes.size() == 1U  ) {
                  assert( isSurfaceElement(*etypes.begin()) );
                  // 0. case INTERIOR_SURFACE where node lies in the middle of an interior surface was already covered above
                  
                  // 1. case PERIMETER_LINE defined by nodes on the perimeter curve of an internal surface
                  // -------------------------------------------------------------------------------------
                  // (criterion: >2 of the surface-elmt edges that the node is part of, must have no neighbor)
                  int edge_elmt_count{0};
                  // 1.1 finding the (neighbor-free) edges of the parent surface elements that the node is part of
                  for ( const auto& eit : nit.second ) {
                       // looping over each elements pfverts
                       for ( auto face{0U}; face<pfverts[eit].size(); face++ )
                         // if the face has no neighbor we check whether the node is a corner node of it
                         if ( pfverts[eit][face] < 0 ) {
                             uint32_t fn0 = CSMP_ElementSpecifications::FaceNodeForElementOfType( ElementType(eit), face, 0U );
                             uint32_t fn1 = CSMP_ElementSpecifications::FaceNodeForElementOfType( ElementType(eit), face, 1U );
                             if ( fn0 == nit.first || fn1 == nit.first )
                               edge_elmt_count++;
                           }
                    }
                  if ( edge_elmt_count >= 2 )
                    BREP_Flag( nit.first, PERIMETER_LINE );
                    
                  // 2. case INTERSECTION_LINE along which multiple surfaces intersect
                  // -----------------------------------------------------------------
                  // (an internal node on a quad-only surface has 4-quad parents)
                  if ( isQuadrilateral(*etypes.begin()) && nit.second.size() > 4U )
                    BREP_Flag( nit.first, INTERSECTION_LINE );
                  // for triangular elements this criterion does not work
                  // so that element orientations must be considered
                  else if ( nit.second.size() > 3U ) {
                       // finding all possible combinations of surface elements
                       deque<vector<size_t> > combinations;
                       const size_t           n_elmts_to_combine{2};
                       vector<size_t>         surf_elmts( nit.second.begin(), nit.second.end() );
                       createUniqueCombinations( surf_elmts, n_elmts_to_combine, combinations );
                       // measuring the angles
                       set<double> angle_between_surf_elmts;
                       for ( const auto& it : combinations ) {
                       angle_between_surf_elmts.insert( AngleBetweenSurfaceElements3D( it[0], it[1] ) );
                       // ideally, a bi-modal distribution would indicate an intersection line
                       // TODO: how can this be detected?
                       // here we use the criterion that 2 element interangles must be greater than 45o
                       int angles_greater45deg{0};
                       for ( const auto& ia : angle_between_surf_elmts )
                         if ( ia > 45. ) angles_greater45deg++;
                       if ( angles_greater45deg >= 2 )
                         BREP_Flag( nit.first, INTERSECTION_LINE );
                    }
                }
           } // end etypes == 1
         
         // 3. case INTERSECTION_POINT between internal curves and surfaces
         //----------------------------------------------------------------
         if ( etypes.size() == 2U ) {
              // nodes that are FLAGGED INTERNAL but not PERIMETER_POINT will become INTERSECTION_POINT
              if ( BREP_Flag( nit.first) != PERIMETER_POINT )
                BREP_Flag( nit.first, INTERSECTION_POINT );
           }
              
       } // end nodes on internal boundaries

     // 4. case EXTERIOR_LINE or EXTERIOR_POINT where internal surfaces touch model boundary
     //-------------------------------------------------------------------------------------
     for ( const auto& nit : node_parent_surf_elmts )
       // for nodes on the model exterior that are connected to internal surfaces
       if ( BFlag(nit.first) != INTERNAL && BFlag(nit.first) != NOT )
         {
            int BREP_flags_assignments_made{0};
            // criteria: only surface elements with nodes on model interior are considered
            // - these are perimeter points if only this node of the parent surface is located on the model boundary
            // - and perimeter lines if two nodes located on the p
            for ( const auto& eit : nit.second ) {
                 // does the element have interior nodes?
                 int elmt_interior_nodes{false};
                 for ( const auto& n : plist[eit] )
                   if ( bflags[n] == NOT || bflags[n] == INTERNAL )
                     elmt_interior_nodes++;
                 // if so, the element gets considered
                 if ( elmt_interior_nodes > 0 ) {
                     for ( auto face{0U}; face<pfverts[eit].size(); face++ )
                       // if the face has no neighbor we check whether the node is a corner node of it
                       if ( pfverts[eit][face] < 0 ) {
                           uint32_t fn0 = CSMP_ElementSpecifications::FaceNodeForElementOfType( ElementType(eit), face, 0U );
                           uint32_t fn1 = CSMP_ElementSpecifications::FaceNodeForElementOfType( ElementType(eit), face, 1U );
                           if ( fn0 == nit.first || fn1 == nit.first ) {
                                BREP_Flag( nit.first, EXTERIOR_LINE );
                                BREP_flags_assignments_made++;
                             }
                         }
                   }
              }
           // if none of the surface elements that the node is part of has an edge on the model boundary,
           // this must be a touching point
           if ( BREP_flags_assignments_made == 0 )
             BREP_Flag( nit.first, EXTERIOR_POINT );
        }
    } // if dim == 3
    
 } // end InitialiseNodeTopologyIdentifiers




/**
Assuming the following numbering:

First tetra
0 1 2 3 - ony one neighbor = nbr

Second tetra
0 4 1 2 -  4 neigbors (inside of model)

The new elements are formed
1)  0 1 4 3 -   neighbors:  2  3  old-tet1-nbor 1   old-tet1-nbor 4
2)  1 2 4 3 -   neighbors:  3  1  old-tet1-nbor 0   old-tet2-nbor 0
3)  0 4 2 3 -   neighbors:  2  old-tet1-nbor 1   1   old-tet2-nbor 3

The corresponding neighbor elements are.

The operations will increase the number of elements by 1 for each corner processed
because 2 elements are replaced by 3.

@attention method returns without doing anything when all the nodes of the 2 tets are at the model boundary
because, in that case, the connecting line between the corner and the extreme node of the inside tetrahedron will not go hrough the tetrahedra.
As a consequence, degenerate tetrahedra would be produced by the splitting operation.

@revision 2, SKM 27/1/2022

@todo does not work yet; too many potential configurations of the neighbor element have to be handled

@todo does not repair any surface elements that are covering the faces of the tetrahedra and would need to be adjusted

*/
void splitCornerTetrahedron( VData& vdata, size_t cnr, size_t nbr )
 {
    assert( cnr < vdata.Elements() );
    assert( nbr < vdata.Elements() );
    
    // getting the volume of the corner element for reference
    array<double,3>  bbox = boundingBox( vdata, cnr );
    const double cnr_vol = bbox[0] * bbox[1] * bbox[2];

    // 1. Establishing the order of the 4 nodes that will be reconnected
    // -----------------------------------------------------------------
    const size_t n_nodes{ vdata.PlistSize(cnr) };
    size_t       n_cnr(UINT_MAX);
    assert( n_nodes == 4 ); // tetrahedra only
    for ( int i{0}; i<n_nodes; ++i )
      if ( vdata.Pfvert(cnr,i) >= 0 ) { // corner element has only one neighbor; all are opposite to nodes
           n_cnr = i;
           break;
        }
    // getting the remaining nodes, numbered in ascending order (assuming that cnr node is 3)
    long n0 = ( n_cnr+1 > 3 ) ? 0 : n_cnr+1;
    long n1 = ( n0   +1 > 3 ) ? 0 : n0+1;
    long n2 = ( n1   +1 > 3 ) ? 0 : n1+1;
    // remembering the neighbors of the original elements
    const long e1n0(vdata.Pfvert(cnr,n0)), e1n1(vdata.Pfvert(cnr,n1)), e1n2(vdata.Pfvert(cnr,n2)), e1n3(vdata.Pfvert(cnr,n_cnr));
    // converting local node numbers to global node numbers
    n0    = vdata.Plist(cnr,n0);
    n1    = vdata.Plist(cnr,n1);
    n2    = vdata.Plist(cnr,n2);
    n_cnr = vdata.Plist(cnr,n_cnr);
    // establishing extra node of opposite tetrahedron 'nbor' (more robust to rely on node numbers only)
    size_t n4{UINT_MAX};
    for ( auto i{0U}; i<n_nodes; ++i )
      if ( vdata.Plist(nbr,i) != n0 &&
           vdata.Plist(nbr,i) != n1 &&
           vdata.Plist(nbr,i) != n2 &&
           vdata.Plist(nbr,i) != n_cnr ) {
           n4 = i;
           break;
        }
    const long e2n3(vdata.Pfvert(nbr,n4) );
    n4 = vdata.Plist(nbr,n4);
    // determining the nodes of the original neighbor element
    //if ( vdata.Plist(nbr,0) == n0
    // remembering the neighbors of the original elements
    //const long e2n0(vdata.Pfvert(nbr,)), e2n1(vdata.Pfvert(nbr,), e2n2(vdata.Pfvert(nbr,), e2n3(vdata.Pfvert(nbr,n4);


    // 2. dealing with some of the special cases where node 4 also is on the outside of the model
    // ------------------------------------------------------------------------------------------
    // if the node of the nbr element that is not shared with 'cnr' also is on the model boundary, the tetra is split into two
    // no new elements have to be created
    if ( vdata.BFlag(n4) < 0 ) {
// debugging
elementToVTK( vdata, cnr, "a_corner_tetrahedron" );
elementToVTK( vdata, nbr, "a_interior_tetrahedron" );
         // new first element (using storage of former corner element)
         // TODO: depedendent on the case, one of the two following options will work, which one?
         // TODO: perhaps this can be resolved if one establishes first which of the 5 nodes lie in the same plane (make this node 1)
         // works on bottom corner case
         vdata.Plist( cnr, 0, n0 ); // <-
         vdata.Plist( cnr, 1, n1 );
         vdata.Plist( cnr, 2, n4 );
         vdata.Plist( cnr, 3, n_cnr );
         bbox = boundingBox( vdata, cnr );
         // works on top corner case
         // here we check whether the first attempt lead to degenerate element with no volume
         if ( bbox[0]*bbox[1]*bbox[2] < cnr_vol* 1.0e-3 ) {
              vdata.Plist( cnr, 0, n2 ); // <-
              vdata.Plist( cnr, 1, n_cnr );
              vdata.Plist( cnr, 2, n4 );
              vdata.Plist( cnr, 3, n0 );
           }
         // neighbors - TODO: may be in this order we are overwriting some neighbors by accident
         vdata.Pfvert( cnr, 0, nbr );
         vdata.Pfvert( cnr, 1, IRREGULAR_OUTSIDE );
         vdata.Pfvert( cnr, 2, vdata.Pfvert( nbr, 2 ) );
         vdata.Pfvert( cnr, 3, IRREGULAR_OUTSIDE ); // TODO?
         // always OK: new second element (using storage of former inside element)
         vdata.Plist( nbr, 0, n1 );
         vdata.Plist( nbr, 1, n4 );
         vdata.Plist( nbr, 2, n_cnr );
         vdata.Plist( nbr, 3, n2 );
         // neighbors
         vdata.Pfvert( nbr, 0, IRREGULAR_OUTSIDE );
         vdata.Pfvert( nbr, 1, cnr );
         vdata.Pfvert( nbr, 2, vdata.Pfvert( nbr, 0 ) );
         vdata.Pfvert( nbr, 3, vdata.Pfvert( cnr, 0 ) );
// debugging
elementToVTK( vdata, cnr, "b_corner_tetrahedron" );
elementToVTK( vdata, nbr, "b_interior_tetrahedron" );
         return;
      }


    // 2. creation of 3 tetrahedra from the corner and inner ones
    // ----------------------------------------------------------
    // 2.1 creating space for one new element (the other 2 will be reused)
    // -------------------------------------------------------------------
    // (and assuming that there are no Faces or InterFaces,
    //  or lower-dimensional elements covering the sides of the tetrahedra)
    assert( vdata.Faces() == 0 );
    assert( vdata.Interfaces() == 0 );
    const size_t n_elements_new{ vdata.Elements() + 1 };
    if ( vdata.HybridElementTypeMesh() ) {
         assert( isTetrahedral( parseFiniteElementTypeEnum(vdata.ElementType(cnr)) ) );
         assert( isTetrahedral( parseFiniteElementTypeEnum(vdata.ElementType(nbr)) ) );
         vdata.ResizeElementTypes( n_elements_new );
         vdata.ElementType( n_elements_new-1U, vdata.ElementType( cnr ) );
      }
    else assert( isTetrahedral( parseFiniteElementTypeEnum(vdata.ElementType(0)) ) );
    
    vdata.ResizePlist( n_elements_new, 4 ); // nodes of tetrahedron
    vdata.ResizePfverts( n_elements_new, 4 ); // nbors of tetrahedron
    
    // 2.2. assignment of nodes (see Pain et al., 2001)
    // ------------------------------------------------
    // 4 nodes involved
    // first new element
    vdata.Plist( n_elements_new-1U, 0, n0 ); // 0
    vdata.Plist( n_elements_new-1U, 1, n1 ); // 1
    vdata.Plist( n_elements_new-1U, 2, n4 ); // 4 (belongs to inside nbor tetra)
    vdata.Plist( n_elements_new-1U, 3, n_cnr ); // 3
    // new second element (using storage of former corner element)
    vdata.Plist( cnr, 0, n1 );
    vdata.Plist( cnr, 1, n2 );
    vdata.Plist( cnr, 2, n4 );
    vdata.Plist( cnr, 3, n_cnr );
    // new third element (using storage of former inside element)
    vdata.Plist( nbr, 0, n0 );
    vdata.Plist( nbr, 1, n4 );
    vdata.Plist( nbr, 2, n2 );
    vdata.Plist( nbr, 3, n_cnr );

    // 3. assignment of neigbors
    // -------------------------
    const long new1(n_elements_new-1U), new2(cnr), new3(nbr),
               old1n0(vdata.Pfvert(cnr,0)), old1n1(vdata.Pfvert(cnr,1)),
               old2n0(vdata.Pfvert(nbr,0)), old2n1(vdata.Pfvert(nbr,1)), old2n3(vdata.Pfvert(nbr,3));
               
    // new element first
    vdata.Pfvert( new1, 0, new2 );
    vdata.Pfvert( new1, 1, new3 );
    vdata.Pfvert( new1, 2, old1n1 ); // not sure
    vdata.Pfvert( new1, 3, old2n3 ); // not sure

    // second element neighbors
    vdata.Pfvert( cnr, 0, new3 );
    vdata.Pfvert( cnr, 1, new1 );
    vdata.Pfvert( cnr, 2, old1n0 ); // not sure
    vdata.Pfvert( cnr, 3, old2n0 ); // not sure

    // third element neighbors
    vdata.Pfvert( nbr, 0, new2 );
    vdata.Pfvert( nbr, 1, old1n1 ); // not sure
    vdata.Pfvert( nbr, 2, new1 );
    vdata.Pfvert( nbr, 3, old2n1 ); // not sure
    
    cout <<"\nsplitCornerTetrahedron: replaced corner "<< cnr;
    cout <<" and its neighbor "<< nbr <<", adding the new tetrahedron "<< new1;

 } // end splitCornerTetrahedron





/// returns the 3D bounding box of the element
array<double,3>  boundingBox( const VData& vdata, size_t elmt )
 {
    //assert( distance(vdata.Pz) > 0 );
    set<double> dx{vdata.Px(vdata.Plist(elmt,0))}, dy{vdata.Py(vdata.Plist(elmt,0))}, dz{vdata.Pz(vdata.Plist(elmt,0))};
    
    for ( auto i{1}; i<vdata.PlistSize(elmt); ++i ) {
         dx.insert( vdata.Px(vdata.Plist(elmt,i)) );
         dy.insert( vdata.Py(vdata.Plist(elmt,i)) );
         dz.insert( vdata.Pz(vdata.Plist(elmt,i)) );
      }
      
    // determining the ranges
    return array<double,3>{ (*dx.rbegin()) - (*dx.begin()), (*dy.rbegin()) - (*dy.begin()), (*dz.rbegin()) - (*dz.begin()) };
 
 } // end boundingBox





/**
       Node numbers, faces etc.
       
       @attention only works for tetrahedra so far
*/
void elementToVTK( const VData& vdata, size_t eidx, const char* outfile )
 {
     // 0. opening data output file in ascii format
     string file_name(outfile);
     file_name += to_string(eidx);
     file_name += ".vtk";

     ofstream ofs( file_name, ios::out|ios::trunc );
     if ( !ofs ) {
           cerr <<"\nelementToVTK(tetrahedron): '"<< outfile <<"'";
           cerr <<"output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // -------------------------
     const string var_name("node_number");
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     const size_t npe = distance(vdata.PlistBegin(eidx),vdata.PlistEnd(eidx));
     assert( npe == 4 );
     ofs <<"POINTS " << npe <<" double"<< endl;
     for ( auto i{0U}; i<npe; i++ )
       ofs << vdata.Px(vdata.Plist(eidx,i)) <<" "<< vdata.Py(vdata.Plist(eidx,i)) <<" "<< vdata.Pz(vdata.Plist(eidx,i)) << endl;
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< 5 << endl;

     ofs << 4 <<" 0 1 2 3" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 10 << endl; // VTK_TETRA
     ofs << endl;

     // 5. writing scalar POINT_DATA (node numbers)
     // ---------------------------------------------------
     ofs <<"POINT_DATA "<< npe << endl;
     // ofs.setf( ios::scientific );
     ofs <<"SCALARS "<< var_name <<" double"<< endl;
     ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
     for ( auto i{0U}; i<npe; i++ ) ofs << vdata.Plist(eidx,i) <<" ";
     ofs << endl;
     ofs.close();
     cout <<"\nelementToVTK: file '"<< file_name <<"' written successfully."<< endl;

 } // end elementToVTK




 
} // end namespace csmp
 
 
