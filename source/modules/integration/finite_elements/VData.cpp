#include "VData.h"
#include "Box.h"
#include "FiniteElement.h"
#include "binaryReadWrite.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "TextFileInterface.h"

using namespace std;

namespace csmp {

/// default constructor: not hybrid, no nodes, nor elements
VData::VData()
 : hybrid_mesh_(false),
   first_face_(0U),
   first_interface_(0U)
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
VData::VData( const deque<size_t>& npes, 
              const deque<size_t>& epes,
              size_t nodes )
  : px(nodes),
    py(nodes),
    pz(nodes),
    pelmt(epes.size()),
    plist(epes.size()),
    pfverts(epes.size()),
    hybrid_mesh_(true),
    first_face_(npes.size()),
    first_interface_(npes.size())
 {
    for ( size_t i=0U; i<epes.size(); i++ )
      plist[i]   = vector<size_t>(npes[i]);
         
    for ( size_t i=0U; i<epes.size(); i++ )
      pfverts[i] = vector<long64>(epes[i]);
 }



/** 
    constructor for meshes with a single type of element

   @note it is assumed that there are no faces nor interfaces
*/
VData::VData( size_t nodes_per_element, size_t nbors_per_element, size_t nodes, size_t elmts )
  : px(nodes), py(nodes), pz(nodes), pelmt(elmts),
    hybrid_mesh_(false),
    first_face_(elmts),
    first_interface_(elmts)
{
    for ( size_t i=0U; i<elmts; i++ )
      plist.push_back( vector<size_t>(nodes_per_element) );
    for ( size_t k=0U; k<elmts; k++ )
      pfverts.push_back( vector<long64>(nbors_per_element) );
}




VData::~VData()
{
}


VData::VData( const VData& vd )
 : px(vd.px), py(vd.py), pz(vd.pz), 
   plist(vd.plist), pfverts(vd.pfverts),
   pelmt(vd.pelmt),
   bflags(vd.bflags),
   hybrid_mesh_(vd.hybrid_mesh_),
   first_face_( vd.first_face_ ),
   first_interface_(vd.first_interface_ )
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
  hybrid_mesh_     = a.hybrid_mesh_;
  first_face_      = a.first_face_;
  first_interface_ = a.first_interface_;

  return *this;
}



void  VData::Px( size_t i, double64 val ) 
{ assert( i<px.size() ); px[i] = val; }

void  VData::Py( size_t i, double64 val ) 
{ assert( i<py.size() ); py[i] = val; }

void  VData::Pz( size_t i, double64 val ) 
{ assert( i<pz.size() ); pz[i] = val; }

double64 VData::Px( size_t i ) const 
{ assert( i<px.size() ); return px[i]; }

double64 VData::Py( size_t i ) const 
{ assert( i<py.size() ); return py[i]; }

double64 VData::Pz( size_t i ) const 
{ assert( i<pz.size() ); return pz[i]; }

void  VData::P( size_t coordinate_axis, size_t i, double64 val )
  {
    if( coordinate_axis == 0U )
    {
        assert( i<px.size() );
        px[i] = val;
    }
    else if( coordinate_axis == 1U )
    {
        assert( i<py.size() );
        py[i] = val;
    }
    else if( coordinate_axis == 2U )
    {
        assert( i<pz.size() );
        pz[i] = val;
    }
    else
        throw std::overflow_error( "VData::P: coordinate axis is out of range");
  }

double64 VData::P( size_t coordinate_axis, size_t i ) const
{
    if( coordinate_axis == 0U )
    {
        assert( i<px.size() );
        return px[i];
    }
    else if( coordinate_axis == 1U )
    {
        assert( i<py.size() );
        return py[i];
    }
    else if( coordinate_axis == 2U )
    {
        assert( i<pz.size() );
        return pz[i];
    }
    else
        throw std::overflow_error( "VData::P: coordinate axis is out of range");
}


size_t VData::Vertices() const { return px.size(); }
    
size_t VData::Elements() const { return plist.size() - (plist.size() - first_face_); }

size_t VData::Faces() const { return first_interface_ - first_face_; }

size_t VData::InterFaces() const { return plist.size() - first_interface_; }

/// the plist contains all: elements, faces and interfaces
size_t VData::Simplices() const { return plist.size(); }

size_t VData::ElementTypes() const { return pelmt.size(); }

size_t VData::ElementNeighbors() const { return pfverts.size(); }

size_t VData::BFlags() const { return bflags.size(); }

bool VData::HybridElementTypeMesh() const { return hybrid_mesh_; }

void VData::HybridElementTypeMesh( bool hybrid_mesh ) { hybrid_mesh_ = hybrid_mesh; }

void VData::AddElementTypes( std::vector<int32>::const_iterator first,
                                    std::vector<int32>::const_iterator last )
 { pelmt.assign( first, last ); }


//  aelement ID's 0...n-1              
std::vector<size_t>::iterator  VData::PlistBegin( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistBegin: input Element ID out of range");
#endif
     return plist[eidx].begin();
  } 
  
std::vector<size_t>::iterator  VData::PlistEnd( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistEnd: input Element ID out of range");
#endif
     return plist[eidx].end();
  }  
     
std::vector<long64>::iterator  VData::PfvertsBegin( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PfvertsBegin: input Element ID out of range");
#endif
     return pfverts[eidx].begin();
  } 
  
std::vector<long64>::iterator  VData::PfvertsEnd( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PfvertsEnd: input Element ID out of range");
#endif
     return pfverts[eidx].end();
  }

// constant versions
std::vector<size_t>::const_iterator  VData::PlistBegin( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistBegin: input Element ID out of range");
#endif
     return plist[eidx].begin();
  } 
  
std::vector<size_t>::const_iterator  VData::PlistEnd( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistEnd: input Element ID out of range");
#endif
     return plist[eidx].end();
  }  
     
std::vector<long64>::const_iterator  VData::PfvertsBegin( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > pfverts.size() )
       throw std::overflow_error( "VData::PfvertsBegin: input Element ID out of range");
#endif
     return pfverts[eidx].begin();
  } 
  
std::vector<long64>::const_iterator  VData::PfvertsEnd( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > pfverts.size() )
       throw std::overflow_error( "VData::PfvertsEnd: input Element ID out of range");
#endif
     return pfverts[eidx].end();
  }


void   VData::ElementType( size_t eidx, int32 etype )
 {
    assert( eidx < pelmt.size() );
    pelmt[ eidx ] = etype;
 }

int32  VData::ElementType( size_t eidx ) const
 {
    assert( eidx < pelmt.size() );
    return pelmt[ eidx ];
 }

// TODO: used to grow VSet cell-by-cell in EclipseInterface, refactor the latter and remove this method
void VData::ResizeElementTypes( size_t new_size )
 {
    pelmt.resize( new_size, UNKNOWN );
 }



// Map Accessors (iterators)

std::deque<std::vector<size_t> >::iterator VData::PlistBegin()
 { return plist.begin(); }

std::deque<std::vector<size_t> >::iterator VData::PlistEnd()
 { return plist.end(); }

std::deque<std::vector<long64> >::iterator VData::PfvertsBegin()
 { return pfverts.begin(); }

std::deque<std::vector<long64> >::iterator VData::PfvertsEnd()
 { return pfverts.end(); }

std::map<size_t,long64>::iterator VData::BFlagsBegin()
 { return bflags.begin(); }

std::map<size_t,long64>::iterator VData::BFlagsEnd()
 { return bflags.end(); }


// uses actual element IDs (1...n)
size_t  VData::PlistSize( size_t eidx ) const
 { return plist[eidx].size(); }


// const iterators
std::vector<int32>::const_iterator  VData::PelmtBegin() const
 { return pelmt.begin(); }
    
std::vector<int32>::const_iterator  VData::PelmtEnd() const
 { return pelmt.end(); }

std::deque<std::vector<size_t> >::const_iterator VData::PlistBegin() const
 { return plist.begin(); }

std::deque<std::vector<size_t> >::const_iterator VData::PlistEnd() const
 { return plist.end(); }

size_t  VData::PfvertsSize( size_t eidx ) const
 { return pfverts[eidx].size(); }

std::deque<std::vector<long64> >::const_iterator VData::PfvertsBegin() const
 { return pfverts.begin(); }

std::deque<std::vector<long64> >::const_iterator VData::PfvertsEnd() const
 { return pfverts.end(); }

std::map<size_t,long64>::const_iterator VData::BFlagsBegin() const
 { return bflags.begin(); }

std::map<size_t,long64>::const_iterator VData::BFlagsEnd() const
 { return bflags.end(); }

void VData::AddBFlag( size_t node_id, long64 bflag )
 { bflags[ node_id ] = bflag; }
 
 
// FACE AND INTERFACE-RELATED ITERATORS

// specific element, face and interface iterators
/// iterator to CSMP finite element type of first face stored in mesh
std::vector<int32>::const_iterator  VData::PelmtFacesBegin() const {
    return std::next( pelmt.begin(), first_face_ );
 }
/// iterator to CSMP finite element type of first interface stored in mesh
std::vector<int32>::const_iterator  VData::PelmtInterfacesEnd() const {
    return std::next( pelmt.begin(), first_interface_ );
 }

/// node iterator for first face in plist (use PlistInterFacesBegin() to find last one)
std::deque<std::vector<size_t> >::const_iterator  VData::PlistFacesBegin() const {
    return std::next( plist.begin(), first_face_ );
 }

/// node iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
std::deque<std::vector<size_t> >::const_iterator  VData::PlistInterFacesBegin() const {
    return std::next( plist.begin(), first_interface_ );
 }

/// neighbor iterator for first face in plist (use PlistInterFacesBegin() to find last one)
std::deque<std::vector<long64> >::const_iterator  VData::PfvertsFacesBegin() const {
    return std::next( pfverts.begin(), first_face_ );
 }
 
/// neighbor iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
std::deque<std::vector<long64> >::const_iterator  VData::PfvertsInterfaceBegin() const {
    return std::next( pfverts.begin(), first_interface_ );
 }
 
 

void VData::Plist( size_t eidx, size_t nidx, size_t val ) 
  { 
     assert( eidx < plist.size() ); 
     assert( nidx < plist[eidx].size() ); 
     plist[eidx][nidx] = val; 
  }

size_t VData::Plist( size_t eidx, size_t nidx ) const 
  { 
     assert( eidx < plist.size() ); 
     assert( nidx < plist[eidx].size() ); 
     return plist[eidx][nidx]; 
  }

/// nidx is neighbour index
void VData::Pfvert( size_t eidx, size_t nidx, int val ) 
  { 
     assert( eidx < pfverts.size() ); 
     assert( nidx < pfverts[eidx].size() ); 
     pfverts[eidx][nidx] = val; 
  }

/// nidx is neighbour index
long64 VData::Pfvert( size_t eidx, size_t nidx ) const
  { 
     assert( eidx < pfverts.size() ); 
     assert( nidx < pfverts[eidx].size() ); 
     return pfverts[eidx][nidx]; 
  }




/**
    Sets the element type of a mono-type element mesh
*/
void VData::SingleElementType( int32 etype )
 {
    if ( pelmt.size() > 1U )
      cout <<"\nVData::SingleElementType: changing VSet from hybrid element to single-element container."<< endl;
    pelmt.clear();
    pelmt.reserve(1U);
    pelmt.push_back( etype );
    hybrid_mesh_ = false;
    // since faces are lower dimensional, at least 2 etypes would be required
    first_face_  = plist.size();
    first_interface_ = plist.size();
 }
    
    
    

/**
    Re-initialises the element-type vector that is associated
    with the VSet.
    
    @note the element type vector contains either:
    
    1. a single element (mono-element type mesh)
 
    2. element type entries for each element of the mesh
       including potential faces and interfaces
*/
void VData::ElementTypes( const vector<int32>& elmt_types )
 {
    assert( !elmt_types.empty() );
    pelmt.clear();
    pelmt.reserve(elmt_types.size());
    pelmt.assign( elmt_types.begin(), elmt_types.end() );
    hybrid_mesh_ = true;
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
                    int32 etype,
                    size_t nodes, size_t elmts )
 {
    px.resize(nodes);  vector<double64>( px ).swap( px );
    py.resize(nodes);  vector<double64>( py ).swap( py );
    pz.resize(nodes);  vector<double64>( pz ).swap( pz );
    pelmt.clear();
    plist.clear();
    pfverts.clear();
    bflags.clear();
    for ( size_t i=0U; i<elmts; ++i ) {
         plist.push_back( vector<size_t>(nodes_per_element) );
      }
    for ( size_t k=0U; k<elmts; ++k )
      pfverts.push_back( vector<long64>(nbors_per_element) );
      
    hybrid_mesh_ = false;
    
    SingleElementType( etype );

    // single element meshes have no faces nor interfaces
    first_face_      = plist.size();
    first_interface_ = plist.size();
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
void VData::Resize( const deque<int32>& etypes,
                    const deque<size_t>& npes, 
                    const deque<size_t>& epes,
                    size_t nodes,
                    size_t faces,
                    size_t interfaces )
 {
    px.resize(nodes);  vector<double64>( px ).swap( px );
    py.resize(nodes);  vector<double64>( py ).swap( py );
    pz.resize(nodes);  vector<double64>( pz ).swap( pz );

    plist.clear();
    pfverts.clear();
    bflags.clear();
    pelmt.clear();
    
    pelmt.assign( etypes.begin(), etypes.end() );
    
    for ( size_t i=0U; i<epes.size(); i++ )
      plist.push_back( vector<size_t>(npes[i]) );

    for ( size_t i=0U; i<epes.size(); i++ )
      pfverts.push_back( vector<long64>(epes[i]) );
      
    if ( etypes.size() > 1U ) hybrid_mesh_ = true;
    else hybrid_mesh_ = false;
   
    // faces
    assert( faces < npes.size() - interfaces );
    first_face_ = plist.size() - faces - interfaces;

    // interfaces
    assert( interfaces < npes.size() - faces );
    first_interface_ = plist.size() - interfaces;
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
       for (  deque<vector<size_t> >::const_iterator
              it=plist.begin(); it!=plist.end(); it++ )
         for ( vector<size_t>::const_iterator
               i=(*it).begin(); i!=(*it).end(); i++ )
           if ( (*i) > node_max ) node_max = (*i);
        
       if ( nodes < node_max )
         csmp_error.notice( CSMP_WARNING, "VData::ResizeNodes", "'plist' contains node numbers larger than desired size" );
      }
      
    if ( px.size() > 0U )
      csmp_error.notice( CSMP_WARNING, "VData::ResizeNodes", "erasing node coordinates in 'vdata'" );
    px.resize( nodes );  vector<double64>( px ).swap( px );
    py.resize( nodes );  vector<double64>( py ).swap( py );
    pz.resize( nodes );  vector<double64>( pz ).swap( pz );
        
 } // end ResizeNodes


/**
   HOPE THIS IS NEVER USED!
*/
//void VData::ResizeElementTypes( size_t elements )
//{
//    pelmt.resize( elements );
//}



/**
     Resizes 'plist' so that it can contain a different number of
     elements.
     
     @note assumes that the number of faces and interfaces is not 
     affected by the change in the number of elements.
*/
void VData::ResizePlist( size_t elements )
 {
    const size_t old_n_elements(plist.size());
  
    vector<size_t>  empty_vec;
    const size_t old_size( plist.size() );
    assert( elements > old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i=0; i<new_elmts; i++ )
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
     Resizes 'plist' so that it can contain a different number of
     elements and neighbors per element.
     
     @attention use only for meshes with a single element type.
     
     @note assumes that the number of faces and interfaces is not 
     affected by the change in the number of elements.
*/
void VData::ResizePlist( size_t elements, size_t nperelmt )
{
    assert( !hybrid_mesh_ );
    vector<size_t>  empty_vec(nperelmt,0);
    const size_t old_size( plist.size() );
    assert( elements > old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i=0; i<new_elmts; i++ )
        plist.push_back( empty_vec );

    // mono-element meshes can have no faces nor interfaces
    first_face_      = elements;
    first_interface_ = elements;
}



/**
    For mono element-type mesh change number of nodes per element.
*/
void VData::ResizeElementNodes( size_t eid, size_t nperelmt )
{
    assert( !hybrid_mesh_ );
    plist[eid].resize( nperelmt, 0 );
}



/**
    Resizes the plist without touching the number of faces or interfaces.
*/
void VData::ResizePlist( const deque<size_t>& mixed_ele_plist )
 {
    const size_t  old_n_elements(plist.size());

    // plist: vector<vector<size_t> >
    if ( !plist.empty() ) {
         cout <<"\nVData::ResizePlist: Warning: Erasing existing 'plist'"<< endl;
         plist.erase( plist.begin(), plist.end() );
      }

    hybrid_mesh_ = false;

    vector<size_t>  empty_vec;
    size_t          n_last(mixed_ele_plist[0]);

    for ( size_t i=0U; i<mixed_ele_plist.size(); ++i )
      {
         if ( !hybrid_mesh_  and  n_last != mixed_ele_plist[i] ) hybrid_mesh_ = true;
         // extending member vectors in place (avoiding copying)
         plist.push_back( empty_vec );
         plist[i].reserve( (n_last=mixed_ele_plist[i]) );
         for ( size_t j=0; j<n_last; ++j ) plist[i].push_back( 0 );
      }

    // if there were no faces or interfaces initially
    if ( first_face_ == old_n_elements or first_interface_ == old_n_elements ) {
         if ( first_face_ == old_n_elements )
           first_face_  = plist.size();
         if ( first_interface_ == old_n_elements )
           first_interface_ = plist.size();
         return;
      }

    // the offsets to first face  and interface are adjusted
    const size_t change_in_n_elements(plist.size() - old_n_elements);
    first_face_      += change_in_n_elements;
    first_interface_ += change_in_n_elements;

 } // end ResizePlist




/**
    changes the length of the 'pfverts' record
*/
void VData::ResizePfverts( size_t elements )
{
    vector<long64>  empty_vec;
    const size_t old_size( plist.size() );
    assert( elements > old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i=0; i<new_elmts; i++ )
        pfverts.push_back( empty_vec );
}


/// change the number of neighbor elements in a single element-type mesh
void VData::ResizePfverts( size_t elements, size_t nperelmt )
{
    assert( !hybrid_mesh_ );
    vector<long64>  empty_vec(nperelmt,0);
    const size_t old_size( plist.size() );
    assert( elements > old_size );
    const size_t new_elmts( elements - old_size );
    for ( size_t i=0; i<new_elmts; i++ )
        pfverts.push_back( empty_vec );
}



/**
     change the size of the storage for the number of neighbors of target element
*/
void VData::ResizeElementNeighbors( size_t eid, size_t nperelmt )
{
    pfverts[eid].resize( nperelmt, 0 );
}




/**
   just the size is changed but no data are transferred
*/
void VData::ResizePfverts( const deque<size_t>& mixed_ele_pfverts )
 {
    vector<long64>  empty_vec;
    const long64    flag_other(0);

    // pfverts: vector<vector<size_t> >
    if ( !pfverts.empty() ) {
         cout <<"\nVData::ResizePfverts: Warning: Erasing existing 'pfverts'"<< endl;
         pfverts.erase( pfverts.begin(), pfverts.end() );
      }

    for ( size_t i=0; i<mixed_ele_pfverts.size(); ++i )
      {
         // extending member vectors in place (avoiding copying)
         pfverts.push_back( empty_vec );
         pfverts[i].reserve( mixed_ele_pfverts[i] );
         for ( size_t j=0U; j<mixed_ele_pfverts[i]; j++ ) pfverts[i].push_back( flag_other );
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
   for ( deque<vector<size_t> >::const_iterator
         pit=plist.begin(); pit!=plist.end(); pit++ )
     for ( vector<size_t>::const_iterator
           it=(*pit).begin(); it!=(*pit).end(); it++ )
       if ( *it > Vertices() )
         {
            cerr <<"\nVData::Check: Wrong node ID in plist: ";
            cerr << *it <<" versus n-nodes: "<< Vertices() << endl;
            ok = false;
         }

   // pfverts
   // -------
   for ( deque<vector<long64> >::iterator
         fpit=pfverts.begin(); fpit!=pfverts.end(); fpit++ )
     for ( vector<long64>::iterator
           fit=(*fpit).begin(); fit!=(*fpit).end(); fit++ ) {
       if ( *fit > static_cast<long64>(Elements()) )
         {
            cerr <<"\nVData::Check: Non-existant high element ID in pfverts: ";
            cerr << *fit <<" versus n-elements: "<< Elements() << endl;
            ok = false;
         }
       if ( *fit < REGION_BOUNDARY )
         {
            if ( ok ) {
                 cerr <<"\n\n\n\nVData::Check: Too low boundary flag ID in pfverts: ";
                 cerr << *fit <<" versus n-elements: "<< Elements() << endl << endl;
              }
            *fit = REGION_BOUNDARY;
            ok = false;
         }
       }  

   // boundary flags
   // --------------
   if ( bflags.empty() ) 
     cout <<"\nVData::Check: No boundary nodes could be identified..."<< endl;
   else {
	    size_t  counter(0U);
	    for ( map<size_t,long64>::iterator
	          bf=bflags.begin(); bf!=bflags.end(); bf++ )
	      if ( (*bf).second < REGION_BOUNDARY ) {
	            (*bf).second = REGION_BOUNDARY;
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
//    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    // test whether nodes are not already numbered 0..1
    const deque<vector<size_t> >::const_iterator constElementsEnd(plist.end());
    for ( deque<vector<size_t> >::const_iterator
          it=plist.begin(); it!=constElementsEnd; it++ )
      for ( vector<size_t>::const_iterator n=(*it).begin(); n!=(*it).end(); n++ )
        if ( *n == 0U ) {
             cerr << "\nVData::EstablishZeroBasedNumbering: Numbering already is 0..n-1 based.\n";
             return;
          }
    // to convert: 
    // plist
    const deque<vector<size_t> >::iterator elementsEnd(plist.end());
    for ( deque<vector<size_t> >::iterator
          it=plist.begin(); it!=elementsEnd; it++ )
      for ( vector<size_t>::iterator n=(*it).begin(); n!=(*it).end(); n++ )
        *n -= 1U;

    // pfverts
    const deque<vector<long64> >::iterator pfvertsEnd(pfverts.end());
    for ( deque<vector<long64> >::iterator
          it=pfverts.begin(); it!=pfvertsEnd; it++ )
      for ( vector<long64>::iterator n=(*it).begin(); n!=(*it).end(); n++ )
        // only the neighbor element ids, not the boundary flags must be decremented
        if ( *n > 0 ) *n -= 1;

    // bconds
    map<size_t,long64>  temp;
    const map<size_t,long64>::const_iterator bflagsEnd(bflags.end());
    for ( map<size_t,long64>::const_iterator
          it=bflags.begin(); it!=bflagsEnd; it++ )
      // TODO: check whether the -1 is still correct
      temp.insert( make_pair( (*it).first-1U, (*it).second ) );
    
    bflags = temp;
    
 } // end EstablishZeroBasedNumbering







/** 
    Writes VData into a binary output file.
*/ 
void VData::OutBinary( FILE* fp ) const
 {
    size_t n0(0), n1(1), records;
    double64*  ptr(0);
    
    // 1. writing whether we are dealing with a mixed mesh
    // ---------------------------------------------------
    if ( hybrid_mesh_ ) fwrite( (void*) &n1, sizeof(size_t), 1, fp );
    else                fwrite( (void*) &n0, sizeof(size_t), 1, fp );
    
    // 2. writing all the p,c arrays or length identifiers = 0
    // -------------------------------------------------------
    // px
    if ( (records=px.size()) > 0 && (ptr=const_cast<double64*>( &(*px.begin()) )) != NULL ) 
      {
         fwrite( (void*) &records, sizeof(size_t), 1, fp );
         fwrite( (void*) ptr, sizeof(double64), records, fp );
      }
    else fwrite( (void*) &n0, sizeof(size_t), 1, fp );
    // py
    if ( (records=py.size()) > 0 && (ptr=const_cast<double64*>( &(*py.begin()) )) != NULL ) 
      {
         fwrite( (void*) &records, sizeof(size_t), 1, fp );
         fwrite( (void*) ptr, sizeof(double64), records, fp );
      }
    else fwrite( (void*) &n0, sizeof(size_t), 1, fp );
    // pz
    if ( (records=pz.size()) > 0 && (ptr=const_cast<double64*>( &(*pz.begin()) )) != NULL ) 
      {
         fwrite( (void*) &records, sizeof(size_t), 1, fp );
         fwrite( (void*) ptr, sizeof(double64), records, fp );
      }
    else fwrite( (void*) &n0, sizeof(size_t), 1, fp );

    // 3. writing pelmt, plist, pfverts, bflags
    // ----------------------------------------
    skm_C_fwrite( fp, pelmt );
    skm_C_fwrite( fp, plist );
    skm_C_fwrite( fp, pfverts );
    skm_C_fwrite( fp, bflags );
   
    // 4. offsets for faces and interfaces
    // -----------------------------------
    fwrite( (void*) &first_face_, sizeof(size_t), 1, fp );
    fwrite( (void*) &first_interface_, sizeof(size_t), 1, fp );
   
    cout <<"\nVData::OutBinary: Mesh has been successfully written to file."<< endl;

 } // end OutBinary
 
 
 
 
 
/** 
     Initialises VData from binary input file.
*/
void VData::InBinary( FILE* fp )
 {
    size_t mixed(0U), records(0U);
    
    // 1. reading whether we are dealing with a mixed mesh
    // ---------------------------------------------------
    fread( (void*) &mixed, sizeof(size_t), 1, fp );
    if ( mixed ) hybrid_mesh_ = true;
    else         hybrid_mesh_ = false;
    
    // 2. reading all the p,c arrays or length identifiers = 0
    // -------------------------------------------------------
    // px
    fread( (void*) &records, sizeof(size_t), 1, fp );
    if ( records > 0U ) {
         px.resize( records );
         vector<double64>( px ).swap( px );
         fread( (void*) &(*px.begin()), sizeof(double64), records, fp );
      }

    // py
    fread( (void*) &records, sizeof(size_t), 1, fp );
    if ( records > 0U ) {
         py.resize( records );
         vector<double64>( py ).swap( py );
         fread( (void*) &(*py.begin()), sizeof(double64), records, fp );
      }

    // pz
    fread( (void*) &records, sizeof(size_t), 1, fp );
    if ( records > 0U ) {
         pz.resize( records );
         vector<double64>( pz ).swap( pz );
         fread( (void*) &(*pz.begin()), sizeof(double64), records, fp );
      }

    // 3. reading pelmt, plist, pfverts, bflags
    // ----------------------------------------
    skm_C_fread( fp, pelmt );
    skm_C_fread( fp, plist );
    skm_C_fread( fp, pfverts );
    skm_C_fread( fp, bflags );
    
    // 4. offsets for faces and interfaces
    // -----------------------------------
    fread( (void*) &first_face_, sizeof(size_t), 1, fp );
    fread( (void*) &first_interface_, sizeof(size_t), 1, fp );
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
     for ( size_t i=0U; i<px.size(); i++ )
       ofs << i <<": \t"<< px[i] <<"\t"<< py[i] <<"\t"<< pz[i] << endl;

     // pelmt
     // -----
     ofs <<"\n'pelmt' finite element types:";
     size_t i(0U);
     for ( vector<int32>::const_iterator
           eit=pelmt.begin(); eit!=pelmt.end(); eit++, i++ )
       ofs <<"\n"<< i <<" = "<< *eit <<" = CSMP type: "<< parseFiniteElementType( static_cast<CSMP_FEM_TYPE>(*eit) );
     ofs << endl;

     // plist
     // -----
     ofs <<"\n'plist' entries of "<< plist.size() <<" elements:"<< endl;
     i = 0U;
     for ( deque<vector<size_t> >::const_iterator
           pt=plist.begin(); pt!=plist.end(); pt++, i++ )
       {
          ofs << i <<": \t";
          for ( size_t j=0U; j<(*pt).size(); j++ ) ofs << (*pt)[j] <<"\t";
          ofs << endl;
       }

     // pfverts
     // -------
     ofs <<"\n'pfverts':"<< endl;
     i = 0U;
     for ( deque<vector<long64> >::const_iterator
           ft=pfverts.begin(); ft!=pfverts.end(); ft++, i++ )
       {
          ofs << i <<": \t";
          for ( size_t j=0U; j<(*ft).size(); j++ ) ofs << (*ft)[j] <<"\t";
          ofs << endl;
       }

     // bflags
     // ------
     if ( !bflags.empty() ) ofs <<"\nBoundary flags 'bflags':"<< endl;
     for ( map<size_t,long64>::const_iterator
           bf=bflags.begin(); bf!=bflags.end(); bf++ )
       ofs << (*bf).first <<": \t"<< (*bf).second << endl;
       
     cout <<"\nVData::OutASCII: ascii file '"<< file_name <<"' written successfully."<< endl;

  } // end Out





/**
    Prints contents of VSet to a stream;
*/
void VData::Out(std::ostream& os) const
  {
     os <<"\nVData::Out: "<< endl;

     // mixed mesh
     // ----------
     if ( hybrid_mesh_ ) os <<"mesh contains different type of elements..."<< endl;
     else os <<"mesh contains only one type of element..."<< endl;

     // faces or interfaces
     // -------------------
     if ( first_face_ < plist.size() and first_face_ < first_interface_  )
       os <<"\tmesh contains "<< first_interface_ - first_face_ <<" descriptors of Face objects."<< endl;
     if ( first_interface_ < plist.size() )
       os <<"\tmesh contains "<< plist.size() - first_interface_ <<" descriptors of InterFace objects."<< endl;

     // px, py, pz
     // ----------
     assert( px.size() == py.size() );
     assert( py.size() == pz.size() );
     os<<"\n'px, py, pz' coordinates of "<< px.size() <<" nodes:"<< endl;
     for ( size_t i=0U; i<px.size(); i++ )
       os << i <<": \t"<< px[i] <<"\t"<< py[i] <<"\t"<< pz[i] << endl;

     // pelmt
     // -----
     os <<"\n'pelmt' finite element types:";
     size_t i(0U);
     for ( vector<int32>::const_iterator
           eit=pelmt.begin(); eit!=pelmt.end(); eit++, i++ )
       os <<"\n"<< i <<" = "<< *eit <<" = CSMP type: "<< parseFiniteElementType( static_cast<CSMP_FEM_TYPE>(*eit) );
     os << endl;

     // plist
     // -----
     os <<"\n'plist' entries of "<< plist.size() <<" elements:"<< endl;
     i = 0U;
     for ( deque<vector<size_t> >::const_iterator
           pt=plist.begin(); pt!=plist.end(); pt++, i++ )
       {
          os << i <<": \t";
          for ( size_t j=0U; j<(*pt).size(); j++ ) os << (*pt)[j] <<"\t";
          os << endl;
       }

     // pfverts
     // -------
     os <<"\n'pfverts':"<< endl;
     i = 0U;
     for ( deque<vector<long64> >::const_iterator
           ft=pfverts.begin(); ft!=pfverts.end(); ft++, i++ )
       {
          os << i <<": \t";
          for ( size_t j=0U; j<(*ft).size(); j++ ) os << (*ft)[j] <<"\t";
          os << endl;
       }

     // bflags
     // ------
     if ( !bflags.empty() ) os <<"\nBoundary flags 'bflags':"<< endl;
     for ( map<size_t,long64>::const_iterator
           bf=bflags.begin(); bf!=bflags.end(); bf++ )
       os << (*bf).first <<": \t"<< (*bf).second << endl;

  } // end Out()





/**
    permits initialisation of VSet from file stream.
    
    @note no faces or interfaces are considered
*/
void VData::InText( std::ifstream& ifs )
 {
     assert( ifs.is_open() );
     
     // erasing any previous records
     Erase();
     
     // 1. read node locations: 1 headlines gives number of nodes and dimensions of model
     // ---------------------------------------------------------------------------------
     char text_line[256];
     const char* const delims =" ,\t,:,\n,\r";
     int  nnodes, dim;
 
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     nnodes = atoi(strtok( text_line, delims ));
     dim    = atoi(strtok(NULL,delims));
     assert( nnodes > 0 );
     assert( dim >= 1 && dim <=3 );
     px.resize( static_cast<uint32>(nnodes) );
     vector<double64>( px ).swap( px );
     if ( dim  > 1 ) {
          py.resize( static_cast<uint32>(nnodes) );
          vector<double64>( py ).swap( py );
       }
     if ( dim == 3 ) {
          pz.resize( static_cast<uint32>(nnodes) );
          vector<double64>( pz ).swap( pz );
       }

     for ( size_t i=0U; i<static_cast<size_t>(nnodes); i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          px[i] = atof(strtok( text_line, delims ));
          if ( static_cast<size_t>(dim)  > 1U ) py[i] = atof(strtok(NULL,delims));
          if ( static_cast<size_t>(dim) == 3U ) pz[i] = atof(strtok(NULL,delims));
       }
       
     // 2. pelmt and plist, mixed mesh indicated by third number  
     //    (number indicates how many element types, thus !1= mixed mesh)
     // -----------------------------------------------------------------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     int nelements = atoi(strtok( text_line, delims ));
     int netypes   = atoi(strtok(NULL,delims));
     assert( nelements > 0 );
     assert( netypes < 20 );
     plist.resize(nelements);
     deque<vector<size_t> >( plist ).swap( plist ); 
     if ( netypes == 1 ) hybrid_mesh_ = false;
     else                hybrid_mesh_ = true;
     if ( hybrid_mesh_ ) pelmt.reserve(nelements);
     
     for ( size_t i=0U; i<static_cast<size_t>(nelements); i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          int etype = atoi(strtok( text_line, delims ));
          int npe   = atoi(strtok(NULL,delims));
          pelmt.push_back( etype );
          plist[i].resize( static_cast<uint32>(npe) );
          vector<size_t>( plist[i] ).swap( plist[i] );
          for ( size_t j=0U; j<static_cast<size_t>(npe); j++ )
            plist[i][j] = static_cast<uint32>(atoi(strtok(NULL,delims)));
       }
     
     
     // 3. pfverts
     // ----------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     int npe = atoi(strtok( text_line, delims ));
     assert( nelements == npe );
     pfverts.resize( static_cast<uint32>(nelements) );
     deque<vector<long64> >( pfverts ).swap( pfverts );
     for ( size_t i=0U; i<static_cast<size_t>(nelements); i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          npe = atoi(strtok( text_line, delims )); // neigbor elements per element (=faces)
          pfverts[i].resize( static_cast<uint32>(npe) );
          vector<long64>( pfverts[i] ).swap( pfverts[i] );
          for ( size_t j=0; j<static_cast<size_t>(npe); j++ )
            pfverts[i][j] = atoi(strtok(NULL,delims));
       }
     
     
     // 4. boundary flags
     // -----------------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     size_t nbnodes = static_cast<uint32>(atoi(strtok( text_line, delims ))); // neigbor elements per element (=faces)
     assert( static_cast<int>(nbnodes) <= nnodes );
     
     for ( size_t i=0U; i<nbnodes; i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          size_t nid   = static_cast<uint32>(atoi(strtok( text_line, delims ))); 
          long64  bflag = atoi(strtok(NULL,delims));
 //         atof(strtok(NULL,delims)); // value is ignored
          bflags.insert( make_pair( nid, bflag ) );
       }
   
    // 5. where the first face or interface - if any start in the records
    // ------------------------------------------------------------------
    do ifs.getline( text_line, 256 );
    while ( (isCommentLine(text_line) && !ifs.eof()) );
    size_t face  = static_cast<uint32>(atoi(strtok( text_line, delims )));
    assert( face <= nelements );
    first_face_ = face;
    face  = static_cast<uint32>(atoi(strtok( text_line, delims )));
    assert( face <= nelements );
    first_interface_ = face;
   
 } // end InText





 
 
/**
    Clear() for VSet
*/
void VData::Erase()
 {
    hybrid_mesh_ = false;
    px.resize(0); vector<double64>().swap( px );
    py.resize(0); vector<double64>().swap( py );
    pz.resize(0); vector<double64>().swap( pz );
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
void VData::CoordinateRange( char coordinate_axis, double64& cmin, double64& cmax ) const
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
        std::cout <<"\nVData::CoordinateRange: Coordinate axis specifier not recognized: ";
        std::cout << coordinate_axis << endl;
        std::cout <<"\n\tNot returning any values."<< endl;
     }
 } 
 
 


/**
    Scales the chosen coordinate range to the supplied minimum - maximum range.
*/
void VData::ScaleCoordinateToRange( char coordinate_axis, double64 cmin, double64 cmax )
 {
    double64 oldmin(-1.0e30), oldmax(1.0e30);
    CoordinateRange( coordinate_axis, oldmin, oldmax );
    
    // test whether we are already O.K.
    if ( cmax == oldmax && cmin == oldmin ) return;
  
    double64  old_range = oldmax - oldmin;
    double64  new_range = cmax   - cmin;  

    if ( coordinate_axis == 'x' || coordinate_axis == 'X' )
      for ( size_t i=0; i<px.size(); i++ )
        px[i] = cmin + ((px[i] - oldmin)/old_range) * new_range;

    else if ( coordinate_axis == 'y' || coordinate_axis == 'Y' )
      for ( size_t i=0; i<py.size(); i++ )
        py[i] = cmin + ((py[i] - oldmin)/old_range) * new_range;

    else if ( coordinate_axis == 'z' || coordinate_axis == 'Z' )
      for ( size_t i=0; i<pz.size(); i++ )
        pz[i] = cmin + ((pz[i] - oldmin)/old_range) * new_range;
 } 
 
 

/**
    comparitor for VSets.
*/
bool  VData::operator==( const VData& vd ) const
 {
    bool return_value(true);
    if ( !(hybrid_mesh_ == vd.hybrid_mesh_) ) return false;
    if ( first_face_      != vd.first_face_ ) return false;
    if ( first_interface_ != vd.first_interface_ ) return false;
    if ( !(px == vd.px) ) { cerr<<"\nVData::operator== failed 'px' comparison."; return_value = false; }
    if ( !(py == vd.py) ) { cerr<<"\nVData::operator== failed 'py' comparison."; return_value = false; }
    if ( !(pz == vd.pz) ) { cerr<<"\nVData::operator== failed 'pz' comparison."; return_value = false; }
    if ( !(pelmt == vd.pelmt) ) { cerr<<"\nVData::operator== failed 'pelmt' comparison."; return_value = false; }
    if ( !(plist == vd.plist) ) { cerr<<"\nVData::operator== failed 'plist' comparison."; return_value = false; }
    if ( !(pfverts == vd.pfverts) ) { cerr<<"\nVData::operator== failed 'pfverts' comparison."; return_value = false; }
    if ( !(bflags == vd.bflags) ) { cerr<<"\nVData::operator== failed 'bflags' comparison."; return_value = false; }
    
    return return_value;
 }




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
    for ( deque<vector<size_t> >::const_iterator
          it=plist.begin(); it!=plist.end(); it++ )
      for ( vector<size_t>::const_iterator 
            vt=(*it).begin(); vt!=(*it).end(); vt++ )
        node_set.insert( *vt );
        
    // if nothing can be done because there are not enough node coordinates
    if ( node_set.size() > px.size() )
      throw csmp::Exception( CSMP_FATAL_ERROR, "VData::CheckForOrphanNodes",
                     "plist contains more nodes node coordinate arrays");
 
    // detecting discontinuities in the node numbering
    if ( (*max_element( node_set.begin(), node_set.end() )) >= node_set.size() )
      throw csmp::Exception( CSMP_FATAL_ERROR, "VData::CheckForOrphanNodes",
                     "there is a discontinuity in the node-numbers stored in 'plist'");
    
    // if there are extra nodes but everything else is OK, these can be removed 
    // swap trick is used to trim extra capacity from vectors                
    if ( node_set.size() < px.size() ) {
          // if the nodes are not stored in the plist they are orphan and can be removed
          if ( eliminate_orphan_nodes ) {
               px.resize( node_set.size() );
               vector<double64>( px ).swap( px );
               if ( !py.empty() ) {
                    py.resize( node_set.size() );
                    vector<double64>( py ).swap( py );
                 }
               if ( !pz.empty() ) {
                    pz.resize( node_set.size() );
                    vector<double64>( pz ).swap( pz );
                 }
               //throw csmp::Exception( CSMP_INFO, "VData::CheckForOrphanNodes",
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
void VData::ReduceTo( const map<size_t,size_t>& o_n_elmt_ids, map<size_t,size_t>& o_n_node_ids )
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
         vector<int32> new_element_types( o_n_elmt_ids.size() );
         for ( map<size_t,size_t>::const_iterator 
               it=o_n_elmt_ids.begin(); it!=o_n_elmt_ids.end(); ++it ) {
              // range-checked access of the 2 containers using at()
              new_element_types.at( (*it).second ) = pelmt.at( (*it).first ); 
           }
	       pelmt.assign( new_element_types.begin(), new_element_types.end() );
      }
    // else nothing needs to be done
   
    // checking wether the mesh is still mixed element type and resizing 
    // the element type deque if not
    set<int32>  n_etypes( pelmt.begin(), pelmt.end() );
    if ( n_etypes.size() == 1U ) {
         pelmt.resize(1U);
         vector<int32>( pelmt ).swap( pelmt );
      }


    // 1. erasing abandonned plist & pfverts entries
    // ---------------------------------------------
    // copying those entries which are to be retained
    // and creating a node correspondance map with the new node IDs
    vector<vector<size_t> >  new_plist( o_n_elmt_ids.size() );
    vector<vector<long64> >  new_pfverts( with_connectivity ? o_n_elmt_ids.size() : 0 );
    size_t                   n_node(0U);
    //  old   & new node numbers
    o_n_node_ids.clear();
    assert( plist.size() == pfverts.size() || !with_connectivity );
    for ( map<size_t,size_t>::const_iterator
          it=o_n_elmt_ids.begin(); it!=o_n_elmt_ids.end(); it++ )
      {
         const size_t  eidx( (*it).first );
         const size_t  neidx( (*it).second );
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
         for ( vector<size_t>::const_iterator
               nit=plist[eidx].begin(); nit!=plist[eidx].end(); nit++ ) {
              pair<map<size_t,size_t>::iterator,bool>
              iit = o_n_node_ids.insert( make_pair( *nit, n_node ) );
              // if a new node was added the node number is incremented
              if ( iit.second == true ) n_node++;
           }
      }
    // verifying that plist and pfverts do not contain empty elements
    for ( vector<vector<size_t> >::const_iterator
          it=new_plist.begin(); it!=new_plist.end(); it++ )
      if ( (*it).empty() )
        throw csmp::Exception( CSMP_FATAL_ERROR, "VData::ReduceTo",
                              "The reduced 'plist' deque contains empty entries. Unable to continue");
    if( with_connectivity )
        for ( vector<vector<long64> >::const_iterator
              it=new_pfverts.begin(); it!=new_pfverts.end(); it++ )
          if ( (*it).empty() )
            throw csmp::Exception( CSMP_FATAL_ERROR, "VData::ReduceTo",
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
        for ( deque<vector<long64> >::iterator it=pfverts.begin(); it!=pfverts.end(); it++ )
          for ( vector<long64>::iterator pit=(*it).begin(); pit!=(*it).end(); pit++ )
            // only if there was a neighboring element before its ID is updated
            if ( (*pit) > 0 ) {
                 map<size_t,size_t>::const_iterator
                   eit=o_n_elmt_ids.find( static_cast<size_t>(*pit) );
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
                 *pit = static_cast<long64>((*eit).second);
              }
    }

    // 3. 'plist': updating node IDs if necessary
    // ------------------------------------------
    if ( px.size() != o_n_node_ids.size() ) {
         for ( deque<vector<size_t> >::iterator it=plist.begin(); it!=plist.end(); it++ )
           for ( vector<size_t>::iterator pit=(*it).begin(); pit!=(*it).end(); pit++ )
             // only if there was a neighboring element before its ID is updated
              {
                  map<size_t,size_t>::const_iterator nit=o_n_node_ids.find( *pit );
                  if ( nit == o_n_node_ids.end() )
                     throw csmp::Exception( CSMP_ERROR, "VData::ReduceTo",
                                                   "'plist' node ID could not be updated");
                  else *pit = (*nit).second;
               }

         // 4. Eliminating orphan nodes
         // ---------------------------
         if( csmp_error.Verbose() )
             cout <<"\nVData::ReduceTo: eliminating orphan nodes..."<< endl;
         // copying the retained nodes over into new vectors
         // px
         vector<double64>  new_px( o_n_node_ids.size() );
         //        old_ID    new_ID
         for ( map<size_t,size_t>::const_iterator
               nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
           new_px[ (*nit).second ] = px[ (*nit).first ];
         // updating px
         px = new_px;
         new_px.clear();
         // py
         if ( !py.empty() ) {
              vector<double64>  new_py( o_n_node_ids.size() );
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
              vector<double64>  new_pz( o_n_node_ids.size() );
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
         map<size_t,long64>  new_bflags;
         for ( map<size_t,long64>::const_iterator bit=bflags.begin(); bit!=bflags.end(); bit++ ) {
              map<size_t,size_t>::const_iterator nit=o_n_node_ids.find( (*bit).first );
              if ( nit!=o_n_node_ids.end() )
                //                             new node idx  old bflag
                new_bflags.insert( make_pair( (*nit).second, (*bit).second ) );
           }
         bflags = new_bflags;
         new_bflags.clear();
      }

    // 6. updating the 'mixed_mesh' boolean variable
    // ---------------------------------------------
    // (this check will not work if there are two element types with the same
    //  number of nodes)
    size_t  nodes_per_element = plist[0].size();

    hybrid_mesh_=false;
    for ( deque<vector<size_t> >::iterator it=plist.begin(); it!=plist.end(); it++ )
      if ( (*it).size() != nodes_per_element ) {
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
             cout <<"\nVData::ReduceTo: Successfully reduced number of simplices by "<< n_elements_eliminated <<" simplex from ";
             cout << old_nelements <<" to "<< plist.size() << endl;
          }
    }

 } // end ReduceTo
 
 
 
 
/**
    Considering all known element types, returns their interpolation order.
    
    @attention assumption is made that all elements in the mesh use the same
    order element interpolation functions.
*/
size_t  VData::OrderOfFiniteElementInterpolationFunctions() const
 {
     for ( vector<int32>::const_iterator 
           it=PelmtBegin(); it!=PelmtEnd(); it++ ) {
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


 
} // end namespace csmp
 
 
