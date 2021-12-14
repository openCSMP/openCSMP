#include "VData.h"
#include "Box.h"
#include "FiniteElement.h"
#include "binaryReadWrite.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "TextFileIO.h"
#include "CSMP_mathUtilities.h"
#include "CSMP_highLevelUtilities.h"
#include "CSMP_ElementSpecifications.h"

using namespace std;

namespace csmp {

/// default constructor: not hybrid, no nodes, nor elements
VData::VData()
 : hybrid_mesh_(false),
   bflags(),
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
VData::VData( const deque<size_t>& npes,
              const deque<size_t>& epes,
              size_t nodes )
  : px(nodes),
    py(nodes),
    pz(nodes),
    bflags(nodes),
    pelmt(epes.size()),
    plist(epes.size()),
    pfverts(epes.size()),
    hybrid_mesh_(true),
    first_face_(npes.size()),
    first_interface_(npes.size())
 {
    for ( size_t i=0U; i<epes.size(); i++ )
      plist[i]   = vector<int64_t>(npes[i]);
         
    for ( size_t i=0U; i<epes.size(); i++ )
      pfverts[i] = vector<int64_t>(epes[i]);
 }




/** 
    constructor for meshes with a single type of element

   @note it is assumed that there are no faces nor interfaces
*/
VData::VData( size_t nodes_per_element, size_t nbors_per_element, size_t nodes, size_t elmts )
  : px(nodes), py(nodes), pz(nodes), bflags(nodes), pelmt(elmts),
    hybrid_mesh_(false),
    first_face_(elmts),
    first_interface_(elmts)
{
    for ( size_t i=0U; i<elmts; i++ )
      plist.push_back( vector<int64_t>(nodes_per_element) );
    for ( size_t k=0U; k<elmts; k++ )
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

void  VData::P( size_t coordinate_axis, size_t i, double val )
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

double VData::P( size_t coordinate_axis, size_t i ) const
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
size_t VData::TotalNumberOfCells() const { return plist.size(); }

size_t VData::ElementTypes() const { return pelmt.size(); }

size_t VData::ElementNeighbors() const { return pfverts.size(); }

size_t VData::BFlags() const { return bflags.size(); }

bool VData::HybridElementTypeMesh() const { return hybrid_mesh_; }

void VData::HybridElementTypeMesh( bool hybrid_mesh ) { hybrid_mesh_ = hybrid_mesh; }

void VData::AddElementTypes( std::vector<int8_t>::const_iterator first,
                             std::vector<int8_t>::const_iterator last )
 { pelmt.assign( first, last ); }

  void VData::AddElementTypes( std::deque<int8_t>::const_iterator first,
                               std::deque<int8_t>::const_iterator last )
  { pelmt.assign( first, last ); }
  

    /// reports whether the model contains only isoparametric element types
bool  VData::IsoparametricElementMesh() const
 {
     for ( auto it : pelmt )
       if ( !CSMP_ElementSpecifications::IsIsoparametric(it) )
         return false;
     return true;
 }



//  aelement ID's 0...n-1              
std::vector<int64_t>::iterator  VData::PlistBegin( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistBegin: input Element ID out of range");
#endif
     return plist[eidx].begin();
  } 
  
std::vector<int64_t>::iterator  VData::PlistEnd( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistEnd: input Element ID out of range");
#endif
     return plist[eidx].end();
  }  
     
std::vector<int64_t>::iterator  VData::PfvertsBegin( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PfvertsBegin: input Element ID out of range");
#endif
     return pfverts[eidx].begin();
  } 
  
std::vector<int64_t>::iterator  VData::PfvertsEnd( size_t eidx ) {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PfvertsEnd: input Element ID out of range");
#endif
     return pfverts[eidx].end();
  }

// constant versions
std::vector<int64_t>::const_iterator  VData::PlistBegin( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistBegin: input Element ID out of range");
#endif
     return plist[eidx].begin();
  } 
  
std::vector<int64_t>::const_iterator  VData::PlistEnd( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > plist.size() )
       throw std::overflow_error( "VData::PlistEnd: input Element ID out of range");
#endif
     return plist[eidx].end();
  }  
     
std::vector<int64_t>::const_iterator  VData::PfvertsBegin( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > pfverts.size() )
       throw std::overflow_error( "VData::PfvertsBegin: input Element ID out of range");
#endif
     return pfverts[eidx].begin();
  } 
  
std::vector<int64_t>::const_iterator  VData::PfvertsEnd( size_t eidx ) const {
#ifndef NDEBUG
     if ( eidx > pfverts.size() )
       throw std::overflow_error( "VData::PfvertsEnd: input Element ID out of range");
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

// TODO: used to grow VSet cell-by-cell in EclipseInterface, refactor the latter and remove this method
/*
void VData::ResizeElementTypes( size_t new_size )
 {
    pelmt.resize( new_size, UNKNOWN );
 }
*/


// Accessors (iterators)

std::deque<std::vector<int64_t> >::iterator VData::PlistBegin()
 { return plist.begin(); }

std::deque<std::vector<int64_t> >::iterator VData::PlistEnd()
 { return plist.end(); }

std::deque<std::vector<int64_t> >::iterator VData::PfvertsBegin()
 { return pfverts.begin(); }

std::deque<std::vector<int64_t> >::iterator VData::PfvertsEnd()
 { return pfverts.end(); }

std::vector<std::int8_t>::iterator VData::BFlagsBegin()
 { return bflags.begin(); }

std::vector<std::int8_t>::iterator VData::BFlagsEnd()
 { return bflags.end(); }


// uses actual element IDs (1...n)
size_t  VData::PlistSize( size_t eidx ) const
 { return plist[eidx].size(); }


// const iterators
std::vector<int8_t>::const_iterator  VData::PelmtBegin() const
 { return pelmt.begin(); }
    
std::vector<int8_t>::const_iterator  VData::PelmtEnd() const
 { return pelmt.end(); }

std::deque<std::vector<int64_t> >::const_iterator VData::PlistBegin() const
 { return plist.begin(); }

std::deque<std::vector<int64_t> >::const_iterator VData::PlistEnd() const
 { return plist.end(); }

size_t  VData::PfvertsSize( size_t eidx ) const
 { return pfverts[eidx].size(); }

std::deque<std::vector<int64_t> >::const_iterator VData::PfvertsBegin() const
 { return pfverts.begin(); }

std::deque<std::vector<int64_t> >::const_iterator VData::PfvertsEnd() const
 { return pfverts.end(); }

std::vector<std::int8_t>::const_iterator VData::BFlagsBegin() const
 { return bflags.begin(); }

std::vector<std::int8_t>::const_iterator VData::BFlagsEnd() const
 { return bflags.end(); }


void VData::AddBFlag( size_t node_id, std::int8_t bflag )
 {
    assert( node_id < bflags.size() );
    BOX_BOUNDARY flag = static_cast<BOX_BOUNDARY>(bflag);
    // if the boundary flag integer value is outside of the range of defined values
    if ( flag < MULTIPLE_BOUNDARIES - 1 ) {
         cerr <<"\nVData::AddBFlag: boundary flag "<< parseBoundary( flag );
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

 
 
/**
      Finds boundary identifier if any.
*/
std::int8_t VData::BoundaryFlag( size_t vertex ) const
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
std::vector<int8_t>::const_iterator  VData::PelmtFacesBegin() const {
    return std::next( pelmt.begin(), first_face_ );
 }
 
/// iterator to CSMP finite element type of first interface stored in mesh
std::vector<int8_t>::const_iterator  VData::PelmtInterfacesBegin() const {
    return std::next( pelmt.begin(), first_interface_ );
 }

std::deque<std::vector<int64_t> >::const_iterator  VData::PlistElmtsBegin() const {
    return plist.begin();
 }

std::deque<std::vector<int64_t> >::const_iterator  VData::PlistElmtsEnd() const {
    return std::next( plist.begin(), first_face_ );
 }

std::deque<std::vector<int64_t> >::const_iterator  VData::PlistFacesBegin() const {
    return std::next( plist.begin(), first_face_ );
 }

std::deque<std::vector<int64_t> >::const_iterator  VData::PlistFacesEnd() const {
    return std::next( plist.begin(), first_interface_ );
 }

std::deque<std::vector<int64_t> >::const_iterator  VData::PlistInterFacesBegin() const {
    return std::next( plist.begin(), first_interface_ );
 }

std::deque<std::vector<int64_t> >::const_iterator  VData::PlistInterFacesEnd() const {
    return plist.end();
 }

/// neighbor iterator for first face in plist (use PlistInterFacesBegin() to find last one)
std::deque<std::vector<int64_t> >::const_iterator  VData::PfvertsFacesBegin() const {
    return std::next( pfverts.begin(), first_face_ );
 }
 
/// neighbor iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
std::deque<std::vector<int64_t> >::const_iterator  VData::PfvertsInterfaceBegin() const {
    return std::next( pfverts.begin(), first_interface_ );
 }
 
 

void VData::Plist( size_t eidx, size_t nidx, size_t val ) 
  { 
     assert( eidx < plist.size() ); 
     assert( nidx < plist[eidx].size() ); 
     plist[eidx][nidx] = val; 
  }

int64_t  VData::Plist( size_t eidx, size_t nidx ) const
  { 
     assert( eidx < plist.size() );
    if (nidx >= plist[eidx].size()) {
      std::cerr << "plist[" << eidx << ".size() = " <<plist[eidx].size() << '\n';
    }
     assert( nidx < plist[eidx].size() ); 
     return plist[eidx][nidx]; 
  }

/// nidx is neighbour index
void VData::Pfvert( size_t eidx, size_t nidx, int64_t  val )
  { 
     assert( eidx < pfverts.size() ); 
     assert( nidx < pfverts[eidx].size() ); 
     pfverts[eidx][nidx] = val; 
  }

/// nidx is neighbour index
int64_t  VData::Pfvert( size_t eidx, size_t nidx ) const
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
    
    pelmt.clear();
    plist.clear();
    pfverts.clear();

    for ( size_t i=0U; i<elmts; ++i ) {
         plist.push_back( vector<int64_t>(nodes_per_element) );
      }
    for ( size_t k=0U; k<elmts; ++k )
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
                      const deque<size_t>& npes,
                      const deque<size_t>& epes,
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
    
    plist.clear();
    pfverts.clear();
    pelmt.clear();
    
    pelmt.assign( etypes.begin(), etypes.end() );
    
    for ( size_t i=0U; i<nrCells; i++ )
      plist.push_back( vector<int64_t>(npes[i]) );
    
    for ( size_t i=0U; i<nrCells; i++ )
      pfverts.push_back( vector<int64_t>(epes[i]) );
    
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
         csmp_error.notice( WARNING, "VData::ResizeNodes", "'plist' contains node numbers larger than desired size" );
      }
      
    if ( px.size() > 0U )
      csmp_error.notice( WARNING, "VData::ResizeNodes", "erasing node coordinates in 'vdata'" );
    px.resize( nodes );  vector<double>( px ).swap( px );
    py.resize( nodes );  vector<double>( py ).swap( py );
    pz.resize( nodes );  vector<double>( pz ).swap( pz );
    ResizeBFlags();
        
 } // end ResizeNodes




/**
      Updates length of vector to that of the node vector, setting potential new flags to zero.
*/
void VData::ResizeBFlags()
 {
    bflags.resize( px.size() );
    vector<std::int8_t>( bflags ).swap( bflags );
    fill( bflags.begin(), bflags.end(), NOT );
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
void VData::ResizePlist( size_t elements, size_t nperelmt )
{
    assert( Faces() == 0 );
    assert( InterFaces() == 0 );
    vector<int64_t>  empty_vec(nperelmt,0);
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

    vector<int64_t>  empty_vec;
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
    for ( size_t i=0; i<new_elmts; i++ )
        pfverts.push_back( empty_vec );
}


/**
    increases the number of neighbor elements in a single element-type mesh
    
        @attention this works only of the size of the pfverts record is increased
 */
void VData::ResizePfverts( size_t elements, size_t nperelmt )
{
    if ( first_face_  != plist.size() || first_interface_ != plist.size() )
      throw logic_error( "VData::ResizePfverts(size_t,size_t): method does not handle VSets with Face and InterFace objects yet");

    vector<int64_t>  empty_vec(nperelmt,0);
    const size_t old_size( pfverts.size() );
    assert( elements >= old_size );
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
    vector<int64_t>  empty_vec;
    const int64_t    flag_other(0);

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
   for ( deque<vector<int64_t> >::const_iterator
         pit=plist.begin(); pit!=plist.end(); pit++ )
     for ( vector<int64_t>::const_iterator
           it=(*pit).begin(); it!=(*pit).end(); it++ )
       if ( *it > Vertices() )
         {
            cerr <<"\nVData::Check: Wrong node ID in plist: ";
            cerr << *it <<" versus n-nodes: "<< Vertices() << endl;
            ok = false;
         }

   // pfverts
   // -------
   for ( deque<vector<int64_t> >::iterator
         fpit=pfverts.begin(); fpit!=pfverts.end(); fpit++ )
     for ( vector<int64_t>::iterator
           fit=(*fpit).begin(); fit!=(*fpit).end(); fit++ ) {
       if ( *fit > static_cast<int64_t>(Elements()) )
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
    const deque<vector<int64_t> >::const_iterator constElementsEnd(plist.end());
    for ( deque<vector<int64_t> >::const_iterator
          it=plist.begin(); it!=constElementsEnd; it++ )
      for ( vector<int64_t>::const_iterator n=(*it).begin(); n!=(*it).end(); n++ )
        if ( *n == 0U ) {
             cerr << "\nVData::EstablishZeroBasedNumbering: Numbering already is 0..n-1 based.\n";
             return;
          }
    // to convert: 
    // plist
    const deque<vector<int64_t> >::iterator elementsEnd(plist.end());
    for ( deque<vector<int64_t> >::iterator
          it=plist.begin(); it!=elementsEnd; it++ )
      for ( vector<int64_t>::iterator n=(*it).begin(); n!=(*it).end(); n++ )
        *n -= 1U;

    // pfverts
    const deque<vector<int64_t> >::iterator pfvertsEnd(pfverts.end());
    for ( deque<vector<int64_t> >::iterator
          it=pfverts.begin(); it!=pfvertsEnd; it++ )
      for ( vector<int64_t>::iterator n=(*it).begin(); n!=(*it).end(); n++ )
        // only the neighbor element ids, not the boundary flags must be decremented
        if ( *n > 0 ) *n -= 1;

    // bflags & bconds - nothing needs to be done
    
 } // end EstablishZeroBasedNumbering







/** 
    Writes VData into a binary output file.
*/ 
void VData::OutBinary( fstream& fp ) const
 {
    static_assert( sizeof(size_t) == sizeof(streamsize), "VData::OutBinary: on this system 'std::streamsize' is not equal to size_t" );
    size_t n0(0), n1(1), records;
    double*  ptr(0);
    
    // 1. writing whether we are dealing with a mixed mesh
    // ---------------------------------------------------
    if ( hybrid_mesh_ ) fp.write( reinterpret_cast<const char*>(&n1), sizeof(size_t) );
    else                fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
    
    // 2. writing all the p,c arrays or length identifiers = 0
    // -------------------------------------------------------
   {
     BinaryFileSectionWrite sect(fp, "VSETCORD");
     // px
    if ( (records=px.size()) > 0 && (ptr=const_cast<double*>( &(*px.begin()) )) != NULL ) 
      {
         fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
         fp.write( reinterpret_cast<const char*>(ptr), sizeof(double) * records );
      }
    else fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
    // py
    if ( (records=py.size()) > 0 && (ptr=const_cast<double*>( &(*py.begin()) )) != NULL ) 
      {
         fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
         fp.write( reinterpret_cast<const char*>(ptr), sizeof(double) * records );
      }
    else fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
    // pz
    if ( (records=pz.size()) > 0 && (ptr=const_cast<double*>( &(*pz.begin()) )) != NULL ) 
      {
         fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );
         fp.write( reinterpret_cast<const char*>(ptr), sizeof(double) * records );
      }
    else fp.write( reinterpret_cast<const char*>(&n0), sizeof(size_t) );
   }

    // 3. writing pelmt, plist, pfverts, bflags
    // ----------------------------------------
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
   
    // 4. offsets for faces and interfaces
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
    static_assert( sizeof(size_t) == sizeof(streamsize), "VData::InBinary: on this system 'std::streamsize' is not equal to size_t" );
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
    
    // 4. offsets for faces and interfaces
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
     for ( size_t i=0U; i<px.size(); i++ )
       ofs << i <<": \t"<< px[i] <<"\t"<< py[i] <<"\t"<< pz[i] << endl;

     // pelmt
     // -----
     ofs <<"\n'pelmt' finite element types:";
     size_t i(0U);
     for ( vector<int8_t>::const_iterator
           eit=pelmt.begin(); eit!=pelmt.end(); eit++, i++ )
       ofs <<"\n"<< i <<" = "<< *eit <<" = CSMP type: "<< parseFiniteElementType( static_cast<CSMP_FEM_TYPE>(*eit) );
     ofs << endl;

     // plist
     // -----
     ofs <<"\n'plist' entries of "<< plist.size() <<" elements:"<< endl;
     i = 0U;
     for ( deque<vector<int64_t> >::const_iterator
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
     for ( deque<vector<int64_t> >::const_iterator
           ft=pfverts.begin(); ft!=pfverts.end(); ft++, i++ )
       {
          ofs << i <<": \t";
          for ( size_t j=0U; j<(*ft).size(); j++ ) ofs << (*ft)[j] <<"\t";
          ofs << endl;
       }

     // bflags
     // ------
     if ( !bflags.empty() ) ofs <<"\nBoundary flags 'bflags':"<< endl;
     size_t n_node(0U);
     for ( auto bf=bflags.begin(); bf!=bflags.end(); bf++ )
       ofs << n_node++ <<": \t"<< (*bf) << endl;
       
     cout <<"\nVData::OutASCII: ascii file '"<< file_name <<"' written successfully."<< endl;

  } // end Out





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
     if ( first_face_ < plist.size() and first_face_ < first_interface_  )
       cout <<"\tmesh contains "<< first_interface_ - first_face_ <<" descriptors of Face objects."<< endl;
     if ( first_interface_ < plist.size() )
       cout <<"\tmesh contains "<< plist.size() - first_interface_ <<" descriptors of InterFace objects."<< endl;

     // px, py, pz
     // ----------
     assert( px.size() == py.size() );
     assert( py.size() == pz.size() );
     cout<<"\n'px, py, pz' coordinates of "<< px.size() <<" nodes:"<< endl;
     for ( size_t i=0U; i<px.size(); i++ )
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
          for ( size_t j=0U; j<(*pt).size(); j++ ) cout << (*pt)[j] <<"\t";
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
          for ( size_t j=0U; j<(*ft).size(); j++ ) {
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
     char text_line[NAME_STRING];
     const char* const delims =" ,\t,:,\n,\r";
 
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     size_t nnodes = atoi(strtok( text_line, delims ));
     size_t dim    = atoi(strtok(NULL,delims));
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

     for ( size_t i=0U; i<nnodes; i++ ) {
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
     
     for ( size_t i=0U; i<nelements; i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          int64_t  etype = atoi(strtok( text_line, delims ));
          assert( etype >= 0 && etype < 50 );
          pelmt.push_back( static_cast<int8_t>(etype) );
          int64_t  npe  = atoi(strtok(NULL,delims));
          assert( npe > 2 );
          plist[i].resize( static_cast<size_t>(npe) );
          vector<int64_t>( plist[i] ).swap( plist[i] );
          for ( size_t j=0U; j<static_cast<size_t>(npe); j++ )
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
     for ( size_t i=0U; i<static_cast<size_t>(nelements); i++ ) {
          do ifs.getline( text_line, 256 );
          while ( (isCommentLine(text_line) && !ifs.eof()) );
          npe = atoi(strtok( text_line, delims )); // neigbor elements per element (=faces)
          pfverts[i].resize( static_cast<size_t>(npe) );
          vector<int64_t>( pfverts[i] ).swap( pfverts[i] );
          for ( size_t j=0; j<static_cast<size_t>(npe); j++ )
            pfverts[i][j] = atoi(strtok(NULL,delims));
       }
     
     
     // 4. boundary flags
     // -----------------
     do ifs.getline( text_line, 256 );
     while ( (isCommentLine(text_line) && !ifs.eof()) );
     size_t nbnodes = atoi(strtok( text_line, delims )); // neigbor elements per element (=faces)
     assert( static_cast<int>(nbnodes) <= nnodes );
     bflags.resize(nnodes,0);
     for ( size_t i=0U; i<nbnodes; i++ ) {
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
         cerr<<"\nVData::operator== failed 'pelmt' comparison.";
         if ( pelmt.size() != vd.pelmt.size() ) cerr <<"\nelement type records have different sizes.\n";
         for ( size_t i{0}; i<pelmt.size(); ++i )
           if ( pelmt[i] != vd.pelmt[i] )
             cerr <<"\n\t\t"<< parseFiniteElementType(pelmt[i]) <<" "<< i <<": "<< static_cast<int>(pelmt[i]) <<" vs "<< static_cast<int>(vd.pelmt[i]);
         cerr << endl << endl;
         return_value = false;
      }
    
    if ( !(plist == vd.plist) ) {
         cerr<<"\nVData::operator== failed 'plist' comparison.";
         for ( size_t i{0}; i<plist.size(); ++i )
           if ( plist[i] != vd.plist[i] ) {
                if ( plist[i].size() != vd.plist[i].size() ) cerr <<"\ncell "<< i <<" nodes-per-element records have different sizes.\n";
                else {
                     cerr <<"\n"<< parseFiniteElementType(pelmt[i]) <<" "<< i <<": 'plist' member comparison: ";
                     for ( size_t j{0}; j<plist[i].size(); ++j )
                       cerr <<"\n\t\t"<< j <<": "<< plist[i][j] <<" vs "<< vd.plist[i][j];
                  }
                cerr << endl << endl;
             }
         return_value = false;
      }
    
    if ( !(pfverts == vd.pfverts) ) {
         cerr<<"\nVData::operator== failed 'pfverts' comparison.\n";
         for ( size_t i{0}; i<pfverts.size(); ++i )
           if ( pfverts[i] != vd.pfverts[i] ) {
                if ( pfverts[i].size() != vd.pfverts[i].size() ) cerr <<"\n\tcell "<< i <<" element-neighbor records have different sizes.\n";
                else {
                     cerr <<"\n"<< parseFiniteElementType(pelmt[i]) <<" "<< i <<": 'pfvert' member comparison: ";
                     for ( size_t j{0}; j<pfverts[i].size(); ++j )
                       cerr <<"\n\t\t"<< j <<": "<< pfverts[i][j] <<" vs "<< vd.pfverts[i][j];
                  }
                cerr << endl << endl;
             }
         return_value = false;
      }
      
    if ( !(bflags == vd.bflags) ) {
         cerr<<"\nVData::operator== failed 'bflags' comparison.";
         if ( bflags.size() != vd.bflags.size() ) cerr <<"\n\tboundary flag records have different sizes.\n";
         for ( size_t i{0}; i<bflags.size(); ++i )
           if ( bflags[i] != vd.bflags[i] )
             cerr <<"\n\t\tnode flag "<< i <<": "<< static_cast<int>(bflags[i]) <<" vs "<< static_cast<int>(vd.bflags[i]);
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
         for ( vector<int64_t>::const_iterator
               nit=plist[eidx].begin(); nit!=plist[eidx].end(); nit++ ) {
              pair<map<size_t,size_t>::iterator,bool>
              iit = o_n_node_ids.insert( make_pair( *nit, n_node ) );
              // if a new node was added the node number is incremented
              if ( iit.second == true ) n_node++;
           }
      }
    // verifying that plist and pfverts do not contain empty elements
    for ( vector<vector<int64_t> >::const_iterator
          it=new_plist.begin(); it!=new_plist.end(); it++ )
      if ( (*it).empty() )
        throw csmp::Exception( FATAL_ERROR, "VData::ReduceTo",
                              "The reduced 'plist' deque contains empty entries. Unable to continue");
    if( with_connectivity )
        for ( vector<vector<int64_t> >::const_iterator
              it=new_pfverts.begin(); it!=new_pfverts.end(); it++ )
          if ( (*it).empty() )
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
         vector<std::int8_t> new_bflags(px.size(),0);
         //        old_ID    new_ID
         for ( map<size_t,size_t>::const_iterator
               nit=o_n_node_ids.begin(); nit!=o_n_node_ids.end(); nit++ )
           new_bflags[ (*nit).second ] = bflags[ (*nit).first ];
      }

    // 6. updating the 'mixed_mesh' boolean variable
    // ---------------------------------------------
    // (this check will not work if there are two element types with the same
    //  number of nodes)
    size_t  nodes_per_element = plist[0].size();

    hybrid_mesh_=false;
    for ( deque<vector<int64_t> >::iterator it=plist.begin(); it!=plist.end(); it++ )
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
size_t  VData::OrderOfFiniteElementInterpolationFunctions() const
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
   
   @attention if the line elements do not know their neighbors already, this method will not work!
   
   @attention assumes node numbering in 'plist' follows the convention corner nodes first, then Interior nodes.
   
   @attention method assumes that all manifolds (connections of more than 2 line elements at a node have been removed
   =disambiguated before. Thus it relies on a valid neighbor connectivity of the line elements, which must stored in pfverts.
   
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
      int64_t  elmt_idx(0U);

       while( eit != end ) {
            // if this is a line element at the beginning of a line element chain
            CSMP_FEM_TYPE etype = parseFiniteElementTypeEnum( (*eit) );
            if ( isLineElement( etype ) ) {
                  // beginnings or endings of polylines inside of model
                  if ( pfverts[elmt_idx][0] < 0 )
                    line_elmts.insert( elmt_idx );
                  // boundary line elements or loop on the boundary
                  if ( bflags[ plist[elmt_idx][0] ] < 0 && bflags[ plist[elmt_idx][1] ] < 0 )
                    boundary_line_elmts.insert( elmt_idx );
              }
            // if this is a surface element  
            else {
                 // looping over the faces (edges), finding those at the boundary
                 // and storing their nodes in correct order for later comparisons
                 if ( isTriangularElement( etype ) ) {
                      const size_t faces{3};
                      for ( size_t i{0}; i < faces; ++i )
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
                 else if ( isQuadrilateralElement( etype ) ) {
                      const size_t faces{4};
                      for ( size_t i{0}; i < faces; ++i )
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
                      // inserting the faces
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
      // ------------------------------------------------------------------------------------------------------
      // root-elmt & numbers of interconnected line elements in discovered chain
      map<int64_t ,deque<int64_t> > polylines;
      set<int64_t>                processed_elmts;
      
      // looping over the line elements that are missing one neighbor, i.e., are at the beginning of a chain
      for ( set<int64_t>::const_iterator it=line_elmts.begin(); it!= line_elmts.end(); ++it )
        {
           // if this is an isolated line segment with no neighbors it is skipped
           elmt_idx = (*it);
           if ( (pfverts[elmt_idx][0] < 0 && pfverts[elmt_idx][1] < 0) ||
                 processed_elmts.find(elmt_idx) != processed_elmts.end() ) continue;

           // if the line element has no first neighbor it must be at the beginning of a chain and correctly oriented
           BOX_BOUNDARY bflag = static_cast<BOX_BOUNDARY>( pfverts[elmt_idx][0] );
           // its neighbor-free side must be at an internal or external boundary
           if ( bflag >= 0 ) {
                cerr <<"\nVData::CreateConsistentLineElementOrientations2D: ";
                cerr <<"'pfvert["<< elmt_idx <<"][0]' entry for neighbour element 1: "<< pfverts[elmt_idx][0];
                cerr <<" not correct: ["<< pfverts[elmt_idx][0] <<","<< pfverts[elmt_idx][1] <<"].";
             }
           assert( bflag >= MULTIPLE );

          // storing the element as the first in the line element sequence
          pair<map<int64_t ,deque<int64_t> >::iterator,bool>
            chain_it=polylines.insert( make_pair( elmt_idx, deque<int64_t>{elmt_idx} ) );
            
          // making sure that the element was indeed inserted (else it is a duplicate)
          assert( chain_it.second == true );
          processed_elmts.insert( elmt_idx );
          
          // traversing the line element chain in the direction of available neighbors, node0 is the one with no neighbor
          // NB: in line elements the element neighbor also is opposite to the node with the same number, e.g.,
          //     [nd0]-line element0-[nd1]-line element1-[nd2] -> elm1 1 is opposite to nd0
          bool end_of_polyline(false);
          while ( end_of_polyline == false )
            {
               // is the neighbor element is correctly oriented its first node will be shared with the second node of the previous line element
               assert( pfverts[elmt_idx][1] >= 0 );
               // if the next neighbor's first neighbor element is element 'elmt_idx', everything is fine and no flip is required,
               bool flip = ( pfverts[ pfverts[elmt_idx][1] ][0] == elmt_idx ) ? false : true;
               // else, we move to this next element,
               elmt_idx = pfverts[elmt_idx][1];
               // and flip its nodes and neighbors.
               if ( flip == true ) {
                    const size_t node0 = plist[elmt_idx][0];
                    const size_t node1 = plist[elmt_idx][1];
                    plist[elmt_idx][0] = node1;
                    plist[elmt_idx][1] = node0;
                    int64_t  swap          = pfverts[elmt_idx][1];
                    pfverts[elmt_idx][1] = pfverts[elmt_idx][0];
                    pfverts[elmt_idx][0] = swap;
                 }
               // Then we store this element in the polyline.
               (*chain_it.first).second.push_back( elmt_idx );
               processed_elmts.insert( elmt_idx );
               // exit condition (if both nodes od line element are o the BOX_BOUNDARY)
               if ( pfverts[elmt_idx][0] < 0 || pfverts[elmt_idx][1] < 0 )
                 end_of_polyline = true;
            }
            
       } // for line_elements
       
      // NOTE: since every chain has a beginning and an end, it would normally stored twice, but:
      //       - of the interior chains, only the ones starting with the lower element number are kept
      //       - for the ones surrounding the model the one consistent with the counter-clockwise numbering of the higher-dim. ele is stored
      //       - the only line elements left now, are those forming part of loops, these must be consistent in their orientation with
      //         the faces of the elements that they surround.
      
      
      // 3. Reordering chains that are located on the model boundary to make them consistent with the counter-clockwise element numbering
      // --------------------------------------------------------------------------------------------------------------------------------
      // using: //  face-nd-ids,      elmt, face   to verify that elements are indeed oriented correctly
      //        map<set<size_t>,pair<size_t,size_t> > surf_elmt_face_nd_ids;
      //
      // getting the surface element deque ready for binary_search
      sort( surf_elmt_face_nd_ids.begin(), surf_elmt_face_nd_ids.end() );
      // looping over the line elements that are missing one neighbor, i.e., are at the beginning of a chain
      for ( set<int64_t>::const_iterator it=boundary_line_elmts.begin(); it!= boundary_line_elmts.end(); ++it )
        {
           // searching for the corresponding face of a higher dimensional element
           // --------------------------------------------------------------------
           // if a surface element face with same node numbering is found the line element is already correctly oriented
           if ( binary_search( surf_elmt_face_nd_ids.begin(), surf_elmt_face_nd_ids.end(), make_pair( plist[*it][0], plist[*it][1]) ) )
             continue;
           // if the face has the opposite orientation, the line element is flipped
           if ( binary_search( surf_elmt_face_nd_ids.begin(), surf_elmt_face_nd_ids.end(), make_pair( plist[*it][1], plist[*it][0]) ) ) {
                // swapping nodes
                int64_t  swap     = plist[*it][1];
                plist[*it][1]   = plist[*it][0];
                plist[*it][0]   = swap;
                // swapping neighbors
                swap            = pfverts[*it][1];
                pfverts[*it][1] = pfverts[*it][0];
                pfverts[*it][0] = swap;
             }
           // this is a line element with no surface element next to it?
           else {
                cerr <<"\nCreateConsistentLineElementOrientations: line element "<< *it <<" at border with the nodes:\n\t\t";
                cerr << plist[*it][0] <<"("<< parseBoundary(intToBOX_BOUNDARY(bflags[plist[*it][0]])) <<"), ";
                cerr << plist[*it][1] <<"("<< parseBoundary(intToBOX_BOUNDARY(bflags[plist[*it][1]])) <<"), ";
                cerr <<" has no higher-dimensional neighbor; its orientation is left untouched.\n";
             }
           
        } // end boundary_line_elmts
              
  } // end CreateConsistentLineElementOrientations





 
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
     
  } // end angleBetweenLineSegments
 
 
 
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
                  pair<map<int64_t ,set<int64_t> >::iterator,bool> // node 0
                    it0 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][0], set{elmt_idx} ) );
                  // if there is already an entry for the node, the elmt id is added to the set
                  if ( !it0.second )
                    (*it0.first).second.insert( elmt_idx );
                  
                  pair<map<int64_t ,set<int64_t> >::iterator,bool> // node 1
                    it1 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][1], set{elmt_idx} ) );
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
        for ( deque<vector<int64_t> >::const_iterator
              pft=pfverts.begin(); pft!=pfverts.end(); ++pft, ++elmt_idx )
          if ( !isLineElement( parseFiniteElementTypeEnum( pelmt[elmt_idx] ) ) )
            {
               size_t boundaries_per_element(0U);
               for ( vector<int64_t>::const_iterator pt=(*pft).begin(); pt!=(*pft).end(); ++pt )
                 // the face is on the boundary
                 if ( (*pt) < 0 ) boundaries_per_element++;
               if ( boundaries_per_element > 1U &&
                   isTriangularElement( parseFiniteElementTypeEnum( pelmt[elmt_idx] ) ) )
                 {
                   cerr <<"\n\n\telement "<< elmt_idx <<" ("<< parseFiniteElementType( pelmt[elmt_idx] ) <<") ";
                   cerr <<" has "<< boundaries_per_element <<" faces on model boundary.\n";
                   csmp_error.notice( WARNING, "VData::EstablishElementConnectivity2D:",
                                     "triangular element with  2 faces on boundary ");
                   triangle_with_all_nodes_on_boundary = true;
                 }
            }
        
        // 3. reconnecting line elements
        // -----------------------------
        // map<size_t,set<size_t> >  line_elmt_that_contain_node;
        for ( const auto& it : line_elmt_that_share_node )
          {
              const size_t n_connections(it.second.size()-1);
              
              // 1. isolated line elements terminating either at an inside node (INTERNAL) or at the BOX_BOUNDARY
              // --------------------------------------------------------------------------------------------------
              if ( n_connections == 0U ) {
                   const size_t elmt = (*it.second.begin());
                   // if there is no neighbor element corresponding to the first node
                   if ( it.first == plist[elmt][0] ) {
                        // identifying the boundary that the missing neighbor is located at
                        pfverts[elmt][0] = (bflags[ it.first ]==0) ? INTERNAL : bflags[ it.first ];
                     }
                   else if ( it.first == plist[elmt][1] ) {
                        // the missing neighbor is located at a boundary
                        pfverts[elmt][1] = (bflags[ plist[elmt][1] ]==0) ? INTERNAL : bflags[ plist[elmt][1] ];
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
                   std::vector<int64_t>    joint_line_elmts( it.second.begin(), it.second.end() ); // actual element ids
                   const size_t            n_elmts_to_combine(2U);
                   deque<vector<int64_t> > combinations;
                   if ( createUniqueCombinations( joint_line_elmts, n_elmts_to_combine, combinations ) == 0 )
                     csmp_error.notice( ERROR, "EstablishElementConnectivity2D", "no combinations between elements available");
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
                            
          } // processing the line element neighbors

     } // if plist empty

   assert( !pfverts.empty() );

   // 4. Elements with all nodes on the boundary can degrade solver convergence for certain boundary conditions
   //    This method eliminates such elements by node swapping
   if ( triangle_with_all_nodes_on_boundary ) SwitchCornerTriangles2D();

   // 5. reorienting line-element chains (done in other method)
   CreateConsistentLineElementOrientations2D();

 } // end EstablishElementConnectivity2D









/**
    Reconnects triangular elements with 3 nodes on the model boundary by switching nodes with their only neighbor; @note needs valid 'pfverts.'
    
    @return the number of replaced triangles;
*/
size_t VData::SwitchCornerTriangles2D()
 {
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
                   for ( size_t i{0}; i<3; ++i ) {
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
                        array elmt1_nds{ plist[elmt_idx][0], plist[elmt_idx][1], plist[elmt_idx][2] };
                        array elmt2_nds{ plist[nb_idx][0],   plist[nb_idx][1],   plist[nb_idx][2] };
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
                        array elmt1_nbors{ pfverts[elmt_idx][0], pfverts[elmt_idx][1], pfverts[elmt_idx][2] };
                        array elmt2_nbors{ pfverts[nb_idx][0],   pfverts[nb_idx][1],   pfverts[nb_idx][2] };
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
         csmp_error.notice( WARNING, "VData::EstablishElementConnectivity3D:", "supplied cell vector is empty; nothing was done." );
         return;
      }
    if ( Faces() > 0 ) {
         csmp_error.notice( WARNING, "VData::EstablishElementConnectivity3D:", "Face objects not handled yet" );
         return;
      }
    if ( InterFaces() > 0 ) {
         csmp_error.notice( WARNING, "VData::EstablishElementConnectivity3D:", "InterFace objects not handled yet" );
         return;
      }
 
    // 1. making separate search vectors of face keys for volume, surface and line elements
    // ------------------------------------------------------------------------------------
    cout << "\nEstablishElementConnectivity3D: Establishing CSMP FE neighbor connectivity...\n";
    cout << "  Building a multimap of the faces of the cells...\n";
    // face-node key  element idx, face number
    map<set<size_t>,map<size_t,size_t> >  volume_neighbor_keys, surface_neighbor_keys;
    //  node key, elements connected at that node
    map<size_t,set<size_t> >  line_elmt_that_share_node;

    const size_t n_elements(plist.size());
    pfverts.resize( plist.size() );

    for ( size_t elmt_idx{0}; elmt_idx < n_elements; ++elmt_idx )
      {
         // getting the element type (unfortunately this is known only at runtime)
         const CSMP_FEM_TYPE etype{ pelmt[elmt_idx] };
         assert( parseFiniteElementTypeEnum( pelmt[elmt_idx] ) != UNKNOWN );
         
          const size_t faces(CSMP_ElementSpecifications::FacesPerElementOfType(etype));
          pfverts[elmt_idx].resize(faces,IRREGULAR);
          for ( size_t face=0U; face<faces; ++face )
            {
               // creating face key of node pointers from indices of face nodes
               set<size_t> key;
               const size_t nodes(CSMP_ElementSpecifications::NodesPerFaceForElementOfType( etype, face ) );
               for ( size_t j=0U; j<nodes; ++j ) {
                    // inserting the global  node numbers into the key
                    const size_t face_node = plist[elmt_idx][ CSMP_ElementSpecifications::FaceNodeForElementOfType( etype, face, j ) ];
                    key.insert( face_node );
                 }
               // inserting newly generated keys into map
               if ( CSMP_ElementSpecifications::VolumeElement( etype ) ) {
                    pair<map<set<size_t>,map<size_t,size_t> >::iterator,bool>
                      vit = volume_neighbor_keys.insert( make_pair( key, map<size_t,size_t>{{elmt_idx,face}} ) );
                    if ( !vit.second ) (*vit.first).second.insert( make_pair(elmt_idx,face) );
                 }
               else if ( CSMP_ElementSpecifications::SurfaceElement( etype ) ) {
                    pair<map<set<size_t>,map<size_t,size_t> >::iterator,bool>
                      vit = surface_neighbor_keys.insert( make_pair( key, map<size_t,size_t>{{elmt_idx,face}} ) );
                    if ( !vit.second ) (*vit.first).second.insert( make_pair(elmt_idx,face) );
                 }
               else {
                    // for all line elements
                    // node 0
                    pair<map<size_t,set<size_t> >::iterator,bool>
                      it0 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][0], set{elmt_idx} ) );
                    // if there is already an entry for the node, the elmt id is added to the set
                    if ( !it0.second ) (*it0.first).second.insert( elmt_idx );
                    // node 1
                    pair<map<size_t,set<size_t> >::iterator,bool>
                      it1 = line_elmt_that_share_node.insert( make_pair( plist[elmt_idx][1], set{elmt_idx} ) );
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
        // map<size_t,set<size_t> >  line_elmt_that_contain_node;
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
                   const size_t            n_elmts_to_combine(2U);
                   deque<vector<int64_t> > combinations;
                   if ( createUniqueCombinations( joint_line_elmts, n_elmts_to_combine, combinations ) == 0 )
                     csmp_error.notice( ERROR, "EstablishElementConnectivity2D", "no combinations between elements available");
                   // finding inter-element angle for all combinations
                   //             angle, combination number
                   vector<pair<double,size_t> > inter_element_angles;
                   inter_element_angles.reserve( combinations.size() );
                   size_t n_combi{0};
                   for ( auto cit : combinations ) {
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
        for ( map<set<size_t>,map<size_t,size_t> >::const_iterator
              it=surface_neighbor_keys.begin(); it!=surface_neighbor_keys.end(); ++it )
          {
             assert( !(*it).second.empty() );
             // 2.2.1 if there is only a single entry, the face is at the model boundary
             // ------------------------------------------------------------------------
             if ( (*it).second.size() == 1U ) {
                  const size_t elmt_idx = (*(*it).second.begin()).first;
                  const size_t face     = (*(*it).second.begin()).second;
                  // finding out on which boundary the edge of the face is
                  // relying on appropriate box-boundary flagging
                  // (using only the two first nodes of the face)
                  const CSMP_FEM_TYPE etype{ pelmt[elmt_idx] };
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
                   std::vector<int64_t>  joint_surf_elmts;
                   joint_surf_elmts.reserve( (*it).second.size() );
                   for ( const auto& i : (*it).second ) joint_surf_elmts.push_back( i.first ); // actual element ids
                   const size_t           n_elmts_to_combine(2U);
                   deque<vector<int64_t> > combinations;
                   if ( createUniqueCombinations( joint_surf_elmts, n_elmts_to_combine, combinations ) == 0 )
                     csmp_error.notice( ERROR, "EstablishElementConnectivity3D", "no combinations between elements available");
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
                         // TODO: one could narrow down which boundary this is, but it will not be used later
                         pfverts[unassigned_elmt][boundary_face] = IRREGULAR;
                      }
               }
          }
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( !volume_neighbor_keys.empty() )
      {
        cout << "\n\t\tvolume elements...\n";
        for ( map<set<size_t>,map<size_t,size_t> >::const_iterator
              it=volume_neighbor_keys.begin(); it!=volume_neighbor_keys.end(); ++it )
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
          
         // 2.3.3 Remeshing tetrahedra that span a corner of the model
         // ----------------------------------------------------------
         // replacing the corner and inner elements with 3 new tetrahedra
         const size_t n_elements{pelmt.size()};
         // if this a tetrahedral only mesh no checks have to be performed
         if ( !HybridElementTypeMesh() && isTetrahedral( parseFiniteElementTypeEnum(ElementType(0)) ) )
           {
              for ( size_t i{0}; i<n_elements; ++i ) {
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
              for ( size_t i{0}; i<n_elements; ++i ) {
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
          
      } // volume elements
 
 } // end EstablishElementConnectivity3D






/**
      Finds the neighbors of each node and returns them into the argument vector..
      
      This method is equivalent to creating a sparsity pattern for matrix accumulation.
      
      @test OK SKM 8/12/21
*/
void VData::EstablishNodeNeighborConnectivity( std::vector<set<size_t>>& pnode ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( plist.empty() ) {
         csmp_error.notice( WARNING, "VData::EstablishNodeNeighborConnectivity:", "plist is empty; nothing was done." );
         return;
      }

    if ( !pnode.empty() ) pnode.clear();
    pnode.resize( px.size() );
    
    bool higher_order_elements{false};
 
    const size_t n_elements{ plist.size() };
    // looping over the elements to get the corner nodes of their segments
    for ( size_t elmt{0}; elmt < n_elements; ++elmt ) {
         // getting the element type
         const auto CSMP_FE_type{ pelmt[elmt] };
         if ( !higher_order_elements && CSMP_ElementSpecifications::InterpolationOrder(CSMP_FE_type) > 1 ) {
              csmp_error.notice( ERROR, "VData::EstablishNodeNeighborConnectivity", "connectivity of midside nodes not tested yet; check!" );
              higher_order_elements = true;
           }
         // for each segment
         const size_t n_segments{ CSMP_ElementSpecifications::SegmentsPerElementOfType( CSMP_FE_type ) };
         for ( size_t segm_id{0}; segm_id < n_segments; ++segm_id ) {
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
             const size_t n_nodes = CSMP_ElementSpecifications::NodesPerElementOfType(CSMP_FE_type);
             for ( size_t segm_id{0}; segm_id < n_segments; ++segm_id ) {
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
             const auto CSMP_FE_type{ pelmt[elmt] };
             if ( CSMP_ElementSpecifications::InterpolationOrder(CSMP_FE_type) == 2 ) {
                 // dealing with quadratic elements that have midside nodes
                 // relying on the numbering convention that midside nodes follow the corner nodes in the same order
                 // and that there is one midside node per segment
                 const size_t n_nodes = CSMP_ElementSpecifications::NodesPerElementOfType(CSMP_FE_type);
                 // for each segment
                 const size_t n_segments{ CSMP_ElementSpecifications::SegmentsPerElementOfType( CSMP_FE_type ) };
                 const size_t first_midside_node = n_nodes - n_segments;
                 // assigning the segment corner nodes as neighbors of the midside node
                 for ( size_t segm_id{0}; segm_id < n_segments; ++segm_id ) {
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
                  csmp_error.notice( ERROR, "VData::EstablishNodeNeighborConnectivity", "connectivity of midside nodes for O>2 meshes not done yet" );
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






/**
    vertex manifolds: pairs of nodes and their INSIDE,OUTSIDE, MIDDLE classifers
    typedef std::deque<std::set<std::pair<size_t,int8_t> > > vertexManifoldIndices;
    
    Checks for collocated vertices into transfer data structure.
*/
bool VData::ExtractNodeManifolds( vertexManifoldIndices& indexes ) const
 {
    indexes.clear();
    
    // if there are no interfaces, method returns false
    if (  first_interface_ == first_face_ ) {
         cerr <<"\nVData::ExtractNodeManifolds: does not contain any manifolds.\n";
         return false;
      }
      
   // creating manifold data from the interface node indices stored in plist
   //for ( auto it=PelmtInterfacesBegin(); it!=PelmtEnd(); ++it ) - element type info
   for ( auto it=PlistInterFacesBegin(); it!=PlistInterFacesEnd(); ++it )
     {
        // the first half of the interface vertices represents inside nodes
        assert( !(*it).empty() );
        const size_t plist_entries((*it).size());
        assert( (plist_entries & 1) == 0 /* even number */ );
        const size_t iface_nodes(plist_entries/2);
        
        // creating a new entry in the manifold map or getting an iterator to an existing one
        // map<size_t,set<pair<size_t,int8_t> > > vertexManifoldIndices
        for ( size_t i=1U; i<iface_nodes; ++i ) {
              const pair<size_t,int8_t> vertex( make_pair( (*it)[i], INSIDE ) );
              // trying the insertion
              auto manif_it = indexes.insert( make_pair( (*it)[i], set<pair<size_t,int8_t> >({vertex}) ) );
              // if this is an existing record the vertex is added to that one
              if ( manif_it.second == false )
                (*manif_it.first).second.insert( vertex );
          }
        for ( size_t i=iface_nodes; i<plist_entries; ++i ) {
               const pair<size_t,int8_t> vertex( make_pair( (*it)[i], OUTSIDE ) );
              // trying the insertion
              auto manif_it = indexes.insert( make_pair( (*it)[i], set<pair<size_t,int8_t> >({vertex}) ) );
              // if this is an existing record the vertex is added to that one
              if ( manif_it.second == false )
                (*manif_it.first).second.insert( vertex );
          }
     }
     
   // now trying to match the vertices of intervening elements, if any
   deque<int64_t> intervening_elmts;
   for ( auto it=next(PfvertsBegin(),first_interface_); it!=PfvertsEnd(); ++it ) {
         assert( (*it).size() > 2U );
         // the last entry in each pfvert record is the index of the intervening element or bflag
         const ssize_t idx = (*it).back();
         if ( idx >= 0 ) // if there is an intervening element
           intervening_elmts.push_back( idx );
     }
     
   // if there are intervening elements, they are added to the manifolds
   if ( !intervening_elmts.empty() ) {
        for ( auto elmt : intervening_elmts )
          for ( auto it=PlistBegin(elmt); it!=PlistEnd(elmt); ++it ) {
              const pair<size_t,int8_t> vertex( make_pair( (*it), MIDDLE ) );
              // trying the insertion
              auto manif_it = indexes.insert( make_pair( (*it), set<pair<size_t,int8_t> >({vertex}) ) );
              // if this is an existing record the vertex is added to that one
              if ( manif_it.second == false )
                (*manif_it.first).second.insert( vertex );
          }
     }
      
   return true;
   
 } // end ExtractNodeManifolds




/**
Assuming the following numbering:

First tetra
0 1 2 3 - ony one neighbor = nbr

Second tetra
0 3 2 4 -  4 neigbors: 0=, 1=, 2=, 3=cnr

The new elements are formed
0 1 2 4 - replacing cnr,             new neighbors: 0=new nbor, 1=nbor1, 2=new, 3=outside
0 3 1 4 - new,                           new neighbors: 0=new nbor, 1=cnr, 2=nbor2, 3=outside
1 3 2 4 - replacing neighbor,    new neighbors: 0=nbor0, 1=nbor1, 2=new, 3=outside

The corresponding neighbor elements are.

*/
void splitCornerTetrahedron( VData& vdata, size_t cnr, size_t nbr )
 {
    assert( cnr < vdata.Elements() );
    assert( nbr < vdata.Elements() );
    
    // 1. creation of space for new element (assuming that there are no Faces of InterFaces)
    // -------------------------------------------------------------------------------------
    assert( vdata.Faces() == 0 );
    assert( vdata.InterFaces() == 0 );
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
    
    // 2. assignment of nodes
    // ----------------------
    // new element first
    vdata.Plist( n_elements_new-1U, 0, vdata.Plist( cnr, 0 ) );
    vdata.Plist( n_elements_new-1U, 1, vdata.Plist( cnr, 3 ) );
    vdata.Plist( n_elements_new-1U, 2, vdata.Plist( cnr, 1 ) );
    vdata.Plist( n_elements_new-1U, 3, vdata.Plist( nbr, 3 ) );
    // temp numbers
    const long cnr_n0{vdata.Plist(cnr,0)}, cnr_n1{vdata.Plist(cnr,3)},
               cnr_n2{vdata.Plist(cnr,1)}, cnr_n3{vdata.Plist(nbr,3)},
               nbr_n0{vdata.Plist(cnr,1)}, nbr_n1{vdata.Plist(cnr,3)},
               nbr_n2{vdata.Plist(cnr,2)}, nbr_n3{vdata.Plist(nbr,3)};
    // reshaped corner element
    vdata.Plist( cnr, 0, cnr_n0 );
    vdata.Plist( cnr, 1, cnr_n1 );
    vdata.Plist( cnr, 2, cnr_n2 );
    vdata.Plist( cnr, 3, cnr_n3 );
    // reshaped neighbor element
    vdata.Plist( nbr, 0, nbr_n0 );
    vdata.Plist( nbr, 1, nbr_n1 );
    vdata.Plist( nbr, 2, nbr_n2 );
    vdata.Plist( nbr, 3, nbr_n3 );

    // 3. assignment of neigbors
    // -------------------------
    // new element first
    vdata.Pfvert( n_elements_new-1U, 0, nbr );
    vdata.Pfvert( n_elements_new-1U, 1, cnr );
    vdata.Pfvert( n_elements_new-1U, 2, vdata.Pfvert( nbr, 2 ) );
    vdata.Pfvert( n_elements_new-1U, 3, vdata.Pfvert( cnr, 3 ) ); // outside
    // temporaries
    const long cnr_p0(nbr), cnr_p1{vdata.Pfvert(nbr,1)},
               cnr_p2(n_elements_new-1), cnr_p3{vdata.Pfvert(nbr,3)},
               nbr_p0{vdata.Pfvert(nbr,0)}, nbr_p1{vdata.Pfvert(nbr,1)},
               nbr_p2(n_elements_new-1), nbr_p3{vdata.Pfvert(cnr,3)};
    // reshaped corner element
    // new neighbors: 0=new nbor, 1=nbor1, 2=new, 3=outside
    vdata.Pfvert( cnr, 0, cnr_p0 );
    vdata.Pfvert( cnr, 1, cnr_p1 );
    vdata.Pfvert( cnr, 2, cnr_p2 );
    vdata.Pfvert( cnr, 3, cnr_p3 ); // outside
    // reshaped neighbor element
    // new neighbors: 0=nbor0, 1=nbor1, 2=new, 3=outside
    vdata.Pfvert( nbr, 0, nbr_p0 );
    vdata.Pfvert( nbr, 1, nbr_p1 );
    vdata.Pfvert( nbr, 2, nbr_p2 );
    vdata.Pfvert( nbr, 3, nbr_p3 ); // outside
    
    cout <<"\nsplitCornerTetrahedron: replaced corner "<< cnr;
    cout <<" and its neighbor "<< nbr <<", adding the new tetrahedron "<< n_elements_new;

 } // end splitCornerTetrahedron

 
} // end namespace csmp
 
 
