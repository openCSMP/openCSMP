#include "BoundaryInterface.h"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"
#include "FaceConstructionData.h"
#include "Face.h"
#include "Boundary.h"
#include "Region.h"
#include "Box.h"
#include "VSet.h"
#include "Exception.h"
#include "SmallSet.h"
#include "ErrorHandler.h"
#include "variableOperations.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundaryInterface()
  {}

template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundaryInterface( const BoundaryInterface& bd )
:faceBoundaryMap_( bd.faceBoundaryMap_ )
  {
  }

template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
BoundaryInterface<dim,BOUNDARY_COMPLEX>::~BoundaryInterface()
  {
  }

template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
typename std::map<std::string,csmp::Boundary<dim> >::iterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesBegin()
  { return faceBoundaryMap_.begin(); }


template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
typename std::map<std::string,csmp::Boundary<dim> >::iterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesEnd()
  { return faceBoundaryMap_.end(); }


template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
typename std::map<std::string,csmp::Boundary<dim> >::const_iterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesBegin() const
  { return faceBoundaryMap_.begin(); }


template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
typename std::map<std::string,csmp::Boundary<dim> >::const_iterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesEnd() const
  { return faceBoundaryMap_.end(); }


template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
typename std::map<std::string,csmp::Boundary<dim> >::iterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary( const csmp::Boundary<dim>& bref )
  { 
    for( boundaryIterator it( BoundariesBegin() ); it != BoundariesEnd(); ++it )
      if( &(it->second) == &bref )
        return it;
    return BoundariesEnd();
  }

template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
size_t  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundaries() const
  { return faceBoundaryMap_.size(); }



template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool  BoundaryInterface<dim,BOUNDARY_COMPLEX>::ContainsBoundary( const std::string& bname ) const
 {
    if ( faceBoundaryMap_.find(bname) != faceBoundaryMap_.end()  )
        return true;
    return false;
 }

/// returns a reference to the boundary 'bname' if it exists
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
Boundary<dim>&  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary( const std::string& bname )
  { 
    boundaryIterator  bit=faceBoundaryMap_.find( bname );
    if( bit != faceBoundaryMap_.end() )
      return bit->second;  
    else
      {
        std::string errMsg("Boundary does not exist!");
        errMsg.append( " (" + bname  + ")" );
        throw csmp::Exception( ERROR,
                               "BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary(std::string&)",
                               errMsg.c_str() );
      }
  } 


template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
const Boundary<dim>&  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary( const std::string& bname ) const
  {    
    boundaryConstIterator  bit=faceBoundaryMap_.find( bname );
    if( bit != faceBoundaryMap_.end() )
      return bit->second;  
    else
      {
        std::string errMsg("Boundary does not exist!");
        errMsg.append( " (" + bname + ")" );
        throw csmp::Exception( ERROR,
                               "BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary(std::string&)",
                               errMsg.c_str() );
      }
  } 


// TODO: SKM: this method appears to make temporary copies of boundaries rather than using move semantics
/**
    @author Philipp Lang
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim,BOUNDARY_COMPLEX>::SwitchBoundaryNames( csmp::Boundary<dim> const& b1, csmp::Boundary<dim> const& b2 )
{
  // since we use a map we have to copy, delete and insert
  csmp::Boundary<dim> b1Copy(b1);
  csmp::Boundary<dim> b2Copy(b2);

  // finding, deleting b1
  boundaryIterator itMatch( faceBoundaryMap_.end() );
  for( boundaryIterator it = faceBoundaryMap_.begin(); it != faceBoundaryMap_.end(); ++it )
      if( &it->second == &b1 )
        itMatch = it;

  if( itMatch == faceBoundaryMap_.end() )
    throw csmp::Exception( ERROR,
                           "BoundaryInterface::SwitchBoundaryNames",
                           "Internal error 1" );

  std::string b1Name = itMatch->first;
  faceBoundaryMap_.erase(itMatch);

  // finding, deleting b1
  itMatch = faceBoundaryMap_.end();
  for( boundaryIterator it = faceBoundaryMap_.begin(); it != faceBoundaryMap_.end(); ++it )
    if( &it->second == &b2 )
      itMatch = it;

  if( itMatch == faceBoundaryMap_.end() )
    throw csmp::Exception( ERROR,
                           "BoundaryInterface::SwitchBoundaryNames",
                           "Internal error 2" );

  std::string b2Name = itMatch->first;
  faceBoundaryMap_.erase(itMatch);
  
  faceBoundaryMap_.insert( make_pair( b2Name, b1Copy) );
  faceBoundaryMap_.insert( make_pair( b1Name, b2Copy) );
} 




/**
     Returns name of the first Boundary (patch) that contains the supplied
     region name strings or any other other strings in an arbitrary order.
     
     if the boundary does not exist, method returns \O string.
     
     @author SKM 
     @date March 2016
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
std::string  BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryName( const std::set<std::string>& intersected_regions ) const
 {
    const BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this) );
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryName:",
                                         "supplied set of substrings is empty; returning '\0'." );
         return std::string("\0");
      }
    // if the model has no boundaries
    if ( boundaryComplex.Boundaries() == 0 ) {
         ErrorHandler::Instance().notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryName:",
                                         "model has no boundaries; returning '\0'." );
         return std::string("\0");
      }
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=boundaryComplex.BoundariesBegin(); it!=boundaryComplex.BoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
            // if the substring is found
            if ( (*it).first.find(*ir) !=std::string::npos ) substrings_found++;
            // when all substrings are contained in the boundary name, it is returned
          if ( substrings_found == substrings_used_in_search )
            return (*it).first;
       }
    return std::string("\0");
 }



/**
    Searches the model for boundaries the name of which contains the search strings
    provided via the first set. The results are returnd into the second set.
    
    @note Use this method, for example, to retrieve multiple boundary patches that were generated
    from a single lower-dimensional regon, like a fault surface.

     @author SKM 
     @date March 2016
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
size_t BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryNames( const set<string>& intersected_regions,
                                                                   std::set<std::string>& region_patches_found ) const
 {
    const BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this) );
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().notice( WARNING, "findBoundaries:", "supplied set of substrings is empty; returning '\0'." );
         return 0U;
      }
    // if the model has no boundaries
    if ( boundaryComplex.Boundaries() == 0 ) {
         ErrorHandler::Instance().notice( WARNING, "findBoundaries:", "model has no boundaries; returning '\0'." );
         return 0U;
      }
    region_patches_found.clear();
   
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=boundaryComplex.BoundariesBegin(); it!=boundaryComplex.BoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
            // if the substring is found
            if ( (*it).first.find(*ir) !=std::string::npos ) substrings_found++;
          // when all substrings are contained in the boundary name, it is returned
          if ( substrings_found == substrings_used_in_search )
            region_patches_found.insert( (*it).first );
       }
    // return how many region patches contain the search string(s)
    return region_patches_found.size();
 }






/**
     Checks whether the name contains any substring that is indicative of a model boundary,
     this includes LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK, IRREGULAR and any name that 
     is preceded by the string BOUNDARY.
     
     @note BOX_BOUNDARY edge and corner names are not checked for
 
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::IsBoundaryName( const std::string& regionName ) const
  {
    // TAG: BOUNDARY should be at the beginning of the region name
    // DON'T confuse with SPLITBOUNDARY
    std::size_t found_position = regionName.find( "BOUNDARY" );
    if ( found_position == 0 && found_position!=std::string::npos )
        return true;
    if( regionName == "BOTTOM" )
      return true;
    if( regionName == "RIGHT" )
      return true;
    if( regionName == "LEFT" )
      return true;
    if( regionName == "TOP" )
      return true;
    if( regionName == "FRONT" )
      return true;
    if( regionName == "BACK" )
      return true;
    if( regionName == "IRREGULAR" )
      return true;
    return false;
  }




/**
     LOCAL METHOD (not known beyond this compilation unit)
 
     creates underscore-separated unique names for the region patches based on the juxtapositions relationships
     across the lower-dimensional regions, the names are composed of:
 
     1. the name of the master region
     
     2. "BOUNDARY"
     
     3. the patch identifier number attached to boundary
     
     3. the name of the inner region, i.e. the region that the lower-dimensional element normals point away from
     
     4. the name of the outer region, i.e. that into which the normals point
     
     @attention  where the boundary just intersects a layer (same material on either side), the layer name appears
     only once. The second instance is replaced by INTERSECTION.
*/
std::string createNameOfInternalBoundaryFrom( const FaceConstructionData& fdata, const std::vector<std::string>& region_names )
 {
     assert( fdata.ElementMaterial() < region_names.size() );
     std::string boundary_name( region_names[ fdata.ElementMaterial() ] );
     boundary_name += "_BOUNDARY";
     boundary_name += to_string(fdata.PatchNumber());
     boundary_name += '_';
     pair<long,long> materials(fdata.Materials());
     assert( materials.first  < region_names.size() );
     assert( materials.second < region_names.size() );
     boundary_name += region_names[ materials.first ];
     boundary_name += '_';
     if ( materials.first == materials.second ) boundary_name +="INTERSECTION";
     else boundary_name += region_names[ materials.second ];
   
     return boundary_name;
 }



/**
    Checks whether a model is box-shaped while
    we can't ask RectangularShapedModel anymore since boundary regions were replaced by boundaries.
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim, BOUNDARY_COMPLEX>::BoxShaped() const
{
	const BOUNDARY_COMPLEX<dim>& boundaryComplex(static_cast<const BOUNDARY_COMPLEX<dim>& >(*this));

	std::cout << "\nBoundaryInterface<" << dim << ">::BoxShaped: checking whether the model contains the boundaries LEFT, RIGHT...";
	size_t boundaries(6);
	if (dim == 2)
		boundaries = 4;
	set<BOX_BOUNDARY> boundariesFound;
	BOX_BOUNDARY currentBoundary(NOT);
	const Region<dim>& region(boundaryComplex.Region("Model"));
	const typename std::vector<Node<dim>*>::const_iterator nodesEnd(region.NodesEnd());
	for (typename std::vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != nodesEnd; ++it)
	{
		currentBoundary = (*it)->AtBoundary();
		if (currentBoundary != LEFT and currentBoundary != RIGHT and currentBoundary != BACK
			and currentBoundary != FRONT and currentBoundary != BOTTOM and currentBoundary != TOP
			and currentBoundary != IRREGULAR)
			continue;
		boundariesFound.insert(currentBoundary);
		if (boundariesFound.size() >= boundaries)
		{
			std::cout << "...true\n";
			return true;
		}
		// TODO: dealing with the case where there might be one IRREGULAR boundary instead of the TOP boundary
	}
	std::cout << "...false\n";
	return false;
}



/**
    Creates boundary from already interconnected faces that also know their parent elements.
    
    @note the supplied Face pointer vector is moved into boundary and will therefore not be 
    accessible anymore after this method has been called.
    
    @author SKM
    @date 1/4/2016
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::AddBoundary( const char* boundary_name,
                                                           typename std::vector<Face<dim>*>::iterator facesBegin,
                                                           typename std::vector<Face<dim>*>::iterator facesEnd,
                                                           BOX_BOUNDARY bflag )
 {
    BOUNDARY_COMPLEX<dim>* const boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>* const>(this) );
    assert( boundaryComplex != nullptr );

    // inserting boundary if it does not existing yet
    std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
        it = faceBoundaryMap_.insert( std::make_pair( boundary_name, move(csmp::Boundary<dim>( boundary_name,
                                                                                               boundaryComplex->Database(),
                                                                                               facesBegin, facesEnd, bflag ) ) ) );
    if ( it.second ) {
         cout << "\nBoundaryInterface<"<< dim <<">::AddBoundary: successfully created boundary '";
         cout << boundary_name <<"' from input faces.";
      }
    else {
         ErrorHandler&  csmp_error( ErrorHandler::Instance() );
         csmp_error.notice( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::AddBoundary:",
                            boundary_name, "Boundary already exists or other problem arose. Nothing was done.");
         return false;
      }

     return true;
   
 } // end AddBoundary




/**
     higherDimensionalNeighbors() - finds the higher-dim neighbor elements of 
     a dim-1 element embedded within the higher-dim mesh.
     
     @attention both neighbors have to be present for this to work.
     
     = LOCAL METHOD ONLY KNOWN TO THIS COMPILATION UNIT
 
     - finds the IDs of the higher dimensional neibors of the current element
 
     - identifies which of the neighbors is on the inside and which on the outside
       as indicated by the normal direction of the lower dimensional element
 
     - identifies the materials on either side
     
     @return two inside-outside pairs of element idx numbers and corresponding materials on either side
     all the data are stored in the returnd FaceConstructionData object.
     
     assumptions
     - assumes that the nodes and elements in the entire model domain are numbered continuously
     
     application
     - use this function for finding neighbors of a surface element that sits on the inside of another region
     
     TODO: output the numbers of matching faces of the discovered the elements so that they can later be connected
 
*/
template<size_t dim>
FaceConstructionData  higherDimensionalNeighbors( const Element<dim>& e, const csmp::Index& mtrl_key )
 {
     assert( e.IsSurfaceElement() );

     // 1. looping over the parent elements of the nodes searching for the faces which are shared with the lower dimensional element
     // -----------------------------------------------------------------------------------------------------------------------------
     // making a set of element nodes to later identify faces by comparison
     set<size_t>   node_set, test_set;
     const size_t  nodes(e.Nodes());
     for ( size_t i=0U; i<nodes; ++i ) node_set.insert(e.N(i)->Idx());
     map<const Element<dim>*,size_t>  nbor_elmts;
     vector<size_t> fnids;
     for ( size_t i=0U; i<nodes; i++ ) {
          const size_t parents(e.N(i)->Parents());
          for ( size_t j=0U; j<parents; ++j ) {
               const Element<dim>* const eptr(e.N(i)->Parent(j));
               const size_t faces(eptr->Faces());
               for ( size_t k=0U; k<faces; ++k ) {
                     eptr->FE()->NodesOfFace( k, fnids );
                     size_t fnodes(fnids.size());
                     for ( size_t l=0U; l<fnodes; ++l )
                       test_set.insert( eptr->N( fnids[l])->Idx() );
                     // if the face is shared the element and its face are recorded
                     if ( node_set == test_set ) {
                          // storing a pointer to this element and its local face number
                          // making sure that no duplicate is received
                          nbor_elmts.insert( make_pair(eptr,k) );
                       }
                     test_set.clear();
                 }
            }
       }
    assert( nbor_elmts.size() == 2U );
   
    // 2. finding which of the neighbors is the inside one by projecting face normals onto lower dim element normal
    // -------------------------------------------------------------------------------------------------------------
    vector<double64>  enrml, fnrml;
    e.UnitNormal( enrml );
    typename map<const Element<dim>*,size_t>::const_iterator  nbit(nbor_elmts.begin());
    bool inside_elmt_found(false);
    bool outside_elmt_found(false);
    pair<size_t,size_t> nbors;
    pair<size_t,size_t> faces;
    pair<long,long>     materials;

    // first element
    // -------------
    assert( (*nbit).first != nullptr );
    faces.first = (*nbit).second;
    (*nbit).first->UnitNormalToFace( faces.first, fnrml );
    double64 dotproduct(0.);
    for ( size_t k=0U; k<dim; ++k )
      dotproduct += enrml[k] * fnrml[k];
   
    // if the projection is negative, the first element lies on the outside
    if ( dotproduct < 0. ) {
         nbors.second       = (*nbit).first->Idx();
         materials.second   = static_cast<long>((*nbit).first->Read( mtrl_key ));
         outside_elmt_found = true;
      }
    else {
         nbors.first        = (*nbit).first->Idx();
         materials.first    = static_cast<long>((*nbit).first->Read( mtrl_key ));
         inside_elmt_found  = true;
      }
    nbit++;


// TODO: SKM: this could be done without 2 normal projections for speedup

    // second element
    // --------------
    assert( (*nbit).first != nullptr );
    faces.second = (*nbit).second;
    (*nbit).first->UnitNormalToFace( faces.second, fnrml );
    dotproduct = 0.;
    for ( size_t k=0U; k<dim; ++k )
      dotproduct += enrml[k] * fnrml[k];

   // if the projection is negative, the second element lies on the outside
    if ( dotproduct < 0. ) {
         // checking that we have no duplication here
         assert( outside_elmt_found == false );
         nbors.second     = (*nbit).first->Idx();
         materials.second = static_cast<long>((*nbit).first->Read( mtrl_key ));
      }
    else {
         assert( inside_elmt_found == false );
         nbors.first      = (*nbit).first->Idx();
         materials.first  = static_cast<long>((*nbit).first->Read( mtrl_key ));
      }
   
    // initialise with nbors, their faces, adjacent materials, and patch numbers
    return FaceConstructionData( e.Idx(), nbors,  faces, materials,
                                 static_cast<long>(e.Read( mtrl_key)) );
   
 } // end higherDimensionalNeighbors





/**
     higherDimensionalNeighbor() - finds a higher-dimensional element, one face of which
     matches  the input supplied lower-dimensional element.
     
     @return pointer to the higher-dimensional adjacent element or NULL when not found.
     
     @return face of the higher-dimensional element that matches the lower-dim element
     
     @return material ID of the higher-dimensional element; if it does not exist, NaN is returned.
     @note an undefined material key is accepted, but in this case NaN will be returned.
 
     The discovered element is considered to be located on the inside of the lower-dim element,
     hence its normal ought to be pointing away from it, else there is an orientation problem.
 
     @attention assumptions
     - assumes that the nodes and elements in the entire model domain are numbered continuously
     
     application
     - use this function for finding the higher-dimensional neighbor of a surface element that sits on the
       outside boundary of the model
 
     @test SKM 22/8/2018 - fixed a bug where element returned had lower spatial dimensional than supplied
     element.
 
*/
template<size_t dim>
const Element<dim>* const higherDimensionalNeighbor( const Element<dim>& e, const csmp::Index& mtrl_key,
                                                     size_t& local_face_number_of_e, double64& material_ID  )
 {
     if      ( dim == 3 ) assert( e.IsSurfaceElement() );
     else if ( dim == 2 ) assert( e.IsLineElement() );
     assert( mtrl_key.type == SCALAR );
     assert( mtrl_key.place == ELEMENT || mtrl_key.place == UNDEFINED );

     // 1. looping over the parent elements of the nodes searching their faces for ones that are shared with the lower dimensional element
     // ----------------------------------------------------------------------------------------------------------------------------------
     // making a set of element nodes to later identify faces by comparison
     set<size_t>   node_set, test_set;
   
     const size_t  nodes(e.Nodes());
     for ( size_t i=0U; i<nodes; ++i ) node_set.insert(e.N(i)->Idx());
   
     const csmp::Element<dim>*  nbor_elmt(nullptr);
     vector<size_t>             fnids;
   
     // since the same element may be discovered by each of the face nodes
     // the loop is stopped after the first discovery
     for ( size_t i=0U; i<nodes; i++ ) {
          const size_t parents(e.N(i)->Parents());
          for ( size_t j=0U; j<parents; ++j ) {
               const Element<dim>* const eptr(e.N(i)->Parent(j));
               // only if the element is not the same and also of a different type
               if ( eptr != &e and
                    eptr->FE_Type() != e.FE_Type() and
                    eptr->Nodes() >= e.Nodes() )
                 {
                   const size_t faces(eptr->Faces());
                   for ( size_t k=0U; k<faces; ++k ) {
                         eptr->FE()->NodesOfFace( k, fnids );
                         size_t fnodes(fnids.size());
                         for ( size_t l=0U; l<fnodes; ++l )
                           test_set.insert( eptr->N( fnids[l] )->Idx() );
                         // if the face is shared the element and its face are recorded
                         if ( node_set == test_set ) {
                              // storing the pointer to this element and its local face number
                              // making sure that no duplicates are received
                              nbor_elmt = eptr;
                              local_face_number_of_e = k;
                              break;
                           }
                         test_set.clear();
                     }
                 }
            }
       }
    assert( nbor_elmt != nullptr );

     if ( dim == 3 ) {
         if (e.IsVolumeElement())
           assert( nbor_elmt->IsVolumeElement() );
         else if (e.IsSurfaceElement())
           assert( nbor_elmt->IsVolumeElement() );
         else if (e.IsLineElement())
           assert( nbor_elmt->IsSurfaceElement() or  nbor_elmt->IsVolumeElement() );
     }
     if ( dim == 2 ) {
         if (e.IsSurfaceElement())
           assert( nbor_elmt->IsSurfaceElement() );
         else if (e.IsLineElement())
           assert( nbor_elmt->IsSurfaceElement() );
     }
   
    // 2. drawing the results
    // -------------------------------------------------------------------------------------------------------------
    material_ID = ( nbor_elmt != nullptr && mtrl_key.place != UNDEFINED )
                  ? nbor_elmt->Read(mtrl_key) : std::numeric_limits<double64>::quiet_NaN();
   
    return move(nbor_elmt);
   
 } // end higherDimensionalNeighbor

// STUB
template<>
const Element<1U>* const higherDimensionalNeighbor( const Element<1U>& e, const csmp::Index&, size_t&, double64& )
 {
    throw logic_error("higherDimensionalNeighbor(in BoundaryInterface: there should be no boundaries in 1D model");
    return &e;
 }


template const Element<2U>* const higherDimensionalNeighbor( const Element<2U>&, const csmp::Index&, size_t&, double64& );
template const Element<3U>* const higherDimensionalNeighbor( const Element<3U>&, const csmp::Index&, size_t&, double64& );

// TESTING
/*
if ( dim == 2 && nbor_elmt->IsLineElement() ) {
     cerr <<"\nhigherDimensionalNeighbor: potentially found duplicate edge elements; candidates are:\n";
     e.Out();
     cerr <<"\n\nand:\n";
     nbor_elmt->Out();
  }
*/   



/** 
     Converts lower-dimensional region into Boundarie(s) of faces, decomposed into patches; 
     region is moved from "Model" to non-unique regions, connectivity is updated.
     
     Uses node-to-parent relationship to find the higher dimensional elements that will 
     share a face with the Face:  usng    higherDimensionalNeighbors()
     
     Algorithmic steps:

     1. Verify input lower-dimensional region object from which the boundary shall be created: 
        must be lower dimensional and must lie inside of model
 
     2. Determine the number of boundary segments (sub-boundaries) that the new boundary will consist of.
        The output of this step will be a map of FaceConstructionData in which the names of the 
        new boundary segments are the keys.
 
     3. In the MeshManager object, 
 
        3.1 create all required face objects
 
        3.2 connect them with one another (neighbors); Boundary::EstablishNeighborConnectivity( std::vector<Face<dim>*>& ); 
            this is important because any ModelSubDomain creation relies on this connectivity during 
            identification of interior and perimeter.
      
     4. Create the Boundary segments, one-by-one from the map<bname, FaceConstructionData>; 
        this involves connecting the faces to their higher dimensional neighbours.
 
     5. For the nodes located on the new boundaries, update / recreate the parent element vectors 
        so that these no longer include neither the elements from which the Boundary was created nor the new faces.
 
     6. Remove the parent region of the new boundary(ies) from the region “Model.” 
        The region is also put into the non-unique region map.
      
     7. Assign BOX_BOUNDARY flags to new nodes if any.
     
     8. By default, but optional removes parent region (including its elements).

*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
size_t BoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateInternalBoundaryFrom( const char* dim_1_region, bool remove_original_region )
 {
    BOUNDARY_COMPLEX<dim>& model( static_cast<BOUNDARY_COMPLEX<dim>&>(*this) );

    cout <<"\nBoundaryInterface<"<< dim <<",BOUNDARY_COMPLEX>::CreateInternalBoundaryFrom: forming boundary(ies) from region: '"<< dim_1_region <<"'...\n";
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    // 1. Verify input lower-dimensional region object from which the boundary shall be created: must be lower dimensional and must lie inside of model
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    // does the parent region exist
    if ( model.ContainsRegion(dim_1_region) == false ) {
          ErrorHandler::Instance().notice( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region, "does not exist; nothing was done." );
          return 0;
      }
    // do such boundaries already exist ?
    const set<string> intersected_regions{dim_1_region};
    set<string> pre_existing_boundaries;
    if ( FindBoundaryNames( intersected_regions, pre_existing_boundaries ) > 0 ) {
          string error_info;
          for ( auto it=pre_existing_boundaries.begin(); it!=pre_existing_boundaries.end(); ++it ) {
               error_info += (*it);
               error_info +=", ";
            }
          ErrorHandler::Instance().notice( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:",
                                           error_info.c_str(), "boundaries are already contained in this model." );
          return 0;
      }
    // does the model contain unique regions
    if ( model.UniqueRegions() < 1 ) {
          ErrorHandler::Instance().notice( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", "model contains no unique regions; cannot proceed." );
          return 0;
      }
    // verifying that we are indeed dealing with a region of surface elements only and that their normals all point into same direction
    Region<dim>&  subdomain(model.Region(dim_1_region));
    if ( checkNeighborNormalsForConsistentOrientation( subdomain ) == false ) {
         ErrorHandler::Instance().notice( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region,
                                         "region appears to have inconstent surface-normal orientations; nothing was done." );
         return 0;
      }
    size_t boundary_elements(0);
    for ( auto eit=subdomain.ElementsBegin(); eit!=subdomain.ElementsEnd(); ++eit ) {
         if ( (*eit)->AtBoundary() != NOT and (*eit)->AtBoundary() != INTERNAL and (*eit)->AtBoundary() != IRREGULAR ) boundary_elements++;
      }
    if ( boundary_elements > 0 ) {
         ErrorHandler::Instance().notice( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region, "region appears to lie at the model boundary; nothing was done." );
         return 0;
      }
    // creating region labels and tagging the regions with unique integer indentifiers
    const string region_tag("region identifier");
    if ( !model.Database().IsDefined(region_tag.c_str()) )
      model.CreateProperty( region_tag.c_str(), "X", SCALAR, ELEMENT );
    const csmp::Index mtrl_key = model.Database().StorageKey(region_tag.c_str());
    vector<string>  region_names;
    const size_t model_regions = model.CountAndLabelRegions( region_tag.c_str(), region_names );
   
    if ( model_regions == 1 )
       ErrorHandler::Instance().notice( INFO, "BoundaryInterface::CreateInternalBoundaryFrom:", region_tag.c_str(), "is single valued; so there is only one patch." );

 
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 2. Determine number of boundary segments (sub-boundaries) that the new boundary will consist of.
    //    The output of this step will be a map of FaceConstructionData in which the names of the new boundary segments are the keys.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // renumbering nodes and elements
    model.Region("Model").UpdateMemberIndexes();
   
    // looping over the region, identifying and recording the juxtaposition relationships
    map<pair<long,long>,size_t>   patches;
    vector<FaceConstructionData>  face_construction_data;
    map<size_t,string>            patch_names;
    string                        patch_name;
    size_t                        n_juxtapositions(0);

    // 2.1 looping over lower dimensional region identifying juxtaposition relationships
    // ----------------------------------------------------------------------------------------
    face_construction_data.reserve(subdomain.Elements());
    for ( auto eit=subdomain.ElementsBegin(); eit!=subdomain.ElementsEnd(); ++eit )
      {
          // 2.1.1 identifying neighbors, facing relations, and juxtaposed materials for current element
          FaceConstructionData  fdata(higherDimensionalNeighbors( *(*eit), mtrl_key ));
        
          // 2.1.2 recording which category of juxtaposition element fall into, naming it and assigning a patch number
          pair<map<pair<long,long>,size_t>::iterator,bool>  it=patches.insert( make_pair(fdata.Materials(),n_juxtapositions) );
          // incrementing number of juxtapositions and corresponding patch names
          if ( it.second == true ) {
               fdata.PatchNumber( (*it.first).second );
               patch_name = createNameOfInternalBoundaryFrom( fdata, region_names );
               patch_names.insert( make_pair(n_juxtapositions,patch_name) );
               n_juxtapositions++;
            }
          fdata.PatchNumber( (*it.first).second );
        
          // 2.1.3 recording the data for the element that will later be used to construct the face from
          face_construction_data.push_back( fdata );
      }
    assert( face_construction_data.size() == subdomain.Elements() );
   
   
    // 2.2 creating labeled boundary patches from the face-defining data
    // -----------------------------------------------------------------
    // 2.2.1 making a map 'patch_numbers' from 'patch_names' to search for patch identifiers
    map<string,size_t>  patch_numbers;
    for ( auto it=patch_names.begin(); it!=patch_names.end(); ++it )
      patch_numbers.insert( make_pair( (*it).second, (*it).first ) );
   
    // 2.2.2 building new map where the patch faces are organised by patch names
    map<string,vector<FaceConstructionData> > patch_simplexes;
    vector<FaceConstructionData>              empty_vec;
    for ( auto it=patch_names.begin(); it!=patch_names.end(); ++it )
      patch_simplexes.insert( make_pair( (*it).second, empty_vec ) );
   
    // 2.2.3 inserting the patch identifiers into the vectors in the map
    for ( auto pit=patch_simplexes.begin(); pit!=patch_simplexes.end(); ++pit )
      {
         assert( patch_numbers.find((*pit).first) != patch_numbers.end() );
         const size_t patch_number((*patch_numbers.find((*pit).first)).second);
         // reserving storage
         (*pit).second.reserve(face_construction_data.size());
         // looping over all face data assigning the ones that are suitable
         for ( auto it=face_construction_data.begin(); it!=face_construction_data.end(); ++it )
           if ( (*it).PatchNumber() == patch_number )
             (*pit).second.push_back( (*it) );
      }
    // 2.2.4 trimming excess storage of the face-data vectors
    for ( auto pit=patch_simplexes.begin(); pit!=patch_simplexes.end(); ++pit )
      vector<FaceConstructionData>( (*pit).second ).swap( (*pit).second );

 
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 3. Getting the MeshManager object to create Face objects for all boundary patches at the same time
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    //  3.1 creating the required face objects
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    Region<dim>&  model_domain(model.Region("Model"));
    const size_t new_faces_required(subdomain.Elements());
    vector<Face<dim>*> face_vector;
    face_vector.reserve(new_faces_required);
    const size_t original_faces(model.Mesh().Faces());

    // establish the storage requirements for face variables
    const LocalVariables             lvsFaces( model.Database().LocalVariablesAt(FACE) );
    const IntegrationPointVariables  lvsIntegrationPoints( model.Database().IntegrationPointVariablesAt(FACE) );
    vector<vector<Face<dim>*> >       face_ptr_per_patch(patch_simplexes.size());
   
    size_t patch_counter(0);	
    for ( map<string,vector<FaceConstructionData> >::const_iterator
          it=patch_simplexes.begin(); it!=patch_simplexes.end(); ++it ) {
        face_ptr_per_patch[patch_counter].reserve( (*it).second.size() );
        // for each of the new patches
        for ( vector<FaceConstructionData>::const_iterator
              pit=(*it).second.begin(); pit!=(*it).second.end(); ++pit ) {
             // creating the faces
             // ------------------
             // storing pointers to the new faces in the vector from which the boundary will be constructed
			Face<dim> newface(*model_domain.E((*pit).Element()), model_domain.E((*pit).InnerElement()), model_domain.E((*pit).OuterElement()), lvsFaces, lvsIntegrationPoints);
			 Face<dim>* faceObj = model.Mesh().Add(newface);
			 
			 // the first face is assigned into the root face of this face group in the mesh			 
			 // the faces should be connected each other, otherwise each root face has only a single face
			 if (face_ptr_per_patch[patch_counter].size() == 1) {
				 model.Mesh().SetRootFace(face_vector.front());
			 }

			 face_vector.push_back(faceObj);

			 // remembering which faces make up the patch
             face_ptr_per_patch[patch_counter].push_back( face_vector.back() );
           }
         patch_counter++;
      }
    patch_simplexes.clear();
    cout << "\n\tAdded "<< model.Mesh().Faces() - original_faces <<" to mesh.\n";

   
    //  3.2 connect them with one another (neighbors); Boundary::EstablishNeighborConnectivity( std::vector<Face<dim>*>& ); this is important because
    //      any ModelSubDomain creation relies on this connectivity during identification of interior and perimeter.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    cout << "\n\tEstablishing neighbor connectivity among faces as it is needed to build the boundaries...\n";
    // 3.2.1 building search map for face neighbors
    // --------------------------------------------
    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<size_t,Face<dim>*> >  surface_neighbor_keys, line_neighbor_keys;
    vector<size_t>  fnids;
    set<Node<dim>*>  key;

    for ( typename vector<Face<dim>*>::const_iterator it=face_vector.begin(); it!= face_vector.end(); ++it )
      for ( size_t face=0U; face<(*it)->Faces(); face++ ) 
        {
           // creating face key from idx's of face
           (*it)->FE()->NodesOfFace( face, fnids ); 
           for ( size_t j=0U; j<fnids.size(); j++ ) key.insert( (*it)->N( fnids[j] ) );
           // inserting newly generated keys into multimap
           if ( (*it)->IsSurfaceElement() )
             surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else // for all line elements
             line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           key.clear();
        }

    // 3.2.2 building face neigborhoods
    // --------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    cout << "\n\tBuilding face neighbor connectivity...";

    // 3.2.2.1 line faces
    // ---------------------
    if ( !line_neighbor_keys.empty() ) {
        Face<dim>* e1Ptr(nullptr);
        Face<dim>* e2Ptr(nullptr);
        cout << "\n\t\tline elements...";
        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Face<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                             it2(line_neighbor_keys.begin());
        it2++;
        while ( it2 != line_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first )
                {
                   assert( (*it1).second.second != nullptr );
                   assert( (*it2).second.second != nullptr );
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment
                   // assigning eachothers faces
                   //     face pointer                   nbor face idx  neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }
              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // line faces
    
    // 3.2.2.2 surface faces
    // ------------------------
    if ( !surface_neighbor_keys.empty() ) {
        Face<dim>* e1Ptr(NULL);
        Face<dim>* e2Ptr(NULL);
        cout << "\n\t\tsurface elements...";
        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Face<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                              it2(surface_neighbor_keys.begin());
        it2++;
        while ( it2 != surface_neighbor_keys.end() ) {
              if ( (*it1).first == (*it2).first )
                {
                   assert( (*it1).second.second != nullptr );
                   assert( (*it2).second.second != nullptr );
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   ++it1;
                   ++it2;
                }
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // etablish neighbors of surface faces
   
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 4. Create the Boundary segments, one-by-one from the map< bname, FaceConstructionData >
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // using map<size_t,string>  patch_names   from above
    for ( size_t i=0U; i<patch_names.size(); ++i )
       // creating the boundary patch
       AddBoundary( patch_names[i].c_str(), face_ptr_per_patch[i].begin(), face_ptr_per_patch[i].end(), INTERNAL );

    face_ptr_per_patch.clear();
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 5. If the Face objects were constructed from lower-dimensional Elements - for nodes located on the new boundary,
    //    update / recreate the parent element vectors of the nodes on the boundary so that these no longer include
    //    neither the elements from which the Boundary was created nor the new faces.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    map<Element<dim>*,size_t>  parents_to_keep;
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         const size_t parent_elements((*nit)->Parents());
         // copying those node parent pointers to the temporary vector which shall be kept
         for ( size_t i=0U; i<parent_elements; ++i ) {
              if ( (*nit)->Parent(i)->IsSurfaceElement() and
                   subdomain.Contains( (*nit)->Parent(i) ) )
                continue;
              else
                parents_to_keep.insert( make_pair( (*nit)->Parent(i), (*nit)->ParentNodeNumber(i) ) );
           }
         // resetting the parent element storage of current node
         (*nit)->EraseParents();
      
         // rebuilding parent element storage of current node
         for ( auto it=parents_to_keep.begin(); it!=parents_to_keep.end(); ++it )
           // parent element node-number, parent element
           (*nit)->Assign( (*it).second, (*it).first  );
         parents_to_keep.clear();
      }
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 6. Assign BOX_BOUNDARY flags to the nodes of each new patch by using the underlying region
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) (*nit)->AtBoundary(INTERNAL);
 
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 7. Remove parent region of the boundary from “Model” and into the non-unique list of regions so that it does not get included into computations.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    model.RemoveFromRegion( "Model", dim_1_region );
    model.MoveToNonUniqueRegions( dim_1_region );
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 8. (optional) remove parent region (including its elements) if no longer required.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    if ( remove_original_region ) {
         // NOTE: works, but then the VTU interface cannot output the boundary anymore
         const bool delete_elements(true);
         model.RemoveRegion( dim_1_region, delete_elements );
      }
   
    return patch_names.size();

 } // end CreateInternalBoundaryFrom







// ================================================================================================================
// OLD CODE - TODO: CHECK WHETHER IT IS STILL RELEVANT
// ================================================================================================================

// TODO: ask the MeshManager to delete the required range of Faces
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim, BOUNDARY_COMPLEX>::RemoveBoundary( csmp::Boundary<dim>& boundary, bool deleteElements )
  {
    assert( deleteElements == false );
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>* >(this) );

    if( deleteElements )
    {
        // Faces removed from all existing boundaries and erased from MeshManager
        std::vector<Face<dim>* > elementsToDelete( boundary.SimplexVector().begin(), boundary.SimplexVector().end() );

        if( !elementsToDelete.empty() )
        {
            // Detach Faces from Neighbors
            // SKM FIX boundary.DetachElementsFromNeighbors();

            // remove redundant Faces from existing Boundaries
            for( boundaryIterator
                 bit = faceBoundaryMap_.begin();
                 bit != faceBoundaryMap_.end();){

                if( removeVectorElements( (*bit).second.SimplexVector(), elementsToDelete ) > 0 )
                {
                    if( !(*bit).second.SimplexVector().empty() )
                    {
                        (*bit).second.CreateNodePointerVector();
                        (*bit).second.IdentifyPerimeter( );
                        ++bit;
                    }
                    else
                        faceBoundaryMap_.erase( bit++ );
                }
                else
                {
                    // we need to reestablish the perimeter because
                    // the order of SimplexVector() was changed by removeVectorElements ()
                    (*bit).second.IdentifyPerimeter( );
                    ++bit;
                }
            }

            // update all indices
            boundaryComplex->UpdateIndices();

            return;
        }

    } // deleting elements

    // finding the boundary
    typename std::map<std::string,csmp::Boundary<dim> >::iterator iterBoundary( faceBoundaryMap_.end() );
    for( boundaryIterator it = faceBoundaryMap_.begin(); it != faceBoundaryMap_.end(); ++it )
      if( &it->second == &boundary )
        iterBoundary = it;

    if( iterBoundary == faceBoundaryMap_.end() )
        throw csmp::Exception( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::RemoveBoundary", "boundary does not exist" );

    // if the boundary was found in the list, it is erased
    faceBoundaryMap_.erase( iterBoundary );
  }




// -----------------------------------------------
// Binary input/output
// -----------------------------------------------



/**
    Writes all the boundaries to binary file, remembering all details of their implementation.
    
    @attention method expects that a unique global nnumbering of Face and Node objects was established prior
    to its application.

    @attention ONLY METHOD that provides complete information about each boundary:
    interior vs. perimeter nodes and element; boudary faces etc.
    
    @author SKM 3/4/2016
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::OutputAllBoundariesToBinary( const char* file_name ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    const BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this) );
    const PropertyDatabase<dim>& database( boundaryComplex.Database() );

    std::string bin_file(file_name);
    fstream fp(bin_file.c_str(), ios::out | ios::binary);
    if ( !fp.is_open() ) {
         csmp_error.notice( ERROR, "BoundaryInterface::OutputAllBoundariesToBinary:",
                            bin_file, "file could not be opened; nothing was done." );
         return false;
      }
     // 0. writing the file header
   {
      BinaryFileSectionWrite sect(fp, "BNDFHEDR");
      std::string heading("BoundaryInterface::OutputAllBoundariesToBinary: ");
      heading +="boundary information for Model '";
      heading += boundaryComplex.Name();
      heading +="' to file: ";
      heading += bin_file;
      heading +="'.";
      skm_C_fwrite( fp, heading.c_str() );
   }

     std::cout <<"\nBoundaryInterface<"<< dim <<">::OutputAllBoundariesToBinary: boundaries written to binary file: ";
     // 1. writing all boundary objects to binary file
   {
     BinaryFileSectionWrite sect(fp, "BOUNDARY");

     const size_t records(this->Boundaries());
     fp.write( (char*) &records, sizeof(size_t) );

     for ( typename std::map<std::string,csmp::Boundary<dim> >::const_iterator
           git=BoundariesBegin(); git!=BoundariesEnd(); ++git )
       {
         BinaryFileSectionWrite hdr(fp, "ONE_BDRY");

          // 1.1 writing the entire connectivity structure to the binary file
          (*git).second.WriteDomainIndexesToBinaryFile( fp );
          // 1.2 writing the boundary flags
          int32 bflag = (*git).second.AtBoundary();
          const size_t record(1U);
		  fp.write((char*) &record, sizeof(size_t));
		  fp.write((char*) &bflag, sizeof(int32));
          // 1.3 writing the stored variables
          domainVariablesOut( fp, (*git).second, database );
          std::cout << (*git).first <<" ";
       }
   }

     // 2. writing the boundary complex variables here
   {
     BinaryFileSectionWrite sect(fp, "BOUNDVAR");

     domainVariablesOut( fp, boundaryComplex, database );
   }

    // 3. footer, and clean up
   {
     BinaryFileSectionWrite sect(fp, "BNDFFOTR");
   }

    fp.close();
    std::cout <<"\nBoundaryInterface<"<< dim <<">::OutputAllBoundariesToBinary: file '";
    std::cout << bin_file <<"' has been successfully written.\n";
   
    return true;

} // end OutputAllBoundariesToBinary




/**
     reads and initialises boundaries from file written by OutputAllBoundariesToBinary()
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim,BOUNDARY_COMPLEX>::InputAllBoundariesFromBinary( const char* file_name )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     std::string bin_file(file_name);
	 fstream fp(bin_file.c_str(), ios::in | ios::binary);
	 if (!fp.is_open()) {
          csmp_error.notice( ERROR, "BoundaryInterface::InputAllBoundariesFromBinary:",
                             bin_file, "file could not be opened; nothing was done." );
          return;
       }
   
     // 0. reading the file header and printing it to screen
   {
     BinaryFileSectionRead sect(fp, "BNDFHEDR");

     char  text[500U];
     skm_C_fread( fp, text );
     std::cout <<"\nBoundaryInterface<"<< dim <<">::InputAllBoundariesFromBinary: Reading file header:\n\t"<< text << std::endl;
   }
     std::cout <<"\n\timporting the boundaries: ";
   
     // 1. reading the boundaries
     const PropertyDatabase<dim>& database( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this).Database() );
     MeshManager<dim>& mesh( static_cast<BOUNDARY_COMPLEX<dim>& >(*this).Mesh() );
   {
     BinaryFileSectionRead sect(fp, "BOUNDARY");
     
     SubDomainInfo  info;
     size_t  records(0);  // region records
     // getting number of unique region records from file
     fp.read( (char*) &records, sizeof(size_t) );
     if ( records > 0 )
        // reading the regions sequentially
        for ( size_t i=0U; i<records; ++i )
          {
            BinaryFileSectionRead hdr(fp, "ONE_BDRY");

             // 1.1 reading name and face indices for each boundaries
             readDomainIndexesFromBinaryFile( dim, fp, info );
            
                  // 1.2 reading BOX boundary flag of the boundary
                  size_t record;
                  fp.read( (char*) &record, sizeof(size_t) );
                  assert( record == 1 );
                  int32 box_boundary_index(IRREGULAR_OUTSIDE);
                  fp.read( (char*) &box_boundary_index, sizeof(int32) );
                  BOX_BOUNDARY bflag = intToBOX_BOUNDARY( box_boundary_index );
                  std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
                    it=faceBoundaryMap_.insert( std::make_pair( info.name.c_str(), csmp::Boundary<dim>(database,mesh,info,bflag) ) );
                  //   ^^^^^^^^^^^^^^^
                  if ( !it.second )
                       throw csmp::Exception( FATAL_ERROR, "BoundaryInterface::InputAllBoundariesFromBinary:",
                                            info.name, "Boundary could not be formed; issue with binary file." );
               
                  // 1.3 reading the variable values
                  domainVariablesIn( fp, (*it.first).second, database );
               
                  // 1.4 reporting out
                  cout <<"\n\t\t"<< (*it.first).first <<" ("<< parseBoundary(bflag) <<", "<< (*it.first).second.Elements() <<" faces).";

				  if (info.name.c_str() == "LEFT")
					  cout << "error after this!!!\n";
         }
     else csmp_error.notice( WARNING, "BoundaryInterface::InputAllBoundariesFromBinary:",
                            bin_file, "does not contain any boundary descriptions; no boundaries were initialised." );
   }

   // 2. read boundary complex variables here
   {
     BinaryFileSectionRead sect(fp, "BOUNDVAR");
     BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>& >(*this) );
     PropertyDatabase<dim>& database( boundaryComplex.Database() );
     domainVariablesIn( fp, boundaryComplex, database );
   }

    // 3. cleaning up
   {
     BinaryFileSectionRead sect(fp, "BNDFFOTR");
   }

    fp.close();
    std::cout <<"\n\nBoundaryInterface<"<< dim <<",BOUNDARY_COMPLEX>::InputAllBoundariesFromBinary: file '";
    std::cout << bin_file <<"' has been read successfully.\n";

} // end InputAllBoundariesFromBinary




// ================================================================================================================
// ================================================================================================================
// ================================================================================================================
// ================================================================================================================
// OLD CODE - TODO: CHECK WHETHER IT IS STILL RELEVANT
// ================================================================================================================
// ================================================================================================================
// ================================================================================================================
// ================================================================================================================







/**
    Method forms a Boundary between the two supplied regions. This will involve the creation
    and connection of Faces.

    @attention:  BoundariesInterface cannot be created in 1D models or between regions which contain 
    one-dimensional elements.

    This method will only work for Regions which are not overlapping.

    This method makes no sense for the region 'Model' as it encompasses all unique regions.
    
    @note the BOX boundary flag will always be internal for this type of boundary because it 
    lies between higher dimensional regions.
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary( const char* group1, const char* group2, bool createRegionBetween )
 {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    std::string  b_normal( std::string( std::string(group1) + std::string("_") + std::string(group2) ) );

    std::string  region1(group1);
    std::string  region2(group2);

    if (region1 == "Model" or region2 == "Model") {
         csmp_error.notice( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(between):", "Region 'Model' not eligible for InsertBoundary.");
         return false;
      }
    if (region1 == region2) {
         csmp_error.notice( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(between):", "Provided Regions are identical.");
         return false;
      }


    if ( !boundaryComplex->IsUnique(group1) or !boundaryComplex->IsUnique(group2) ) {
         csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(between):",
                           "This method is intended for the creation of boundaries between unique Regions");
 
         // checking for a potential overlap of the regions, if the regions are non-unique
         if ( boundaryComplex->RegionIntersection( group1,  group2, "groupintersection" ) ) {
               csmp_error.notice( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(between):",
                                 "one of the supplied regions is not unique and they overlap",
                                 "It was therefore impossible to insert a boundary");
               boundaryComplex->RemoveRegion( "groupintersection", false );
               return false;
        }
      }
   
    const csmp::Region<dim>&  gref1(boundaryComplex->Region(group1));
    const csmp::Region<dim>&  gref2(boundaryComplex->Region(group2));
    
    // checking whether the two regions share some nodes (these will mark their common boundary)
    const size_t  shared_nodes(sharedNodes( gref1, gref2 ));
    
    bool succeeded(false);
    if ( shared_nodes == 0U )
      csmp_error.notice( WARNING,
                         "Model<dim,BOUNDARY_COMPLEX>::InsertBoundary(between):",
                        "The regions of interest do not share any nodes; trying to create a boundary");
   
    // attempt to create a regular (Face-based) boundary  
    else {
         std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
             it = faceBoundaryMap_.insert( std::make_pair( b_normal, csmp::Boundary<dim>( b_normal, boundaryComplex->Database(), INTERNAL ) ) );
         if ( it.second ) 
           {
             std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary(between): creating boundary between " << group1 << " and " << group2 << std::endl;
             boundaryComplex->UpdateIndices();
             if( createRegionBetween )
             {
                 boundaryComplex->RegionBetween( group1, group2, std::string( std::string(group1) + std::string("_") + std::string(group2) ).c_str() );
                 const csmp::Region<dim>&  gref( boundaryComplex->Region( std::string( std::string(group1) + std::string("_") + std::string(group2) ).c_str() ) );
                 succeeded = (*it.first).second.CreateFrom( boundaryComplex->Mesh(), gref, csmp::Index(), INTERNAL );
             }
             else
             {
                 succeeded = (*it.first).second.CreateBetween( boundaryComplex->Mesh(), boundaryComplex->FE_Manager(), gref1, gref2 );
             }
             boundaryComplex->UpdateIndices();
             std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary(between): created boundary between " << group1 << " and " << group2 << std::endl;
             return succeeded;
           }
         else throw csmp::Exception( INFO, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(between):",
                                     b_normal.c_str(),
                                    "boundary already exists. Nothing was done.");
      }
    return false; 
 } // end InsertBoundary





/**

  Method forms a Boundary (ModelSubDomain<Face>) from an existing Region<Element> of lower
  dimensional representation.

  @author P. Lang Aug 2011
  
  @todo !!! SKM: method fails when the same material is on either side of the boundary.
  in this case we get all kinds of normal orientations.
  
  SKM fix: Rewrite as follows:
  
  0. Test that the boundary actually consists of lower-dimensional elements and make sure that all normals 
     point in the same direction
  
  1. Normal orientations
     - The normals of the lower dimensional region are already all pointing in the same direction
      (this was verified when the model was built)
     - We therefore leave the normals alone and just name the boundaries consistent with the inside-outside relations
  
  2. Naming conventions
     - it is not sufficient to just use the names of the neighbouring regions:
     - where the material on the inside is the same as on the outside, we also need INSIDE and OUTSIDE as new keywords
     - if we stick with a single normal convention, we also need the name of the original region to get a unique name
       (think through a salt diapir surface example)
     - we also want the name of the original region so that we know which boundary patches belong together
       (for instance if we want to process entire faults turned into boundaries)
    suggestion: start the names of boundaries with that of the original region, e.g., fault_sand_shale, fault_INSIDE_sand
    
  3. Boundaries stored in the BoundaryInterface
    - we should now have unique names that can be searched for either by the name of the original region or the
      combination of lithologies that border each other
    - by starting the name with the name of the parent region, we ascertain that all boundary segments are stored sequentially
      in the boundary map
    - in addition or alternatively we might keep sets of the boundary patches that were created by InsertBoundary(), 
      at least we should output them to a set<boundaryName> for the user to check
      
   NB: I am not sure whether the side identification algorithm works correctly at layer boundaries; this is worth a check! 
 
*/
 /*
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary( const char* region, BOX_BOUNDARY boxBoundary, bool deleteRegionAndItsElements )
  {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // if region appears to be valid reference is created
    if( !boundaryComplex->ContainsRegion( region ) )
      csmp_error.notice( ERROR, "Model<dim,BOUNDARY_COMPLEX>::InsertBoundary(from dim-1 region):", region, "region does not exist." );
    // referencing the region of interest
    const csmp::Region<dim>&  rref( boundaryComplex->Region(region) );    
    // asserting the case where a 3D boundary is attempted to be built for a model containg volume elements
    if( dim == 3U and containsVolumeElements( boundaryComplex->Region( "Model" ) ) and !containsSurfaceElements( rref ) )
      csmp_error.notice( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(from dim-1 region):",
                                "Attempting to create a line element boundary for a mesh containing volumetric elements!" );
    // establishing name following convention "Face-RegionName", using BOX_BOUNDARY if supplied
    std::string regionName( region );
    std::string bName( regionName );
    // inserting boundary if not existing yet
    std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
        it = faceBoundaryMap_.insert( std::make_pair( bName, csmp::Boundary<dim>( bName, boundaryComplex->Database(), boxBoundary ) ) );
    if ( it.second )
      {
        std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary(from dim-1 region): creating boundary from Region "<< region << "\n";
        boundaryComplex->UpdateIndices();
        bool succeeded( (*it.first).second.CreateFrom( boundaryComplex->Mesh(), rref, csmp::Index(), boxBoundary ) );
        boundaryComplex->UpdateIndices();
        if( deleteRegionAndItsElements ) {
          // SKM FIX boundaryComplex->RemoveRegion( region );
             boundaryComplex->RemoveFromRegion( "Model", region );
             boundaryComplex->MoveToNonUniqueRegions( region );
          }
          
        std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary(from dim-1 region): created boundary from Region "<< region << "\n";
        return succeeded;
      }
    else csmp_error.notice( INFO, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(from dim-1 region):",
                            bName.c_str(), "Boundary already exists. Nothing was done.");
    // shouldn't get here
    return false;
  } // InsertBoundary

  */


/** Method forms a Boundary (ModelSubDomain<Face>) from nodes and elements flagged
    as box boundary as provided in parameter.

  @author P. Lang Aug 2011
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary( BOX_BOUNDARY boxBoundary, const char* region )
  {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // if the region appears to be valid a reference to it is created
    if ( !boundaryComplex->ContainsRegion(region) )
      csmp_error.notice( ERROR, "Model<dim,BOUNDARY_COMPLEX>::InsertBoundary(from dim-1 region):", region, "region does not exist." );
    const csmp::Region<dim>&  rref( boundaryComplex->Region(region) );
    // establishing name following convention "Face-RegionName", using BOX_BOUNDARY if supplied
    std::string regionName( parseBoundary(boxBoundary) );
    std::string bName( regionName );
    // inserting boundary if not existing yet
    std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
        it = faceBoundaryMap_.insert( std::make_pair( bName, csmp::Boundary<dim>( bName, boundaryComplex->Database(), boxBoundary ) ) );
    if ( it.second )
      {
        std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary(from dim-1 region): creating boundary " <<  parseBoundary( boxBoundary ) << std::endl;
        bool succeeded( (*it.first).second.CreateFrom( boundaryComplex->Mesh(), rref, csmp::Index(), boxBoundary ) );
        boundaryComplex->UpdateIndices();
        std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary(from dim-1 region): created boundary " <<  parseBoundary( boxBoundary ) << std::endl;
        return succeeded;
      }
    else csmp_error.notice( INFO, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary(from dim-1 region):",
                            bName.c_str(), "Boundary already exists. Nothing was done.");
    // shouldn't get here
    return false;
    
  } // end InsertBoundary



	/**
	Method forms a Boundary (ModelSubDomain<Face>) around a named existing Region<Element>.
	The string "BOUNDARY" is appended to the region name to identify the newly generated boundary.

	@author P. Lang Aug 2011
	*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim, BOUNDARY_COMPLEX>::AddFaces(const char* region)
{
	BOUNDARY_COMPLEX<dim>* boundaryComplex(static_cast<BOUNDARY_COMPLEX<dim>*>(this));
	ErrorHandler&  csmp_error(ErrorHandler::Instance());

	// if region appears to be valid reference is created
	if (!boundaryComplex->ContainsRegion(region))
		csmp_error.notice(ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::AddFaces", region, "Region does not exist.");
	const csmp::Region<dim>&  rref(boundaryComplex->Region(region));
	// establishing name following convention "RegionName_BOUNDARY"
	std::string bName(region);
	BOX_BOUNDARY bflag = (bName == "Model") ? IRREGULAR : INTERNAL;	
	bName += "_BOUNDARY";

    // inserting boundary if not existing yet
    std:: pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
        it = faceBoundaryMap_.insert( std::make_pair( bName, csmp::Boundary<dim>( bName, boundaryComplex->Database(), bflag ) ) );
    if ( it.second )
      {
        std::cout << "\nBoundaryInterface<"<< dim <<">::AddFaces: creating boundary around " << region << std::endl;
        //boundaryComplex->UpdateIndices(region);
        //                                 FACE & BOUNDARY CREATION
        bool succeeded( (*it.first).second.CreateAround( boundaryComplex->Mesh(), boundaryComplex->FE_Manager(), rref ) );
		assert( succeeded == true );
        //boundaryComplex->UpdateIndices(region);
        std::cout << "\nBoundaryInterface<"<< dim <<">::AddFaces: created boundary around " << region << std::endl;
        return succeeded;
      }
    else csmp_error.notice( INFO, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary",
                            bName.c_str(),
                            "Boundary already exists. Nothing was done.");
    return false;
  
  } // AddFaces




/** 
  Method forms a Boundary (ModelSubDomain<Face>) from supplied set of faces and name. 
 
  @attention This could create non-unique boundaries

  @author P. Lang
  @date April 2013
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::InsertBoundary( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                                                              const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                                                              const std::string& bName )
  {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );

    // inserting boundary if not existing yet
    std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
        it = faceBoundaryMap_.insert( std::make_pair( bName, csmp::Boundary<dim>( bName, boundaryComplex->Database(), parseBoundary(bName) ) ) );
    if ( it.second )
      {
        std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary creating boundary " <<  bName << std::endl;
        boundaryComplex->UpdateIndices();
        bool succeeded = (*it.first).second.CreateFrom( facesBegin, facesEnd );
        boundaryComplex->UpdateIndices();
        std::cout << "\nBoundaryInterface<"<< dim <<">::InsertBoundary created boundary " <<  bName << std::endl;
        return succeeded;
      }
    else throw csmp::Exception( ERROR, "Model<dim,BOUNDARY_COMPLEX>::InsertBoundary",
                                bName.c_str(),
                                "Boundary already exists. Nothing was done.");
    // shouldn't get here
    return false;
  } // end InsertBoundary



template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::DivideBoundary( typename std::map<std::string,csmp::Boundary<dim> >::iterator boundary,
                                                              const typename std::map<std::string,Region<dim> >::const_iterator subRegion )
  {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
	std::string  bName( subRegion->first );
	
	// rename it due to discontiguous model
	auto& meshMgr = boundaryComplex->Mesh();
	if( meshMgr.NodeGroups() > 1 )
		bName = (*boundary).second.Name() + "_" + subRegion->first;    

	const BOX_BOUNDARY bFlag((*boundary).second.AtBoundary());
	std::pair<typename std::map<std::string, csmp::Boundary<dim> >::iterator, bool>
		subBoundary = faceBoundaryMap_.insert(std::make_pair(bName, csmp::Boundary<dim>(bName, boundaryComplex->Database(), bFlag)));
	if (subBoundary.second)
	{
		if (boundary->second.Divide(subRegion->second, subBoundary.first->second)) {
			auto& face_vec = subBoundary.first->second.FaceVector();
			if (face_vec.size() > 0) boundaryComplex->Mesh().SetRootFace(face_vec[0]);
			std::cout << "\nBoundaryInterface<" << dim << ">::DivideBoundary created boundary " << bName << std::endl;
			return true;
		}
		else {
			faceBoundaryMap_.erase(subBoundary.first);
			csmp_error.notice(INFO, "BoundaryInterface::DivideBoundary", bName.c_str(), "Nothing to divide. Unable to form boundary.");			
		}
	}

    return false;    
  }




/**
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::DivideBoundary( typename std::map<std::string,csmp::Boundary<dim> >::iterator boundary,
                                                              const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                                                              const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                                                              const std::string& bName )
  {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const BOX_BOUNDARY bFlag( (*boundary).second.AtBoundary() );
    std::pair<typename std::map<std::string,csmp::Boundary<dim> >::iterator,bool>
      subBoundary = faceBoundaryMap_.insert( std::make_pair( bName, csmp::Boundary<dim>( bName, boundaryComplex->Database(), bFlag ) ) );
    if( subBoundary.second )
    {
      if( boundary->second.Divide( facesBegin, facesEnd, subBoundary.first->second ) )
      {
		auto& face_vec = subBoundary.first->second.FaceVector();
		if (face_vec.size() > 0) boundaryComplex->Mesh().SetRootFace(face_vec[0]);
        std::cout << "\nBoundaryInterface<"<< dim <<">::DivideBoundary created boundary " << bName << std::endl;
        return true;
      }
      else
        throw csmp::Exception( ERROR, "BoundaryInterface::DivideBoundary", "Unable to form boundary" );
    }
    else csmp_error.notice( INFO, "BoundaryInterface::DivideBoundary",
                            bName.c_str(),
                            "Boundary already exists. Nothing was done.");

    return false;    
  }



/**
    loops over the perimeter of the boundary, making keys from the node-pointers of the boundary face
    nodes and recording the boundary elements for the faces of which the keys were made
     
     @attention makes sense in 3D only
*/
void createPerimeterKeysFor( const Boundary<3U>& boundary, map<set<csmp::Node<3U>*>,Face<3U>*>& perimeter_keys )
 {
    if ( !perimeter_keys.empty() ) perimeter_keys.clear();
   
    // creating keys for the perimeter element faces of boundary1
    vector<size_t>  fnids;
    for ( size_t i=boundary.InteriorElements(); i<boundary.Elements(); ++i )
      for ( size_t j=0U; j<boundary.PerimeterFaces(i); ++j )
        {
            const size_t pface = boundary.PerimeterFace(i,j);
            boundary.E(i)->FE()->NodesOfFace( pface, fnids );
            set<Node<3U>*>  key;
            for ( size_t k=0U; k<fnids.size(); ++k )
              key.insert( boundary.E(i)->N( fnids[k] ) );
            perimeter_keys.insert( make_pair(key,boundary.E(i)) );
              
        }
   
 } // end createPerimeterKeysFor




/**
    Connects the line faces representing the edge with their neighbors    
    logic: where the line elements share a node they are connected
    
    @note makes sense only in 3D.
    
    @todo TODO: deal with the special case of manifolds. At the corner
    of a box-shaped model, for instance, 3 line-element faces contact each other.
    
    @attention where multiple line element faces connect at a Node,
    this method assigns a nullptr to that face.
    
    @attention this method assumes that the first 2 nodes of the line-element
    Face are the end-point nodes
*/
void createLineFaceConnectivity( std::vector<Face<3U>*>& line_faces )
 {
    if ( line_faces.empty() ) return;
   
    // 1. making a map of the parent faces that each node is connected to
    //  key       faces that are connected to the node (should be 2 at most)
    unordered_map<Node<3U>*,set<Face<3U>*> > parent_faces;
   
    for ( const auto it : line_faces ) {
         set<Face<3U>*> parents({it});
         for ( size_t i=0U; i<it->Nodes(); ++i ) {
              // inserting a new set or inserting a face pointer into the set if the node key already exists
              auto nit = parent_faces.insert( make_pair( it->N(i), parents ) );
              if ( !nit.second )
                (*nit.first).second.insert(it);
           }
      }
   
   // 2. connecting the faces with one another
   //   (the assumption is that each face has 1-2 equidimensional neighbors that coincide with its corner nodes!=midside nodes if any)
   ErrorHandler& csmp_error(ErrorHandler::Instance());

   for ( auto it : parent_faces ) {
        // 2.1 nodes/faces at the end-points of the edge = perimeter faces
        // (there may only be one neighbor or a manifold interpreted as endpoint)
        if ( it.second.size() != 2 ) {
          if ( it.second.size() > 2 ) {
               it.first->Out();
               csmp_error.notice( ERROR, "creatLineFaceConnectivity(Face):",
                                "edge node is connected to more than 2 line Faces;\
                                 don't know how to deal with this manifold.");
            }
           /* NOTHING NEEDS TO BE DONE BECAUSE FACE POINTERS ALREADY ARE NULLPTRs
             // this neighbor of the face is set to nullptr=no neighbor
             for ( auto fit=it.second.begin(); fit!=it.second.end(); ++fit )
               for ( size_t node=0U; node<(*fit)->Nodes(); ++node )
                 if ( it.first == (*fit)->N(node) )
                   (*fit)->Assign( node, static_cast<Face<3U>*>(nullptr) );
          */
          }
        // 2.2 when both nodes/faces have line-element neighbors
        else {
             Face<3U>* edge1 = (*it.second.begin());
             Face<3U>* edge2 = (*it.second.rbegin());
             // finding the node in the first Face = line element
             for ( size_t i=0U; i<edge1->Nodes(); ++i )
               if ( it.first == edge1->N(i) ) {
                    // assigning the opposite neighbor
                    edge1->Assign( i, edge2 );
                    break;
                 }
             // finding node number in second Face
             for ( size_t i=0U; i<edge2->Nodes(); ++i )
               if ( it.first == edge2->N(i) ) {
                    // assigning neighbors
                    edge2->Assign( i, edge1 );
                    break;
                 }
          }
     }
   
 } // create line element neighbor connectivity

bool createBoundaryFromSharedEdge( Model<1U>&, const Boundary<1U>&, const Boundary<1U>&, std::vector<Face<1U>*>& ) { throw logic_error("createBoundaryFromSharedEdge(1D)"); }
bool createBoundaryFromSharedEdge( Model<2U>&, const Boundary<2U>&, const Boundary<2U>&, std::vector<Face<2U>*>& ) { throw logic_error("createBoundaryFromSharedEdge(2D)"); }




/**
    @param shared_faces is a vector of the Face objects that were created in the MeshManager in the boundary creation process
 
    @return if the operation was successful
*/
bool createBoundaryFromSharedEdge( Model<3U>& model, const Boundary<3U>& boundary1, const Boundary<3U>& boundary2, std::vector<Face<3U>*>&  shared_faces )
 {
    if ( !shared_faces.empty() ) shared_faces.clear();
   
    // creating keys for the perimeter element faces of boundary1
    map<set<csmp::Node<3U>*>,Face<3U>*>  perimeter_keys1, perimeter_keys2;
    createPerimeterKeysFor( boundary1, perimeter_keys1 );
    createPerimeterKeysFor( boundary2, perimeter_keys2 );

    // LVS
    const LocalVariables lvsFaces( model.Database().LocalVariablesAt(FACE) );
    const IntegrationPointVariables lvsIntegrationPoints( model.Database().IntegrationPointVariablesAt(FACE) );

    shared_faces.reserve( perimeter_keys1.size() );	
    for ( typename map<set<csmp::Node<3U>*>,Face<3U>*>::const_iterator
          it1=perimeter_keys1.begin(); it1!=perimeter_keys1.end(); ++it1 )
      {
         // finding pairs of matching faces along the boundary
         typename map<set<csmp::Node<3U>*>,Face<3U>*>::const_iterator it2(perimeter_keys2.find( (*it1).first ));
         // if there is a matching face pair of adjacent Faces is captured
         // and their parent elements are discovered
         if ( it2 != perimeter_keys2.end() ) {
               vector<Node<3U>*>  segment_nodes;
               // finding the parent elements of the neighbors faces as these will
               // be used to create the higher dimensional neighbors of the edge Faces ?!
               // convention: the element in the first map is the inner element
               //             it is taken as the edge created from boundary which precedes the other in the enumeration
               Element<3U>*   inner   = (*it1).second->InnerParent();
               Element<3U>*   outer   = (*it2).second->InnerParent();
               // the same finite element type represents all edges in volumetric elements
               FiniteElement* fem_ptr = model.FE_Manager().E( inner->FE()->ElementTypeOfSegment(0) );
               assert( fem_ptr != nullptr );
               // establishing the face-node sequence of the inner element face that will be shared with the new Face object
               for ( size_t i=0U; i<inner->Segments(); ++i ) {
                    vector<size_t> snids; // local segment node ids
                    inner->FE()->NodesOfSegment( i, snids );
                    set<csmp::Node<3U>*> nset;
                    for ( size_t k=0U; k<snids.size(); ++k ) nset.insert( inner->N(snids[k]) );
                    if ( nset == (*it1).first ) {
                          // establishing the finite-element type associated with Segment (of
                          // fem_ptr = model.FE_Manager().E( inner->FE()->ElementTypeOfSegment(i) );
                          // assert( fem_ptr != nullptr );
                          // capturing the segment nodes for the construction of the Face object
                          segment_nodes.reserve( snids.size() );
                          for ( size_t j=0U; j<snids.size(); ++j )
                            segment_nodes.push_back( inner->N(snids[j]) );
                          assert( segment_nodes.size() >= 2U );
                          break;
                      }
                 }
              // creating the faces
			  Face<3U> new_face(fem_ptr, inner, outer, segment_nodes, lvsFaces, lvsIntegrationPoints);
			  Face<3U>* faceObj = model.Mesh().Add(new_face);
			  
			  // the first face is assigned into the root face of this face group in the mesh
			  // the faces should be connected each other, otherwise each root face has only a single face
			  if (shared_faces.size()==1) {
				  model.Mesh().SetRootFace(shared_faces.front());
			  }
			  
			  shared_faces.push_back(faceObj);
           }
      }
   
    // assigning the equidimensional neighbors to the newly created faces
    createLineFaceConnectivity( shared_faces );
		
    // false if no shared faces could be detected
    return ( !shared_faces.empty() );

 } // end createBoundaryFromSharedEdge



  /// High Level Functions to create Boundaries from provided Region names


/** creates edge Boundary objects (of dim-2 Face objects) for box-shaped model from side boundaries
   
   1. verifies that that the side boundaries of the box-shaped model are there
   
   @attention EXCEPTION: the edge curves are represented by Face objects, however these 
   are dim-2 entities. Following CSMP's rules, these edges will be connected to their face neighbors
   AND their higher-dimensional parent elements with which they share edges.

*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel()
 {
    ErrorHandler& csmp_error(ErrorHandler::Instance());

    if ( dim != 3 ) {
         cout <<"\nBoundaryInterface<"<< dim <<">::EstablishEdgeBoundariesOfBoxShapedModel: edges only are required in 3D; nothing was done.\n";
         return false;
      }
    // getting references to all necessary side boundaries of the box-shaped model
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    csmp::Boundary<dim>&  bottom(boundaryComplex->Boundary("BOTTOM"));
    csmp::Boundary<dim>&  right(boundaryComplex->Boundary("RIGHT"));
    csmp::Boundary<dim>&  top(boundaryComplex->Boundary("TOP"));
    csmp::Boundary<dim>&  left(boundaryComplex->Boundary("LEFT"));
    csmp::Boundary<dim>&  back(boundaryComplex->Boundary("BACK"));
    csmp::Boundary<dim>&  front(boundaryComplex->Boundary("FRONT"));

    // edge Faces that were created and are needed to create edge boundary
    vector<Face<dim>*>  shared_faces;
    bool                return_value(true);
   
    // EDGE1 = BACK_BOTTOM
    // -------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, bottom, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE1 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE1", shared_faces.begin(), shared_faces.end(), EDGE1 );
    if ( shared_faces.empty() ) return_value=false;
 
    // EDGE2 = BACK_RIGHT
    // ------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, right, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE2 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE2", shared_faces.begin(), shared_faces.end(), EDGE2 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE3 = BACK_TOP
    // ----------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, top, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE3 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE3", shared_faces.begin(), shared_faces.end(), EDGE3 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE4 = BACK_LEFT
    // -----------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, left, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE4 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE4", shared_faces.begin(), shared_faces.end(), EDGE4 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE5 = BOTTOM_LEFT
    // -------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, bottom, left, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE5 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE5", shared_faces.begin(), shared_faces.end(), EDGE5 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE6 = BOTTOM_RIGHT
    // --------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, bottom, right, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE6 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE6", shared_faces.begin(), shared_faces.end(), EDGE6 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE7 = TOP_RIGHT
    // -----------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, top, right, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE7 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE7", shared_faces.begin(), shared_faces.end(), EDGE7 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE8 = TOP_LEFT
    // ----------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, top, left, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE8 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE8", shared_faces.begin(), shared_faces.end(), EDGE8 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE9 = FRONT_BOTTOM
    // --------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, bottom, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE9 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE9", shared_faces.begin(), shared_faces.end(), EDGE9 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE10 = FRONT_RIGHT
    // --------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, right, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE10 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE10", shared_faces.begin(), shared_faces.end(), EDGE10 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE11 = FRONT_TOP
    // ------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, top, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE11 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE11", shared_faces.begin(), shared_faces.end(), EDGE11 );
    if ( shared_faces.empty() ) return_value=false;

    // EDGE12 = FRONT_LEFT
    // -------------------
    if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, left, shared_faces ) )
      csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                        "Boundary EDGE12 could not be created because no shared edge was found." );

    else AddBoundary( "EDGE12", shared_faces.begin(), shared_faces.end(), EDGE12 );
    if ( shared_faces.empty() ) return_value=false;
   
    return return_value; // whether all edges could be established
   
 } // end EstablishEdgeBoundariesOfBoxShapedModel





  /**
       Searching the VSet for line elements identified by their id (0..n-1) 
       to form Region objects of line elements making up the EDGES of the box shaped model.
       These edge elements are then used to flag the boundary nodes of the
       model.
       
       @attention for this method to work, line elements must be present on the edges of the input model.
       
       @attention the presence of such line elements changes the properties of the model
       if they are kept in the subdomain "Model"
       
       @attention unless specific line-element regions were created by the model builder,
       2D  (rectangular) models have no edges, but only BOX_BOUNDARY-flagged boundary points;
       their sides are called boundaries.

       @attention model corners are not considered.

  */
  template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
  bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeRegionsOfBoxShapedModel( const VSet<dim>& vset )
    {
      BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );

      std::cout << "\nBoundaryInterface<"<< dim <<">::EstablishEdgeRegionsOfBoxShapedModel: Establishing Edge Regions for " << dim << " dimensional box shaped model...";
      ErrorHandler& csmp_error(ErrorHandler::Instance());

      // Forming edge regions
      size_t number_of_edges_to_be_found( 12 );
      if( dim == 2 )
          number_of_edges_to_be_found = 4;

      size_t number_of_edges( 12 );

      std::map<size_t,int32>          bflags(vset.BFlagsBegin(),vset.BFlagsEnd());
      std::vector<bool>               edge( number_of_edges,false);
      std::vector<BOX_BOUNDARY>       edge_flag(number_of_edges,NOT);
      std::vector<std::set<size_t> >  edge_elmt_idx(number_of_edges);

      edge_flag[0]  = EDGE1;
      edge_flag[1]  = EDGE2;
      edge_flag[2]  = EDGE3;
      edge_flag[3]  = EDGE4;
      edge_flag[4]  = EDGE5;
      edge_flag[5]  = EDGE6;
      edge_flag[6]  = EDGE7;
      edge_flag[7]  = EDGE8;
      edge_flag[8]  = EDGE9;
      edge_flag[9]  = EDGE10;
      edge_flag[10] = EDGE11;
      edge_flag[11] = EDGE12;

      // 1. collecting the ids of the line elements that sit on the model edges also remembering their flags
      std::cout <<"\nBoundaryInterface<"<< dim <<">::EstablishEdgeRegionsOfBoxShapedModel: Assigning line elements to model edges..." << std::endl;
      // (only if there are multiple element types in the VSet)
      if ( vset.ElementTypes() > 1 )
        for ( size_t e=0U; e<vset.Elements(); ++e )
          // if we are dealing with a line element
          if ( vset.ElementType(e) == ISOPARAMETRIC_LINEAR_BAR || vset.ElementType(e) == ISOPARAMETRIC_QUADRATIC_BAR )
            {
                size_t edge_idx( 0 );
                bool   belongs_to_edge( false );

                // check nodes to see which boundary they belong to
                for ( size_t j=0U; j < vset.PlistSize( e ); ++j )
                  {
                      // get boundary flag
                      std::map<size_t,int32>::const_iterator it = bflags.find( vset.Plist( e, j ) );
                      if ( it == bflags.end() ) break;

                      // checking whether we are on an edge
                      const BOX_BOUNDARY b = intToBOX_BOUNDARY( (*it).second );

                      if ( !belongs_to_edge )
                          for( size_t i=0; i<number_of_edges; i++ )
                          {
                              if( belongsToEdge( edge_flag[ i ],  b ) )
                              {
                                  belongs_to_edge = true;
                                  edge_idx = i;
                                  break;
                              }
                          }
                      else if( !belongsToEdge( edge_flag[ edge_idx ],  b ) )
                          belongs_to_edge = false;

                      if ( !belongs_to_edge ) break;
                  }

                if ( belongs_to_edge ) {
                     edge[ edge_idx ] = true;
                     edge_elmt_idx[ edge_idx ].insert( e );
                     break;
                  }
            }

      // 2. Forming edge regions from the collected line elements
      std::vector<size_t>  elmt_indexes;
      for( size_t i=0; i<number_of_edges; i++ )
          if( edge[i] )
          {
              // collecting the element ids into a vector from which the region will be formed
              elmt_indexes.assign( edge_elmt_idx[i].begin(), edge_elmt_idx[i].end() );
              // forming region with the edge name
              boundaryComplex->FormRegionFrom( parseBoundary( edge_flag[ i ] ).c_str(), elmt_indexes );
              elmt_indexes.clear();
          }

      // 3. checking that all edges were found
      size_t number_of_found_edges(0);
      for( size_t i = 0 ; i< number_of_edges_to_be_found; i++ )
          if( edge[ i ] )
              number_of_found_edges++;
      bool all_edges_were_found( number_of_edges_to_be_found == number_of_found_edges );
    
      // SKM: in 2D models I do not expect to see edges
      if ( !all_edges_were_found and dim == 3 )
        {
            for( size_t i=0; i<number_of_edges_to_be_found; i++ )
               std::cerr<<"\n EDGE[ "<<i+1<<" ], Flag = "<<parseBoundary( edge_flag[ i ] )<<":\t"<<( edge[i] ? "exist":"does not exist.");
            std::cerr<<std::endl;
            csmp_error.notice( INFO, "BoundaryInterface::Initialize (box-shaped model):",
                              "Not all of the edges could be formed");
        }
       std::cout << "\n...done!\n";

   return true;
      
 } // end EstablishEdgeRegionsOfBoxShapedModel
 
 



/**
     Forms boundaries of CSMP box-shaped model if corresponding regions are present.
     For 3D models, method also creates edge regions where the side boundaries intersect.
     These regions are given the standard names and are flagged correspondingly.  
     
     @note The boundary creation itself does not deal with the generation of BOX_BOUNDARY flags for the model
     corners. This is accomplished subsequently (in this method) by calling recreateBoxBoundaryFlags().
     
     @author refactored by SKM 2016
     @author refactored by SKM 2018
*/
  template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
  bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishBoxBoundaries()
   {
      BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );	  
      std::cout << "\nBoundaryInterface<"<< dim <<">::EstablishBoxBoundaries: Establishing Box-object boundaries for " << dim << " dimensional box shaped model...";

	  // TOP may be missing, if there is an IRREGULAR boundary instead
	  if (boundaryComplex->ContainsRegion("TOP"))
		  boundaryComplex->InsertBoundary(TOP, "TOP");
      boundaryComplex->InsertBoundary( BOTTOM, "BOTTOM" );
      boundaryComplex->InsertBoundary( RIGHT, "RIGHT" );
      boundaryComplex->InsertBoundary( LEFT, "LEFT" );
      if( dim == 3 ) {
           boundaryComplex->InsertBoundary( FRONT, "FRONT" );
           boundaryComplex->InsertBoundary( BACK, "BACK" );
        }
      // potential irregular model outside boundaries, like for instance in a box with topography on top
      if ( boundaryComplex->ContainsRegion("IRREGULAR") )
        boundaryComplex->InsertBoundary( IRREGULAR, "IRREGULAR" );
     
      // creating the edges needed in a three-dimensional model
      if( dim == 3 ) EstablishEdgeBoundariesOfBoxShapedModel();

      std::set<std::string> regionsToRemove;
      for ( typename std::map<std::string,csmp::Region<dim> >::iterator
          it = boundaryComplex->UniqueRegionsBegin(); it != boundaryComplex->UniqueRegionsEnd(); ++it )
        if ( isSide( parseBoundary((*it).first) ) || isEdge( parseBoundary((*it).first) ) )
          regionsToRemove.insert(it->first);

     // doing the removal of all box boundaries in one go
      set<csmp::Element<dim>*>  elmts_to_remove;
      if ( !regionsToRemove.empty() )
        std::cout << "\n\nBoundaryInterface<"<< dim <<">::EstablishBoxBoundaries: Removing regions since they were transformed to Boundaries:\n\t";
      for ( auto rit=regionsToRemove.begin(); rit!=regionsToRemove.end(); ++rit )
        if ( (*rit) != "Model" ) {
             std::cout << (*rit) << " ";
             const Region<dim>& subdomain(boundaryComplex->Region(*rit));
             for ( auto eit=subdomain.ElementsBegin(); eit!=subdomain.ElementsEnd(); ++eit )
               elmts_to_remove.insert(*eit);
          }
      std::cout << std::endl;
     
      // the regions that were converted into boundaries are removed from the region 'Model'	  
      boundaryComplex->RemoveFromRegion( "Model", elmts_to_remove );
      // then the regions and their elements are removed
      const bool remove_elements(true);
	  for (auto rit = regionsToRemove.begin(); rit != regionsToRemove.end(); ++rit)
	  {
		  boundaryComplex->RemoveRegion((*rit).c_str(), remove_elements);
		  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoxBoundaries: Removing the region " << (*rit).c_str();
	  }
	        
      std::cout << "\n\n EstablishBoxBoundaries: done!\n";
      return true;
  }




  /**
  Creates boundary around "Model" region determining whether lower-dimensional
  regions shall become a boundary segments on the basis of whether their name contains
  the string "BOUNDARY" or is one of the reserved boundary names.

  When a boundary (of Face objects) is created, the underpinning lower-dimensional (parent)
  region is removed from the region "Model".
  It is also moved from the unique to the non-unique region map if it was
  stored there originally.

  @attention by contrast to EstablishBoxBoundaries(),
  EstablishBoundaries() does not require to be applied to box-shaped models, but is
  more general.

  @test updated by SKM 2016
  */
  template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
  bool BoundaryInterface<dim, BOUNDARY_COMPLEX>::EstablishBoundariesFromRegions(bool remove_original_lower_dimensional_regions)
  {

	  BOUNDARY_COMPLEX<dim>* boundaryComplex(static_cast<BOUNDARY_COMPLEX<dim>*>(this));
	  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundaries: searching for eligible boundary domains...\n";

	  // generate hull of Face objects around model calling it "Model_BOUNDARY"
	  AddFaces("Model");
	  csmp::Boundary<dim>& modelBoundary(Boundary(std::string("Model_BOUNDARY")));
	  boundaryIterator     modelBoundaryIt(Boundary(modelBoundary));

	  std::vector<std::string> eligibleRegions;
	  eligibleRegions.reserve(boundaryComplex->UniqueRegions());

	  for (typename std::map<std::string, csmp::Region<dim> >::iterator
		  it = boundaryComplex->UniqueRegionsBegin(); it != boundaryComplex->UniqueRegionsEnd(); ++it)
	  {
		  if (it->first == "Model_BOUNDARY")
			  continue;
		  if (!IsBoundaryName(it->first))
			  continue;
		  if (!isOfLowerDimensionalRepresentation(it->second))
			  continue;
		  if (dim == 3U and containsVolumeElements(boundaryComplex->Region("Model")) and !containsSurfaceElements(it->second))
			  continue;
		  eligibleRegions.push_back(it->first);
		  DivideBoundary(modelBoundaryIt, it);
		  if (modelBoundary.Elements() == 0)
		  {
			  boundaryComplex->RemoveBoundary(modelBoundary);
			  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundariesFromRegions: Boundary 'Model_BOUNDARY' entirely replaced by sub boundaries.\n";
			  break;
		  }
	  }

	  // if the boundary called 'Model_BOUNDARY' that was created by AddFaces() is removed, i.e. it is a leftover that is no-longer needed
	  if (count(eligibleRegions.begin(), eligibleRegions.end(), "Model_BOUNDARY") == 0) {
		  faceBoundaryMap_.erase("Model_BOUNDARY");		  
	  }

	  // moving the original regions from which the boundaries were created from model and into the non-unique regions map
	  for (size_t i = 0; i < eligibleRegions.size(); ++i) {
		  auto& region_name = eligibleRegions[i];
		  if (region_name != "Model_BOUNDARY") {
			  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundariesFromRegions: Removing region '";
			  std::cout << region_name << "' from 'Model' since it was transformed into Boundary...";
			  boundaryComplex->RemoveFromRegion("Model", region_name.c_str());
			  if (remove_original_lower_dimensional_regions) {
				  boundaryComplex->RemoveRegion(region_name.c_str(), true);
			  }
			  else {
				  boundaryComplex->MoveToNonUniqueRegions(region_name.c_str());
			  }
		  }
	  }

	  // changing all BOX_BOUNDARY flags on the outside of the model to IRREGULAR, unless a box-boundary name is recognised
	  // (edges are not considered)
	  for (auto it = boundaryComplex->BoundariesBegin(); it != boundaryComplex->BoundariesEnd(); ++it) {
		  BOX_BOUNDARY bflag(IRREGULAR);
		  if (isDiagnosticBoxBoundaryClassifier((*it).first)) bflag = parseBoundary((*it).first);
		  for (auto nit = (*it).second.NodesBegin(); nit != (*it).second.NodesEnd(); ++nit)
			  (*nit)->AtBoundary(bflag);
	  }

	  std::cout << "\n\nBoundaryInterface::EstablishBoundariesFromRegions: done!\n";
	  return true;

  } // end EstablishBoundariesFromRegions

	/**
	Creates boundaries around discontiguous model regions determining whether lower-dimensional
	regions shall become a boundary segments on the basis of whether their name contains
	the string "BOUNDARY" or is one of the reserved boundary names.

	When a boundary (of Face objects) is created, the underpinning lower-dimensional (parent)
	region is removed from the region "Model".
	It is also moved from the unique to the non-unique region map if it was
	stored there originally.

	@attention by contrast to EstablishBoxBoundaries(),
	EstablishBoundariesFromDiscontiguousModel() does not require to be applied to box-shaped models, but is
	more general, and inserts irregular csmp::Boundary for all eligible regions in the discontiguous model
	*/
  template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
  bool BoundaryInterface<dim, BOUNDARY_COMPLEX>::EstablishBoundariesFromDiscontiguousModel(bool remove_original_lower_dimensional_regions)
  {

	  BOUNDARY_COMPLEX<dim>* boundaryComplex(static_cast<BOUNDARY_COMPLEX<dim>*>(this));
	  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundariesFromDiscontiguousModel: searching for eligible boundary domains...\n";

	  // generate hull of Face objects around model calling it "Model_#n_BOUNDARY"
	  // ie. if there are 2 subdomains, one is "Model_0" and second is "Model_1".
	  //     their boundary names are "Model_0_BOUNDARY" and "Model_1_BOUNDARY"
	  auto& meshMgr = boundaryComplex->Mesh();

	  std::vector<std::string> subdomains;	  
	  std::string model_name = "Model";	  
	  for (size_t i = 0; i < meshMgr.NodeGroups(); i++) {
		  if (i > 0) model_name = "Model_" + std::to_string(i); // first model name is 'Model', ann then Model_1, Model_2 and so on.
		  subdomains.push_back(model_name);
	  }

	  bool ret = true;
	  for (auto model_name : subdomains)	  
		  AddFaces(model_name.c_str());
	  
	  std::set<std::string> eligibleRegions;
	  for (typename std::map<std::string, csmp::Region<dim> >::iterator
		  it = boundaryComplex->UniqueRegionsBegin(); it != boundaryComplex->UniqueRegionsEnd(); ++it)
	  {
		  for (auto model_name : subdomains){
			  std::string boundary_name = model_name + "_BOUNDARY";

			  csmp::Boundary<dim>& modelBoundary(Boundary(boundary_name));
			  boundaryIterator     modelBoundaryIt(Boundary(modelBoundary));

			  if (it->first == boundary_name)
				  continue;
			  if (!IsBoundaryName(it->first))
				  continue;
			  if (!isOfLowerDimensionalRepresentation(it->second))
				  continue;
			  if (dim == 3U and containsVolumeElements(boundaryComplex->Region(model_name.c_str())) and !containsSurfaceElements(it->second))
				  continue;
			  eligibleRegions.insert(it->first);
			  DivideBoundary(modelBoundaryIt, it);
			  if (modelBoundary.Elements() == 0)
			  {
				  boundaryComplex->RemoveBoundary(modelBoundary);
				  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundariesFromDiscontiguousModel: Boundary '" << boundary_name << "' entirely replaced by sub boundaries.\n";
				  continue;
			  }
		  }
	  }

	  // if the boundary called 'Model_#n_BOUNDARY' that was created by AddFaces() is removed, i.e. it is a leftover that is no-longer needed
	  for (auto model_name : subdomains) {
		  std::string boundary_name = model_name + "_BOUNDARY";
		  if (count(eligibleRegions.begin(), eligibleRegions.end(), boundary_name) == 0) {
			  faceBoundaryMap_.erase(boundary_name);			
		  }
	  }

	  // moving the original regions from which the boundaries were created from model and into the non-unique regions map
	  for (auto region_name : eligibleRegions) {		  
		  for (size_t i = 0; i < meshMgr.NodeGroups(); i++) {
			  std::string model_name = "Model_" + std::to_string(i);
			  if (i == 0) model_name = "Model"; // first model name is 'Model', ann then Model_1, Model_2 and so on.
			  std::string boundary_name = model_name + "_BOUNDARY";
			  if (region_name != boundary_name) {			  
				  std::cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundariesFromDiscontiguousModel: Removing region '";
				  std::cout << region_name << "' from '"<< model_name << "' since it was transformed into Boundary...";
				  boundaryComplex->RemoveFromRegion(model_name.c_str(), region_name.c_str());
			  }
		  }
	  }

	  for (auto region_name : eligibleRegions) {
		  if (remove_original_lower_dimensional_regions) {
			  boundaryComplex->RemoveRegion(region_name.c_str(), true);
		  }
		  else {
			  boundaryComplex->MoveToNonUniqueRegions(region_name.c_str());
		  }
	  }
  
	  // changing all BOX_BOUNDARY flags on the outside of the model to IRREGULAR, unless a box-boundary name is recognised
	  // (edges are not considered)
	  for (auto it = boundaryComplex->BoundariesBegin(); it != boundaryComplex->BoundariesEnd(); ++it) {
		  BOX_BOUNDARY bflag(IRREGULAR);
		  if (isDiagnosticBoxBoundaryClassifier((*it).first)) bflag = parseBoundary((*it).first);
		  for (auto nit = (*it).second.NodesBegin(); nit != (*it).second.NodesEnd(); ++nit)
			  (*nit)->AtBoundary(bflag);
	  }

	  if (ret)
		  std::cout << "\n\nBoundaryInterface::EstablishBoundariesFromDiscontiguousModel: established successfully!\n";

	  return ret;

  } // end EstablishBoundariesFromRegions
	
  


// little utility
template<size_t dim> void setToVector( const set<Face<dim>*>& input, vector<Face<dim>*>& output )
  {
     if ( !output.empty() ) output.clear();
     output.reserve( input.size() );
     for ( auto sit : input )
       output.push_back( sit );
  }

/**
    Tries to partition and replace general boundary 'Model' with more computationally useful model patches
    such as TOP, BOTTOM, INTERNAL, IRREGULAR, VERTICAL_SIDE etc.
    
    @attention this method expects that there is already a Boundary 'Model' around the model domain.
    @attention method assumes that model perimeter correctly captures the outside faces of the model.
    @attention this method was designed primarily for three-dimensional models.
    
    @author SKM
    @date 21/8/2018
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishRegularities()
  {

    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    std::cout << "\nBoundaryInterface<"<< dim <<">::EstablishRegularities: searching for eligible boundary domains...\n";

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( !boundaryComplex->ContainsBoundary("Model") ) {
         csmp_error.notice( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::AddBoundary:",
                            "Model", "Boundary was missing. Creating it now.");

         AddFaces("Model");
      }

    csmp::Boundary<dim>& modelBoundary( Boundary( std::string("Model_BOUNDARY") ) );

    std::vector<std::string> eligibleRegions;
    eligibleRegions.reserve( boundaryComplex->UniqueRegions() );
    
    // 1. Grouping pointers to Face objects of 'Model' boundary according to their facing direction
    // -------------------------------------------------------------------------------------------------
    vector<double64>      nrml, nrml_right, nrml_left, nrml_top, nrml_bottom, nrml_front, nrml_back;
    Box                   box;
    double64              minLength(0.71); // dot-product of 2 unit vectors at an angle >=45 degrees
   
    box.UnitNormalTo( BOTTOM, dim, nrml_bottom );
    box.UnitNormalTo( TOP,    dim, nrml_top );
    box.UnitNormalTo( LEFT,   dim, nrml_left );
    box.UnitNormalTo( RIGHT,  dim, nrml_right );
    box.UnitNormalTo( FRONT,  dim, nrml_front );
    box.UnitNormalTo( BACK,   dim, nrml_back );
    
    set<Face<dim>*>  top_faces, bottom_faces, left_faces, right_faces, front_faces, back_faces, irregular_faces;
    
    // for all Face objects on the model boundary
    for ( auto fit=modelBoundary.ElementsBegin(); fit!=modelBoundary.ElementsEnd(); ++fit )
      {
         // getting the (outward) pointing unit normal to face
         (*fit)->UnitNormal( nrml );
         // classifying the faces in terms of their facing direction
         // (projecting: perfect alignment would give dot-product equal 1, inclinations up to 37 degrees cos(37)~0.8 are tolerated)
         if      ( dotProduct<dim>(nrml,nrml_bottom) >= minLength ) bottom_faces.insert(*fit); // BOTTOM
         else if ( dotProduct<dim>(nrml,nrml_top)    >= minLength ) top_faces.insert(*fit);    // TOP
         else if ( dotProduct<dim>(nrml,nrml_left)   >= minLength ) left_faces.insert(*fit);   // LEFT
         else if ( dotProduct<dim>(nrml,nrml_right)  >= minLength ) right_faces.insert(*fit);  // RIGHT
         else if ( dotProduct<dim>(nrml,nrml_front)  >= minLength ) front_faces.insert(*fit);  // FRONT
         else if ( dotProduct<dim>(nrml,nrml_back)   >= minLength ) back_faces.insert(*fit);   // BACK
         // deal with the remaining cases, distinguishing sides etc.
         else {
              irregular_faces.insert(*fit);
           }
      } // end perimeter faces
    
    // 2. Creating boundaries from the non-empty sets of faces
    // -------------------------------------------------------------------------------------------------
    vector<Face<dim>*> boundary_faces;
    // BOTTOM
    if ( !bottom_faces.empty() ) {
         setToVector( bottom_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "BOTTOM" );
      } // TOP
    if ( !top_faces.empty() ) {
         setToVector( top_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "TOP" );
      } // LEFT
    if ( !left_faces.empty() ) {
         setToVector( left_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "LEFT" );
      } // RIGHT
    if ( !right_faces.empty() ) {
         setToVector( right_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "RIGHT" );
      } // FRONT
    if ( !front_faces.empty() ) {
         setToVector( front_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "FRONT" );
      } // BACK
    if ( !back_faces.empty() ) {
         setToVector( back_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "BACK" );
      } // IRREGULAR
    if ( !irregular_faces.empty() ) {
         setToVector( irregular_faces, boundary_faces );
         boundaryComplex->InsertBoundary( boundary_faces.begin(), boundary_faces.end(), "IRREGULAR" );
      }

    
    // 4. changing BOX_BOUNDARY flags on the outside of the model so that TOP and BOTTOM are recognised; else they were set to irregular
    // (edges are not considered)
    // -------------------------------------------------------------------------------------------------
    for ( auto fit=modelBoundary.NodesBegin(); fit!=modelBoundary.NodesEnd(); ++fit )
      if ( (*fit)->AtBoundary() == NOT )
        (*fit)->AtBoundary( IRREGULAR );
 
    // the boundary called 'Model' that was created by AddFaces() is removed, i.e. it is a leftover that is no-longer needed
    faceBoundaryMap_.erase("Model_BOUNDARY");

    // for remaining (new) TOP and BOTTOM boundaries
  	for ( auto it = boundaryComplex->BoundariesBegin(); it != boundaryComplex->BoundariesEnd(); ++it )
      if ( (*it).first == "TOP" or (*it).first == "BOTTOM" )
        {
           BOX_BOUNDARY bflag(IRREGULAR);
           if ( isDiagnosticBoxBoundaryClassifier( (*it).first ) ) bflag = parseBoundary( (*it).first );
           for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.NodesEnd(); ++nit )
             (*nit)->AtBoundary( bflag );
        }

    std::cout << "\n\nBoundaryInterface::EstablishRegularities: done!\n";
    return true;
    
} // end EstablishRegularities










  template class BoundaryInterface<1U, Model>;
  template class BoundaryInterface<2U, Model>;
  template class BoundaryInterface<3U, Model>;

} // csmp

