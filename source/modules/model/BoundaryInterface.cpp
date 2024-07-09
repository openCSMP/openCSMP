#include "BoundaryInterface.h"
#include "ModelTopology.h"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"
#include "MeshManagementUtilities.h"
#include "FaceConstructionData.h"
#include "Element.h"
#include "Face.h"
#include "Boundary.h"
#include "Region.h"
#include "Box.h"
#include "VSet.h"
#include "Exception.h"
#include "SmallSet.h"
#include "ErrorHandler.h"
#include "variableOperations.h"
#include "binaryReadWrite.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
typename BoundaryInterface<dim,BOUNDARY_COMPLEX>::boundaryIterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesBegin()
  { return boundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
typename BoundaryInterface<dim,BOUNDARY_COMPLEX>::boundaryIterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesEnd()
  { return boundaryMap_.end(); }


template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
typename BoundaryInterface<dim,BOUNDARY_COMPLEX>::boundaryConstIterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesBegin() const
  { return boundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
typename BoundaryInterface<dim,BOUNDARY_COMPLEX>::boundaryConstIterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::BoundariesEnd() const
  { return boundaryMap_.end(); }


template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
typename BoundaryInterface<dim,BOUNDARY_COMPLEX>::boundaryIterator  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary( const csmp::Boundary<dim>& bref )
  { 
    for( boundaryIterator it( BoundariesBegin() ); it != BoundariesEnd(); ++it )
      if( &(it->second) == &bref )
        return it;
    return BoundariesEnd();
  }

template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
size_t  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundaries() const
  { return boundaryMap_.size(); }



template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool  BoundaryInterface<dim,BOUNDARY_COMPLEX>::ContainsBoundary( const string& bname ) const
 {
    if ( boundaryMap_.find(bname) != boundaryMap_.end()  ) return true;
    return false;
 }



/// returns a reference to the boundary 'bname' if it exists
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
Boundary<dim>&  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary( const string& bname )
  { 
    boundaryIterator  bit=boundaryMap_.find( bname );
    if( bit != boundaryMap_.end() )
      return bit->second;  
    else
      {
        string errMsg("Boundary does not exist!");
        errMsg.append( " (" + bname  + ")" );
        throw csmp::Exception( ERROR,
                               "BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary(string&)",
                               errMsg.c_str() );
      }
  } 



template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
const Boundary<dim>&  BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary( const string& bname ) const
  {    
    boundaryConstIterator  bit=boundaryMap_.find( bname );
    if( bit != boundaryMap_.end() )
      return bit->second;  
    else
      {
        string errMsg("Boundary does not exist!");
        errMsg.append( " (" + bname + ")" );
        throw csmp::Exception( ERROR,
                               "BoundaryInterface<dim,BOUNDARY_COMPLEX>::Boundary(string&)",
                               errMsg.c_str() );
      }
  } 




/**
    Creates a name for a new boundary that shall be created between two touching unique (non-overlapping) regions
    
    Name will be 'BOUNDARY' + inner + "_" + outer contacting regions. If the name already exists, a number will be appended to boundary
    to make it unique.
    
    @author SKM
    @date 17/4/22
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
string  BoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateBoundaryName( const std::string& inside_region, const std::string& outside_region )
 {
    string boundary_name = "BOUNDARY_" + inside_region + "_" + outside_region;
    
    int boundary_count{1};
    while( boundaryMap_.count( boundary_name ) > 0 ) {
         boundary_name = "BOUNDARY" + to_string( boundary_count );
         boundary_name +="_" + inside_region + "_" + outside_region;
         boundary_count++;
      }
      
    return boundary_name;
    
 } // end CreateBoundaryName




/**
     Returns name of the first Boundary (patch) that contains the supplied
     region name strings or any other other strings in an arbitrary order.
     
     if the boundary does not exist, method returns \O string.
     
     @author SKM 
     @date March 2016
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
string  BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryByName( const set<string>& intersected_regions ) const
 {
    const BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this) );
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryByName:",
                                         "supplied set of substrings is empty; returning '\0'." );
         return string("\0");
      }
    // if the model has no boundaries
    if ( boundaryComplex.Boundaries() == 0 ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryByName:",
                                         "model has no boundaries; returning '\0'." );
         return string("\0");
      }
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=boundaryComplex.BoundariesBegin(); it!=boundaryComplex.BoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
            // if the substring is found
            if ( (*it).first.find(*ir) !=string::npos ) substrings_found++;
            // when all substrings are contained in the boundary name, it is returned
          if ( substrings_found == substrings_used_in_search )
            return (*it).first;
       }
    return string("\0");
 }



/**
    Searches the model for boundaries the name of which contains the search strings
    provided via the first set. The results are returnd into the second set.
    
    @note Use this method, for example, to retrieve multiple boundary patches that were generated
    from a single lower-dimensional regon, like a fault surface.

     @author SKM 
     @date March 2016
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
size_t BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryByNames( const set<string>& intersected_regions,
                                                                     set<string>& region_patches_found ) const
 {
    const BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this) );
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryByNames:",
                                         "supplied set of substrings is empty; returning '\0'." );
         return 0U;
      }
    // if the model has no boundaries
    if ( boundaryComplex.Boundaries() == 0 ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::FindBoundaryByNames:",
                                         "model has no boundaries; returning '\0'." );
         return 0U;
      }
    region_patches_found.clear();
   
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=boundaryComplex.BoundariesBegin(); it!=boundaryComplex.BoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
            // if the substring is found
            if ( (*it).first.find(*ir) !=string::npos ) substrings_found++;
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
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::IsBoundaryName( const string& regionName ) const
  {
    // TAG: BOUNDARY should be at the beginning of the region name
    // DON'T confuse with SPLITBOUNDARY
    size_t found_position = regionName.find( "BOUNDARY" );
    if ( found_position == 0 && found_position!=string::npos )
        return true;
        
    if( regionName == "INTERNAL" || regionName == "IRREGULAR" )
      return true;
    if( isDiagnosticBoxBoundaryClassifier(regionName) )
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
     
     4. the name of the inner region, i.e. the region that the lower-dimensional element normals point away from
     
     5. the name of the outer region, i.e. that into which the normals point
     
     @attention  where the boundary just intersects a single layer (same material on either side), the layer name appears
     only once. The second instance is replaced by INTERSECTION.
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
string BoundaryInterface<dim, BOUNDARY_COMPLEX>::CreateBoundaryNameFrom( const FaceConstructionData<dim>& fdata,
                                                                         const vector<string>& region_names ) const
 {
     assert( fdata.ElementMaterial() < region_names.size() );
     string boundary_name( region_names[ fdata.ElementMaterial() ] );
     boundary_name += "_BOUNDARY";
     boundary_name += to_string(fdata.PatchNumber());
     boundary_name += '_';
     // inside/outside
     pair<long,long> materials{ fdata.Materials() };
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
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim, BOUNDARY_COMPLEX>::BoxShaped() const
{
	const BOUNDARY_COMPLEX<dim>& boundaryComplex(static_cast<const BOUNDARY_COMPLEX<dim>& >(*this));

//	cout << "\nBoundaryInterface<" << dim << ">::BoxShaped: checking whether the model contains the boundaries LEFT, RIGHT...";
	size_t boundaries(6);
	if constexpr (dim == 2) boundaries = 4;
	set<BOX_BOUNDARY> boundariesFound;
	BOX_BOUNDARY currentBoundary(NOT);
	const Region<dim>& region(boundaryComplex.Region("Model"));
	const auto nodesEnd(region.NodesEnd());
	for ( auto it = region.NodesBegin(); it != nodesEnd; ++it )
    {
      currentBoundary = (*it)->AtBoundary();
      // also covers the case where there is an IRREGULAR boundary instead of the TOP boundary
      if (currentBoundary != LEFT and currentBoundary != RIGHT and currentBoundary != BACK
        and currentBoundary != FRONT and currentBoundary != BOTTOM and currentBoundary != TOP
        and currentBoundary != IRREGULAR)
        continue;
      boundariesFound.insert(currentBoundary);
      if ( boundariesFound.size() >= boundaries )
        return true;
    }
	return false;
  
} // end BoxShaped




/**
  Assuming that the necessary Face objects already exist, this method creates Boundaries interpreting the Cell 'idx'
  stored in the model topology object as Face 'idx' ranging between n_elements..interfaces-1
 */
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
size_t BoundaryInterface<dim,BOUNDARY_COMPLEX>::FormBoundariesFrom( const ModelTopology& topo )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. getting the names of the regions
  list<string> boundaries;
  topo.OutputBoundaries( boundaries );

  // 2. assigning the regions to groups in the Model
  cout << "\nBoundaryInterface::FormBoundariesFrom: Forming the boundaries: ";

  size_t new_boundaries{0};
  for ( const auto& lit : boundaries )
    {
      string domain_name( lit );
      const auto boundary_flag = parseBoundary( domain_name );
      auto it = boundaryMap_.insert( make_pair( domain_name, csmp::Boundary<dim>( domain_name,
                                                static_cast<BOUNDARY_COMPLEX<dim>*>(this)->Database(), boundary_flag ) ) );
      // if the region was successfully inserted
      if ( it.second )
        {
          // making a list of the element numbers
          vector<size_t>  cell_ids;
          cell_ids.reserve( topo.CellsWithinDomain( lit.c_str() ) );
          copy( topo.CellsOfDomainBegin( lit.c_str() ),
                topo.CellsOfDomainEnd( lit.c_str() ),
                back_inserter( cell_ids ) );

          // retrieving the faces by their 'idx' numbers and assigning them to the Boundary
          (*it.first).second.AccumulateByNumber( static_cast<BOUNDARY_COMPLEX<dim>*>(this)->Mesh(), cell_ids );
          cell_ids.erase( cell_ids.begin(), cell_ids.end() );

          // removing the Boundary again if it contains no elements
          if ( (*it.first).second.Cells() == 0U ) {
              boundaryMap_.erase( it.first );
              csmp_error.Note( WARNING, "BoundaryInterface::FormBoundariesFrom",
                                 "Boundary could not be formed", lit.c_str() );
            }
          else {
               // reporting the name of the newly generated region
               cout << domain_name << " ";
               new_boundaries++;
            }
        }
      else
        throw csmp::Exception( ERROR, "BoundaryInterface::FormBoundariesFrom",
                               "Boundary could not be formed. Does this region already exist?", lit.c_str() );
    }
  cout << endl;

  return new_boundaries;

} // end FormBoundariesFrom



/**
    Creates boundary from Faces that already know their parent elements and are interconnected to one another!
    
    The connectivity between the Faces is (re)established.
    
    @note the supplied Face pointer vector is moved into boundary and will therefore not be 
    accessible anymore after this method has been called.
    
    @author SKM
    @date 1/4/2016
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::AddBoundary( const char* boundary_name,
                                                           typename vector<Face<dim>*>::iterator facesBegin,
                                                           typename vector<Face<dim>*>::iterator facesEnd,
                                                           BOX_BOUNDARY bflag )
 {
    BOUNDARY_COMPLEX<dim>* const boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>* const>(this) );
    assert( boundaryComplex != nullptr );

    // inserting boundary if it does not existing yet
    auto it = boundaryMap_.insert( make_pair( boundary_name, csmp::Boundary<dim>( boundary_name,
                                                                                  boundaryComplex->Database(),
                                                                                  facesBegin, facesEnd, bflag ) ) );
    if ( it.second ) {
         cout << "\nBoundaryInterface<"<< dim <<">::AddBoundary: successfully created boundary '";
         cout << boundary_name <<"' from input faces.";
      }
    else {
         ErrorHandler&  csmp_error( ErrorHandler::Instance() );
         csmp_error.Note( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::AddBoundary:",
                            boundary_name, "Boundary already exists or other problem arose. Nothing was done.");
         return false;
      }

     return true;
   
 } // end AddBoundary








/** 
     Converts lower-dimensional region into Boundarie(s) of faces, decomposed into patches. The MeshManager updates the connectivity.
     The original region is removed (including removal from "Model") to non-unique regions.
     
     Uses node-to-parent relationship to find the higher dimensional elements that will 
     share a face with the Face:  usng    higherDimensionalNeighbors()
     
     @return pair of the set oof boundary patches that were created from the region and boolean that tells whether the operation was completely successful.
     
     Algorithmic steps:

     1. Verify input lower-dimensional region object from which the boundary shall be created: 
        must be lower dimensional and must lie inside of model
 
     2. Determine the number of boundary segments (sub-boundaries) that the new boundary will consist of.
        The output of this step will be a map of FaceConstructionData in which the names of the 
        new boundary segments are the keys.
 
     3. In the MeshManager object, 
 
        3.1 create all required face objects
 
        3.2 connect them with one another (neighbors); Boundary::EstablishNeighborConnectivity( vector<Face<dim>*>& );
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
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
pair<set<string>,bool>   BoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateInternalBoundaryFrom( const char* dim_1_region )
 {
    BOUNDARY_COMPLEX<dim>& model( static_cast<BOUNDARY_COMPLEX<dim>&>(*this) );

    cout <<"\nBoundaryInterface<"<< dim <<",BOUNDARY_COMPLEX>::CreateInternalBoundaryFrom: forming boundary(ies) from region: '"<< dim_1_region <<"'...\n";
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    // 1. Verify input lower-dimensional region object from which the boundary shall be created: must be lower dimensional and must lie inside of model
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    const string creation_failed( string("CreateInternalBoundaryFrom(") +  dim_1_region +") failed");
    // does the parent region exist
    if ( model.ContainsRegion(dim_1_region) == false ) {
          ErrorHandler::Instance().Note( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region, "does not exist; nothing was done." );
          return make_pair( set<string>({creation_failed}), false );
      }
    // do such boundaries already exist ?
    const set<string> intersected_regions{dim_1_region};
    set<string> pre_existing_boundaries;
    if ( FindBoundaryByNames( intersected_regions, pre_existing_boundaries ) > 0 ) {
          string error_info;
          for ( auto it=pre_existing_boundaries.begin(); it!=pre_existing_boundaries.end(); ++it ) {
               error_info += (*it);
               error_info +=", ";
            }
          ErrorHandler::Instance().Note( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:",
                                           error_info.c_str(), "boundaries are already contained in this model." );
          return make_pair( set<string>({creation_failed}), false );
      }
    // does the model contain unique regions
    if ( model.UniqueRegions() < 1 ) {
          ErrorHandler::Instance().Note( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", "model contains no unique regions; cannot proceed." );
          return make_pair( set<string>({creation_failed}), false );
      }
    // verifying that we are indeed dealing with a region of surface or line elements only and that their normals all point into same direction
    Region<dim>&  subdomain(model.Region(dim_1_region));
    if ( checkNeighborNormalsForConsistentOrientation( subdomain ) == false ) {
         ErrorHandler::Instance().Note( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region,
                                         "region appears to have inconsistent surface-normal orientations; nothing was done." );
         return make_pair( set<string>({creation_failed}), false );
      }
    // checking that the region is not already an internal boundary
    size_t nodes_flagged_internal_boundary{0U}, nodes_flagged_external_boundary{0U};
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         if ( (*nit)->AtBoundary() == INTERNAL ) nodes_flagged_internal_boundary++;
         else if ( (*nit)->AtBoundary() != NOT ) nodes_flagged_external_boundary++;
      }
    if ( nodes_flagged_internal_boundary >= subdomain.Nodes() - nodes_flagged_external_boundary ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region,
                                                   "region may already be a boundary; nothing was done." );
         return make_pair( set<string>({creation_failed}), false );
      }
    
    // checking that the region is not located at the model boundary
    size_t boundary_elements{0U};
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
      for ( auto i{0U}; i<(*eit)->Neighbors(); ++i ) {
           const BOX_BOUNDARY bflag = (*eit)->AtBoundary(i);
           if ( bflag != NOT and bflag != INTERNAL and bflag != IRREGULAR ) boundary_elements++;
        }
    if ( boundary_elements == subdomain.Cells() ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region,
                                                   "region appears to lie at the model boundary; nothing was done." );
         return make_pair( set<string>({creation_failed}), false );
      }
      
    // creating region labels and tagging the regions with unique integer indentifiers
    // if "region identifier is already defined it is assumed that the unique regions have already been labeled correctly
    const string region_tag("region identifier");
    vector<string>  region_names;
    if ( !model.Database().IsDefined(region_tag.c_str()) ) {
         model.CreateProperty( region_tag.c_str(), "rid", "uint", SCALAR, ELEMENT );
      }
    const csmp::Index mtrl_key = model.Database().StorageKey(region_tag.c_str());
    // needs to be done everytime because the number of unique regions may have changed
    const size_t model_regions = model.CountAndLabelUniqueRegions( region_tag.c_str(), region_names );
    if ( model_regions == 1 )
      ErrorHandler::Instance().Note( INFO, "BoundaryInterface::CreateInternalBoundaryFrom:", region_tag.c_str(),
                                           "is single valued; so there is only one patch." );
    // assigning region names
    region_names.reserve( distance( model.UniqueRegionsBegin(),model.UniqueRegionsEnd()) );
    for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit )
      region_names.push_back( (*rit).first );
 
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 2. Determine number of boundary segments (sub-boundaries) that the new boundary will consist of.
    //    The output of this step will be a map of FaceConstructionData in which the names of the new boundary segments are the keys.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // looping over the region, identifying and recording the juxtaposition relationships
    map<pair<long,long>,uint32_t>      patches;
    vector<FaceConstructionData<dim>>  face_construction_data;
    map<long,string>                   patch_names;
    string                             patch_name;
    uint32_t                           n_juxtapositions(0);

    // 2.1 looping over lower dimensional region identifying juxtaposition relationships
    // ----------------------------------------------------------------------------------------
    face_construction_data.reserve(subdomain.Cells());
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
      {
          // 2.1.1 identifying neighbors, facing relations, and juxtaposed materials for current element
          FaceConstructionData  fdata( higherDimensionalNeighbors( *(*eit), mtrl_key ) );
        
          // 2.1.2 recording which category of juxtaposition element fall into, naming it and assigning a patch number
          pair<map<pair<long,long>,uint32_t>::iterator,bool>  it=patches.insert( make_pair(fdata.Materials(),n_juxtapositions) );
          // incrementing number of juxtapositions and corresponding patch names
          if ( it.second == true ) {
               fdata.PatchNumber( (*it.first).second );
               patch_name = CreateBoundaryNameFrom( fdata, region_names );
               patch_names.insert( make_pair(n_juxtapositions,patch_name) );
               n_juxtapositions++;
            }
          fdata.PatchNumber( (*it.first).second );
        
          // 2.1.3 recording the data for the element that will later be used to construct the face from
          face_construction_data.push_back( fdata );
      }
    assert( face_construction_data.size() == subdomain.Cells() );

   
    // 2.2 creating labeled boundary patches from the face-defining data
    // -----------------------------------------------------------------
    // 2.2.1 making a map 'patch_numbers' from 'patch_names' to search for patch identifiers
    map<string,size_t>  patch_numbers;
    for ( const auto& it : patch_names )
      patch_numbers.insert( make_pair( it.second, static_cast<uint32_t>(it.first) ) );
   
    // 2.2.2 building new map where the patch faces are organised by patch names
    map<string,vector<FaceConstructionData<dim>> > patch_data;
    vector<FaceConstructionData<dim>>              empty_vec;
    for ( const auto& it : patch_names )
      patch_data.insert( make_pair( it.second, empty_vec ) );
   
    // 2.2.3 inserting the patch identifiers into the vectors in the map
    for ( auto& pit : patch_data )
      {
         assert( patch_numbers.find(pit.first) != patch_numbers.end() );
         const size_t patch_number((*patch_numbers.find(pit.first)).second);
         // reserving storage
         pit.second.reserve(face_construction_data.size());
         // looping over all face data assigning the ones that are suitable
         for ( const auto& it : face_construction_data )
           if ( it.PatchNumber() == patch_number )
             pit.second.push_back( it );
      }
    // 2.2.4 trimming excess storage of the face-data vectors
    for ( auto& pit : patch_data ) pit.second.shrink_to_fit();

 
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 3. Getting the MeshManager object to create Face objects for all boundary patches at the same time
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    //  3.1 creating the required face objects
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    const size_t       new_faces_required(subdomain.Cells());
    vector<Face<dim>*> face_vector;
    face_vector.reserve(new_faces_required);
#ifdef DEBUG
    const size_t n_original_faces(model.Mesh().Faces());
    const size_t n_original_elmts(model.Mesh().Elements());
#endif
    // establish the storage requirements for face variables
    const LocalVariables             lvsFaces( model.Database().LocalVariablesAt(FACE) );
    const IntegrationPointVariables  lvsIntegrationPoints( model.Database().IntegrationPointVariablesAt(FACE) );
    vector<vector<Face<dim>*> >      face_ptr_per_patch(patch_data.size());
   
    size_t patch_counter{0U};
    for ( auto& it : patch_data ) // not const because elements will be deleted
      {
         face_ptr_per_patch[patch_counter].reserve( it.second.size() );

         // for each of the new patches
         for ( auto& pit : it.second )
           {
              // creating the faces
              // ------------------
              // storing pointers to the new faces in the vector from which the boundary will be constructed
              face_vector.push_back( model.Mesh().ReplaceElementByFace( pit.LowerDimElement(),
                                                                        pit.InnerElement(),
                                                                        pit.OuterElement(),
                                                                        pit.InnerElementFace(),
                                                                        pit.OuterElementFace(),
                                                                        lvsFaces, lvsIntegrationPoints ) );
              // remembering which faces make up the patch
              face_ptr_per_patch[patch_counter].push_back( face_vector.back() );
           }
         patch_counter++;
      }
    patch_data.clear();
    
#ifdef DEBUG
    cout <<"\n\nBoundaryInterface<"<< dim <<">::CreateInternalBoundaryFrom:";
    cout << "\n\t\t"<<"Added "<< model.Mesh().Faces() - n_original_faces <<" faces to mesh.";
    cout << "\n\t\t"<<"Removed "<< n_original_elmts - model.Mesh().Elements()  <<" elements from the mesh."<< endl;
#endif
   
    //  3.2 connect them with one another (neighbors); Boundary::EstablishNeighborConnectivity( vector<Face<dim>*>& ); this is important because
    //      any ModelSubDomain creation relies on this connectivity during identification of interior and perimeter.
    //      - this method also updates node to parent element connectivity
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // TODO: these are global changes! - not sure how to improve this because so many regions are affected
    model.Mesh().UpdateConnectivity();
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 4. Create Boundary objects for each of the mesh patches established above
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    for ( auto i{0U}; i<patch_names.size(); ++i )
      AddBoundary( patch_names[i].c_str(), face_ptr_per_patch[i].begin(), face_ptr_per_patch[i].end(), INTERNAL );
      
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 5. Assign BOX_BOUNDARY flags to the nodes of each new patch by using the underlying region
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit )
      if ( (*nit)->AtBoundary() == NOT )
        (*nit)->AtBoundary(INTERNAL);
 
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 6. remove lower-dimensional input region (its elements were already removed above).
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    const bool remove_elmts{ false };
    model.RemoveRegion( dim_1_region, remove_elmts );
    model.UpdateRegions();
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 7. extra diagnostics and output of boundary names
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    if ( patch_names.empty() ) {
         ErrorHandler::Instance().Note( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region, "no Boundary patches could be created." );
         return make_pair( set<string>({}), false );
      }
      
    set<string> boundary_names;
    for ( const auto& it : patch_names ) boundary_names.insert( it.second );

    if ( patch_names.size() > boundary_names.size() ) {
         ErrorHandler::Instance().Note( ERROR, "BoundaryInterface::CreateInternalBoundaryFrom:", dim_1_region, "not all of the created patchnames are unique." );
         return make_pair( set<string>({}), false );
      }
    
    return make_pair( boundary_names, true );

 } // end CreateInternalBoundaryFrom


/* CUT OUT
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 5. Since the Face objects were constructed from lower-dimensional Elements - for nodes located on that new boundary,
    //    update / recreate the parent element vectors of the nodes on the boundary so that these no longer include
    //    neither the elements from which the Boundary was created nor the new faces.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    map<Element<dim>*,uint32_t>  parents_to_keep;
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         const auto parent_elements((*nit)->Parents());
         // copying those node parent pointers to the temporary vector which shall be kept
         for ( uint32_t i{0U}; i<parent_elements; ++i ) {
              if ( (*nit)->Parent(i)->IsSurface() and subdomain.Contains( (*nit)->Parent(i) ) )
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
*/










/**
   Removes the boundary from the map, prompting the MeshManager to delete the corresponding Face objects if so required.
   Interior nodes that had a boundary specific BOX_BOUNDARY  flag are set back to NOT.
   
   @param boundary_name either one of the box boundaries or a name that contains the string BOUNDARY (case insentivite)
   
   @attention this does not remove the Face objects associated with the boundary; call MeshManager to do this.
   
       @author SKM (refactored - since design was flawed)
       @date 15/8/2020
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim, BOUNDARY_COMPLEX>::RemoveBoundary( const char* boundary_name, bool erase_faces )
  {
    if ( !ContainsBoundary(boundary_name) ) {
         ErrorHandler::Instance().Note( WARNING, "BoundaryInterface::RemoveBoundary",
                                          boundary_name, "no boundary with this name found; nothing was done");
         return;
      }

    // for interior nodes with a boundary specific flag this flag is set to NOT
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    csmp::Boundary<dim>& boundary = boundaryComplex->Boundary( boundary_name );
    for ( auto nit=boundary.NodesBegin(); nit!=boundary.PerimeterNodesBegin(); ++nit )
      if ( (*nit)->AtBoundary() == boundary.AtBoundary() )
        (*nit)->AtBoundary( NOT );
    
     // erasing the faces
     if ( erase_faces ) {
         // getting the mesh manager to delete faces and nodes and fix up the connectivity
         boundaryComplex->Mesh().DeleteCellsAndRepairConnnectivity( boundary.CellVector().begin(), boundary.CellVector().end() );
       }

    // erasing the boundary
    boundaryMap_.erase( boundary_name );
    
  } // end RemoveBoundary



/**
   Removes the boundary from the map, prompting the MeshManager to delete the corresponding Face objects if so required.
   Interior nodes that had a boundary specific BOX_BOUNDARY  flag are set back to NOT.
         
   @param boundary reference to valid boundary in the model

   @attention this does not remove the Face objects associated with the boundary; call MeshManager to do this.
   
       @author SKM (refactored - since design was flawed)
       @date 15/8/2020
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim, BOUNDARY_COMPLEX>::RemoveBoundary( csmp::Boundary<dim>& boundary, bool erase_faces )
  {
    // for interior nodes with a boundary specific flag this flag is set to NOT
     for ( auto nit=boundary.NodesBegin(); nit!=boundary.PerimeterNodesBegin(); ++nit )
       if ( (*nit)->AtBoundary() == boundary.AtBoundary() )
         (*nit)->AtBoundary( NOT );

     // erasing the faces
     if ( erase_faces ) {
         BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
         // getting the mesh manager to delete faces and nodes and fix up the connectivity
         boundaryComplex->Mesh().DeleteCellsAndRepairConnnectivity( boundary.CellVector().begin(), boundary.CellVector().end() );
       }

     // deleting the Boundary
     boundaryMap_.erase( boundary.Name() );

   } // end RemoveBoundary





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
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::OutputBoundariesToBinary( const char* file_name ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    const BOUNDARY_COMPLEX<dim>& boundaryComplex( static_cast<const BOUNDARY_COMPLEX<dim>& >(*this) );
    
    string  bin_file(file_name);
    fstream fp(bin_file.c_str(), ios::out | ios::binary);
    if ( !fp.is_open() ) {
         csmp_error.Note( ERROR, "BoundaryInterface::OutputBoundariesToBinary:",
                            bin_file, "file could not be opened; nothing was done." );
         return false;
      }
     // 0. writing file header
     {
        BinaryFileSectionWrite sect(fp, "BNDFHEDR");
        string heading("BoundaryInterface::OutputBoundariesToBinary: ");
        heading +="boundary information for Model '";
        heading += boundaryComplex.Name();
        heading +="' to file: ";
        heading += bin_file;
        heading +="'.";
        binaryFileWrite( fp, heading.c_str() );
     }

     cout <<"\nBoundaryInterface<"<< dim <<">::OutputBoundariesToBinary: boundaries written to binary file: ";
     // 1. writing all boundary objects to binary file
     {
       BinaryFileSectionWrite sect(fp, "BOUNDARY");

       const uint64_t  records(this->Boundaries());
       fp.write( reinterpret_cast<const char*>(&records), sizeof(uint64_t ) );

       for ( typename map<string,csmp::Boundary<dim> >::const_iterator
             git=BoundariesBegin(); git!=BoundariesEnd(); ++git )
         {
            BinaryFileSectionWrite hdr(fp, "ONE_BDRY");
            cout <<"'"<< (*git).first <<"' ";
            cout.flush();
            // 1.1 writing the entire connectivity structure to the binary file
            (*git).second.WriteDomainIndexesToBinaryFile( fp );
            // 1.2 writing the boundary flags
            auto bflag = (*git).second.AtBoundary();
            const int64_t  record(1U);
            fp.write( reinterpret_cast<const char*>(&record), sizeof( int64_t  ) );
            fp.write( reinterpret_cast<const char*>(&bflag), sizeof( int8_t ) );
            // 1.3 writing the stored variables
            domainVariablesOut( fp, (*git).second, boundaryComplex.Database() );
         }
     }

    // 2. footer, and clean up
    BinaryFileSectionWrite sect(fp, "BNDFFOTR");

    fp.close();
    cout <<"\n\nBoundaryInterface<"<< dim <<">::OutputBoundariesToBinary: file '";
    cout << bin_file <<"' has been successfully written.\n";
   
    return true;

} // end OutputBoundariesToBinary




/**
     reads and initialises boundaries from file written by OutputBoundariesToBinary()
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim,BOUNDARY_COMPLEX>::InputBoundariesFromBinary( const char* file_name,
                                                                         const set<string>& subset_variables )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );

   string  bin_file(file_name);
	 fstream fp(bin_file.c_str(), ios::in | ios::binary);
	 if (!fp.is_open()) {
          csmp_error.Note( ERROR, "BoundaryInterface::InputBoundariesFromBinary:",
                             bin_file, "file could not be opened; nothing was done." );
          return;
       }
   
     // 0. reading the file header and printing it to screen
     {
       BinaryFileSectionRead sect(fp, "BNDFHEDR");

       char  text[INFO_STRING];
       binaryFileRead( fp, text );
       cout <<"\nBoundaryInterface<"<< dim <<">::InputBoundariesFromBinary: Reading file header:\n\t"<< text << endl;
     }
     cout <<"\n\timporting the boundaries: ";
   
     // 1. reading the boundaries
     BOUNDARY_COMPLEX<dim>&       boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>& >(*this) );
     const PropertyDatabase<dim>& database( boundaryComplex.Database() );
     MeshManager<dim>&            mesh( boundaryComplex.Mesh() );

     {
       BinaryFileSectionRead sect(fp, "BOUNDARY");
       
       SubDomainInfo  info;
       uint64_t   records(0);  // region records
       // getting number of unique region records from file
       fp.read( reinterpret_cast<char*>(&records), sizeof(uint64_t ) );
       if ( records > 0 )
          // reading the regions sequentially
          for ( auto i{0U}; i<records; ++i )
            {
               BinaryFileSectionRead hdr(fp, "ONE_BDRY");

               // 1.1 reading name and face indices for each boundaries
               readDomainIndexesFromBinaryFile( dim, fp, info );
              
               // 1.2 reading BOX boundary flag of the boundary
               int64_t  record;
               fp.read( reinterpret_cast<char*>(&record), sizeof(int64_t ) );
               assert( record == 1 );
               int8_t box_boundary_index(IRREGULAR_OUTSIDE);
               fp.read( reinterpret_cast<char*>(&box_boundary_index), sizeof(int8_t) );
               BOX_BOUNDARY bflag = intToBOX_BOUNDARY( box_boundary_index );
               
               // 1.3 creating the boundary objects
               pair<typename map<string,csmp::Boundary<dim> >::iterator,bool>
                 it=boundaryMap_.insert( make_pair( info.name, csmp::Boundary<dim>(database, mesh, info, bflag) ) );
              
               if ( !it.second )
                 throw csmp::Exception( FATAL_ERROR, "BoundaryInterface::InputBoundariesFromBinary:",
                                        info.name, "Boundary could not be formed; issue with binary file." );
           
               // 1.4 reading the variable values
               if ( subset_variables.empty() ) domainVariablesIn( fp, (*it.first).second, database );
               else selectedDomainVariablesIn( fp, (*it.first).second, database, subset_variables );
           
               // 1.5 reporting to stdout
               cout <<"\n\t\t"<< (*it.first).first <<" ("<< parseBoundary(bflag) <<", "<< (*it.first).second.Cells() <<" faces).";
           }
    }

    // 2. cleaning up
    BinaryFileSectionRead sect(fp, "BNDFFOTR");

    fp.close();
    cout <<"\n\nBoundaryInterface<"<< dim <<",Face>::InputBoundariesFromBinary: file '";
    cout << bin_file <<"' has been read successfully.\n";

} // end InputBoundariesFromBinary





/**
    Method forms a Boundary between the two supplied regions. This will involve the creation
    and connection of Faces.

    This method will only work for unique Regions which are not overlapping.

    This method makes no sense for the region 'Model' as it encompasses all unique regions.
    
    @return returns pair singnalling the success of the operation and giving the name of the boundary.
    
    @attention:  BoundariesInterface cannot be created in 1D models or between regions which contain
    one-dimensional elements.

    @note the BOX boundary flag will always be internal for this type of boundary because it
    lies between higher dimensional regions.
    
    NAMING THE NEW BOUNDARY
    
    (1. name of the master region from which the boundary was created if any)
    
    2. "BOUNDARY"
    
    3. the patch identifier number attached to boundary (making the boundary patch unique, avoiding duplications)
    
    4. the name of the inner region, i.e. the region that the lower-dimensional element normals point away from
    
    5. the name of the outer region, i.e. that into which the normals point
    
    @note Where the boundary just intersects a layer (same material on either side), the layer name appears
    only once. The second instance is replaced by INTERSECTION.
    
    @section implementation Implementation

    Attempts to create a lower-dimensional region between higher dimensional ones.
    This is done in the following steps:

    0. Checks:
    - do the input regions exist
    - are they higher dimensional
    - is there not already a region that has the name that the new region will get?

    1. Using the Perimeter faces of the candiate regions, matching faces are found and the inside and outside
    elements are determined as well as recording their face numbers.

    2. The MeshManager is instructed to create the required Face objects AND connect them with one another.
    (no objects need to be deleted because it is assumed that there is no line element region at this boundary) @todo check

    3. The Boundary is constructed from the Face objects.

    @note the Region will be oriented such that the elements of the first region will be on the inside (normals pointing outward from this region).

    @note this region is not necessarily contiguous

    @author SKM
    @date refactored 2/4/2022
    
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
pair<string,bool>  BoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateBoundaryBetween( const char* group1, const char* group2 )
 {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    string  region1(group1);
    string  region2(group2);

    // 1. preliminary checks
    // ---------------------
    if (region1 == "Model" or region2 == "Model") {
         csmp_error.Note( ERROR, "BoundaryInterface::CreateBetween:", "Region 'Model' not eligible for InsertBoundary.");
         return make_pair("boundary not created",false);
      }
    if (region1 == region2) {
         csmp_error.Note( ERROR, "BoundaryInterface::CreateBetween:", "Provided Regions are identical.");
         return make_pair("boundary not created",false);
      }
    if ( !boundaryComplex->IsUnique(group1) || !boundaryComplex->IsUnique(group2) ) {
         csmp_error.Note( ERROR, "BoundaryInterface::CreateBetween:",
                           "This method is intended for the creation of unique non-overlapping Regions");
         return make_pair("boundary not created",false);
      }
    // TODO: test whether cases are handled correctly where one of the regions is a surface or a line while the other is a volume?
   
   
    // 2. do the two regions share part of (are in contact with eachother) their perimeter
    // -----------------------------------------------------------------------------------
    // (this also considers the case where one of the regions has no usable interface because it is at the model boundary)
    const csmp::Region<dim>&  gref1(boundaryComplex->Region(group1));
    const csmp::Region<dim>&  gref2(boundaryComplex->Region(group2));
    vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >  matching_element_faces;
    const size_t contacting_faces = sharedPerimeterCells( gref1, gref2, matching_element_faces );

    if ( contacting_faces == 0U ) {
         csmp_error.Note( ERROR, "BoundaryInterface::CreateBoundaryBetween:",
                           "the supplied regions have no contacting faces; boundary could not be created");
         return make_pair("boundary not created",false);
      }

    // 3. creating a boundary from the faces of the shared perimeter elements
    // ----------------------------------------------------------------------
    string boundary_name = CreateBoundaryName( region1, region2 );
    
    // 4. creating the Face objects for the Boundary
    // ---------------------------------------------
    // (includes the update of the node to parent connectivity etc.)
    vector<Face<dim>*> face_ptrs = boundaryComplex->Mesh().CreateFacesBetweenNodeSharingElements( boundaryComplex->Database(), matching_element_faces );

    // 5. creating the the Boundary from the Face objects
    // --------------------------------------------------
    cout << "\nBoundaryInterface<"<< dim <<">::CreateBoundaryBetween: creating boundary between ";
    cout << group1 << " and " << group2 << endl;
    AddBoundary( boundary_name.c_str(), face_ptrs.begin(), face_ptrs.end(), INTERNAL );

    return make_pair( boundary_name, true );
    
 } // end CreateBoundaryBetween








/** Creates Boundary around the supplied region ignoring pre-exisiting Boundaries of SplitBoundaries
 identified by nodes flagged as INTERNAL or the presence of NodeManifolds.
 Also ignores external model boundaries.
 
 The resulting Boundary is flagged as INTERNAL.
 
 @attention By default, there is no partitioning of the boundary into patches based in the juxtaposition
 of different layers and the region. Instead, the boundary is named BOUNDARY_"region"_HULL.
 
 @attention the detection of existing boundaries relies on their nodes being flagged as !NOT, if these are
 not present the method might duplicate existing boundaries.

 @todo should HULL be partitioned for consistency with boundary naming conventions? - this would require
 the identification of the names of the adjacent regions. It could be made optional / dependent on the presence
 of an initialised region identifier.

  @author SKM
  @date 27/10/21 refactored
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateBoundaryAround( const char* region_name )
  {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 1. initial checks
    // -----------------
    if ( !boundaryComplex->ContainsRegion(region_name) ) {
         csmp_error.Note( ERROR, "BoundaryInterFace::CreateBoundaryAround:", region_name, "region does not exist. Nothing could be done" );
         return false;
      }
    const Region<dim>& region( boundaryComplex->Region(region_name) );
      
    // 2. Creation of unique boundary name
    // -----------------------------------
    // (BOUNDARY_region_HULL)
    int count{0};
    string boundaryName = string("BOUNDARY_") + region_name + "_HULL";
    while ( boundaryMap_.count( boundaryName ) > 0 ) {
         boundaryName  = "BOUNDARY" + to_string( count++ );
         boundaryName += string("_") + region_name + "_HULL";
      }

    // 3. Finding the perimeter faces from which the boundary shall be created, including their higher dimensional neighbors
    // ---------------------------------------------------------------------------------------------------------------------
    vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >  hull_elmt_faces;
    hull_elmt_faces.reserve( region.PerimeterCells() );
    
    // looping over the perimeter faces of the region
    for ( size_t i{ region.InteriorCells() }; i<region.Cells(); i++ )
      for ( auto j{0U}; j<region.PerimeterFaces(i); j++ ) {
           auto inside_face = region.PerimeterFace(i,j);
           assert( !region.E(i)->IsLine() );
           assert( inside_face < region.E(i)->Neighbors() );
           // checking whether the corner nodes of the face lie on a boundary (in this case the face is discarded)
           bool is_boundary_face{true};
           auto fnodes = region.E(i)->CornerNodesOfFace( inside_face );
           for ( const auto& nit : fnodes )
             if ( nit->AtBoundary() == NOT && nit->IsManifold() == false ) {
                  is_boundary_face = false;
                  break;
               }
           // if this is an interior face that is not already a boundary or a split boundary,
           // the element on the outside of the region is searched for
           if ( !is_boundary_face ) {
                auto outside_eptr = region.E(i)->Neighbor( inside_face );
                if ( outside_eptr != nullptr ) {
                     // this is a face eligible for hull creation and we find the number of the face of the outside
                     // element that matches it
                     uint32_t outside_face{ UINT_MAX };
                     for ( auto face{0U}; face<outside_eptr->Neighbors(); face++ )
                       if ( outside_eptr->Neighbor(face) == region.E(i) ) {
                            outside_face = face;
                            break;
                         }
                     assert( outside_face < outside_eptr->Faces() );
                     
                     // now all data is there to create an entry into
                     hull_elmt_faces.push_back( make_pair( make_pair( region.E(i), inside_face ),
                                                           make_pair( region.E(i)->Neighbor(inside_face), outside_face ) ) );
                  }
             }
           
        }
    
    if ( hull_elmt_faces.empty() ) {
         csmp_error.Note( ERROR, "BoundaryInterFace::CreateBoundaryAround:", region_name, "No eligible outside faces were found." );
         return false;
      }

    // 4. getting the MeshManager to create and interconnect the missing Face objects
    // ------------------------------------------------------------------------------
    vector<Face<dim>*> face_ptrs = boundaryComplex->Mesh().CreateFacesBetweenNodeSharingElements( boundaryComplex->Database(), hull_elmt_faces );
    

    // 5. inserting a corresponding boundary if not existing yet
    // ---------------------------------------------------------
    AddBoundary( boundaryName.c_str(), face_ptrs.begin(), face_ptrs.end(), INTERNAL );

    // shouldn't get here
    return true;
    
  } // end CreateAround









/**
    Loops over the perimeter of the boundary, making keys from the node-pointers of the boundary face
    nodes and recording the boundary elements for the faces of which the keys were made
     
     @attention makes sense in 3D only
*/
static void createPerimeterKeysFor( const Boundary<3U>& boundary, map<set<csmp::Node<3U>*>,Face<3U>*>& perimeter_keys )
 {
    if ( !perimeter_keys.empty() ) perimeter_keys.clear();
   
    // looping over the Face edges on the boundary, creating the keys from sets of node pointers
    for ( size_t i=boundary.InteriorCells(); i<boundary.Cells(); ++i )
      for ( auto j{0U}; j<boundary.PerimeterFaces(i); ++j )
        {
            const auto      pface = boundary.PerimeterFace(i,j);
            set<Node<3U>*>  key;
            for ( const auto& k : boundary.E(i)->FE()->NodesOfFace(pface) )
              key.insert( boundary.E(i)->N(k) );
            perimeter_keys.insert( make_pair(key,boundary.E(i)) );
              
        }
   
 } // end createPerimeterKeysFor




/**
    Connects the line faces representing the edge with their neighbors.
    logic: where the line elements share a node they are connected.
    
    @note makes sense only in 3D.
    
    @todo TODO: deal with the special case of manifolds. At the corner
    of a box-shaped model, for instance, 3 line-element faces contact each other.
    
    @attention where multiple line element faces connect at a Node,
    this method assigns a nullptr to that face.
    
    @attention this method assumes that the first 2 nodes of the line-element
    Face are the end-point nodes
*/
static void createLineFaceConnectivity( vector<Face<3U>*>& line_faces )
 {
    if ( line_faces.empty() ) return;
   
    // 1. making a map of the parent faces that each node is connected to
    //  key       faces that are connected to the node (should be 2 at most)
    map<Node<3U>*,set<Face<3U>*> > parent_faces;
   
    for ( auto& it : line_faces ) {
         set<Face<3U>*> parents{it};
         for ( auto i{0U}; i<it->Nodes(); ++i ) {
              // inserting a new set or inserting a face pointer into the set if the node key already exists
              auto nit = parent_faces.insert( make_pair( it->N(i), parents ) );
              if ( !nit.second )
                (*nit.first).second.insert(it);
           }
      }
   
   // 2. connecting the faces with one another
   //   (the assumption is that each face has 1-2 equidimensional neighbors that
   //    coincide with its corner nodes!=midside nodes if any)
   ErrorHandler& csmp_error(ErrorHandler::Instance());

   for ( auto& it : parent_faces ) {
        // 2.1 nodes/faces at the end-points of the edge = perimeter faces
        // (there may only be one neighbor or a manifold interpreted as endpoint)
        if ( it.second.size() != 2 ) {
          if ( it.second.size() > 2 ) {
               // line faces only have a role on model edges
               // ony those edge node parents are kept in the map which have two nodes flagged as edge
               for ( auto iit=it.second.begin(); iit!=it.second.end(); ++iit ) {
                    for ( auto j{0U}; j<(*iit)->Nodes(); ++j )
                      if ( !isEdge( (*iit)->N(j)->AtBoundary() ) ) {
                           // returns a valid iterator to the set after the erasure
                           iit = it.second.erase( iit );
                           break;
                        }
                    if ( iit == it.second.end() ) break;
                 }
               if ( it.second.size() != 2 )
                 csmp_error.Note( ERROR, "creatLineFaceConnectivity(Face):",
                                   "edge node is connected to more than 2 line Faces;\
                                    don't know how to deal with this manifold.");
            }
           /* NOTHING EXTRA NEEDS TO BE DONE BECAUSE FACE POINTERS ALREADY ARE NULLPTRs
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
             for ( auto i{0U}; i<edge1->Nodes(); ++i )
               if ( it.first == edge1->N(i) ) {
                    // assigning the opposite neighbor
                    edge1->Assign( i, edge2 );
                    break;
                 }
             // finding node number in second Face
             for ( auto i{0U}; i<edge2->Nodes(); ++i )
               if ( it.first == edge2->N(i) ) {
                    // assigning neighbors
                    edge2->Assign( i, edge1 );
                    break;
                 }
          }
     }
   
 } // create line element neighbor connectivity



/**
    Creates a model boundary EDGE#  of Face objects that are line elements and have the boundary surface elements as higher dimensional neighbors.
     
     The method does not require the presence of line elements on the edge of interest but finds it by matching the edges of boundary perimeter faces.

    @param shared_faces is a vector of the Face objects that were created in the MeshManager in the boundary creation process
 
    @return if the operation was successful
*/
static bool createBoundaryFromSharedEdge( Model<3U>& model,
                                          const Boundary<3U>& boundary1, const Boundary<3U>& boundary2,
                                          vector<Face<3U>*>&  shared_faces )
 {
    if ( !shared_faces.empty() ) shared_faces.clear();
 
    // 1. creating keys for the perimeter element faces of boundary1
    // ------------------------------------------------------------------
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
               // 2. finding parent elements of Faces and their segments that coincide with the new line Face
               // -------------------------------------------------------------------------------------------
               Element<3U>*  parent1 = (*it1).second->InnerParent();
               Element<3U>*  parent2 = (*it2).second->InnerParent();
               assert( parent1 != nullptr );
               if ( parent1 == parent2 ) parent2 = nullptr;
               vector<Node<3U>*>  segment_nodes;
               uint32_t segment_id_parent1 = UNSPECIFIED;
               uint32_t segment_id_parent2 = UNSPECIFIED;
               // 2.1 parent element of Face 1
               // ----------------------------
               // establishing the face-node sequence of the inner element face that will be shared with the new Face object
               for ( auto i{0U}; i<parent1->Segments(); ++i ) {
                    vector<uint32_t> snids; // local segment node ids
                    parent1->FE()->NodesOfSegment( i, snids );
                    set<csmp::Node<3U>*> nset;
                    for ( size_t k{0U}; k<snids.size(); ++k ) nset.insert( parent1->N(snids[k]) );
                    if ( nset == (*it1).first ) {
                          segment_id_parent1 = i;
                          // capturing the segment nodes for the construction of the Face object
                          segment_nodes.reserve( snids.size() );
                          for ( size_t j{0U}; j<snids.size(); ++j )
                            segment_nodes.push_back( parent1->N(snids[j]) );
                          assert( segment_nodes.size() >= 2U );
                          // 2.2 parent of Face 2
                          // --------------------
                          // finding which segment this in parent2 if it exists
                          if ( parent2 != nullptr ) {
                               for ( auto j{0U}; j<parent2->Segments(); ++j ) {
                                    vector<uint32_t> snids2; // local segment node ids
                                    parent2->FE()->NodesOfSegment( j, snids2 );
                                    set<csmp::Node<3U>*> nset2; // search set of node pointers
                                    for ( auto n=0U; n<snids2.size(); ++n ) nset2.insert( parent2->N(snids2[n]) );
                                    if ( nset == nset2 ) {
                                         segment_id_parent2 = j;
                                         break;
                                      }
                                 }
                            }
                          break;
                      }
                 }
        
               // creating the face in MeshManager
			         Face<3U>* const faceObj = model.Mesh().AddEdgeFace( (*it1).second, segment_id_parent1,
                                                                   (*it2).second, segment_id_parent2,
                                                                   lvsFaces, lvsIntegrationPoints,
                                                                   segment_nodes );
			        shared_faces.push_back(faceObj);
           }
      }
   
    // assigning the equidimensional neighbors to the newly created faces
    // createLineFaceConnectivity( shared_faces ); // works fine but this is the task of the MeshManager
    
    // interconnecting the edges with one-another
    model.Mesh().BuildLineConnectivity<Face>( shared_faces.begin(), shared_faces.end() );
		
    // false if no shared faces could be detected
    return ( !shared_faces.empty() );

 } // end createBoundaryFromSharedEdge





  // High Level Functions to create Boundaries from provided Region names





/** creates edge Boundary objects (of dim-2 Face objects) for box-shaped model from side boundaries
   
   1. verifies that that the side boundaries of the box-shaped model are there
   
   @attention EXCEPTION: the edge curves are represented by Face objects, however these 
   are dim-2 entities. Following CSMP's rules, these edges will be connected to their face neighbors
   AND their higher-dimensional parent elements with which they share edges.
   
   @test OK - for BoxHalfs3D

*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel()
 {
    if constexpr ( dim != 3U )
      throw csmp::Exception( ERROR, "BoundaryInterface<::EstablishEdgeBoundariesOfBoxShapedModel: ","edges required only in 3D; nothing was done" );
     
    if constexpr ( dim == 3U ) {
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
     
      ErrorHandler& csmp_error(ErrorHandler::Instance());

      // EDGE1 = BACK_BOTTOM
      // -------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, bottom, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE1 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE1", shared_faces.begin(), shared_faces.end(), EDGE1 );
      if ( shared_faces.empty() ) return_value=false;
   
      // EDGE2 = BACK_RIGHT
      // ------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, right, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE2 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE2", shared_faces.begin(), shared_faces.end(), EDGE2 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE3 = BACK_TOP
      // ----------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, top, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE3 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE3", shared_faces.begin(), shared_faces.end(), EDGE3 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE4 = BACK_LEFT
      // -----------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, back, left, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE4 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE4", shared_faces.begin(), shared_faces.end(), EDGE4 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE5 = BOTTOM_LEFT
      // -------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, bottom, left, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE5 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE5", shared_faces.begin(), shared_faces.end(), EDGE5 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE6 = BOTTOM_RIGHT
      // --------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, bottom, right, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE6 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE6", shared_faces.begin(), shared_faces.end(), EDGE6 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE7 = TOP_RIGHT
      // -----------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, top, right, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE7 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE7", shared_faces.begin(), shared_faces.end(), EDGE7 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE8 = TOP_LEFT
      // ----------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, top, left, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE8 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE8", shared_faces.begin(), shared_faces.end(), EDGE8 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE9 = FRONT_BOTTOM
      // --------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, bottom, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE9 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE9", shared_faces.begin(), shared_faces.end(), EDGE9 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE10 = FRONT_RIGHT
      // --------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, right, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE10 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE10", shared_faces.begin(), shared_faces.end(), EDGE10 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE11 = FRONT_TOP
      // ------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, top, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE11 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE11", shared_faces.begin(), shared_faces.end(), EDGE11 );
      if ( shared_faces.empty() ) return_value=false;

      // EDGE12 = FRONT_LEFT
      // -------------------
      if ( !createBoundaryFromSharedEdge( *boundaryComplex, front, left, shared_faces ) )
        csmp_error.Note( WARNING, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishEdgeBoundariesOfBoxShapedModel:",
                          "Boundary EDGE12 could not be created because no shared edge was found." );

      else AddBoundary( "EDGE12", shared_faces.begin(), shared_faces.end(), EDGE12 );
      if ( shared_faces.empty() ) return_value=false;
    
      return return_value; // whether all edges could be established
   }
   
 return false; // whether all edges could be established
   
} // end EstablishEdgeBoundariesOfBoxShapedModel






/** Helper function for method EstablishBoxBoundaries() below
    returns index of first and last element of the checked region
*/
template<uint32_t dim>
static pair<size_t,size_t>  collectLowerDimensionalElementsFrom( Model<dim>& model, const char* region_name,
                                                                 vector<Element<dim>*>& elements )
 {
    ErrorHandler& csmp_error(ErrorHandler::Instance());
    
    const size_t first_index{ elements.size() };
    size_t       last_index{ first_index };
    
    if ( !model.ContainsRegion( region_name ) ) {
         csmp_error.Note( WARNING, "collectLowerDimensionalElementsFrom", region_name,
                           "does not exist and could therefore not be added vector");
                           
         return make_pair(first_index,first_index);
      }
     
    Region<dim>& domain(model.Region(region_name));
    if ( elements.empty() ) elements.reserve( domain.Cells() );
   
    // indexing and storing the cells for later identification
    for ( auto& it : domain.CellVector() ) {
        // checking that we are indeed dealing with a lower-dimensional element
        if ( (dim == 3 && !it->IsSurface()) || (dim == 2 && !it->IsLine()) ) {
            cerr <<"\n\t"<< parseFiniteElementType( it->FE_Type() ) <<": idx: "<< it->Idx();
            csmp_error.Note( ERROR, "collectLowerDimensionalElementsFrom", "element is not lower dimensional");
          }
        else {
             elements.push_back( it );
             last_index++;
          }
     }

   return make_pair( first_index, last_index );
   
 } // end collectLowerDimensionalElements




/**
     Forms boundaries considering names of Box.h-defined strings only, e.g.,
     BACK, BOTTOM, RIGHT, TOP, FRONT, IRREGULAR and INTERNAL. Regions with these names are converted to boundaries
     removing them and their Elements that are replaced by Face objects.
     The new boundaries are given the standard names and are flagged correspondingly.
     
     @remark input regions must be lower-dimensional (surfaces in 3D and lines in 1D)
     @remark input regions must be unique (space exclusive)
     @remark connectivity of Face objects must match / preserve that of previous Element objects
     @remark when elements become faces they need to be tagged so that one can detect which region they were created from (idx?)
     
     @note if the TOP region is missing, but a IRREGULAR region is there in stead, this is converted into the corresponding boundary.
     In this case, the model is still regarded as BOX_SHAPED.
     
     @note Boundary creation itself does not deal with the generation of BOX_BOUNDARY flags for the model edges and
     corners. This is accomplished subsequently (in this method) by calling recreateBoxBoundaryFlags().
     
     @attention this method does not take care of the updating of the non-unique regions that are affected by the conversion
     of Regions into boundaries. This is done afterwards by methods of the Model class.
     
     @author refactored by SKM 2016
     @author refactored by SKM 2018
     @author SKM 2021
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishBoxBoundaries()
   {
      ErrorHandler& csmp_error(ErrorHandler::Instance());

      BOUNDARY_COMPLEX<dim>* model( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
      
      cout << "\nBoundaryInterface<"<< dim <<">::EstablishBoxBoundaries: Establishing Box-object boundaries for ";
      cout << dim << " dimensional box shaped model...";
     
      // 1. do diagnostics, creating vector of iterators for the elements that shall be replaced by Faces
      // ------------------------------------------------------------------------------------------------
      vector<Element<dim>*> elmts_to_become_faces;
      elmts_to_become_faces.reserve( model->Mesh().Elements() );

      // 2. adding the input regions
      // ---------------------------
	    // TOP may be missing, but if so, there must be an IRREGULAR boundary instead
	    if ( !model->ContainsRegion("TOP") && model->ContainsRegion("IRREGULAR") ) {
           csmp_error.Note( ERROR, "BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishBoxBoundaries",
                             "model is not box-shaped, nothing was done");
           return false;
        }
      pair<size_t,size_t> top{0,0}, irregular{0,0}, front{0,0}, back{0,0};
      
      if ( model->ContainsRegion("TOP") )
        top = collectLowerDimensionalElementsFrom( *model, "TOP", elmts_to_become_faces );

      if ( model->ContainsRegion("IRREGULAR") )
        irregular = collectLowerDimensionalElementsFrom( *model, "IRREGULAR", elmts_to_become_faces );
        
      auto bottom = collectLowerDimensionalElementsFrom( *model, "BOTTOM", elmts_to_become_faces );
      auto right  = collectLowerDimensionalElementsFrom( *model, "RIGHT", elmts_to_become_faces );
      auto left   = collectLowerDimensionalElementsFrom( *model, "LEFT", elmts_to_become_faces );
      
      if constexpr ( dim == 3U ) {
           front = collectLowerDimensionalElementsFrom( *model, "FRONT", elmts_to_become_faces );
           back  = collectLowerDimensionalElementsFrom( *model, "BACK", elmts_to_become_faces );
        }

// TODO: include edges into the transformations if there are any in the input (-regions.txt) file

      // 3. getting MeshManager to create faces and delete pre-cursor elements
      // ---------------------------------------------------------------------
      vector<Face<dim>*> faces = model->Mesh().ReplaceBoundaryElementsByFaces( model->Database(),
                                                                               elmts_to_become_faces.begin(),
                                                                               elmts_to_become_faces.end() );
      // 4. creating the Boundaries from the faces
      // -----------------------------------------
      typename vector<Face<dim>*>::iterator fit{ faces.begin() };

      AddBoundary( "BOTTOM", next(fit,bottom.first), next(fit,bottom.second), BOTTOM );
      AddBoundary( "RIGHT",  next(fit,right.first),  next(fit,right.second), RIGHT );
      AddBoundary( "LEFT",   next(fit,left.first),   next(fit,left.second), LEFT );
      if ( top.first != top.second )
        AddBoundary( "TOP", next(fit,top.first), next(fit,top.second), TOP );
      if ( irregular.first != irregular.second )
        AddBoundary( "IRREGULAR", next(fit,irregular.first), next(fit,irregular.second), IRREGULAR );
      if constexpr( dim == 3U ) {
          AddBoundary( "BACK",  next(fit,back.first),  next(fit,back.second), BACK );
          AddBoundary( "FRONT", next(fit,front.first), next(fit,front.second), FRONT );
       }
      
      // 5. removing input regions and updating other regions
      // ------------------------------------------------------------------------------------
      // (no flagging for rebuilt of regions is necessary as they will be completely removed)
      const bool erase_elements{ false }; // this was already done above
      if ( top.first != top.second ) model->RemoveRegion( "TOP", erase_elements );
      if ( irregular.first != irregular.second ) model->RemoveRegion( "IRREGULAR", erase_elements );
      model->RemoveRegion( "BOTTOM", erase_elements );
      model->RemoveRegion( "RIGHT", erase_elements );
      model->RemoveRegion( "LEFT", erase_elements );
      if constexpr( dim == 3U ) {
           model->RemoveRegion( "BACK", erase_elements );
           model->RemoveRegion( "FRONT", erase_elements );
        }
        
      model->UpdateRegions();
      cout << "\n\n EstablishBoxBoundaries: done!\n";
      return true;
      
  } // EstablishBoxBoundaries








/**
    converts lower-dimensional Region on the outside of the model into a Boundary; returns whether this conversion was successful as well as the boundary name
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
pair<string,bool>  BoundaryInterface<dim, BOUNDARY_COMPLEX>::CreateExternalBoundaryFrom( const char* dimension_minus1_region, bool check_topo_attributes_of_nodes )
  {
	  BOUNDARY_COMPLEX<dim>* model(static_cast<BOUNDARY_COMPLEX<dim>*>(this));

    // 0. does the Region exist and is it lower dimensional?
    // -----------------------------------------------------
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    // is the region there
    if ( !model->ContainsRegion(dimension_minus1_region) ) {
          csmp_error.Note( WARNING, "BoundaryInterface::CreateExternalBoundaryFrom",
                                     "unable to find lower-dimensional regions to create Boundary from");
         return make_pair( string(dimension_minus1_region) + " not found", false );
      }
    Region<dim>& subdomain = model->Region(dimension_minus1_region);
    // is it lower dimensional?
    if ( !hasLowerDimensionalRepresentation(subdomain) ) {
         csmp_error.Note( WARNING, "BoundaryInterface::CreateExternalBoundaryFrom",
                          dimension_minus1_region, "is not lower-dimensional (dim-1)");
         return make_pair( string(dimension_minus1_region) + " not lower dimensional", false );
      }
    // is the region located at the outer boundary of the model
    if ( !formsPartOfExternalBoundary(subdomain,check_topo_attributes_of_nodes) ) {
         csmp_error.Note( WARNING, "BoundaryInterface::CreateExternalBoundaryFrom",
                          dimension_minus1_region, "is not entirely located on an external boundary of the model");
         return make_pair( string(dimension_minus1_region) + " not on external boundary", false );
       }
    
    // 1. making a map of the region's elements that will be converted to faces
    // -----------------------------------------------------------------------
    vector<Element<dim>*>  elmts_to_become_faces;
    elmts_to_become_faces.reserve( subdomain.Cells() );
    elmts_to_become_faces.insert( elmts_to_become_faces.end(),
                                  subdomain.CellsBegin(), subdomain.CellsEnd() );

    // 2. replacing the elements by Faces (input elements are deleted and nullptrs returned)
    // -------------------------------------------------------------------------------------
    assert( connectivityCheck<dim>( elmts_to_become_faces.begin(), elmts_to_become_faces.end() ) == 0 ); // zero means that there are no issues
    vector<Face<dim>*> faces = model->Mesh().ReplaceBoundaryElementsByFaces( model->Database(),
                                                                             elmts_to_become_faces.begin(),
                                                                             elmts_to_become_faces.end() );
    // 3. creating the Boundary from the faces
    // ---------------------------------------
    // create boundary name by appending '_BOUNDARY' to the original name of the region
    const string boundary_name( string(dimension_minus1_region) + "_BOUNDARY" );
    
    // create boundary and trying to find suitable  BOX_BOUNDARY flag for it
    BOX_BOUNDARY boundary_flag = parseBoundary( boundary_name.c_str() );
    if ( boundary_flag == MULTIPLE ) boundary_flag = IRREGULAR; // outside
    if ( !AddBoundary( boundary_name.c_str(), faces.begin(), faces.end(), boundary_flag ) ) {
         csmp_error.Note( WARNING, "BoundaryInterface::CreateExternalBoundaryFrom",
                          boundary_name.c_str(), "could not be created");
         return make_pair( string(dimension_minus1_region) + " could not be created", false );
      }
    cout << "\n\nBoundaryInterface::CreateExternalBoundaryFrom: created external boundary '";
    cout <<" "<< boundary_name <<"' successfully.";
    cout << endl;

	  // 4. removing the original regions from which the boundaries were created
    // ------------------------------------------------------------------------------------
    // (no flagging for rebuilt of regions is necessary as they will be completely removed)
    const bool erase_elements{ true };
    model->RemoveRegion( dimension_minus1_region, erase_elements );
    model->UpdateRegions();

	  cout <<"\n\n"<<"BoundaryInterface::CreateExternalBoundaryFrom: removed input region."<< endl;

	  return make_pair( boundary_name, true );

  } // end CreateExternalBoundaryFrom










  /**
  Creates boundary around "Model" region determining whether lower-dimensional
  regions shall become a boundary segments on the basis of whether their name contains
  the string "BOUNDARY" or is one of the reserved boundary names.

  When a boundary (of Face objects) is created, the underpinning lower-dimensional (parent)
  region is removed from the region "Model".
  It is also moved from the unique to the non-unique region map if it was
  stored there originally.
  
  @remark builds boundaries only considering  unique, lower-dimensional regions (line-element boundaries will not be created).
  
  @return if all the created Face objects can be assigned to boundaries (as identified by the string BOUNDARY) the method returns true.
   If not, the remaining Face objects will be retained as a Boundary object named "Model_BOUNDARY and the method returns false.

  @attention by contrast to EstablishBoxBoundaries(),
  EstablishBoundaries() does not require to be applied to box-shaped models, but is
  more general.

  @test updated by SKM 2016
  */
  template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
  set<string>  BoundaryInterface<dim, BOUNDARY_COMPLEX>::EstablishBoundariesFromRegions()
  {
	  BOUNDARY_COMPLEX<dim>* model(static_cast<BOUNDARY_COMPLEX<dim>*>(this));
	  cout << "\nBoundaryInterface<" << dim << ">::EstablishBoundaries: searching for eligible boundary domains...\n";

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 1. compiling the unique regions that will be used as input for boundary creation
    //    and making a map of their elements that will be converted to faces
    // ---------------------------------------------------------------------
	  map<string,pair<size_t,size_t>>  eligibleRegions; // first and last element idx for input region
    vector<Element<dim>*>            elmts_to_become_faces;
    elmts_to_become_faces.reserve( model->Mesh().Elements() );
 
	  for ( typename map<string,csmp::Region<dim> >::iterator
		      it = model->UniqueRegionsBegin(); it != model->UniqueRegionsEnd(); ++it )
      {
        if ( !IsBoundaryName(it->first) )
          continue;
        if ( !hasLowerDimensionalRepresentation(it->second) )
          continue;
        // accumulating the elements of the eligible regions
        const pair<size_t,size_t> elmt_range{ elmts_to_become_faces.size(), elmts_to_become_faces.size() + (*it).second.Cells() };
        elmts_to_become_faces.insert( elmts_to_become_faces.end(),
                                     (*it).second.CellsBegin(), (*it).second.CellsEnd() );
                                     
        eligibleRegions.insert( make_pair( it->first, elmt_range ) );
      }

	  if ( eligibleRegions.empty() ) {
         csmp_error.Note( WARNING, "BoundaryInterface::EstablishBoundariesFromRegions",
                                     "unable to find eligible lower-dimensional regions to create Boundary objects from");
         return set<string>{};
      }


    // 2. replacing the elements by Faces (input elements are deleted and nullptrs returned)
    // -------------------------------------------------------------------------------------
    assert( connectivityCheck<dim>( elmts_to_become_faces.begin(), elmts_to_become_faces.end() ) == 0 );
    vector<Face<dim>*> faces = model->Mesh().ReplaceBoundaryElementsByFaces( model->Database(),
                                                                             elmts_to_become_faces.begin(),
                                                                             elmts_to_become_faces.end() );
    // 3. creating the Boundaries from the faces
    // -----------------------------------------
    typename vector<Face<dim>*>::iterator fit{ faces.begin() };
    set<string>  boundaries_created;

    for ( auto& it : eligibleRegions ) {
         // TODO: boundary names may have to be adjusted to meet CSMP conventions
         BOX_BOUNDARY boundary_flag = parseBoundary( it.first );
         if ( boundary_flag == MULTIPLE ) boundary_flag = IRREGULAR; // outside
         if ( AddBoundary( it.first.c_str(), next(fit,it.second.first), next(fit,it.second.second), boundary_flag ) )
           boundaries_created.insert( it.first );
      }
    
    if ( !boundaries_created.empty() ) {
         cout << "\n\nBoundaryInterface::EstablishBoundariesFromRegions: successfully created the external boundaries:\n\t";
         for ( auto bit : boundaries_created )
           cout <<" "<< bit;
         cout << endl;
      }

	  // 4. removing the original regions from which the boundaries were created
    // ------------------------------------------------------------------------------------
    // (no flagging for rebuilt of regions is necessary as they will be completely removed)
    const bool erase_elements{ false }; // this was already done above
    for ( auto& it : eligibleRegions ) {
         model->RemoveRegion( it.first.c_str(), erase_elements );
      }
    model->UpdateRegions();

	  cout << "\n\nBoundaryInterface::EstablishBoundariesFromRegions: done!\n";
    // if there are some unattributed faces left the method returs false
	  return boundaries_created;

} // end EstablishBoundariesFromRegions



 // DEBUGGING - checked that there are no duplicates or nullptrs in the 'elmts_to_become_faces' vector
// checking input vector for duplicates (OK for prism_test
/*
sort( elmts_to_become_faces.begin(), elmts_to_become_faces.end() );
bool hasDuplicates = adjacent_find( elmts_to_become_faces.begin(), elmts_to_become_faces.end()) !=
                                                                                 elmts_to_become_faces.end();
// contains null pointers ?
bool hasNullPointer = find( elmts_to_become_faces.begin(), elmts_to_become_faces.end(), nullptr) !=
                                                                         elmts_to_become_faces.end();
*/
  
  
// DEBUGGING - checked that there are no duplicates or nullptrs in the 'elmts_to_become_faces' vector
//           - no elements numbered as follows are among the memory violations:  for ( auto& it : elmts_to_become_faces ) it->Idx(999);
// - no issues uo to here here
/*
 {
    csmp::Region<dim>& fracture_domain(model->Region("FRACTURE"));
    for ( auto it=fracture_domain.CellsBegin(); it!=fracture_domain.CellsEnd(); ++it ) {
          assert( (*it)->Idx() >= 0 );
          for ( auto i{0U}; i<(*it)->Neighbors(); ++i )
            if ( (*it)->Neighbor(i) != nullptr )
              cerr << (*it)->Neighbor(i)->Idx() <<" ";
      }
 }
// are there shared elements or neighbors? - NO!
csmp::Region<dim>& frac_domain(model->Region("FRACTURE"));
for ( auto& it : elmts_to_become_faces ) {
     if ( frac_domain.Contains( it ) )
       cerr <<"region overlap at "<< it->Idx();
       for ( auto i{0U}; i<it->Neighbors(); ++i )
         if ( it->Neighbor(i) != nullptr )
           if ( frac_domain.Contains( it->Neighbor(i) ) )
             cerr <<"neighbor overlap at "<< it->Neighbor(i)->Idx();
  }
*/




/**
    Tries to replace bondary surface elements with Faces and assign these to Box boundaries.
    TOP, BOTTOM, INTERNAL, IRREGULAR, VERTICAL_SIDE etc. Edge boundaries are not created.
    
    @attention method assumes that the  perimeter of region "Model" correctly identifies the outside faces of the elements of the model.
    @attention this method was designed primarily for three-dimensional models.
    
    @author SKM
    @date 21/8/2018
*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
bool BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishBoxBoundariesFromOrientation()
  {

    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    MeshManager<dim>& mesh = boundaryComplex->Mesh();
    const size_t n_initial_faces = mesh.Faces();

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    const LocalVariables lvsFaces( boundaryComplex->Database().LocalVariablesAt(FACE) );
    const IntegrationPointVariables lvsIntegrationPoints( boundaryComplex->Database().IntegrationPointVariablesAt(FACE_INTEGRATION_POINT) );

    // creating faces on the outside of the model
    const csmp::Region<dim>& model_domain( boundaryComplex->Region("Model") );
    assert( model_domain.PerimeterCells() > 0 );
    for ( size_t i=model_domain.InteriorCells(); i<model_domain.Cells(); ++i )
      for ( auto j{0U}; j < model_domain.PerimeterFaces(i); ++j ) {
           // ascertaining that we are indeed at the model boundary
           assert( model_domain.E(i)->Neighbor( model_domain.PerimeterFace(i,j) ) == nullptr );
           // creating the boundary face
           mesh.AddBoundaryFace( model_domain.E(i), model_domain.PerimeterFace(i,j), lvsFaces, lvsIntegrationPoints );
        }

    vector<string> eligibleRegions;
    eligibleRegions.reserve( boundaryComplex->UniqueRegions() );
    
    // 1. Grouping pointers to Face objects of 'Model' boundary according to their facing direction
    // -------------------------------------------------------------------------------------------------
    vector<double>  nrml, nrml_right, nrml_left, nrml_top, nrml_bottom, nrml_front, nrml_back;
    Box             box;
    double          minLength(0.71); // dot-product of 2 unit vectors at an angle >=45 degrees
   
    box.UnitNormalTo( BOTTOM, dim, nrml_bottom );
    box.UnitNormalTo( TOP,    dim, nrml_top );
    box.UnitNormalTo( LEFT,   dim, nrml_left );
    box.UnitNormalTo( RIGHT,  dim, nrml_right );
    box.UnitNormalTo( FRONT,  dim, nrml_front );
    box.UnitNormalTo( BACK,   dim, nrml_back );
    
    set<Face<dim>*>  top_faces, bottom_faces, left_faces, right_faces, front_faces, back_faces, irregular_faces;
    
    // for all Face objects on the model boundary (initial faces would normally be zero)
    for ( auto fit=next(mesh.FacesBegin(),static_cast<long>(n_initial_faces)); fit!=mesh.FacesEnd(); ++fit )
      {
         // getting the (outward) pointing unit normal to face
         (*fit).UnitNormal( nrml );
         // classifying the faces in terms of their facing direction
         // (projecting: perfect alignment would give dot-product equal 1, inclinations up to 37 degrees cos(37)~0.8 are tolerated)
         if      ( dotProduct<dim>(nrml,nrml_bottom) >= minLength ) bottom_faces.insert(&(*fit)); // BOTTOM
         else if ( dotProduct<dim>(nrml,nrml_top)    >= minLength ) top_faces.insert(&(*fit));    // TOP
         else if ( dotProduct<dim>(nrml,nrml_left)   >= minLength ) left_faces.insert(&(*fit));   // LEFT
         else if ( dotProduct<dim>(nrml,nrml_right)  >= minLength ) right_faces.insert(&(*fit));  // RIGHT
         else if ( dotProduct<dim>(nrml,nrml_front)  >= minLength ) front_faces.insert(&(*fit));  // FRONT
         else if ( dotProduct<dim>(nrml,nrml_back)   >= minLength ) back_faces.insert(&(*fit));   // BACK
         // deal with the remaining cases, distinguishing sides etc.
         else irregular_faces.insert(&(*fit));
      } // end perimeter faces
    
    // 2. Creating boundaries from the non-empty sets of faces
    // -------------------------------------------------------------------------------------------------
    // BOTTOM
    if ( !bottom_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( bottom_faces.begin(), bottom_faces.end() );
         boundaryComplex->AddBoundary( "BOTTOM", boundary_faces.begin(), boundary_faces.end(), BOTTOM );
      } // TOP
    if ( !top_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( top_faces.begin(), top_faces.end() );
         boundaryComplex->AddBoundary( "TOP", boundary_faces.begin(), boundary_faces.end(), TOP );
      } // LEFT
    if ( !left_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( left_faces.begin(), left_faces.end() );
         boundaryComplex->AddBoundary( "LEFT", boundary_faces.begin(), boundary_faces.end(), LEFT );
      } // RIGHT
    if ( !right_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( right_faces.begin(), right_faces.end() );
         boundaryComplex->AddBoundary( "RIGHT", boundary_faces.begin(), boundary_faces.end(), RIGHT );
      } // FRONT
    if ( !front_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( front_faces.begin(), front_faces.end() );
         boundaryComplex->AddBoundary( "FRONT", boundary_faces.begin(), boundary_faces.end(), FRONT );
      } // BACK
    if ( !back_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( back_faces.begin(), back_faces.end() );
         boundaryComplex->AddBoundary( "BACK", boundary_faces.begin(), boundary_faces.end(), BACK );
      } // IRREGULAR (left-over faces)
    if ( !irregular_faces.empty() ) {
         vector<Face<dim>*> boundary_faces( irregular_faces.begin(), irregular_faces.end() );
         boundaryComplex->AddBoundary( "IRREGULAR", boundary_faces.begin(), boundary_faces.end(), IRREGULAR );
      }

    
    // 4. checking whether faces remain that could not be assigned
    // -------------------------------------------------------------------------------------------------
    size_t n_faces_assigned = Boundary("BOTTOM").Cells() + Boundary("RIGHT").Cells() +
                              Boundary("TOP").Cells() + Boundary("LEFT").Cells();
    if ( !irregular_faces.empty() ) n_faces_assigned += Boundary("IRREGULAR").Cells() ;
    
    if constexpr ( dim == 3 )
      n_faces_assigned += Boundary("FRONT").Cells() + Boundary("BACK").Cells();
      
      if ( n_faces_assigned != mesh.Faces() - n_initial_faces ) {
        csmp_error.Note( WARNING, "BoundaryInterFace::EstablishBoundaryFlagsFromOrientation",
                           "Not all boundary Face faces could be assigned to standard boundaries");
         return false;
      }
    
    
    // 5. adjusting node flags for TOP or BOTTOM boundary nodes that were flagged as irregular
    // ---------------------------------------------------------------------------------------
    for ( auto it = boundaryComplex->BoundariesBegin(); it != boundaryComplex->BoundariesEnd(); ++it )
      if ( (*it).first == "TOP" or (*it).first == "BOTTOM" )
        {
           for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.NodesEnd(); ++nit )
             if ( (*nit)->AtBoundary() == IRREGULAR )
               (*nit)->AtBoundary( parseBoundary( (*it).first ) );
        }

    cout << "\n\nBoundaryInterface::EstablishBoxBoundariesFromOrientation: done!\n";
    return true;
    
} // end EstablishBoundaryFlagsFromOrientation




/** Helper function for method below
     
    Creates desired range of faces and interconnects them with one-another.
*/
template<uint32_t dim>
void createBoundaryFaces( MeshManager<dim>& mesh, const PropertyDatabase<dim>& dbase,
                          const set<pair<Element<dim>*,uint32_t> >& face_set, vector<Face<dim>*>& boundary_faces )
 {
    assert( !face_set.empty() );
    
    boundary_faces.clear();
    boundary_faces.reserve( face_set.size() );
    const LocalVariables             lvars( dbase.LocalVariablesAt(FACE) );
    const IntegrationPointVariables  ivars(dbase.IntegrationPointVariablesAt(FACE) );

    // creating the faces
    for ( auto& it : face_set ) {
         // building the face info structure
         boundary_faces.push_back( mesh.AddBoundaryFace( it.first, it.second, lvars, ivars ) );
      }
    // connecting the faces
    mesh.template BuildConnectivity<Face>( boundary_faces.begin(), boundary_faces.end() );

} // end createBoundaryFaces

//template void createBoundaryFaces( MeshManager<3>&, const PropertyDatabase<3>&, const set<pair<Element<3>*,size_t> >&, vector<Face<3>*>& );
 

/**
 
Creates BOX boundaries using the node flags to identify sides, edges, and corners; use for simple models where corresponding lines or surfaces are missing.

The result are Boundary objects with standard names such as TOP, BOTTOM, INTERNAL, IRREGULAR, EDGE1, CNR2 etc.

@attention method assumes that boundary flags were  assigned correctly to the nodes, else a method which does this is called.

Procedure
- using the node flags of the perimeter elements faces, boundary by boundary VData are created that are subsequently employed to create the Face objects of the new boundaries.
Then these are created.

*/
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim,BOUNDARY_COMPLEX>::EstablishBoxBoundariesFromNodeFlags( bool recreate_box_boundary_flags_before )
 {
    BOUNDARY_COMPLEX<dim>* boundaryComplex( static_cast<BOUNDARY_COMPLEX<dim>*>(this) );
    if ( recreate_box_boundary_flags_before ) recreateBoxBoundaryFlags( *boundaryComplex );
    
 
  cout << "\nBoundaryInterface<" << dim << ">::EstablishBoxBoundariesFromNodeFlags: searching for eligible boundary domains...\n";

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !boundaryMap_.empty() ) {
       csmp_error.Note( WARNING, "BoundaryInterface::EstablishBoxBoundariesFromNodeFlags:",
                                   "model already contains Boundary objects; nothing was done.");
       return;
    }

  const csmp::Region<dim>& model_domain = boundaryComplex->Region("Model");

  set<pair<Element<dim>*,uint32_t> > top_faces, bottom_faces, left_faces, right_faces, front_faces, back_faces, irregular_faces, internal_faces;

  // for all element faces on the model boundary
  for ( size_t eid{model_domain.InteriorCells()}; eid < model_domain.Cells(); ++eid )
    for ( auto j{0U}; j < model_domain.PerimeterFaces(eid); ++j ) {
         const auto face_id{ model_domain.PerimeterFace(eid,j) };
         // getting the boundary flag of the face
         const BOX_BOUNDARY bflag = atBoundary( model_domain.E(eid), face_id );
         // storing the FaceInfo for the correct boundary
         switch ( bflag ) {
              case BOTTOM:
                  bottom_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case RIGHT:
                  right_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case TOP:
                  top_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case LEFT:
                  left_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case IRREGULAR:
                  irregular_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case INTERNAL:
                  internal_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case BACK:
                  back_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              case FRONT:
                  front_faces.insert( make_pair( model_domain.E(eid), face_id ) );
                break;
              default:
                csmp_error.Note( ERROR, "BoundaryInterface::EstablishBoxBoundariesFromNodeFlags:"
                                  "boundary flag not recognised as valid option:", parseBoundary(bflag) );
           }
      }
 
 
 
  // 2. Creating boundaries from the non-empty sets of faces
  // -------------------------------------------------------------------------------------------------
  vector<Face<dim>*> boundary_faces;
  // BOTTOM
  if ( !bottom_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), bottom_faces, boundary_faces );
      boundaryComplex->AddBoundary( "BOTTOM", boundary_faces.begin(), boundary_faces.end(), BOTTOM );
    } // TOP
  if ( !top_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), top_faces, boundary_faces );
      boundaryComplex->AddBoundary( "TOP", boundary_faces.begin(), boundary_faces.end(), TOP );
    } // LEFT
  if ( !left_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), left_faces, boundary_faces );
      boundaryComplex->AddBoundary( "LEFT", boundary_faces.begin(), boundary_faces.end(), LEFT );
    } // RIGHT
  if ( !right_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), right_faces, boundary_faces );
      boundaryComplex->AddBoundary( "RIGHT", boundary_faces.begin(), boundary_faces.end(), RIGHT );
    } // FRONT
  if ( !front_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), front_faces, boundary_faces );
      boundaryComplex->AddBoundary( "FRONT", boundary_faces.begin(), boundary_faces.end(), FRONT );
    } // BACK
  if ( !back_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), back_faces, boundary_faces );
      boundaryComplex->AddBoundary( "BACK", boundary_faces.begin(), boundary_faces.end(), BACK );
    } // IRREGULAR
  if ( !irregular_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), irregular_faces, boundary_faces );
      boundaryComplex->AddBoundary( "IRREGULAR", boundary_faces.begin(), boundary_faces.end(), IRREGULAR );
    } // INTERNAL
  if ( !internal_faces.empty() ) {
      createBoundaryFaces( boundaryComplex->Mesh(), boundaryComplex->Database(), internal_faces, boundary_faces );
      boundaryComplex->AddBoundary( "INTERNAL", boundary_faces.begin(), boundary_faces.end(), INTERNAL );
    }


  // 3. Creating lower-dimensional edge faces
  // -------------------------------------------------------------------------------------------------
  EstablishEdgeBoundariesOfBoxShapedModel();

  cout << "\n\nBoundaryInterface::EstablishBoxBoundariesFromFlags: done!\n";
  
} // end EstablishBoxBoundariesFromNodeFlags







template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
void BoundaryInterface<dim, BOUNDARY_COMPLEX>::BoundariesOut() const
 {
     cout <<"\nBoundaryInterface<"<< dim <<",Boundary<Face>>::BoundariesOut: ";
     if ( Boundaries() == 0U ) {
          cout <<"\tmodel does not contain any boundaries.\n\n";
          return;
       }
     cout <<"boundaries of ";
     if ( BoxShaped() ) cout <<"box-shaped model:\n";
     else cout <<"irregularly-shaped model:\n";
     for ( auto bit=BoundariesBegin(); bit!=BoundariesEnd(); ++bit ) {
          cout <<"\n\t"<< (*bit).first <<", box-flag: "<< parseBoundary( (*bit).second.AtBoundary() );
          cout <<" "<< (*bit).second.Cells() <<" faces, ";
          // in 3D a boudary is a surface
           if constexpr ( dim == 3U ) {
                cout <<"area (m2): "<< (*bit).second.Area();
                cout <<", perimeter length (m): "<< (*bit).second.Perimeter();
             }
           if constexpr ( dim == 2U )
             cout <<" length (m): "<< (*bit).second.Area();
       }
     cout << endl << endl;
     cout.flush();
 }




  template class BoundaryInterface<1U, Model>;
  template class BoundaryInterface<2U, Model>;
  template class BoundaryInterface<3U, Model>;

} // csmp

