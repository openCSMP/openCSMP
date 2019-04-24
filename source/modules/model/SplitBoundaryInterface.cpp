#include "SplitBoundaryInterface.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

#define SPLIT_BOUNDARY_DEBUG

namespace csmp {

// JC: check wehether default constructors/destructors are actually neeeded?
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundaryInterface()
{}

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundaryInterface( const SplitBoundaryInterface& bd )
  : splitBoundaryMap_( bd.splitBoundaryMap_ )
{}

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::~SplitBoundaryInterface()
{}


/**
Creates a name like SPLITBOUNDARY_region1_region2 adding a number if this is necessary to make the
name unique.

@attention the facing relationships of the boundary are preserved so that region1 is the first in the
argument pair.

@author SKM
@date January 2018
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
std::string SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryName( const std::pair<std::string, std::string>& juxtaposed_regions ) const
{
  string split_boundary_name( "SPLITBOUNDARY_" + juxtaposed_regions.first + '_' + juxtaposed_regions.second );

  // if the substring set is empty
  if ( ContainsSplitBoundary( split_boundary_name.c_str() ) ) {
    ErrorHandler::Instance().notice( WARNING, "SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateSplitBoundaryName:",
                                     split_boundary_name.c_str(), "already exists, try other name.'\0'." );
    return std::string( "\0" );
  }

  // making a set of boundary names
  /*
  const size_t substrings_used_in_search(intersected_regions.size());
  for ( auto it=SplitBoundariesBegin(); it!=SplitBoundariesEnd(); ++it ) {
  size_t substrings_found(0U);
  for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
  // if the substring is found
  if ( (*it).first.find(*ir) !=std::string::npos ) substrings_found++;
  // when all substrings are contained in the boundary name, it is returned
  if ( substrings_found == substrings_used_in_search )
  return (*it).first;
  }
  return std::string("\0");
  */
  return split_boundary_name;

} // end CreateSplitBoundaryName


/**
tag SPLITBOUNDARY should be at the beginning of the region name
DON'T confuse with BOUNDARY

@todo remove if this is not used anywhere
*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, BOUNDARY_COMPLEX>::IsEligibleSplitBoundaryRegionName( const std::string& spbNname ) const
{
  std::size_t found_position = spbNname.find( "SPLITBOUNDARY" );
  if ( found_position == 0 && found_position != std::string::npos )
    return true;
  return false;
}


/**
returns reference to SplitBoundary
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundary<dim>&  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundary( const std::string& spbName )
{
  splitBoundaryIterator  sbit = splitBoundaryMap_.find( spbName );
  if ( sbit != splitBoundaryMap_.end() )
    return (*sbit).second;
  else {
    std::string errMsg( "SplitBoundary does not exist!" );
    errMsg.append( " (" + spbName + ")" );
    throw csmp::Exception( ERROR,
                           "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundary",
                           errMsg.c_str() );
  }
}


/**
returns const reference to SplitBoundary
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
const SplitBoundary<dim>&  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundary( const std::string& spbName ) const
{
  splitBoundaryConstIterator  sbit = splitBoundaryMap_.find( spbName );
  if ( sbit != splitBoundaryMap_.end() )
    return (*sbit).second;
  else {
    std::string errMsg( "SplitBoundary does not exist!" );
    errMsg.append( " (" + spbName + ")" );
    throw csmp::Exception( ERROR,
                           "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundary",
                           errMsg.c_str() );
  }
}


/**
Reports whether a SplitBoundary with the corresponding name exists inside the model.
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::ContainsSplitBoundary( const std::string& bname ) const
{
  if ( splitBoundaryMap_.find( bname ) != splitBoundaryMap_.end() ) return true;
  return false;
}


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin()
{ return splitBoundaryMap_.begin(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd()
{ return splitBoundaryMap_.end(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin() const
{ return splitBoundaryMap_.begin(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd() const
{ return splitBoundaryMap_.end(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
size_t  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundaries() const
{ return splitBoundaryMap_.size(); }


/**
Removes the SplitBoundary object with the given name if it exists inside the model; else returns
with a warning.

User has the option to trigger the MeshManager to also delete the corresponding InterFace elements
from the model.
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary,
                                                                                bool deleteElements )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsSplitBoundary( splitboundary.Name().c_str() ) ) {
    csmp_error.notice( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary",
                       splitboundary.Name().c_str(), "Does not exist; nothing was done." );
    return;
  }

  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>* >(this) );

  if ( deleteElements )
  {
    // Faces removed from all existing boundaries and erased from MeshManager
    std::vector<InterFace<dim>* > elementsToDelete( splitboundary.SimplexVector().begin(), splitboundary.SimplexVector().end() );

    if ( !elementsToDelete.empty() )
    {
      // Detach InterFaces from Neighbors
      splitboundary.DetachElementsFromNeighbors();

      // remove redundant InterFaces from existing SplitBoundaries
      for ( splitBoundaryIterator
            sbit = splitBoundaryMap_.begin();
            sbit != splitBoundaryMap_.end();) {

        if ( removeVectorElements( (*sbit).second.SimplexVector(), elementsToDelete ) > 0 )
        {
          if ( !(*sbit).second.SimplexVector().empty() )
          {
            (*sbit).second.CreateNodePointerVector();
            (*sbit).second.IdentifyPerimeter();
            ++sbit;
          }
          else
            splitBoundaryMap_.erase( sbit++ );
        }
        else
        {
          // we need to reestablish the perimeter because
          // the order of SimplexVector() was changed by removeVectorElements ()
          (*sbit).second.IdentifyPerimeter();
          ++sbit;
        }
      }

      // update all indices
      splitboundaryComplex->UpdateIndices();

      return;
    }

  } // deleting elements

    // finding the boundary
  auto iterSplitBoundary( splitBoundaryMap_.end() );
  for ( splitBoundaryIterator it = splitBoundaryMap_.begin(); it != splitBoundaryMap_.end(); ++it )
    if ( &it->second == &splitboundary )
      iterSplitBoundary = it;

  if ( iterSplitBoundary == splitBoundaryMap_.end() )
    throw csmp::Exception( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary", "splitboundary does not exist" );

  // if the boundary was found in the list, it is erased
  splitBoundaryMap_.erase( iterSplitBoundary );

} // end RemoveSplitBoundary


// -----------------------------------------------
// Binary input/output
// -----------------------------------------------
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::OutputSplitBoundariesToBinary( const char* file_name ) const
{
  if ( this->SplitBoundaries() == 0 ) {
    std::cout << "\nSplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary: There are no split boundaries. Nothing to write.'";
    return false;
  }

  std::string bin_file( file_name );
  std::string heading( "Binary splitboundary list for CSMP model, file '" );
  heading += bin_file;
  heading += "'.";

  fstream fp( bin_file.c_str(), ios::out | ios::binary );
  if ( !fp.is_open() ) {
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary:", "Binary file could not be created." );
    return false;
  }
  // writing the file header
  skm_C_fwrite( fp, heading.c_str() );

  // writing number of splitboundaries
  size_t records( this->SplitBoundaries() );
  fp.write( (char*)&records, sizeof( size_t ) );

  for ( auto bit( this->SplitBoundariesBegin() ); bit != this->SplitBoundariesEnd(); bit++ )
  {
    // writing name of splitboundary
    skm_C_fwrite( fp, bit->first.c_str() );
    // splitboundary
    if ( !bit->second.Out( fp ) )
      throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary:", "SplitBoundary could not be stored." );
  }

  // cleaning up
  fp.close();
  std::cout << "\nSplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary: Split boundaries have been successfully written to: '";
  std::cout << bin_file << "'" << std::endl;
  return true;

} // end OutputSplitBoundariesToBinary


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InputSplitBoundariesFromBinary( const char* file_name )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  std::string bin_file( file_name );
  fstream fp( bin_file.c_str(), ios::in | ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.notice( ERROR, "SplitBoundaryInterface::InputSplitBoundariesFromBinary:",
                       bin_file, "file could not be opened; nothing was done." );
    return false;
  }

  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  const csmp::Region<dim>&  rref( splitboundaryComplex->Region( "Model" ) );

  // reading the file header
  char csCache[255];
  skm_C_fread( fp, csCache );

  // reading number of splitboundaries
  size_t records( 0 );
  fp.read( (char*)&records, sizeof( size_t ) );

  std::cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: reading " << csCache << " containing "
    << records << " splitboundaries\n" << std::endl;

  for ( size_t i( 0 ); i < records; ++i )
  {
    std::string bName;
    // reading name of splitboundary
    skm_C_fread( fp, csCache );
    bName = csCache;
    // inserting splitboundary if not existing yet
    std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
      bit = splitBoundaryMap_.insert( std::make_pair( bName, csmp::SplitBoundary<dim>( bName, splitboundaryComplex->Database() ) ) );
    // splitboundary
    if ( !bit.second )
      return false;
    if ( !bit.first->second.In( splitboundaryComplex->Mesh(), splitboundaryComplex->FE_Manager(), rref, fp ) )
      return false;
    std::cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: read boundary " << bName << " successfully.\n";
  }
  fp.close();
  return true;

} // end InputSplitBoundariesFromBinary


/// container of juxtaposed element pairs for SplitBoundary creation:
template<size_t dim>
struct SplitBoundaryElementSets : public
  //       key consisting out the names of regions juxtaposed at the SplitBoundary
  std::map<std::pair<std::string, std::string>,
  // set that stores pointers to the pairs of elements juxtaposed across boundary
  // size_t parameter gives the local number of the element face that sits at the split boundary
  std::set<std::pair<std::pair<Element<dim>*, size_t>,
  std::pair<Element<dim>*, size_t> > > > {
};


/**
Creates SplitBoundaries detecting, them in the (ANSYS) input model
as node-matched interfaces, connecting such disconnected perimeter element faces
in mesh; these are detected and grouped by bordering regions and turned into SplitBoundary objects
with names following the same conventions as for Boundary objects.

@author SKM 20/03/2018
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries()
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&               csmp_error( ErrorHandler::Instance() );

  // 1. detecting potential SplitBoundaries
  // --------------------------------------
  // reports local ids of shared face on either side of the split boundary
  set<pair<pair<Element<dim>*, size_t>, pair<Element<dim>*, size_t> > >  interface_elmt_pairs;
  if ( !findSplitInterfaceElements( splitboundaryComplex->Region( "Model" ), interface_elmt_pairs ) ) {
    csmp_error.notice( WARNING, "SplitBoundaryInterface::DetectAndCreateSplitBoundaries:",
                       "node-coordinate matched faces / internal boundaries could not be detected; nothing was done." );
    return false;
  }

  // 2. classifying the detected interfaces in terms of the regions that they juxtapose
  // ----------------------------------------------------------------------------------
  // 2.1 assigning integer keys to the elements of the model regions so that their name can be identified from the perimeter elements
  //     names by their index integer
  vector<string>  region_names;
  if ( !splitboundaryComplex->Database().IsDefined( "region number" ) ) {
    splitboundaryComplex->CreateProperty( "region number", "-", SCALAR, ELEMENT );
    splitboundaryComplex->UpdateIndices();
  }
  splitboundaryComplex->CountAndLabelRegions( "region number", region_names );

  // 2.2 grouping interface element pairs into ones that juxtapose specific regions against one another
  //     these will later become specific split boundaries
  typedef set<pair<pair<Element<dim>*, size_t>, pair<Element<dim>*, size_t> > > INTERFACE_ELEMENT_PAIRS;
  map<pair<string, string>, INTERFACE_ELEMENT_PAIRS>  split_boundary_map;
  const csmp::Index reg_key( splitboundaryComplex->Database().StorageKey( "region number" ) );

  // for all the SplitBoundary objects supplied as sets of pairs of Element pointers and interface idx values
  for ( auto iit : interface_elmt_pairs ) {
    // extracting region names from the name-integer vector
    pair<string, string> key = make_pair( region_names[static_cast<long>(iit.first.first->Read( reg_key ))],
                                          region_names[static_cast<long>(iit.second.first->Read( reg_key ))] );
    // storing interfaces in split boundary maps
    auto eit = split_boundary_map.insert( make_pair( key, INTERFACE_ELEMENT_PAIRS( { iit } ) ) );
    // if no insertion could be performed, the element pair is added to an existing set
    if ( !eit.second ) (*eit.first).second.insert( iit );
  }

  // echoing the map to the screen
#ifdef SPLIT_BOUNDARY_DEBUG
  cerr << "\nSplitBoundaryInterface::DetectAndCreateSplitBoundaries: interface region pairs found:\n";
  for ( auto i : split_boundary_map )
    cerr << i.first.first << "," << i.first.second << "\n";
  cerr << endl;
#endif

  // 3. Creating splitboundaries for each of the discovered juxtapositions of regions
  // --------------------------------------------------------------------------------
  for ( auto it : split_boundary_map ) {
    // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
    string bname = CreateSplitBoundaryName( it.first );

    // extracting the split boundary neighbor elements into interface element pairs
    InterFaceSet<dim>  ifset( it.second );

    // asking the MeshManager to create the required number of InterFace objects
    // and creating the SplitBoundary
    std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
      bit = splitBoundaryMap_.insert( std::make_pair( bname,
                                      csmp::SplitBoundary<dim>( bname, splitboundaryComplex->Database(),
                                      splitboundaryComplex->FE_Manager(),
                                      splitboundaryComplex->Mesh(), ifset ) ) );
    if ( bit.second ) {
      std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries: ";
      std::cout << " boundary created successfully.\n";
      splitboundaryComplex->UpdateIndices();
    }
    else
      throw csmp::Exception( INFO,
                             "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries:",
                             bname.c_str(),
                             "boundary already exists. Nothing was done." );
  }

  return true;

} // end DetectAndCreateSplitBoundaries


/**
Creates SplitBoundarty around a unique volumetric(3D) / surface(2D) region; except for model boundary.

This will involve the creation and connection of InterFace elements by the MeshManager.
This will only work for unique Regions which are not overlapping.

@attention This method makes no sense for the region 'Model' as it encompasses all unique regions.

@attention:  SplitBoundaryInterface objects cannot be created in 1D models or between regions which contain
one-dimensional elements.

Option deleteRegionAndItsElements in fact provides a choice for keeping fracture\fault or any other low dimensional Region
or use such Region just to introduce the interface around it and delete it afterwards
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const string& region )
{
  bool succeeded( false );
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( region == "Model" ) {
    csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                       "Region 'Model' not eligible for InsertSplitBoundary." );
    return false;
  }

  if ( !splitboundaryComplex->IsUnique( region.c_str() ) )
    csmp_error.notice( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                       "This method is intended for the creation of splitboundaries around unique Regions" );

  std::string  splitboundaryName( region );
  splitboundaryComplex->InsertBoundary( IRREGULAR, region.c_str() );
  Boundary<dim>& boundary( splitboundaryComplex->Boundary( splitboundaryName ) );

  // attempt to create a regular (InterFace-based) splitboundary
  std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
    it = splitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName, splitboundaryComplex->Database() ) ) );
  if ( it.second )
  {
    std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary around: " << region << std::endl;
    splitboundaryComplex->UpdateIndices();
    succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );

  }
  else
    throw csmp::Exception( INFO, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                           splitboundaryName.c_str(),
                           "boundary already exists. Nothing was done." );

  // remove temporal boundary
  splitboundaryComplex->RemoveBoundary( boundary );


  // reestablish containers
  // unique regions
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity();
    rit->second.IdentifyPerimeter();
  }
  // non-unique regions
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity();
    rit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
    bit->second.CreateNodePointerVector();
    bit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
    bit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
    sbit->second.CreateNodePointerVector();
    sbit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
    sbit->second.IdentifyPerimeter();
  }

  // update indexes
  splitboundaryComplex->UpdateIndices();

  if ( succeeded )
    std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary created splitboundary around " << region << std::endl;
  else
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

  return true;

} // end InsertSplitBoundary


/**
Method forms a SplitBoundary between the two supplied regions. This will involve the creation
and connection of InterFaces.

@attention:  BoundariesInterface cannot be created in 1D models or between regions which contain
one-dimensional elements.

This method will only work for Regions which are not overlapping.

This method makes no sense for the region 'Model' as it encompasses all unique regions.
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const std::string& group1, const std::string& group2, bool createRegionBetween )
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( group1 == "Model" or group2 == "Model" ) {
    csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Region 'Model' not eligible for InsertSplitBoundary." );
    return false;
  }
  if ( group1 == group2 ) {
    csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Provided Regions are identical." );
    return false;
  }


  if ( !splitboundaryComplex->IsUnique( group1.c_str() ) or !splitboundaryComplex->IsUnique( group2.c_str() ) ) {
    csmp_error.notice( WARNING,
                       "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                       "This method is intended for the creation of splitboundaries between unique Regions" );

    // checking for a potential overlap of the regions, if the regions are non-unique
    if ( splitboundaryComplex->RegionIntersection( group1.c_str(), group2.c_str(), "groupintersection" ) ) {
      csmp_error.notice( ERROR,
                         "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                         "one of the supplied regions is not unique and they overlap",
                         "It was therefore impossible to insert a boundary" );
      splitboundaryComplex->RemoveRegion( "groupintersection", false );
      return false;
    }
  }

  std::string  splitboundaryName( std::string( group1 ) + std::string( "_" ) + std::string( group2 ) );

  splitboundaryComplex->InsertBoundary( group1.c_str(), group2.c_str(), createRegionBetween );
  Boundary<dim>& boundary( splitboundaryComplex->Boundary( splitboundaryName ) );

  const csmp::Region<dim>&  gref1( splitboundaryComplex->Region( group1.c_str() ) );
  const csmp::Region<dim>&  gref2( splitboundaryComplex->Region( group2.c_str() ) );

  // checking whether the two regions share some nodes (these will mark their common boundary)
  const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
  bool succeeded( false );

  if ( shared_nodes == 0U )
    csmp_error.notice( WARNING,
                       "Model<dim,SPLITBOUNDARY_COMPLEX>::InsertBoundary",
                       "The regions of interest do not share any nodes; trying to create a split boundary" );


  // attempt to create a regular (InterFace-based) splitboundary
  else {

    std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
      it = splitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,
                                     splitboundaryComplex->Database() ) ) );
    if ( it.second )
    {
      std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary between " << group1 << " and " << group2 << std::endl;
      splitboundaryComplex->UpdateIndices();
      if ( createRegionBetween )
        splitboundaryComplex->RegionBetween( group1.c_str(), group2.c_str(), std::string( group1 + std::string( "_" ) + group2 ).c_str() );
      succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );
    }
    else
      throw csmp::Exception( INFO,
                             "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                             splitboundaryName.c_str(),
                             "boundary already exists. Nothing was done." );
  }

  // remove temporal boundary
  splitboundaryComplex->RemoveBoundary( boundary );

  // reestablish containers
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity();
    rit->second.IdentifyPerimeter();
  }
  // SKM FIX also the non-unique regions need to be rebuilt
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity();
    rit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
    bit->second.CreateNodePointerVector();
    bit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
    bit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
    sbit->second.CreateNodePointerVector();
    sbit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
    sbit->second.IdentifyPerimeter();
  }

  // update indexes
  splitboundaryComplex->UpdateIndices();

  if ( succeeded )
    std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary created splitboundary between " << group1 << " and " << group2 << std::endl;
  else
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

  return true;

} // end InsertSplitBoundary


/**
Prints current SplitBoundaries to screen
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesOut() const
{
  cout << "\nSplitBoundaryInterface::SplitBoundariesOut: current split boundaries in the model: " << splitBoundaryMap_.size();
  // std::map<std::string,csmp::SplitBoundary<dim> >  splitBoundaryMap_
  for ( auto it = splitBoundaryMap_.begin(); it != splitBoundaryMap_.end(); ++it ) {
    cout << "\n\n\tSplitBoundary: " << (*it).first << "\n";
    (*it).second.Out();
  }

} // end Out


template class SplitBoundaryInterface<1U, Model>;
template class SplitBoundaryInterface<2U, Model>;
template class SplitBoundaryInterface<3U, Model>;

} // csmp
