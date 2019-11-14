#include "RegionInterface.h"
#include "Model.h"
#include "UnionFind.h"

namespace csmp {

template<size_t dim, template<size_t> class REGION_COMPLEX>
RegionInterface<dim, REGION_COMPLEX>::RegionInterface()
{}

template<size_t dim, template<size_t> class REGION_COMPLEX>
RegionInterface<dim, REGION_COMPLEX>::RegionInterface( const RegionInterface& re )
  : uniqueGroupMap_( re.uniqueGroupMap_ ),
  groupMap_( re.groupMap_ )
{}

template<size_t dim, template<size_t> class REGION_COMPLEX>
RegionInterface<dim, REGION_COMPLEX>::~RegionInterface()
{}


template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::iterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsBegin()
{ return uniqueGroupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::iterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsEnd()
{ return uniqueGroupMap_.end(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::iterator  RegionInterface<dim, REGION_COMPLEX>::RegionsBegin()
{ return groupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::iterator  RegionInterface<dim, REGION_COMPLEX>::RegionsEnd()
{ return groupMap_.end(); }


template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::const_iterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsBegin() const
{ return uniqueGroupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::const_iterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsEnd() const
{ return uniqueGroupMap_.end(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::const_iterator  RegionInterface<dim, REGION_COMPLEX>::RegionsBegin() const
{ return groupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename std::map<std::string, csmp::Region<dim> >::const_iterator  RegionInterface<dim, REGION_COMPLEX>::RegionsEnd() const
{ return groupMap_.end(); }


template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::Regions() const
{ return uniqueGroupMap_.size() + groupMap_.size(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::UniqueRegions() const
{ return uniqueGroupMap_.size(); }


/**
Method first searches the Region in the unique region map, then in the
non-unique region map list where regions may overlap. If the desired region
can not be found, an ERROR csmp::Exception is thrown

In the case of failure, a reference to the Model Region may be returned.
This region is guaranteed to be
there always.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
const Region<dim>&  RegionInterface<dim, REGION_COMPLEX>::Region( const std::string& region_name ) const
{
  if ( region_name.empty() )
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region: ",
                           "regions search string is empty" );

  // first a look in the unique group list
  typename std::map<std::string, csmp::Region<dim> >::const_iterator  iter( uniqueGroupMap_.find( region_name ) );
  if ( iter != uniqueGroupMap_.end() )
    return (*iter).second;

  // now a look at the generic group list
  iter = groupMap_.find( std::string( region_name ) );

  if ( iter != groupMap_.end() )
    return (*iter).second;
  else
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region: ",
                           (std::string( "region does not exist: " ) + region_name) );

  return Region( "Model" );
}



/**
Method first searches the Region in the unique region map, then in the
non-unique region map list where regions may overlap. If the desired region
can not be found, an ERROR csmp::Exception is thrown

In the case of failure, a reference to the Model Region may be returned.
This region is guaranteed to be
there always.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
Region<dim>&  RegionInterface<dim, REGION_COMPLEX>::Region( const std::string& region_name )
{
  if ( region_name.empty() )
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region: ",
                           "regions search string is empty" );

  typename std::map<std::string, csmp::Region<dim> >::iterator  iter( uniqueGroupMap_.find( region_name ) );
  if ( iter != uniqueGroupMap_.end() )
    return (*iter).second;

  iter = groupMap_.find( region_name );

  if ( iter != groupMap_.end() )
    return (*iter).second;
  else
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region:",
                           (std::string( "region does not exist: " ) + region_name) );

  return Region( "Model" );
}


template<size_t dim, template<size_t> class REGION_COMPLEX>
bool  RegionInterface<dim, REGION_COMPLEX>::IsUnique( const std::string& region_name ) const
{
  typename std::map<std::string, csmp::Region<dim> >::const_iterator  iter( uniqueGroupMap_.find( region_name ) );
  if ( iter != uniqueGroupMap_.end() )
    return true;
  return false;
}



template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::ContainsRegion( const std::string& region_name ) const
{
  if ( uniqueGroupMap_.find( region_name ) != uniqueGroupMap_.end() )
    return true;
  if ( groupMap_.find( region_name ) != groupMap_.end() )
    return true;

  return false;
} // end ContainsRegion




/**
Forms non-unique user-defined region by graph traversal, relying only on node-to-parent element connections.

@return element and node numbers are compared with MeshManager entries to verify that
all elements and nodes were discovered. If so method returns true, else false

@param reestablishNeighborConnectivity will prompt CSMP to recreate element neighbor connectivity.

@author SKM
@date 5/4/2016
*/
/*
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim,REGION_COMPLEX>::CreateNonUniqueMasterRegionFromRootNode( bool reestablishNeighborConnectivity )
{
ErrorHandler&  csmp_error( ErrorHandler::Instance() );

// does this region already exist
if ( ContainsRegion(masterRegion_) ) {
csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::CreateNonUniqueMasterRegionFromRootNode:",
masterRegion_.c_str(), " already exists; nothing was done." );
return false;
}

std::pair<typename std::map<std::string,csmp::Region<dim> >::iterator,bool>
newRegion = groupMap_.insert( std::make_pair( masterRegion_, csmp::Region<dim>( masterRegion_,
static_cast<REGION_COMPLEX<dim>*>(this)->Database()) ) );

// if region could be inserted successfully
if ( newRegion.second ) {
const size_t elements = (*newRegion.first).second.AccumulateAll( &(static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().RootNode()),
reestablishNeighborConnectivity );
if ( elements < static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() ) {
std::cerr <<"\n\n\tdiscovered only "<< elements <<" versus "<< static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() <<" elements.\n\n";
csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::CreateNonUniqueMasterRegionFromRootNode:",
"mesh tree travel discovered less elements than there are in the model; is the mesh disconnected? - is there stand-alone mesh?");
}
}

else {
csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::CreateNonUniqueMasterRegionFromRootNode:",
masterRegion_, "region could not be formed.");
return false;
}

// checking that all elements and nodes were discovered
if ( (*newRegion.first).second.Elements() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() or
(*newRegion.first).second.Nodes() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Nodes() )
return false;

return true;
}
*/

/**
Forms contiguous multiple domains by graph traversal of all reachable elements in the mesh without expecting element-to-neighbor connections
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::CreateRegions( bool is_unique, bool reestablishNeighborConnectivity )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) );
  auto& meshMgr = regionComplex->Mesh();

  bool ret = true;
  for ( size_t i = 0; i < meshMgr.NodeGroups(); i++ ) {
    std::string regionname = "Model_" + std::to_string( i );
    if ( i == 0 ) regionname = "Model"; // first region name is 'Model', ann then Model_1, Model_2 and so on.

                                        // does this region already exist
    if ( ContainsRegion( regionname ) ) {
      RemoveRegion( regionname.c_str(), false );
    }

    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      newRegion = (is_unique) ?
      uniqueGroupMap_.insert( std::make_pair( regionname, csmp::Region<dim>( regionname,
                              static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
      :
      groupMap_.insert( std::make_pair( regionname, csmp::Region<dim>( regionname,
                        static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

    // if region was inserted successfully
    if ( newRegion.second )
      (*newRegion.first).second.AccumulateAll( meshMgr.RootNode( i ), reestablishNeighborConnectivity );

    else {
      csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::CreateRegions:",
                         regionname, "region could not be formed." );
      ret = false;
    }
  }

  return ret;
}


template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::CreateRegionFromRootNode( const char* regionname, bool is_unique, bool reestablishNeighborConnectivity )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  // does this region already exist
  if ( ContainsRegion( regionname ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::CreateRegionFromRootNode:",
                       regionname, "region already exists; nothing was done." );
    return false;
  }

  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
    newRegion = (is_unique) ?
    uniqueGroupMap_.insert( std::make_pair( regionname, csmp::Region<dim>( regionname,
                            static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
    :
    groupMap_.insert( std::make_pair( regionname, csmp::Region<dim>( regionname,
                      static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  // if region was inserted successfully	
  if ( newRegion.second )
    (*newRegion.first).second.AccumulateAll( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().RootNode( 0 ),
                                             reestablishNeighborConnectivity );

  else {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::CreateRegionFromRootNode:",
                       regionname, "region could not be formed." );
    return false;
  }

  // checking that all elements and nodes were discovered
  if ( (*newRegion.first).second.Elements() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() or
       (*newRegion.first).second.Nodes() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Nodes() )
    return false;

  return true;
}




/**
SKM trying to make sense of Andrew Bromage's undocumented code:
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::CreateRegionFromLargestComponent( const char* regionname, bool is_unique, bool reestablishNeighborConnectivity )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  // does this region already exist
  if ( ContainsRegion( regionname ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::CreateRegionFromRootNode:",
                       regionname, "region already exists; nothing was done." );
    return false;
  }

  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
    newRegion = (is_unique) ?
    uniqueGroupMap_.insert( std::make_pair( regionname, csmp::Region<dim>( regionname,
                            static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
    :
    groupMap_.insert( std::make_pair( regionname, csmp::Region<dim>( regionname,
                      static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  // if region was inserted successfully
  if ( newRegion.second ) {
    (*newRegion.first).second.FromLargestComponent( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh(), reestablishNeighborConnectivity );
  }
  else {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::CreateRegionFromRootNode:",
                       regionname, "region could not be formed." );
    return false;
  }

  // checking that all elements and nodes were discovered
  if ( (*newRegion.first).second.Elements() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() or
       (*newRegion.first).second.Nodes() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Nodes() )
    return false;

  return true;
}






/**
Removes a named Region object from the map of regions which is stored
inside of the Model object. If the region was unique and element
region IDs were set to it, these are reset to ULONG_MAX.

@param  regionName The name of the group object which shall be removed.


@section messages Messages

@note: this can trigger removal of elements

@note If the region which shall be removed does not exist, the method reports a warning.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::RemoveRegion( const char* regionName, bool delete_elements )
{
  // check whether region exists (should be a notice only, nothrow)
  if ( !ContainsRegion( regionName ) )
    throw csmp::Exception( WARNING,
                           "RegionsInterface<dim,REGION_COMPLEX>::RemoveRegion",
                           "region did not exist: ",
                           regionName );

  // finding the region in the corresponding map
  typename std::map<std::string, csmp::Region<dim> >::iterator
    iterRegion( groupMap_.find( std::string( regionName ) ) ),
    iterUniqueRegion( uniqueGroupMap_.find( std::string( regionName ) ) );

  // alerting user that other regions may be accidentally damaged by deleting non-unique regions
  if ( delete_elements ) {
    if ( iterRegion != groupMap_.end() ) {
      ErrorHandler::Instance().notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::RemoveRegion:",
                                       "Deleting the elements of a non-unique region: ", regionName );
    }
    else {
      auto& subdomain = iterUniqueRegion->second;

      REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) );
      auto& meshMgr = regionComplex->Mesh();
      auto spatialDimensions = subdomain.ElementSpatialDimensions();
      auto& elementVector = subdomain.ElementVector();

      // 1. Delete elements              
      for ( size_t i = 0U; i < elementVector.size(); i++ ) {
        // 1.1 Remove this element from its neighbour's connections
        Element<dim>* e = elementVector[i];
        auto& neighbourVector = e->NeighborElementVector();
        for ( size_t j = 0U; j < neighbourVector.size(); j++ ) {
          if ( neighbourVector[j] == NULL ) continue;
          auto& nnVector = neighbourVector[j]->NeighborElementVector();
          for ( size_t k = 0U; k < nnVector.size(); k++ ) {
            auto& nn = nnVector[k];
            if ( nn == e )
              nn = NULL;
          }
        }

        // 1.2. Remove it				
        meshMgr.Erase( e );
        elementVector.erase( elementVector.begin() + i );
        elementVector.swap( elementVector );
        i--;
      }

      // 2. Rebuild node connections if necessary			  
      if ( spatialDimensions.second == dim ) {
        // This is a region whose dimension is dim, so interior
        // nodes must be removed.
        for ( auto nit = subdomain.InteriorNodesBegin(); nit != subdomain.InteriorNodesEnd(); ++nit ) {
          meshMgr.Erase( *nit );
        }

        // Update node connections on the region's perimeter nodes that were retained.
        meshMgr.RebuildParentRelationships( subdomain.PerimeterNodesBegin(), subdomain.PerimeterNodesEnd() );
      }
      else {
        // This is a region whose dimension is less than dim (i.e. a boundary or split boundary). Just update nodes.
        meshMgr.RebuildParentRelationships( subdomain.NodesBegin(), subdomain.NodesEnd() );
      }
    }
  }


  // if the region was found in the list, it is erased
  if ( iterRegion != groupMap_.end() )
    groupMap_.erase( std::string( regionName ) );
  if ( iterUniqueRegion != uniqueGroupMap_.end() )
    uniqueGroupMap_.erase( std::string( regionName ) );

} // end RemoveRegion


/**
Removes specific elements according to the set of elmt_numbers.

@param  region_name The name of region which includes the elements which shall be removed.
@param  elmt_numbers the element numbers which shall be removed.
@return the number of the removed elements is returned.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::RemoveElements( const char* region_name, const std::set<long>& elmt_numbers )
{
  // check whether region exists (should be a notice only, nothrow)
  if ( !ContainsRegion( region_name ) )
    throw csmp::Exception( WARNING,
                           "RegionsInterface<dim,REGION_COMPLEX>::RemoveRegion",
                           "region did not exist: ",
                           region_name );

  // finding the region in the corresponding map
  typename std::map<std::string, csmp::Region<dim> >::iterator iterUniqueRegion( uniqueGroupMap_.find( std::string( region_name ) ) );

  // deleting the elements according to elmt_numbers
  auto& subdomain = iterUniqueRegion->second;
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) ) ;
  auto& meshMgr = regionComplex->Mesh();  
  auto& elementVector = subdomain.ElementVector();
  const csmp::Index key_enr = regionComplex->Database().StorageKey( "element number" );

  long removed_elements( 0 );
  std::set<Node<dim>*> candidate_nodes;  
  for ( size_t i = 0U; i < elementVector.size(); i++ ) {
    // 1.1 Remove this element from its neighbour's connections
    Element<dim>* e = elementVector[i];
    if ( e->AtBoundary() != NOT && e->AtBoundary() != IRREGULAR ) continue;
    long idx = static_cast<long>(e->Read( key_enr ));
    if ( elmt_numbers.find( idx ) == elmt_numbers.end() ) continue;    
    
    // 1.2. find candidate nodes belonging to this element
    for ( size_t j = 0U; j < e->Nodes(); j++ ) {
      Node<dim>* n = e->N( j );
      candidate_nodes.insert( n );
    }    
    meshMgr.Erase( e );
    elementVector.erase( elementVector.begin() + i );
    elementVector.swap( elementVector );
    removed_elements++;
    i--;
  }
  for ( size_t i = 0U; i < elementVector.size(); i++ ) {
    Element<dim>* e = elementVector[i];
    // remove this element from its neighbour's connections
    auto& neighbourVector = e->NeighborElementVector();
    for ( size_t j = 0U; j < neighbourVector.size(); j++ ) {
      if ( neighbourVector[j] == NULL ) continue;
      auto& nnVector = neighbourVector[j]->NeighborElementVector();
      for ( size_t k = 0U; k < nnVector.size(); k++ ) {
        if ( nnVector[k] == NULL ) continue;
        auto& nn = nnVector[k];        
        if ( nn->Nodes() == 0 ) {
          nn = NULL;
        }
      }
    }
  }
  for(Node<dim>* n: candidate_nodes )
  {
    if ( n->Parents() == 0 )
      meshMgr.Erase( n );
  }
  subdomain.CreateNodePointerVector();
  subdomain.EstablishNeighborConnectivity();
  subdomain.IdentifyPerimeter();

  // remove those elements from the default region 'Model'
  const char* default_model = "Model";
  if ( ContainsRegion( default_model ) ) {
    // finding the region in the corresponding map
    typename std::map<std::string, csmp::Region<dim> >::iterator iterRegion( groupMap_.find( std::string( default_model ) ) );

    auto& defulat_domain = iterRegion->second;
    auto& elementVector = defulat_domain.ElementVector();

    for ( size_t i = 0U; i < elementVector.size(); i++ ) {
      // remove this element 
      Element<dim>* e = elementVector[i];
      if ( e->Nodes() == 0 ) {
        elementVector.erase( elementVector.begin() + i );
        elementVector.swap( elementVector );
        i--;
      }
    }
    for ( size_t i = 0U; i < elementVector.size(); i++ ) {
      Element<dim>* e = elementVector[i];
      // remove this element from its neighbour's connections
      auto& neighbourVector = e->NeighborElementVector();
      for ( size_t j = 0U; j < neighbourVector.size(); j++ ) {
        auto& ne = neighbourVector[j];
        if ( ne == NULL ) continue;
        if ( ne->Nodes() == 0 ) {
          ne = NULL;
        }
      }
    }
    defulat_domain.CreateNodePointerVector();
    defulat_domain.EstablishNeighborConnectivity();
    defulat_domain.IdentifyPerimeter();
  }

  meshMgr.RebuildParentRelationships( subdomain.NodesBegin(), subdomain.NodesEnd() );
  
  return removed_elements;
} // end RemoveElements


// -----------------------------------------------
// Binary input/output
// -----------------------------------------------

/**
Writes all the unique and non-unique regions to binary file.
The unique regions are written first.

@attention ONLY METHOD that provides complete information about each region:
interior vs. perimeter nodes and element; boudary faces etc.

@note the element and node Idx indices must have a global unique numbering for this method to work.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::OutputAllRegionsToBinary( const char* file_name ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( regionComplex.Database() );

  std::string bin_file( file_name );
  std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::OutputAllRegionsToBinary:",
                       bin_file, "file could not be opened; nothing was done." );
    return;
  }

  // -----------------------------
  // 1. File header
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "REGFHEDR" );

    std::string heading( "RegionInterface::OutputAllRegionsToBinary: " );
    heading += "region information for Model '";
    heading += regionComplex.Name();
    heading += "' to file: ";
    heading += bin_file;
    heading += "'.";

    skm_C_fwrite( fp, heading.c_str() );
  }

  std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::OutputAllRegionsToBinary: regions written to binary file: ";
  std::cout << "Unique regions: ";

  // -----------------------------
  // 2. writing the unique regions
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "UNIQREGN" );

    size_t records = this->UniqueRegions();

    // writing number of unique regions
    fp.write( (char*)&records, sizeof( size_t ) );

    for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
          git = UniqueRegionsBegin(); git != UniqueRegionsEnd(); ++git )
    {
      BinaryFileSectionWrite hdr( fp, "ONE_REGN" );
      (*git).second.WriteDomainIndexesToBinaryFile( fp );
      domainVariablesOut( fp, (*git).second, database );
      std::cout << (*git).first << " ";
    }
  }
  std::cout << ", non-unique (potentially overlapping) regions: ";

  // -----------------------------
  // 3. writing non-unique regions
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "NONUREGN" );

    size_t records = this->Regions() - this->UniqueRegions();

    fp.write( (char*)&records, sizeof( size_t ) );

    for ( auto git = RegionsBegin(); git != RegionsEnd(); git++ )
    {
      BinaryFileSectionWrite hdr( fp, "ONE_REGN" );
      (*git).second.WriteDomainIndexesToBinaryFile( fp );            
      if ( (*git).second.InteriorElementsBegin() != (*git).second.InteriorElementsEnd() ) {
        domainVariablesOut( fp, (*git).second, database );
      }
      std::cout << (*git).first << " ";
    }
    std::cout << std::endl;
  }

  // -----------------------------
  // 3. writing model variables
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "MODLVARS" );

    // writing the Model variables here TODO: check for correctness
    domainVariablesOut( fp, regionComplex, database );
  }

  // -----------------------------
  // 4. writing file footer
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "REGFFOTR" );
  }

  // 5. clean up
  fp.close();
  std::cout << "\nRegionInterface<" << dim << ",REGION_COMPLEX>::OutputAllRegionsToBinary: file '";
  std::cout << bin_file << "' has been successfully written.\n";

} // end OutputAllRegionsToBinary






/**
Reads all regions stored by OutputAllRegionsToBinary() into Interface,
constructing them using the SubDomainInfo data, thereby avoiding the costly
re-initialisation.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::InputAllRegionsFromBinary( const char* file_name, const std::set<std::string>* subset_variables )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );

  std::string bin_file( file_name );
  std::fstream fp( bin_file.c_str(), std::ios::in | std::ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputAllRegionsFromBinary:",
                       bin_file, "file could not be opened; nothing was done." );
    return;
  }

  // 1. File header
  {
    BinaryFileSectionRead hdr( fp, "REGFHEDR" );
    char  text[500U];

    skm_C_fread( fp, text );
    std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::InputAllRegionsFromBinary: Reading file header:\n\t" << text << std::endl;
  }
  std::cout << "\n\timporting the regions: ";

  const PropertyDatabase<dim>& database( static_cast<const REGION_COMPLEX<dim>& >(*this).Database() );

  // getting a reference to this newly created master region
  std::cout << "\n\tunique regions: ";

  // -----------------------------
  // 2. unique regions
  // -----------------------------

  std::deque<csmp::Node<dim>*>	 nodes;
  std::deque<csmp::Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( &regionComplex.Mesh(), nodes, elmts );
  std::sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  std::sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  SubDomainInfo  info;

  {
    BinaryFileSectionRead hdr( fp, "UNIQREGN" );
    size_t  records( 0 );  // region records

                           // number of regions
    fp.read( (char*)&records, sizeof( size_t ) );
    if ( records > 0 )
      // reading the regions sequentially
      for ( size_t i = 0U; i<records; i++ )
      {
        BinaryFileSectionRead hdr( fp, "ONE_REGN" );
        // reading name and element indices for each unique region
        readDomainIndexesFromBinaryFile( dim, fp, info );
        // reconstruct the region
        std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
          it = uniqueGroupMap_.insert( std::make_pair( info.name, csmp::Region<dim>( database, nodes, elmts, info ) ) );

        if ( !it.second )
          throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputAllRegionsFromBinary:",
                                 info.name, "Region could not be formed; issue with binary file." );

        std::cout << "\n\t\t'" << (*it.first).first << "'(" << (*it.first).second.Elements() << " elements).";
        domainVariablesIn( fp, (*it.first).second, database );
      }
    else csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputAllRegionsFromBinary:",
                            bin_file, "does not contain any unique region descriptions; no regions were initialised." );
  }

  // -----------------------------
  // 3. reading non-unique regions
  // -----------------------------
  std::cout << "\n\tnon-unique regions: ";
  {
    BinaryFileSectionRead hdr( fp, "NONUREGN" );

    // getting number of non-unique region records from file
    size_t records( 0 );
    fp.read( (char*)&records, sizeof( size_t ) );
    if ( records > 0 )
      // reading the regions sequentially
      for ( size_t i = 0U; i<records; i++ )
      {
        BinaryFileSectionRead hdr( fp, "ONE_REGN" );

        // reading name and element indices for each unique region
        readDomainIndexesFromBinaryFile( dim, fp, info );
        // if the region info record is not empty the region is reconstructed
        std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
          it = groupMap_.insert( std::make_pair( info.name, csmp::Region<dim>( database, nodes, elmts, info ) ) );
        if ( !info.interior_elmts.empty() ) {
          if ( !it.second )
            throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputAllRegionsFromBinary:",
                                   info.name, "Region could not be formed; issue with binary file." );

          std::cout << "\n\t\t'" << (*it.first).first << "'(" << (*it.first).second.Elements() << " elements).";
          domainVariablesIn( fp, (*it.first).second, database );
        }
      }
    else csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputAllRegionsFromBinary:",
                            bin_file, "does not contain any non-unique region descriptions; no regions were initialised." );
  }

  // -----------------------------
  // 4. reading model vars
  // -----------------------------

  {
    BinaryFileSectionRead hdr( fp, "MODLVARS" );

    // reading the Model variables here TODO: check for correctness
    domainVariablesIn( fp, regionComplex, database );
  }


  // -----------------------------
  // 5. file footer
  // -----------------------------
  {
    BinaryFileSectionRead hdr( fp, "REGFFOTR" );
  }

  // -----------------------------
  // 6. cleanup
  // -----------------------------

  fp.close();
  std::cout << "\n\nRegionInterface<" << dim << ",REGION_COMPLEX>::InputAllRegionsFromBinary: file '";
  std::cout << bin_file << "' has been read successfully.\n";

} // end InputAllRegionsFromBinary






/**
Creates a new binary file into which the region is stored.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::OutputRegionToBinary( const char* region_name, const char* file_name ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  // does this region already exist
  if ( !ContainsRegion( region_name ) ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::OutputRegionToBinary:",
                       region_name, "region does not exist; nothing was done." );
    return;
  }

  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( regionComplex.Database() );

  std::string bin_file( file_name );
  std::string heading( "RegionInterface<dim,REGION_COMPLEX>::OutputRegionToBinary: output region '" );
  heading += region_name;
  heading += "' of Model '";
  heading += regionComplex.Name();
  heading += "' to file: ";
  heading += bin_file;
  heading += "'.";

  std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::OutputRegionToBinary:",
                       bin_file, "file could not be opened; nothing was done." );
    return;
  }

  // 1. writing the file header
  skm_C_fwrite( fp, heading.c_str() );

  // 2. writing number of regions=1
  const size_t records( 1U );
  fp.write( (char*)&records, sizeof( size_t ) );

  // 3. writing the name of region
  skm_C_fwrite( fp, region_name );

  // 4. writing the element records of the region
  const csmp::Region<dim>& subdomain( Region( region_name ) );
  std::vector<size_t>  elmtIDs;
  subdomain.MemberElementIndexes( elmtIDs );
  skm_C_fwrite( fp, elmtIDs );

  // 5. writing the variable values associated with the region to file
  domainVariablesOut( fp, subdomain, database );

  // 6. cleaning up
  fp.close();
  std::cout << "\nRegionInterface<" << dim << ",REGION_COMPLEX>::OutputRegionToBinary:: region '";
  std::cout << region_name << "' has been successfully written to: '";
  std::cout << bin_file << "'" << std::endl;

} // end OutputRegionToBinary




/**
For all regions in the current model, this method writes name and
contained elements into a binary file with the name 'file_name'. The
extension '.dat' is appended.

The data in the file are organised as follows:

1. header line
2. number of region
3. for each region, name followed by array of the indexes of the elements
stored in this region (correct=unique numbering is expected).
end of file

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::AppendRegionsToBinary( const char* file_name ) const
{
  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( regionComplex.Database() );

  std::string bin_file( file_name );
  std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::app | std::ios::binary );
  if ( !fp.is_open() ) {
    std::cerr << "\nRegionInterface<dim,REGION_COMPLEX>::AppendRegionsToBinary: file: '" << bin_file;
    std::cerr << "' could not be opened." << std::endl;
    return;
  }

  std::vector<size_t>  elmtIDs;

  // --------------------------
  // writing the unique regions
  // --------------------------
  // writing number of unique regions
  size_t records( this->UniqueRegions() );
  fp.write( (char*)&records, sizeof( size_t ) );

  for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
        git = UniqueRegionsBegin(); git != UniqueRegionsEnd(); ++git )
  {
    // writing name of the region
    skm_C_fwrite( fp, (*git).first.c_str() );
    // writing the element records of the region
    (*git).second.MemberElementIndexes( elmtIDs );
    skm_C_fwrite( fp, elmtIDs );
    // writing the values ofthe variables associated with the region
    domainVariablesOut( fp, (*git).second, database );
  }

  // --------------------------
  // writing non-unique regions
  // --------------------------
  records = this->Regions() - this->UniqueRegions();
  fp.write( (char*)&records, sizeof( size_t ) );

  for ( auto git = RegionsBegin(); git != RegionsEnd(); git++ ) {
    // avoiding the region which is the master region since it was already written before
    skm_C_fwrite( fp, (*git).first.c_str() );
    (*git).second.MemberElementIndexes( elmtIDs );
    skm_C_fwrite( fp, elmtIDs );
    domainVariablesOut( fp, (*git).second, database );
  }

  // writing the Model variables here TODO: check for correctness
  domainVariablesOut( fp, regionComplex, database );

  // cleaning up
  fp.close();
  std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::AppendRegionsToBinary: regions have been successfully written to: '";
  std::cout << bin_file << "'" << std::endl;

} // end AppendRegionsToBinary




/**
Reads a single non-unique region from the supplied binary file.

@author SKM
@date 5/4/2016
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::InputRegionFromBinary( const char* region_name, bool is_unique, const char* file_name )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  std::string bin_file( file_name );

  std::fstream fp( bin_file.c_str(), std::ios::in | std::ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionFromBinary:",
                       bin_file, "file could not be opened; nothing was done." );
    return;
  }

  // reading the file header and printing it to screen
  char  text[500U];
  skm_C_fread( fp, text );
  std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::InputRegionFromBinary: Reading file header:\n\t" << text << std::endl;

  // reading how many records will follow (should be 1 for this method)
  size_t  records( 0 );  // region records
  fp.read( (char*)&records, sizeof( size_t ) );

  // there should just be the region of interest, else something is wrong
  if ( records == 1 ) {
    // region name
    skm_C_fread( fp, text );
    assert( std::strcmp( text, region_name ) == 0 );
    // element ids
    std::vector<size_t>  elmtIDs; // unsigned integer element identifiers to be read
    skm_C_fread( fp, elmtIDs );
    // creating the region
    const bool reestablishNeighborConnectivity( true );
    CreateRegionFromRootNode( region_name, is_unique, reestablishNeighborConnectivity );
    // reading the associated variable values
    const PropertyDatabase<dim>& database( static_cast<const REGION_COMPLEX<dim>& >(*this).Database() );
    csmp::Region<dim> binRegion( text, database );
    // elements were already accumulated so only the region variables need to be read
    domainVariablesIn( fp, this->Region( region_name ), database );
  }
  else csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionFromBinary:",
                          region_name, "inconsistent binary record for region; no data could be read." );

} // end InputRegionFromBinary





/**
Adds regions to an existing Model using the information from file.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::InputRegionsFromBinary( const char* file_name )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( static_cast<const REGION_COMPLEX<dim>& >(*this).Database() );

  std::string bin_file( file_name );
  size_t               records( 0 );
  std::vector<size_t>  elmtIDs;

  std::fstream fp( bin_file.c_str(), std::ios::in | std::ios::binary );
  if ( !fp.is_open() ) {
    std::cerr << "\nRegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary: file: '" << bin_file;
    std::cerr << "' could not be opened." << std::endl;
    return;
  }

  // skipping header and the first "master region"
  char  text[256U];
  skm_C_fread( fp, text );
  fp.read( (char*)&records, sizeof( size_t ) ); // record
  assert( records == 1 );
  skm_C_fread( fp, text );    // name of region
  skm_C_fread( fp, elmtIDs ); // element indices
                              // TODO: deal with this redundant step although it does not affect many variable values

                              // XXX AJB FIXME
  domainVariablesIn( fp, this->Region( "Model" ), database );

  auto& mesh = regionComplex.Mesh();

  // ------------------
  // unique regions
  // ------------------
  fp.read( (char*)&records, sizeof( size_t ) );
  if ( records > 0 )
    // reading the regions sequentially
    for ( size_t i = 0U; i<records; i++ )
    {
      // reading name and element indices for each unique region
      skm_C_fread( fp, text );
      skm_C_fread( fp, elmtIDs );
      csmp::Region<dim> binRegion( text, database );
      domainVariablesIn( fp, binRegion, database );
      // if the region is not empty it is assembled
      if ( !elmtIDs.empty() ) {
        // adding the regions to Model object
        std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
          it = uniqueGroupMap_.insert( std::make_pair( text, binRegion ) );
        //   ^^^^^^^^^^^^^^^
        if ( it.second )
          (*it.first).second.AccumulateByNumber( mesh, elmtIDs );

        if ( (*it.first).second.Elements() == 0U ) {
          uniqueGroupMap_.erase( it.first );
          throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                                 text, "Region could not be formed; dataset contained no elements." );
        }
      }
    }
  else csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                          bin_file, "does not contain any unique region descriptions; none were initialised." );

  // ------------------
  // non-unique regions
  // ------------------
  fp.read( (char*)&records, sizeof( size_t ) );
  if ( records > 0U )
    for ( size_t i = 0U; i<records; i++ )
    {
      skm_C_fread( fp, text );
      skm_C_fread( fp, elmtIDs );
      csmp::Region<dim> binRegion( text, database );
      // the master region should not have been stored to disk via the corresponding append to binary function
      domainVariablesIn( fp, binRegion, database );

      if ( !elmtIDs.empty() ) {
        std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
          it = groupMap_.insert( std::make_pair( text, binRegion ) );
        //   ^^^^^^^^^
        if ( it.second )
          (*it.first).second.AccumulateByNumber( mesh, elmtIDs );

        if ( (*it.first).second.Elements() == 0U ) {
          this->groupMap_.erase( it.first );
          throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                                 text, "Region could not be formed; dataset contained no elements." );
        }
      }
    }
  else csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                          bin_file, "does not contain any non-unique region descriptions; none were initialised." );

  /// @todo (3-D) We misuse the regions file here to store Model variables
  domainVariablesIn( fp, regionComplex, database );

  fp.close();

  if ( Regions() == records ) {
    std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary: " << records;
    if ( UniqueRegions() == 1U ) std::cout << " unique region read successfully from '";
    else std::cout << " unique regions read successfully from '";
    std::cout << bin_file << "'" << std::endl;
    std::cout.flush();
  }

} // end inputRegionsFromBinary







/**
RegionsFromPropertyValues() defines one group for each value of the target
property.  The user will be prompted to assign a name to each group as these
are being created. The names of successfully created groups are returned
into the list argument.

@param prop the property on the basis of whose variations
the group regions will be defined.

The names of the newly created groups are returned into an STL set
which uses the less<> functional to order the names alphabetically. If
the set is not empty, it will be erased before the group names are
stored within it.

@section implementation Implementation

Thus far, the method only handles scalar variables. Since there is always
just one property value per element, the groups will be unique.

@section application Application

RegionsFromPropertyValues() is useful if the desired subregions of a model
coincide with changes of a specific property. Clearly, one does not want
to apply this method on a continuously changing property, since it might
produce as many Regions as there are finite-elements in the mesh.

@section messages Messages

The method is interactive and will prompt the user for the names of the
groups which are created in the course of its execution.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues( const char* prop,
                                                                          std::set<std::string>& group_names )
{
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  if ( !group_names.empty() )
    group_names.erase( group_names.begin(), group_names.end() );

  csmp::Index  prop_key = regionComplex->Database().StorageKey( prop );

  if ( prop_key.type != SCALAR )
    throw csmp::Exception( INFO, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromPropertyValues",
                           "method can only be applied to SCALAR variables",
                           "no assignments were made" );

  auto& mesh = regionComplex->Mesh();
  std::map<double64, std::string>  groups;
  ScalarVariable                  sc;

  // 1. Making a map with one entry for each property value
  // ------------------------------------------------------

  // traversal of the existing mesh nodes to find all its elements	
  std::deque<csmp::Element<dim>*>	elmts;
  std::deque<csmp::Node<dim>*>	nodes;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  switch ( prop_key.place )
  {
    case NODE:
      for ( auto nit : nodes ) {
        nit->Read( prop_key, sc );
        groups[sc()] = "undefined";
      }
      break;
    case ELEMENT_INTEGRATION_POINT:
      for ( auto eit : elmts ) {
        for ( size_t i = 0U; i < eit->IntegrationPoints(); i++ )
        {
          eit->Read( i, prop_key, sc );
          groups[sc()] = "undefined";
        }
      }
      break;
    case ELEMENT:
      for ( auto eit : elmts ) {
        eit->Read( prop_key, sc );
        groups[sc()] = "undefined";
      }
      break;
    default:
      throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromPropertyValues",
                             "Property placement not recognized; presumable placement REGION is not allowed in this context. " );
  }

  // 2. Prompting user for group names and making groups
  // ---------------------------------------------------
  std::cout << "\nModel<" << dim << ">::FormRegionsFromPropertyValues: Generating region names";
  std::cout << " corresponding to unique values of the property: '" << prop << "' ";
  std::cout << "[" << regionComplex->Database().Unit( prop ) << "]." << std::endl;
  std::string  gname( "region_" );
  gname += prop;
  size_t  group_idx( 0U );
  char    num[30U];

  for ( typename std::map<double64, std::string>::iterator
        it = groups.begin(); it != groups.end(); it++ )
  {
    std::cout << "\nProperty value: " << (*it).first;
    sprintf( num, "%lu", group_idx++ );
    (*it).second = std::string( gname + num );

    if ( (*it).second != "undefined" )
    {
      assert( !ContainsRegion( (*it).second.c_str() ) );
      FormRegionFrom( (*it).second.c_str(), prop,
                      (*it).first - std::numeric_limits<double64>::epsilon(),
                      (*it).first + std::numeric_limits<double64>::epsilon(), true );

      group_names.insert( (*it).second );
    }
    else
      std::cout << "\nModel<" << dim << ">::FormRegionsFromPropertyValues: check predefined region name." << std::endl;
  }

} // end RegionsFromPropertyValues



template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionName, const csmp::Region<dim>& region, bool unique/*=false */ )
{
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool> it;
  if ( unique )
    it = this->uniqueGroupMap_.insert( std::make_pair( regionName, region ) );
  else
    it = this->groupMap_.insert( std::make_pair( regionName, region ) );
  return it.second;
}



/**
FormAndAddRegion() used combined constraints suppplied in the form of
a constraints object to determine which elements shall be used to
form a (unique/non-unique) group with the target name.

@param groupname the name of the group that shall be formed and a reference to
the initialized PropertyConstraints object which must contain the
ranges of the variables that shall be used to discriminate elements
that shall be accumulated into the group.

@section application Application

To form groups from dynamic criteria like the combined pressure and
temperature ranges that get computed via Algorithms or other criteria
that are not known at the onset of a computation.

@section messages Messages

The method reports an INFO if no values fall into the target ranges and
a FATAL_ERROR if the group cannot be created because the name is already
in use.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* groupname,
                                                           PropertyConstraints& constraints,
                                                           bool unique_group )
{
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  std::string output_region( groupname );
  if ( ContainsRegion( groupname ) ) {
    std::cout << "\nModel<" << dim << ">::FormRegionFrom: WARNING: region '" << groupname;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  // the 'bool' member of pair indicates whether insertion into map worked or not
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>  it;
  if ( unique_group )
    it = uniqueGroupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
  else
    it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( it.second )
  {
    constraints.InitializePropertyIndices( regionComplex->Database() );
    auto& mesh = regionComplex->Mesh();
    (*it.first).second.AccumulateWithinRange( mesh, constraints );

    // removing the group if it contains no elements
    if ( (*it.first).second.Elements() == 0U ) {
      groupMap_.erase( it.first );
      csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                         "Region could not be formed", output_region.c_str() );
      return false;
    }
  }
  else {
    csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                       "Region could not be formed", output_region.c_str() );
    return false;
  }

  return true;

} // end FormRegionFrom




/**
Forms a group of finite elements who host property values which lie within
a user-specified range. The new group is added to the Model group map.

The selection process implies that this type of region can include several
regions inside your problem domain. For instance, all shale horizons as
identified by a clay content between 50 to 100% may constitute a new group.

When a new group is formed a group internal flag will be assigned to
each node, constraint point and element.
There are two group-internal object flags, PLAIN and BOUNDARY. When a
new group is formed, the Region method IdentifyBoundaryAs() assigns the
group-internal object flags, depending on whether the nodes, constraint
points or elements in the group lie at the group boundary or inside of
the group. Elements are assigned a boundary flag if at least one of
their faces coincides with the group boundary.

@section arguments Input Arguments

The newly formed group of finite elements will have a name specified
by the first method argument. The selection criterion is that values of
the physical variable identified by the second argument, are within the
open interval given by [min,max].

@section implementation Implementation

Once all member elements of the new group have been identified, the group
analyzes which elements, nodes, and constraint points lie at its boundary.
Subsequently the group is added to a map of groups that is contained in
the Model object.

Inn this and other group-forming methods Regions are first added to the
group map before they are filled, avoiding a costly copy construction
of non-empty groups.

@section application Application

To form groups on the basis of characteristic material properties.

@section messages Messages

A warning is issued if a group with the same name already exists or if no
elements with the desired properties were found.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* groupname,
                                                           const char* prop,
                                                           double64 min, double64 max,
                                                           bool unique_group )
{
  std::string output_region( groupname );
  if ( ContainsRegion( groupname ) ) {
    std::cout << "\nModel<" << dim << ">::FormRegionFrom: WARNING: region '" << groupname;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // the 'bool' member of pair indicates whether insertion into map worked or not
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>  it;
  if ( unique_group )
    it = uniqueGroupMap_.insert( std::make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
  else
    it = groupMap_.insert( std::make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  if ( it.second )
  {
    REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
    auto& mesh = regionComplex->Mesh();
    (*it.first).second.AccumulateWithinRange( mesh, prop, min, max );

    // removing the group if it contains no elements
    if ( (*it.first).second.Elements() == 0U ) {
      groupMap_.erase( it.first );
      csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                         "Region could not be formed", output_region.c_str() );
      return false;
    }
  }
  else {
    csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                       "Region could not be formed. ", output_region.c_str() );
    return false;
  }

  return true;

} // end FormRegionFrom



/// Forms group from the supplied vector of global Element IDs (0..n-1)
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* groupname, std::vector<size_t>& element_ids, bool unique_region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( element_ids.empty() )
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom (element numbers)",
                           "Region could not be formed because input element number vector is empty", groupname );

  std::string output_region( groupname );
  if ( ContainsRegion( groupname ) ) {
    std::cout << "\nModel<" << dim << ">::FormRegionFrom (element numbers): WARNING: region '" << groupname;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
    it = (unique_region == true) ?
    uniqueGroupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) ) :
    groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  if ( it.second ) {
    REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
    auto& mesh = regionComplex->Mesh();
    (*it.first).second.AccumulateByNumber( mesh, element_ids );

    // removing the group if it contains no elements
    if ( (*it.first).second.Elements() == 0U ) {
      if ( unique_region ) uniqueGroupMap_.erase( it.first );
      else groupMap_.erase( it.first );
      csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom (element numbers)",
                         "Region could not be formed", output_region.c_str() );
      return false;
    }
  }
  else {
    csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom (element numbers)",
                       "Region could not be formed. Does this region already exist?", output_region.c_str() );
    return false;
  }

  return true;

} // FormRegionFrom


/**
Forms a Region from a set of region names by using the method MergeRegions.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionname, const std::set<std::string>& region_names )
{
  std::string output_region( regionname );
  if ( ContainsRegion( regionname ) ) {
    std::cout << "\nModel<" << dim << ">::FormRegionFrom (element numbers): WARNING: region '" << regionname;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  this->MergeRegions( region_names, output_region.c_str() );
}



/**
@author P Lang

This forms a region from elements eligible as reported from elementComp.

ElementComp is a model of binary predicate, i.e.
@code
template<size_t dim>
struct ElementsLessX
{
bool operator () ( Element<dim> const* ePtr ) const
{
return ...is ePtr eligible?;
}
};
@endcode
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
template<template<size_t> class ElementComp>
bool RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* newRegionName, ElementComp<dim> const& elementComp, const char* hostRegion )
{
  static_cast<REGION_COMPLEX<dim>*>(this)->UpdateIndices();
  csmp::Region<dim> const& rref( this->Region( hostRegion ) );
  std::vector<size_t> elementIds;
  elementIds.reserve( rref.Elements() );
  for ( typename csmp::Region<dim>::SimplexContainer::const_iterator it( rref.ElementsBegin() ); it != rref.ElementsEnd(); ++it )
    if ( elementComp( (*it) ) )
      elementIds.push_back( (*it)->Idx() );
  std::vector<size_t>( elementIds ).swap( elementIds );
  if ( elementIds.empty() )
    return false;
  this->FormRegionFrom( newRegionName, elementIds );
  return true;
}



/**
Forms regions using the element ID containers stored in the model topology
object. The regions are numbered in their alphabetical order and these numbers are assigned to the material ID of the element class.

@note the material IDs may later be overwritten by rocktypes .

@section arguments Input Arguments

The method takes the information on the basis of which the regions
are formed from a model topology object via const reference.

@section implementation Implementation

The assumption is made that a master Region already exists.
The method uses the region interface 'AccumulateElements()'.

@section application Application

When a model is build and external topological information is available,
such as subvolume names etc., this method provides a convenient and
rapid way to form computational domains in a Model on the basis
of the supplied information.

@section messages Messages

The method reports which regions are being formed.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::FormRegionsFrom( const ModelTopology& topo )
{
  // 1. getting the names of the regions
  std::list<std::string> regions;
  topo.Out( regions );

  // 2. assigning the regions to groups in the Model
  std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::FormRegionsFrom: Forming the regions: ";
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  auto& mesh = regionComplex->Mesh();

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  int32 new_regions( 0U );
  for ( typename std::list<std::string>::const_iterator lit = regions.begin(); lit != regions.end(); lit++ )
  {
    std::string group_name( *lit );
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = uniqueGroupMap_.insert( make_pair( group_name, csmp::Region<dim>( group_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    // if the region was successfully inserted
    if ( it.second )
      {
        // making a list of the element numbers
        std::vector<size_t>  element_ids;
        element_ids.reserve( topo.ElementsOfRegion( (*lit).c_str() ) );
        copy( topo.ElementsOfRegionBegin( (*lit).c_str() ),
              topo.ElementsOfRegionEnd( (*lit).c_str() ),
              back_inserter( element_ids ) );

        // assigning the element IDs to the group & cleaning up
        (*it.first).second.AccumulateByNumber( mesh, element_ids );
        element_ids.erase( element_ids.begin(), element_ids.end() );

        // removing the group if it contains no elements
        if ( (*it.first).second.Elements() == 0U ) {
            uniqueGroupMap_.erase( it.first );
            csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFrom",
                               "Region could not be formed", (*lit).c_str() );
          }
        else {
            // assigning unique material IDs to the element members of the region
            for ( auto eit=(*it.first).second.ElementsBegin(); eit!=(*it.first).second.ElementsEnd(); ++eit )
              (*eit)->Material_ID( new_regions );
            // reporting the name of the newly generated region
            std::cout << group_name << " ";
            new_regions++;
          }
      }
    else
      throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFrom",
                             "Region could not be formed. Does this region already exist?", (*lit).c_str() );
  }
  std::cout << std::endl;

  if ( new_regions == topo.ModelRegions() ) return true;
  return false;

} // end FormRegionsFrom


/**
Tests whether the elements of the region are connected to each-other.

@attention This test cannot be performed if the region has elements of different spatial
dimensions since these are not interconnected. Therefore, this method returns false if
the region consists of elements from different spatial dimensions.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::IsContiguous( const std::string& region_name ) const
{
  const csmp::Region<dim>& mref = Region( region_name );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  std::pair<int32, int32>  dimensionality = mref.ElementSpatialDimensions();
  if ( dimensionality.first > 1U ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::IsContiguous:",
                       "method can determine contiguity only for regions which consist only of same spatial dimension elements; returned false." );
    return false;
  }

  std::set<Element<dim>*>  elmts( mref.ElementsBegin(), mref.ElementsEnd() ),
    contiguous_elmts;

  floodFill( (*elmts.begin()), contiguous_elmts );
  if ( contiguous_elmts.size() != mref.Elements() ) return false;

  return true;

} // end IsContiguous





/**
Breaks non-contiguous Regions into contiguous subregions that carry the name
of the region but have a number as suffix to their name to denote the
partition.

@return The method returns the number of subgregions that were created.

@attention The master region that was successfully partitioned is removed.

@param group The name of the region that may be non-contiguous.
If so, new sbregions will be created to the name of which integers
will be appended that correspond to the number of subdomains
that are created in this process.

@return The method returns the number of contiguous subdomains which it
created.

@section implementation Implementation

The method uses the union-find algorithm, using neighbours to determine
components.

@section application Application

To automatically partition groups that consist of a multitude of
non-contiguous model subdomains so that the latter can be addressed
individually in computations.

@section messages Messages

The method will report if the group is already contiguous in which
case no changes are made.

@attention this method cannot be applied to the region model or the master region
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions( const char* group )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  const std::string  group_name( group );

  if ( group_name == "Model" ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions:",
                       "this operation is not allowed for region 'Model' or the master region." );
    return 0U;
  }

  // renumbering elements and nodes of model
  csmp::Region<dim>&  gref( Region( group ) );
  const bool          unique_group( IsUnique( group ) );

  // Perform union-find
  UnionFind<Element<dim>*> unionFind;
  auto eend = gref.ElementsEnd();
  for ( auto eit = gref.ElementsBegin(); eit != eend; eit++ ) {
    auto e = *eit;
    const size_t  neighbors( e->Neighbors() );
    for ( size_t i = 0U; i<neighbors; i++ )
      if ( e->Neighbor( i ) != NULL )
        unionFind.SameComponent( e, e->Neighbor( i ) );
  }

  std::deque<std::pair<size_t, Element<dim>*>> components;
  unionFind.Components( components );
  if ( components.size() <= 1 ) {
    std::cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions: ";
    std::cout << "region '" << group << "' is already contiguous, nothing was done." << std::endl;
    return 0U;
  }

  // Sort all elements by component
  std::vector<std::pair<Element<dim>*, Element<dim>*>> componentMemberships;
  componentMemberships.reserve( gref.Elements() );
  for ( auto eit = gref.ElementsBegin(); eit != eend; eit++ ) {
    componentMemberships.emplace_back( unionFind.resolve( *eit ), *eit );
  }
  std::sort( componentMemberships.begin(), componentMemberships.end() );

  // Turn components into regions
  std::cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions: ";
  std::cout << "region '" << group << "' is divided into the subregion(s):\n";

  size_t subgroupNum( 0 );
  auto cmCBegin = componentMemberships.begin();
  auto cmEnd = componentMemberships.end();
  while ( cmCBegin != cmEnd ) {
    auto cmCEnd = cmCBegin;
    size_t subgroupSize = 0;
    while ( cmCEnd != cmEnd && cmCBegin->first == cmCEnd->first ) {
      ++cmCEnd;
      ++subgroupSize;
    }

    ++subgroupNum;
    std::string  subgroup_name;
    char         num[128];
    sprintf( num, "%lu", subgroupNum );
    subgroup_name = group_name + num;
    std::cout << "\t\t\t'" << subgroup_name << "'";
    std::cout << " (" << subgroupSize << " elmts)" << std::endl;

    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = (unique_group) ? uniqueGroupMap_.insert( make_pair( subgroup_name, csmp::Region<dim>( subgroup_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
      : groupMap_.insert( make_pair( subgroup_name, csmp::Region<dim>( subgroup_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( !it.second )
      throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions",
                             subgroup_name.c_str(), "region could not be formed (name is probably not unique)" );
    else {
      // accumulating the subregion
      std::vector<Element<dim>*> subgroup;
      subgroup.reserve( subgroupSize );
      for ( auto it = cmCBegin; it != cmCEnd; ++it ) {
        subgroup.push_back( it->second );
      }
      (*it.first).second.Accumulate( subgroup.begin(), subgroup.end() );

      // copy all values of region properties from parent to child region
      (*it.first).second.LVS( gref.LVS() );
    }
    cmCBegin = cmCEnd;
  }

  // if the region has been partitioned succesfully and its name is not model, it will be removed
  if ( IsUnique( group ) )
    RemoveRegion( group, false );

  return subgroupNum;

} // end partitionRegionIntoContiguousSubRegions







/**
As previous method but using element Idx information to find the elements.

@attention The master region that was successfully partitioned is removed.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegionsByIdx( const char* group )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  if ( std::string( "Model" ) == group ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegionsByIdx:",
                       "this operation is not allowed for region 'Model'." );
    return 0U;
  }
  // the region gref remains intact
  csmp::Region<dim>&  gref( Region( group ) );
  const bool          unique_group( IsUnique( group ) );
  // renumbering elements and nodes of model
  gref.UpdateMemberIndexes();

  // getting a set of the element numbers of the target group
  std::set<size_t>  elements;
  for ( typename std::vector<csmp::Element<dim>*>::const_iterator
        eit = gref.ElementsBegin(); eit != gref.ElementsEnd(); eit++ )
    elements.insert( (*eit)->Idx() );

  // detecting via a flood-fill whether the group can be partitioned, else nothing is done
  std::set<size_t>  elements_contiguous_subset;
  floodFillViaIndexes( gref, (*elements.begin()), elements_contiguous_subset );
  // if the first flood-fill reached all elements of the group it is contiguous
  if ( elements.size() == elements_contiguous_subset.size() ) {
    std::cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions: ";
    std::cout << "region '" << group << "' is already contiguous, nothing was done." << std::endl;
    return 0U;
  }

  // else partitions can be created 
  std::string  group_name( group );
  std::string  subgroup_name;
  char         num[128];
  size_t       n_subgroups( 1 );

  // creating new contiguous group from the element subset
  while ( !elements.empty() )
  {
    // creating name of contiguous subgroup
    sprintf( num, "%lu", n_subgroups );
    subgroup_name = group_name + num;
    if ( n_subgroups == 1U ) {
      std::cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions: ";
      std::cout << " forming new subregion(s):\n";
    }
    std::cout << "\t\t\t'" << subgroup_name << "'";
    std::cout << " (" << elements_contiguous_subset.size() << " elmts)" << std::endl;

    // creating either a unique or non-unique group depending on uniqueness of original region
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = (unique_group) ? uniqueGroupMap_.insert( make_pair( subgroup_name, csmp::Region<dim>( subgroup_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
      : groupMap_.insert( make_pair( subgroup_name, csmp::Region<dim>( subgroup_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( !it.second )
      throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions",
                             subgroup_name.c_str(), "region could not be formed (name is probably not unique)" );
    else {
      std::vector<size_t>  elementsIdx( elements_contiguous_subset.begin(), elements_contiguous_subset.end() );
      // accumulating the subregion
      (*it.first).second.AccumulateByNumber( gref.ElementsBegin(), gref.ElementsEnd(), elementsIdx );

      // copy all values of region properties from parent to child region
      (*it.first).second.LVS( gref.LVS() );
    }

    // subtracting the elements that constitute the new group from the remaining element list
    for ( typename std::set<size_t>::const_iterator
          sit = elements_contiguous_subset.begin(); sit != elements_contiguous_subset.end(); ++sit )
      elements.erase( (*sit) );

    // computing the next subset
    if ( elements.empty() ) break;
    else floodFillViaIndexes( gref, (*elements.begin()), elements_contiguous_subset );

    n_subgroups++;
  }

  if ( IsUnique( group ) && !(strncmp( group, "Model", NAME_STRING ) == 0) )
    RemoveRegion( group, false );

  return n_subgroups;

} // end partitionRegionIntoContiguousSubRegionsByIdx (by numbers)






/**
Removes subregions of the region identified by name. The subregions are
defined as regions that have the same name as the aforementioned region,
but with numbers appended, e.g., 'fractures' and 'fractures1'.

If the user supplies the string 'all subgroups' as method argument
all subgroups in the current Model object will be removed.

@param group The name of the master region, the subregions of which shall be removed.

@section application Application

To save memory by removing excessive region subdivisions as can be created
via the method PartitionRegionIntoContiguousSubRegions().

@section messages Messages

The method will report the names of the subregions that were removed.

@return the method will return the number of subregions that it removed.

@attention the master region is not touched or recreated; to achieve this use MergeRegions() before

@todo (1) SKM: logic of this method seems to be broken and it does not always work. Refactor!
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::RemoveRegionPartitionsFor( const char* group )
{
  std::string            target( group );
  std::set<std::string>  group_names, groups_to_remove;

  // 1. making a set of all region names
  for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
        grit = UniqueRegionsBegin(); grit != UniqueRegionsEnd(); grit++ )
    group_names.insert( (*grit).first );

  // 2. For all regions whose name does not contain any numbers,
  //    find subgroups identified by numbers attached to their names
  for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
        grit = UniqueRegionsBegin(); grit != UniqueRegionsEnd(); grit++ )
    if ( (*grit).first.find( target ) == std::string::npos )
    {
      // checking whether the name contains a number
      bool hasnumber = false;
      for ( std::string::const_iterator sit = (*grit).first.begin(); sit != (*grit).first.end(); sit++ )
        if ( isdigit( *sit ) ) { hasnumber = true; break; }
      // if not it is assumed that this is a primary group
      if ( !hasnumber ) {
        // and subgroups are searched for in the group_name set
        for ( std::set<std::string>::const_iterator it = group_names.begin(); it != group_names.end(); it++ )
          // if the groupname contains the search string, this is a subgroup to be deleted
          if ( (*it).find( target ) != std::string::npos )
            groups_to_remove.insert( (*it) );

        // the original group, however is kept by removing its name from the deletion list
        groups_to_remove.erase( (*grit).first );
      }
    }


  // Any region can only be unique or non unique. If no unique region was found the non-unique ones are searched
  if ( groups_to_remove.empty() ) {
    for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
          grit = RegionsBegin(); grit != RegionsEnd(); grit++ )
      group_names.insert( (*grit).first );

    // 2. For all groups whose names do not contain any numbers,
    //    find subgroups identified by numbers attached to their names
    for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
          grit = RegionsBegin(); grit != RegionsEnd(); grit++ )
      if ( (*grit).first.find( target ) == std::string::npos )
      {
        // checking whether the name contains a number
        bool hasnumber = false;
        for ( std::string::const_iterator sit = (*grit).first.begin(); sit != (*grit).first.end(); sit++ )
          if ( isdigit( *sit ) ) { hasnumber = true; break; }
        // if not it is assumed that this is a primary group
        if ( !hasnumber ) {
          // and subgroups are searched for in the group_name set
          for ( std::set<std::string>::const_iterator it = group_names.begin(); it != group_names.end(); it++ )
            // if the groupname contains the search string, this is a subgroup to be deleted
            if ( (*it).find( target ) != std::string::npos )
              groups_to_remove.insert( (*it) );

          // the original group, however is kept by removing its name from the deletion list
          groups_to_remove.erase( (*grit).first );
        }
      }
  }

  // removing the original region from the removal list, just in case it is contained therein
  groups_to_remove.erase( target );

  // 3. removing the subgroups and extra groups
  size_t  groups_removed( 0 );
  if ( !groups_to_remove.empty() ) std::cout << "\nRegionsInterface<dim,REGION_COMPLEX>::RemoveRegionPartitionsFor: removing region(s): ";
  for ( std::set<std::string>::const_iterator it = groups_to_remove.begin(); it != groups_to_remove.end(); it++ ) {
    std::cout << "'" << (*it) << "' ";
    RemoveRegion( (*it).c_str(), false );
    groups_removed++;
  }
  if ( !groups_to_remove.empty() ) std::cout << std::endl << std::endl;

  return groups_removed;

} // end RemoveRegionPartitionsFor






/**
FormAndAddRectangularRegion() builds a non-unique group from those elements
whose node coordinates lie within a user-specified bounding box. This box is
defined by its lower left and upper right corners.

When a new group is formed a group internal flag will be assigned to
each node, constraint point and element.
There are two group-internal object flags, PLAIN and BOUNDARY. When a
new group is formed, the Region method IdentifyBoundaryAs() assigns the
group-internal object flags, depending on whether the nodes, constraint
points or elements in the group lie at the group boundary or inside of
the group. Elements are assigned a boundary flag if at least one of
their faces coincides with the group boundary.

@section arguments Input Arguments

Assigns the finite elements that are located in a rectangular region that
is given by the double64 arguments (x = horizontal, y = vertical) as
members of a new group with the name char* (first argument).

@section implementation Implementation

The method works only for two-dimensional models and it calls the method
AccumulateRectangularRegion() of the Region class.

@section application Application

With FormAndAddRectangularRegion() regions to monitor changing properties
in a model can be defined at exact positions. Using property minimum
and maximum values in such regions, one can also measure property
gradients at positions of interest.

@section messages Messages

If the group cannot be formed because no elements in the desired
coordinate range can be found or because another group with the same
name already exists, the method will terminate the program, by reporting
a fatal error.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::FormRectangularRegion( const char* groupname,
                                                                  const Point<dim>& min_xyz,
                                                                  const Point<dim>& max_xyz )
{
  std::string output_region( groupname );
  if ( ContainsRegion( groupname ) ) {
    std::cout << "\nModel<" << dim << ">::FormRectangularRegion: WARNING: region '" << groupname;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) );

  // the 'bool' member of pair indicates whether insertion into map worked or not
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
    it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, regionComplex->Database() ) ) );
  if ( it.second )
  {
    auto& meshMgr = regionComplex->Mesh();
    (*it.first).second.AccumulateRectangularRegion( meshMgr, min_xyz, max_xyz );

    // removing the group if it contains no elements
    if ( (*it.first).second.Elements() == 0U ) {
      groupMap_.erase( it.first );
      csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRectangularRegion",
                         "Region could not be formed", output_region.c_str() );
      return false;
    }
  }
  else {
    csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRectangularRegion",
                       "Region could not be formed. Does this region already exist?", output_region.c_str() );
    return false;
  }

  // assigning new region name
  (*it.first).second.Name( groupname );

  return true;

} // end FormRectangularRegion 



/**
Makes a copy of an existing group and stores it under a new name.
NB: The new group will not be unique.

@param existing_group The names of the existing and the new group that shall be created.

@attention since this region will be a copy of an existing one, it will be non-unique

TODO: check whether this copying process is as fast as could be.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void  RegionInterface<dim, REGION_COMPLEX>::CopyRegion( const char* existing_group, const char* new_copied_group, bool unique_group )
{
  if ( !ContainsRegion( existing_group ) )
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::CopyRegion",
                           "region to make a copy of could not be found",
                           existing_group );

  std::string output_region( new_copied_group );
  if ( ContainsRegion( new_copied_group ) ) {
    std::cout << "\nModel<" << dim << ">::CopyRegion: WARNING: region '" << new_copied_group;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  const csmp::Region<dim>&  gr_ref( Region( existing_group ) );

  // irrespective of whether the original region was unique or non-unique its copy
  // will not be unique because it overlaps with the original region
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
    it = (unique_group) ? uniqueGroupMap_.insert( make_pair( output_region, csmp::Region<dim>( gr_ref ) ) ) :
    groupMap_.insert( make_pair( output_region, csmp::Region<dim>( gr_ref ) ) );
  if ( !it.second )
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::CopyRegion",
                           "region to copy to- could not be formed",
                           output_region.c_str() );

  // assigning new region name
  (*it.first).second.Name( new_copied_group );

} // end copyRegion





/**
Adds the first group to the second group so that the latter contains
both groups after the operation is complete.

@param group_to_add The names of the first group that shall be added to the second group.

@section messages Messages

ERRORs are reported if either the first or the second group does
not exist.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::AssimilateRegion( const char* group_to_add, const char* group_to_be_added_to )
{
  // finding the group in the group list
  try {
    const csmp::Region<dim>&  gref_to_add( Region( group_to_add ) );
    csmp::Region<dim>&        gref_to_be_added_to( Region( group_to_be_added_to ) );

    gref_to_be_added_to.Add( gref_to_add );
  }
  catch ( Exception& e ) {
    std::cout << "\nRegionsInterface<dim,REGION_COMPLEX>::AssimilateRegion: nothing was done; handled exception: " << std::endl;
    e.Out();
  }

} // end AssimilateRegion







/** Merges regions supplied as argument set into a single region with a new name.

@param input_groups is a set of unique names of the regions that
shall be merged into the
@param ensemble_group result region (second argument).

@section application Application

To join dynamically created regions during a simulation.

@section messages Messages

If the set of regions is empty or if one of the specified regions does not
exist, an error or an info message is reported.

@todo (3) Rewrite so that element ids are not longer required
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::MergeRegions( const std::set<std::string>& input_groups, const char* ensemble_group )
{
  std::string output_region( ensemble_group );

  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  if ( ContainsRegion( ensemble_group ) ) {
    csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                       ensemble_group, "output region already exists, adding an underscore to its name." );
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  if ( input_groups.empty() ) {
    csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                       "No input regions were specified; merge could not be performed; no new region.",
                       ensemble_group );
    return;
  }

  typename std::map<std::string, csmp::Region<dim> >::const_iterator  iter;

  // renumbering elements and nodes of model
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) );
  auto& meshMgr = regionComplex->Mesh();
  regionComplex->UpdateIndices();

  // collecting element indexes from input groups into set for output
  std::vector<size_t>  element_ids;
  for ( std::set<std::string>::const_iterator
        it = input_groups.begin(); it != input_groups.end(); it++ ) {
    // finding the group in the group list
    if ( (iter = groupMap_.find( *it )) != groupMap_.end() or
         (iter = uniqueGroupMap_.find( *it )) != uniqueGroupMap_.end() ) {
      //  outputting the ids of the member elements of the group
      element_ids.reserve( element_ids.size() + (*iter).second.Elements() );
      for ( typename std::vector<csmp::Element<dim>*>::const_iterator
            eit = (*iter).second.ElementsBegin(); eit != (*iter).second.ElementsEnd(); eit++ )
        element_ids.push_back( (*eit)->Idx() );
    }
    else csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                            (*it).c_str(), "region does not exist and was therefore not considered." );
  }

  if ( !element_ids.empty() ) {
    // making a non-unique new region
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region,
                             static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( !it.second )
      throw csmp::Exception( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                             output_region.c_str(), "region could not be formed." );

    else (*it.first).second.AccumulateByNumber( meshMgr, element_ids );
  }
  else
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions",
                           "No elements in target list; merged region could not be build",
                           output_region.c_str() );

  std::cout << "\nRegionsInterface<dim,REGION_COMPLEX>::MergeRegions: new region '" << output_region << "' formed successfully.\n";

} // end MergeRegions






template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::MergeRegions( const char* region_name_tag, const char* ensemble_group )
{
  std::set<std::string> regions;
  for ( regionConstIterator rit( UniqueRegionsBegin() ); rit != UniqueRegionsEnd(); ++rit )
    if ( rit->first.find( region_name_tag ) != std::string::npos )
      regions.insert( rit->first );
  if ( !regions.empty() )
    MergeRegions( regions, ensemble_group );
  return regions.size();
}







/**
If the first region a contains all the elements of the second region
b, the method will return the boolean variable 'true'.

@param groupa  a first group which is tested for whether it
contains all elements of a second group.

@return The method return either 'true' or 'false' depending on whether the first
group contains the second group or not.

@section application Application

RegionIncludes() is one of a set of boolean algebraic operations which are
supported for Region objects.

@section messages Messages

An error will be reported, if one of the evaluated groups does not
exist.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionIncludes( const char* groupa, const char* groupb ) const
{
  try {
    // finding the groups in the group list
    const csmp::Region<dim>&  a = Region( groupa );
    const csmp::Region<dim>&  b = Region( groupb );
    return a.Includes( b );
  }
  catch ( Exception& e ) {
    std::cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionIncludes: returning false; handled exception: " << std::endl;
    e.Out();
  }

  return false;

} // end RegionIncludes







/**
RegionUnion() will create a new region, defined as the union of two existing
groups.

@section arguments Input Arguments

The names of the two existing groups are supplied as first and second
method arguments. The third argument specifies the name of the new group
which will contain both the first and the second group.

@section application Application

To build groups on the basis of complex criteria, for instance, one could
form the union of a group which includes all elements that are hotter than
670oC with a group representing a granite melt in the model.

@section messages Messages

If either one of the groups is empty or does not exist or if the target
group cannot be formed an error is reported.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionUnion( const char* groupa, const char* groupb,
                                                        const char* groupunion )
{
  std::string output_region( groupunion );
  if ( ContainsRegion( groupunion ) ) {
    std::cout << "\nModel<" << dim << ">::RegionUnion: WARNING: region '" << groupunion;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  // finding the groups in the group list
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding new group
    // the 'bool' member of pair indicates whether insertion into map worked or not
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second )
    {
      groupUnion( itera, iterb, (*it.first).second );

      // removing the group if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Elements() == 0U ) {
        groupMap_.erase( it.first );
        return false;
      }
      // assigning new region name
      (*it.first).second.Name( groupunion );
    }
    else {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::RegionUnion",
                         "union of regions cannot be build", output_region.c_str() );
      return false;
    }

    return true;
  }
  catch ( Exception& e ) {
    std::cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionUnion: returning false; handled exception: " << std::endl;
    e.Out();
  }

  return false;

} // end RegionUnion








/**
RegionIntersection() finds those finite-elements which belong both to a
first and a second group and forms a third group from them.

Method returns true if there is an intersection and false if none
can be found.

@section arguments Input Arguments

The names of the two existing groups are supplied as first and second
method arguments. The third argument specifies the name of the new group
which will contain some members of the first and of the second group.

@section application Application

RegionIntersection() allows the user to define groups on the basis of
multiple property criteria. For instance, the user may form a group of
elements which represent granite above a temperature of 500oC. This would
be achieved by intersecting a 'granite' group with an 'above-500oC'
group.

@section messages Messages

If either the first or the second group does not exist or if the
desired intersection group would have zero elements, an error is
reported.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionIntersection( const char* groupa,
                                                               const char* groupb,
                                                               const char* groupintersection )
{
  std::string output_region( groupintersection );
  if ( ContainsRegion( groupintersection ) ) {
    std::cout << "\nModel<" << dim << ">::RegionIntersection: WARNING: region '" << groupintersection;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  // finding the groups in the group list
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding new group
    // the 'bool' member of pair indicates whether insertion into map worked or not
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second ) {
      intersection( itera, iterb, (*it.first).second );
      // removing the group if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Elements() == 0U ) {
        groupMap_.erase( it.first );
        return false;
      }
      // assigning new region name
      (*it.first).second.Name( groupintersection );
    }
    else return false;

    return true;
  }
  catch ( Exception& e ) {
    std::cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionIntersection: returning false; handled exception: " << std::endl;
    e.Out();
  }

  return false;

} // end RegionIntersection







/**
RegionDifference() forms a new group which will contain those finite
elements of the first group, which are not contained in the second
group.

If the difference does not exist false is returned.

@section arguments Input Arguments

The names of the two existing groups are supplied as first and second
method arguments. The third argument specifies the name of the new group.


@section application Application

To identify for instance what distinguishes a first from a second
group.

@section messages Messages

If either the first or the second group does not exist or if the
desired distinction group would have zero elements, an error is
reported.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionDifference( const char* groupa, const char* groupb,
                                                             const char* groupdiff )
{
  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  std::string output_region( groupdiff );
  if ( ContainsRegion( groupdiff ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RegionDifference: '", groupdiff,
                       "' already exists, adding an underscore at end of name." );
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  // finding the regions in the map
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding the new region
    // the 'bool' member of pair indicates whether insertion into map worked or not
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region,
                             static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second )
    {
      size_t elements_of_new_region = difference( itera, iterb, (*it.first).second );
      if ( elements_of_new_region == 0 )
        csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RegionDifference:", groupdiff,
                           "to be formed. No elements found that belong only to one of the 2 input regions." );

      // removing the group if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Elements() == 0 ) {
        groupMap_.erase( it.first );
        return false;
      }

      // assigning new region name
      (*it.first).second.Name( groupdiff );
    }

    return true;
  }
  catch ( Exception& e ) {
    std::cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionDifference: returning false; handled exception: " << std::endl;
    e.Out();
  }

  return false;

} // end RegionDifference









/**
RegionSymmetricDifference() forms a new region which will contain those finite
elements of a first and a second region, which are not contained in both
groups.

@section arguments Input Arguments

The names of the two existing groups are supplied as first and second
method arguments. The third argument specifies the name of the new group.


@section application Application

To identify for instance what two groups do not have in common.

@section messages Messages

If either the first or the second group does not exist or if the
desired distinction group would have zero elements, an error is
reported.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionSymmetricDifference( const char* groupa,
                                                                      const char* groupb,
                                                                      const char* groupsymdiff )
{
  std::string output_region( groupsymdiff );
  if ( ContainsRegion( groupsymdiff ) ) {
    std::cout << "\nModel<" << dim << ">::RegionSymmetricDifference: WARNING: region '" << groupsymdiff;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  // finding the groups in the group list
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding new group
    // the 'bool' member of pair indicates whether insertion into map worked or not
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second )
    {
      symmetricDifference( itera, iterb, (*it.first).second );
      // removing the group if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Elements() == 0 ) {
        groupMap_.erase( it.first );
        return false;
      }
      // assigning new region name
      (*it.first).second.Name( groupsymdiff );
    }
    else return false;

    return true;
  }
  catch ( Exception& e ) {
    std::cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionSymmetricDifference: returning false; handled exception: " << std::endl;
    e.Out();
  }

  return false;

} // end RegionSymmetricDifference






/**
Attempts to create a lower-dimensional region between higher dimensional ones.

TODO: @todo deprecate, better approach avoids search over the entire region, but looks at the perimeter
faces only.

@author R. Manasipov (2014)
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionBetween( const char* group1, const char* group2, const char* region_between )
{
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  std::string  region_name_between( region_between );

  std::string region1( group1 );
  std::string region2( group2 );

  if ( (region1 == "Model") || (region2 == "Model") ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween", "Region 'Model' not eligible for RegionBetween." );
    return false;
  }
  if ( group1 == group2 ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween", "Provided Regions are identical." );
    return false;
  }


  if ( !regionComplex->IsUnique( group1 ) or !regionComplex->IsUnique( group2 ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween",
                       "This method is intended for the creation of layer between unique Regions" );

    // checking for a potential overlap of the regions, if the regions are non-unique
    if ( regionComplex->RegionIntersection( group1, group2, "groupintersection" ) ) {
      csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween",
                         "one of the supplied regions is not unique and they overlap",
                         "It was therefore impossible to insert a boundary" );
      regionComplex->RemoveRegion( "groupintersection", false );
      return false;
    }
  }

  const csmp::Region<dim>&  gref1( regionComplex->Region( group1 ) );
  const csmp::Region<dim>&  gref2( regionComplex->Region( group2 ) );

  // checking whether the two regions share some nodes (these will mark their common boundary)
  const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
  bool  succeeded( false );

  if ( shared_nodes == 0U )
    csmp_error.notice( WARNING, "Model<dim,REGION_COMPLEX>::RegionBetween",
                       "The regions of interest do not share any nodes; trying to create a region between" );


  // attempt to create a regular (Face-based) boundary
  else {
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
      it = uniqueGroupMap_.insert( std::make_pair( region_name_between, csmp::Region<dim>( region_name_between, regionComplex->Database() ) ) );
    if ( it.second )
    {
      std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::RegionBetween creating Region between " << group1;
      std::cout << " and " << group2 << std::endl;
      regionComplex->UpdateIndices();
      succeeded = (*it.first).second.CreateBetween( regionComplex->Mesh(), regionComplex->FE_Manager(), gref1, gref2 );
      regionComplex->UpdateIndices();
      // assigning new region name
      (*it.first).second.Name( region_between );
    }
    else throw csmp::Exception( INFO, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween",
                                static_cast<std::string>(region_name_between).c_str(),
                                "region already exists. Nothing was done." );
  }

  if ( succeeded )
    std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::RegionBetweeen created region between " << group1 << " and " << group2 << std::endl;
  else
    throw csmp::Exception( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetweeen", "Creating Region failed!" );

  return true;
}



/**
Removes the elements shared with 'region-to-subtract' from the current (non-unique) region.

@note if the target region is unique, this implies that it does not overlap with any other region
except for non-unique regions.

@todo check (for volumetric regions) whether the removal of these elements has thus-far undiscovered consequences

@author SKM
@date 30/3/2016
@test OK
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RemoveFromRegion( const char* region, const char* region_to_subtract )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion( region ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                       region, "does not exist; nothing was done." );
    return false;
  }
  if ( !ContainsRegion( region_to_subtract ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                       region_to_subtract, "does not exist; nothing was done." );
    return false;
  }

  // finding the elements that are shared among the 2 regions
  csmp::Region<dim>&        r1_ref( Region( region ) );
  const csmp::Region<dim>&  r2_ref( Region( region_to_subtract ) );
  const size_t r1_elements( r1_ref.Elements() );
  const size_t r2_elements( r2_ref.Elements() );

  // new_region1 = region1 - region2
  std::vector<csmp::Element<dim>*>  new_region1;
  new_region1.reserve( r1_ref.Elements() );
  // if the element is not contained in region-to-subtract, it is kept
  for ( auto it = r1_ref.ElementsBegin(); it != r1_ref.ElementsEnd(); ++it )
    if ( !r2_ref.Contains( *it ) )
      new_region1.push_back( *it );

  // trimming excess capacity from new_region1
  std::vector<csmp::Element<dim>*>( new_region1 ).swap( new_region1 );

  // rebuilding the decimated region
  r1_ref.ElementVector() = std::move( new_region1 );
  r1_ref.EstablishNeighborConnectivity( false ); // TODO: needed, but this connectivity should have been established long ago !
  r1_ref.CreateNodePointerVector();
  r1_ref.IdentifyPerimeter();

  // reporting
  std::cout << "\nRegionInterface::RemoveFromRegion: removed region'" << region_to_subtract << "' (" << r2_elements << ")";
  std::cout << " from region '" << region << "' (" << r1_elements << ").\n";
  std::cout << "\t" << r1_ref.Elements() << " elements remaining in '" << region << "'.\n";

  return true;

} // RemoveFromRegion




/**
This version removes the supplied elements (as pointed to) from the target region.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RemoveFromRegion( const char* region, const std::set<Element<dim>*>& elmt_set )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion( region ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                       region, "does not exist; nothing was done." );
    return false;
  }
  if ( elmt_set.empty() ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                       "supplied element set was empty; nothing was done." );
    return false;
  }

  // finding the elements that are shared among the 2 regions
  csmp::Region<dim>&  r1_ref( Region( region ) );
  const size_t r1_elements( r1_ref.Elements() );

  // new_region1 = region1 - region2
  std::vector<csmp::Element<dim>*>  new_region1;
  new_region1.reserve( r1_ref.Elements() );
  for ( auto it = r1_ref.ElementsBegin(); it != r1_ref.ElementsEnd(); ++it )
    // if the element is not contained supplied set, it is kept
    if ( elmt_set.find( *it ) == elmt_set.end() )
      new_region1.push_back( *it );

  // trimming excess capacity from new_region1
  std::vector<csmp::Element<dim>*>( new_region1 ).swap( new_region1 );

  // rebuilding the decimated region
  r1_ref.ElementVector() = std::move( new_region1 );
  r1_ref.EstablishNeighborConnectivity( false ); // TODO: needed, but this connectivity should have been established long ago !
  r1_ref.CreateNodePointerVector();
  r1_ref.IdentifyPerimeter();

  // reporting
  std::cout << "\nRegionInterface::RemoveFromRegion: removed " << r1_elements - r1_ref.Elements();
  std::cout << " elements from region '" << region << "'.\n";

  return true;

} // RemoveFromRegion





/**
Moves region to from the unique- to the non-unique regions map.

@return whether the operation was performed successfully (also true when the region is already
non unique.

@author SKM
@date 30/3/2016
@test OK
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::MoveToNonUniqueRegions( const char* unique_region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion( unique_region ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::MoveToNonUniqueRegions:",
                       unique_region, "does not exist; nothing was done." );
    return false;
  }
  if ( !IsUnique( unique_region ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::MoveToNonUniqueRegions:",
                       unique_region, "is already a non-unique region; nothing was done." );
    return true;
  }

  // performing the move
  typename std::map<std::string, csmp::Region<dim> >::iterator  iterRegion( uniqueGroupMap_.find( std::string( unique_region ) ) );
  groupMap_.insert( std::make_pair( (*iterRegion).first, std::move( (*iterRegion).second ) ) );
  uniqueGroupMap_.erase( std::string( unique_region ) );

  return true;

} // end MoveToNonUniqueRegions



/**
Numbers the unique regions of the models, labeling their elements with the region number
as "region identifier".

@parameter vector of region names to retrieve them from the integer keys
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::CountAndLabelRegions( const char* region_identifier, std::vector<std::string>& region_names )
{
  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );

  // setting element variable up to identify all unique regions - uniquely
  const std::string rvariable( region_identifier );
  if ( !regionComplex.Database().IsDefined( rvariable.c_str() ) )
    regionComplex.CreateProperty( rvariable.c_str(), "X", SCALAR, ELEMENT );

  // counting the regions and initialising them with the unique identifiers
  region_names.resize( UniqueRegions() );
  size_t regions( 0 );
  for ( auto it = UniqueRegionsBegin(); it != UniqueRegionsEnd(); it++ ) {
    (*it).second.InputPropertyValue( rvariable.c_str(), makeScalar( PLAIN, static_cast<double64>(regions) ) );
    region_names[regions] = (*it).first;
    regions++;
  }

  assert( regions == region_names.size() );
  return regions;

} // end countAndLabelRegions


/**
Finds the contact area between regions a and b, logging pairs of element pointers and face numbers; @return number of shared faces

Region a will be the inner region whose perimeter Element pointers will be the first in the pairs.

@author SKM
@date 18/3/2017
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::SharedPerimeterFaces( const char* region_a, const char* region_b,
                                                                   std::vector<std::tuple<Element<dim>*, ///< inner element
                                                                   Element<dim>*, ///< outer element
                                                                   size_t,        ///< inner face
                                                                   size_t> >&     ///< outer face
                                                                   shared ) const
{
  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !IsUnique( region_a ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::SharedPerimeterFaces:",
                       region_a, "is not a unique region; cannot proceed." );
    return true;
  }
  if ( !IsUnique( region_b ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::SharedPerimeterFaces:",
                       region_b, "is not a unique region; cannot proceed." );
    return true;
  }

  const csmp::Region<dim>& subdomain_a( regionComplex.Region( region_a ) );
  const csmp::Region<dim>& subdomain_b( regionComplex.Region( region_b ) );

  if ( !shared.empty() ) shared.clear();
  shared.reserve( subdomain_a.PerimeterElements() );

  // 1. searching for shared faces between the regions
  for ( size_t eid = subdomain_a.InteriorElements(); eid<subdomain_a.Elements(); ++eid )
    for ( size_t face = 0U; face<subdomain_a.PerimeterFaces( eid ); ++face ) {
      // checking whether the neighbor of the perimeter face is in region b
      // ------------------------------------------------------------------
      const size_t pface = subdomain_a.PerimeterFace( eid, face );
      Element<dim>* nptr = subdomain_a.E( eid )->Neighbor( pface );
      // avoiding searches for neighbors that do not exist because one is at the model boundary
      if ( nptr != nullptr and subdomain_b.IsPerimeterElement( nptr ) ) {
        // finding which face is the perimeter face in the neighbor element
        size_t opposite_pface( UINT_MAX );
        for ( size_t i = 0U; i<nptr->Faces(); ++i )
          if ( nptr->Neighbor( i ) == subdomain_a.E( eid ) ) {
            opposite_pface = i;
            break;
          }
        assert( opposite_pface < nptr->Faces() );
        // creating connection record
        std::tuple<Element<dim>*, Element<dim>*, size_t, size_t>
          connection( subdomain_a.E( eid ), nptr, pface, opposite_pface );
        // storing the information
        shared.emplace_back( connection );
      }
    }

  // return how many shared faces were found
  return shared.size();

} // end SharedPerimeterFaces



template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::RegionsOut() const
 {
     std::cout <<"\n\nRegionInterface<"<< dim <<",Region<Element>>::RegionsOut:\n";
     std::cout <<"\n\tUnique regions of model:\n";
     for ( auto rit=UniqueRegionsBegin(); rit!=UniqueRegionsEnd(); ++rit ) {
          std::cout <<"\t\t"<< (*rit).first;
          std::cout <<" "<< (*rit).second.Elements() <<" elements,";
          std::cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          std::pair<int32, int32> rdim = (*rit).second.ElementSpatialDimensions();
          std::cout <<", range of spatial dimensions: "<< rdim.first <<" to "<< rdim.second << std::endl;
       }
     std::cout <<"\n\tNon-unique regions of model:\n";
     for ( auto rit=RegionsBegin(); rit!=RegionsEnd(); ++rit ) {
          std::cout <<"\t\t"<< (*rit).first;
          std::cout <<" "<< (*rit).second.Elements() <<" elements,"<< (*rit).second.Volume();
          std::cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          std::pair<int32, int32> rdim = (*rit).second.ElementSpatialDimensions();
          std::cout <<", range of spatial dimensions: "<< rdim.first <<" to "<< rdim.second << std::endl;
       }
     std::cout << std::endl << std::endl;
 }




// EXPLICIT TEMPLATE INSTANTIATION FOR REGION_INTERFACE
template class RegionInterface<1U, Model>;
template class RegionInterface<2U, Model>;
template class RegionInterface<3U, Model>;


} // csmp
