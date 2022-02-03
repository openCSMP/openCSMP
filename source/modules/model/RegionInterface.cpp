#include "RegionInterface.h"
#include "Region.h"
#include "Model.h"
#include "PropertyConstraints.h"
#include "ModelTopology.h"
#include "MeshManagementUtilities.h"
#include "CSMP_highLevelUtilities.h"
#include "UnionFind.h"
#include "binaryReadWrite.h"
#include "ErrorHandler.h"


using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsBegin()
{ return uniqueGroupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsEnd()
{ return uniqueGroupMap_.end(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::RegionsBegin()
{ return groupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::RegionsEnd()
{ return groupMap_.end(); }


template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsBegin() const
{ return uniqueGroupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsEnd() const
{ return uniqueGroupMap_.end(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator RegionInterface<dim, REGION_COMPLEX>::RegionsBegin() const
{ return groupMap_.begin(); }

template<size_t dim, template<size_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator RegionInterface<dim, REGION_COMPLEX>::RegionsEnd() const
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




// TODO: remove contiguity requirement in the presence of SplitBoundaries
  /// checks that there is a model region and that it contains elements
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::HasValidModelRegion() const
 {
    // does the region model exist?
    if ( !this->ContainsRegion("Model") ) return false;
    
    // is it contiguous?
    if ( !this->IsContiguous("Model") ) return false;
    
    return true;
 }




/**

Takes all elements from the mesh manager and creates a region that contains them.

Uses the highest order of elements (dim == model dimension) to perform a flood fill to determine whether the region is contiguous.

@return the number of elements contained in the region "Model" and whether this region is contiguous (true) or discontiguous (false).

@author SKM
@date 4/9/21

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormModelRegion( bool is_unique )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
  // 0. initial checks
  // -----------------
  const string regionname("Model");
  if ( ContainsRegion( "Model" ) ) {
      csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                         regionname, "'Model' region already exists; nothing was done." );
      return false;
    }

  // 1. building the 'Model' region
  // ------------------------------
  std::pair<typename map<string,csmp::Region<dim> >::iterator, bool>
    newRegion = (is_unique) ?
    uniqueGroupMap_.insert( make_pair( regionname, csmp::Region<dim>( regionname,
                            static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
    :
    groupMap_.insert( make_pair( regionname, csmp::Region<dim>( regionname,
                      static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  // if region was inserted successfully
  if ( newRegion.second ) {
       size_t elmts = (*newRegion.first).second.AccumulateAll( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh() );
       if ( elmts == 0 )
         csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                            regionname, "region could not be formed." );
    }
  else {
      csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                         regionname, "region could not be formed." );
      return 0U;
    }

  // checking that all elements and nodes were discovered
  if ( (*newRegion.first).second.Elements() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() or
       (*newRegion.first).second.Nodes() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Nodes() )
    {
       csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                          regionname, "'Model' not all elements were incorporated into the new 'Model' region." );
    }

  return (*newRegion.first).second.Elements();
  
} // end FormModelRegion




/**
      Forms unique regions called MATERIAL1..n from the material IDs assigned to the elements.
      
      @attention this method assumes that the unique numbers of elements, faces, and interfaces via the MeshManage
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromMaterialIDs( bool reestablishNeighborConnectivity )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
    MeshManager<dim>&    mesh = regionComplex->Mesh();
    size_t               n_regions(0U);

    // checking the existence of valid material ID values (mapping the element ids)
    map<int32_t,vector<Element<dim>*> >  mtrl_ids;
    const typename plf::colony<Element<dim>>::iterator elmts_end(mesh.ElementsEnd());
    typename plf::colony<Element<dim>>::iterator       eit(mesh.ElementsBegin());
    string                                             region_name("undefined");
    
    // collecting the elements that make up the different materials of the model
    while ( eit != elmts_end ) {
         pair<typename map<int32_t,vector<Element<dim>*> >::iterator,bool>
           it = mtrl_ids.insert( make_pair( (*eit).Material_ID(),
                                             vector<Element<dim>*>{ 1, &(*eit) } ) );
         if ( it.second == false )
           (*it.first).second.push_back( &(*eit) );
         eit++;
      }

    if ( mtrl_ids.size() <= 1 ) {
         csmp_error.notice( ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromMaterialIDs:",
                           "Material IDs do not appear to have been initialized; nothing was done.");
         return 0u;
      }

    // creating the regions from the material identifiers
    for ( const auto& mit : mtrl_ids )
      {
         region_name = "MATERIAL" + to_string( mit.first );
         pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
           it = uniqueGroupMap_.insert( make_pair( region_name, csmp::Region<dim>( region_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
         // if the region was successfully inserted
         if ( it.second )
           {
             // retrieving the elements by their IDs and assigning them  to the region
             (*it.first).second.Accumulate( mit.second.begin(), mit.second.end() );

             // removing the group if it contains no elements
             if ( (*it.first).second.Elements() == 0U ) {
                 uniqueGroupMap_.erase( it.first );
                 csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromMaterialIDs",
                                    region_name, "could not be formed" );
               }
             else n_regions++;
          }
        else
          throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromMaterialIDs",
                                 region_name, "could not be formed. Does this region already exist?" );
        
      } // end materials loop

    return n_regions;
    
 } // end FormRegionsFromMaterialIDs




/**
Removes named Region object and its elements and nodes.

@param  regionName The name of the group object which shall be removed.

@section messages Messages

@note: this can trigger removal of elements

@note If the region which shall be removed does not exist, the method reports a warning.

@attention to delete underlying elements, the MeshManager had to be called upon.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::RemoveRegion( const char* regionName )
{
  // check whether region exists (should be a notice only, nothrow)
  if ( !ContainsRegion( regionName ) )
    throw csmp::Exception( WARNING,
                           "RegionsInterface<dim,Model>::RemoveRegion",
                           "region did not exist: ",
                           regionName );

  // finding the region in the corresponding map
  typename map<string, csmp::Region<dim> >::iterator
    iterRegion( groupMap_.find( string( regionName ) ) ),
    iterUniqueRegion( uniqueGroupMap_.find( string( regionName ) ) );

  // if the region was found in the respective map, it is erased
  if ( iterRegion != groupMap_.end() )
    groupMap_.erase( std::string( regionName ) );
  if ( iterUniqueRegion != uniqueGroupMap_.end() )
    uniqueGroupMap_.erase( std::string( regionName ) );

} // end RemoveRegion


// DEBUG - check element vector for duplicates (OK)
//set<const csmp::Element<dim>*> elmts( subdomain.ElementsBegin(), subdomain.ElementsEnd() );
//assert( elmts.size() == subdomain.Elements() );
// DEBUG - check element vector for dead elements (OK - no corrupt elements)
//cerr <<"\nRegionInterface<dim, REGION_COMPLEX>::RemoveRegion: "<< regionName <<"\n";
//for ( auto it=subdomain.ElementsBegin(); it!=subdomain.ElementsEnd(); ++it )
//  cerr <<" "<< parseFiniteElementType( (*it)->FE_Type() ) <<":"<< (*it)->Idx();
//cerr << endl;
// DEBUG - check MeshManager vector for dead elements (OK - no corrupt elements)
//cerr <<"\nRegionInterface<dim, REGION_COMPLEX>::RemoveRegion: "<< regionName <<"\n";
//for ( auto it=meshMgr.ElementsBegin(); it!=meshMgr.ElementsEnd(); ++it )
//  cerr <<" "<< parseFiniteElementType( (*it).FE_Type() ) <<":"<< (*it).Idx();
//cerr << endl;

 







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
void RegionInterface<dim, REGION_COMPLEX>::OutputRegionsToBinary( const char* file_name ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( regionComplex.Database() );

  std::string bin_file( file_name );
  std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::binary );
  if ( !fp.is_open() ) {
      csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::OutputRegionsToBinary:",
                         bin_file, "file could not be opened; nothing was done." );
      return;
    }

  // -----------------------------
  // 1. File header
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "REGFHEDR" );

    std::string heading( "RegionInterface::OutputRegionsToBinary: " );
    heading += "region information for Model '";
    heading += regionComplex.Name();
    heading += "' written to file: ";
    heading += bin_file;
    heading += "'.";

    binaryFileWrite( fp, heading.c_str() );
  }

  std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::OutputRegionsToBinary: regions written to binary file: ";
  std::cout << "\n\n\tUnique regions: ";

  // -----------------------------
  // 2. writing the unique regions
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "UNIQREGN" );

    int64_t  records = this->UniqueRegions();
    
    // writing number of unique regions
    fp.write( reinterpret_cast<const char*>(&records), sizeof( int64_t  ) );

    for ( typename std::map<std::string, csmp::Region<dim> >::const_iterator
          git = UniqueRegionsBegin(); git != UniqueRegionsEnd(); ++git )
      {
        BinaryFileSectionWrite hdr( fp, "ONE_REGN" );
        std::cout << (*git).first << " ";
        // refer to ModelSubDomain for following method
        (*git).second.WriteDomainIndexesToBinaryFile( fp );
        domainVariablesOut( fp, (*git).second, database );
      }
  }

  // -----------------------------
  // 3. writing non-unique regions
  // -----------------------------
  {
    cout << "\n\n\tNon-unique regions overlapping unique ones and potentially each other: ";
    BinaryFileSectionWrite hdr( fp, "NONUREGN" );

    int64_t  records = this->Regions() - this->UniqueRegions();

    fp.write( reinterpret_cast<const char*>(&records), sizeof( int64_t  ) );

    for ( auto git = RegionsBegin(); git != RegionsEnd(); git++ )
      {
        BinaryFileSectionWrite hdr( fp, "ONE_REGN" );
        std::cout << (*git).first << " ";
        (*git).second.WriteDomainIndexesToBinaryFile( fp );
        domainVariablesOut( fp, (*git).second, database );
      }
    std::cout << std::endl;
  }

  // -----------------------------
  // 3. writing model variables
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "MODLVARS" );
    // writing the Model variables here
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
  std::cout << "\nRegionInterface<" << dim << ",REGION_COMPLEX>::OutputRegionsToBinary: file '";
  std::cout << bin_file << "' has been successfully written.\n";

} // end OutputRegionsToBinary






/**
Reads all regions stored by OutputRegionsToBinary() into Interface,
constructing them using the SubDomainInfo data, thereby avoiding the costly
re-initialisation.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::InputRegionsFromBinary( const char* file_name,
                                                                   const set<string>& subset_variables )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );
  MeshManager<dim>& mesh( static_cast<REGION_COMPLEX<dim>& >(*this).Mesh() );

  std::string bin_file( file_name );
  std::fstream fp( bin_file.c_str(), std::ios::in | std::ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                       bin_file, "file could not be opened; nothing was done." );
    return;
  }

  // -----------------------------
  // 1. Reading file header
  // -----------------------------
    {
      BinaryFileSectionRead hdr( fp, "REGFHEDR" );
      char                  text[INFO_STRING];
      binaryFileRead( fp, text );
      std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary: Reading file header:\n\t" << text << std::endl;
    }
  std::cout << "\n\timporting the regions: ";

  const PropertyDatabase<dim>& database( static_cast<const REGION_COMPLEX<dim>& >(*this).Database() );

  // -----------------------------
  // 2. unique regions
  // -----------------------------
  std::cout << "\n\tunique regions: ";
    {
      BinaryFileSectionRead hdr( fp, "UNIQREGN" );
      // getting number of unique region records from file
      int64_t   records( 0 );  // region records
      fp.read( reinterpret_cast<char*>(&records), sizeof( int64_t  ) );
      if ( records > 0 )
          // reading the regions sequentially
          for ( size_t i = 0U; i<records; i++ )
            {
              BinaryFileSectionRead hdr( fp, "ONE_REGN" );
              // reading name and element indices for each unique region
              SubDomainInfo  info;
              readDomainIndexesFromBinaryFile( dim, fp, info );
              // reconstruct the region
              std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
                it = uniqueGroupMap_.insert( std::make_pair( info.name, csmp::Region<dim>( database, mesh, info ) ) );

              if ( !it.second )
                throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                                       info.name, "Region could not be formed; issue with binary file." );

              std::cout << "\n\t\t'" << (*it.first).first << "'(" << (*it.first).second.Elements() << " elements).";
              
              // reading variables stored on the regions
              if ( subset_variables.empty() ) domainVariablesIn( fp, (*it.first).second, database );
              else selectedDomainVariablesIn( fp, (*it.first).second, database, subset_variables );
            }
       else csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                               bin_file, "does not contain any unique region descriptions; no regions were initialised." );
    }

  // -----------------------------
  // 3. reading non-unique regions
  // -----------------------------
  std::cout << "\n\tnon-unique regions: ";
  {
    BinaryFileSectionRead hdr( fp, "NONUREGN" );
    // getting number of non-unique region records from file
    int64_t   records( 0 );
    fp.read( reinterpret_cast<char*>(&records), sizeof( int64_t  ) );
    if ( records > 0 )
      // reading the regions sequentially
      for ( size_t i = 0U; i<records; i++ )
        {
          BinaryFileSectionRead hdr( fp, "ONE_REGN" );
          // reading name and element indices for each unique region
          SubDomainInfo  info;
          readDomainIndexesFromBinaryFile( dim, fp, info );
          // if the region info record is not empty the region is reconstructed
          pair<typename map<string, csmp::Region<dim> >::iterator, bool>
            it = groupMap_.insert( std::make_pair( info.name, csmp::Region<dim>( database, mesh, info ) ) );
          if ( !info.interior_elmts.empty() ) {
              if ( !it.second )
                throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                                       info.name, "Region could not be formed; issue with binary file." );

              cout << "\n\t\t'" << (*it.first).first << "'(" << (*it.first).second.Elements() << " elements).";

              // reading variables stored on the regions
              if ( subset_variables.empty() ) domainVariablesIn( fp, (*it.first).second, database );
              else selectedDomainVariablesIn( fp, (*it.first).second, database, subset_variables );
            }
        }
    else csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                            bin_file, "does not contain any non-unique region descriptions; no regions were initialised." );
  }

  // -----------------------------
  // 4. reading model vars
  // -----------------------------
  {
    BinaryFileSectionRead hdr( fp, "MODLVARS" );
    // reading variables placed on the model here
    domainVariablesIn( fp, regionComplex, database );
  }

  // -----------------------------
  // 5. file footer
  // -----------------------------
  BinaryFileSectionRead hdr( fp, "REGFFOTR" );

  // -----------------------------
  // 6. cleanup
  // -----------------------------

  fp.close();
  std::cout << "\n\nRegionInterface<" << dim << ",REGION_COMPLEX>::InputRegionsFromBinary: file '";
  std::cout << bin_file << "' has been read successfully.\n";

} // end InputRegionsFromBinary











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
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues( const char* prop,
                                                                            std::set<std::string>& group_names )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.notice( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  if ( !group_names.empty() )
    group_names.erase( group_names.begin(), group_names.end() );

  csmp::Index  prop_key = regionComplex->Database().StorageKey( prop );

  if ( prop_key.type != SCALAR )
    throw csmp::Exception( INFO, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromPropertyValues",
                           "method can only be applied to SCALAR variables",
                           "no assignments were made" );

  if ( prop_key.place != NODE && prop_key.place != ELEMENT && prop_key.place != ELEMENT_INTEGRATION_POINT )
    throw csmp::Exception( INFO, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromPropertyValues",
                           "method can only be applied to Node, Element, and Element-IP variables",
                           "no assignments were made" );

  map<double,string>  groups;

  // 1. Making a map with one entry for each property value
  // ------------------------------------------------------
  const csmp::Region<dim>& model_domain(Region("Model"));

  switch ( prop_key.place )
  {
    case NODE:
      for ( auto& nit : model_domain.NodeVector() ) {
          double sc = nit->Read( prop_key );
          groups[sc] = "undefined";
        }
      break;
    case ELEMENT_INTEGRATION_POINT:
      for ( auto eit : model_domain.CellVector() ) {
        for ( size_t i = 0U; i < eit->IntegrationPoints(); i++ )
          {
            double sc = eit->Read( i, prop_key );
            groups[sc] = "undefined";
          }
      }
      break;
    case ELEMENT:
      for ( auto eit : model_domain.CellVector() ) {
          double sc = eit->Read( prop_key );
          groups[sc] = "undefined";
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

  for ( typename std::map<double, std::string>::iterator
        it = groups.begin(); it != groups.end(); it++ )
  {
    std::cout << "\nProperty value: " << (*it).first;
    sprintf( num, "%lu", group_idx++ );
    (*it).second = std::string( gname + num );

    if ( (*it).second != "undefined" )
      {
        assert( !ContainsRegion( (*it).second.c_str() ) );
        FormRegionFrom( (*it).second.c_str(), prop,
                        (*it).first - std::numeric_limits<double>::epsilon(),
                        (*it).first + std::numeric_limits<double>::epsilon(), true );

        group_names.insert( (*it).second );
      }
    else
      std::cout << "\nModel<" << dim << ">::FormRegionsFromPropertyValues: check predefined region name." << std::endl;
  }
  
  return groups.size();

} // end RegionsFromPropertyValues






/**
     Uses the region model, to form a region from its elements with the corresponding ID numbers.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionName, vector<size_t>& elmt_ids, bool unique/*=false */ )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.notice( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  csmp::Region<dim>& model_domain(Region("Model"));
  
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool> it;
  if ( unique ) {
        it = this->uniqueGroupMap_.insert( make_pair( regionName, csmp::Region<dim>( regionName, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
        if ( it.second )
          (*it.first).second.AccumulateByNumber( model_domain.ElementsBegin(), model_domain.ElementsEnd(), elmt_ids );
    }
  else {
        it = this->groupMap_.insert( std::make_pair( regionName, csmp::Region<dim>( regionName, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
        if ( it.second )
          (*it.first).second.AccumulateByNumber( model_domain.ElementsBegin(), model_domain.ElementsEnd(), elmt_ids );
    }

  return (*it.first).second.Elements();
}







  /// forms unique or non-unique region from range of elements; returns reference to it
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionname,
                                                             typename vector<Element<dim>*>::iterator first,
                                                             typename vector<Element<dim>*>::iterator last,
                                                             bool unique )
 {
    string output_region( regionname );
    if ( ContainsRegion( regionname ) ) {
      cerr << "\nModel<" << dim << ">::FormRegionFrom: WARNING: region '" << regionname;
      cerr << "' already exists, adding an underscore at end of name: ";
      output_region += "_";
      cout << output_region << endl;
    }

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // the 'bool' member of pair indicates whether insertion into map worked or not
    std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>  it;
    if ( unique )
      it = uniqueGroupMap_.insert( std::make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    else
      it = groupMap_.insert( std::make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

    if ( it.second ) {
        (*it.first).second.Accumulate( first, last );

        // removing the group if it contains no elements
        if ( (*it.first).second.Elements() == 0U ) {
          groupMap_.erase( it.first );
          csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                             "Region could not be formed", output_region.c_str() );
          return 0U;
        }
      }
    else {
        csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                           "Region could not be formed. ", output_region.c_str() );
        return 0U;
      }

    return (*it.first).second.Elements();
    
 } // end FormRegionFrom (range of element pointers)



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
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* groupname,
                                                             PropertyConstraints& constraints,
                                                             bool unique_group )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.notice( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  string output_region( groupname );
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

  if ( it.second )
    {
      REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
      constraints.InitializePropertyIndices( regionComplex->Database() );
      csmp::Region<dim>& model_domain(Region("Model"));
      (*it.first).second.AccumulateWithinRange( regionComplex->Mesh(), constraints );

      // removing the group if it contains no elements
      if ( (*it.first).second.Elements() == 0U ) {
        groupMap_.erase( it.first );
        csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                           "Region could not be formed", output_region.c_str() );
        return 0U;
      }
    }
  else {
    csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                       "Region could not be formed", output_region.c_str() );
    return 0U;
  }

  return (*it.first).second.Elements();

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
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* groupname,
                                                             const char* prop,
                                                             double min, double max,
                                                             bool unique_group )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.notice( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  string output_region( groupname );
  if ( ContainsRegion( groupname ) ) {
    std::cout << "\nModel<" << dim << ">::FormRegionFrom: WARNING: region '" << groupname;
    std::cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    std::cout << output_region << std::endl;
  }

  // the 'bool' member of pair indicates whether insertion into map worked or not
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>  it;
  if ( unique_group )
    it = uniqueGroupMap_.insert( std::make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
  else
    it = groupMap_.insert( std::make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  if ( it.second )
    {
      (*it.first).second.AccumulateWithinRange( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh(), prop, min, max );

      // removing the group if it contains no elements
      if ( (*it.first).second.Elements() == 0U ) {
        groupMap_.erase( it.first );
        csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                           "Region could not be formed", output_region.c_str() );
        return 0U;
      }
    }
  else {
    csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                       "Region could not be formed. ", output_region.c_str() );
    return 0U;
  }

  return (*it.first).second.Elements();

} // end FormRegionFrom











/**
Forms a Region from a set of region names by using the method MergeRegions.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionname, const std::set<std::string>& region_names )
{
    std::string output_region( regionname );
    if ( ContainsRegion( regionname ) ) {
      std::cout << "\nModel<" << dim << ">::FormRegionFrom (element numbers): WARNING: region '" << regionname;
      std::cout << "' already exists, adding an underscore at end of name: ";
      output_region += "_";
      std::cout << output_region << std::endl;
    }

    this->MergeRegions( region_names, output_region.c_str() );
    
    return this->Region(output_region.c_str()).Elements();
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
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* newRegionName, ElementComp<dim> const& elementComp, const char* hostRegion )
  {
    static_cast<REGION_COMPLEX<dim>*>(this)->UpdateIndices();
    csmp::Region<dim> const& rref( this->Region( hostRegion ) );
    rref.UpdateMemberIndexes();
    std::vector<size_t> elementIds;
    elementIds.reserve( rref.Elements() );
    for ( typename csmp::Region<dim>::SimplexContainer::const_iterator it( rref.ElementsBegin() ); it != rref.ElementsEnd(); ++it )
      if ( elementComp( (*it) ) )
        elementIds.push_back( (*it)->Idx() );
    std::vector<size_t>( elementIds ).swap( elementIds );
    if ( elementIds.empty() )
      return 0U;
    
    return this->FormRegionFrom( newRegionName, elementIds );
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

@todo for the debug version, put in a check that verifies that all element ids stored in the model topology are actually contained in the mesh.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionsFrom( const ModelTopology& topo )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. getting the names of the regions
  std::list<std::string> regions;
  topo.Out( regions );

  // 2. assigning the regions to groups in the Model
  std::cout << "\nRegionInterface<dim,REGION_COMPLEX>::FormRegionsFrom: Forming the regions: ";

  int32_t new_regions( 0U );
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

          // retrieving the elements by their IDs and assigning them  to the region
          (*it.first).second.AccumulateByNumber( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh(), element_ids );
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

  if ( new_regions == topo.ModelRegions() ) return regions.size();
  return 0U;

} // end FormRegionsFrom





/**
Tests whether the elements of the region are connected to each-other.

@attention This test cannot be performed if the region has elements of different spatial
dimensions since these are not interconnected. Therefore, this method returns false if
the region consists of elements from different spatial dimensions.

@note This method is based on the floofFill() algorithm implemented in CSMP. 

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::IsContiguous( const std::string& region_name ) const
{
  const csmp::Region<dim>& mref = Region( region_name );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  std::pair<int32_t, int32_t>  dimensionality = mref.ElementSpatialDimensions();
  if ( dimensionality.first > 1U ) {
      csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::IsContiguous:",
                         "method can determine contiguity only for regions which consist only of same spatial dimension elements; returned false." );
      return false;
    }

  set<Element<dim>*> contiguous_elmts;

  floodFill( const_cast<Element<dim>* const>(*mref.ElementsBegin()), contiguous_elmts );
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

      The method uses the element neighbor information to perform floodfills
      as required until all member elements are partitioned.

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
size_t  RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions( const char* group )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      if ( std::string("Model") == group ) {
         csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions:",
                           "this operation is not allowed for region 'Model' or the master region." );
            return 0U;
        }
     
      // renumbering elements and nodes of model
      csmp::Region<dim>&  gref(Region(group));
      const bool          unique_group(IsUnique(group));

      // getting a set of the element numbers of the target group
      set<Element<dim>*>  elements( gref.ElementsBegin(), gref.ElementsEnd() );

      // detecting via a flood-fill whether the group can be partitioned, else nothing is done
      set<Element<dim>*>  elements_contiguous_subset;
      floodFill( gref.E(0), elements_contiguous_subset );
      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( elements.size() <= elements_contiguous_subset.size() ) {
           std::cout <<"\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions: ";
           std::cout <<"region '"<< group <<"' is already contiguous, nothing was done."<< std::endl;
           return 0U;
        }

      // else partitions can be created
      std::string  group_name(group);
      std::string  subgroup_name;
      char         num[128];
      size_t       n_subgroups(1);
    
      // creating new contiguous group from the element subset
      while ( !elements.empty() )
        {
           // creating name of contiguous subgroup
           sprintf( num, "%lu", n_subgroups );
           subgroup_name = group_name + num;
           if ( n_subgroups == 1U ) {
                 std::cout <<"\nModel<"<< dim <<">::PartitionRegionIntoContiguousSubRegions: ";
                 std::cout <<"region '"<< group <<"' is divided into the subregion(s):\n";
             }
           std::cout <<"\t\t\t'"<< subgroup_name <<"'";
           std::cout <<" ("<< elements_contiguous_subset.size() <<" elmts)"<< std::endl;
         
           // creating either a unique or non-unique group depending on uniqueness of original region
           std::pair<typename map<string,csmp::Region<dim> >::iterator,bool>
             it = ( unique_group ) ? uniqueGroupMap_.insert( make_pair( subgroup_name, csmp::Region<dim>( subgroup_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database()) ) )
                                   : groupMap_.insert( make_pair( subgroup_name, csmp::Region<dim>( subgroup_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database()) ) );
           if ( !it.second )
             throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions",
                                    subgroup_name.c_str(), "region could not be formed (name is probably not unique)" );
           else {
                // accumulating the subregion
                (*it.first).second.Accumulate( elements_contiguous_subset.begin(),
                                               elements_contiguous_subset.end() );
                                               
                // copy all values of region properties from parent to child region
                (*it.first).second.LVS( gref.LVS() );
             }
         
           // subtracting the elements that constitute the new group from the remaining element list
           for ( typename std::set<Element<dim>*>::const_iterator
                 sit=elements_contiguous_subset.begin(); sit!=elements_contiguous_subset.end(); ++sit )
             elements.erase( (*sit) );

           // computing the next subset
           if ( elements.empty() ) break;
           else floodFill( (*elements.begin()), elements_contiguous_subset );
      
           n_subgroups++;
        }

      // if the region has been partitioned succesfully and its name is not model, it will be removed
      if ( IsUnique(group) )
        RemoveRegion( group );

      return n_subgroups;
    
   } // end partitionRegionIntoContiguousSubRegions






/**

ANDREW BROMAGE VERSION

Breaks non-contiguous Region into contiguous subregions that carry the name
of the original region, but have a number as a suffix to their name to distinguish each partition.
All elements and the original region are retained, but the original region is moved into the 
non-unique regions storage.

@return The method returns the number of subgregions that were created, if any.

@param group The name of the region that may be non-contiguous.

@section implementation Implementation

The method uses the union-find algorithm, using element neighbours to determine
which parts of it are contiguous.

@attention this method is implemented using UnionFind  as opposed to standard floodfill operations. This might cause performance issues.

@section application Application

To  partition regions  into contiguous subdomains such as fractures so that these 
can be processed one-by-one by various algorithms. 

@section messages Messages

The method will report if the region is  contiguous to start with. In this case no subregions will be created.

@attention it makes no sense to apply this method to region "Model" or complex non-unique regions.

TODO: legacy of Andrew Bromage, not sure whether it makes any sense, needs testing!
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions_Bromage( const char* group )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  const std::string  group_name( group );

  if ( group_name == "Model" ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions_Bromage:",
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

  vector<pair<size_t, Element<dim>*> > components;
  unionFind.Components( components );
  if ( components.size() <= 1 ) {
    std::cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions_Bromage: ";
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
  std::cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions_Bromage: ";
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
      throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions_Bromage",
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
    MoveToNonUniqueRegions( group );

  return subgroupNum;

} // end partitionRegionIntoContiguousSubRegions













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
    RemoveRegion( (*it).c_str() );
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
is given by the double arguments (x = horizontal, y = vertical) as
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
size_t RegionInterface<dim, REGION_COMPLEX>::FormRectangularRegion( const char* groupname,
                                                                    const Point<dim>& min_xyz,
                                                                    const Point<dim>& max_xyz )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.notice( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  std::string output_region( groupname );
  if ( ContainsRegion( groupname ) ) {
      std::cout << "\nModel<" << dim << ">::FormRectangularRegion: WARNING: region '" << groupname;
      std::cout << "' already exists, adding an underscore at end of name: ";
      output_region += "_";
      std::cout << output_region << std::endl;
    }

  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) );
  // the 'bool' member of pair indicates whether insertion into map worked or not
  std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
    it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, regionComplex->Database() ) ) );
  if ( it.second )
    {
       (*it.first).second.AccumulateRectangularRegion( regionComplex->Mesh(), min_xyz, max_xyz );

       // removing the group if it contains no elements
       if ( (*it.first).second.Elements() == 0U ) {
          groupMap_.erase( it.first );
          csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRectangularRegion",
                             "Region could not be formed", output_region.c_str() );
          return 0U;
        }
    }
  else {
    csmp_error.notice( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRectangularRegion",
                       "Region could not be formed. Does this region already exist?", output_region.c_str() );
    return 0U;
  }

  // assigning new region name
  (*it.first).second.Name( groupname );

  return (*it.first).second.Elements();

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

  // collecting element indexes from input groups into set for output
  std::vector<Element<dim>*>  element_ptrs;
  for ( auto it = input_groups.begin(); it != input_groups.end(); it++ ) {
    // finding the group in the group list
    if ( (iter = groupMap_.find( *it )) != groupMap_.end() or
         (iter = uniqueGroupMap_.find( *it )) != uniqueGroupMap_.end() ) {
      //  outputting the ids of the member elements of the group
      element_ptrs.reserve( element_ptrs.size() + (*iter).second.Elements() );
      for ( auto eit = (*iter).second.ElementsBegin(); eit != (*iter).second.ElementsEnd(); eit++ )
        element_ptrs.push_back( const_cast<Element<dim>*>(*eit) );
    }
    else csmp_error.notice( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                            (*it).c_str(), "region does not exist and was therefore not considered." );
  }

  if ( !element_ptrs.empty() )
    {
      // eliminating duplicate entries from pointer vector
      sort( element_ptrs.begin(), element_ptrs.end() );
      element_ptrs.erase( unique( element_ptrs.begin(), element_ptrs.end() ), element_ptrs.end() );
      
      // making a non-unique new region
      std::pair<typename std::map<std::string, csmp::Region<dim> >::iterator, bool>
        it = groupMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region,
                               static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
      if ( !it.second )
        throw csmp::Exception( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                               output_region.c_str(), "region could not be formed." );

      else (*it.first).second.Accumulate( element_ptrs.begin(), element_ptrs.end() );
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
size_t RegionInterface<dim, REGION_COMPLEX>::RegionBetween( const char* group1, const char* group2,
                                                            const char* region_between, int32_t material_id )
{
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  std::string  region_name_between( region_between );

  std::string region1( group1 );
  std::string region2( group2 );

  if ( (region1 == "Model") || (region2 == "Model") ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween", "Region 'Model' not eligible for RegionBetween." );
    return 0U;
  }
  if ( group1 == group2 ) {
    csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween", "Provided Regions are identical." );
    return 0U;
  }


  if ( !regionComplex->IsUnique( group1 ) or !regionComplex->IsUnique( group2 ) ) {
    csmp_error.notice( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween",
                       "This method is intended for the creation of layer between unique Regions" );

    // checking for a potential overlap of the regions, if the regions are non-unique
    if ( regionComplex->RegionIntersection( group1, group2, "groupintersection" ) ) {
      csmp_error.notice( ERROR, "RegionInterface<dim,REGION_COMPLEX>::RegionBetween",
                         "one of the supplied regions is not unique and they overlap",
                         "It was therefore impossible to insert a boundary" );
      // regionComplex->RemoveRegion( "groupintersection", false );
      return 0U;
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
        succeeded = (*it.first).second.CreateBetween( regionComplex->Mesh(), regionComplex->FE_Manager(), gref1, gref2, material_id );
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

  return Region(region_between).Elements();
  
} // end CreateBetween






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
  r1_ref.CellVector() = std::move( new_region1 );
// NOT AFFECTED  r1_ref.EstablishNeighborConnectivity( false ); // TODO: needed, but this connectivity should have been established long ago !
  r1_ref.CreateNodePointerVector();
  r1_ref.IdentifyPerimeter();

  // reporting
  std::cout << "\nRegionInterface::RemoveFromRegion: removed region'" << region_to_subtract << "' (" << r2_elements << ")";
  std::cout << " from region '" << region << "' (" << r1_elements << ").\n";
  std::cout << "\t" << r1_ref.Elements() << " elements remaining in '" << region << "'.\n";

  return true;

} // RemoveFromRegion




/**
    This version removes the pointers to the supplied elements from the target region.
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RemoveFromRegion( const char* region, const set<Element<dim>* const>& elmt_set )
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
  csmp::Region<dim>&   subdomain( Region(region) );
  vector<Element<dim>* const> elmts_to_remove( elmt_set.begin(), elmt_set.end() );
  
  const size_t elmts_removed = subdomain.RemoveRange( elmts_to_remove.begin(), elmts_to_remove.end() );
  
  // reporting
  std::cout << "\nRegionInterface::RemoveFromRegion: removed " << elmts_removed;
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
    (*it).second.InputPropertyValue( rvariable.c_str(), makeScalar( PLAIN, static_cast<double>(regions) ) );
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



/**
      if the mesh changed this brute-force method rebuild the node and element vectors of all regions
      TODO: find way to do this more selectively
      
      @attention the assumption is made the inter-element connectivity has was updated before
*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::RebuildRegions()
 {
    // since this region may now contain a different number of elements
    RemoveRegion("Model" );

     // rebuilding the region 'Model'
     const bool is_unique = ( distance(UniqueRegionsBegin(), UniqueRegionsEnd()) > 0 ) ? false : true;
     FormModelRegion( is_unique );

     for ( auto rit=RegionsBegin(); rit!=RegionsEnd(); ++rit ) {
           rit->second.CreateNodePointerVector();
           rit->second.IdentifyPerimeter();
        }

/* unique regions have already been dealt with
     for ( auto rit=UniqueRegionsBegin(); rit!=UniqueRegionsEnd(); ++rit ) {
           rit->second.CreateNodePointerVector();
           rit->second.IdentifyPerimeter(); // calls PartitionCellVector
        }
*/
       
 } // end RebuildRegions





template<size_t dim, template<size_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::RegionsOut() const
 {
     std::cout <<"\n\nRegionInterface<"<< dim <<",Region<Element>>::RegionsOut:\n";
     std::cout <<"\n\tUnique regions of model:\n";
     for ( auto rit=UniqueRegionsBegin(); rit!=UniqueRegionsEnd(); ++rit ) {
          std::cout <<"\t\t"<< (*rit).first;
          std::cout <<" "<< (*rit).second.Elements() <<" elements,";
          std::pair<int32_t, int32_t> rdim = (*rit).second.ElementSpatialDimensions();
          if ( rdim.second == 3 )
            std::cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 2 )
            std::cout <<" surface area (m2): "<< (*rit).second.Volume() <<", perimeter length (m): "<< (*rit).second.SurfaceArea();
          std::cout <<", range of spatial dimensions: "<< rdim.first <<", highest spatial dimension "<< rdim.second << std::endl;
       }
     std::cout <<"\n\tNon-unique regions of model:\n";
     for ( auto rit=RegionsBegin(); rit!=RegionsEnd(); ++rit ) {
          std::cout <<"\t\t"<< (*rit).first;
          std::cout <<" "<< (*rit).second.Elements() <<" elements,"<< (*rit).second.Volume();
          std::pair<int32_t, int32_t> rdim = (*rit).second.ElementSpatialDimensions();
          if ( rdim.second == 3 )
            std::cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 2 )
            std::cout <<" surface area (m2): "<< (*rit).second.Volume() <<", perimeter length (m): "<< (*rit).second.SurfaceArea();
          std::cout <<", range of spatial dimensions: "<< rdim.first <<", highest spatial dimension "<< rdim.second << std::endl;
       }
     std::cout << std::endl << std::endl;
 }




// EXPLICIT TEMPLATE INSTANTIATION FOR REGION_INTERFACE
template class RegionInterface<1U, Model>;
template class RegionInterface<2U, Model>;
template class RegionInterface<3U, Model>;


} // csmp
