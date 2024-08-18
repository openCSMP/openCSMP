#include "RegionInterface.h"
#include "Region.h"
#include "Model.h"
#include "PropertyConstraints.h"
#include "ModelTopology.h"
#include "MeshManagementUtilities.h"
#include "UnionFind.h"
#include "binaryReadWrite.h"
#include "ErrorHandler.h"


using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsBegin()
{ return uniqueRegionMap_.begin(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsEnd()
{ return uniqueRegionMap_.end(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::RegionsBegin()
{ return regionMap_.begin(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionIterator  RegionInterface<dim, REGION_COMPLEX>::RegionsEnd()
{ return regionMap_.end(); }


template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsBegin() const
{ return uniqueRegionMap_.begin(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator  RegionInterface<dim, REGION_COMPLEX>::UniqueRegionsEnd() const
{ return uniqueRegionMap_.end(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator RegionInterface<dim, REGION_COMPLEX>::RegionsBegin() const
{ return regionMap_.begin(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
typename RegionInterface<dim, REGION_COMPLEX>::regionConstIterator RegionInterface<dim, REGION_COMPLEX>::RegionsEnd() const
{ return regionMap_.end(); }


template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::Regions() const
{ return uniqueRegionMap_.size() + regionMap_.size(); }

template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::UniqueRegions() const
{ return uniqueRegionMap_.size(); }






/**
    Changes the subdomain name and the search key in the region map with deletion or copying of elements
    @return reports on whether the name change was successful.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool  RegionInterface<dim, REGION_COMPLEX>::RenameRegion( const string& old_name, const string& new_name )
 {
    // does the region to be renamed exist
    const auto uniqueRegionIterator = uniqueRegionMap_.find(old_name);
    const auto regionIterator       = regionMap_.find(old_name);
    if ( uniqueRegionIterator == uniqueRegionMap_.end() && regionIterator == regionMap_.end() ) return false;
    
    // renaming the corresponding model subdomain
    if ( uniqueRegionIterator != uniqueRegionMap_.end() )
      (*uniqueRegionIterator).second.Name(new_name );
    else
      (*regionIterator).second.Name(new_name );
         
    // moving the unique region according to its new key
    if ( uniqueRegionMap_.find(old_name) != uniqueRegionMap_.end() ) {
         auto regionHandler = uniqueRegionMap_.extract(old_name);
         regionHandler.key() = new_name;
         uniqueRegionMap_.insert(std::move(regionHandler));
      }

    // moving the region according to its new key
    if ( regionMap_.find(old_name) != regionMap_.end() ) {
         auto regionHandler = regionMap_.extract(old_name);
         regionHandler.key() = new_name;
         regionMap_.insert(std::move(regionHandler));
      }

    return true;
 }




/**
Method first searches the Region in the unique region map, then in the
non-unique region map list where regions may overlap. If the desired region
can not be found, an ERROR csmp::Exception is thrown

In the case of failure, a reference to the Model Region may be returned.
This region is guaranteed to be
there always.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
const Region<dim>&  RegionInterface<dim, REGION_COMPLEX>::Region( const string& region_name ) const
{
  if ( region_name.empty() )
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region: ",
                           "regions search string is empty" );

  // first a look in the unique region list
  typename map<string, csmp::Region<dim> >::const_iterator  iter( uniqueRegionMap_.find( region_name ) );
  if ( iter != uniqueRegionMap_.end() )
    return (*iter).second;

  // now a look at the generic region list
  iter = regionMap_.find( string( region_name ) );

  if ( iter != regionMap_.end() )
    return (*iter).second;
  else
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region: ",
                           (string( "region does not exist: " ) + region_name) );

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
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
Region<dim>&  RegionInterface<dim, REGION_COMPLEX>::Region( const string& region_name )
{
  if ( region_name.empty() )
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region: ",
                           "regions search string is empty" );

  typename map<string, csmp::Region<dim> >::iterator  iter( uniqueRegionMap_.find( region_name ) );
  if ( iter != uniqueRegionMap_.end() )
    return (*iter).second;

  // see whether this is perhaps a non-unique region
  iter = regionMap_.find( region_name );

  if ( iter != regionMap_.end() )
    return (*iter).second;
  else
    throw csmp::Exception( ERROR,
                           "RegionsInterface<dim,REGION_COMPLEX>::Region:",
                           (string( "region does not exist: " ) + region_name) );

  return Region( "Model" );
}



template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool  RegionInterface<dim, REGION_COMPLEX>::IsUnique( const string& region_name ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion(region_name) )
    csmp_error.Note( ERROR, "RegionInterface<dim, REGION_COMPLEX>::IsUnique:",
                       region_name, "does not exist.");

  typename map<string, csmp::Region<dim> >::const_iterator  iter( uniqueRegionMap_.find( region_name ) );
  if ( iter != uniqueRegionMap_.end() )
    return true;
  return false;
}



template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::ContainsRegion( const string& region_name ) const
{
  if ( uniqueRegionMap_.find( region_name ) != uniqueRegionMap_.end() )
    return true;
  if ( regionMap_.find( region_name ) != regionMap_.end() )
    return true;

  return false;
} // end ContainsRegion




/**
     Finds Region by its unique region ID (=domain idx) assigned to any ModelSubDomain upon creation, this method searches for the corresponding unique region
 */
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
Region<dim>&  RegionInterface<dim, REGION_COMPLEX>::RegionByDomainIndex( int32_t domain_index )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( uniqueRegionMap_.empty() )
      throw csmp::Exception( ERROR,
                             "RegionsInterface<dim,REGION_COMPLEX>::RegionByDomainIndex:",
                             "model does not contain any unique regions" );
    // linear search
    assert( domain_index > 0 );
    for ( auto& rit : uniqueRegionMap_ ) {
         if ( rit.second.DomainIndex() == domain_index )
           return rit.second;
      }

    csmp_error.Note( ERROR, "RegionInterface<dim, REGION_COMPLEX>::RegionByDomainIndex:",
                         to_string( domain_index ), "region with this domain index does not exist; returning 'Model'.");
    
    return Region("Model");
 }



// TODO: remove contiguity requirement in the presence of SplitBoundaries
  /// checks that there is a model region and that it contains elements
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::HasValidModelRegion() const
 {
    // does the region model exist?
    if ( !this->ContainsRegion("Model") ) return false;
    
    // is it contiguous?
    //if ( !this->Region("Model").IsContiguous() ) return false;
    
    return true;
 }




/**

Takes all elements from the mesh manager and creates a region that contains them.

Uses the highest order of elements (dim == model dimension) to perform a flood fill to determine whether the region is contiguous.

@return the number of elements contained in the region "Model" and whether this region is contiguous (true) or discontiguous (false).

@author SKM
@date 4/9/21

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormModelRegion( bool is_unique )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
  // 0. initial checks
  // -----------------
  const string regionname("Model");
  if ( ContainsRegion( "Model" ) ) {
      csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                         regionname, "'Model' region already exists; nothing was done." );
      return false;
    }

  // 1. building the 'Model' region
  // ------------------------------
  pair<typename map<string,csmp::Region<dim> >::iterator, bool>
    newRegion = (is_unique) ?
    uniqueRegionMap_.insert( make_pair( regionname, csmp::Region<dim>( regionname,
                            static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
    :
    regionMap_.insert( make_pair( regionname, csmp::Region<dim>( regionname,
                      static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  // if region was inserted successfully
  if ( newRegion.second ) {
       size_t elmts = (*newRegion.first).second.AccumulateAll( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh() );
       if ( elmts == 0 ) //                     ^^^^^^^^^^^^^
         csmp_error.Note( ERROR, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                            regionname, "region could not be formed." );
       // do not renumber
       // (*newRegion.first).second.UpdateMemberIndexes();
    }
  else {
      csmp_error.Note( ERROR, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                         regionname, "region could not be formed." );
      return 0U;
    }

  // checking that all elements and nodes were discovered
  if ( (*newRegion.first).second.Cells() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Elements() or
       (*newRegion.first).second.Nodes() != static_cast<REGION_COMPLEX<dim>*>(this)->Mesh().Nodes() )
    {
       csmp_error.Note( ERROR, "RegionInterface<dim,REGION_COMPLEX>::FormModelRegion:",
                          regionname, "'Model' not all elements were incorporated into the new 'Model' region." );
    }
    
  // setting the region ID
  if ( is_unique )
    (*newRegion.first).second.SetRegion_ID();

  return (*newRegion.first).second.Cells();
  
} // end FormModelRegion




/**
      Forms unique regions called MATERIAL1..n from the material IDs assigned to the elements.
      
      @attention this method assumes that the unique numbers of elements, faces, and interfaces via the MeshManager
      
      TODO: add a PropertyConstraint here
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromMaterialIDs()
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
         csmp_error.Note( ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromMaterialIDs:",
                           "Material IDs do not appear to have been initialized; nothing was done.");
         return 0u;
      }

    // creating the regions from the material identifiers
    for ( const auto& mit : mtrl_ids )
      {
         region_name = "MATERIAL" + to_string( mit.first );
         pair<typename map<string, csmp::Region<dim> >::iterator, bool>
           it = uniqueRegionMap_.insert( make_pair( region_name, csmp::Region<dim>( region_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
         // if the region was successfully inserted
         if ( it.second )
           {
             // retrieving the elements by their IDs and assigning them  to the region
             (*it.first).second.Accumulate( mit.second.begin(), mit.second.end() );

             // removing the region if it contains no elements
             if ( (*it.first).second.Cells() == 0U ) {
                 uniqueRegionMap_.erase( it.first );
                 csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromMaterialIDs",
                                    region_name, "could not be formed" );
               }
             else n_regions++;
          }
        else
          throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionsFromMaterialIDs",
                                 region_name, "could not be formed. Does this region already exist?" );

        // setting region IDs
        (*it.first).second.SetRegion_ID();

      } // end materials loop

    return n_regions;
    
 } // end FormRegionsFromMaterialIDs




/**

Removes named Region object and (optionally) its elements and nodes.

@param  regionName The name of the region object which shall be removed.
@param  erase_elmts_and_update_connectivity  gets MeshManager to delete elements and nodes and rebuilt local connectivity.

@note If the region which shall be removed does not exist, the method reports a warning.

@attention to delete underlying elements, the MeshManager had to be called upon.

@attention this call in itself is not enough to remove the elements defining a regiion, but its perimeter elements need to be disconnected from their
neighbors and deleted. These tasks are done by the MeshManager:

Mesh().DetachOutsideNeighborsAlongPerimeter(  ModelSubdomain&  );

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::RemoveRegion( const char* regionName,
                                                         bool erase_elmts_and_update_connectivity )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 0. check whether region exists (should be a notice only, nothrow)
  if ( !ContainsRegion( regionName ) ) {
       csmp_error.Note( WARNING, "RegionsInterface<dim,Model>::RemoveRegion",
                         "region did not exist: ", regionName );
       return;
    }
    
  if ( erase_elmts_and_update_connectivity && !IsUnique(regionName) ) {
       csmp_error.Note( ERROR, "RegionsInterface<dim,Model>::RemoveRegion",
                          regionName, "is potentially overlapping other regions; this case is not handled yet" );
       return;
    }

  REGION_COMPLEX<dim>&  regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );
  csmp::Region<dim>&    region = this->Region( regionName );

  if ( erase_elmts_and_update_connectivity )
    // get mesh manager to delete elements and nodes and fix up the connectivity
    regionComplex.Mesh().DeleteElementsAndRepairConnnectivity( region.CellVector().begin(), region.CellVector().end() );

  // finding the region in the corresponding map
  typename map<string, csmp::Region<dim> >::iterator
    iterRegion( regionMap_.find( string( regionName ) ) ),
    iterUniqueRegion( uniqueRegionMap_.find( string( regionName ) ) );

  // if the region was found in the respective map, it is erased
  if ( iterRegion != regionMap_.end() ) {
       regionMap_.erase( string( regionName ) );
    }
  if ( iterUniqueRegion != uniqueRegionMap_.end() ) {
       uniqueRegionMap_.erase( string( regionName ) );
    }

} // end RemoveRegion



// DEBUG - check element vector for duplicates (OK)
//set<const csmp::Element<dim>*> elmts( subdomain.CellsBegin(), subdomain.CellsEnd() );
//assert( elmts.size() == subdomain.Cells() );
// DEBUG - check element vector for dead elements (OK - no corrupt elements)
//cerr <<"\nRegionInterface<dim, REGION_COMPLEX>::RemoveRegion: "<< regionName <<"\n";
//for ( auto it=subdomain.CellsBegin(); it!=subdomain.CellsEnd(); ++it )
//  cerr <<" "<< parseFiniteElementType( (*it)->FE_Type() ) <<":"<< (*it)->Idx();
//cerr << endl;
// DEBUG - check MeshManager vector for dead elements (OK - no corrupt elements)
//cerr <<"\nRegionInterface<dim, REGION_COMPLEX>::RemoveRegion: "<< regionName <<"\n";
//for ( auto it=meshMgr.CellsBegin(); it!=meshMgr.CellsEnd(); ++it )
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
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::OutputRegionsToBinary( const char* file_name ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( regionComplex.Database() );

  string bin_file( file_name );
  fstream fp( bin_file.c_str(), ios::out | ios::binary );
  if ( !fp.is_open() ) {
      csmp_error.Note( ERROR, "RegionInterface<dim,REGION_COMPLEX>::OutputRegionsToBinary:",
                         bin_file, "file could not be opened; nothing was done." );
      return;
    }

  // -----------------------------
  // 1. File header
  // -----------------------------
  {
    BinaryFileSectionWrite hdr( fp, "REGFHEDR" );

    string heading( "RegionInterface::OutputRegionsToBinary: " );
    heading += "region information for Model '";
    heading += regionComplex.Name();
    heading += "' written to file: ";
    heading += bin_file;
    heading += "'.";

    binaryFileWrite( fp, heading.c_str() );
  }

  cout << "\nRegionInterface<dim,REGION_COMPLEX>::OutputRegionsToBinary: regions written to binary file: ";
  cout << "\n\n\tUnique regions: ";

  // -----------------------------
  // 2. writing the unique regions
  // -----------------------------
  {
    BinaryFileSectionWrite uhdr( fp, "UNIQREGN" );

    int64_t  records = static_cast<int64_t>(this->UniqueRegions());
    
    // writing number of unique regions
    fp.write( reinterpret_cast<const char*>(&records), sizeof( int64_t  ) );

    for ( typename map<string, csmp::Region<dim> >::const_iterator
          git = UniqueRegionsBegin(); git != UniqueRegionsEnd(); ++git )
      {
        BinaryFileSectionWrite hdr( fp, "ONE_REGN" );
        cout <<"'"<< (*git).first <<"' ";
        cout.flush();
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
    BinaryFileSectionWrite nhdr( fp, "NONUREGN" );

    int64_t  records = static_cast<int64_t>(this->Regions() - this->UniqueRegions());

    fp.write( reinterpret_cast<const char*>(&records), sizeof( int64_t  ) );

    for ( auto git = RegionsBegin(); git != RegionsEnd(); git++ )
      {
        BinaryFileSectionWrite hdr( fp, "ONE_REGN" );
        cout << (*git).first << " ";
        (*git).second.WriteDomainIndexesToBinaryFile( fp );
        domainVariablesOut( fp, (*git).second, database );
      }
    cout << endl;
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
  cout << "\nRegionInterface<" << dim << ",REGION_COMPLEX>::OutputRegionsToBinary: file '";
  cout << bin_file << "' has been successfully written.\n";

} // end OutputRegionsToBinary






/**
Reads all regions stored by OutputRegionsToBinary() into Interface,
constructing them using the SubDomainInfo data, thereby avoiding costly
re-initialisation.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::InputRegionsFromBinary( const char* file_name,
                                                                   const set<string>& subset_variables )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );
  MeshManager<dim>& mesh( static_cast<REGION_COMPLEX<dim>& >(*this).Mesh() );

  string bin_file( file_name );
  fstream fp( bin_file.c_str(), ios::in | ios::binary );
  if ( !fp.is_open() ) {
    csmp_error.Note( ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
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
      cout << "\nRegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary: Reading file header:\n\t" << text << endl;
    }
  cout << "\n\timporting the regions: ";

  const PropertyDatabase<dim>& database( static_cast<const REGION_COMPLEX<dim>& >(*this).Database() );

  // -----------------------------
  // 2. unique regions
  // -----------------------------
  cout << "\n\tunique regions: ";
    {
      BinaryFileSectionRead uhdr( fp, "UNIQREGN" );
      // getting number of unique region records from file
      int64_t   records( 0 );  // region records
      fp.read( reinterpret_cast<char*>(&records), sizeof( int64_t  ) );
      if ( records > 0 )
          // reading the regions sequentially
          for ( int64_t i{0UL}; i<records; i++ )
            {
              BinaryFileSectionRead hdr( fp, "ONE_REGN" );
              // reading name and element indices for each unique region
              SubDomainInfo  info;
              readDomainIndexesFromBinaryFile( dim, fp, info );
              // reconstruct the region
              pair<typename map<string, csmp::Region<dim> >::iterator, bool>
                it = uniqueRegionMap_.insert( make_pair( info.name, csmp::Region<dim>( database, mesh, info ) ) );

              if ( !it.second )
                throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                                       info.name, "Region could not be formed; issue with binary file." );
              else
                // new domain indices matching the ref-count of regions must be assigned
                (*it.first).second.SetRegion_ID();

              cout << "\n\t\t'" << (*it.first).first << "'(" << (*it.first).second.Cells() << " elements).";
              
              // reading variables stored on the regions
              if ( subset_variables.empty() ) domainVariablesIn( fp, (*it.first).second, database );
              else selectedDomainVariablesIn( fp, (*it.first).second, database, subset_variables );
            }
       else csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                               bin_file, "does not contain any unique region descriptions; no regions were initialised." );
    }

  // -----------------------------
  // 3. reading non-unique regions
  // -----------------------------
  cout << "\n\tnon-unique regions: ";
  {
    BinaryFileSectionRead nhdr( fp, "NONUREGN" );
    // getting number of non-unique region records from file
    int64_t   records( 0 );
    fp.read( reinterpret_cast<char*>(&records), sizeof( int64_t  ) );
    if ( records > 0 )
      // reading the regions sequentially
      for ( auto i{0U}; i<records; i++ )
        {
          BinaryFileSectionRead hdr( fp, "ONE_REGN" );
          // reading name and element indices for each unique region
          SubDomainInfo  info;
          readDomainIndexesFromBinaryFile( dim, fp, info );
          // if the region info record is not empty the region is reconstructed
          pair<typename map<string, csmp::Region<dim> >::iterator, bool>
            it = regionMap_.insert( make_pair( info.name, csmp::Region<dim>( database, mesh, info ) ) );
          if ( !info.interior_elmts.empty() ) {
              if ( !it.second )
                throw csmp::Exception( FATAL_ERROR, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
                                       info.name, "Region could not be formed; issue with binary file." );

              cout << "\n\t\t'" << (*it.first).first << "'(" << (*it.first).second.Cells() << " elements).";

              // reading variables stored on the regions
              if ( subset_variables.empty() ) domainVariablesIn( fp, (*it.first).second, database );
              else selectedDomainVariablesIn( fp, (*it.first).second, database, subset_variables );
            }
        }
    else csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::InputRegionsFromBinary:",
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
  cout << "\n\nRegionInterface<" << dim << ",REGION_COMPLEX>::InputRegionsFromBinary: file '";
  cout << bin_file << "' has been read successfully.\n";

} // end InputRegionsFromBinary











/**
RegionsFromPropertyValues() defines one region for each value of the target
property.  The user will be prompted to assign a name to each region as these
are being created. The names of successfully created regions are returned
into the list argument.

@param prop the property on the basis of whose variations
the region regions will be defined.

The names of the newly created regions are returned into an STL set
which uses the less<> functional to order the names alphabetically. If
the set is not empty, it will be erased before the region names are
stored within it.

@section implementation Implementation

Thus far, the method only handles scalar variables. Since there is always
just one property value per element, the regions will be unique.

@section application Application

RegionsFromPropertyValues() is useful if the desired subregions of a model
coincide with changes of a specific property. Clearly, one does not want
to apply this method on a continuously changing property, since it might
produce as many Regions as there are finite-elements in the mesh.

@section messages Messages

The method is interactive and will prompt the user for the names of the
regions which are created in the course of its execution.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues( const char* prop,
                                                                            set<string>& region_names )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion("Model") )
    csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
  if ( !region_names.empty() )
    region_names.erase( region_names.begin(), region_names.end() );

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
        for ( auto i{0U}; i < eit->IntegrationPoints(); i++ )
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

  // 2. Prompting user for region names and making regions
  // -----------------------------------------------------
  cout << "\nModel<" << dim << ">::FormRegionsFromPropertyValues: Generating region names";
  cout << " corresponding to unique values of the property: '" << prop << "' ";
  cout << "[" << regionComplex->Database().Unit( prop ) << "]." << endl;
  string  gname( "region_" );
  gname += prop;
  size_t  group_idx( 0U );
  char    num[30U];

  for ( auto it = groups.begin(); it != groups.end(); it++ )
  {
    cout << "\nProperty value: " << (*it).first;
    snprintf( num, sizeof(num), "%lu", group_idx++ );
    (*it).second = string( gname + num );

    if ( (*it).second != "undefined" )
      {
        assert( !ContainsRegion( (*it).second.c_str() ) );
        FormRegionFrom( (*it).second.c_str(), prop,
                        (*it).first - numeric_limits<double>::epsilon(),
                        (*it).first + numeric_limits<double>::epsilon(), true );

        region_names.insert( (*it).second );
      }
    else
      cout << "\nModel<" << dim << ">::FormRegionsFromPropertyValues: check predefined region name." << endl;
  }
  
  return groups.size();

} // end RegionsFromPropertyValues






/**
     Uses the region model, to form a region from its elements with the corresponding ID numbers.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionName, vector<size_t>& elmt_ids, bool unique/*=false */ )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion("Model") )
    csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  csmp::Region<dim>& model_domain(Region("Model"));
  
  pair<typename map<string, csmp::Region<dim> >::iterator, bool> it;
  if ( unique ) {
        it = this->uniqueRegionMap_.insert( make_pair( regionName, csmp::Region<dim>( regionName, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
        if ( it.second )
          (*it.first).second.AccumulateByNumber( model_domain.CellsBegin(), model_domain.CellsEnd(), elmt_ids );

        // setting the region ID
        (*it.first).second.SetRegion_ID();
    }
  else {
        it = this->regionMap_.insert( make_pair( regionName, csmp::Region<dim>( regionName, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
        if ( it.second )
          (*it.first).second.AccumulateByNumber( model_domain.CellsBegin(), model_domain.CellsEnd(), elmt_ids );
    }

  return (*it.first).second.Cells();
}




/**
  Forms a region on the basis of its node coordinates.
  
  @param regionName unique name for the region that shall be formed
  @param cnr_min minimum x,y,z of bounding box
  @param cnr_max maximum coordinate of bounding box
  @param all_nodes_must_be_within to choose whether just a single node of the targeted elements must be within the bbox or all of them
  @return the number of elements found inside of the bounding box added to the target region
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionInBoundingBox( const char* regionName,
                                                                      const Point<dim>& cnr_min, const Point<dim>& cnr_max,
                                                                      bool all_nodes_must_be_within )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion("Model") )
    csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionInBoundingBox:",
                      "method relies on the existence of region 'Model', which does not exist");

  csmp::Region<dim>&     model_domain(Region("Model"));
  vector<Element<dim>*>  inside_elmts;
  
  if ( all_nodes_must_be_within ) {
       for ( const auto& it : model_domain.CellVector() ) {
            bool elmt_is_inside{ true };
            for ( auto nit=it->NodesBegin(); nit!=it->NodesEnd(); ++nit )
              if ( !isWithinBoundingBox<dim>( cnr_min, cnr_max, (*nit) ) ) {
                   elmt_is_inside = false;
                   break;
                }
            if ( elmt_is_inside ) inside_elmts.push_back( it );
         }
    }
  // if only a single element node is required to lie withing the bounding box
  // we add all parent elements for nodes that are contained in it
  else {
       for ( const auto& nit : model_domain.NodeVector() )
         if ( !isWithinBoundingBox<dim>( cnr_min, cnr_max, nit ) )
           for ( uint32_t i{0u}; i<nit->Parents(); ++ i )
             if ( nit->Parent(i) )
               inside_elmts.push_back( nit->Parent(i) );
       // since this initialisation of the vector leads to duplicates it has to be sorted and made unique
       sort( inside_elmts.begin(), inside_elmts.end() );
       inside_elmts.erase( unique( inside_elmts.begin(), inside_elmts.end() ), inside_elmts.end() );
    }

  // this region will always overlap other ones
  const bool unique{ false };
  return FormRegionFrom( regionName, inside_elmts.begin(), inside_elmts.end(), unique );
}







/** Assuming that the supplied n-nodes are in sequence of a polyline, method forms n-1 line elements putting them into a now limer-dimensional region
 */
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionName, vector<Node<dim>*>& nodes, bool is_unique )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( !ContainsRegion("Model") )
      csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                        "method relies on the existence of region 'Model', which does not exist");

    REGION_COMPLEX<dim>* model( static_cast<REGION_COMPLEX<dim>*>(this) );
    MeshManager<dim>&    mesh = model->Mesh();

    // getting the MeshManager to create the required number of line elements
    assert( mesh.FiniteElements().InterpolationOrder() == 1 );
    const CSMP_FEM_TYPE fe_type = (mesh.FiniteElements().UsesElementsWithLocalCoordinateSystem()==true) ? ISOPARAMETRIC_LINEAR_BAR : LINEAR_BAR;
    const LocalVariables            lvsElementVars( model->Database().LocalVariablesAt(ELEMENT) );
    const IntegrationPointVariables lvsIntegrationPointVars( model->Database().IntegrationPointVariablesAt(ELEMENT) );
    int32_t                         material_ID=numeric_limits<int32_t>::max(); // TODO: should be the domain index
    
    // forming element vector (the elements get their material from the domain index
    vector<Element<dim>*> elmts;
    elmts.reserve( nodes.size() - 1u );
    for ( auto nit1=nodes.begin(), nit2=next(nodes.begin(),1); nit2!=nodes.end(); ++nit1, ++nit2 )
      elmts.push_back( mesh.AddElement( fe_type, lvsElementVars, lvsIntegrationPointVars,
                                        vector<Node<dim>*>{ (*nit1), (*nit2) }, material_ID ) );
                                        
    // connectiong the elements with one another
    mesh.template BuildLineConnectivity<Element>( elmts.begin(), elmts.end() );
 
    // forming the region from the new line element vector
    return FormRegionFrom( regionName, elmts.begin(), elmts.end(), is_unique );

 } // end FormRegionFrom(nodes)







  /// forms unique or non-unique region from range of elements; returns reference to it
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
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
    pair<typename map<string, csmp::Region<dim> >::iterator, bool>  it;
    if ( unique )
      it = uniqueRegionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    else
      it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

    if ( it.second ) {
        (*it.first).second.Accumulate( first, last );

        // removing the region if it contains no elements
        if ( (*it.first).second.Cells() == 0U ) {
          regionMap_.erase( it.first );
          csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                             "Region could not be formed", output_region.c_str() );
          return 0U;
        }
      }
    else {
        csmp_error.Note( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                           "Region could not be formed. ", output_region.c_str() );
        return 0U;
      }

    if ( unique ) (*it.first).second.SetRegion_ID();

    return (*it.first).second.Cells();
    
 } // end FormRegionFrom (range of element pointers)



/**
FormAndAddDomain() used combined constraints suppplied in the form of
a constraints object to determine which elements shall be used to
form a (unique/non-unique) region with the target name.

@param region_name the name of the region that shall be formed and a reference to
the initialized PropertyConstraints object which must contain the
ranges of the variables that shall be used to discriminate elements
that shall be accumulated into the region.

@section application Application

To form regions from dynamic criteria like the combined pressure and
temperature ranges that get computed via Algorithms or other criteria
that are not known at the onset of a computation.

@section messages Messages

The method reports an INFO if no values fall into the target ranges and
a FATAL_ERROR if the region cannot be created because the name is already
in use.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* region_name,
                                                             PropertyConstraints& constraints,
                                                             bool unique_region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion("Model") )
    csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  string output_region( region_name );
  if ( ContainsRegion( region_name ) ) {
    cout << "\nModel<" << dim << ">::FormRegionFrom: WARNING: region '" << region_name;
    cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    cout << output_region << endl;
  }

  // the 'bool' member of pair indicates whether insertion into map worked or not
  pair<typename map<string, csmp::Region<dim> >::iterator, bool>  it;
  if ( unique_region )
    it = uniqueRegionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
  else
    it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  if ( it.second )
    {
      REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>*>(this) );
      constraints.InitializePropertyIndices( regionComplex->Database() );
      regionTraits_.insert( make_pair( region_name, constraints ) );
      (*it.first).second.AccumulateWithinRange( regionComplex->Mesh(), constraints );

      // removing the region if it contains no elements
      if ( (*it.first).second.Cells() == 0U ) {
        regionMap_.erase( it.first );
        csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                           "Region could not be formed", output_region.c_str() );
        return 0U;
      }
    }
  else {
    csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                       "Region could not be formed", output_region.c_str() );
    return 0U;
  }

  if ( unique_region ) (*it.first).second.SetRegion_ID();

  return (*it.first).second.Cells();

} // end FormRegionFrom







/**
Forms a region of finite elements who host property values which lie within
a user-specified range. The new region is added to the Model region map.

The selection process implies that this type of region can include several
regions inside your problem domain. For instance, all shale horizons as
identified by a clay content between 50 to 100% may constitute a new region.

When a new region is formed a region internal flag will be assigned to
each node, constraint point and element.
There are two region-internal object flags, PLAIN and BOUNDARY. When a
new region is formed, the Region method IdentifyBoundaryAs() assigns the
region-internal object flags, depending on whether the nodes, constraint
points or elements in the region lie at the region boundary or inside of
the region. Elements are assigned a boundary flag if at least one of
their faces coincides with the region boundary.

@section arguments Input Arguments

The newly formed region of finite elements will have a name specified
by the first method argument. The selection criterion is that values of
the physical variable identified by the second argument, are within the
open interval given by [min,max].

@section implementation Implementation

Once all member elements of the new region have been identified, the region
analyzes which elements, nodes, and constraint points lie at its boundary.
Subsequently the region is added to a map of regions that is contained in
the Model object.

Inn this and other region-forming methods Regions are first added to the
region map before they are filled, avoiding a costly copy construction
of non-empty regions.

@section application Application

To form regions on the basis of characteristic material properties.

@section messages Messages

A warning is issued if a region with the same name already exists or if no
elements with the desired properties were found.

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* region_name,
                                                             const char* prop,
                                                             double min, double max,
                                                             bool unique_region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  string output_region( region_name );
  if ( ContainsRegion( region_name ) ) {
    cout << "\nModel<" << dim << ">::FormRegionFrom: WARNING: region '" << region_name;
    cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    cout << output_region << endl;
  }

  // the 'bool' member of pair indicates whether insertion into map worked or not
  pair<typename map<string, csmp::Region<dim> >::iterator, bool>  it;
  if ( unique_region )
    it = uniqueRegionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
  else
    it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );

  if ( it.second )
    {
      (*it.first).second.AccumulateWithinRange( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh(), prop, min, max );

      // removing the region if it contains no elements
      if ( (*it.first).second.Cells() == 0U ) {
        regionMap_.erase( it.first );
        csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                           "Region could not be formed", output_region.c_str() );
        return 0U;
      }
    }
  else {
    csmp_error.Note( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRegionFrom",
                       "Region could not be formed. ", output_region.c_str() );
    return 0U;
  }

  if ( unique_region ) (*it.first).second.SetRegion_ID();

  return (*it.first).second.Cells();

} // end FormRegionFrom











/**
Forms a Region from a set of region names by using the method MergeRegions.

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* regionname, const set<string>& region_names )
{
    string output_region( regionname );
    if ( ContainsRegion( regionname ) ) {
      cout << "\nModel<" << dim << ">::FormRegionFrom (element numbers): WARNING: region '" << regionname;
      cout << "' already exists, adding an underscore at end of name: ";
      output_region += "_";
      cout << output_region << endl;
    }

    this->MergeRegions( region_names, output_region.c_str() );
    
    return this->Region(output_region.c_str()).Cells();
}



/**
@author P Lang

This forms a region from elements eligible as reported by elementComp.

ElementComp is a model of binary predicates, i.e.
@code
template<uint32_t dim>
struct ElementsLessX
{
bool operator () ( Element<dim> const* ePtr ) const
{
return ...is ePtr eligible?;
}
};
@endcode
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
template<template<uint32_t> class ElementComp>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionFrom( const char* newRegionName, ElementComp<dim> const& elementComp, const char* hostRegion )
  {
    static_cast<REGION_COMPLEX<dim>*>(this)->UpdateIndices();
    csmp::Region<dim>& rref( this->Region( hostRegion ) );
    rref.UpdateMemberIndexes();
    vector<uint32_t> elementIds;
    elementIds.reserve( rref.Elements() );
    for ( typename csmp::Region<dim>::SimplexContainer::const_iterator it( rref.CellsBegin() ); it != rref.CellsEnd(); ++it )
      if ( elementComp( (*it) ) )
        elementIds.push_back( (*it)->Idx() );
    vector<uint32_t>( elementIds ).swap( elementIds );
    if ( elementIds.empty() )
      return 0U;
    
    return this->FormRegionFrom( newRegionName, elementIds, false );
  }





/**
Forms regions using the element ID containers stored in the model topology
object. The regions are numbered in their alphabetical order and corresponding numbers
are assigned to the material ID of the element class.

@param ignore_domain_type_identifiers when true all topological entities are turned into regions
even if the contain boundary names etc.

@note the material IDs may later be overwritten by rocktypes.

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
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRegionsFrom( const ModelTopology& topo, bool ignore_domain_type_identifiers )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. getting the names of the regions
  list<string> regions;
  if ( ignore_domain_type_identifiers ) topo.OutputAll( regions );
  else topo.OutputRegions( regions );

  // 2. assigning the elements to regions in the Model
  cout << "\nRegionInterface::FormRegionsFrom: Forming the regions: ";

  uint32_t new_regions( 0U );
  for ( typename list<string>::const_iterator lit = regions.begin(); lit != regions.end(); lit++ )
    {
      string region_name( *lit );
      pair<typename map<string, csmp::Region<dim> >::iterator, bool>
        it = uniqueRegionMap_.insert( make_pair( region_name, csmp::Region<dim>( region_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
      // if the region was successfully inserted
      if ( it.second )
        {
          // making a list of the element numbers
          vector<size_t>  element_ids;
          element_ids.reserve( topo.CellsWithinDomain( (*lit).c_str() ) );
          copy( topo.CellsOfDomainBegin( (*lit).c_str() ),
                topo.CellsOfDomainEnd( (*lit).c_str() ),
                back_inserter( element_ids ) );

          // retrieving the elements by their IDs and assigning them  to the region
          (*it.first).second.AccumulateByNumber( static_cast<REGION_COMPLEX<dim>*>(this)->Mesh(), element_ids );
          element_ids.erase( element_ids.begin(), element_ids.end() );

          // removing the region if it contains no elements
          if ( (*it.first).second.Cells() == 0U ) {
              uniqueRegionMap_.erase( it.first );
              csmp_error.Note( WARNING, "RegionsInterface::FormRegionsFrom",
                                 "Region could not be formed", (*lit).c_str() );
            }
          else {
               // assigning unique material IDs and region identifiers to the element members of the region
               for ( auto eit=(*it.first).second.CellsBegin(); eit!=(*it.first).second.CellsEnd(); ++eit ) {
                    (*eit)->Material_ID( new_regions );
                    (*eit)->Region_ID( (*it.first).second.DomainIndex() );
                 }
               // reporting the name of the newly generated region
               cout << region_name << " ";
               
               (*it.first).second.SetRegion_ID();

               new_regions++;
            }
        }
      else
        throw csmp::Exception( ERROR, "RegionsInterface::FormRegionsFrom",
                               "Region could not be formed. Does this region already exist?", (*lit).c_str() );
    }
  cout << endl;

  return new_regions;

} // end FormRegionsFrom










  /**
      Breaks non-contiguous Regions into contiguous subregions that carry the name
      of the region but have a number as suffix to their name to denote the
      partition.
      
      @return The method returns the number of subgregions that were created.
      
      @attention The master region that was successfully partitioned is removed.

      @param region The name of the region that may be non-contiguous.
      If so, new sbregions will be created to the name of which integers
      will be appended that correspond to the number of subdomains
      that are created in this process.

      @return The method returns the number of contiguous subdomains which it
      created.

      @section implementation Implementation

      The method uses the element neighbor information to perform floodfills
      as required until all member elements are partitioned.

      @section application Application

      To automatically partition regions that consist of a multitude of
      non-contiguous model subdomains so that the latter can be addressed
      individually in computations.

      @section messages Messages

      The method will report if the region is already contiguous in which
      case no changes are made.
      
      @attention this method cannot be applied to the region model or the master region
  */
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t  RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions( const char* region )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      if ( string("Model") == region ) {
         csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions:",
                           "this operation is not allowed for region 'Model' or the master region." );
            return 0U;
        }
     
      // renumbering elements and nodes of model
      csmp::Region<dim>&  gref(Region(region));
      const bool          unique_region(IsUnique(region));

      // getting a set of the element numbers of the target region
      set<Element<dim>*>  elements( gref.CellsBegin(), gref.CellsEnd() );

      // detecting via a flood-fill whether the region can be partitioned, else nothing is done
      set<Element<dim>*>  elements_contiguous_subset;
      floodFill( gref.E(0), elements_contiguous_subset );
      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( elements.size() <= elements_contiguous_subset.size() ) {
           cout <<"\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions: ";
           cout <<"region '"<< region <<"' is already contiguous, nothing was done."<< endl;
           return 0U;
        }

      // else partitions can be created
      string  region_name(region);
      string  subregion_name;
      char    num[128];
      size_t  n_subgroups(1);
    
      // creating new contiguous region from the element subset
      while ( !elements.empty() )
        {
           // creating name of contiguous subregion
           snprintf( num, sizeof(num), "%lu", n_subgroups );
           subregion_name = region_name + num;
           if ( n_subgroups == 1U ) {
                 cout <<"\nModel<"<< dim <<">::PartitionRegionIntoContiguousSubRegions: ";
                 cout <<"region '"<< region <<"' is divided into the subregion(s):\n";
             }
           cout <<"\t\t\t'"<< subregion_name <<"'";
           cout <<" ("<< elements_contiguous_subset.size() <<" elmts)"<< endl;
         
           // creating either a unique or non-unique region depending on uniqueness of original region
           pair<typename map<string,csmp::Region<dim> >::iterator,bool>
             it = ( unique_region ) ? uniqueRegionMap_.insert( make_pair( subregion_name, csmp::Region<dim>( subregion_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database()) ) )
                                   : regionMap_.insert( make_pair( subregion_name, csmp::Region<dim>( subregion_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database()) ) );
           if ( !it.second )
             throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions",
                                    subregion_name.c_str(), "region could not be formed (name is probably not unique)" );
           else {
                // accumulating the subregion
                (*it.first).second.Accumulate( elements_contiguous_subset.begin(),
                                               elements_contiguous_subset.end() );
                                               
                // copy all values of region properties from parent to child region
                (*it.first).second.LVS( gref.LVS() );
                
                // assign a region ID to the elements of the new region
                (*it.first).second.SetRegion_ID();
             }
         
           // subtracting the elements that constitute the new region from the remaining element list
           for ( typename set<Element<dim>*>::const_iterator
                 sit=elements_contiguous_subset.begin(); sit!=elements_contiguous_subset.end(); ++sit )
             elements.erase( (*sit) );

           // computing the next subset
           if ( elements.empty() ) break;
           else floodFill( (*elements.begin()), elements_contiguous_subset );
      
           n_subgroups++;
        }

      // if the region has been partitioned successfully and its name is not model, it will be removed
      if ( IsUnique( region ) ) {
           const bool also_remove_elmts{ false };
           RemoveRegion( region, also_remove_elmts );
        }

      return n_subgroups;
    
   } // end partitionRegionIntoContiguousSubRegions






/**

ANDREW BROMAGE VERSION

Breaks non-contiguous Region into contiguous subregions that carry the name
of the original region, but have a number as a suffix to their name to distinguish each partition.
All elements and the original region are retained, but the original region is moved into the 
non-unique regions storage.

@return The method returns the number of subgregions that were created, if any.

@param region The name of the region that may be non-contiguous.

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
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t  RegionInterface<dim, REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions_Bromage( const char* region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  const string  region_name( region );

  if ( region_name == "Model" ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions_Bromage:",
                       "this operation is not allowed for region 'Model' or the master region." );
    return 0U;
  }

  // renumbering elements and nodes of model
  csmp::Region<dim>&  gref( Region( region ) );
  const bool          unique_region( IsUnique( region ) );

  // Perform union-find
  UnionFind<Element<dim>*> unionFind;
  auto eend = gref.CellsEnd();
  for ( auto eit = gref.CellsBegin(); eit != eend; eit++ ) {
    auto e = *eit;
    const size_t  neighbors( e->Neighbors() );
    for ( auto i{0U}; i<neighbors; i++ )
      if ( e->Neighbor( i ) != NULL )
        unionFind.SameComponent( e, e->Neighbor( i ) );
  }

  vector<pair<size_t, Element<dim>*> > components;
  unionFind.Components( components );
  if ( components.size() <= 1 ) {
    cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions_Bromage: ";
    cout << "region '" << region << "' is already contiguous, nothing was done." << endl;
    return 0U;
  }

  // Sort all elements by component
  vector<pair<Element<dim>*, Element<dim>*>> componentMemberships;
  componentMemberships.reserve( gref.Cells() );
  for ( auto eit = gref.CellsBegin(); eit != eend; eit++ ) {
    componentMemberships.emplace_back( unionFind.resolve( *eit ), *eit );
  }
  sort( componentMemberships.begin(), componentMemberships.end() );

  // Turn components into regions
  cout << "\nModel<" << dim << ">::PartitionRegionIntoContiguousSubRegions_Bromage: ";
  cout << "region '" << region << "' is divided into the subregion(s):\n";

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
    string  subregion_name;
    char    num[128];
    snprintf( num, sizeof(num), "%lu", subgroupNum );
    subregion_name = region_name + num;
    cout << "\t\t\t'" << subregion_name << "'";
    cout << " (" << subgroupSize << " elmts)" << endl;

    pair<typename map<string, csmp::Region<dim> >::iterator, bool>
      it = (unique_region) ? uniqueRegionMap_.insert( make_pair( subregion_name, csmp::Region<dim>( subregion_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) )
      : regionMap_.insert( make_pair( subregion_name, csmp::Region<dim>( subregion_name, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( !it.second )
      throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::PartitionRegionIntoContiguousSubRegions_Bromage",
                             subregion_name.c_str(), "region could not be formed (name is probably not unique)" );
    else {
      // accumulating the subregion
      vector<Element<dim>*> subgroup;
      subgroup.reserve( subgroupSize );
      for ( auto rit = cmCBegin; rit != cmCEnd; ++rit )
        subgroup.push_back( rit->second );

        (*it.first).second.Accumulate( subgroup.begin(), subgroup.end() );

      // copy all values of region properties from parent to child region
      (*it.first).second.LVS( gref.LVS() );
    }
    cmCBegin = cmCEnd;
  }

  // if the region has been partitioned succesfully and its name is not model, it will be removed
  if ( IsUnique( region ) )
    MoveToNonUniqueRegions( region );

  return subgroupNum;

} // end partitionRegionIntoContiguousSubRegions













/**
Removes subregions of the region identified by name. The subregions are
defined as regions that have the same name as the aforementioned region,
but with numbers appended, e.g., 'fractures' and 'fractures1'.

If the user supplies the string 'all subregions' as method argument
all subregions in the current Model object will be removed.

@param region The name of the master region, the subregions of which shall be removed.

@section application Application

To save memory by removing excessive region subdivisions as can be created
via the method PartitionRegionIntoContiguousSubRegions().

@section messages Messages

The method will report the names of the subregions that were removed.

@return the method will return the number of subregions that it removed.

@attention the master region is not touched or recreated; to achieve this use MergeRegions() before

@todo (1) SKM: logic of this method seems to be broken and it does not always work. Refactor!
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::RemoveRegionPartitionsFor( const char* region )
{
  string            target( region );
  set<string>  region_names, groups_to_remove;

  // 1. making a set of all region names
  for ( typename map<string, csmp::Region<dim> >::const_iterator
        grit = UniqueRegionsBegin(); grit != UniqueRegionsEnd(); grit++ )
    region_names.insert( (*grit).first );

  // 2. For all regions whose name does not contain any numbers,
  //    find subregions identified by numbers attached to their names
  for ( typename map<string, csmp::Region<dim> >::const_iterator
        grit = UniqueRegionsBegin(); grit != UniqueRegionsEnd(); grit++ )
    if ( (*grit).first.find( target ) == string::npos )
    {
      // checking whether the name contains a number
      bool hasnumber = false;
      for ( string::const_iterator sit = (*grit).first.begin(); sit != (*grit).first.end(); sit++ )
        if ( isdigit( *sit ) ) { hasnumber = true; break; }
      // if not it is assumed that this is a primary region
      if ( !hasnumber ) {
        // and subregions are searched for in the region_name set
        for ( set<string>::const_iterator it = region_names.begin(); it != region_names.end(); it++ )
          // if the region_name contains the search string, this is a subregion to be deleted
          if ( (*it).find( target ) != string::npos )
            groups_to_remove.insert( (*it) );

        // the original region, however is kept by removing its name from the deletion list
        groups_to_remove.erase( (*grit).first );
      }
    }


  // Any region can only be unique or non unique. If no unique region was found the non-unique ones are searched
  if ( groups_to_remove.empty() ) {
    for ( typename map<string, csmp::Region<dim> >::const_iterator
          grit = RegionsBegin(); grit != RegionsEnd(); grit++ )
      region_names.insert( (*grit).first );

    // 2. For all regions whose names do not contain any numbers,
    //    find subregions identified by numbers attached to their names
    for ( typename map<string, csmp::Region<dim> >::const_iterator
          grit = RegionsBegin(); grit != RegionsEnd(); grit++ )
      if ( (*grit).first.find( target ) == string::npos )
      {
        // checking whether the name contains a number
        bool hasnumber = false;
        for ( string::const_iterator sit = (*grit).first.begin(); sit != (*grit).first.end(); sit++ )
          if ( isdigit( *sit ) ) { hasnumber = true; break; }
        // if not it is assumed that this is a primary region
        if ( !hasnumber ) {
          // and subregions are searched for in the region_name set
          for ( set<string>::const_iterator it = region_names.begin(); it != region_names.end(); it++ )
            // if the region_name contains the search string, this is a subregion to be deleted
            if ( (*it).find( target ) != string::npos )
              groups_to_remove.insert( (*it) );

          // the original region, however is kept by removing its name from the deletion list
          groups_to_remove.erase( (*grit).first );
        }
      }
  }

  // removing the original region from the removal list, just in case it is contained therein
  groups_to_remove.erase( target );

  // 3. removing the subregions and extra regions
  size_t  groups_removed( 0 );
  if ( !groups_to_remove.empty() ) cout << "\nRegionsInterface<dim,REGION_COMPLEX>::RemoveRegionPartitionsFor: removing region(s): ";
  for ( set<string>::const_iterator it = groups_to_remove.begin(); it != groups_to_remove.end(); it++ ) {
       cout << "'" << (*it) << "' ";
       const bool also_remove_elmts{ false };
       RemoveRegion( (*it).c_str(), also_remove_elmts );
       groups_removed++;
    }
  if ( !groups_to_remove.empty() ) cout << endl << endl;

  return groups_removed;

} // end RemoveRegionPartitionsFor






/**
FormAndAddRectangularRegion() builds a non-unique region from those elements
whose node coordinates lie within a user-specified bounding box. This box is
defined by its lower left and upper right corners.

When a new region is formed a region internal flag will be assigned to
each node, constraint point and element.
There are two region-internal object flags, PLAIN and BOUNDARY. When a
new region is formed, the Region method IdentifyBoundaryAs() assigns the
region-internal object flags, depending on whether the nodes, constraint
points or elements in the region lie at the region boundary or inside of
the region. Elements are assigned a boundary flag if at least one of
their faces coincides with the region boundary.

@section arguments Input Arguments

Assigns the finite elements that are located in a rectangular region that
is given by the double arguments (x = horizontal, y = vertical) as
members of a new region with the name char* (first argument).

@section implementation Implementation

The method works only for two-dimensional models and it calls the method
AccumulateRectangularRegion() of the Region class.

@section application Application

With FormAndAddRectangularRegion() regions to monitor changing properties
in a model can be defined at exact positions. Using property minimum
and maximum values in such regions, one can also measure property
gradients at positions of interest.

@section messages Messages

If the region cannot be formed because no elements in the desired
coordinate range can be found or because another region with the same
name already exists, the method will terminate the program, by reporting
a fatal error.

@attention this Region is always non-unique because it overlaps other regions

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::FormRectangularRegion( const char* region_name,
                                                                    const Point<dim>& min_xyz,
                                                                    const Point<dim>& max_xyz )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !HasValidModelRegion() )
    csmp_error.Note( FATAL_ERROR, "RegionInterface<dim, REGION_COMPLEX>::FormRegionsFromPropertyValues:",
                      "method relies on the existence of region 'Model', which does not exist");

  string output_region( region_name );
  if ( ContainsRegion( region_name ) ) {
      cout << "\nModel<" << dim << ">::FormRectangularRegion: WARNING: region '" << region_name;
      cout << "' already exists, adding an underscore at end of name: ";
      output_region += "_";
      cout << output_region << endl;
    }

  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) );
  // the 'bool' member of pair indicates whether insertion into map worked or not
  pair<typename map<string, csmp::Region<dim> >::iterator, bool>
    it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, regionComplex->Database() ) ) );
  if ( it.second )
    {
       (*it.first).second.AccumulateRectangularRegion( regionComplex->Mesh(), min_xyz, max_xyz );

       // removing the region if it contains no elements
       if ( (*it.first).second.Cells() == 0U ) {
          regionMap_.erase( it.first );
          csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::FormRectangularRegion",
                             "Region could not be formed", output_region.c_str() );
          return 0U;
        }
    }
  else {
    csmp_error.Note( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::FormRectangularRegion",
                       "Region could not be formed. Does this region already exist?", output_region.c_str() );
    return 0U;
  }

  // assigning new region name
  (*it.first).second.Name( region_name );

  return (*it.first).second.Cells();

} // end FormRectangularRegion 



/**
Makes a copy of an existing region and stores it under a new name.
NB: The new region will not be unique.

@param existing_region The names of the existing and the new region that shall be created.

@attention since this region will be a copy of an existing one, it will be non-unique

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void  RegionInterface<dim, REGION_COMPLEX>::CopyRegion( const char* existing_region, const char* new_copied_region, bool unique_region )
{
  if ( !ContainsRegion( existing_region ) )
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::CopyRegion",
                           "region to make a copy of could not be found",
                           existing_region );

  string output_region( new_copied_region );
  if ( ContainsRegion( new_copied_region ) ) {
    cout << "\nModel<" << dim << ">::CopyRegion: WARNING: region '" << new_copied_region;
    cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    cout << output_region << endl;
  }

  const csmp::Region<dim>&  gr_ref{ Region( existing_region ) };

  // irrespective of whether the original region was unique or non-unique its copy
  // will not be unique because it overlaps with the original region
  pair<typename map<string, csmp::Region<dim> >::iterator, bool>
    it = (unique_region) ? uniqueRegionMap_.insert( make_pair( output_region, csmp::Region<dim>( gr_ref ) ) ) :
    regionMap_.insert( make_pair( output_region, std::move( csmp::Region<dim>( gr_ref ) ) ) );
  if ( !it.second )
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::CopyRegion",
                           "region to copy to- could not be formed",
                           output_region.c_str() );

  // assigning new region name
  (*it.first).second.Name( new_copied_region );

} // end copyRegion





/**
Adds the first region to the second region so that the latter contains
both regions after the operation is complete.

@param region_to_add The names of the first region that shall be added to the second region.

@section messages Messages

ERRORs are reported if either the first or the second region does
not exist.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::AssimilateRegion( const char* region_to_add, const char* region_to_be_added_to )
{
  // finding the region in the region list
  try {
    const csmp::Region<dim>&  gref_to_add( Region( region_to_add ) );
    csmp::Region<dim>&        gref_to_be_added_to( Region( region_to_be_added_to ) );

    gref_to_be_added_to.Add( gref_to_add );
    
    // if the region was unique the elements added are assigned the same unique region ID
    if ( IsUnique(region_to_be_added_to) )
      gref_to_be_added_to.SetRegion_ID();
  }
  catch ( Exception& e ) {
    cout << "\nRegionsInterface<dim,REGION_COMPLEX>::AssimilateRegion: nothing was done; handled exception: " << endl;
    e.Out();
  }

} // end AssimilateRegion







/** Merges regions supplied as argument set into a single region with a new name.

@param input_regions is a set of unique names of the regions that
shall be merged into the
@param ensemble_region result region (second argument).

@section application Application

To join dynamically created regions during a simulation.

@section messages Messages

If the set of regions is empty or if one of the specified regions does not
exist, an error or an info message is reported.

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::MergeRegions( const set<string>& input_regions, const char* ensemble_region )
{
  string output_region( ensemble_region );

  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  if ( ContainsRegion( ensemble_region ) ) {
    csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                       ensemble_region, "output region already exists, adding an underscore to its name." );
    output_region += "_";
    cout << output_region << endl;
  }

  if ( input_regions.empty() ) {
    csmp_error.Note( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                       "No input regions were specified; merge could not be performed; no new region.",
                       ensemble_region );
    return;
  }

  typename map<string, csmp::Region<dim> >::const_iterator  iter;

  // collecting element indexes from input regions into set for output
  bool all_regions_are_unique{ true };
  vector<Element<dim>*>  element_ptrs;
  for ( auto it = input_regions.begin(); it != input_regions.end(); it++ ) {
    // finding the region in the region list
    if ( (iter = regionMap_.find( *it )) != regionMap_.end() or
         (iter = uniqueRegionMap_.find( *it )) != uniqueRegionMap_.end() ) {
      // outputting the ids of the member elements of the region
      element_ptrs.reserve( element_ptrs.size() + (*iter).second.Cells() );
      for ( auto eit = (*iter).second.CellsBegin(); eit != (*iter).second.CellsEnd(); eit++ )
        element_ptrs.push_back( const_cast<Element<dim>*>(*eit) );
      // are all input regions unique?
      if ( regionMap_.find(*it) != regionMap_.end() )
        all_regions_are_unique = false;
    }
    else csmp_error.Note( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                            (*it).c_str(), "region does not exist and was therefore not considered." );
  }

  if ( !element_ptrs.empty() )
    {
      // eliminating duplicate entries from pointer vector
      sort( element_ptrs.begin(), element_ptrs.end() );
      element_ptrs.erase( unique( element_ptrs.begin(), element_ptrs.end() ), element_ptrs.end() );
      
      // rebuilding the inter-element connectivity
      //REGION_COMPLEX<dim>& model( static_cast<REGION_COMPLEX<dim>&>(*this) );
      //model.Mesh().template BuildConnectivity<Element>( element_ptrs.begin(), element_ptrs.end() );
      
      // making a non-unique new region
      pair<typename map<string, csmp::Region<dim> >::iterator, bool>
        it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region,
                               static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
      if ( !it.second )
        throw csmp::Exception( WARNING, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions:",
                               output_region.c_str(), "region could not be formed." );

      else (*it.first).second.Accumulate( element_ptrs.begin(), element_ptrs.end() );

      // if all of the merged regions are unique the elements of the new region is assigned a new unique region ID
      if ( all_regions_are_unique )
        (*it.first).second.SetRegion_ID();
   }
  else
    throw csmp::Exception( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::MergeRegions",
                           "No elements in target list; merged region could not be build",
                           output_region.c_str() );

  cout << "\nRegionsInterface<dim,REGION_COMPLEX>::MergeRegions: new region '" << output_region << "' formed successfully.\n";

} // end MergeRegions






template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::MergeRegions( const char* region_name_tag, const char* ensemble_region )
{
  set<string> regions;
  for ( regionConstIterator rit( UniqueRegionsBegin() ); rit != UniqueRegionsEnd(); ++rit )
    if ( rit->first.find( region_name_tag ) != string::npos )
      regions.insert( rit->first );
  if ( !regions.empty() )
    MergeRegions( regions, ensemble_region );
  return regions.size();
}







/**
If the first region a contains all the elements of the second region
b, the method will return the boolean variable 'true'.

@param groupa  a first region which is tested for whether it
contains all elements of a second region.

@return The method return either 'true' or 'false' depending on whether the first
region contains the second region or not.

@section application Application

RegionIncludes() is one of a set of boolean algebraic operations which are
supported for Region objects.

@section messages Messages

An error will be reported, if one of the evaluated regions does not
exist.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionIncludes( const char* groupa, const char* groupb ) const
{
  try {
    // finding the groups in the group list
    const csmp::Region<dim>&  a = Region( groupa );
    const csmp::Region<dim>&  b = Region( groupb );
    return a.Includes( b );
  }
  catch ( Exception& e ) {
    cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionIncludes: returning false; handled exception: " << endl;
    e.Out();
  }

  return false;

} // end RegionIncludes







/**
RegionUnion() will create a new region, defined as the union of two existing
regions.

@section arguments Input Arguments

The names of the two existing regions are supplied as first and second
method arguments. The third argument specifies the name of the new region
which will contain both the first and the second region.

@section application Application

To build regions on the basis of complex criteria, for instance, one could
form the union of a region which includes all elements that are hotter than
670oC with a region representing a granite melt in the model.

@section messages Messages

If either one of the regions is empty or does not exist or if the target
region cannot be formed an error is reported.

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionUnion( const char* groupa, const char* groupb,
                                                        const char* groupunion )
{
  string output_region( groupunion );
  if ( ContainsRegion( groupunion ) ) {
    cout << "\nModel<" << dim << ">::RegionUnion: WARNING: region '" << groupunion;
    cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    cout << output_region << endl;
  }

  // finding the regions in the region list
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding new region
    // the 'bool' member of pair indicates whether insertion into map worked or not
    pair<typename map<string, csmp::Region<dim> >::iterator, bool>
      it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second )
    {
      groupUnion( itera, iterb, (*it.first).second );

      // removing the region if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Cells() == 0U ) {
        regionMap_.erase( it.first );
        return false;
      }
      // assigning new region name
      (*it.first).second.Name( groupunion );
      
      // new region ID is assigned to the ensemble of elements
      if ( IsUnique(groupa) && IsUnique(groupb) )
        (*it.first).second.SetRegion_ID();
    }
    else {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      csmp_error.Note( ERROR, "RegionsInterface<dim,REGION_COMPLEX>::RegionUnion",
                         "union of regions cannot be build", output_region.c_str() );
      return false;
    }

    return true;
  }
  catch ( Exception& e ) {
    cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionUnion: returning false; handled exception: " << endl;
    e.Out();
  }

  return false;

} // end RegionUnion








/**
RegionIntersection() finds those finite-elements which belong both to a
first and a second region and forms a third region from them.

Method returns true if there is an intersection and false if none
can be found.

@section arguments Input Arguments

The names of the two existing regions are supplied as first and second
method arguments. The third argument specifies the name of the new region
which will contain some members of the first and of the second region.

@section application Application

RegionIntersection() allows the user to define regions on the basis of
multiple property criteria. For instance, the user may form a region of
elements which represent granite above a temperature of 500oC. This would
be achieved by intersecting a 'granite' region with an 'above-500oC'
region.

@section messages Messages

If either the first or the second region does not exist or if the
desired intersection region would have zero elements, an error is
reported.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionIntersection( const char* groupa,
                                                               const char* groupb,
                                                               const char* groupintersection )
{
  string output_region( groupintersection );
  if ( ContainsRegion( groupintersection ) ) {
    cout << "\nModel<" << dim << ">::RegionIntersection: WARNING: region '" << groupintersection;
    cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    cout << output_region << endl;
  }

  // finding the regions in the region list
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding new region
    // the 'bool' member of pair indicates whether insertion into map worked or not
    pair<typename map<string, csmp::Region<dim> >::iterator, bool>
      it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second ) {
      intersection( itera, iterb, (*it.first).second );
      // removing the region if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Cells() == 0U ) {
        regionMap_.erase( it.first );
        return false;
      }
      // assigning new region name
      (*it.first).second.Name( groupintersection );
    }
    else return false;

    return true;
  }
  catch ( Exception& e ) {
    cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionIntersection: returning false; handled exception: " << endl;
    e.Out();
  }

  return false;

} // end RegionIntersection







/**
RegionDifference() forms a new region which will contain those finite
elements of the first region, which are not contained in the second
region.

If the difference does not exist false is returned.

@section arguments Input Arguments

The names of the two existing regions are supplied as first and second
method arguments. The third argument specifies the name of the new region.


@section application Application

To identify for instance what distinguishes a first from a second
region.

@section messages Messages

If either the first or the second region does not exist or if the
desired distinction region would have zero elements, an error is
reported.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionDifference( const char* groupa, const char* groupb,
                                                             const char* groupdiff )
{
  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  string output_region( groupdiff );
  if ( ContainsRegion( groupdiff ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RegionDifference: '", groupdiff,
                       "' already exists, adding an underscore at end of name." );
    output_region += "_";
    cout << output_region << endl;
  }

  // finding the regions in the map
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding the new region
    // the 'bool' member of pair indicates whether insertion into map worked or not
    pair<typename map<string, csmp::Region<dim> >::iterator, bool>
      it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region,
                             static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second )
    {
      size_t elements_of_new_region = difference( itera, iterb, (*it.first).second );
      if ( elements_of_new_region == 0 )
        csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RegionDifference:", groupdiff,
                           "to be formed. No elements found that belong only to one of the 2 input regions." );

      // removing the region if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Cells() == 0 ) {
        regionMap_.erase( it.first );
        return false;
      }

      // assigning new region name
      (*it.first).second.Name( groupdiff );
    }

    return true;
  }
  catch ( Exception& e ) {
    cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionDifference: returning false; handled exception: " << endl;
    e.Out();
  }

  return false;

} // end RegionDifference









/**
RegionSymmetricDifference() forms a new region which will contain those finite
elements of a first and a second region, which are not contained in both
regions.

@section arguments Input Arguments

The names of the two existing regions are supplied as first and second
method arguments. The third argument specifies the name of the new region.


@section application Application

To identify for instance what two regions do not have in common.

@section messages Messages

If either the first or the second region does not exist or if the
desired distinction region would have zero elements, an error is
reported.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RegionSymmetricDifference( const char* groupa,
                                                                      const char* groupb,
                                                                      const char* groupsymdiff )
{
  string output_region( groupsymdiff );
  if ( ContainsRegion( groupsymdiff ) ) {
    cout << "\nModel<" << dim << ">::RegionSymmetricDifference: WARNING: region '" << groupsymdiff;
    cout << "' already exists, adding an underscore at end of name: ";
    output_region += "_";
    cout << output_region << endl;
  }

  // finding the regions in the region list
  try {
    const csmp::Region<dim>&  itera = Region( groupa );
    const csmp::Region<dim>&  iterb = Region( groupb );

    // adding new region
    // the 'bool' member of pair indicates whether insertion into map worked or not
    pair<typename map<string, csmp::Region<dim> >::iterator, bool>
      it = regionMap_.insert( make_pair( output_region, csmp::Region<dim>( output_region, static_cast<REGION_COMPLEX<dim>*>(this)->Database() ) ) );
    if ( it.second )
    {
      symmetricDifference( itera, iterb, (*it.first).second );
      // removing the region if it contains no elements (extra error message is generated in function)
      if ( (*it.first).second.Cells() == 0 ) {
        regionMap_.erase( it.first );
        return false;
      }
      // assigning new region name
      (*it.first).second.Name( groupsymdiff );
    }
    else return false;

    return true;
  }
  catch ( Exception& e ) {
    cerr << "\nRegionsInterface<dim,REGION_COMPLEX>::RegionSymmetricDifference: returning false; handled exception: " << endl;
    e.Out();
  }

  return false;

} // end RegionSymmetricDifference







/**
Removes the elements shared with 'region-to-subtract' from the current (non-unique) region.

@note if the target region is unique, this implies that it does not overlap with any other region
except for non-unique regions.

@todo check (for volumetric regions) whether the removal of these elements has thus-far undiscovered consequences

@author SKM
@date 30/3/2016
@test OK
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RemoveFromRegion( const char* region, const char* region_to_subtract )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion( region ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                       region, "does not exist; nothing was done." );
    return false;
  }
  if ( !ContainsRegion( region_to_subtract ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                       region_to_subtract, "does not exist; nothing was done." );
    return false;
  }

  // finding the elements that are shared among the 2 regions
  csmp::Region<dim>&        r1_ref( Region( region ) );
  const csmp::Region<dim>&  r2_ref( Region( region_to_subtract ) );
  const size_t r1_elements( r1_ref.Cells() );
  const size_t r2_elements( r2_ref.Cells() );

  // new_region1 = region1 - region2
  vector<csmp::Element<dim>*>  new_region1;
  new_region1.reserve( r1_ref.Cells() );
  // if the element is not contained in region-to-subtract, it is kept
  for ( auto it = r1_ref.CellsBegin(); it != r1_ref.CellsEnd(); ++it )
    if ( !r2_ref.Contains( *it ) )
      new_region1.push_back( *it );

  // trimming excess capacity from new_region1
  vector<csmp::Element<dim>*>( new_region1 ).swap( new_region1 );

  // rebuilding the decimated region
  r1_ref.CellVector() = std::move( new_region1 );
  r1_ref.CreateNodePointerVector();
  r1_ref.IdentifyPerimeter();

  // reporting
  cout << "\nRegionInterface::RemoveFromRegion: removed region'" << region_to_subtract << "' (" << r2_elements << ")";
  cout << " from region '" << region << "' (" << r1_elements << ").\n";
  cout << "\t" << r1_ref.Cells() << " elements remaining in '" << region << "'.\n";

  return true;

} // RemoveFromRegion




/**
    This version removes the pointers to the supplied elements from the target region.
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::RemoveFromRegion( const char* region, const set<Element<dim>*>& elmt_set )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion( region ) ) {
      csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                         region, "does not exist; nothing was done." );
      return false;
    }
  if ( elmt_set.empty() ) {
      csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::RemoveFromRegion:",
                         "supplied element set was empty; nothing was done." );
      return false;
    }

  // finding the elements that are shared among the 2 regions
  csmp::Region<dim>&    subdomain( Region(region) );
  vector<Element<dim>*> elmts_to_remove( elmt_set.begin(), elmt_set.end() );
  
  const size_t elmts_removed = subdomain.RemoveRange( elmts_to_remove.begin(), elmts_to_remove.end() );
  
  // reporting
  cout << "\nRegionInterface::RemoveFromRegion: removed " << elmts_removed;
  cout << " elements from region '" << region << "'.\n";

  return true;

} // RemoveFromRegion





/**
Moves region to from the unique- to the non-unique regions map.

@return whether the operation was performed successfully (also true when the region is already
non unique.

@author SKM
@date 30/3/2016, 26/9/2022
@test OK
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::MoveToNonUniqueRegions( const char* unique_region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !ContainsRegion( unique_region ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::MoveToNonUniqueRegions:",
                       unique_region, "does not exist; nothing was done." );
    return false;
  }
  if ( !IsUnique( unique_region ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::MoveToNonUniqueRegions:",
                       unique_region, "is already a non-unique region; nothing was done." );
    return true;
  }

  // performing the move
  auto region_handle = uniqueRegionMap_.extract( unique_region );
  regionMap_.insert( std::move(region_handle) );

  return true;

} // end MoveToNonUniqueRegions



/**
Numbers the unique regions of the models, labeling their elements with the region number
as "region identifier". If the supplied variable does not exist, it is created by this method.

@param region_identifier scalar element property the unique value of which shall be used to distinguish the unique regions
@param region_names vector of region names to retrieve them from the integer keys

@attention this functionality is superseded by the unique RegionID stored on each element of each unique region

*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::CountAndLabelUniqueRegions( const char* region_identifier, vector<string>& region_names )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );

  // setting element variable up to identify all unique regions - uniquely
  if ( regionComplex.Database().IsDefined( region_identifier ) ) {
       if ( regionComplex.Database().Placement(region_identifier) != ELEMENT ||
            regionComplex.Database().Type(region_identifier)      != SCALAR  )
         csmp_error.Note( ERROR, "RegionInterface<dim,REGION_COMPLEX>::CountAndLabelUniqueRegions:",
                          region_identifier, "this property must be a scalar placed on the element." );
    }
  else
    regionComplex.CreateProperty( region_identifier, "rid", "uint", SCALAR, ELEMENT );

  // counting the regions and initialising them with the unique identifiers
  region_names.clear();
  region_names.resize( UniqueRegions() );
  size_t regions{ 0U };
  for ( auto it = UniqueRegionsBegin(); it != UniqueRegionsEnd(); it++ ) {
    (*it).second.InputPropertyValue( region_identifier, makeScalar( PLAIN, static_cast<double>(regions) ) );
    region_names[regions] = (*it).first;
    regions++;
  }

  assert( regions == region_names.size() );
  return regions;

} // end countAndLabelRegions






template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
bool RegionInterface<dim, REGION_COMPLEX>::CheckRegionIdentifierIsUpToDate( const char* region_identifier )
{
  REGION_COMPLEX<dim>& regionComplex( static_cast<REGION_COMPLEX<dim>& >(*this) );

  //if region_identifier not yet defined
  if (regionComplex.Database().IsDefined( region_identifier ) == false )
    return false;

  //if region_identifier is defined, we need to check its values are up to date with the number of regions in the model
  ScalarVariable highest_material_id;
  (*(uniqueRegionMap_.rbegin()->second.CellsBegin()))->Read( regionComplex.Database().StorageKey( region_identifier), highest_material_id );

  //Highest material id should be smaller than the number of regions by 1
  bool approximatelyEqual = std::fabs(  static_cast<double>(this->UniqueRegions()) -  highest_material_id() - 1 ) < 0.1 ;

  return  approximatelyEqual;
}




/**
Finds the contact area between regions a and b, logging pairs of element pointers and face numbers; @return number of shared faces

Region a will be the inner region whose perimeter Element pointers will be the first in the pairs.

@atention there is similar functionality inside the MeshManager

@author SKM
@date 18/3/2017
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::SharedPerimeterFaces( const char* region_a, const char* region_b,
                                                                   vector<tuple<Element<dim>*, ///< inner element
                                                                   Element<dim>*, ///< outer element
                                                                   size_t,        ///< inner face
                                                                   size_t> >&     ///< outer face
                                                                   shared ) const
{
  const REGION_COMPLEX<dim>& regionComplex( static_cast<const REGION_COMPLEX<dim>& >(*this) );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !IsUnique( region_a ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::SharedPerimeterFaces:",
                       region_a, "is not a unique region; cannot proceed." );
    return true;
  }
  if ( !IsUnique( region_b ) ) {
    csmp_error.Note( WARNING, "RegionInterface<dim,REGION_COMPLEX>::SharedPerimeterFaces:",
                       region_b, "is not a unique region; cannot proceed." );
    return true;
  }

  const csmp::Region<dim>& subdomain_a( regionComplex.Region( region_a ) );
  const csmp::Region<dim>& subdomain_b( regionComplex.Region( region_b ) );

  if ( !shared.empty() ) shared.clear();
  shared.reserve( subdomain_a.PerimeterCells() );

  // 1. searching for shared faces between the regions
  for ( size_t eid = subdomain_a.InteriorCells(); eid<subdomain_a.Cells(); ++eid )
    for ( auto face = 0U; face<subdomain_a.PerimeterFaces( eid ); ++face ) {
      // checking whether the neighbor of the perimeter face is in region b
      // ------------------------------------------------------------------
      const auto   pface = subdomain_a.PerimeterFace( eid, face );
      Element<dim>* nptr = subdomain_a.E( eid )->Neighbor( pface );
      // avoiding searches for neighbors that do not exist because one is at the model boundary
      if ( nptr != nullptr and subdomain_b.IsPerimeterCell( nptr ) ) {
        // finding which face is the perimeter face in the neighbor element
        size_t opposite_pface( numeric_limits<size_t>::max() );
        for ( auto i{0U}; i<nptr->Faces(); ++i )
          if ( nptr->Neighbor( i ) == subdomain_a.E( eid ) ) {
            opposite_pface = i;
            break;
          }
        assert( opposite_pface < nptr->Faces() );
        // creating connection record
        tuple<Element<dim>*, Element<dim>*, size_t, size_t>
          connection( subdomain_a.E( eid ), nptr, pface, opposite_pface );
        // storing the information
        shared.emplace_back( connection );
      }
    }

  // return how many shared faces were found
  return shared.size();

} // end SharedPerimeterFaces






/**
     Rebuilds modified regions. Region 'Model' will always be rebuilt.
     Any model subdomain that has been ScheduledForRebuilt()  will also be rebuilt, whether unique or non-unique.
      
      @attention the assumption is made the inter-element connectivity has was updated before
*/
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::UpdateRegions()
 {
     // since this region may now contain a different number of elements
     const bool also_remove_elmts{ false };
     RemoveRegion("Model", also_remove_elmts );

     // 1. rebuilding the region 'Model'
     // --------------------------------
     const bool is_unique = ( distance(UniqueRegionsBegin(), UniqueRegionsEnd()) > 0 ) ? false : true;
     FormModelRegion( is_unique );

     csmp::Region<dim>&  model_domain = RegionInterface<dim,REGION_COMPLEX>::Region("Model");

     // 2. non-unique, potentially overlapping regions
     // ----------------------------------------------
     //    (they are rebuilt using original creation constraints)
     for ( auto rit=RegionsBegin(); rit!=RegionsEnd(); ++rit )
       if ( (*rit).second.NeedsRebuilt() && (*rit).first != "Model" ) {
             auto crit = regionTraits_.find( (*rit).first );
             PropertyConstraints region_traits = ( crit == regionTraits_.end() )
                                                    ? PropertyConstraints("permeability", 1e-21,1e-5) : (*crit).second;
                                                    
             (*rit).second.UpdateCellMembershipApplyingConstraints( model_domain.CellsBegin(), model_domain.CellsEnd(), region_traits );
             
             // assuming the the element neighbor connectivity was updated before by the MeshManager
             (*rit).second.RebuildSubDomainAfterChangeOfCellVector();
          }

     // 3. unique regions: only get modified if they have been ScheduledForRebuilt()
     // ----------------------------------------------------------------------------
     for ( auto rit=UniqueRegionsBegin(); rit!=UniqueRegionsEnd(); ++rit )
       if ( (*rit).second.NeedsRebuilt() ){
         // assuming the the element neighbor connectivity was updated before by the MeshManager
         (*rit).second.RebuildSubDomainAfterChangeOfCellVector();
        std::cout << "Rebuilt Subdomain " << (*rit).first << " -> element and node vector now up to date" << std::endl;
       }
       
 } // end UpdateRegions






template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::RegionsOut() const
 {
     cout <<"\n\nRegionInterface<"<< dim <<",Region<Element>>::RegionsOut:\n";
     cout <<"\n\tUnique regions of model:\n";
     for ( auto rit=UniqueRegionsBegin(); rit!=UniqueRegionsEnd(); ++rit ) {
          cout <<"\t\t"<< (*rit).first <<":"<< (*rit).second.DomainIndex();
          cout <<" "<< (*rit).second.Cells() <<" elements,";
          pair<int32_t, int32_t> rdim = (*rit).second.ElementSpatialDimensions();
          if ( rdim.second == 3U )
            cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 2U )
            cout <<" surface area (m2): "<< (*rit).second.Volume() <<", perimeter length (m): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 1U )
            cout <<" length (m): "<< (*rit).second.Volume();
          cout <<", range of element spatial dimensions: "<< rdim.first <<", highest spatial dimension "<< rdim.second << endl;
       }
     cout <<"\n\tNon-unique regions of model:\n";
     for ( auto rit=RegionsBegin(); rit!=RegionsEnd(); ++rit ) {
          cout <<"\t\t"<< (*rit).first;
          cout <<" "<< (*rit).second.Cells() <<" elements,";
          pair<int32_t, int32_t> rdim = (*rit).second.ElementSpatialDimensions();
          if ( rdim.second == 3U )
            cout <<" volume (m3): "<< (*rit).second.Volume() <<", surface area (m2): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 2U )
            cout <<" surface area (m2): "<< (*rit).second.Volume() <<", perimeter length (m): "<< (*rit).second.SurfaceArea();
          else if ( rdim.second == 1U )
            cout <<" length (m): "<< (*rit).second.Volume();
          cout <<", range of element spatial dimensions: "<< rdim.first <<", highest spatial dimension "<< rdim.second << endl;
       }
     cout << endl << endl;
     
     return distance( RegionsBegin(), RegionsEnd() );
 }



    /// prints ModelSubDomain::domain_idx_ assigned automatically using reference count mechanism
template<uint32_t dim, template<uint32_t> class REGION_COMPLEX>
void RegionInterface<dim, REGION_COMPLEX>::PrintDomainIndices() const
 {
     cout <<"\n"<<"RegionInterface<dim, REGION_COMPLEX>::PrintDomainIndices:"<< endl;
     cout <<"\t"<<"indices of unique regions:"<< endl;
     for ( auto rit=UniqueRegionsBegin(); rit!=UniqueRegionsEnd(); ++rit )
       cout <<"\t\t"<< (*rit).first <<": "<< (*rit).second.DomainIndex() << endl;
     cout <<"\t"<<"indices of non-unique regions:"<< endl;
     for ( auto rit=RegionsBegin(); rit!=RegionsEnd(); ++rit )
       cout <<"\t\t"<< (*rit).first <<": "<< (*rit).second.DomainIndex() << endl;
     cout << endl;
 }




// EXPLICIT TEMPLATE INSTANTIATION FOR REGION_INTERFACE
template class RegionInterface<1U, Model>;
template class RegionInterface<2U, Model>;
template class RegionInterface<3U, Model>;




/**
   In preparation for the generation of a SplitBoundary at a region boundary,
   search the elements outside that do not have a face but some nodes on its perimeter.
   
   @param perimeter_nodes input  that can be computed with sharedPerimeterNodes()
   @param perimeter_cells elements on the perimeter (inside) of the region, and those matching on the outside
   @param touching_elmts those outside elements that contain the perimeter node but do not have a face on the boundary
   
   @return how many of such nodes were found.
   
   TODO: not used at the moment, remove?
   TODO: generalise to all types of ModelSubDomains
  
*/
template<uint32_t dim>
size_t outsideElementsWithNodesTouchingPerimeter( const Region<dim>& subdomain,
                                                  const vector<Node<dim>*>& perimeter_nodes,
                                                  const vector<pair<pair<Element<dim>*,uint32_t>,
                                                               pair<Element<dim>*,uint32_t> > >& perimeter_cells,
                                                  map<Node<dim>*,map<Element<dim>*,uint32_t>>& touching_elmts )
 {
    touching_elmts.clear();
    
    // creating a search vector with outside cells that have a face on the regions perimeter
    // (since these are already unique no check is required)
    vector<Element<dim>*> outside_perimeter_cells;
    outside_perimeter_cells.reserve( perimeter_cells.size() );
    for ( const auto& it : perimeter_cells )
      outside_perimeter_cells.push_back( it.second.first );
    sort( outside_perimeter_cells.begin(), outside_perimeter_cells.end() );
    
    // search perimeter nodes for parent elements that are not perimeter elements
    // (as determined by searching the elements that share perimeter face)
    for ( const auto& nit : perimeter_nodes )
      for ( uint32_t i{0u}; i<nit->Parents(); ++i )
        // record ‘eptr’ and local node number (from with node-pointer can be deduced for later lookup)
        if ( !binary_search( outside_perimeter_cells.begin(), outside_perimeter_cells.end(), nit->Parent(i) ) &&
             !subdomain.Contains( nit->Parent(i) )  ) {
             auto insert_it = touching_elmts.insert( make_pair( nit,
                                                                map<Element<dim>*,uint32_t>{ {nit->Parent(i),
                                                                                              nit->ParentNodeNumber(i)} } ) );
             // if there already is an entry for this node another map entry is inserted for it
             if ( !insert_it.second )
               (*insert_it.first).second.insert( make_pair( nit->Parent(i), nit->ParentNodeNumber(i) ) );
          }

    // return how many of such outside elements were found
    size_t n_extra_outside_elmts{0ul};
    for ( const auto& nit : touching_elmts )
      n_extra_outside_elmts += nit.second.size();
      
    return n_extra_outside_elmts;
 
 } // outsideElementsWithNodesTouchingPerimeter

template size_t outsideElementsWithNodesTouchingPerimeter( const Region<1U>&,
                                                           const vector<Node<1U>*>&,
                                                           const vector<pair<pair<Element<1U>*,uint32_t>,
                                                                             pair<Element<1U>*,uint32_t> > >&,
                                                           map<Node<1U>*,map<Element<1U>*,uint32_t>>& );

template size_t outsideElementsWithNodesTouchingPerimeter( const Region<2U>&,
                                                           const vector<Node<2U>*>&,
                                                           const vector<pair<pair<Element<2U>*,uint32_t>,
                                                                             pair<Element<2U>*,uint32_t> > >&,
                                                           map<Node<2U>*,map<Element<2U>*,uint32_t>>& );

template size_t outsideElementsWithNodesTouchingPerimeter( const Region<3U>&,
                                                           const vector<Node<3U>*>&,
                                                           const vector<pair<pair<Element<3U>*,uint32_t>,
                                                                             pair<Element<3U>*,uint32_t> > >&,
                                                           map<Node<3U>*,map<Element<3U>*,uint32_t>>& );





/**
      Returns those Elements from outside Region 2, which share a face or touch (with a node) the perimeter of the inside Region 1.
      All these elements are returned into an unordered (hash) map for fast access.
      
      @param region1 unique model subdomain that is supposed to have  a shared border with subdomain2
      @param region2 unique model subdomain bordering subdomain 2
      @param halo_elmts  into which the cells forming a halo to subdomain1 will be stored
      
      @attention this method only works for Region because nodes do not store Face or InterFace parents
      
      @note to find the outside elements ModelSubDomain<>::Contains() is used; this is expensive but always works
*/
template<uint32_t dim>
size_t haloElements( const Region<dim>& region1, const Region<dim>& region2, unordered_set<Element<dim>*>& halo_elmts )
 {
    // finding the nodes that are shared between the subdomains which will be on their perimeter
    vector<Node<dim>*> perimeter_nodes;
    if ( sharedPerimeterNodes( region1, region2, perimeter_nodes ) == 0 )
      return 0ul;
    
    // finding the halo cells by searching the parent elements of the perimeter nodes
    halo_elmts.clear();
    for ( const auto& nit : perimeter_nodes )
      for ( uint32_t i{0u}; i<nit->Parents(); ++i ) {
            if ( nit->Parent(i) == nullptr ) {
                 cerr <<" nd:"<< nit->Idx() <<" parent:"<< i <<"null";
                 continue;
              }
            if ( region2.Contains( nit->Parent(i) ) )
              halo_elmts.insert( nit->Parent(i) );
        }
      
    return halo_elmts.size();
 
 } // haloCells
 
template size_t haloElements( const Region<1U>&, const Region<1U>&, unordered_set<Element<1U>*>& );
template size_t haloElements( const Region<2U>&, const Region<2U>&, unordered_set<Element<2U>*>& );
template size_t haloElements( const Region<3U>&, const Region<3U>&, unordered_set<Element<3U>*>& );





} // csmp
