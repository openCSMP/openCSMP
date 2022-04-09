#include "SplitBoundaryInterface.h"
#include "ModelTopology.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"
#include "MeshManagementUtilities.h"
#include "smoothElementData.h"
#include "MeshManager.h"
#include "Node.h"
#include "NodeManifold.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "IsoparametricLinearLineElement.h" 
#include "IsoparametricLinearTriangle.h" 
#include "IsoparametricLinearQuadrilateral.h"
#include "binaryReadWrite.h"

#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {


/**
returns reference to SplitBoundary
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
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
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
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
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::ContainsSplitBoundary( const std::string& bname ) const
{
  if ( splitBoundaryMap_.find( bname ) != splitBoundaryMap_.end() ) return true;
  return false;
}


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin()
{ return splitBoundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd()
{ return splitBoundaryMap_.end(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin() const
{ return splitBoundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd() const
{ return splitBoundaryMap_.end(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
size_t  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundaries() const
{ return splitBoundaryMap_.size(); }




/**
Removes the SplitBoundary object with the given name if it exists inside the model; else returns
with a warning.

Removing a SplitBoundary will also prompt the MeshManager to  delete the corresponding InterFace elements
from the model, potentially turning it into a disconnected group of mesh patches.

@note he mesh is not fused together again.

@todo Provide the option to fuse the mesh back together again.

       @author SKM (refactored - since design was flawed)
       @date 15/8/2020
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( const char* splitboundary )
 {
    if ( !ContainsSplitBoundary(splitboundary) ) {
         ErrorHandler::Instance().notice( WARNING, "SplitBoundaryInterface::RemoveSplitBoundary",
                                          splitboundary, "split boundary is not contained in model; nothing was done");
         return;
      }
      
    // locating the boundary in the split boundary map
    splitBoundaryIterator  spit = splitBoundaryMap_.find( splitboundary );
    if ( spit != splitBoundaryMap_.end() ) {
         // deleting the split boundary
         splitBoundaryMap_.erase( (*spit).first );
         
         // TODO: remove duplicated nodes again
         cout <<"\n"<<"SplitBoundaryInterface<"<< dim <<">::RemoveSplitBoundary: removed split boundary named '";
         cout << splitboundary <<"' albeing retaining duplicated nodes."<< endl;
         return;
      }

    // if the split boundary was not contained in the map, feedback is given
    ErrorHandler::Instance().notice( WARNING, "SplitBoundaryInterface::RemoveSplitBoundary", splitboundary,
                                    "split boundary was not contained in SplitBoundary map");
  } // end RemoveSplitBoundary


/**
Removes the SplitBoundary object with the given name if it exists inside the model; else returns
with a warning.

Removing a SplitBoundary will also prompt the MeshManager to  delete the corresponding InterFace elements
from the model, potentially turning it into a disconnected group of mesh patches.

@note he mesh is not fused together again.

@todo Provide the option to fuse the mesh back together again.

       @author SKM (refactored - since design was flawed)
       @date 15/8/2020
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary )
 {
    // locating the boundary in the split boundary map
    splitBoundaryIterator  spit = splitBoundaryMap_.find( splitboundary.Name() );
    // if the addresses are the same
    if ( spit !=  splitBoundaryMap_.end() ) {
         // if the split boundary is contained in the map, it is erased
         splitBoundaryMap_.erase( (*spit).first );
         return;
      }

    ErrorHandler::Instance().notice( WARNING, "SplitBoundaryInterface<dim,Model>::RemoveBoundary", "boundary does not exist" );
    
  } // end RemoveSplitBoundary


 
 


// -----------------------------------------------
// Binary input/output
// -----------------------------------------------

/**
   Writes SplitBoundary objects to file using the section tags:
   
        SBDFHEDR
        SPLITBDRY
        ONE_BDRY
        SBDRYVAR
        SBDFFOTR
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::OutputSplitBoundariesToBinary( const char* file_name ) const
{
  if ( this->SplitBoundaries() == 0 ) return false;

  const SPLITBOUNDARY_COMPLEX<dim>& splitBoundaryComplex( static_cast<const SPLITBOUNDARY_COMPLEX<dim>& >(*this) );

  string  bin_file( file_name );
  fstream fp( bin_file.c_str(), ios::out | ios::binary );
  if ( !fp.is_open() ) {
      throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary:", "Binary file could not be created." );
      return false;
    }
  // 0. writing file header
  {
     BinaryFileSectionWrite sect(fp, "SBDFHEDR");
     string heading("SplitBoundaryInterface::OutputSplitBoundariesToBinary: ");
     heading +="split boundary information for Model '";
     heading += splitBoundaryComplex.Name();
     heading +="' to file: ";
     heading += bin_file;
     heading +="'.";
     binaryFileWrite( fp, heading.c_str() );
  }
 cout <<"\nSplitBoundaryInterface<"<< dim <<">::OutputSplitBoundariesToBinary: split boundaries written to binary file: ";
     
  // 1. writing number of splitboundaries
  {
     BinaryFileSectionWrite sect(fp, "SPLITBDR");
  
     const uint64_t  records( this->SplitBoundaries() );
     fp.write( reinterpret_cast<const char*>(&records), sizeof(uint64_t) );

     for ( auto bit{ SplitBoundariesBegin() }; bit != SplitBoundariesEnd(); ++bit )
       {
          BinaryFileSectionWrite hdr(fp, "ONE_BDRY");
          cout <<"'"<< (*bit).first <<"' ";
          cout.flush();
          (*bit).second.WriteDomainIndexesToBinaryFile( fp );
          // NB: splitboundary objects have no BOX_BOUNDARY flag values because these always default to INTERNAL.
          // writing the stored variables
          domainVariablesOut( fp, (*bit).second, splitBoundaryComplex.Database() );
       }
   }

  // 2. footer and clean up
  BinaryFileSectionWrite sect(fp, "SBDFFOTR");

  fp.close();
  cout << "\n\nSplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary: Split boundaries have been successfully written to: '";
  cout << bin_file << "'" << std::endl;
  return true;

} // end OutputSplitBoundariesToBinary




/**
    Reads all split boundary objects in a model from binary file, expecting that interfaces and their connectivity are already present in MeshManager.
    @author SKM
    @date 6/9/2021
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InputSplitBoundariesFromBinary( const char* file_name,
                                                                                         const set<string>& subset_variables )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
  string bin_file( file_name );
  fstream fp( bin_file.c_str(), ios::in | ios::binary );
	 if (!fp.is_open()) {
          csmp_error.notice( ERROR, "SplitBoundaryInterface::InputSplitBoundariesFromBinary:",
                             bin_file, "file could not be opened; nothing was done." );
          return false;
       }

  // 0. reading the file header and printing it to the screen
  {
     BinaryFileSectionRead sect(fp, "SBDFHEDR");

     char  text[INFO_STRING];
     binaryFileRead( fp, text );
     cout <<"\nSplitBoundaryInterface<"<< dim <<">::InputSplitBoundariesFromBinary: Reading file header:\n\t"<< text << std::endl;
  }
  cout <<"\n\timporting the split boundaries: ";

  // 1. reading the spli boundaries
  SPLITBOUNDARY_COMPLEX<dim>& splitBoundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>& >(*this) );
  const PropertyDatabase<dim>& database( splitBoundaryComplex.Database() );
  MeshManager<dim>& mesh( splitBoundaryComplex.Mesh() );

  {
    BinaryFileSectionRead sect(fp, "SPLITBDR");
   
    uint64_t  records(0);  // region records
    // getting number of unique region records from file
    fp.read( reinterpret_cast<char*>(&records), sizeof(uint64_t ) );
    if ( records > 0 )
      // reading the regions sequentially
      for ( auto i{0U}; i<records; ++i )
        {
           BinaryFileSectionRead hdr(fp, "ONE_BDRY");

           // 1.1 reading name and face indices for each boundaries
           SubDomainInfo  info;
           readDomainIndexesFromBinaryFile( dim, fp, info );
          
           // 1.3 reading the split boundary objects
           pair<typename map<string,csmp::SplitBoundary<dim> >::iterator,bool>
             it=splitBoundaryMap_.insert( make_pair( info.name, csmp::SplitBoundary<dim>( database, mesh, info ) ) );
          
           if ( !it.second )
             throw csmp::Exception( FATAL_ERROR, "BoundaryInterface::InputBoundariesFromBinary:",
                                    info.name, "Boundary could not be formed; issue with binary file." );
       
           // 1.4 reading the values of variables placed on the SplitBoundary
           if ( subset_variables.empty() ) domainVariablesIn( fp, (*it.first).second, database );
           else selectedDomainVariablesIn( fp, (*it.first).second, database, subset_variables );
       
           // 1.5 reporting out
           cout <<"\n\t\t"<< (*it.first).first << (*it.first).second.Cells() <<" faces).";
       }
   }

  // 2. cleaning up
  BinaryFileSectionRead sect(fp, "SBDFFOTR");

  fp.close();
  cout <<"\n\nSplitBoundaryInterface<"<< dim <<"InterFace>::InputSplitBoundariesFromBinary: file '";
  cout << bin_file <<"' has been read successfully.\n";

  return true;

} // end InputSplitBoundariesFromBinary


/*
  // 0. reading number of splitboundaries
  size_t records( 0 );
  fp.read( (char*)&records, sizeof( size_t ) );

  std::cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: reading " << csCache << " containing "
    << records << " splitboundaries\n" << std::endl;

  for ( size_t i( 0 ); i < records; ++i )
  {
    std::string bName;
    // reading name of splitboundary
    binaryFileRead( fp, csCache );
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

*/


/// container of juxtaposed element pairs for SplitBoundary creation:
template<uint32_t dim>
struct SplitBoundaryElementSets : public
  //       key consisting out the names of regions juxtaposed at the SplitBoundary
  std::map<std::pair<std::string, std::string>,
  // set that stores pointers to the pairs of elements juxtaposed across boundary
  // size_t parameter gives the local number of the element face that sits at the split boundary
  std::set<std::pair<std::pair<Element<dim>*, size_t>,
  std::pair<Element<dim>*, size_t> > > > {
};





/**
Creates a name like SPLITBOUNDARY_region1_region2 adding a number if this is necessary to make the
name unique.

@attention the facing relationships of the boundary are preserved so that region1 is the first in the
argument pair.

@author SKM
@date January 2018
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
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
    Creates SplitBoundaries detecting, them in the (ANSYS) input model
    as node-matched interfaces, connecting such disconnected perimeter element faces
    in mesh; these are detected and grouped by bordering regions and turned into SplitBoundary objects
    with names following the same conventions as for Boundary objects.
    
    @attention this method does not prompt the generation of any nodes because they are expected to be already there in ANSYS.
    
    @attention this method relies on the existing numbering of nodes, elements, faces and interfaces.

    @author SKM
    @date 20/03/2018 updated 5/9/2021

*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
std::pair<std::set<std::string>,bool> SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries()
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&               csmp_error( ErrorHandler::Instance() );

  // 1. detecting potential SplitBoundaries
  // --------------------------------------
  // finding the (local) ids of the element faces on either side of the split boundary
  set<pair<pair<Element<dim>*, size_t>, pair<Element<dim>*, size_t> > >  interface_elmt_pairs;
  if ( !findSplitInterfaceElements( splitboundaryComplex->Region( "Model" ), interface_elmt_pairs ) ) {
      csmp_error.notice( WARNING, "SplitBoundaryInterface::DetectAndCreateSplitBoundaries:",
                         "node-coordinate matched faces / internal boundaries could not be detected; nothing was done." );
      return make_pair(set<string>(),false);
    }

  // 2. classifying the detected interfaces in terms of the regions that they juxtapose
  // ----------------------------------------------------------------------------------
  // 2.1 assigning integer keys to the elements of the model regions so that their name can be identified from the perimeter elements
  //     names by their index integer
  vector<string>  region_names;
  if ( !splitboundaryComplex->Database().IsDefined( "region number" ) )
    splitboundaryComplex->CreateProperty( "region number", "-", SCALAR, ELEMENT );
  splitboundaryComplex->CountAndLabelRegions( "region number", region_names ); // method in RegionInterFace

  // 2.2 grouping interface element pairs into ones that juxtapose specific regions against one another
  //     these will later become specific split boundaries
  typedef set<pair<pair<Element<dim>*, uint32_t>, pair<Element<dim>*, uint32_t> > > INTERFACE_ELEMENT_PAIRS;
  map<pair<string, string>, INTERFACE_ELEMENT_PAIRS>  split_boundary_map;
  const csmp::Index reg_key( splitboundaryComplex->Database().StorageKey( "region number" ) );

  // for all the SplitBoundary objects supplied as sets of pairs of Element pointers and interface idx values
  for ( auto& iit : interface_elmt_pairs ) {
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
  for ( auto& i : split_boundary_map )
    cerr << i.first.first << "," << i.first.second << "\n";
  cerr << endl;
#endif

  // 3. Creating splitboundaries for each of the discovered juxtapositions of regions
  // --------------------------------------------------------------------------------
  pair<set<string>,bool>  splitBoundaryNames; 
  splitBoundaryNames.second = true;
  for ( auto& it : split_boundary_map ) {
      // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
      string bname = CreateSplitBoundaryName( it.first );
      splitBoundaryNames.first.insert( bname );
      
      // extracting the split boundary neighbor elements into interface element pairs
      InterFaceSet<dim>  ifset( it.second );

      // creating the SplitBoundary asking the MeshManager to create the required number of InterFace objects
      std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
        bit = splitBoundaryMap_.insert( std::make_pair( bname,
                                                        csmp::SplitBoundary<dim>( bname, splitboundaryComplex->Database(),
                                                                                  splitboundaryComplex->FE_Manager(),
                                                                                  splitboundaryComplex->Mesh(), ifset ) ) );
      if ( bit.second ) {
          std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries: '"<< bname;
          std::cout << "' created successfully.\n";
        }
      else
        throw csmp::Exception( INFO,
                               "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries:",
                               bname.c_str(),
                               "boundary already exists. Nothing was done." );
    }

  // reporting the names of the split boundaries that were created
  return splitBoundaryNames;

} // end DetectAndCreateSplitBoundaries




template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
size_t SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::FormSplitBoundariesFrom( const ModelTopology& topo )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. getting the names of the regions
  std::list<std::string> split_boundaries;
  topo.OutputSplitBoundaries( split_boundaries );

  // 2. assigning the regions to groups in the Model
  std::cout << "\nSplitBoundaryInterface::FormSplitBoundariesFrom: Forming the split boundaries: ";

  size_t new_split_boundaries{0};
  for ( auto lit = split_boundaries.begin(); lit != split_boundaries.end(); lit++ )
    {
      std::string domain_name( *lit );
      auto it = splitBoundaryMap_.insert( make_pair( domain_name, csmp::SplitBoundary<dim>( domain_name, static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this)->Database() ) ) );
      // if the region was successfully inserted
      if ( it.second )
        {
          // making a list of the element numbers
          std::vector<size_t>  cell_ids;
          cell_ids.reserve( topo.CellsWithinDomain( (*lit).c_str() ) );
          copy( topo.CellsOfDomainBegin( (*lit).c_str() ),
                topo.CellsOfDomainEnd( (*lit).c_str() ),
                back_inserter( cell_ids ) );

          // retrieving the elements by their IDs and assigning them  to the region
          (*it.first).second.AccumulateByNumber( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this)->Mesh(), cell_ids );
          cell_ids.erase( cell_ids.begin(), cell_ids.end() );

          // removing the group if it contains no elements
          if ( (*it.first).second.Cells() == 0U ) {
              splitBoundaryMap_.erase( it.first );
              csmp_error.notice( WARNING, "SplitBoundaryInterface::FormSplitBoundariesFrom",
                                 "SplitBoundary could not be formed", (*lit).c_str() );
            }
          else {
               // reporting the name of the newly generated region
               std::cout << domain_name << " ";
               new_split_boundaries++;
            }
        }
      else
        throw csmp::Exception( ERROR, "SplitBoundaryInterface::FormSplitBoundariesFrom",
                              "SplitBoundary could not be formed. Does this region already exist?", (*lit).c_str() );
    }
  std::cout << std::endl;

  return new_split_boundaries;

} // end FormSplitBoundariesFrom





/**
    Creates SplitBoundar(ies) from lower dimensional region without the need for a user to create Boundary objects first. The underlying region is removed in the process.
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
pair<set<string>,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom( const char* dim_1_region )
 {
    SPLITBOUNDARY_COMPLEX<dim>* model( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  
    // 1. Converting the lower dimensional region into a single Boundary or multiple Boundaries (patches of juxtaposed rocks)
    const bool remove_original_region(true);
    //    CreateInternalBoundaryFrom checks whether dim_1_region actually exists
    pair<set<string>,bool> boundary_names = model->CreateInternalBoundaryFrom( dim_1_region, remove_original_region );
    if ( boundary_names.second == false ) 
      return boundary_names;
    
    set<string>  split_boundary_names;
    for ( set<string>::const_iterator it=boundary_names.first.begin(); it!=boundary_names.first.end(); ++it ) {
          Boundary<dim>& boundary = model->Boundary( (*it) );
          // CreateSplitBoundaryFrom removes the boundary from which the split boundary was created
          split_boundary_names.insert( CreateSplitBoundaryFrom( boundary ).first );
       }
       
    return make_pair( split_boundary_names, true );
       
 } // end CreateSplitBoundaryFrom




/**
     CreateSplitBoundaryFrom( const Boundary<dim>& );
     
     @attention the input boundary is removed in the process.
     
     @return boolean indicating whether the method was able to create a singe Split boundary and its name
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
pair<string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom( Boundary<dim>& boundary )
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // checking whether boundary is external to the model in which case a SplitBoundary cannot be buid
  if ( boundary.IsExternal() ) {
       csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom", "Region 'Model' not eligible for InsertSplitBoundary." );
       return make_pair( "no SplitBoundary was created", false );
    }

  // creating boundary name by replacing BOUNDARY with SPLIT_BOUNDARY
  std::string   splitboundaryName( boundary.Name() );
  const size_t  str_length(string("BOUNDARY").length());
  splitboundaryName.replace( splitboundaryName.find("BOUNDARY"), str_length, "SPLIT_BOUNDARY" );

  // attempt to create a splitboundary
  if ( ContainsSplitBoundary(splitboundaryName) ) {
       csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom", splitboundaryName,
                        "a SplitBoundary with this name already exists; nothing was done." );
       return make_pair( "no SplitBoundary was created", false );
    }
  
  pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
    it = splitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,
                                                                                                splitboundaryComplex->Database() ) ) );
  if ( it.second ) {
      cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom:";
      cout <<" creating SplitBoundary from 'Boundary' "<< boundary.Name() << endl;
      (*it.first).second.CreateFrom( splitboundaryComplex->Database(), splitboundaryComplex->Mesh(), boundary );
    }

  // removes boundary also deleting its interface objects
  splitboundaryComplex->RemoveBoundary( boundary );

  cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom: created splitboundary: '";
  cout << splitboundaryName <<"' successfully.\n\n";

  return make_pair( splitboundaryName, false );

} // end CreateSplitBoundaryFrom( Boundary )









/**

Creates SplitBoundary via the creation of a boundary between region.
The original boundary gets removed.

@author SKM
@date 25/1/20

 */
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
pair<string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryBetween( const char* region1, const char* region2 )
 {
   SPLITBOUNDARY_COMPLEX<dim>*  modelComplex(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
   
   pair<string,bool> result = modelComplex->BoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateBoundaryBetween( region1, region2 );
   
   assert( result.second == true );
   return CreateSplitBoundaryFrom( modelComplex->Boundary(result.first) );

 } // end InsertSplitBoundary





/**
    @author E.P
    @date 07/2/2020

    @brief
    Creates lower dimensional elements and adds them to
    1) interface object as middle element
    2) Region of all the middle elements in the split boundary
    3) Model Region

    @note This method can be used either
    1) Create an entirely new lower dimensional Region within the SplitBoundary with its name followed by the string "_REGION".
    or
    2) Extends an already existing lower dimensional middle region of the split boundary to include any new interface objects which have the Middle_Element = nullptr.

    //E.P TODO: Setting boundary flags to Elements which have faces on a boundary must be done for 3D case!
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
std::pair<std::string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertRegionIntoSplitBoundary( const char* split_boundary,
                                                                                                                int32_t material_id )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( !ContainsSplitBoundary( split_boundary ) ) {
         csmp_error.notice( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertRegionIntoSplitBoundary",
                            split_boundary, "Does not exist; nothing was done." );
         return make_pair("no Region created",false);
       }
     
     SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
     csmp::SplitBoundary<dim>&    splitBoundary( model->SplitBoundary(split_boundary) );
     csmp::MeshManager<dim>&      mesh( model->Mesh() );
     
     // 1. creating unique set of nodes matching those on the inside of the SplitBoundary in position
     // --------------------------------------------------------------------------------------------------------------------------
     // NB: the use of Duplicate() ascertains that the MeshManager registers the new collocated nodes with the NodeManifoldManager
     set<Node<dim>*>     unique_new_nodes;
     vector<Node<dim>*>  node_pointers; // linear array of the nodes of one face after another
     node_pointers.reserve( splitBoundary.Nodes() );
     for ( auto it=splitBoundary.CellsBegin(); it!=splitBoundary.CellsEnd(); ++it ) {
          const size_t n_nodes((*it)->Nodes());
          // looping over the nodes on the inside of the interface which must be manifolds
          for ( auto i{0U}; i<n_nodes; ++i ) {
               assert( (*it)->N(i)->IsManifold() );
               pair<typename set<Node<dim>*>::iterator,bool> nit=unique_new_nodes.insert( (*it)->N(i) );
               // if the node is not yet contained in 'unique_new_nodes' it is created and added, but only if
               if ( nit.second ) {
                    // the corresponding NodeManifold does not already contain an isolated central node from another SB with intervening elements
                    // (this also covers the case of a SplitBoundary intersection where there should only be one intervening node in the middle)
                    if (  (*it)->N(i)->Manifold()->Branches() > 2U ) {
                         // finding intervening node, if any, and inserting it into the node pointer
                         Node<dim>* nptr = (*it)->N(i)->Manifold()->MIDDLE_Node();
                         if ( nptr != nullptr )
                           node_pointers.push_back( nptr );
                         else // a new node has to be generated
                           node_pointers.push_back( mesh.Duplicate( (*it)->N(i), MIDDLE ) );
                      }
                    //                            inserts the new node it creates into the corresponding NodeManifold
                    else node_pointers.push_back( mesh.Duplicate( (*it)->N(i), MIDDLE ) );
                 }
               // else the already created new node is added
               else node_pointers.push_back( (*nit.first) );
            }
       }
     unique_new_nodes.clear();
     
     
     // 2. creating elements within InterFace objects with node-numbering matching that of corresponding INNER parent element face
     // --------------------------------------------------------------------------------------------------------------------------
     const LocalVariables             element_props           = model->Database().LocalVariablesAt( ELEMENT );
     const IntegrationPointVariables  integration_point_props = model->Database().IntegrationPointVariablesAt( ELEMENT );
     vector<Element<dim>*>            new_elmts;
     new_elmts.reserve( splitBoundary.Cells() );
     size_t node_offset(0U); // for moving through the node-pointer vector
     
     for ( auto it=splitBoundary.CellsBegin(); it!=splitBoundary.CellsEnd(); ++it )
       {
          assert( (*it)->Parent(MIDDLE) == nullptr );
          // extracting element-node subvector
          const size_t n_nodes( (*it)->FE()->Nodes() );
          vector<Node<dim>*> nodes( &node_pointers[node_offset], &node_pointers[node_offset+n_nodes] );

          // construct the new element
          new_elmts.push_back(  model->Mesh().AddInterveningElement( (*it),
                                                                     element_props, integration_point_props,
                                                                     nodes, material_id ) );
          node_offset += n_nodes;
       }
       
     // 3. establishing neighbor connectivity among the new elements
     // ------------------------------------------------------------
     establishNeighborConnectivity( new_elmts, false, false );


     // 4. construct the new unique region between the interface elements in the model
     //    given it the same name as the split boundary but calling it region instead
     // -----------------------------------------------------------------------------
     string       region_name( splitBoundary.Name() );
     const size_t str_length( string("SPLIT_BOUNDARY").length() );
     region_name.replace( region_name.find("SPLIT_BOUNDARY"), str_length, "REGION" );
     
     const bool  unique_map(true);
     model->FormRegionFrom( region_name.c_str(), new_elmts.begin(), new_elmts.end(), unique_map );
     // add new unique region to model region
     model->Region("Model").Add( model->Region(region_name) );
     
     return make_pair( region_name, true ); 
     
 } // InsertRegionIntoSplitBoundary
 









/**
Insert a lower-dimensional element mesh into the model, which is fully disconnected from the model,
to represent a fracture that coincides with a SplitBoundary according to
the desire to use SplitBoundary objects for lower-dimensional fractures
for the modelling of both, fluid flow and geomechanics.
It only shares the locations of the nodes (if the SplitBoundary has not been modified) or
lies in the symmetry plane of the 2 sides if the nodes have been separated by deformation.

@author JCK updated 26/08/2019

TODO: like Boundary patches, the split boundary patches should have unique names that reflect juxtaposition relationships of regions

*/
/*
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SingleRegionFromAllSplitBoundaries( const char* region_name, int32_t material_id )
{
  SPLITBOUNDARY_COMPLEX<dim>*       splitboundaryComplex(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
  MeshManager<dim>&                 mesh    = splitboundaryComplex->Mesh();
  FiniteElementManager&             fem_mgr = splitboundaryComplex->FE_Manager();
  FiniteVolumeStencilManager<dim>*  fvm_mgr = splitboundaryComplex->FV_Manager();

  // 1. collect the nodes and the elements that belong to split boundaries
  set<Node<dim>*>         ifnodes;
  vector<InterFace<dim>*> ifelmts;
  
  for ( typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator
        it = splitboundaryComplex->SplitBoundariesBegin(); it != splitboundaryComplex->SplitBoundariesEnd(); ++it )
    {
      for (auto eit = (*it).second.CellsBegin(); eit != (*it).second.CellsEnd(); eit++) {
          csmp::InterFace<dim>* pInterFace = (*eit);
          ifelmts.push_back(pInterFace);
          for (size_t j = 0; j < pInterFace->Nodes() / 2; j++)
            ifnodes.insert(pInterFace->N(j, INSIDE));
        }
    }

  // 2. create new nodes and elements for a new mesh which coincides with these split boundaries
  size_t node_idx(mesh.Nodes());
  size_t elmt_idx(mesh.Elements());

  map<const Node<dim>* const,Node<dim>*>  node_pairs;

  for ( auto& ifnode : ifnodes )
    {
      // TODO: check interface to determine whether the nodes are indeed collocated
      Node<dim>* const new_node = mesh.Duplicate( ifnode );
      if (new_node != nullptr) {
          new_node->Idx(node_idx++);
          // storing pointers to the interface nodes as keys to retrieve the new nodes
          node_pairs.insert(make_pair(ifnode, new_node));
        }
    }

  const LocalVariables&             ifvars( splitboundaryComplex->Database().LocalVariablesAt( ELEMENT ) );
  const IntegrationPointVariables&  if_ip_vars( splitboundaryComplex->Database().IntegrationPointVariablesAt( ELEMENT_INTEGRATION_POINT ) );

  // creating the new lower-dimensional elements between the sides of collected InterFace objects and assigning them to middle-element pointer
  vector<uint32_t>  element_numbers;
  element_numbers.reserve( ifelmts.size() );

  for ( auto& ifelmt : ifelmts )
    {
      CSMP_FEM_TYPE csmp_elmt = ifelmt->FE_Type();
      Element<dim>* const new_elmt = mesh.AddElement( fem_mgr.E(csmp_elmt),
                                                      fvm_mgr->Stencil(csmp_elmt),
                                                      ifvars, if_ip_vars,
                                                      node_vec, neighbor_vec,
                                                      material_id );
      new_elmt->Idx(elmt_idx++);
      element_numbers.push_back( new_elmt->Idx() );
      // assign the new element to the intervening element pointer of the InterFace element in the Splitboundary
      ifelmt->Assign( new_elmt ); 
   }

  // 3. Construct the connections between the new nodes, new elements and interfaces.
  //    Assign a root node and a root element for new nodes elements respectively
  vector<Element<dim>*>  cell_vec_to_establish_nbor_connectivity;
  cell_vec_to_establish_nbor_connectivity.reserve( ifelmts.size() );
  
  for ( auto& ifelmt : ifelmts )
    {
      Element<dim>* const new_elmt = ifelmt->InterveningElement();

      for (auto i = 0; i < ifelmt->Nodes() / 2; i++) {
          // finding the new node using the matching node of the inside of the interface as search key
          Node<dim>* new_node = node_pairs[ ifelmt->N(i, INSIDE) ];
          assert( new_node != nullptr );
          new_elmt->Assign(i, new_node);
          new_node->ResizeParentStorage(new_node->Parents() + 1);
          new_node->Assign(new_node->Parents(), new_elmt);
        }

      cell_vec_to_establish_nbor_connectivity.push_back(new_elmt);
    }

  // update the neighbor connectivity of new region
  // TODO: potential manifolds have to be disambiguated 
  establishNeighborConnectivity( cell_vec_to_establish_nbor_connectivity, false, false ); 

  // for each contiguous patch of the new region supply a pointer any of its elements into mesh manager 
  sort( cell_vec_to_establish_nbor_connectivity.begin(), cell_vec_to_establish_nbor_connectivity.end() );
  set<Element<dim>*>    contiguous_subset;
  vector<Element<dim>*> leftovers;
  
  while ( !cell_vec_to_establish_nbor_connectivity.empty() ) 
    {
       // finding contiguous element patch and inserting its first element into the root element vector
       floodFill( (*cell_vec_to_establish_nbor_connectivity.begin()), contiguous_subset );
       // if the split boundary is already contiguous
       if ( contiguous_subset.size() == cell_vec_to_establish_nbor_connectivity.size() ) break;
       // removing the pointers to the recovered elements from 'cell_vec_to_establish_nbor_connectivity'
       leftovers.reserve( cell_vec_to_establish_nbor_connectivity.size() - contiguous_subset.size() );
       for ( typename vector<Element<dim>*>::const_iterator 
             it=cell_vec_to_establish_nbor_connectivity.begin(); it!=cell_vec_to_establish_nbor_connectivity.end(); ++it )
         // copy remaining elements into the leftover vector   
         if ( contiguous_subset.count( (*it) ) == 0 )
           leftovers.push_back( (*it) ); 

       // assigning result vector to repeat operation                
       cell_vec_to_establish_nbor_connectivity = leftovers;
       // contiguous_subset.clear(); - done in floodfill
       leftovers.clear();
    }

  // 4. constuct a new unique region with the interface elements in the model
  const bool unique_map(true);
  splitboundaryComplex->FormRegionFrom( region_name, element_numbers, unique_map );

  return true;

} // end SingleRegionFromAllSplitBoundaries

*/ // MAYBE BRING THIS BACK LATER in rewritten form





/**
    In all SplitBoundaries of the model, lower dimensional Element regions are inserted.

       As many lower-dimensional element regions are created as there SplitBoundaries in the model
       @author SKM
       @date 26/1/2020
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
set<string>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertLowerDimensionalRegionsIntoSplitBoundaries( int32_t material_id )
{
  ErrorHandler&                csmp_error(ErrorHandler::Instance());
  SPLITBOUNDARY_COMPLEX<dim>*  splitboundaryComplex(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
  set<string>                  new_regions;
  
  for ( typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator
        it = splitboundaryComplex->SplitBoundariesBegin(); it != splitboundaryComplex->SplitBoundariesEnd(); ++it )
    {
       pair<string,bool>  result = splitboundaryComplex->InsertRegionIntoSplitBoundary( (*it).first.c_str(), material_id );
       if ( result.second == false )
         csmp_error.notice( ERROR, "SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertLowerDimensionalRegionsIntoSplitBoundaries", 
                            (*it).first, "unable to create Region from this boundary" );
       else {
            new_regions.insert( result.first );
            cout <<"\nSplitBoundaryInterface<"<< dim <<",Model>::InsertLowerDimensionalRegionsIntoSplitBoundaries: created new region '";
            cout << result.first <<"' from SplitBoundary."<< endl;
         }
    }
  
  return new_regions;

} // end InsertLowerDimensionalRegionsIntoSplitBoundaries







/**
Prints current SplitBoundaries to screen
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesOut() const
  {
     const SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<const SPLITBOUNDARY_COMPLEX<dim>*>(this));

     cout <<"\nSplitBoundaryInterface<"<< dim <<",SplitBoundary<InterFace>>::SplitBoundariesOut: ";
     if ( SplitBoundaries() == 0 ) {
          cout <<"\tmodel does not contain any split boundaries.\n\n";
          return;
       }
     cout <<"split boundaries of ";
     if ( model->BoxShaped() ) cout <<"box-shaped model:\n";
     else cout <<"irregularly-shaped model:\n";
     for ( auto bit=SplitBoundariesBegin(); bit!=SplitBoundariesEnd(); ++bit ) {
          cout <<"\n\t"<< (*bit).first <<", box-flag: "<< parseBoundary( (*bit).second.AtBoundary() );
          cout <<" "<< (*bit).second.Cells() <<" interfaces, ";
          // in 3D a boudary is a surface
           if constexpr ( dim == 3 ) {
                cout <<"area (m2): "<< (*bit).second.Area();
                cout <<", perimeter length (m): "<< (*bit).second.Perimeter();
             }
           if constexpr ( dim == 2 )
             cout <<" length (m): "<< (*bit).second.Area();
       }
     cout << endl << endl;
     cout.flush();

} // end Out





template class SplitBoundaryInterface<1U, Model>;
template class SplitBoundaryInterface<2U, Model>;
template class SplitBoundaryInterface<3U, Model>;

} // csmp







// DEPRECATED - SKM 25/1/2020

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
 /*
 template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
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

   std::string  splitboundaryName( "SPLITBOUNDARY_" + region );
   splitboundaryComplex->InsertBoundary( IRREGULAR, region.c_str() );
   Boundary<dim>& boundary( splitboundaryComplex->Boundary( region.c_str() ) );

   splitboundaryComplex->RemoveRegion( region.c_str(), false );

   // attempt to create a regular (InterFace-based) splitboundary
   std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
     it = splitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName, splitboundaryComplex->Database() ) ) );
   if ( it.second )
   {
     //std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary around: " << region << std::endl;
     splitboundaryComplex->UpdateIndices();
     succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );

   }
   else
     throw csmp::Exception( INFO, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                            splitboundaryName.c_str(),
                            "boundary already exists. Nothing was done." );

   // remove temporal boundary
   splitboundaryComplex->RemoveBoundary( boundary );

   // update indexes
   splitboundaryComplex->UpdateIndices();

   UpdateSplitBoundaryComplex();

   if ( succeeded )
     std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary created splitboundary around " << region << std::endl;
   else
     throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

   return true;

 } // end InsertSplitBoundary

*/





/**
Method forms a SplitBoundary between the two supplied regions. This will involve the creation
and connection of InterFaces by the MeshManager.

@attention:  BoundariesInterface cannot be created in 1D models or between regions which contain
one-dimensional elements.

This method will only work for Regions which are not overlapping.

This method makes no sense for the region 'Model' as it encompasses all unique regions.
*/
/*
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const std::string& group1, const std::string& group2, bool createRegionBetween )
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( dim == 1 ) {
    csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "There are no SplitBoundaries in 1D models." );
    return false;
  }

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

  std::string  boundaryName( string(group1) + string("_") + string(group2) );
  pair<string, string> key = make_pair( group1, group2 );

  splitboundaryComplex->InsertBoundary( group1.c_str(), group2.c_str(), createRegionBetween );
  Boundary<dim>& boundary( splitboundaryComplex->Boundary( boundaryName ) );

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
    // TODO: make name consistent with name of boundary (replace string BOUNDARY with SPLIT_BOUNDARY)
    std::string  splitboundaryName( CreateSplitBoundaryName( key ) );
    boundaryName = splitboundaryName;
    std::pair<typename std::map<std::string, csmp::SplitBoundary<dim> >::iterator, bool>
      it = splitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,
                                     splitboundaryComplex->Database() ) ) );
    if ( it.second )
    {
      //std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary between " << group1 << " and " << group2 << std::endl;
      splitboundaryComplex->UpdateIndices();
      succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );
    }
    else
      throw csmp::Exception( INFO,
                             "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                             splitboundaryName.c_str(),
                             "boundary already exists. Nothing was done." );
  }

  splitboundaryComplex->UpdateIndices();

  // remove temporal boundary
  splitboundaryComplex->RemoveBoundary( boundary );

  // update indexes
  splitboundaryComplex->UpdateIndices();

  UpdateSplitBoundaryComplex();

  if ( succeeded ) {
       cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary: created split boundary '";
       cout << boundaryName <<"' between " << group1 << " and " << group2 << std::endl;
    }
  else
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

  return true;

} // end InsertSplitBoundary



DEPRECATED END */
