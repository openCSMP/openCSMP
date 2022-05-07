#include "SplitBoundaryInterface.h"
#include "ModelTopology.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Model.h"
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

#define SPLIT_BOUNDARY_DEBUG

using namespace std;

namespace csmp {

/**
returns reference to SplitBoundary
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundary<dim>&  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundary( const string& spbName )
{
  splitBoundaryIterator  sbit = splitBoundaryMap_.find( spbName );
  if ( sbit != splitBoundaryMap_.end() )
    return (*sbit).second;
  else {
    string errMsg( "SplitBoundary does not exist!" );
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
const SplitBoundary<dim>&  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundary( const string& spbName ) const
{
  splitBoundaryConstIterator  sbit = splitBoundaryMap_.find( spbName );
  if ( sbit != splitBoundaryMap_.end() )
    return (*sbit).second;
  else {
    string errMsg( "SplitBoundary does not exist!" );
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
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::ContainsSplitBoundary( const string& bname ) const
{
  if ( splitBoundaryMap_.find( bname ) != splitBoundaryMap_.end() ) return true;
  return false;
}


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename map<string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin()
{ return splitBoundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename map<string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd()
{ return splitBoundaryMap_.end(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename map<string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin() const
{ return splitBoundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
typename map<string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd() const
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
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( const char* splitboundary,
                                                                              bool erase_interfaces )
 {
     if ( !ContainsSplitBoundary(splitboundary) ) {
         ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface::RemoveSplitBoundary",
                                          splitboundary, "split boundary is not contained in model; nothing was done");
         return;
      }
    
    SPLITBOUNDARY_COMPLEX<dim>&  splitBoundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>& >(*this) );
    csmp::SplitBoundary<dim>&    split_boundary = splitBoundaryComplex.SplitBoundary( splitboundary );

     // erasing the faces
     if ( erase_interfaces ) {
         // getting the mesh manager to delete faces and nodes and fix up the connectivity
         splitBoundaryComplex.Mesh().DeleteAndRepairConnnectivity( split_boundary.CellVector().begin(), split_boundary.CellVector().end() );
       }

     // deleting the split boundary
     splitBoundaryMap_.erase( splitboundary );
     
     // TODO: should we remove duplicated nodes again?
     cout <<"\n"<<"SplitBoundaryInterface<"<< dim <<">::RemoveSplitBoundary: removed split boundary named '";
     cout << splitboundary <<"' albeing retaining duplicated nodes."<< endl;

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
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary,
                                                                              bool erase_interfaces )
 {
     // erasing the faces
     if ( erase_interfaces ) {
         SPLITBOUNDARY_COMPLEX<dim>&  splitBoundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>&>(*this) );
         // getting the mesh manager to delete faces and nodes and fix up the connectivity
         splitBoundaryComplex.Mesh().DeleteAndRepairConnnectivity( splitboundary.CellVector().begin(), splitboundary.CellVector().end() );
       }

     splitBoundaryMap_.erase( splitboundary.Name() );
    
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
  cout << bin_file << "'" << endl;
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
          csmp_error.Note( ERROR, "SplitBoundaryInterface::InputSplitBoundariesFromBinary:",
                             bin_file, "file could not be opened; nothing was done." );
          return false;
       }

  // 0. reading the file header and printing it to the screen
  {
     BinaryFileSectionRead sect(fp, "SBDFHEDR");

     char  text[INFO_STRING];
     binaryFileRead( fp, text );
     cout <<"\nSplitBoundaryInterface<"<< dim <<">::InputSplitBoundariesFromBinary: Reading file header:\n\t"<< text << endl;
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

  cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: reading " << csCache << " containing "
    << records << " splitboundaries\n" << endl;

  for ( size_t i( 0 ); i < records; ++i )
  {
    string bName;
    // reading name of splitboundary
    binaryFileRead( fp, csCache );
    bName = csCache;
    // inserting splitboundary if not existing yet
    pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
      bit = splitBoundaryMap_.insert( make_pair( bName, csmp::SplitBoundary<dim>( bName, splitboundaryComplex->Database() ) ) );
    // splitboundary
    if ( !bit.second )
      return false;
    if ( !bit.first->second.In( splitboundaryComplex->Mesh(), splitboundaryComplex->FE_Manager(), rref, fp ) )
      return false;
    cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: read boundary " << bName << " successfully.\n";
  }
  fp.close();

*/


/// container of juxtaposed element pairs for SplitBoundary creation:
template<uint32_t dim>
struct SplitBoundaryElementSets : public
  //       key consisting out the names of regions juxtaposed at the SplitBoundary
  map<pair<string, string>,
  // set that stores pointers to the pairs of elements juxtaposed across boundary
  // size_t parameter gives the local number of the element face that sits at the split boundary
  set<pair<pair<Element<dim>*, size_t>,
  pair<Element<dim>*, size_t> > > > {
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
string SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryName( const pair<string, string>& juxtaposed_regions ) const
{
  string split_boundary_name( "SPLITBOUNDARY_" + juxtaposed_regions.first + '_' + juxtaposed_regions.second );

  // if the substring set is empty
  if ( ContainsSplitBoundary( split_boundary_name.c_str() ) ) {
    ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateSplitBoundaryName:",
                                     split_boundary_name.c_str(), "already exists, try other name.'\0'." );
    return string( "\0" );
  }

  // making a set of boundary names
  /*
  const size_t substrings_used_in_search(intersected_regions.size());
  for ( auto it=SplitBoundariesBegin(); it!=SplitBoundariesEnd(); ++it ) {
  size_t substrings_found(0U);
  for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
  // if the substring is found
  if ( (*it).first.find(*ir) !=string::npos ) substrings_found++;
  // when all substrings are contained in the boundary name, it is returned
  if ( substrings_found == substrings_used_in_search )
  return (*it).first;
  }
  return string("\0");
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
pair<set<string>,bool> SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries()
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&               csmp_error( ErrorHandler::Instance() );
  
  // 0. Is the model contigous? - if so, the self-intersection problem needs to be addressed during boundary creation
  // ----------------------------------------------------------------------------------------------------------------
#ifdef DEBUG
  map<string,vector<Element<dim>*> >  mesh_patches;
  size_t n_mesh_patches = findStandAloneMeshPatches( splitboundaryComplex->Mesh().ElementsBegin(),
                                                     splitboundaryComplex->Mesh().ElementsEnd(),
                                                     mesh_patches );
  // checking the dimensionality of the patches
  size_t volume_patches{0U}, surface_patches{0U}, line_patches{0U}, equidimensional_patches{0};
  for ( const auto& it : mesh_patches ) {
       // only looking at the first element of each patch
       if ( it.second[0]->IsEquidimensional() ) equidimensional_patches++;
       if      ( it.second[0]->IsVolume() ) volume_patches++;
       else if ( it.second[0]->IsSurface() ) surface_patches++;
       else if ( it.second[0]->IsLine() ) line_patches++;
    }
  cout <<"\n"<<"SplitBoundaryInterface<"<< dim << ">::DetectAndCreateSplitBoundaries: model contains: "<< n_mesh_patches;
  cout <<" element patches. "<< equidimensional_patches <<" with same dimension as model."<< endl;
  cout << splitboundaryComplex->Name() <<" contains "<< volume_patches <<" volume, ";
  cout << surface_patches <<" surface, and "<< line_patches <<" line element patches."<< endl;
#endif

  // 0. creating region labels and tagging the regions with unique integer indentifiers
  // ----------------------------------------------------------------------------------
  const string region_tag("region identifier");
  if ( !splitboundaryComplex->Database().IsDefined(region_tag.c_str()) )
    splitboundaryComplex->CreateProperty( region_tag.c_str(), "X", SCALAR, ELEMENT );
  const csmp::Index reg_key = splitboundaryComplex->Database().StorageKey(region_tag.c_str());
  vector<string>  region_names;
  const size_t model_regions = splitboundaryComplex->CountAndLabelRegions( region_tag.c_str(), region_names );
  assert( model_regions > 1U );


  // 1. detecting potential SplitBoundaries
  // --------------------------------------
  // finding the (local) ids of the element faces on either side of the split boundary
  vector<pair<pair<Element<dim>*,uint32_t>, pair<Element<dim>*,uint32_t> > >  interface_elmt_pairs;
  if ( !findSplitInterfaceElements( *splitboundaryComplex, region_tag, interface_elmt_pairs ) ) {
      csmp_error.Note( WARNING, "SplitBoundaryInterface::DetectAndCreateSplitBoundaries:",
                         "node-coordinate matched faces / internal boundaries could not be detected; nothing was done." );
      return make_pair(set<string>(),false);
    }

  // 2. classifying the detected interfaces in terms of the regions that they juxtapose
  // ----------------------------------------------------------------------------------
  // 2.1 grouping interface element pairs into ones that juxtapose specific regions against one another
  //     these will later become specific split boundaries
  typedef vector<pair<pair<Element<dim>*, uint32_t>, pair<Element<dim>*, uint32_t> > > INTERFACE_ELEMENT_PAIRS;
  map<set<string>, INTERFACE_ELEMENT_PAIRS>  split_boundary_map;

  // for all the SplitBoundary objects supplied as sets of pairs of Element pointers and interface idx values
  for ( const auto& iit : interface_elmt_pairs ) {
      // extracting region names from the name-integer vector
      set<string> key{ region_names[static_cast<long>(iit.first.first->Read( reg_key ))],
                       region_names[static_cast<long>(iit.second.first->Read( reg_key ))] };
      // storing interfaces in split boundary maps
      auto eit = split_boundary_map.insert( make_pair( key, INTERFACE_ELEMENT_PAIRS( { iit } ) ) );
      // if no insertion could be performed, the element pair is added to an existing set
      if ( !eit.second ) (*eit.first).second.push_back( iit );
    }

  // echoing the map to the screen
#ifdef SPLIT_BOUNDARY_DEBUG
  cerr << "\nSplitBoundaryInterface::DetectAndCreateSplitBoundaries: interface region pairs found:\n";
  for ( const auto& split : split_boundary_map ) {
       for ( const auto& j : split.first ) cout << j <<",";
       cout << endl;
    }
  cerr << endl;
#endif

  // 3. Creating splitboundaries for each of the discovered juxtapositions of regions
  // --------------------------------------------------------------------------------
  pair<set<string>,bool>  splitBoundaryNames; 
  splitBoundaryNames.second = true;
  for ( auto& it : split_boundary_map ) {
      // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
      string bname = CreateSplitBoundaryName( make_pair( *(it.first.begin()),  *(next(it.first.begin(),1)) ) );
      splitBoundaryNames.first.insert( bname );
      
      // extracting the split boundary neighbor elements into interface element pairs
      InterFaceParentElements<dim>  ifset( it.second );
      
      // establishing the InterFace and NodeManifold objects and connecting them all
      vector<InterFace<dim>*>  iface_ptrs = splitboundaryComplex->Mesh().CreateInterfacesBetweenNodeMatchingElements(
                                                                                    splitboundaryComplex->Database(),
                                                                                    ifset );

      // creating the SplitBoundary asking the MeshManager to create the required number of InterFace objects
      pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
        bit = splitBoundaryMap_.insert( make_pair( bname, csmp::SplitBoundary<dim>( bname, splitboundaryComplex->Database() ) ) );
      if ( bit.second ) {
           (*bit.first).second.CreateFrom( iface_ptrs.begin(), iface_ptrs.end() );
           cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries: '"<< bname;
           cout << "' created successfully.\n";
           // NOTE: no rebuild of regions is needed here because no new nodes were created
           //       and none of the existing regions were changed!
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





/** Form consistently named split boundaries between all the unique regions in the current model, while detecting and honouring already existing
    Boundary and SplitBoundary objects identified by INTERNAL boundary flags. These lower-dimensional model subdomains will not be touched.
    Only those unique regions are considered which have the same spatial dimension as the model.
    
    @attention all perimeter nodes / contact points between split boundaries in 2D will become NodeManifolds connecting as many Node objects
    as regions are in contact with one another at that split boundary.
 */
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
size_t SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SeparateUniqueRegionsBySplitBoundaries()
  {
    SPLITBOUNDARY_COMPLEX<dim>* model( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
    //ErrorHandler&               csmp_error( ErrorHandler::Instance() );
    
    cout << "\nStart creating splitBoundaries between the unique regions (and respecting the existing boundaries:\n"<<endl;
    // creating region labels and tagging the regions with unique integer indentifiers
    const string region_tag("region identifier");
    if ( !model->Database().IsDefined(region_tag.c_str()) )
      model->CreateProperty( region_tag.c_str(), "X", SCALAR, ELEMENT );
    // TODO: use ElementMaterial_ID here rather than relying on a new variable
    const csmp::Index mtrl_key = model->Database().StorageKey(region_tag.c_str());
    vector<string>  region_names;
    const size_t model_regions = model->CountAndLabelRegions( region_tag.c_str(), region_names );

    // search the existing unique sub-regions
    set<pair<string, string>>     discovered;
    deque<string>                 current_regions;
    vector<pair<string, string>>  region_final_pairs;
    
    for ( const auto& region_name : region_names )
      {
         // starting at the first region
         current_regions.push_back( region_name );
        
         while ( !current_regions.empty() )
           {
             string current_region( *current_regions.begin() );
             set<string> neighbors;
             // for all neighbor sub-regions of the current region
             const csmp::Region<dim>&  region1( model->Region( current_region ) );
             for ( auto j{0U}; j < region_names.size(); j++ ) {
                  if ( current_region == region_names[j] ) continue;
                  const csmp::Region<dim>&  region2( model->Region( region_names[j] ) );
                  bool same_material(false);
                  if ( static_cast<int32_t>(region1.E(0)->Read(mtrl_key)) ==
                       static_cast<int32_t>(region2.E(0)->Read(mtrl_key)) ) same_material = true;
                  // only perimeter nodes may be shared in the case of unique (non-overlapping) regions
                  const size_t shared_nodes( sharedPerimeterNodes( region1, region2 ) );
                  if ( shared_nodes > 1 && !same_material)
                    neighbors.insert( region_names[j] );
               }
             for ( auto neighbour_region : neighbors ) {
                  // if this neighbor is new one
                  auto new_region = discovered.insert( make_pair( current_region, neighbour_region ) );
                  if ( new_region.second ) {
                      discovered.insert( make_pair( neighbour_region, current_region ) );
                      current_regions.push_back( neighbour_region );
                      region_final_pairs.push_back( make_pair( current_region, neighbour_region ) );
                    }
               }
            // removing the sub-region from the discovered (but not yet explored) deque
            current_regions.pop_front();
          }
      }

    // CREATION OF SPLITBOUNDARY BETWEEN 2 REGIONS
    // -------------------------------------------
    // do this sequentially according to neighbours, otherwise boundaries are not assgiend properly
    for ( auto it : region_final_pairs ) // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
      model->CreateSplitBoundaryBetween( it.first.c_str(), it.second.c_str() );

    // screen output
    // cout <<"\n\n"<<"SplitBoundaryInterface::SeparateUniqueRegionsBySplitBoundaries: split boundaries after separation.\n";
    // model->SplitBoundariesOut();
  
    return SplitBoundaries();
  
} // end SeparateUniqueRegionsBySplitBoundaries
    
/* there are no duplicates
    // removing potential duplicates
    for ( auto& pit : region_final_pairs ) if ( pit.first > pit.second ) swap( pit.first, pit.second );
    region_final_pairs.erase( unique( region_final_pairs.begin(), region_final_pairs.end() ), region_final_pairs.end() );
*/






template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
size_t SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::FormSplitBoundariesFrom( const ModelTopology& topo )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. getting the names of the regions
  list<string> split_boundaries;
  topo.OutputSplitBoundaries( split_boundaries );

  // 2. assigning the regions to groups in the Model
  cout << "\nSplitBoundaryInterface::FormSplitBoundariesFrom: Forming the split boundaries: ";

  size_t new_split_boundaries{0};
  for ( auto lit = split_boundaries.begin(); lit != split_boundaries.end(); lit++ )
    {
      string domain_name( *lit );
      auto it = splitBoundaryMap_.insert( make_pair( domain_name, csmp::SplitBoundary<dim>( domain_name, static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this)->Database() ) ) );
      // if the region was successfully inserted
      if ( it.second )
        {
          // making a list of the element numbers
          vector<size_t>  cell_ids;
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
              csmp_error.Note( WARNING, "SplitBoundaryInterface::FormSplitBoundariesFrom",
                                 "SplitBoundary could not be formed", (*lit).c_str() );
            }
          else {
               // reporting the name of the newly generated region
               cout << domain_name << " ";
               new_split_boundaries++;
            }
        }
      else
        throw csmp::Exception( ERROR, "SplitBoundaryInterface::FormSplitBoundariesFrom",
                              "SplitBoundary could not be formed. Does this region already exist?", (*lit).c_str() );
    }
  cout << endl;

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
    //    CreateInternalBoundaryFrom checks whether dim_1_region actually exists
    pair<set<string>,bool> boundary_names = model->CreateInternalBoundaryFrom( dim_1_region );
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
     
     @attention the input boundary is removed in the process and its Face objects will be erased.
     
     @return boolean indicating whether the method was able to create a singe Split boundary and its name
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
pair<string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom( Boundary<dim>& boundary )
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // checking whether boundary is external to the model in which case a SplitBoundary cannot be buid
  if ( boundary.IsExternal() ) {
       csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom", "Region 'Model' not eligible for InsertSplitBoundary." );
       return make_pair( "no SplitBoundary was created", false );
    }

  // creating boundary name by replacing BOUNDARY with SPLIT_BOUNDARY
  string   splitboundaryName( boundary.Name() );
  const size_t  str_length(string("BOUNDARY").length());
  splitboundaryName.replace( splitboundaryName.find("BOUNDARY"), str_length, "SPLIT_BOUNDARY" );

  // attempt to create a splitboundary
  if ( ContainsSplitBoundary(splitboundaryName) ) {
       csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom", splitboundaryName,
                        "a SplitBoundary with this name already exists; nothing was done." );
       return make_pair( "no SplitBoundary was created", false );
    }
  
  pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
    it = splitBoundaryMap_.insert( make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,
                                                                                           splitboundaryComplex->Database() ) ) );
  if ( it.second ) {
      cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom:";
      cout <<" creating SplitBoundary from 'Boundary' "<< boundary.Name() << endl;
      (*it.first).second.CreateFrom( splitboundaryComplex->Database(), splitboundaryComplex->Mesh(), boundary );
    }

  // removes boundary, noting that its interface objects were already deleted by CreateFrom()
  const bool erase_faces{ false };
  splitboundaryComplex->RemoveBoundary( boundary, erase_faces );

  cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom: created splitboundary: '";
  cout << splitboundaryName <<"' successfully.\n\n";

  return make_pair( splitboundaryName, false );

} // end CreateSplitBoundaryFrom( Boundary )









/**

Creates SplitBoundary between two unique regions that share nodes along their perimeter.
The first region will be placed on the inside of the new split boundary.

@param region1_name name of the unique region that shall be on the inside of the new split boundary
@param region2_name of the unique region that touches (shares nodes with)  region1 and shall become the outside region of the new split boundary

@note Method will add new Node and InterFace objects to the model

@note method does not create a temporary Boundary on the way to the SplitBoundary construction

@note method ignores Elements that are lines

@author SKM
@date 25/1/20

 */
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
pair<string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryBetween( const char* region1_name, const char* region2_name )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    SPLITBOUNDARY_COMPLEX<dim>*  modelComplex(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
  
    // 1. initial diagnostics: verifying inputs and shared faces between the two regions
    // ---------------------------------------------------------------------------------
    if ( string{region1_name} == region2_name ) {
         csmp_error.Note( ERROR, "BoundaryInterface::CreateSplitBoundaryBetween:", "Provided Regions are the same.");
         return make_pair("split boundary not created",false);
      }
    if ( !modelComplex->IsUnique(region1_name) || !modelComplex->IsUnique(region2_name) ) {
         csmp_error.Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryBetween:", "This method only works for unique Region objects");
         return make_pair("split boundary not created",false);
      }
   
    const csmp::Region<dim>&  region1(modelComplex->Region(region1_name));
    csmp::Region<dim>&        region2(modelComplex->Region(region2_name)); // not const because must be rebuilt
    
    // 2. Is there a shared interface? - if so InterFace objects are created along it
    // ------------------------------------------------------------------------------
    vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > > matched_elmts;
    sharedPerimeterCells( region1, region2, matched_elmts );
    
    if ( matched_elmts.empty() ) {
         string message( string(" input regions '") + region1_name + "' and '" + region2_name +"'");
         csmp_error.Note( WARNING, "SplitBoundaryInterface::CreateSplitBoundaryBetween:",
                          message, "may share a node but do not share any faces");
         return make_pair("split boundary not created",false);
      }
        
    // creating the necessary interfaces and nodes, and updates the connectivity of the mesh
    const bool create_manifolds_on_perimeter{true};
    vector<InterFace<dim>*>  interfaces = modelComplex->Mesh().CreateInterfacesBetweenNodeSharingElements( modelComplex->Database(),
                                                                                                           matched_elmts,
                                                                                                           create_manifolds_on_perimeter );
    // 3. creation of the new SplitBoundary
    // ------------------------------------
    string split_boundary_name = CreateSplitBoundaryName( make_pair( region1_name, region2_name ) );
    pair<typename map<string,csmp::SplitBoundary<dim> >::iterator,bool>
      it = splitBoundaryMap_.insert( make_pair( split_boundary_name, csmp::SplitBoundary<dim>( split_boundary_name, modelComplex->Database() ) ) );

    if ( it.second ) {
        cout << "\nSplitBoundaryInterface<"<< dim <<">::CreateSplitBoundaryBetween: creating split boundary between '";
        cout << region1_name << "' and '" << region2_name <<"'"<< endl;
        bool succeeded = (*it.first).second.CreateFrom( interfaces.begin(), interfaces.end() );
 
        if ( !succeeded )
          csmp_error.Note( ERROR, "SplitBoundaryInterFace::CreateSplitBoundaryBetween:", "creation failed");
        else {
             cout << "\nSplitBoundaryInterface<"<< dim <<">::CreateSplitBoundaryBetween: split boundary '";
             cout << split_boundary_name << "' created successfully."<< endl;
             // flagging outside region for rebuild because it has new nodes
             region2.ScheduleForRebuilt();
             modelComplex->UpdateRegions();
             // all done
             return make_pair(split_boundary_name,true);
          }
      }
    else throw csmp::Exception( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryBetween:",
                                split_boundary_name, "split boundary already exists. Nothing was done.");

    return make_pair("split boundary not created",false);

 } // end CreateSplitBoundaryBetween





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
pair<string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertRegionIntoSplitBoundary( const char* split_boundary,
                                                                                                      int32_t material_id )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( !ContainsSplitBoundary( split_boundary ) ) {
         csmp_error.Note( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertRegionIntoSplitBoundary",
                            split_boundary, "Does not exist; nothing was done." );
         return make_pair("no Region created",false);
       }
     
     SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
     csmp::SplitBoundary<dim>&    splitBoundary( model->SplitBoundary(split_boundary) );
     csmp::MeshManager<dim>&      mesh( model->Mesh() );
     const LocalVariables&        lvars( model->Database().LocalVariablesAt( NODE ) );
     
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
                           node_pointers.push_back( mesh.Duplicate( (*it)->N(i), MIDDLE, lvars ) );
                      }
                    //                            inserts the new node it creates into the corresponding NodeManifold
                    else node_pointers.push_back( mesh.Duplicate( (*it)->N(i), MIDDLE, lvars ) );
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
          const auto n_nodes( (*it)->FE()->Nodes() );
          vector<Node<dim>*> nodes( &node_pointers[node_offset], &node_pointers[node_offset+n_nodes] );

          // construct the new element
          new_elmts.push_back(  model->Mesh().AddInterveningElement( (*it),
                                                                     element_props, integration_point_props,
                                                                     nodes, material_id ) );
          node_offset += n_nodes;
       }
       
     // 3. establishing neighbor connectivity among the new elements
     // ------------------------------------------------------------
     if constexpr ( dim == 2U )
       mesh.template BuildLineConnectivity<Element>( new_elmts.begin(), new_elmts.end() );
     if constexpr ( dim == 3U )
       mesh.template BuildSurfaceConnectivity<Element>( new_elmts.begin(), new_elmts.end() );


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
  
  for ( typename map<string, csmp::SplitBoundary<dim> >::const_iterator
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
  
  for ( auto it = splitboundaryComplex->SplitBoundariesBegin(); it != splitboundaryComplex->SplitBoundariesEnd(); ++it )
    {
       pair<string,bool>  result = splitboundaryComplex->InsertRegionIntoSplitBoundary( (*it).first.c_str(), material_id );
       if ( result.second == false )
         csmp_error.Note( ERROR, "SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertLowerDimensionalRegionsIntoSplitBoundaries", 
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
     if ( SplitBoundaries() == 0U ) {
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
           if constexpr ( dim == 3U ) {
                cout <<"area (m2): "<< (*bit).second.Area();
                cout <<", perimeter length (m): "<< (*bit).second.Perimeter();
             }
           if constexpr ( dim == 2U )
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
     csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                        "Region 'Model' not eligible for InsertSplitBoundary." );
     return false;
   }

   if ( !splitboundaryComplex->IsUnique( region.c_str() ) )
     csmp_error.Note( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                        "This method is intended for the creation of splitboundaries around unique Regions" );

   string  splitboundaryName( "SPLITBOUNDARY_" + region );
   splitboundaryComplex->InsertBoundary( IRREGULAR, region.c_str() );
   Boundary<dim>& boundary( splitboundaryComplex->Boundary( region.c_str() ) );

   splitboundaryComplex->RemoveRegion( region.c_str(), false );

   // attempt to create a regular (InterFace-based) splitboundary
   pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
     it = splitBoundaryMap_.insert( make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName, splitboundaryComplex->Database() ) ) );
   if ( it.second )
   {
     //cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary around: " << region << endl;
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
     cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary created splitboundary around " << region << endl;
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
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const string& group1, const string& group2, bool createRegionBetween )
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( dim == 1 ) {
    csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "There are no SplitBoundaries in 1D models." );
    return false;
  }

  if ( group1 == "Model" or group2 == "Model" ) {
    csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Region 'Model' not eligible for InsertSplitBoundary." );
    return false;
  }
  if ( group1 == group2 ) {
    csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Provided Regions are identical." );
    return false;
  }


  if ( !splitboundaryComplex->IsUnique( group1.c_str() ) or !splitboundaryComplex->IsUnique( group2.c_str() ) ) {
    csmp_error.Note( WARNING,
                       "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                       "This method is intended for the creation of splitboundaries between unique Regions" );

    // checking for a potential overlap of the regions, if the regions are non-unique
    if ( splitboundaryComplex->RegionIntersection( group1.c_str(), group2.c_str(), "groupintersection" ) ) {
      csmp_error.Note( ERROR,
                         "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                         "one of the supplied regions is not unique and they overlap",
                         "It was therefore impossible to insert a boundary" );
      splitboundaryComplex->RemoveRegion( "groupintersection", false );
      return false;
    }
  }

  string  boundaryName( string(group1) + string("_") + string(group2) );
  pair<string, string> key = make_pair( group1, group2 );

  splitboundaryComplex->InsertBoundary( group1.c_str(), group2.c_str(), createRegionBetween );
  Boundary<dim>& boundary( splitboundaryComplex->Boundary( boundaryName ) );

  const csmp::Region<dim>&  gref1( splitboundaryComplex->Region( group1.c_str() ) );
  const csmp::Region<dim>&  gref2( splitboundaryComplex->Region( group2.c_str() ) );

  // checking whether the two regions share some nodes (these will mark their common boundary)
  const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
  bool succeeded( false );

  if ( shared_nodes == 0U )
    csmp_error.Note( WARNING,
                       "Model<dim,SPLITBOUNDARY_COMPLEX>::InsertBoundary",
                       "The regions of interest do not share any nodes; trying to create a split boundary" );


  // attempt to create a regular (InterFace-based) splitboundary
  else {
    // TODO: make name consistent with name of boundary (replace string BOUNDARY with SPLIT_BOUNDARY)
    string  splitboundaryName( CreateSplitBoundaryName( key ) );
    boundaryName = splitboundaryName;
    pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
      it = splitBoundaryMap_.insert( make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,
                                     splitboundaryComplex->Database() ) ) );
    if ( it.second )
    {
      //cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary between " << group1 << " and " << group2 << endl;
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
       cout << boundaryName <<"' between " << group1 << " and " << group2 << endl;
    }
  else
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

  return true;

} // end InsertSplitBoundary



DEPRECATED END */
