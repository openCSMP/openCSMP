// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SplitBoundaryInterface.h"

#include <ranges>

#include "ModelTopology.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "meshManagementUtilities.h"
#include "MeshManager.h"
#include "FaceConstructionData.h"
#include "Node.h"
#include "binaryReadWrite.h"

#include "Exception.h"
#include "ErrorHandler.h"

#define CSMP_SPLIT_BOUNDARY_DEBUG

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
                           errMsg );
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
                           errMsg );
  }
}


/**
Reports whether a SplitBoundary with the corresponding name exists inside the model.
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::ContainsSplitBoundary( const string& bname ) const
{
  if ( splitBoundaryMap_.contains(bname) ) return true;
  return false;
}


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
 map<string, typename csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin()
{ return splitBoundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
 map<string, csmp::SplitBoundary<dim> >::iterator
SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd()
{ return splitBoundaryMap_.end(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
 map<string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin() const
{ return splitBoundaryMap_.begin(); }


template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
 map<string, csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd() const
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

     // erasing the Interfaces
     if ( erase_interfaces ) {
         // finding the region(s) on the outside of the split boundary which will need to be updated after split boundary removal
         set<int32_t> domain_indices_of_outside_regions;
         for ( const auto& it : split_boundary.CellVector() )
           domain_indices_of_outside_regions.insert( it->OuterParent()->Region_ID() );
         // getting MeshManager to delete interfaces and nodes and fix up the connectivity
         splitBoundaryComplex.Mesh().DeleteInterfacesAndRepairConnnectivity( split_boundary.CellVector().begin(),
                                                                             split_boundary.CellVector().end() );
         // updating outside region(s) and 'Model'
         assert( !domain_indices_of_outside_regions.empty() );
         for ( auto& it : domain_indices_of_outside_regions ) {
              Region<dim>& region = splitBoundaryComplex.RegionByDomainIndex( it );
              region.RebuildSubDomainAfterChangeOfNodeVector();
           }
         splitBoundaryComplex.Region("Model").RebuildSubDomainAfterChangeOfNodeVector();
       }

     // deleting the split boundary
     splitBoundaryMap_.erase( splitboundary );
     
     cout <<"\n"<<"SplitBoundaryInterface<"<< dim <<">::RemoveSplitBoundary: removed split boundary named '"<< splitboundary <<"'"<< endl;

  } // end RemoveSplitBoundary

// testing
//splitBoundaryComplex.PrintDomainIndices();





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
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( csmp::SplitBoundary<dim>& split_boundary,
                                                                              bool erase_interfaces )
 {
     // erasing the faces
     if ( erase_interfaces ) {
         SPLITBOUNDARY_COMPLEX<dim>&  splitBoundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>&>(*this) );
         // finding the region(s) on the outside of the domain which will need to be updated after split boundary removal
         set<int32_t> domain_indices_of_outside_regions;
         for ( const auto& it : split_boundary.CellVector() )
           domain_indices_of_outside_regions.insert( it->OuterParent()->Region_ID() );
         // getting MeshManager to delete interfaces and nodes and fix up the connectivity
         splitBoundaryComplex.Mesh().DeleteInterfacesAndRepairConnnectivity( split_boundary.CellVector().begin(),
                                                                             split_boundary.CellVector().end() );
         // updating outside region(s)
         assert( !domain_indices_of_outside_regions.empty() );
         for ( auto& it : domain_indices_of_outside_regions ) {
              Region<dim>& region = splitBoundaryComplex.RegionByDomainIndex( it );
              region.RebuildSubDomainAfterChangeOfNodeVector();
           }
       }

     splitBoundaryMap_.erase( split_boundary.Name() );
    
  } // end RemoveSplitBoundary


 

/**
   Merges splitboundary objects and gives resultant Splitboundary a new name, original Splitboundaries are removed.
   Original interfaces are not deleted.
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
string SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::MergeSplitBoundaries( const char* new_sb_name, const set<string>& splitboundaries )
 {
      vector<InterFace<dim>*> iface_ptrs;
      for ( auto& sb_name : splitboundaries){
        //get splitboundary
        csmp::SplitBoundary<dim>& sb = this->SplitBoundary( sb_name );
        iface_ptrs.reserve( iface_ptrs.size() + sb.Cells() );
        for ( auto& ifit : sb.CellVector() )
          iface_ptrs.push_back(ifit);
      } //end of iface collection

      //Create merged splitboundary
      string new_name = string(new_sb_name) + "_SPLITBOUNDARY";
      AddSplitBoundary( new_name.c_str(), iface_ptrs.begin(), iface_ptrs.end() );

      //Delete original splitboundaries (but not interfaces) after merged splitboundary has been created
      for (auto& sb_name : splitboundaries){
        //Remove original splitboundary
        bool delete_interfaces = false;
        RemoveSplitBoundary( sb_name.c_str() , delete_interfaces ); //We do not want to delete interfaces, just moving them to a new splitboundary
      }

      cout << "SplitBoundaryInterface::MergeSplitBoundary() -> Successfully merged splitboundaries to " << new_name << endl;

      return new_name;

  } // end MergeSplitBoundaries




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
  
     const size_t records = this->SplitBoundaries();
     fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t) );

     for ( auto bit{ SplitBoundariesBegin() }; bit != SplitBoundariesEnd(); ++bit )
       {
          BinaryFileSectionWrite hdr(fp, "ONE_BDRY");
          cout <<"'"<< (*bit).first <<"' ";
          cout.flush();
          (*bit).second.WriteSplitBoundaryIndexesToBinaryFile( fp );
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
   
    size_t records(0);  // region records
    // getting number of unique region records from file
    fp.read( reinterpret_cast<char*>(&records), sizeof(size_t) );
    if ( records > 0 )
      // reading the regions sequentially
      for ( size_t i{0U}; i<records; ++i )
        {
           BinaryFileSectionRead hdr(fp, "ONE_BDRY");

           // 1.1 reading name and face indices for each boundaries
           SubDomainInfo  info;
           readSplitBoundaryIndexesFromBinaryFile( fp, info );
          
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
  if ( ContainsSplitBoundary( split_boundary_name ) ) {
    ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::CreateSplitBoundaryName:",
                                     split_boundary_name, "already exists, try other name" );
    return "\0";
  }

  return split_boundary_name;

} // end CreateSplitBoundaryName

  // TESTING making a set of boundary names
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





/**
     protected METHOD (not known beyond this compilation unit)
 
     creates underscore-separated unique names for the region patches based on the juxtapositions relationships
     across the lower-dimensional regions, the names are composed of:
 
     1. the name of the master region
     
     2. "SPLITBOUNDARY"
     
     3. the patch identifier number attached to boundary
     
     4. the name of the inner region, i.e. the region that the lower-dimensional element normals point away from
     
     5. the name of the outer region, i.e. that into which the normals point
     
     @attention  where the boundary just intersects a single layer (same material on either side), the layer name appears
     only once. The second instance is replaced by INTERSECTION.
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
string SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryNameFrom( const FaceConstructionData<dim>& fdata,
                                                                                       const vector<string>& region_names ) const
 {
     const long n_region_names = static_cast<long>(region_names.size());
     assert( fdata.ElementMaterial() < n_region_names );
     string boundary_name( region_names[ static_cast<size_t>(fdata.ElementMaterial()) ] );
     boundary_name += "_SPLITBOUNDARY";
     boundary_name += to_string(fdata.PatchNumber());
     boundary_name += '_';
     pair<long,long> materials(fdata.Materials()); // material identifiers can be negative
     assert( materials.first  < n_region_names );
     assert( materials.second < n_region_names );
     boundary_name += region_names[ static_cast<size_t>(materials.first) ];
     boundary_name += '_';
     if ( materials.first == materials.second ) boundary_name +="INTERSECTION";
     else boundary_name += region_names[ static_cast<size_t>(materials.second) ];
   
     return boundary_name;
 }







/**
    Searches the model for splitboundaries the name of which contains the search strings
    provided via the first set. The results are returnd into the second set.
    
    @note Use this method, for example, to retrieve multiple splitboundary patches that were generated
    from a single lower-dimensional regon, like a fault surface.

     @author SKM
     @date March 2016
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
size_t SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::FindSplitBoundaryByNames( const set<string>& intersected_regions,
                                                                                    set<string>& region_patches_found ) const
 {
    const SPLITBOUNDARY_COMPLEX<dim>& splitBoundaryComplex( static_cast<const SPLITBOUNDARY_COMPLEX<dim>& >(*this) );
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::FindSplitBoundaryByNames:",
                                         "supplied set of substrings is empty; returning '\0'." );
         return 0U;
      }
    // if the model has no boundaries
    if ( splitBoundaryComplex.Boundaries() == 0 ) {
         ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::FindSplitBoundaryByNames:",
                                         "model has no boundaries; returning '\0'." );
         return 0U;
      }
    region_patches_found.clear();
   
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=splitBoundaryComplex.SplitBoundariesBegin(); it!=splitBoundaryComplex.SplitBoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for (const auto& intersected_region : intersected_regions)
            // if the substring is found
            if ( (*it).first.find(intersected_region) !=string::npos ) substrings_found++;
          // when all substrings are contained in the boundary name, it is returned
          if ( substrings_found == substrings_used_in_search )
            region_patches_found.insert( (*it).first );
       }
    // return how many region patches contain the search string(s)
    return region_patches_found.size();
    
 } // end FindSplitBoundaryByNames






/**
  Expects a unique range of valid interfaces that were aready connected with one another (=neighbors assigned)
 */
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::AddSplitBoundary( const char* split_boundary_name,
                                                                          typename vector<InterFace<dim>*>::iterator ifacesBegin,
                                                                          typename vector<InterFace<dim>*>::iterator ifacesEnd )
 {
    SPLITBOUNDARY_COMPLEX<dim>* const splitBoundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>* const>(this) );
    assert( splitBoundaryComplex != nullptr );
    assert( ifacesBegin != ifacesEnd );
    
    // inserting boundary if it does not existing yet
    auto it = splitBoundaryMap_.insert( make_pair( split_boundary_name, csmp::SplitBoundary<dim>( split_boundary_name,
                                                                                                  splitBoundaryComplex->Database() ) ) );
    if ( it.second ) {
         // at this point, the interfaces have already been connected with each-other (SKM verified)
         (*it.first).second.CreateFrom( ifacesBegin, ifacesEnd );
         cout << "\nSplitBoundaryInterface<"<< dim <<">::AddSplitBoundary: successfully created split boundary '";
         cout << split_boundary_name <<"' from input faces.";
         (*it.first).second.Name(split_boundary_name);
      }
    else {
         ErrorHandler&  csmp_error( ErrorHandler::Instance() );
         csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::AddSplitBoundary:",
                            split_boundary_name, "SplitBoundary already exists or other problem arose. Nothing was done.");
         return false;
      }

     return true;
   
 } // end AddSplitBoundary






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
  size_t n_mesh_patches = findContiguousMeshPatches( splitboundaryComplex->Mesh().ElementsBegin(),
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
    splitboundaryComplex->CreateProperty( region_tag.c_str(), "rid", "uint", SCALAR, ELEMENT );
  const csmp::Index reg_key = splitboundaryComplex->Database().StorageKey(region_tag.c_str());
  vector<string>  region_names;
  const size_t model_regions = splitboundaryComplex->CountAndLabelUniqueRegions( region_tag.c_str(), region_names );
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
      set<string> key{ region_names[static_cast<size_t>(iit.first.first->Read( reg_key ))],
                       region_names[static_cast<size_t>(iit.second.first->Read( reg_key ))] };
      // storing interfaces in split boundary maps
      auto eit = split_boundary_map.insert( make_pair( key, INTERFACE_ELEMENT_PAIRS( { iit } ) ) );
      // if no insertion could be performed, the element pair is added to an existing set
      if ( !eit.second ) (*eit.first).second.push_back( iit );
    }

  // echoing the map to the screen
#if defined(DEBUG) && defined(CSMP_SPLIT_BOUNDARY_DEBUG)
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
                               bname,
                               "boundary already exists. Nothing was done." );
    }

  // update region on the outside of the new SplitBoundary because it needs a new node vector
  for ( auto rit=splitboundaryComplex->UniqueRegionsBegin(); rit!=splitboundaryComplex->UniqueRegionsEnd(); ++rit )
   (*rit).second.ScheduleForRebuild();

  //update outside regions
  // RebuildSubDomainAfterChangeOfCellVector() is not sufficient because the model needs to be rebuild as well
  splitboundaryComplex->UpdateRegions();

  // reporting the names of the split boundaries that were created
  return splitBoundaryNames;

} // end DetectAndCreateSplitBoundaries





/**
    Form consistently named split boundaries between all the unique regions of the current model,  detecting and honouring already existing
    Boundary and SplitBoundary objects identified by INTERNAL boundary flags. These lower-dimensional model subdomains will not be touched.
    Only those unique regions are considered which have the same spatial dimension as the model.
    
    @attention all perimeter nodes / contact points between split boundaries in 2D will become NodeManifolds connecting as many Node objects
    as regions are in contact with one another at that split boundary.
 */
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
size_t SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SeparateUniqueRegionsBySplitBoundaries()
  {
    SPLITBOUNDARY_COMPLEX<dim>* model( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
    //ErrorHandler&             csmp_error( ErrorHandler::Instance() );
    
    cout << "\nStart creating splitBoundaries between the unique regions (and respecting the existing boundaries:\n"<<endl;
    // creating region labels and tagging the regions with unique integer indentifiers
    const string    region_tag("region identifier");
    vector<string>  region_names;
    if ( !model->Database().IsDefined(region_tag.c_str()) ) {
          model->CreateProperty( region_tag.c_str(), "rid", "uint", SCALAR, ELEMENT );
         model->CountAndLabelUniqueRegions( region_tag.c_str(), region_names );
      }
    else { // assigning region names
         region_names.reserve( static_cast<size_t>(distance( model->UniqueRegionsBegin(),model->UniqueRegionsEnd())) );
         for ( auto rit=model->UniqueRegionsBegin(); rit!=model->UniqueRegionsEnd(); ++rit )
           region_names.push_back( (*rit).first );
      }
    const csmp::Index mtrl_key = model->Database().StorageKey(region_tag.c_str());

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
             for ( size_t j{0U}; j < region_names.size(); j++ ) {
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
                      region_final_pairs.emplace_back( current_region, neighbour_region );
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
    // the new interfaces are already connected with one another

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
  if ( !split_boundaries.empty() ) cout << "\nSplitBoundaryInterface::FormSplitBoundariesFrom: Forming the split boundaries: ";
  else return 0;

  size_t new_split_boundaries{0};
  for (auto& split_boundarie : split_boundaries)
    {
      string domain_name( split_boundarie );
      auto it = splitBoundaryMap_.insert( make_pair( domain_name, csmp::SplitBoundary<dim>( domain_name, static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this)->Database() ) ) );
      // if the region was successfully inserted
      if ( it.second )
        {
          // making a list of the element numbers
          vector<size_t>  cell_ids;
          cell_ids.reserve( topo.CellsWithinDomain( split_boundarie.c_str() ) );
          copy( topo.CellsOfDomainBegin( split_boundarie.c_str() ),
                topo.CellsOfDomainEnd( split_boundarie.c_str() ),
                back_inserter( cell_ids ) );

          // retrieving the interfaces by their IDs and assigning them  to the region
          (*it.first).second.AccumulateByNumber( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this)->Mesh(), cell_ids );
          cell_ids.erase( cell_ids.begin(), cell_ids.end() );

          // removing the split boundary if it contains no interfaces
          if ( (*it.first).second.Cells() == 0U ) {
              splitBoundaryMap_.erase( it.first );
              csmp_error.Note( WARNING, "SplitBoundaryInterface::FormSplitBoundariesFrom",
                                 "SplitBoundary could not be formed", split_boundarie );
            }
          else {
               // reporting the name of the newly generated region
               cout << domain_name << " ";
               new_split_boundaries++;
            }
        }
      else
        throw csmp::Exception( ERROR, "SplitBoundaryInterface::FormSplitBoundariesFrom",
                              "SplitBoundary could not be formed. Does this region already exist?", split_boundarie );
    }
  cout << endl;

  return new_split_boundaries;

} // end FormSplitBoundariesFrom





/**
        Eddy Pezulli's method for postprocessing SplitBoundary intersections.
        
        In essence, his code performs a topological disambiguation.
        When a new SplitBoundary (a 2D surface in 3D, or line in  1D in 2D) intersects an existing one,
        the nodes at the intersection are already "split" into a NodeManifold.
        This code's job is to look at each lower-dimensional element in the new boundary and ensure it points to the specific node in the manifold
        that "belongs" to the same mesh patch.
        
        To achieve this goal, the method goes over all manifolds in the model,
        checks them, adding Nodes / NodeManifolds, trying to rectify them in the Face of the newly emerging Topology.

        @attention ASSUMPTION: each split boundaries terminates at the intersection with another.
    
        Algorithm: If we detect a FE-node manifold along the surface mesh that shall become the splitboundary, the surface
        is intersecting a pre-existing splitboundary.
       
        If so, we need to determine whether the lower dimensional element is attached on the inside or outside of this interface.
        Once determined, the lower dimensional element can be replaced by an InterFace element with correct OUTSIDE or INSIDE node.
        Now splitting can proceed like before (the lack of neighbor connectivity of the higher dimensional elements on either side of
        the InterFace (contributing to the SplitBoundary) will ensure that outside neighbor searches will not reach INSIDE elements
        and parent element searches for the nodes on the inside will not find elements on the outside.
       
        The algorithmic approach taken is not to view the intersection as two interfaces which cut through a splitboundary,
        but to view an intersection on the scale of a single interface which touches the node of a pre-existing interface.
        This is important because it allows T intersections and X intersections to be handled with the same logic.
       
        Algorithm

        1) Identify all split-boundary nodes on the surface (lower-dim mesh to be converted into SplitBoundary)
       
        2) For each Element node pair, disambiguate which node the lower dim InterFace should be assigned assigned
          (relying on the nodes dim-dimensional parent having the correct node).
          
        3) assign the correct node to the lower-dim element, etc., before continueing with the splitting process.
        
        @param subdomain must be the model region unless SplitBoundary objects only exist in a subregion of Model
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::ResolveInconsistentNodeManifolds( Region<dim>& subdomain )
 {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    SPLITBOUNDARY_COMPLEX<dim>& model( static_cast<SPLITBOUNDARY_COMPLEX<dim>&>(*this) );

    // 1. Loop over all elements of the lower-dimensional finite-element input region (mesh)
    // 1.1 For all nodes, if they have been duplicated (=split) indicated by IsManifold()==true
    //     store in map<Element, ( vec<Manifold Node>, vec<node ids> ) > .
    //     note: this handles the case where the perimeter has been split by another splitboundary (adjacent mesh patch with collocated nodes)
    map<Element<dim>*, pair<vector<Node<dim>*>, vector<uint32_t> > > elements_with_manifold_nodes;
    set<Node<dim>*>                                                  manifold_on_perimeter;
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit ) {
      const uint32_t     nodes = (*eit)->Nodes();
      vector<Node<dim>*> manifold_nodes;
      vector<uint32_t>   node_ids;
      for (uint32_t i{0U}; i<nodes; ++i){         //loop over number of nodes of element
        if ( (*eit)->N(i)->IsManifold() ){
          manifold_nodes.push_back( (*eit)->N(i) );  //insert manifold node
          node_ids.push_back(i);
          if (subdomain.IsPerimeterNode((*eit)->N(i))){
              manifold_on_perimeter.insert((*eit)->N(i));
          }
        }
      } //end of node loop
      if (!manifold_nodes.empty()){      //if found manifolds
        //add element to map if manifold nodes were found.
        elements_with_manifold_nodes.insert(make_pair(*eit, make_pair( manifold_nodes, node_ids )) );
      }
    }//end of element loop


    // 2. Loop over element map
    for (auto& it : elements_with_manifold_nodes) {
        Element<dim>* lower_dim_elem = it.first;
        auto& man_nodes_info = it.second;

        Node<dim>* reference_node = nullptr;
        // Attempt to find a "stable" node (not manifold, not perimeter)
        for (auto nit = lower_dim_elem->NodesBegin(); nit != lower_dim_elem->NodesEnd(); ++nit) {
            if (!(*nit)->IsManifold()) {
                if ( isPerimeterNode<dim>((*nit)->Attribute()) == false ) {
                    reference_node = *nit;
                    break;
                }
                reference_node = *nit; // Weak fallback
            }
        }

        set<Element<dim>*> candidate_parents;

        if (reference_node != nullptr) {
            // CASE A: Standard logic using a stable reference node
            for (uint32_t p = 0; p < reference_node->Parents(); ++p) {
                Element<dim>* parent = reference_node->Parent(p);
                if (parent->IsEquidimensional()) {
                    candidate_parents.insert(parent);
                }
            }
        } else {
            /* CASE B: All-Manifold Element Fallback
               We find the higher-dimensional parent by intersecting the parent sets
               of all nodes in the lower-dim element. Only the correct mesh patch
               volume cell will contain all (the correct branches of) these nodes. */
            
            bool first_node = true;
            for (auto nit = lower_dim_elem->NodesBegin(); nit != lower_dim_elem->NodesEnd(); ++nit) {
                set<Element<dim>*> current_node_parents;
                for (uint32_t p = 0; p < (*nit)->Parents(); ++p) {
                    Element<dim>* parent = (*nit)->Parent(p);
                    if (parent->IsEquidimensional()) {
                        current_node_parents.insert(parent);
                    }
                }

                if (first_node) {
                    candidate_parents = current_node_parents;
                    first_node = false;
                } else {
                    // Intersection: candidate_parents = candidate_parents ∩ current_node_parents
                    for (auto it_p = candidate_parents.begin(); it_p != candidate_parents.end(); ) {
                        if (current_node_parents.find(*it_p) == current_node_parents.end()) {
                            it_p = candidate_parents.erase(it_p);
                        } else {
                            ++it_p;
                        }
                    }
                }
            }
        }

    if ( candidate_parents.empty() ) {
        csmp_error.Note(ERROR, "ManifoldResolution", "No common higher-dim parent found for element " + to_string(lower_dim_elem->Idx()));
        return;
    }

    // 3. Checking the node manifolds (split nodes) to make sure that their parent elements belong to different regions
    //    (the check uses the Region IDs of each nodes parents )
    //    - if not, the duplicated node (and if necessary the manifold) is removed
    for (size_t m = 0; m < man_nodes_info.first.size(); ++m )
      {
          Node<dim>* current_man_node = man_nodes_info.first[m];
          uint32_t          local_idx = man_nodes_info.second[m];
          NodeManifold<dim>* manifold = current_man_node->Manifold();
          
          set<Node<dim>*> nodes_to_assign;

          // For each potential volume cell parents (candidate_parents)
          for ( auto* vol_parent : candidate_parents ) {
                for ( uint32_t b{0}; b < manifold->Branches(); ++b ) {
                    // Check if this specific node branch is part of this volume cell
                    if (manifold->N(b)->IsParent(vol_parent)) {
                        nodes_to_assign.insert(manifold->N(b));
                    }
                }
            }

          // If nodes_to_assign > 1, the lower-dim element might be an interface between
          // two patches that haven't been fully separated yet.
          // If size == 0, the topology is disconnected; nothing needs to be done.
          if ( nodes_to_assign.empty() )
            ;
          // If nodes_to_assign = 1, the lower-dim element sits at an interface between
          // two modelsubdomains that haven't been fully separated yet.
          else if ( nodes_to_assign.size() == 1 ) {
               // disambiguating manifold node
               lower_dim_elem->Assign(local_idx, *nodes_to_assign.begin());
            }
          // If nodes_to_assign > 1 more diagnostics are needed
          else
            {
               // Verifying that the Nodes in the Manifold belong to different ModelSubDomains
               set<int32_t> subdomains_in_manifold;
               // rebuiding the nodes parent element vectors
               for ( auto* nptr : nodes_to_assign )
                 {
                    // find current parents
                    set<Element<dim>*> actual_parents;
                    for ( uint32_t i{0}; i<nptr->Parents(); ++i )
                      if ( nptr->Parent(i) && nptr->Parent(i)->IsEquidimensional() )
                        for ( uint32_t j{0}; j<nptr->Parent(i)->Nodes(); ++j )
                          if ( nptr->Parent(i)->N(j) == nptr ) {
                               actual_parents.insert( nptr->Parent(i) );
                               break;
                            }
                    // are there any?
                    if ( actual_parents.empty() ) {
                         csmp_error.Note( ERROR, "SplitBoundaryInterface::ResolveInconsistentNodeManifolds",
                                         "Parent element reassignment failed for " + to_string(lower_dim_elem->Idx()));
                         continue;
                      }
                    // updating the parent elements
                    nptr->EraseParents();
                    for ( auto* parent_ptr : actual_parents )
                      for ( uint32_t i{0}; i<parent_ptr->Nodes(); ++i )
                        if ( parent_ptr->N(i) == nptr ) {
                             nptr->Assign( i, parent_ptr );
                             break;
                          }
                    // verifying that the parents only belong to a single modelsubdomain
                    set<int32_t> parentRegion_IDs;
                    for ( uint32_t i{0}; i<nptr->Parents(); ++i )
                      parentRegion_IDs.insert( nptr->Parent(i)->Region_ID() );
                    // the nodes parents must only belong to a single region
                    if ( parentRegion_IDs.size() > 1 ) {
                         cout <<"\n"<<"Node "<< nptr->Idx() <<": "<< parseBoundary(nptr->AtBoundary()) <<": "<< parseTopology(nptr->Attribute());
                         cout <<": parents belong to multiple subdomains: ";
                         for ( const auto& region : parentRegion_IDs ) cout << model.RegionByDomainIndex(region).Name() <<" ("<< region <<") ";
                         cout << endl << endl;
                         csmp_error.Note( WARNING, "SplitBoundaryInterface::ResolveInconsistentNodeManifolds",
                                         "parent elements of Manifold Node belong to different modelsubdomains" );
                         //TODO: this node must be split
                      }
                    // recording regions avoiding duplicates
                    subdomains_in_manifold.insert( parentRegion_IDs.begin(), parentRegion_IDs.end() );
                    
                 }
                 
               if ( subdomains_in_manifold.size() != nodes_to_assign.size() ) {
                    // TODO: if the extra node cannot be separated, it must be deleted
                    // TODO: if there are more domains than nodes extra entries must be created for manifold
                    throw csmp::Exception(ERROR, "SplitBoundaryInterface::ResolveInconsistentNodeManifolds",
                                    "Ambiguity in node assignment at element " + to_string(lower_dim_elem->Idx()));
                
                 }
          } // else
      } // for
  }
 
 } // end SplitBoundaryInterface::ResolveInconsistentNodeManifolds





/* EDDY'S ORIGINAL CODE

    // 1. Loop over all elements of the lower-dimensional finite-element input region (mesh)
    // 1.1 For all nodes, if they have been duplicated (=split) indicated by IsManifold()==true
    //     store in map<Element, ( vec<Manifold Node>, vec<node ids> ) > .
    //     note: this handles the case where the perimeter has been split by another splitboundary (adjacent mesh patch with collocated nodes)
    map<Element<dim>*, pair<vector<Node<dim>*>, vector<uint32_t> > > elements_with_manifold_nodes;
    set<Node<dim>*>                                                  manifold_on_perimeter;
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit ) {
      const uint32_t     nodes = (*eit)->Nodes();
      vector<Node<dim>*> manifold_nodes;
      vector<uint32_t>   node_ids;
      for (uint32_t i{0U}; i<nodes; ++i){         //loop over number of nodes of element
        if ( (*eit)->N(i)->IsManifold() ){
          manifold_nodes.push_back( (*eit)->N(i) );  //insert manifold node
          node_ids.push_back(i);
          if (subdomain.IsPerimeterNode((*eit)->N(i))){
              manifold_on_perimeter.insert((*eit)->N(i));
          }
        }
      } //end of node loop
      if (!manifold_nodes.empty()){      //if found manifolds
        //add element to map if manifold nodes were found.
        elements_with_manifold_nodes.insert(make_pair(*eit, make_pair( manifold_nodes, node_ids )) );
        assert(manifold_nodes.size() != nodes );
      }
    }//end of element loop


    // 2. Loop over element map
    // ------------------------
    //  2.1 Find non-manifold node of element (must exist, otherwise we have split an existing splitboundary...)
    //      Loop over higher dimensional parents of non-manifold node
    //        2.2 Ask each node node in the node manifold, if they have the higher dim element as a parent
    //        2.3 insert node with matching parent into map of nodes to assign. assert(set.size() = number_manifold_nodes of lower_dim_elmt).
    //        2.4 IMPORTANT -> If manifold node is also a perimeter node, we must add it to perimeter nodes of model subdomain
    set<Node<dim>*> extra_perimeter_nodes;  //nodes that were on perimeter but split by another splitboundary
    for ( auto& it : elements_with_manifold_nodes )
      {
        // 2.1 find non-manifold node (Not manifold, and not perimeter either).  (Has Manifold()==nullptr : not even perimeter)
        //     This non-manifold node will have higher dim parents with correct nodes assigned
        Node<dim>* non_manifold_node = nullptr;
        for (auto nit=it.first->NodesBegin(); nit != it.first->NodesEnd(); ++nit ){
            if ( (*nit)->IsManifold() == false && ( (*nit)->Attribute() != PERIMETER_POINT && (*nit)->Attribute() != PERIMETER_LINE ) ) { //Should not be on perimeter (since this means node shares parents with both inside and out)
                non_manifold_node = *nit; //we found non-manifold node
                break;
              }
          }//end of search
        assert( non_manifold_node != nullptr );

        // 2.2 For each manifold node, start search for correct manifold node to assign to lower dim element
        const uint32_t parents = non_manifold_node->Parents();
        uint32_t m{0U};
        for ( auto& man_node : it.second.first ) {
              set<Node<dim>*> node_to_assign;
              NodeManifold<dim>* manifold = man_node->Manifold();
              const uint32_t branches = manifold->Branches();

              //Check if manifold node was also on perimeter
              bool add_to_perimeter_node = false;
              if (manifold_on_perimeter.find(man_node) != manifold_on_perimeter.end())
                add_to_perimeter_node = true;

              //2.3 Search over all parents of non-manifold node to see if the manifold node also shares the parent
              for (uint32_t j{0U}; j<parents;++j){ //loop over all higher dim parents
                if ( non_manifold_node->Parent(j)->IsEquidimensional() ){
                  //Find which branch of manifold shares parent
                  for (uint32_t i{0U}; i<branches; ++i){
                    if ( manifold->N(i)->IsParent( non_manifold_node->Parent(j) ) ){
                      node_to_assign.insert(manifold->N(i));            //this manifold node is on correct side of lower dim region
                      if (add_to_perimeter_node)
                        extra_perimeter_nodes.insert(manifold->N(i)); //This manifold node is also on perimeter of lower dim region
                    }
                  }//end of branch search
                }
              }//end of parent search

              // DEBUGGING: parent elements are not shared across SplitBoundaries, i.e.,
              // different manifold nodes cannot share element parents
              if ( node_to_assign.size() != 1 ) {
                   for ( const auto& node : node_to_assign ) {
                        cout <<"\n"<<"Node "<< node->Idx() <<": "<< node->Coordinate();
                        node->Manifold()->Out();
                        cout <<"\n"<<"Element to which node should be assigned "<< (*it.first).Idx() << endl;;
                     }
                   throw csmp::Exception( ERROR, "::CreateSplitBoundaryFrom:", "when searching outside Node to be assigned to dim-1 element");
                }
              // We have now found the node which should be on the lower dim element
              // 2.4 Assign node to the lower-dim element
              it.first->Assign( it.second.second[m], *node_to_assign.begin() );
              m++; //increment manifold index (allows to get correct node id information)
           }

    }//end of elements_with_manifold_nodes loop


*/





/**
    Creates SplitBoundar(ies) from lower dimensional region without the need for a user to create Boundary objects first. The underlying region is removed in the process.

    @brief 1) Checks dim-1 region exists
           2) Checks no Nodes on the dim-1 region have the label INTERNAL before proceeding
           3) Counts Unique Regions and creates property "region identifyier" for all elements within each region with a region-specific integer
              (to be identified after node duplication for region updating)
           4) For each Element in the dim-1 Region, it constructs FaceConstructionData, which stores:
                   i)    Inner/Outer Parents of dim-1 element,
                   ii)   Faces of inner/outer parent elements that are sharing the dim-1 element.
                   iii)  Region_identifier (Element property) of the Inside and Outside Parent Element -> (materials_)
                              (THIS HOWEVER IS NOT ALL THE ELEMENTS THAT SHARE THE NODE)
                   iii)  Region_Identifier of the dim-1 element (material_)
                   iv)   PatchNumber which is the unique identifier of the subregions determined by the n_juxtoposition
                            (Patch is the same if you have the same inside/outside materials_ (region identifier))

                4.2) Breaks dim-1 region into smaller subregions based on t
           5) Breaks dim-1 Region into smaller subregions, and gives each one a unique name

*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
pair<set<string>,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom( const char* dim_1_region,
                                                                                                     bool retain_elmts_as_intervening_elements )
 {
    SPLITBOUNDARY_COMPLEX<dim>& model( static_cast<SPLITBOUNDARY_COMPLEX<dim>&>(*this) );

    cout <<"\nSplitBoundaryInterface<"<< dim <<",Model>::CreateSplitBoundaryFrom: forming splitboundary(ies) from region: '";
    cout << dim_1_region <<"'...\n";
    
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    // 1. Verify input lower-dimensional region object from which splitboundary shall be created: must be lower dimensional and lie on inside of model
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    const string creation_failed( string("CreateSplitBoundaryFrom(") +  dim_1_region +") failed");
    // does the parent region exist
    if ( model.ContainsRegion(dim_1_region) == false ) {
          ErrorHandler::Instance().Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", dim_1_region, "does not exist; nothing was done." );
          return make_pair( set<string>({creation_failed}), false );
      }
    // 1.1 do such split boundaries already exist ?
    const set<string> intersected_regions{dim_1_region};
    set<string>       pre_existing_splitboundaries;
    if ( FindSplitBoundaryByNames( intersected_regions, pre_existing_splitboundaries ) > 0 ) {
          string error_info;
          for ( auto& it : pre_existing_splitboundaries ) {
               error_info += it;
               error_info +=", ";
            }
          ErrorHandler::Instance().Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryFrom:",
                                           error_info, "boundaries are already contained in this model." );
          return make_pair( set<string>({creation_failed}), false );
      }
    // 1.2 does the model contain unique regions
    if ( model.UniqueRegions() < 1U ) {
          ErrorHandler::Instance().Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", "model contains no unique regions; cannot proceed." );
          return make_pair( set<string>({creation_failed}), false );
      }
    // 1.3 verifying that we are indeed dealing with a region of surface or line elements only and that their normals all point into same direction
    Region<dim>&  subdomain(model.Region(dim_1_region));
    if ( checkNeighborNormalsForConsistentOrientation( subdomain ) == false ) {
         ErrorHandler::Instance().Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", dim_1_region,
                                       "region appears to have inconsistent surface-normal orientations; nothing was done." );
         return make_pair( set<string>({creation_failed}), false );
      }
    // 1.4 checking that the region is not already an internal boundary
    size_t nodes_flagged_internal_boundary{0U}, nodes_flagged_external_boundary{0U};
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         if ( (*nit)->AtBoundary() == INTERNAL ) nodes_flagged_internal_boundary++;
         else if ( (*nit)->AtBoundary() != NOT ) nodes_flagged_external_boundary++;
      }
    if ( nodes_flagged_internal_boundary >= subdomain.Nodes() - nodes_flagged_external_boundary ) {
         ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", dim_1_region,
                                                   "region may already be a boundary; nothing was done." );
         return make_pair( set<string>({creation_failed}), false );
      }
    
    // 1.5 checking that the region is not located at the model boundary
    size_t boundary_elements{0U};
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
      for ( uint32_t i{0U}; i<(*eit)->Neighbors(); ++i ) {
           const BOX_BOUNDARY bflag = (*eit)->AtBoundary(i);
           if ( bflag != NOT and bflag != INTERNAL and bflag != IRREGULAR ) boundary_elements++;
        }
    if ( boundary_elements == subdomain.Cells() ) {
         ErrorHandler::Instance().Note( WARNING, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", dim_1_region,
                                                   "region appears to lie at the model boundary; nothing could be done." );
         return make_pair( set<string>({creation_failed}), false );
      }
      
    // TODO: use domain index tags instead of creating a new variable
    // 1.6 creating region labels and tagging the regions with unique integer identifiers
    //     if "region identifier" is already defined it is assumed that it has already been initialised as well
    const string    region_tag("region identifier");  //Note: this name is hard coded in MeshManager::ReplaceElementsByInterface()
    vector<string>  region_names;
    if ( !model.CheckRegionIdentifierIsUpToDate(region_tag.c_str()) ) {
         // for each of labels created (0..regions-1), region_names remembers which region the label refers to
         if ( model.CountAndLabelUniqueRegions( region_tag.c_str(), region_names ) == 1U )
            ErrorHandler::Instance().Note( INFO, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", region_tag,
                                                  "is single valued; so there is only one patch." );
      }
    else { // assigning region names
         region_names.reserve( static_cast<size_t>(abs(distance( model.UniqueRegionsBegin(),model.UniqueRegionsEnd()))) );
         for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit )
           region_names.push_back( (*rit).first );
      }
    const csmp::Index mtrl_key = model.Database().StorageKey(region_tag.c_str());
         
         


    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 2. Preprocessing to FaceConstructionData and getting MeshManager object to create InterFace objects for all elements in lower-dim region
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    //  2.1 creating the required objects
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    const size_t                       new_interfaces_required(subdomain.Cells());
    vector<FaceConstructionData<dim>>  iface_construction_vector;
    iface_construction_vector.reserve(new_interfaces_required);
    const size_t n_original_faces(model.Mesh().Interfaces());
    const size_t n_original_elmts(model.Mesh().Elements());

    // looping over the region, identifying and recording the juxtaposition relationships
    map<pair<long,long>,uint32_t>      patches;
    map<long,string>                   patch_names;
    string                             patch_name;
    uint32_t                           n_juxtapositions{0}; // among regions

    // 2.2 Preprocessing the entire lower dimensional region collecting the Data needed for InterFace construction
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
      {
        //making FaceConstructionData
        FaceConstructionData  fdata( higherDimensionalNeighbors( *(*eit), mtrl_key ) );

        // 2.1.2 recording which category of juxtaposition element fall into, naming it and assigning a patch number
        auto it = patches.insert( make_pair(fdata.Materials(),n_juxtapositions) );
        // incrementing number of juxtapositions and corresponding patch names
        if ( it.second == true ) {                                          //If new element was inserted, then
             fdata.PatchNumber( (*it.first).second );                       //Set patch number to n_juxtapositions (this should be unique...)
             patch_name = CreateSplitBoundaryNameFrom( fdata, region_names ); //Make a name based on fdata( materials, juxta, region name)
             patch_names.insert( make_pair(n_juxtapositions,patch_name) );
             n_juxtapositions++;
          }
        fdata.PatchNumber( (*it.first).second );          //name set t

        //Adding to vector
        iface_construction_vector.push_back(fdata);
      }
    assert( iface_construction_vector.size() == subdomain.Cells() );


    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 2.3 Creating Interfaces from vector of FaceConstructionDataInterFace using MeshManager (expects FaceConstructionData with well defined perimeter).
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    set<size_t>     region_material_ids;   // needed for selective rebuilding of regions with OUTSIDE elements (vector is configured by MeshManager)
    set<Node<dim>*> extra_perimeter_nodes; // hopefully there are none

    // disconnecting any equidimensional elements surrounding the element patch even if it will not be deleted
    updateHaloElementConnectivity<dim>( subdomain.CellVector().begin(), subdomain.CellVector().end() );
    // Creating interfaces for the entire input region, considering its perimeter to avoid node duplication where it terminates inside another region
    // split boundary intersections are also handled .... Hopefully
    vector<InterFace<dim>*> iface_vector = model.Mesh().ReplaceElementsByInterFaces( model.Database(), // <<<<<<< INTERFACE CREATION
                                                             iface_construction_vector.begin(),
                                                             iface_construction_vector.end(),
                                                             subdomain.PerimeterNodesBegin(),
                                                             subdomain.NodesEnd(),
                                                             extra_perimeter_nodes, // will contain old perim-node - manifold node pairs (only need the manifold node).
                                                             region_material_ids,
                                                             retain_elmts_as_intervening_elements ); // <<<<<<<<<< KEEP ORIGINAL ELEMENTS

    assert(iface_vector.size() == iface_construction_vector.size());
    assert( (*iface_vector.begin())->InnerParent() == iface_construction_vector.begin()->InnerElement()) ;
    assert( (*iface_vector.back()).InnerParent() == iface_construction_vector.back().InnerElement()) ;
    
    // while the new interfaces were already connected with one another by ReplaceElementsByInterFaces this updates their neighborhood
    model.Mesh().UpdateConnectivity( iface_vector.begin(), iface_vector.end() );
    // TOPOTYPE perimeter flagging is not done at this stage because the interfaces have to be partitioned
    // into Split(sub)Boundaries first

    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 3. Determine number of mesh patches that will be turned into new splitboundary objects.
    //    The output of this step will be a map of FaceConstructionData in which the names of new Split(sub)boundaries are the keys.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 3.2 creating labeled boundary patches from the interface-defining data
    // ----------------------------------------------------------------------
    // 3.2.1 making a map 'patch_numbers' from 'patch_names' to search for patch identifiers
    map<string,size_t>  patch_numbers;
    for ( const auto& it : patch_names )
      patch_numbers.insert( make_pair( it.second, static_cast<uint32_t>(it.first) ) );

    // 3.2.2 building new map of InterFace objects grouped by patch names
    map<string,vector<InterFace<dim>*>>  patch_data;
    vector<InterFace<dim>*>              empty_vec;
    for (const auto& val: patch_names | views::values)
      patch_data.insert( make_pair( val, empty_vec ) );

    // 3.2.3 inserting the patch identifiers into the vectors in the map
    for ( auto& pit : patch_data )
      {
         assert( patch_numbers.find(pit.first) != patch_numbers.end() );
         const size_t patch_number((*patch_numbers.find(pit.first)).second);

         // looping over all face data assigning the ones that are suitable
         const size_t ifaces = iface_construction_vector.size();
         // reserving storage
         pit.second.reserve(ifaces);

         for ( size_t i{0U}; i < ifaces; ++i )
           if ( iface_construction_vector[i].PatchNumber() == patch_number )
             pit.second.push_back( iface_vector[i] );  //assumes that iface_vector[i] has InterFace which corresponds with iface_construction_vector[i]
      }
    // 3.2.4 trimming excess storage of the face-data vectors
    for ( auto& pit : patch_data ) pit.second.shrink_to_fit();


    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 4. Creating SplitBoundary objects for each of the mesh patches established above
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // (the interfaces are already connected but this step identifies patch interior and perimeters)
    for ( long i{0U}; i < static_cast<long>(patch_names.size()); ++i )
      AddSplitBoundary( patch_names[i].c_str(), patch_data[ patch_names[i] ].begin(), patch_data[ patch_names[i] ].end() );
      
#ifdef DEBUG
    cout <<"\n\n"<<"SplitBoundaryInterface<"<< dim <<">::CreateSplitBoundaryFrom:";
    cout << "\n\t\t"<<"Added "<< model.Mesh().Interfaces() - n_original_faces <<" interfaces to mesh.";
#endif

    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 5. Removing lower-dimensional input region if so requested.
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    const bool remove_elmts{ !retain_elmts_as_intervening_elements };
    // disconnecting the intervening element pointers in the new SplitBoundary from them elements that will be deleted
    if ( remove_elmts ) {
         for ( long i{0U}; i < static_cast<long>(patch_names.size()); ++i ) {
              auto& split_boundary = SplitBoundary( patch_names[i] );
              for ( auto& it : split_boundary.CellVector() )
                 it->UnAssignInterveningElement();
           }
         // removing the elements and the region that contained them
         model.RemoveRegion( dim_1_region, retain_elmts_as_intervening_elements );
      }

#ifdef DEBUG
    if ( remove_elmts ) cout << "\n\t\t"<<"Removed "<< n_original_elmts - model.Mesh().Elements()  <<" elements from the mesh."<< endl;
#endif

    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 6. Flagging the regions on the outside of the new split boundaries for update of their connectivity
    // because they will contain new nodes
    for (size_t i : region_material_ids )
      model.Region( region_names[ i ] ).ScheduleForRebuild();

    //update outside regions
    model.Region("Model").ScheduleForRebuild();
    model.UpdateRegions();
   
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    // 7. extra diagnostics and output of boundary names
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    if ( patch_names.empty() ) {
         ErrorHandler::Instance().Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", dim_1_region, "no SplitBoundary patches could be created." );
         return make_pair( set<string>({}), false );
      }
      
    set<string> split_boundary_names;
    for ( const auto& it : patch_names ) split_boundary_names.insert( it.second );

    if ( patch_names.size() > split_boundary_names.size() ) {
         ErrorHandler::Instance().Note( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryFrom:", dim_1_region, "not all of the created patchnames are unique." );
         return make_pair( set<string>({}), false );
      }

    return make_pair( split_boundary_names, true );
       
 } // end CreateSplitBoundaryFrom








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
   
    const csmp::Region<dim>&  region1(modelComplex->Region(region1_name)); // inside region
    csmp::Region<dim>&        region2(modelComplex->Region(region2_name)); // not const because outside region must be rebuilt
    
    // 2. Is there a shared interface? - if so InterFace objects are created along it
    // ------------------------------------------------------------------------------
    vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > > matched_elmts;
    sharedPerimeterCells( region1, region2, matched_elmts );
    
    if ( matched_elmts.empty() ) {
         string message( string(" input regions '") + region1_name + "' and '" + region2_name +"'");
         csmp_error.Note( WARNING, "SplitBoundaryInterface::CreateSplitBoundaryBetween:",
                          message, "may share a node but do not share any faces; is there already a SplitBoundary between them?");
         return make_pair("split boundary not created",false);
      }
        
    // SKM: new: recording elements touching the perimeter from the outside so that they will need to be rebuilt
    unordered_set<Element<dim>*>  halo_elements;
    haloElements( region1, region2, halo_elements );

    // creating the necessary interfaces and nodes, and establishing the connectivity of the faces while updating the mesh
    const bool create_manifolds_on_perimeter{true};
    vector<InterFace<dim>*>  interfaces = modelComplex->Mesh().CreateInterfacesBetweenNodeSharingElements( modelComplex->Database(),
                                                                                                           matched_elmts,
                                                                                                           create_manifolds_on_perimeter,
                                                                                                           halo_elements );
    // Flagging the topology of interface nodes
    updateSplitBoundaryTopoFlags<dim>( interfaces );
    
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
             region2.ScheduleForRebuild();
             modelComplex->UpdateRegions();
             // all done
             return make_pair(split_boundary_name,true);
          }
      }
    else throw csmp::Exception( ERROR, "SplitBoundaryInterface::CreateSplitBoundaryBetween:",
                                split_boundary_name, "split boundary already exists. Nothing was done.");
    
    // update region outside of new SplitBoundary because its nodes have changed
    region2.RebuildSubDomainAfterChangeOfNodeVector();

    //update outside regions
    // RebuildSubDomainAfterChangeOfCellVector() is not sufficient because the model needs to be rebuild as well
    modelComplex->UpdateRegions();

    return make_pair("split boundary not created",false);

 } // end CreateSplitBoundaryBetween




/**
 * Replaces specific boundary keywords with "REGION" or appends "_REGION".
 * Avoids string::contains() for MSVC compatibility.
 */
static void regionNameFromSplitBoundaryName( string& name ) {
    // Define the target substrings
    const std::vector<std::string_view> targets = {
        "SPLITBOUNDARY",
        "SPLIT_BOUNDARY",
        "split_boundary",
        "split-boundary"
    };

    bool replaced = false;

    for (const auto& target : targets) {
        size_t pos = name.find(target);
        
        // While loop in case the substring appears multiple times
        while (pos != std::string::npos) {
            name.replace(pos, target.length(), "REGION");
            replaced = true;
            
            // Advance position to avoid infinite loop if "REGION" contained the target
            pos = name.find(target, pos + 6); // 6 is length of "REGION"
        }
    }

    // If none of the targets were found, append the suffix
    if (!replaced) {
        name += "_REGION";
    }
}




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
     
    SPLITBOUNDARY_COMPLEX<dim>*     model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
    csmp::SplitBoundary<dim>&       splitBoundary( model->SplitBoundary(split_boundary) );
    csmp::MeshManager<dim>&         mesh( model->Mesh() );
    const LocalVariables&           nlvars( model->Database().LocalVariablesAt( NODE ) );
    const LocalVariables&           lvars( model->Database().LocalVariablesAt( ELEMENT ) );
    const IntegrationPointVariables& ivars( model->Database().IntegrationPointVariablesAt( ELEMENT ) );
     
    // 1. Create a map to associate the original inside nodes with the newly duplicated nodes
    // ---------------------------------------------------------------------------------------------
    pair<vector<Node<dim>*>,size_t>  inside_nodes = splitBoundary.InsideNodes();
     
    // Use an unordered_map instead of overwriting the global Node::Idx()
    std::unordered_map<const Node<dim>*, Node<dim>*> old_to_new_node_map;
    old_to_new_node_map.reserve(inside_nodes.first.size());

    // Duplicate nodes and populate the map
    for ( auto* old_node : inside_nodes.first ) {
        // NB: here new manifolds are generated or the new nodes are inserted into existing manifolds
        Node<dim>* new_node = mesh.Duplicate( old_node, nlvars );
        old_to_new_node_map[old_node] = new_node;
    }
      
    // 2. Creating elements within InterFace objects with node-numbering matching the INNER parent
    // --------------------------------------------------------------------------------------------------------------------------
    vector<Element<dim>*> elmt_pointers; // new elements
    elmt_pointers.reserve(splitBoundary.CellVector().size());

    for ( auto& it : splitBoundary.CellVector() ) {
        // Guard clause: Only add elements to interfaces that don't already have one
        // (Assuming a method like HasInterveningElement() or similar exists in your API)
        if ( it->HasInterveningElement() ) {
            continue; 
        }

        vector<Node<dim>*>  nodes;  
        nodes.reserve(it->FE()->Nodes());
        
        for ( uint32_t i{0U}; i<it->FE()->Nodes(); i++ ) {
            const Node<dim>* target_inside_node = it->N(i, INSIDE);
            
            // Safely look up the new node without relying on corrupted global IDs
            auto map_it = old_to_new_node_map.find(target_inside_node);
            if (map_it != old_to_new_node_map.end()) {
                nodes.push_back( map_it->second );
            } else {
                csmp_error.Note( FATAL_ERROR, "InsertRegionIntoSplitBoundary",
                                 "InterFace INSIDE node not found in unique InsideNodes list!" );
            }
        }
        
        // New interior and perimeter elements get connected to their nodes and the middle element
        elmt_pointers.push_back( mesh.AddInterveningElement( it, lvars, ivars, nodes, material_id ) );
    } 
       
    // If no new elements were added (e.g., region was already fully populated), exit cleanly
    if (elmt_pointers.empty()) {
        return make_pair(splitBoundary.Name(), false);
    }

    // 3. Establishing neighbor connectivity among the new elements
    // ------------------------------------------------------------
    mesh.template BuildConnectivity<Element>( elmt_pointers.begin(), elmt_pointers.end() );
    mesh.ConnectNodesToParentsAndNeighbors( elmt_pointers.begin(), elmt_pointers.end() );

    // 4. Construct the new unique region between the interface elements in the model
    // -----------------------------------------------------------------------------
    string region_name = splitBoundary.Name();
    regionNameFromSplitBoundaryName( region_name ); 
     
    const bool  unique_map(true);
    model->FormRegionFrom( region_name.c_str(), elmt_pointers.begin(), elmt_pointers.end(), unique_map );
     
    // Grab the newly formed region to utilize its sorted perimeter nodes
    auto& new_region = model->Region(region_name);

    // 4.1 Flag the perimeter nodes of this new lower-dimensional region
    for (auto nit = new_region.PerimeterNodesBegin(); nit != new_region.NodesEnd(); ++nit) {
        (*nit)->AtBoundary(INTERNAL);
        if constexpr (dim == 3) {
            (*nit)->Attribute(PERIMETER_LINE);  // Edge of the internal surface patch
        } else if constexpr (dim == 2) {
            (*nit)->Attribute(PERIMETER_POINT); // Tip of the internal line polyline
        }
    }

    // 4.2 Flag the interior nodes of this new lower-dimensional region
    for (auto nit = new_region.NodesBegin(); nit != new_region.PerimeterNodesBegin(); ++nit) {
        (*nit)->AtBoundary(INTERNAL);
        if constexpr (dim == 3) {
            (*nit)->Attribute(INTERIOR_SURFACE); // Nodes moving within the surface plane
        } else if constexpr (dim == 2) {
            (*nit)->Attribute(INTERIOR_LINE);    // Nodes sliding along the internal line
        }
    }

    // Add new unique region to model region
    model->Region("Model").Add( new_region );
     
    return make_pair( region_name, true );    
} // end InsertRegionIntoSplitBoundary


/* REPLACED Idx() based version (not reliable)

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
     
     SPLITBOUNDARY_COMPLEX<dim>*      model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
     csmp::SplitBoundary<dim>&        splitBoundary( model->SplitBoundary(split_boundary) );
     csmp::MeshManager<dim>&          mesh( model->Mesh() );
     const LocalVariables&            nlvars( model->Database().LocalVariablesAt( NODE ) );
     const LocalVariables&            lvars( model->Database().LocalVariablesAt( ELEMENT ) );
     const IntegrationPointVariables& ivars( model->Database().IntegrationPointVariablesAt( ELEMENT ) );
     
     // 1. creating unique set of nodes matching those on the inside of the SplitBoundary in position
     // ---------------------------------------------------------------------------------------------
     
     // 1.1 to start with, a unique set of inside nodes is created for duplication
     pair<vector<Node<dim>*>,size_t>  inside_nodes = splitBoundary.InsideNodes();
     
     // 1.2 the indices of these nodes are numbered so that the new duplicated node vector can be accessed for assignments
     size_t counter{0U};
     for ( auto& nit : inside_nodes.first )
       nit->Idx( counter++ );
     
     // 1.3 now all nodes are duplicated AND inserted into the corresponding manifolds
     //     (NB: even the tip nodes now become manifolds)
     vector<Node<dim>*>    node_pointers; // to the new nodes
     for ( size_t i{0U}; i<inside_nodes.first.size(); i++ )
       // NB: here new manifolds are generated or the new nodes are inserted into existing manifolds
       node_pointers.push_back( mesh.Duplicate( inside_nodes.first[i], nlvars ) );
       //                       NODE CREATION ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
     // 2. creating elements within InterFace objects with node-numbering matching that of corresponding INNER parent element face
     // --------------------------------------------------------------------------------------------------------------------------
     vector<Element<dim>*> elmt_pointers; // new elements
     counter = 0U;
     for ( auto& it : splitBoundary.CellVector() ) {
           // nodes
           vector<Node<dim>*>  nodes;  nodes.reserve(4);
           for ( uint32_t i{0U}; i<it->FE()->Nodes(); i++ )
             nodes.push_back( node_pointers[ it->N(i,INSIDE)->Idx() ] );
           // new interior and perimeter elements get connected to their nodes and the middle element
           elmt_pointers.push_back( mesh.AddInterveningElement( it, lvars, ivars, nodes, material_id ) );
       }   //         ELEMENT GENERATION ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       
     // 3. establishing neighbor connectivity among the new elements
     // ------------------------------------------------------------
     // 3.1 connect elements with their neighbors
     mesh.template BuildConnectivity<Element>( elmt_pointers.begin(), elmt_pointers.end() );

     // 3.2 establish node-to-node and parent connectivity
     mesh.ConnectNodesToParentsAndNeighbors( elmt_pointers.begin(), elmt_pointers.end() );

     // 4. construct the new unique region between the interface elements in the model
     //    given it the same name as the split boundary but calling it region instead
     // -----------------------------------------------------------------------------
     string region_name = splitBoundary.Name();
     regionNameFromSplitBoundaryName( region_name );
     
     const bool  unique_map(true);
     model->FormRegionFrom( region_name.c_str(), elmt_pointers.begin(), elmt_pointers.end(), unique_map );
     
     // add new unique region to model region
     model->Region("Model").Add( model->Region(region_name) );
     
     return make_pair( region_name, true ); 
     
 } // InsertRegionIntoSplitBoundary
 */





    /// Creates a single lower-dimensional mesh region taking into account all split boundaries objects; returning its name and whether this operation was successful
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFromInterfaces( const char* splitboundary_name )
 {
     SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));

     // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
     string bname = splitboundary_name;
     bname       += "_SPLITBOUNDARY";

     // collecting all interfaces into a vector
     vector<InterFace<dim>*> interface_ptrs;
     interface_ptrs.reserve( model->Mesh().Interfaces() );
     for ( auto itf=model->Mesh().InterfacesBegin(); itf!=model->Mesh().InterfacesEnd(); ++itf )
       interface_ptrs.push_back( &(*itf) );
      
     return AddSplitBoundary( bname.c_str(), interface_ptrs.begin(), interface_ptrs.end() );

 } // end CreateSplitBoundaryFromAllInterfaces








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
  
  // what will become the first new Element generated
  size_t n_elmts_before = splitboundaryComplex->Mesh().Elements();
  
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
    
  // rebuilding the connectivity among all the new surface elements
  size_t n_elmts_after = splitboundaryComplex->Mesh().Elements();
  size_t n_new_elmts = n_elmts_after - n_elmts_before;
  if ( n_new_elmts > 0 ) {
       cout <<"\n"<<"SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertLowerDimensionalRegionsIntoSplitBoundaries: ";
       cout <<" created "<< n_new_elmts << " lower-dimensional intervening elements forming the regions:\n\t";
       // getting pointers to the new elements created in the regions
       vector<Element<dim>*> new_elmts;
       for ( const auto& it : new_regions ) {
            cout << it <<" ";
            auto& subdomain = splitboundaryComplex->Region(it);
            // appending the cells of the region
            new_elmts.insert( new_elmts.end(), subdomain.CellsBegin(), subdomain.CellsEnd() );
         }
       cout << endl << endl;
       splitboundaryComplex->Mesh().template BuildConnectivity<Element>( new_elmts.begin(), new_elmts.end() );
    }
  else csmp_error.Note( ERROR, "SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertLowerDimensionalRegionsIntoSplitBoundaries",
                       "no new intervening elements could be created" );
    
  return new_regions;

} // end InsertLowerDimensionalRegionsIntoSplitBoundaries







template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateNonUniqueSplitBoundaryGroup( const set<string>& input_split_boundaries,
                                                                                            const char* ensemble_split_boundary )
{
  string output_domain( ensemble_split_boundary );

  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  if ( ContainsSplitBoundary( ensemble_split_boundary ) ) {
    csmp_error.Note( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateNonUniqueSplitBoundaryGroup:",
                     ensemble_split_boundary, "output domain already exists, adding an underscore to its name." );
    output_domain += "_";
    cout << output_domain << endl;
  }

  if ( input_split_boundaries.empty() ) {
    csmp_error.Note( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateNonUniqueSplitBoundaryGroup:",
                    "No input split boundaries were specified; merge could not be performed; no new SplitBoundary.",
                    ensemble_split_boundary );
    return;
  }

  typename map<string, csmp::SplitBoundary<dim> >::const_iterator  iter;

  // collecting element indexes from input regions into set for output
  vector<InterFace<dim>*>  interface_ptrs;
  for ( const auto& it : input_split_boundaries ) {
      // finding the region in the region list
      if ( (iter = splitBoundaryMap_.find( it )) != splitBoundaryMap_.end() ) {
          // outputting the ids of the member elements of the region
          interface_ptrs.reserve( interface_ptrs.size() + (*iter).second.Cells() );
          for ( auto eit = (*iter).second.CellsBegin(); eit != (*iter).second.CellsEnd(); eit++ )
            interface_ptrs.push_back( const_cast<InterFace<dim>*>(*eit) );
        }
      else csmp_error.Note( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateNonUniqueSplitBoundaryGroup:",
                            it, "SplitBoundary does not exist and was therefore not considered." );
    }

  if ( !interface_ptrs.empty() )
    {
      // eliminating duplicate entries from pointer vector
      sort( interface_ptrs.begin(), interface_ptrs.end() );
      interface_ptrs.erase( unique( interface_ptrs.begin(), interface_ptrs.end() ), interface_ptrs.end() );
      
      // rebuilding the interface connectivity
      SPLITBOUNDARY_COMPLEX<dim>& model( static_cast<SPLITBOUNDARY_COMPLEX<dim>&>(*this) );
      model.Mesh().template BuildConnectivity<InterFace>( interface_ptrs.begin(), interface_ptrs.end() );
      
      AddSplitBoundary( ensemble_split_boundary, interface_ptrs.begin(), interface_ptrs.end() );
   }
  else
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateNonUniqueSplitBoundaryGroup",
                          "No InterFace objects in input; SplitBoundary group could not be assembled",
                           output_domain.c_str() );

  cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateNonUniqueSplitBoundaryGroup: new split boundary '" << output_domain << "' formed successfully.\n";

} // end CreateNonUniqueSplitBoundaryGroup





/**
Prints current SplitBoundaries to screen
*/
template<uint32_t dim, template<uint32_t> class SPLITBOUNDARY_COMPLEX>
long SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitBoundariesOut() const
  {
     const SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<const SPLITBOUNDARY_COMPLEX<dim>*>(this));

     cout <<"\nSplitBoundaryInterface<"<< dim <<",SplitBoundary<InterFace>>::SplitBoundariesOut: ";
     if ( SplitBoundaries() == 0U ) {
          cout <<"\tmodel does not contain any split boundaries.\n\n";
          return 0U;
       }
     cout <<"split boundaries of ";
     if ( model->BoxShaped() ) cout <<"box-shaped model:\n";
     else cout <<"irregularly-shaped model:\n";
     for ( auto bit=SplitBoundariesBegin(); bit!=SplitBoundariesEnd(); ++bit ) {
          cout <<"\n\t"<< (*bit).first <<": d-idx:"<< (*bit).second.DomainIndex() <<", box-flag: "<< parseBoundary( (*bit).second.AtBoundary() );
          cout <<" "<< (*bit).second.Cells() <<" interfaces, ";
          // in 3D a boudary is a surface
           if constexpr ( dim == 3U ) {
                const pair<CELL_SHAPE,bool> cell_type = (*bit).second.SingleCellShapeDomain();
                assert( cell_type.second == true );
                if ( cell_type.first == SURFACE ) {
                     cout <<"area (m2): "<< (*bit).second.Area();
                     cout <<", perimeter length (m): "<< (*bit).second.Perimeter();
                  }
                else cout <<"length (m): "<< (*bit).second.Area();
             }
           if constexpr ( dim == 2U )
             cout <<"length (m): "<< (*bit).second.Length();
       }
     cout << endl << endl;
     cout.flush();

     return distance( SplitBoundariesBegin(), SplitBoundariesEnd() );

} // end Out





template class SplitBoundaryInterface<1U, Model>;
template class SplitBoundaryInterface<2U, Model>;
template class SplitBoundaryInterface<3U, Model>;

} // csmp


