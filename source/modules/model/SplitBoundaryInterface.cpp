#include "SplitBoundaryInterface.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"
#include "MeshManager.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "IsoparametricLinearLineElement.h" 
#include "IsoparametricLinearTriangle.h" 
#include "IsoparametricLinearQuadrilateral.h" 

using namespace std;

namespace csmp {

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
/*
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

  return split_boundary_name;

} // end CreateSplitBoundaryName
*/





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
    std::vector<InterFace<dim>* > elementsToDelete( splitboundary.CellVector().begin(), splitboundary.CellVector().end() );

    if ( !elementsToDelete.empty() )
    {
      // remove redundant InterFaces from existing SplitBoundaries
      for ( splitBoundaryIterator
            sbit = splitBoundaryMap_.begin();
            sbit != splitBoundaryMap_.end();) {

        if ( removeVectorElements( (*sbit).second.CellVector(), elementsToDelete ) > 0 )
        {
          if ( !(*sbit).second.CellVector().empty() )
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
  if ( this->SplitBoundaries() == 0 ) return false;

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

  std::cout << "\nSplitBoundaryInterface<" << dim << ">::OutputSplitBoundariesToBinary: split boundaries are being written to binary file...\n";
  for ( auto bit( this->SplitBoundariesBegin() ); bit != this->SplitBoundariesEnd(); bit++ )
  {
    // writing name of splitboundary
    skm_C_fwrite( fp, bit->first.c_str() );
    std::cout << (*bit).first << " ";
    // splitboundary
    if ( !bit->second.Out( fp ) )
      throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary:", "SplitBoundary could not be stored." );
  }

  // cleaning up
  fp.close();
  std::cout << "\n\nSplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary: Split boundaries have been successfully written to: '";
  std::cout << bin_file << "'" << std::endl;
  return true;

} // end OutputSplitBoundariesToBinary





template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InputSplitBoundariesFromBinary( const char* file_name, const std::set<std::string>* subset_variables )
{
  std::string bin_file( file_name );
  fstream fp( bin_file.c_str(), ios::in | ios::binary );
  if ( !fp.is_open() ) return false;

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
    @author JCK updated 15/05/2019
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
std::pair<std::set<std::string>,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries()
{
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  ErrorHandler&               csmp_error( ErrorHandler::Instance() );

  // 1. detecting potential SplitBoundaries
  // --------------------------------------
  set<pair<string, string>>  region_pairs;
  std::vector<std::string> regions;
  regions.reserve( splitboundaryComplex->UniqueRegions() );

  const pair<int32, int32>  model_dim = splitboundaryComplex->Region( "Model" ).SpatialDimensions();
  for ( typename std::map<std::string, csmp::Region<dim> >::iterator
        it = splitboundaryComplex->UniqueRegionsBegin(); it != splitboundaryComplex->UniqueRegionsEnd(); ++it ) {
    const pair<int32, int32>  sub_dim = (*it).second.SpatialDimensions();
    if ( sub_dim.second == model_dim.second ) // check whether the highest dimension of the region is equal to the highest dimension of the model
      regions.push_back( (*it).second.Name() );
  }

  if ( regions.size() > 1 )
    for ( size_t i = 0; i < regions.size(); i++ ) {
      const csmp::Region<dim>&  gref1( splitboundaryComplex->Region( regions[i] ) );
      for ( size_t j = i + 1; j < regions.size(); j++ ) {
        const csmp::Region<dim>&  gref2( splitboundaryComplex->Region( regions[j] ) );
        if ( regions[i] == regions[j] ) continue;
        const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
        if ( shared_nodes > 0 )
          region_pairs.insert( make_pair( regions[i], regions[j] ) );
      }
    }

  if ( region_pairs.size() == 0 ) {
    csmp_error.notice( WARNING, "SplitBoundaryInterface::DetectAndCreateSplitBoundaries:",
                        "node-coordinate matched faces / internal boundaries could not be detected; nothing was done." );
    return make_pair( set<string>(), false );
  }

  // backup the original nodes and elements since the model might not be contiguous anymore after spliting
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( &splitboundaryComplex->Mesh(), nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // 2. creating splitboundaries for each of the discovered regions
  // --------------------------------------------------------------------------------
  set<std::string>  split_boundary_names;
   
  for ( auto it : region_pairs ) {
      pair<string,bool>  split_boundary = splitboundaryComplex->InsertSplitBoundary( it.first.c_str(), it.second.c_str() );
      // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
      if ( split_boundary.second == true ) {
           std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries: ";
           std::cout << " splitboundary created successfully.\n";
           splitboundaryComplex->Mesh().Update( nodes, elmts );
           split_boundary_names.insert( split_boundary.first ); 
        }
      else {
          csmp_error.notice( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::DetectAndCreateSplitBoundaries:",
                                      "splitboundary already exists. Nothing was done." );
          return make_pair( set<string>(), false );
        }
    }

  return make_pair( split_boundary_names, true );

} // end DetectAndCreateSplitBoundaries







/**
    Creates SplitBoundar(ies) from lower dimensional region without the need for a user to create Boundary objects first. The underlying region is removed in the process.
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
pair<set<string>,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom( const char* dim_1_region )
 {
    SPLITBOUNDARY_COMPLEX<dim>* model( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  
    // 1. Converting the lower dimensional region into a single Boundary or multiple Boundaries (patches of juxtaposed rocks)
    const bool remove_original_region(true);
    pair<set<string>,bool> boundary_names = model->CreateInternalBoundaryFrom( dim_1_region, remove_original_region );
    if ( boundary_names.second == false ) 
      return boundary_names;
    
    set<string>  split_boundary_names;
    for ( set<string>::const_iterator it=boundary_names.first.begin(); it!=boundary_names.first.end(); ++it ) {
          Boundary<dim>& boundary = model->Boundary( (*it) );
          split_boundary_names.insert( CreateSplitBoundaryFrom( boundary ).first );
       }
       
    return make_pair( split_boundary_names, true );
       
 } // end CreateSplitBoundaryFrom




/**
     CreateSplitBoundaryFrom( const Boundary<dim>& );
     
     @attention the input boundary is removed in the process.
     
     @return boolean indicating whether the method was able to create a singe Split boundary and its name
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
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
  std::string  splitboundaryName( boundary.Name() );
  splitboundaryName.replace(splitboundaryName.find("BOUNDARY"),splitboundaryName.length(),"SPLIT_BOUNDARY");

  // attempt to create a splitboundary
  bool succeeded(false);
  pair<typename map<string, csmp::SplitBoundary<dim> >::iterator, bool>
    it = splitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,
                                                                                                splitboundaryComplex->Database() ) ) );
  if ( it.second ) {
      cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom:";
      cout <<" creating SplitBoundary from 'Boundary' "<< boundary.Name() << endl;
      splitboundaryComplex->UpdateIndices();
      succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );
      // SKM FIX - 
      if ( succeeded ) {
           SplitNodes( boundary, (*it.first).second ); 
           // TODO: no need to identify perimeter again because it exists already in Boundary
           (*it.first).second.IdentifyPerimeter();
        }
    }
  else
    throw csmp::Exception( WARNING,
                           "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom",
                           splitboundaryName.c_str(),
                           "boundary already exists. Nothing was done." );

  splitboundaryComplex->UpdateIndices();

  // removes boundary also deleting its elements
  splitboundaryComplex->RemoveBoundary( boundary, true );

  // update indexes
  splitboundaryComplex->UpdateIndices();

// TODO: put this repeated code into a private method that can be called separately
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity(false);
    rit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity(false);
    rit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
    bit->second.CreateNodePointerVector();
    bit->second.EstablishNeighborConnectivity(false);
    bit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
    sbit->second.CreateNodePointerVector();
    sbit->second.EstablishNeighborConnectivity(false);
    sbit->second.IdentifyPerimeter();
  }

  if ( succeeded ) {
        cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::CreateSplitBoundaryFrom: created splitboundary: '";
        cout << splitboundaryName <<"' successfully.\n\n";
    }

  return make_pair( splitboundaryName, false );

} // end CreateSplitBoundaryFrom( Boundary )




/**
Creates SplitBoundary via the creation of a boundary between regions.
The original boundary gets removed.

@author SKM
@date 25/1/20

 */
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
pair<string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const char* region1, const char* region2, bool createRegionBetween )
 {
   SPLITBOUNDARY_COMPLEX<dim>*  splitboundaryComplex(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
   
   pair<string,bool> result = splitboundaryComplex->InsertBoundary( region1, region2, createRegionBetween );
   
   return CreateSplitBoundaryFrom( splitboundaryComplex->Boundary(result.first) );

 } // end InsertSplitBoundary





/**
    Creates lower dimensional Region bisecting the SplitBoundary with its name followed by the string "_REGION".
    
    @author SKM
    @date 25/1/2020
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
std::pair<std::string,bool>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::InsertRegionIntoSplitBoundary( const char* split_boundary )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( !ContainsSplitBoundary( split_boundary ) ) {
         csmp_error.notice( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertRegionIntoSplitBoundary",
                            split_boundary, "Does not exist; nothing was done." );
         return make_pair("no Region created",false);
       }
     
     SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
     csmp::SplitBoundary<dim>&    splitBoundary( model->SplitBoundary(split_boundary) );
     
     // creating elements in the InterFace objects of the SplitBoundary with the node-numbering order matching that of the corresponding INNER parent element face
     const LocalVariables             element_props           = model->Database().LocalVariablesAt( ELEMENT );
     const LocalVariables             node_props              = model->Database().LocalVariablesAt( NODE );
     const IntegrationPointVariables  integration_point_props = model->Database().IntegrationPointVariablesAt( ELEMENT );
     size_t                           elmt_idx = model->Mesh().Elements(), node_idx = model->Mesh().Nodes();
     vector<size_t>                   element_numbers;
     const bool                       with_finite_volumes = (model->FV_Manager() == nullptr) ? false : true;
     vector<Element<dim>*>            new_elmts;
     
     element_numbers.reserve( splitBoundary.Elements() );
     new_elmts.reserve( splitBoundary.Elements() );
     
     for ( typename vector<InterFace<dim>*>::iterator 
           it=splitBoundary.ElementsBegin(); it!=splitBoundary.ElementsEnd(); ++it ) 
       {
           // construct the new element 
           element_numbers.push_back( elmt_idx );
           csmp::Element<dim>* new_elmt  = model->Mesh().Add( Element<dim>( elmt_idx++, (*it)->FE(), 
                                                                            element_props, integration_point_props, INTERNAL ) );
           if ( with_finite_volumes ) new_elmt->AssignFiniteVolume( (*it)->FV() );
           
           // assign the new element to the intervening element pointer of the InterFace element in the Splitboundary
           (*it)->Assign( new_elmt );
           new_elmts.push_back( new_elmt ); 
           
           // create new nodes and assign them to new element 
           // ASSUMING that the nodes are collocated
           for ( size_t n=0U; n<(*it)->Nodes()/2; ++n ) {
                csmp::Node<dim>* nptr = model->Mesh().Add( Node<dim>( node_idx++, (*it)->N(n,INSIDE)->Coordinate(), node_props, (*it)->N(n,INSIDE)->AtBoundary() ) );
                // assignment
                assert( nptr != nullptr );
                new_elmt->Assign( n, nptr );
                nptr->ResizeParentStorage( nptr->Parents() + 1U );
                nptr->Assign( nptr->Parents(), new_elmt );
                // TODO: only one node per patch needed 
             }
       }
       
     // update the neighbor connectivity of new elements
     // TODO: potential manifolds have to be disambiguated 
     establishNeighborConnectivity( new_elmts, false, false ); 

     // for each contiguous patch of the new region supply a pointer any of its elements into mesh manager 
     sort( new_elmts.begin(), new_elmts.end() );
     set<Element<dim>*>    contiguous_subset;
     vector<Element<dim>*> leftovers;
     
     while ( !new_elmts.empty() ) 
       {
          // finding contiguous element patch and inserting its first element into the root element vector
          floodFill( (*new_elmts.begin()), contiguous_subset );
          model->Mesh().SetRootElement( (*contiguous_subset.begin()) );
          model->Mesh().SetRootNode( (*contiguous_subset.begin())->N(0) );
          // if the split boundary is already contiguous
          if ( contiguous_subset.size() == new_elmts.size() ) break;
          // removing the pointers to the recovered elements from 'elmt_vec_to_establish_nbor_connectivity'
          leftovers.reserve( new_elmts.size() - contiguous_subset.size() );
          for ( typename vector<Element<dim>*>::const_iterator 
                it=new_elmts.begin(); it!=new_elmts.end(); ++it )
            // copy remaining elements into the leftover vector   
            if ( contiguous_subset.count( (*it) ) == 0 )
              leftovers.push_back( (*it) ); 

          // assigning result vector to repeat operation                
          new_elmts = leftovers;
          // contiguous_subset.clear(); - done in floodfill
          leftovers.clear();
       }

     // 4. constuct the new unique region with the interface elements in the model
     const string region_name( string(splitBoundary.Name()) + "_REGION" );
     const bool   unique_map(true);
     model->FormRegionFrom( region_name.c_str(), element_numbers, unique_map );
     
     return make_pair( region_name, true ); 
     
 } // InsertRegionIntoSplitBoundary
 








/** 
     Surt's code to efficiently remove a single element from a sorted vector, without preserving sorted order
     
     https://stackoverflow.com/questions/26719144/how-to-erase-a-value-efficiently-from-a-sorted-vector/26720032
 */
//inline void erase_v4(std::vector<int> &vec, int value)
template<size_t dim>
void eraseElementPointerFromVector( vector<csmp::Element<dim>*>& vec, const Element<dim>* eptr )
 {
    // get the range in 2*log2(N), N=vec.size()
    auto bounds = std::equal_range(vec.begin(), vec.end(), eptr );

    // calculate the index of the first to be deleted O(1)
    auto last = vec.end() - std::distance(bounds.first, bounds.second);

    // swap the 2 ranges O(equals) , equal = std::distance(bounds.first, bounds.last)
    std::swap_ranges(bounds.first, bounds.second, last);

    // erase the victims O(equals)
    vec.erase(last, vec.end());
}

template void eraseElementPointerFromVector( vector<csmp::Element<3U>*>&, const Element<3U>* );



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

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SingleRegionFromAllSplitBoundaries( const char* region_name )
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
      for (auto eit = (*it).second.ElementsBegin(); eit != (*it).second.ElementsEnd(); eit++) {
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

  for ( auto ifnode : ifnodes )
    {
      // TODO: check interface to determine whether the nodes are indeed collocated
      Node<dim>* new_node = mesh.Duplicate( *ifnode );
      if (new_node) {
          new_node->Idx(node_idx++);
          // storing pointers to the interface nodes as keys to retrieve the new nodes
          node_pairs.insert(make_pair(ifnode, new_node));
        }
    }

  const LocalVariables&             ifvars( splitboundaryComplex->Database().LocalVariablesAt( ELEMENT ) );
  const IntegrationPointVariables&  if_ip_vars( splitboundaryComplex->Database().IntegrationPointVariablesAt( ELEMENT_INTEGRATION_POINT ) );

  // creating the new lower-dimensional elements between the sides of collected InterFace objects and assigning them to middle-element pointer
  vector<size_t>  element_numbers;
  element_numbers.reserve( ifelmts.size() );

  for ( auto ifelmt : ifelmts )
    {
      CSMP_FEM_TYPE csmp_elmt = ifelmt->FE_Type();
      Element<dim>* new_elmt  = mesh.Add( Element<dim>( fem_mgr.E(csmp_elmt), 
                                                        const_cast<FiniteVolumeStencil<dim>*>(fvm_mgr->Stencil(csmp_elmt)), 
                                                        ifvars, if_ip_vars ) );
      new_elmt->Idx(elmt_idx++);
      element_numbers.push_back( new_elmt->Idx() );
      // assign the new element to the intervening element pointer of the InterFace element in the Splitboundary
      ifelmt->Assign( new_elmt ); 
   }

  // 3. Construct the connections between the new nodes, new elements and interfaces.
  //    Assign a root node and a root element for new nodes elements respectively
  vector<Element<dim>*>  elmt_vec_to_establish_nbor_connectivity;
  elmt_vec_to_establish_nbor_connectivity.reserve( ifelmts.size() );
  
  for ( auto ifelmt : ifelmts )
    {
      Element<dim>* const new_elmt = ifelmt->InterveningElement();

      for (size_t i = 0; i < ifelmt->Nodes() / 2; i++) {
          // finding the new node using the matching node of the inside of the interface as search key
          Node<dim>* new_node = node_pairs[ ifelmt->N(i, INSIDE) ];
          assert( new_node != nullptr );
          new_elmt->Assign(i, new_node);
          new_node->ResizeParentStorage(new_node->Parents() + 1);
          new_node->Assign(new_node->Parents(), new_elmt);
        }

      elmt_vec_to_establish_nbor_connectivity.push_back(new_elmt);
    }

  // update the neighbor connectivity of new region
  // TODO: potential manifolds have to be disambiguated 
  establishNeighborConnectivity( elmt_vec_to_establish_nbor_connectivity, false, false ); 

  // for each contiguous patch of the new region supply a pointer any of its elements into mesh manager 
  sort( elmt_vec_to_establish_nbor_connectivity.begin(), elmt_vec_to_establish_nbor_connectivity.end() );
  set<Element<dim>*>    contiguous_subset;
  vector<Element<dim>*> leftovers;
  
  while ( !elmt_vec_to_establish_nbor_connectivity.empty() ) 
    {
       // finding contiguous element patch and inserting its first element into the root element vector
       floodFill( (*elmt_vec_to_establish_nbor_connectivity.begin()), contiguous_subset );
       mesh.SetRootElement( (*contiguous_subset.begin()) );
       mesh.SetRootNode( (*contiguous_subset.begin())->N(0) );
       // if the split boundary is already contiguous
       if ( contiguous_subset.size() == elmt_vec_to_establish_nbor_connectivity.size() ) break;
       // removing the pointers to the recovered elements from 'elmt_vec_to_establish_nbor_connectivity'
       leftovers.reserve( elmt_vec_to_establish_nbor_connectivity.size() - contiguous_subset.size() );
       for ( typename vector<Element<dim>*>::const_iterator 
             it=elmt_vec_to_establish_nbor_connectivity.begin(); it!=elmt_vec_to_establish_nbor_connectivity.end(); ++it )
         // copy remaining elements into the leftover vector   
         if ( contiguous_subset.count( (*it) ) == 0 )
           leftovers.push_back( (*it) ); 

       // assigning result vector to repeat operation                
       elmt_vec_to_establish_nbor_connectivity = leftovers;
       // contiguous_subset.clear(); - done in floodfill
       leftovers.clear();
    }

  // 4. constuct a new unique region with the interface elements in the model
  const bool unique_map(true);
  splitboundaryComplex->FormRegionFrom( region_name, element_numbers, unique_map );

  return true;

} // end SingleRegionFromAllSplitBoundaries







/**
    In all SplitBoundaries of the model, lower dimensional Element regions are inserted.

       As many lower-dimensional element regions are created as there SplitBoundaries in the model
       @author SKM
       @date 26/1/2020
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
set<string>  SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RegionsFromSplitBoundaries()
{
  ErrorHandler&                csmp_error(ErrorHandler::Instance());
  SPLITBOUNDARY_COMPLEX<dim>*  splitboundaryComplex(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));
  set<string>                  new_regions;
  
  for ( typename std::map<std::string, csmp::SplitBoundary<dim> >::const_iterator
        it = splitboundaryComplex->SplitBoundariesBegin(); it != splitboundaryComplex->SplitBoundariesEnd(); ++it )
    {
       pair<string,bool>  result = splitboundaryComplex->InsertRegionIntoSplitBoundary( (*it).first.c_str() );
       if ( result.second == false )
         csmp_error.notice( ERROR, "SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RegionsFromSplitBoundaries", 
                            (*it).first, "unable to create Region from this boundary" );
    }

  return new_regions;

} // end RegionsFromSplitBoundaries




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




/**
    Replace all old nodes with duplicated ones and updates connectivity
    @author JC 1/7/2019
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::SplitNodes(  const Boundary<dim>& boundary, 
                                                                      csmp::SplitBoundary<dim>& splitboundary  )
{
  SPLITBOUNDARY_COMPLEX<dim>*  model(static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this));

  // Backup pointers to the nodes of the boundary that will be duplicated into a new set
  std::set<Node<dim>*> nodesToDuplicate;
  for ( typename vector<Node<dim>*>::const_iterator 
        nit( boundary.NodesBegin() ); nit != boundary.NodesEnd(); ++nit )
    nodesToDuplicate.insert( (*nit) );
  
  std::set<Node<dim>*>             outsideElementNodes, insideElementNodes;
  std::map<Node<dim>*,Node<dim>*>  manyfoldNodes;
  
  for ( typename std::vector<InterFace<dim>*>::const_iterator 
        ifit( splitboundary.ElementsBegin() ); ifit != splitboundary.ElementsEnd(); ++ifit )
  {
    Element<dim>* eit = (*ifit)->OuterParent();    
    for ( size_t en( 0 ); en < eit->Nodes(); ++en )
    {
      bool found( false );
      for ( size_t ifn( 0 ); ifn < (*ifit)->Nodes(); ++ifn )
      {
        if ( (*ifit)->N( ifn )->Idx() == eit->N( en )->Idx() )
          found = true;
      }
      if( found )
        outsideElementNodes.insert( eit->N( en ) );
    }
  }

  for ( typename std::vector<InterFace<dim>*>::const_iterator 
        ifit( splitboundary.ElementsBegin() ); ifit != splitboundary.ElementsEnd(); ++ifit )
  {
    Element<dim>* eit = (*ifit)->InnerParent();    
    for ( size_t en( 0 ); en < eit->Nodes(); ++en )
    {
      bool found( false );
      for ( size_t ifn( 0 ); ifn < (*ifit)->Nodes(); ++ifn )
      {
        if ( (*ifit)->N( ifn )->Idx() == eit->N( en )->Idx() )
          found = true;
      }
      if ( found )
      insideElementNodes.insert( eit->N( en ) );
    }
  }

  // prompting mesh manager to create new nodes
  for ( auto oen : outsideElementNodes ) {
      Node<dim>* duplicatedNode = model->Mesh().Add( Node<dim>(*oen) );
      //                          =====================================
      duplicatedNode->Idx( model->Mesh().Nodes() );
      manyfoldNodes.insert( make_pair( oen, duplicatedNode ) );
    }

  Region<dim>* outer_region = nullptr;
  for ( typename std::vector<InterFace<dim>*>::const_iterator 
        ifit( splitboundary.ElementsBegin() ); ifit != splitboundary.ElementsEnd(); ++ifit )
  {
    Element<dim>* oeit = (*ifit)->OuterParent();
    for ( typename std::map<std::string, csmp::Region<dim> >::iterator
          it = model->UniqueRegionsBegin(); it != model->UniqueRegionsEnd(); ++it ) {
      if ( (*it).second.Contains( oeit ) ) {
        outer_region = &(*it).second;
        break;
      }
    }
  }
  
  for ( typename std::vector<InterFace<dim>*>::const_iterator 
        ifit( splitboundary.ElementsBegin() ); ifit != splitboundary.ElementsEnd(); ++ifit )
  {
    Element<dim>* oeit = (*ifit)->OuterParent();
    for ( size_t en( 0 ); en < oeit->Nodes(); ++en ) {
      if ( manyfoldNodes.find( oeit->N( en ) ) != manyfoldNodes.end() ) {
        oeit->Assign( en, manyfoldNodes[oeit->N( en )] );
      }
    }

  // assigning manifold nodes to elements of the model 
  Region<dim>&  mref( model->Region( "Model" ) );
  for ( typename vector<Element<dim>*>::iterator eit( mref.ElementsBegin() ); eit != mref.ElementsEnd(); ++eit )
    for ( size_t en( 0 ); en < (*eit)->Nodes(); ++en )
      if ( outer_region->Contains( (*eit) ) )
        if ( manyfoldNodes.find( (*eit)->N( en ) ) != manyfoldNodes.end() )
          (*eit)->Assign( en, manyfoldNodes[(*eit)->N( en )] );
    }
  
} // end SplitNodes




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

   // reestablish containers
   // unique regions
   for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
     rit->second.CreateNodePointerVector();
     rit->second.EstablishNeighborConnectivity(false);
     rit->second.IdentifyPerimeter();
   }
   // non-unique regions
   for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
     rit->second.CreateNodePointerVector();
     rit->second.EstablishNeighborConnectivity(false);
     rit->second.IdentifyPerimeter();
   }
   for ( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
     bit->second.CreateNodePointerVector();
     bit->second.EstablishNeighborConnectivity(false);
     bit->second.IdentifyPerimeter();
   }
   for ( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
     sbit->second.CreateNodePointerVector();
     sbit->second.EstablishNeighborConnectivity(false);
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
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
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

// TODO: put this repeated code into a private method that can be called separately
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity(false);
    rit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
    rit->second.CreateNodePointerVector();
    rit->second.EstablishNeighborConnectivity(false);
    rit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
    bit->second.CreateNodePointerVector();
    bit->second.EstablishNeighborConnectivity(false);
    bit->second.IdentifyPerimeter();
  }
  for ( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
    sbit->second.CreateNodePointerVector();
    sbit->second.EstablishNeighborConnectivity(false);
    sbit->second.IdentifyPerimeter();
  }

  if ( succeeded ) {
       cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary: created split boundary '";
       cout << boundaryName <<"' between " << group1 << " and " << group2 << std::endl;
    }
  else
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

  return true;

} // end InsertSplitBoundary



DEPRECATED END */
