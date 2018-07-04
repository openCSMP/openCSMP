#include "SplitBoundaryInterface.h"
#include "Boundary.h"
#include "Model.h"

namespace csmp
{

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundaryInterface()
{}

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundaryInterface( const SplitBoundaryInterface& bd )
:interFaceSplitBoundaryMap_( bd.interFaceSplitBoundaryMap_ )
{

}

/// Deletes InterFaces
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::~SplitBoundaryInterface()
{
}

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
void SplitBoundaryInterface<dim, SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary, bool deleteElements )
  {

    SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>* >(this) );

    if( deleteElements )
    {
        // Faces removed from all existing boundaries and erased from MeshManager
        std::vector<InterFace<dim>* > elementsToDelete( splitboundary.SimplexVector().begin(), splitboundary.SimplexVector().end() );

        if( !elementsToDelete.empty() )
        {
            // Detach InterFaces from Neighbors
            splitboundary.DetachElementsFromNeighbors();

            // remove redundant InterFaces from existing SplitBoundaries
            for( splitBoundaryIterator
                 sbit  = interFaceSplitBoundaryMap_.begin();
                 sbit != interFaceSplitBoundaryMap_.end();){

                if( removeVectorElements( (*sbit).second.SimplexVector(), elementsToDelete ) > 0 )
                {
                    if( !(*sbit).second.SimplexVector().empty() )
                    {
                        (*sbit).second.CreateNodePointerVector();
                        (*sbit).second.IdentifyPerimeter();
                        ++sbit;
                    }
                    else
                        interFaceSplitBoundaryMap_.erase( sbit++ );
                }
                else
                {
                    // we need to reestablish the perimeter because
                    // the order of SimplexVector() was changed by removeVectorElements ()
                    (*sbit).second.IdentifyPerimeter( );
                    ++sbit;
                }
            }

            // update all indices
            splitboundaryComplex->UpdateIndices();

            return;
        }

    } // deleting elements

    // finding the boundary
    typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator iterSplitBoundary( interFaceSplitBoundaryMap_.end() );
    for( splitBoundaryIterator it = interFaceSplitBoundaryMap_.begin(); it != interFaceSplitBoundaryMap_.end(); ++it )
        if( &it->second == &splitboundary )
            iterSplitBoundary = it;

    if( iterSplitBoundary == interFaceSplitBoundaryMap_.end() )
        throw csmp::Exception( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::RemoveSplitBoundary", "splitboundary does not exist" );

    // if the boundary was found in the list, it is erased
    interFaceSplitBoundaryMap_.erase( iterSplitBoundary );
  }

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin()
  { return interFaceSplitBoundaryMap_.begin(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd()
  { return interFaceSplitBoundaryMap_.end(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string,csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundariesBegin() const
  { return interFaceSplitBoundaryMap_.begin(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
typename std::map<std::string,csmp::SplitBoundary<dim> >::const_iterator  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundariesEnd() const
  { return interFaceSplitBoundaryMap_.end(); }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
size_t  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundaries() const
  { return interFaceSplitBoundaryMap_.size(); }















template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::ContainsSplitBoundary( const std::string& bname ) const
 {
    std::string  b_normal(bname);
    if ( interFaceSplitBoundaryMap_.find(b_normal) != interFaceSplitBoundaryMap_.end() ) return true;
    return false;
 }


template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim,BOUNDARY_COMPLEX>::IsEligibleSplitBoundaryRegionName( const std::string& spbNname ) const
  {
    // TAG: SPLITBOUNDARY should be at the beginning of the region name
    // DON'T confuse with BOUNDARY
    std::size_t found_position = spbNname.find( "SPLITBOUNDARY" );
    if ( found_position == 0 && found_position!=std::string::npos )
        return true;
    return false;
  }


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
SplitBoundary<dim>&  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundary( const std::string& spbName )
  { 
    splitBoundaryIterator  sbit = interFaceSplitBoundaryMap_.find( spbName );
    if ( sbit != interFaceSplitBoundaryMap_.end() )
      return (*sbit).second;
    else
    {
        std::string errMsg("SplitBoundary does not exist!");
        errMsg.append( " (" + spbName + ")" );
        throw csmp::Exception( ERROR,
                             "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundary",
                             errMsg.c_str() );
    }
  } 


template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
const SplitBoundary<dim>&  SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundary( const std::string& spbName ) const
 {
   splitBoundaryConstIterator  sbit = interFaceSplitBoundaryMap_.find( spbName );
   if ( sbit != interFaceSplitBoundaryMap_.end() )
     return (*sbit).second;
   else
   {
     std::string errMsg("SplitBoundary does not exist!");
     errMsg.append( " (" + spbName + ")" );
     throw csmp::Exception( ERROR,
                            "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::SplitBoundary",
                            errMsg.c_str() );
   }
 } 














// -----------------------------------------------
// Binary input/output
// -----------------------------------------------

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::OutputSplitBoundariesToBinary( const char* file_name ) const
  {
  std::string bin_file(file_name);
  std::string heading("Binary splitboundary list for CSMP model, file '");
  heading += bin_file;
  heading +="'.";

  FILE*  fp(0);
  if ( (fp=fopen( bin_file.c_str(), "wb")) == NULL ) {
    throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary:", "Binary file could not be created." );
    return false;
    }
  // writing the file header
  skm_C_fwrite( fp, heading.c_str() );

  // writing number of splitboundaries
  size_t records(this->SplitBoundaries());
  fwrite( (void*) &records, sizeof(size_t), 1, fp );

  for ( typename std::map<std::string,csmp::SplitBoundary<dim> >::const_iterator
        bit( this->SplitBoundariesBegin() ); bit != this->SplitBoundariesEnd(); bit++ )
    {
    // writing name of splitboundary
    skm_C_fwrite( fp, bit->first.c_str() );
    // splitboundary
    if( !bit->second.Out(fp) )
      throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary:", "SplitBoundary could not be stored." );
    }

  // cleaning up
  fclose( fp );
  std::cout <<"\nSplitBoundaryInterface<dim>::OutputSplitBoundariesToBinary: split boundaries have been successfully written to: '";
  std::cout << bin_file <<"'"<< std::endl;
  return true;
  } // OutputSplitBoundariesToBinary




template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InputSplitBoundariesFromBinary( const char* file_name )
  {
  SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
  const csmp::Region<dim>&  rref( splitboundaryComplex->Region("Model") );

  std::string bin_file(file_name);
  FILE*  fp(0);
  if ( (fp=fopen( bin_file.c_str(), "rb")) == NULL ) {
    throw csmp::Exception( ERROR, "BoundaryInterface<dim>::InputSplitBoundariesFromBinary:", "Binary file could not be opened." );
    return false;
    }

  // reading the file header
  char csCache[255];
  skm_C_fread( fp, csCache );

  // reading number of splitboundaries
  size_t records(0);
  fread( (void*) &records, sizeof(size_t), 1, fp );

  std::cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: reading " << csCache << " containing "
            << records << " splitboundaries\n" <<  std::endl;

  for ( size_t i(0); i < records; ++i )
    {
    std::string bName;
    // reading name of splitboundary
    skm_C_fread( fp, csCache );
    bName = csCache;
    // inserting splitboundary if not existing yet
    std::pair<typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator,bool>
      bit = interFaceSplitBoundaryMap_.insert( std::make_pair( bName, csmp::SplitBoundary<dim>( bName, splitboundaryComplex->Database() ) ) );
    // splitboundary
    if( !bit.second )
      return false;
    if( !bit.first->second.In( splitboundaryComplex->Mesh(), splitboundaryComplex->FE_Manager(), rref, fp) )
      return false;
    std::cout << "\nSplitBoundaryInterface<dim>::InputSplitBoundariesFromBinary: read boundary " << bName << " successfully.\n";
    }
  fclose(fp);
  return true;
  } // InsertSplitBoundaryFromBinary





















/**
Method forms a SplitBoundary between the two supplied regions. This will involve the creation
and connection of InterFaces.

@attention:  SplitBoundaryInterface cannot be created in 1D models or between regions which contain
one-dimensional elements.

This method will only work for Regions which are not overlapping.

This method makes no sense for the region 'Model' as it encompasses all unique regions.

Option deleteRegionAndItsElements in fact provides a choice for keeping fracture\fault or any other low dimensional Region
or use such Region just to introduce the interface around it and delete it afterwards
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
bool SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const std::string& group, bool deleteRegionAndItsElements )
 {
    bool succeeded( false );
    SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( group == "Model" ) {
         csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Region 'Model' not eligible for InsertSplitBoundary.");
         return false;
      }

    if ( !splitboundaryComplex->IsUnique(group.c_str()) )
         csmp_error.notice( WARNING, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                           "This method is intended for the creation of splitboundaries around unique Regions");

    std::string  splitboundaryName( group );
    splitboundaryComplex->InsertBoundary( group.c_str(), IRREGULAR, false /* do not delete region*/ );
    Boundary<dim>& boundary( splitboundaryComplex->Boundary( splitboundaryName ) );

    // attempt to create a regular (InterFace-based) splitboundary
     std::pair<typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator,bool>
         it = interFaceSplitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName, splitboundaryComplex->Database() ) ) );
     if ( it.second )
       {
         std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary around " << group << std::endl;
         splitboundaryComplex->UpdateIndices();
         succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );
         if( deleteRegionAndItsElements ) {
              // SKM FIX
              splitboundaryComplex->RemoveRegion( group.c_str(), false );
              splitboundaryComplex->MoveToNonUniqueRegions( group.c_str() );
           }
       }
     else
         throw csmp::Exception( INFO, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                                      splitboundaryName.c_str(),
                                      "boundary already exists. Nothing was done.");

     // remove temporal boundary
     splitboundaryComplex->RemoveBoundary( boundary );


     // reestablish containers
     // unique regions
     for( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
         rit->second.CreateNodePointerVector();
         rit->second.EstablishNeighborConnectivity();
         rit->second.IdentifyPerimeter(  );
       }
     // non-unique regions
     for( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
         rit->second.CreateNodePointerVector();
         rit->second.EstablishNeighborConnectivity();
         rit->second.IdentifyPerimeter(  );
       }
     for( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
         bit->second.CreateNodePointerVector();
         bit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
         bit->second.IdentifyPerimeter(  );
       }
     for( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
         sbit->second.CreateNodePointerVector();
         sbit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
         sbit->second.IdentifyPerimeter(  );
       }

    // update indexes
    splitboundaryComplex->UpdateIndices();

    if( succeeded )
        std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary created splitboundary around " << group << std::endl;
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
bool SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary( const std::string& group1, const std::string& group2, bool createRegionBetween )
 {
    SPLITBOUNDARY_COMPLEX<dim>* splitboundaryComplex( static_cast<SPLITBOUNDARY_COMPLEX<dim>*>(this) );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( group1 == "Model" or group2 == "Model" ) {
         csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Region 'Model' not eligible for InsertSplitBoundary.");
         return false;
      }
    if ( group1 == group2 ) {
         csmp_error.notice( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Provided Regions are identical.");
         return false;
      }


    if ( !splitboundaryComplex->IsUnique(group1.c_str()) or !splitboundaryComplex->IsUnique(group2.c_str()) ) {
         csmp_error.notice( WARNING,
                            "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                            "This method is intended for the creation of splitboundaries between unique Regions");

         // checking for a potential overlap of the regions, if the regions are non-unique
         if ( splitboundaryComplex->RegionIntersection( group1.c_str(),  group2.c_str(), "groupintersection" ) ) {
               csmp_error.notice( ERROR,
                                  "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                                  "one of the supplied regions is not unique and they overlap",
                                  "It was therefore impossible to insert a boundary");
               splitboundaryComplex->RemoveRegion( "groupintersection", false );
               return false;
         }
      }

    std::string  splitboundaryName( std::string(group1) + std::string("_") + std::string(group2) );

    splitboundaryComplex->InsertBoundary( group1.c_str(), group2.c_str(), createRegionBetween );
    Boundary<dim>& boundary( splitboundaryComplex->Boundary( splitboundaryName ) );

    const csmp::Region<dim>&  gref1( splitboundaryComplex->Region(group1.c_str()) );
    const csmp::Region<dim>&  gref2( splitboundaryComplex->Region(group2.c_str()) );

    // checking whether the two regions share some nodes (these will mark their common boundary)
    const size_t  shared_nodes(sharedNodes( gref1, gref2 ));
    bool succeeded( false );

    if ( shared_nodes == 0U )
      csmp_error.notice( WARNING,
                         "Model<dim,SPLITBOUNDARY_COMPLEX>::InsertBoundary",
                        "The regions of interest do not share any nodes; trying to create a split boundary");


    // attempt to create a regular (InterFace-based) splitboundary
    else {

         std::pair<typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator,bool>
             it = interFaceSplitBoundaryMap_.insert( std::make_pair( splitboundaryName, csmp::SplitBoundary<dim>( splitboundaryName,  splitboundaryComplex->Database() ) ) );
         if ( it.second )
           {
             std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary creating splitboundary between " << group1 << " and " << group2 << std::endl;
             splitboundaryComplex->UpdateIndices();
             if( createRegionBetween )
                 splitboundaryComplex->RegionBetween( group1.c_str(), group2.c_str(), std::string( group1 + std::string("_") + group2 ).c_str() );
             succeeded = (*it.first).second.CreateFrom( *splitboundaryComplex, boundary );
           }
         else
             throw csmp::Exception( INFO,
                                    "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary",
                                    splitboundaryName.c_str(),
                                    "boundary already exists. Nothing was done.");
      }

    // remove temporal boundary
    splitboundaryComplex->RemoveBoundary( boundary );

    // reestablish containers
    for( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->UniqueRegionsBegin(); rit != splitboundaryComplex->UniqueRegionsEnd(); ++rit ) {
        rit->second.CreateNodePointerVector();
        rit->second.EstablishNeighborConnectivity();
        rit->second.IdentifyPerimeter();
      }
    // SKM FIX also the non-unique regions need to be rebuilt
    for( typename SPLITBOUNDARY_COMPLEX<dim>::regionIterator rit = splitboundaryComplex->RegionsBegin(); rit != splitboundaryComplex->RegionsEnd(); ++rit ) {
        rit->second.CreateNodePointerVector();
        rit->second.EstablishNeighborConnectivity();
        rit->second.IdentifyPerimeter();
      }
    for( typename SPLITBOUNDARY_COMPLEX<dim>::boundaryIterator bit = splitboundaryComplex->BoundariesBegin(); bit != splitboundaryComplex->BoundariesEnd(); ++bit ) {
        bit->second.CreateNodePointerVector();
        bit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
        bit->second.IdentifyPerimeter();
      }
    for( typename SPLITBOUNDARY_COMPLEX<dim>::splitBoundaryIterator sbit = splitboundaryComplex->SplitBoundariesBegin(); sbit != splitboundaryComplex->SplitBoundariesEnd(); ++sbit ) {
        sbit->second.CreateNodePointerVector();
        sbit->second.EstablishNeighborConnectivity(); // TODO: unassign neighbors outside
        sbit->second.IdentifyPerimeter();
      }

    // update indexes
    splitboundaryComplex->UpdateIndices();

    if( succeeded )
        std::cout << "\nSplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary created splitboundary between " << group1 << " and " << group2 << std::endl;
    else
      throw csmp::Exception( ERROR, "SplitBoundaryInterface<dim,SPLITBOUNDARY_COMPLEX>::InsertSplitBoundary", "Splitting failed!" );

    return true;

 } // end InsertSplitBoundary


 template class SplitBoundaryInterface<1U,Model>;
 template class SplitBoundaryInterface<2U,Model>;
 template class SplitBoundaryInterface<3U,Model>;

} // csmp
