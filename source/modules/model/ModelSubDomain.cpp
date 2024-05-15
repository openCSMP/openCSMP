/*
 *  ModelSubDomain.cpp
 *
 *  Created by Stephan Matthai on 7/18/10.
 *  Copyright 2010 SKM. All rights reserved.
 *
 */
#include <type_traits>
#include "writeVariableIf.h"
#include "ModelSubDomain.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Box.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "FiniteElementManager.h"
#include "MeshManager.h"
#include "FiniteVolumeStencilManager.h"

#include "Visitor.h"
#include "PropertyDatabase.h"
#include "Interrelation.h"
#include "CopyReplaceVisitor.h"
#include "PropertyConstraints.h"

#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include "binaryReadWrite.h"

#include "MeshManagementUtilities.h"

//#define MODEL_SUBDOMAIN_DEBUG

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
ModelSubDomain<dim, CELL>::ModelSubDomain( const string& subdomain_name, const PropertyDatabase<dim>& pref )
	: pref_(pref),
	  subdomain_name_(subdomain_name),
    domain_idx_(++domain_count_)
 {
    if ( verbose_ ) cout <<"\nModelSubDomain(idx="<< domain_idx_ <<"): called custom constructor.\n";
 }


template<uint32_t dim, template<uint32_t> class CELL>
ModelSubDomain<dim,CELL>::ModelSubDomain( const ModelSubDomain& ed )
  : pref_(ed.pref_),
    cell_vec_(ed.cell_vec_),
    node_vec_(ed.node_vec_),
    first_bd_node_(ed.first_bd_node_),
    bd_face_vec_(ed.bd_face_vec_),
    subdomain_name_(ed.subdomain_name_),
    rebuilt_needed_(ed.rebuilt_needed_),
    domain_idx_(++domain_count_)
 {
    if ( verbose_ ) cout <<"\nModelSubDomain(idx="<< domain_idx_ <<"): called copy constructor.\n";
 }


/// move constructor; @attention remove verbose output after testing
template<uint32_t dim, template<uint32_t> class CELL>
ModelSubDomain<dim,CELL>::ModelSubDomain( ModelSubDomain&& ed )
 : pref_( std::move(ed.pref_) ),
   cell_vec_( std::move(ed.cell_vec_) ),
   node_vec_( std::move(ed.node_vec_) ),
   first_bd_node_( std::move(ed.first_bd_node_) ),
   bd_face_vec_( std::move(ed.bd_face_vec_) ),
   subdomain_name_( std::move(ed.subdomain_name_) ),
   rebuilt_needed_(std::move(ed.rebuilt_needed_) ),
   domain_idx_( std::move(ed.domain_idx_) ) // since argument object gets destroyed there is no incrementation of domain_idx_
 {
    domain_count_++; // needed because when destructor is called on 'ed' the object count will be decremented!
    if ( verbose_ ) cout <<"\nModelSubDomain(idx="<< domain_idx_ <<"): called move constructor.\n";
 }




template<uint32_t dim, template<uint32_t> class CELL>
ModelSubDomain<dim,CELL>::~ModelSubDomain()
 {
    if ( verbose_ ) cout <<"\nModelSubDomain(idx="<< domain_idx_ <<"): called destructor.\n";
    domain_count_--;
 }



template<uint32_t dim, template<uint32_t> class CELL>
ModelSubDomain<dim,CELL>&  ModelSubDomain<dim,CELL>::operator=( const ModelSubDomain& ed )
 {
     if ( &ed != this ) {
          cell_vec_       = ed.cell_vec_;
          node_vec_       = ed.node_vec_;
          first_bd_node_  = ed.first_bd_node_;
          //domain_idx_     = ed.domain_idx_; - keep the domain index unique!
          bd_face_vec_    = ed.bd_face_vec_;
          subdomain_name_ = ed.subdomain_name_;
          rebuilt_needed_ = ed.rebuilt_needed_;
          if ( verbose_ ) cout <<"\nModelSubDomain(idx="<< domain_idx_ <<"): called assignment operator.\n";
       }
     return *this;
 }


template<uint32_t dim, template<uint32_t> class CELL>
ModelSubDomain<dim,CELL>&  ModelSubDomain<dim,CELL>::operator=( ModelSubDomain&& ed )
 {
     if ( &ed != this ) {
         cell_vec_       = std::move( ed.cell_vec_ );
         node_vec_       = std::move( ed.node_vec_ );
         first_bd_node_  = std::move( ed.first_bd_node_ );
         domain_idx_     = std::move( ed.domain_idx_ );
         bd_face_vec_    = std::move( ed.bd_face_vec_ );
         subdomain_name_ = std::move( ed.subdomain_name_ );
         rebuilt_needed_ = std::move( ed.rebuilt_needed_ );
         if ( verbose_ ) cout <<"\nModelSubDomain(idx="<< domain_idx_ <<"): called move assignment operator.\n";
       }
     return *this;
 }


template<uint32_t dim, template<uint32_t> class CELL>
string  ModelSubDomain<dim,CELL>::Name() const
 {
    return subdomain_name_;
 }


template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::Name( const string& name )
 {
    subdomain_name_ = name;
 }


template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::ScheduleForRebuilt()
 {
    rebuilt_needed_ = true;
 }

template<uint32_t dim, template<uint32_t> class CELL>
bool  ModelSubDomain<dim,CELL>::NeedsRebuilt() const
 {
    return rebuilt_needed_;
 }



template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::Nodes() const
  {
     return node_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::InteriorNodes() const
  {
     return first_bd_node_;
  }

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::PerimeterNodes() const
  {
     return node_vec_.size() - InteriorNodes();
  }

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::Cells() const
  {
     return cell_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::InteriorCells() const
  {
     return cell_vec_.size() - bd_face_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::PerimeterCells() const
  {
     return bd_face_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
bool ModelSubDomain<dim,CELL>::Empty() const
  {
     return cell_vec_.empty();
  }



/**
    To loop over the faces of the perimeter cells which lie on the subdomain boundary.
    
    @param cell_idx marks the location of the cell in  bd_face_vec_  and has to be in the range of InteriorCells() and total number of Cells in the subdomain-1U.
    
    To loop over all the boundary faces of a subdomain, use the following code snippet:
    
    @code
    for ( size_t i{ subdomain.InteriorCells() }; i<subdomain.Cells(); i++ )
      for ( auto j{0U}; j<subdomain.PerimeterFaces(i); j++ ) {
           auto perim_face = subdomain.PerimeterFace(i,j);
           // get a normal to the cell face
           Point<dim> unrml = subdomain.E(i)->FE()->Point<dim> UnitNormalToFace(j);
           ...
        }
    @endcode
        
    @return returns how many faces of the target cell lie on the subdomain boundary
 
    @attention the cell index that is supplied as a method argument has to
    range between e = interior cells and cells-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  ModelSubDomain<dim,CELL>::PerimeterFaces( size_t cell_idx ) const
 {
    assert( cell_idx >= InteriorCells() );
    assert( cell_idx < cell_vec_.size() );
    return static_cast<uint32_t>( bd_face_vec_[cell_idx - InteriorCells()].size() );
 }


/**
    @return returns cells local cell face number (0..faces-1) for
    the n'th face that is on the subdomain boundary.

    @attention the cell index that is supplied as a method argument has to
    range between e = interior cells and cells-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  ModelSubDomain<dim,CELL>::PerimeterFace( size_t e, uint32_t face ) const
 {
    assert( e >= InteriorCells() );
    assert( e < cell_vec_.size() );
    assert( face < PerimeterFaces(e) );
    return static_cast<uint32_t>(bd_face_vec_[e-InteriorCells()][face]);
 }


template<uint32_t dim, template<uint32_t> class CELL>
csmp::Node<dim>*  ModelSubDomain<dim,CELL>::N( size_t nd ) const
 { assert( nd < node_vec_.size() ); return node_vec_[nd]; }


template<uint32_t dim, template<uint32_t> class CELL>
CELL<dim>*  ModelSubDomain<dim,CELL>::E( size_t e ) const
 { assert( e < cell_vec_.size() ); return cell_vec_[e]; }


template<uint32_t dim, template<uint32_t> class CELL>
const typename std::vector<CELL<dim>*>&  ModelSubDomain<dim,CELL>::CellVector() const
 { return cell_vec_; }

template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<CELL<dim>*>&  ModelSubDomain<dim,CELL>::CellVector()
 { return cell_vec_; }

template<uint32_t dim, template<uint32_t> class CELL>
const typename std::vector<Node<dim>*>&  ModelSubDomain<dim,CELL>::NodeVector() const
  { return node_vec_; }



template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::NodesBegin() const
 { return node_vec_.begin(); }

template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::NodesEnd() const
 { return node_vec_.end(); }

template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::PerimeterNodesBegin() const
 { return std::next( node_vec_.begin(), InteriorNodes() ); }

template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<CELL<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::CellsBegin() const
 { return cell_vec_.begin(); }

template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<CELL<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::CellsEnd() const
 { return cell_vec_.end(); }

template<uint32_t dim, template<uint32_t> class CELL>
typename std::vector<CELL<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::PerimeterCellsBegin() const
  { return std::next( cell_vec_.begin(), InteriorCells() ); }
   




/**
Tests whether the cells of the region are connected to each-other.

@attention This test cannot be performed if the region has cells of different spatial
dimensions since these are not interconnected. Therefore, this method returns false if
the region consists of cells from different spatial dimensions.

@note This method is based on the floofFill() algorithm implemented in CSMP.

*/
template<uint32_t dim, template<uint32_t> class CELL>
bool  ModelSubDomain<dim,CELL>::IsContiguous() const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  pair<int32_t, int32_t>  dimensionality = SpatialDimensions();
  if ( dimensionality.first > 1 ) {
      csmp_error.Note( WARNING, "ModelSubDomain<dim,CELL>::IsContiguous:",
                         "method can determine contiguity only for domains which consist only of same shape cells (line, surface or volume type); returned false." );
      return false;
    }

  set<CELL<dim>*> contiguous_elmts;

  floodFill( (*cell_vec_.begin()), contiguous_elmts );
  if ( contiguous_elmts.size() != cell_vec_.size() ) return false;

  return true;

} // end IsContiguous




/**
    Detects of how many spatial dimensions cell types are contained in model.
    
    @return method returns a pair: first value gives number of different spatial dimensions contained,
    second value returns the highest spatial dimension contained.

    @author SKM 1/11/2013

*/
template<uint32_t dim, template<uint32_t> class CELL>
pair<int32_t,int32_t>  ModelSubDomain<dim,CELL>::SpatialDimensions() const
 {
    bool with_volume_cells(false);
    bool with_surface_cells(false);
    bool with_line_cells(false);

    for( const auto& it : cell_vec_ ) {
         if      ( it->IsLine() )    with_line_cells = true;
         else if ( it->IsSurface() ) with_surface_cells = true;
         else if ( it->IsVolume() )  with_volume_cells = true;
      }

    int32_t counter(0);
    if ( with_volume_cells )  counter++;
    if ( with_surface_cells ) counter++;
    if ( with_line_cells )    counter++;
    int32_t highest_spatial_dim(LINE);
    if      ( with_volume_cells )  highest_spatial_dim = VOLUME;
    else if ( with_surface_cells ) highest_spatial_dim = SURFACE;

    return make_pair( counter, highest_spatial_dim );

 } // end ElementSpatialDimensions




template<uint32_t dim, template<uint32_t> class CELL>
pair<CELL_SHAPE,bool>  ModelSubDomain<dim,CELL>::SingleCellShapeDomain() const
 {
    const auto cell_dim = SpatialDimensions();
    const bool single_elmt_domain = ( cell_dim.first == 1 ) ? true : false;
    if ( !single_elmt_domain ) return make_pair( static_cast<CELL_SHAPE>(UNSPECIFIED), single_elmt_domain );
    
    // cell shape
    switch ( cell_dim.second ) {
        case 1:
          return make_pair( LINE, single_elmt_domain );
        case 2:
          return make_pair( SURFACE, single_elmt_domain );
        case 3:
          return make_pair( VOLUME, single_elmt_domain );
      }
    
    return make_pair( static_cast<CELL_SHAPE>(dim), false );
    
 } // end SingleCellShapeDomain













/**
    In order to use the region node/cell flags 'INTERIOR' or 'PERIMETER',
    boundary nodes and cells must be identified first.
    This is done every time when a new region is formed.

    @attention convention: cells are considered 'PERIMETER' cells
    of a Region only if they have at least one edge or face at the
    region boundary. This is the case irrespective of the dimensionality
    of such cells. In the extreme case of line cells inside a volume,
    they are flagged boundary, if they have one node on the perimeter of
    the volume region.
    Another way of looking at this is to consider all cells boundary
    cells, that share a face with the boundary.

    @attention convention: If a Region consists of cells of different
    dimensionality, the highest dimensional ones define the position of the boundary, i.e.
    all lower-dimensional cells that stick out of the region (and including all their
    nodes) are flagged 'PERIMETER'.
    Lower-dimensional cells inside a higher dimensional region are considered
    'INTERIOR'.

    @section application Application

    IdentifyPerimeter() is called as part of the region-forming process.
    The perimeter flagging is used to access boundaries selectively.

    @todo (1) Does not work for quadratic regions

   @section implementation Implementation

   1. Perimeter cells and nodes are found looking for the absence
      of cell neighbors, but only the highest dimensionality cells
      inside the Region

      - all the nodes of the highest-dimensionality cells are stored
        in a set.

   2. Lower-dimensional cells are found that do not share all of their nodes
      with the highest dimensional cells.

      - These cells nodes must be outside of the highest dimensional region.
      - These cells and their nodes are flagged 'PERIMETER'.

   3. All lower dimensional cells that share all of their nodes with the
      highest dimensional cells AND have at least nodes flagged 'PERIMETER'
      on one edge (a single node for a line cell)
      will be flagged 'PERIMETER' cells.

   It follows that all lower-dimensional mesh that 'sticks out' of a higher
   dimensional region has cells and nodes flagged as PERIMETER.

*/
template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::IdentifyPerimeter()
 {
    // 0. recreate Pfverts information for testing
    //EstablishNeighborConnectivity();
   
    // 1. putting the perimeter cells at the end of the elmt_vec and sorting
    //    interior and perimeter cell ranges subsequently
    // ----------------------------------------------------
    const size_t first_bd_elmt = PartitionCellVector();
    assert( first_bd_elmt <= cell_vec_.size() );

 } // end IdentifyPerimeter






/**
     Method for reconstruction of model  subdomains from binary file.
     
         @author  SKM
         @date 11/10/2019
*/
template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::BuildPerimeterFaceVector( int64_t interior_cells )
  {
     ErrorHandler&  csmp_error(ErrorHandler::Instance());

     // verification of suitable model state
     if ( cell_vec_.empty() )
       csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::BuildPerimeterFaceVector:",
                          Name(), "model subdomain: CELL vector not initialised yet.");
                              
     if ( interior_cells > static_cast<int64_t>(cell_vec_.size()) ) {
           cerr <<"\ninterior cells: "<< interior_cells;
           csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::BuildPerimeterFaceVector:",
                              Name(), "model subdomain: less CELLs in CELL vector than interior cells specified.");
       }

     if ( !is_sorted( cell_vec_.begin(), next(cell_vec_.begin(),interior_cells)) ) {
          cerr <<"\ninterior cells: "<< interior_cells;
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::BuildPerimeterFaceVector:",
                             Name(), "model subdomain: supplied CELL vector not sorted.");
       }
       
     // (re)setting the CELL face vector; note: size must already be correct else Contains(eptr) function used below will fail
     const int64_t n_elmts_with_perimeter_faces = static_cast<int64_t>(cell_vec_.size()) - interior_cells;
     assert( n_elmts_with_perimeter_faces > 0 );
     bd_face_vec_.resize( static_cast<size_t>(n_elmts_with_perimeter_faces) );

     const typename vector<CELL<dim>*>::const_iterator perimeterElementsBegin( next(cell_vec_.begin(),interior_cells) );
     const typename vector<CELL<dim>*>::const_iterator cellsEnd( cell_vec_.end() );
     vector<uint32_t>  boundary_faces;
     size_t            counter{0U};

     // for all faces of CELLs that are located on the subdomain boundary
     for ( auto it = perimeterElementsBegin; it != cellsEnd; ++it )
       {
           // line elements cannot have a face on the outside of the model
           if ( (*it)->IsLine() == false )
             {
                assert( (*it)->Faces() == (*it)->Neighbors() );
                const uint32_t faces( (*it)->Faces() );
                boundary_faces.reserve( faces );
                for ( uint32_t face{0U}; face<faces; ++face )
                  //   located on model boundary          or  neighbor is not contained in this subdomain
                  if ( (*it)->Neighbor( face ) == nullptr || !(Contains((*it)->Neighbor(face))) )
                    boundary_faces.push_back( face );
                // storing the boundary face vector for the current cell
                if ( boundary_faces.size() == faces ) {
                     if ( (dim == 2 && (*it)->IsSurface()) ||
                          (dim == 3 && (*it)->IsVolume()) ) {
                          cout <<"\n\tINFO, ModelSubDomain<dim,CELL>::BuildPerimeterFaceVector: subdomain '"<< Name();
                          cout <<"', cell: "<< (*it)->Idx() <<"("<< parseFiniteElementType((*it)->FE_Type()) <<")";
                          cout <<" is a stand-alone cell in this subdomain.";
                       }
                  }
                bd_face_vec_[counter] = boundary_faces;
                boundary_faces.clear();
             }
           counter++;
         }
      
  } // end BuildPerimeterFaceVector







/**
   Distinguishes 'interior' from 'perimeter' cells of the model subdomain by checking for each cell face whether
   this face is located at a model boundary or has a neighbor that does not belong to the current region (=cell range of ModelSubdomain).
   
   Method also creates the bd_face_vector   enlisting all faces of cells that are located on the region boundary.

  Remarks on binary search:

  - The sorted source range referenced must be valid; all pointers must be dereferenceable and,
    within the sequence, the last position must be reachable from the first by incrementation.

  - The sorted range must each be arranged as a precondition to the application of the binary_search
    algorithm in accordance with the same ordering as is to be used by the algorithm
    to sort the combined ranges.

    @attention any lower-dimensional cells and their nodes that stick outside of a higher
    dimensional region will be flagged as boundary.
    
    @todo due to the searching, this method is a speed bottleneck; can it be improved.
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  ModelSubDomain<dim,CELL>::PartitionCellVector()
 {
    ErrorHandler&  csmp_error(ErrorHandler::Instance());

    if ( cell_vec_.empty() )
      throw logic_error( (string("ModelSubDomain<dim>::PartitionCellVector: method called on empty region: ") + Name()).c_str() );

    sort( cell_vec_.begin(), cell_vec_.end() );

    // --------------------------------------------------------------------
    // 0. detecting whether this region contains lower-dimensional cells
    // --------------------------------------------------------------------
    // pair contains: 1) number of spatial cell dimensions in region, 2)  highest contained dimension
    const pair<int32_t,int32_t>  elmt_dim = SpatialDimensions();

    // -----------------------------------------------------------------
    // 1. distinguishing boundary from interior cells, same for nodes
    //   (at this point the cells and nodes are already known)
    // -----------------------------------------------------------------
    set<CELL<dim>*> interior_elmts, boundary_elmts;
    set<Node<dim>*>                 boundary_nodes;
    set<pair<CELL<dim>*,uint32_t> > boundary_faces;

    // 1.1 If all cells have the same spatial dimension
    // ---------------------------------------------------
    if ( elmt_dim.first == 1U )
      {
        for ( auto& eit : cell_vec_ )
          {
            // identifying the boundary faces and their nodes
            // (each face potentially has a neighbor cell)
            auto nbors_that_belong_to_group{ eit->Neighbors() };
            const auto n_faces{ eit->Faces() };
            for ( auto i{0U}; i<n_faces; i++ )
              // if the face is at a model boundary or has a neighbor that does not belong to the region
              if ( !eit->Neighbor(i) || !binary_search( cell_vec_.begin(), cell_vec_.end(), eit->Neighbor(i)) )
                {
                  // boundary faces
                  boundary_faces.insert( make_pair( eit, i ) );
                  // boundary nodes
                  assert( eit->FE() != nullptr );
                  for ( const auto& j : eit->FE()->NodesOfFace(i) ) {
                       assert( eit->N(j) );
                       boundary_nodes.insert( eit->N(j) );
                    }
                  // counting neighbors
                  nbors_that_belong_to_group--;
                }

            // storing the distinguished cells in the respective vectors
            // ------------------------------------------------------------
            // interior cells
              if ( nbors_that_belong_to_group == eit->Neighbors() )
                interior_elmts.insert( eit );
              // cells with at least one face on the region boundary
              else
                boundary_elmts.insert( eit );
          }
     }
    
    // 1.2 If there are cells with different spatial dimensions
    // -----------------------------------------------------------
    //     the ones with highest dimensions are used to define perimeter
    //     all lower dimensional mesh that sticks out is flagged as perimeter as well.
    else
      {
        // a. identify the boundary cells among the highest dimensional cells,
        //    also collecting all their node pointers into a set.
        set<CELL<dim>*> lesser_dim_elmts;
        set<Node<dim>*> highest_dim_elmt_nodes;

        for ( auto& eit : cell_vec_ )
          {
              // cells of the highest spatial dimension are used to define the boundary
              assert( parseFiniteElementDimension( eit->FE_Type() ) != 0 );
              if ( parseFiniteElementDimension( eit->FE_Type() ) == elmt_dim.second )
                {
                   // creating a subset with their nodes
                   const auto n_nodes{ eit->Nodes() };
                   for ( auto i{0U}; i<n_nodes; ++i ) {
                        assert( eit->N(i) != nullptr );
                        highest_dim_elmt_nodes.insert( eit->N(i) );
                     }
                   // if the cell has faces that lie on the region boundary
                   // it is considered a boudary cell
                   long  nbors_that_belong_to_group(eit->Neighbors());
                   const auto n_faces{ eit->Faces() };
                   for ( auto i{0U}; i<n_faces; ++i )
                     // 1) the cell is on model boundary  or  2) one of its neighbors does not belong to its parent region
                     if ( !eit->Neighbor(i) || !binary_search( cell_vec_.begin(), cell_vec_.end(), eit->Neighbor(i) ) )
                       {
                          // the cell pointer and the face number are used to create a unique key for the discovered boundary face
                          boundary_faces.insert( make_pair( eit, i ) );
                          // recording the nodes of the boundary face as boundary nodes
                          for ( const auto& j : eit->FE()->NodesOfFace(i) )
                            boundary_nodes.insert( eit->N(j) );

                          nbors_that_belong_to_group--;
                       }
                   if ( nbors_that_belong_to_group == eit->Neighbors() ) interior_elmts.insert( eit );
                   else boundary_elmts.insert( eit );
                }
               else lesser_dim_elmts.insert( eit );
            }
        assert( /* all cells are accounted for */ cell_vec_.size() == interior_elmts.size() + boundary_elmts.size() + lesser_dim_elmts.size() );

#ifdef MODEL_SUBDOMAIN_DEBUG
cout <<"\nModelSubDomain<" << dim <<">::PartitionCellVector: '"<< Name() <<"' highest subdomain cell dim: "<< elmt_dim.second <<"\n";
cout <<"\n\tlesser-dim cells: "<< lesser_dim_elmts.size() <<", boundary cells: "<< boundary_elmts.size();
cout <<"\n\ttotal nodes: "<< node_vec_.size() <<", highest-dim cell nodes: "<< highest_dim_elmt_nodes.size() << endl;
cout.flush();
#endif

         // ---------------------------------------------------------------------------
         // 1.3  processing lower dimensional cells and their nodes in the subdomain
         // ---------------------------------------------------------------------------
         // 1.3.1 nodes that are not contained in the higher-dimensional cell subset are identified as extra boundary node
         for ( auto& it : lesser_dim_elmts ) {
              const auto n_nodes{ it->Nodes() };
              for ( auto i{0U}; i<n_nodes; ++i )
                if ( highest_dim_elmt_nodes.find( it->N(i) ) == highest_dim_elmt_nodes.end() )
                  boundary_nodes.insert( it->N(i) );
           }

         // 1.3.2 finding the lesser dimensional cells on the region boundary
         set<CELL<dim>*> lesser_dim_elmts_detached; // to distinguish stand-alone lower dimensional mesh

         for ( auto it=lesser_dim_elmts.begin(); it!=lesser_dim_elmts.end(); ++it )
           {
              // a) lower-dim cells sticking out
              // ----------------------------------
              // lower-dimensional cells with nodes that do not belong to the node set of the
              // higher dimensional cells must be boundary cells
              size_t  exterior_nodes(0U);
              for ( auto i{0U}; i<(*it)->Nodes(); ++i )
                if ( highest_dim_elmt_nodes.find( (*it)->N(i) ) == highest_dim_elmt_nodes.end() )
                  exterior_nodes++;

              // if individual nodes stick out the parent cell sticks out as well.
              if ( exterior_nodes >= 1U ) {
                   // adding boundary cells and boundary faces
                   const auto n_faces{ (*it)->Faces() };
                   for ( auto i{0U}; i<n_faces; ++i )
                     if ( !(*it)->Neighbor(i) || lesser_dim_elmts.find( static_cast<CELL<dim>*>((*it)->Neighbor(i)) ) == lesser_dim_elmts.end() ) {
                          // the cell is a boundary cell that sticks out of the region
                          boundary_elmts.insert( (*it) );
                          boundary_faces.insert( make_pair( (*it), i ) );
                       }
                   // if the entire cell sticks out of the higher dimensional domain, it is saved for a warning issued later
                   if ( exterior_nodes == (*it)->Nodes() )
                     lesser_dim_elmts_detached.insert( (*it) );
                }
              // b) lower-dim cell that shares all their nodes with the higher dimensional ones
              // ---------------------------------------------------------------------------------
              // are boundary cells if they have at least one face on the subdomain boundary:
              // - for line cells this means at least one node
              // - for surface cells this means at least one edge
              else {
                   // line cells (assuming that the faces correspond to the nodes)
                   if ( (*it)->IsLine() ) {
                        const auto n_nodes{ (*it)->Nodes() };
                        for ( auto i{0U}; i<n_nodes; ++i )
                          // if the node is a boundary noode
                          if ( boundary_nodes.find( (*it)->N(i) ) != boundary_nodes.end() ) {
                               boundary_elmts.insert( (*it) );
                               // boundary faces for line cells
                               boundary_faces.insert( make_pair( (*it), i ) );
                            }
                     }
                   // surface cells (this will only be possible in a 3D model) are on the boundary
                   // if they share at least one face with it
                   else {
                        const auto n_faces{ (*it)->Faces() };
                        for ( auto i{0U}; i<n_faces; ++i ) {
                            size_t  bnodes{0U};
                            for ( const auto& j : (*it)->FE()->NodesOfFace(i) )
                              if ( boundary_nodes.find( (*it)->N(j) ) != boundary_nodes.end() )
                                bnodes++;
                            // if all the nodes of at least one face lie at the boundary, so does the cell
                            if ( bnodes == (*it)->FE()->NodesPerFace(i) ) {
                                 boundary_elmts.insert( (*it) );
                                 // boundary faces
                                 boundary_faces.insert( make_pair( (*it), i ) );
                              }
                          }
                     }
                }
           }

         // adding lesser-dimensional cells that are not at the boundary to the interior domain
         for ( auto& it : lesser_dim_elmts )
           if ( boundary_elmts.find( it ) == boundary_elmts.end() )
             interior_elmts.insert( it );

         if ( !lesser_dim_elmts_detached.empty() ) {             
              csmp_error.Note( WARNING, "ModelSubdomain<dim,CELL>::PartitionCellVector:", Name().c_str(),
                                "subdomain contains lower-dimensional cells detached from higher dimensional domain; these will be treated as boundary.");
              // do some additional diagnostics on these cells
              // ------------------------------------------------
              cerr <<"\n\tdetached cells: "<< lesser_dim_elmts_detached.size() <<":";
              //for ( typename set<CELL<dim>*>::const_iterator
              //      it=lesser_dim_elmts_detached.begin(); it!=lesser_dim_elmts_detached.end(); ++it ) cerr <<" "<< (*it)->Idx();
              cerr << endl;
           }

#ifdef MODEL_SUBDOMAIN_DEBUG
cout <<"\n\ttotal cells: "<< cell_vec_.size() <<", interior ones: "<< interior_elmts.size() <<", boundary cells: "<< boundary_elmts.size();
cout <<", stand-alone lower-dim cells: "<< lesser_dim_elmts_detached.size();
cout <<" (sum="<< interior_elmts.size() + boundary_elmts.size() <<").";
cout <<"\n\ttotal nodes: "<< node_vec_.size() <<", boundary nodes: "<< boundary_nodes.size() << endl << endl;
cout.flush();
#endif

       } // end multi-dim cell region


    // --------------------------------------------------
    // 2. rebuilding the cell vector
    // --------------------------------------------------
    assert( interior_elmts.size() + boundary_elmts.size() == cell_vec_.size() );
    // appending the boundary cell vector<double> to the interior cell vector
    cell_vec_.assign( interior_elmts.begin(), interior_elmts.end() );
    back_insert_iterator<vector<CELL<dim>*> >  back_it(cell_vec_);
    copy( boundary_elmts.begin(), boundary_elmts.end(), back_it );

    // swap trick to trim excess memory from end of vector
    vector<CELL<dim>*>( cell_vec_ ).swap( cell_vec_ );


#ifdef MODEL_SUBDOMAIN_DEBUG
// TESTING - is there an cell with a boundary face that is not in the boundary cell vector and vice versa
bool no_error_yet(true);
set<CELL<dim>*> elmts_with_bfaces;
for ( auto it=boundary_faces.begin(); it!=boundary_faces.end(); ++it ) {
      if ( boundary_elmts.find( (*it).first ) == boundary_elmts.end() ) {
           if ( no_error_yet ) {
                cerr <<"\nboundary face parent cells vs. boundary cells:\n";
                no_error_yet=false;
             }
           cerr <<" "<< (*it).first->Idx() <<": "<< parseFiniteElementType( (*it).first->FE_Type() );
        }
      elmts_with_bfaces.insert( (*it).first );
    }
if ( elmts_with_bfaces.size() != boundary_elmts.size() )
  cerr <<"\n\tmismatch between boundary face-  and boundary cell set: "<< elmts_with_bfaces.size() <<" vs. "<< boundary_elmts.size() << endl;
assert( elmts_with_bfaces.size() == boundary_elmts.size() );
#endif

    // --------------------------------------------------
    // 3. creating the boundary face vector
    // --------------------------------------------------
    if ( !bd_face_vec_.empty() ) bd_face_vec_.clear();
    // bd_face_vec_ is only available if there are boundary cells
    if (boundary_elmts.size() > 0)
      {
        bd_face_vec_.reserve(cell_vec_.size() - boundary_elmts.size());
        //       parent cell of face, face
        typename set<pair<CELL<dim>*,uint32_t> >::const_iterator  bfit(boundary_faces.begin());
        typename set<pair<CELL<dim>*,uint32_t> >::const_iterator  ffit(boundary_faces.begin());
        vector<uint32_t>  bface_data;
        size_t            counter(0U);

        while ( bfit != boundary_faces.end() ) {
            assert((*ffit).first == cell_vec_[counter + interior_elmts.size()]);
            bface_data.reserve(3);
            // as long as we considering faces of the same cell
            while ( (*bfit).first == (*ffit).first ) {
                bface_data.push_back((*bfit).second);
                bfit++;
                if ( bfit == boundary_faces.end() ) break;
              }
            ffit = bfit;
            assert( bface_data.size() > 0 );
            bd_face_vec_.push_back(bface_data);
            bface_data.clear();
            counter++;
          }
        bd_face_vec_.shrink_to_fit();
        // debug checks
        for ( auto it=bd_face_vec_.begin(); it!=bd_face_vec_.end(); ++it )
          assert( (*it).size() >= 1 );
        assert(bd_face_vec_.size() == Cells() - InteriorCells());
     }


    // --------------------------------------------------
    // 4. building and partitioning the node vector
    //   (which SplitBoundary objects do not have)
    // --------------------------------------------------
    if ( !node_vec_.empty() )
      {
        assert( node_vec_.size() >= boundary_nodes.size() );
        first_bd_node_ = node_vec_.size() - boundary_nodes.size();
        
        // rebuilding and sorting the node vector (noting that set nodes are already sorted)
        // ---------------------------------------------------------------------------------
        sort( node_vec_.begin(), node_vec_.end() );
        // a temporary new node vector is created
        vector<csmp::Node<dim>*>  temp;
        temp.reserve(node_vec_.size());
        // all nodes that are not in the boundary node vector are considered as interior nodes
        for ( const auto& nit : node_vec_ ) {
              assert( nit );
              if ( boundary_nodes.find(nit) == boundary_nodes.end() )
                temp.push_back( nit );
          }
        // second, the already sorted perimeter nodes are appended
        for ( const auto& nit : boundary_nodes )
          temp.push_back( nit );
        // now the temporary vector is assigned to the permanent one
        assert( temp.size() == node_vec_.size() );
        node_vec_ = std::move( temp );
        node_vec_.shrink_to_fit();
        
        assert( first_bd_node_ <= node_vec_.size() );
      }
#ifdef MODEL_SUBDOMAIN_DEBUG
cout <<"\nModelSubDomain<dim,CELL>::PartitionCellVector: '"<< Name() <<"': of the ";
cout << cell_vec_.size() <<" cells, "<< boundary_elmts.size() <<" lie at the domain boundary."<< endl;
cout.flush();
#endif

    return cell_vec_.size() - boundary_elmts.size();

 } // end PartitionCellVector






/**
    Uses vector to create unique node vector by pushing back all node of the elements including duplicates,
    then sorting it and eliminating the duplicates.
    
    TODO: speed critical function. Perhaps refactor with unordered set as intermediate container for unique nodes because there will be so many duplicates.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::CreateNodePointerVector()
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( cell_vec_.empty() ) {
        csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::CreateNodePointerVector",
                                  "cannot create 'node_vec_', current cell vector is empty");
        return;
     }

  if ( !node_vec_.empty() )
    node_vec_.clear();

  // creating the node pointer vector
  node_vec_.reserve( cell_vec_.size() * 4 );
  for ( auto& it : cell_vec_ ) {
       const auto nodes{ it->Nodes() };
       for ( auto i{0U}; i<nodes; i++ ) {
            assert( it->N( i ) != nullptr );
            node_vec_.push_back( it->N( i ) );
         }
    }

  // removing duplicates and trimming excess memory from node vector
  sort( node_vec_.begin(), node_vec_.end() );
  node_vec_.erase( unique( node_vec_.begin(), node_vec_.end() ), node_vec_.end() );
  node_vec_.shrink_to_fit();
}






    /// creates node pointer vector from the shared face nodes of the supplied range of contacting cells
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::CreateNodePointerVector( std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& contacting_cells )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( contacting_cells.empty() ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::CreateNodePointerVector",
                                    "provided argument vector does not contain any cells");
          return;
       }

  if ( !node_vec_.empty() ) node_vec_.clear();
  node_vec_.reserve( contacting_cells.size() ); // rough guess of n-nodes
   
  // creating the node pointer vector
  for ( const auto& it : contacting_cells ) {
       // taking the nodes from the inside of the cell pairs
       for ( const auto& i : it.first.first->FE()->NodesOfFace( it.first.second ) )
         node_vec_.push_back( it.first.first->N(i) );
    }

  // removing duplicates and trimming excess memory from node vector
  sort( node_vec_.begin(), node_vec_.end() );
  node_vec_.erase( unique( node_vec_.begin(), node_vec_.end() ), node_vec_.end() );
  node_vec_.shrink_to_fit();

 } // end CreateNodePointerVector





/**
      SKM - September 2019.
      Jan 2020, updated for the case where there are no interior nodes or cells
 */
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::SortVectors( size_t interior_cells, size_t interior_nodes )
  {
      // verification of suitable model state
      if ( cell_vec_.empty() )
        throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::SortVectors",
                               Name(), "model subdomain: CELL vector not initialised yet.");
                               
      if ( interior_cells > cell_vec_.size() ) {
            cerr <<"\ninterior cells: "<< interior_cells;
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::SortVectors",
                                   Name(), "model subdomain: less CELLs in CELL vector than interior cells specified.");
        }
     if ( node_vec_.empty() )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::SortVectors",
                              Name(), "model subdomain: node vector not initialised yet.");
                              
     if ( interior_nodes > node_vec_.size() ) {
           cerr <<"\ninterior cells: "<< interior_nodes;
           throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::SortVectors",
                                  Name(), "model subdomain: less nodes in node vector than interior nodes specified.");
       }
       
     // sorting
     // -------
     // dealing with the case where all cells and nodes are located on the regions perimeter,
     // i.e. there are no interior nodes or cells so that there is only a single sorted range
     if ( interior_nodes == 0 ) {
          sort( node_vec_.begin(), node_vec_.end() );
       }
     else {
         sort( node_vec_.begin(), next(node_vec_.begin(),interior_nodes) );
         sort( next(node_vec_.begin(),interior_nodes), node_vec_.end() );
       }
     
     if ( interior_cells == 0 ) {
          sort( cell_vec_.begin(), cell_vec_.end() );
       }
     else {
         sort( cell_vec_.begin(), next(cell_vec_.begin(),interior_cells) );
         sort( next(cell_vec_.begin(),interior_cells), cell_vec_.end() );
       }

  } // end SortVectors






// PLACEMENT OF PROPERTY ASSIGNMENT

SUBDOMAIN_PART parseSubdomainPart( const char* subdomain )
{
    std::string ssubdomain(subdomain);

    std::transform(ssubdomain.begin(),ssubdomain.end(),ssubdomain.begin(),::toupper);

    if (  ssubdomain == "COMPLETE"  )  return COMPLETE;
    if (  ssubdomain == "INTERIOR"  )  return INTERIOR;
    if (  ssubdomain == "PERIMETER" )  return PERIMETER;
    return COMPLETE;
}


std::string parseSubdomainPart( SUBDOMAIN_PART ssubdomain )
 {
    if (  ssubdomain == COMPLETE  )  return std::string("COMPLETE");
    if (  ssubdomain == INTERIOR  )  return std::string("INTERIOR");
    if (  ssubdomain == PERIMETER )  return std::string("PERIMETER");
    return std::string("COMPLETE");
}


/*
template<uint32_t dim, template<uint32_t> class CELL>
PLACEMENT ModelSubDomain<dim,CELL>::Placement() const 
 { 
    throw csmp::Exception( ERROR, "ModelSubdomain<>::Placement:", "erratic call"); 
    return UNDEFINED; 
 }
*/

template<uint32_t dim, template<uint32_t> class CELL>
bool ModelSubDomain<dim,CELL>::ValidVariable( const char* variableName ) const
  {
    const PLACEMENT p( pref_.Placement( variableName ) );
    if ( p == NODE || p == ELEMENT )
      return true;
    return false;
  }


// INTERRELATIONS INTERFACE

/** applies Interrelation to ModelSubDomain
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::Apply( Interrelation<dim>& relation )
 {
    relation.Apply( *this );

 } // end Apply


// VISITORS INTERFACE


/// virtual function stub that will be overwritten by base classes
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::Accept( csmp::Visitor<dim>& )
 {
    throw csmp::Exception( ERROR,
                           "ModelSubDomain<dim,CELL>::Accept",
                           "Subclass method should be called; nothing was done");

 } // end Accept


// INTEGRATION POINTS

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::IntegrationPoints() const
 {
    if ( !(*cell_vec_.begin())->FE()->UsesLocalCoordinates() ) return 0U;

    size_t  cpoints(0U);
    for( typename vector<CELL<dim>*>::const_iterator
         it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
      cpoints += (*it)->IntegrationPoints();

    return cpoints;
 }


template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::SectorIntegrationPoints() const
  {
  if ( !(*cell_vec_.begin())->FE()->UsesLocalCoordinates() ) return 0U;

  size_t  cpoints(0U);
  for( typename vector<CELL<dim>*>::const_iterator
    it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
    cpoints += (*it)->Sectors()*(*it)->IntegrationPointsPerSector();

  return cpoints;
  }

template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::FacetIntegrationPoints() const
  {
  if ( !(*cell_vec_.begin())->FE()->UsesLocalCoordinates() ) return 0U;

  size_t  cpoints(0U);
  for( typename vector<CELL<dim>*>::const_iterator
    it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
    cpoints += (*it)->Facets()*(*it)->IntegrationPointsPerFacet();

  return cpoints;
  }




// INDEXES

template<uint32_t dim, template<uint32_t> class CELL>
int32_t  ModelSubDomain<dim,CELL>::DomainIndex() const
{
  return domain_idx_;
}

/**

Returns a vector<double> with the ID numbers of the Elements which belong
to the Region.
*/
template<uint32_t dim, template<uint32_t> class CELL>
vector<size_t>  ModelSubDomain<dim,CELL>::MemberCellIndexes() const
 {
    vector<size_t> ids;
    ids.reserve( cell_vec_.size() );
    for ( const auto& it : cell_vec_ ) ids.push_back( it->Idx() );
    
    return ids;
 }



/** Renumbers nodes from 0 to n-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::RenumberNodes() const
 {
    size_t  counter(0U);

    for ( auto& it : node_vec_ ) it->Idx( counter++ );

    return counter;
 }



/** Renumbers cells from 0 to n-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::RenumberCells() const
 {
    size_t counter(0U);

    for( auto& it : cell_vec_ ) it->Idx(counter++);
    cell_vec_[0]->FE()->CurrentID( numeric_limits<size_t>::max() );

    return counter;
    
 } // end RenumberCells



/** Renumbers cells and nodes from 0 to n-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::UpdateMemberIndexes() const
 {
    RenumberNodes();
    RenumberCells();
   
 } // end UpdateRegionMemberIndexes





// ACCESSORY

/**
    Uses binary_search on both ranges of the sorted cell vector to find the cell in question.
    returns true or false.
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool ModelSubDomain<dim,CELL>::Contains( const CELL<dim>* const eptr ) const
 {
    assert( eptr != nullptr );

    // search perimeter first
    if ( binary_search( next(cell_vec_.begin(),InteriorCells()), cell_vec_.end(), eptr ) )
       return true;

    if ( binary_search( cell_vec_.begin(), next(cell_vec_.begin(),InteriorCells()), eptr ) )
       return true;

    return false;

 } // end


/**
    Uses binary_search on both ranges of the sorted node vector to find the node in question.
    returns true or false.
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool ModelSubDomain<dim,CELL>::Contains( const Node<dim>* const nptr ) const
 {
    assert( nptr != nullptr );

    if ( binary_search( next(node_vec_.begin(), InteriorNodes()), node_vec_.end(), nptr ) )
       return true;

    if ( binary_search( node_vec_.begin(), next(node_vec_.begin(), InteriorNodes()), nptr ) )
       return true;

    return false;

 } // end

template<uint32_t dim, template<uint32_t> class CELL>
bool ModelSubDomain<dim,CELL>::IsPerimeterNode( const csmp::Node<dim>* const nd_ptr ) const
 {
    assert( nd_ptr != nullptr );
    return std::binary_search( PerimeterNodesBegin(), NodesEnd(), nd_ptr );
 }


template<uint32_t dim, template<uint32_t> class CELL>
bool  ModelSubDomain<dim,CELL>::IsPerimeterCell( const CELL<dim>* const e_ptr ) const
 {
    assert( e_ptr != nullptr );
    return std::binary_search( PerimeterCellsBegin(), CellsEnd(), e_ptr );
 }


// still used by legacy NodeCenteredFiniteVolumeTransport
template<uint32_t dim, template<uint32_t> class CELL>
bool  ModelSubDomain<dim,CELL>::IsPerimeterNode( size_t i ) const
 {
    if ( i >= Nodes() ) return false;
    return ( i >=  first_bd_node_ );
 }

/**
    Returns true when the queried cell index (0..n-1)
    is among those of the cells located on the boundary of the region;
    else false.
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool  ModelSubDomain<dim,CELL>::IsPerimeterCell( size_t e ) const
 {
    if ( e >= cell_vec_.size() ) return false;
    return (e < InteriorCells()) ? false : true;
 }



// GEOMETRY


template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const
 {
    xyz_min = xyz_max = (*node_vec_.begin())->Coordinate();

    // only nodes at the group boundary have to be checked
    for ( auto bit=PerimeterNodesBegin(); bit!=NodesEnd(); bit++ )
      {
         csmp::Point<dim> p = (*bit)->Coordinate();
         xyz_min[0] = std::min( p[0], xyz_min[0] );
         xyz_max[0] = std::max( p[0], xyz_max[0] );
         if constexpr( dim != 1U ) {
             xyz_min[1] = std::min( p[1], xyz_min[1] );
             xyz_max[1] = std::max( p[1], xyz_max[1] );
           }
         if constexpr( dim == 3U ) {
             xyz_min[2] = std::min( p[2], xyz_min[2] );
             xyz_max[2] = std::max( p[2], xyz_max[2] );
           }
      }

 } // end MinMaxCoordinates

/**

AssignCellCharacteristicsTo() allows to assign a number of Element
characteristics as identified by strings (second argument) to scalar physical
variables. These characteristics are:

"volume" (3D)
"inner radius"
"aspect ratio" (longest / shortest segment)

@section arguments Input Arguments

The cell characteristic indicated by the first string argument is
assigned to a property of the user's choice specified by the second method
argument and as identified by its name in the variable database.

@section implementation Implementation

The method was implemented mainly for two-dimensional calculations.
Therefore, height and width only work in two-dimensional calculations.
The characteristics are obtained from the FiniteElement which is
bridged to the current Element that is queried.

@section application Application

The characteristic "inner radius" is commonly used to find the
appropriate resolution for an advection or visualization grid on which a
variable is to be mapped on.

AssignCellCharacteristicsTo() can be used to:
- test the shape of cells in a mesh for their suitability for a
  computation (aspect ratio).
- testing how skewed the cells got by deformation
- integrate cell properties in specific calculations, for instance,
  if a fracture is just one cell wide, the minimum height of
  these fracture cells may be equivalent to the fracture aperture.

@section messages Messages

Firstly, the method will report an error and return, if the target
property to which the cell characteristic shall be assigned to is not
an cell variable.

Further errors may be reported if the user tries to use the characteristic
area in a 3D computation, or volume in a 2D computation. Also, cell
height and width are thus far only available in 2D.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo( const char* characteristic,
                                                            const char* var )
 {
     ErrorHandler&     csmp_error( ErrorHandler::Instance() );
     const csmp::Index prop_key = pref_.StorageKey(var);
     ScalarVariable    sc;

     if ( prop_key.place != ELEMENT and prop_key.type != SCALAR ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo",
                                    "method assigns only to scalar cell properties");
          return;
       }
     if ( dim != 1U and !strcmp( characteristic, "length" ) )
       csmp_error.Note( WARNING, "ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo",
                                   "'length' will only be assigned to line cells, nothing done to others");

     if ( dim != 2U and !strcmp( characteristic, "area" ) )
       csmp_error.Note( WARNING, "ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo",
                                   "'area' will only be assigned to surface cells, nothing done to others");

     if ( dim != 3U and !strcmp( characteristic, "volume" ) ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo",
                                    "'volume' can only be assigned to volume cells. Nothing was done");
          return;
       }

     if ( !strcmp( characteristic, "length" ) )
       for ( auto& eit : cell_vec_ )
         {
             if ( eit->IsLine() ) {
                  sc = fabs( eit->Volume() );
                  eit->Store( prop_key, sc );
               }
         }
     else if ( !strcmp( characteristic, "area" ) )
       for ( auto& eit : cell_vec_ )
         {
             if ( eit->IsSurface() ) {
                  sc = fabs( eit->Volume() );
                  eit->Store( prop_key, sc );
               }
         }
     else if ( !strcmp( characteristic, "volume" ) )
       for ( auto& eit : cell_vec_ )
         {
             if ( eit->IsVolume() ) {
                  sc = fabs( eit->Volume() );
                  eit->Store( prop_key, sc );
               }
         }
     else if ( !strcmp( characteristic, "aspect ratio" ) )
       {
          for ( auto& eit : cell_vec_ )
            {
               sc = eit->AspectRatio();
               eit->Store( prop_key, sc );
            }
      }
     else if ( !strcmp( characteristic, "inner radius" ) )
       {
          for ( auto& eit : cell_vec_ )
            {
               sc = eit->InnerRadius();
               eit->Store( prop_key, sc );
            }
      }
     else {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo",
                             characteristic,  " entered as characteristic was not identified; nothing done");
          cout <<"\nYour options are: "<< endl;
          cout <<"\n\t inner radius"   << endl;
          cout <<"\t aspect ratio"     << endl;
          cout <<"\t length" << endl;
          cout <<"\t area (>=2D models)" << endl;
          cout <<"\t volume (3D models only)" << endl;
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::AssignCellCharacteristicsTo", "Now exciting");
      }

 } // end





/**

Assigns node coordinates to the supplied node variable which may be either
of scalar or vector<double> type. If either the X, Y, or Z coordinate are to be
assigned to a scalar variable an index from 0...dim-1 must be supplied
to indicate which coordinate axis shall be mapped.

@section arguments Input Arguments

The name of the physical variable to which the node coordinate shall
be assigned to, and, if a scalar property shall be initialized with
coordinate values, an int index to indicate which coordinate axis is
desired (0=X, 1=Y, 2=Z).

@section implementation Implementation

Method uses the Coordinates() and X(), Y(), and Z() interfaces of the
Node class.

@section application Application

To arrive at coordinate values for a model where topography is
important etc.

@section messages Messages

The method will report a fatal error and terminate the computation, if
the target variable is not placed on the node or if it is not of the
correct type.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo(  const char* vector_variable )
 {
     const csmp::Index  prop_key = pref_.StorageKey(vector_variable);
     VectorVariable<dim>  vc;

      if ( prop_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "Node coordinates can only be assigned to a vector<double> variable");

     if ( prop_key.place != NODE )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "Node coordinates can only be assigned to a node variable");

     for ( auto& nit : node_vec_ ) {
          vc = nit->Coordinate();
          nit->Store( prop_key, vc );
       }

 } // end AssignNodeCoordinatesTo



template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo( const char* scalar_variable, char c )
 {
     const csmp::Index prop_key = pref_.StorageKey(scalar_variable);
     ScalarVariable  sc;

     if ( c != 'x' && c != 'y' && c != 'z' &&
          c != 'X' && c != 'Y' && c != 'Z')
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "coordinate index letter is not valid",
                                    "should be either of 'x,y,z,X,Y,Z'");

     if ( prop_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "Individual node coordinates can only be assigned to a scalars");

     if ( prop_key.place != NODE )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "Node coordinates can only be assigned to a node variable");

     if ( dim == 2U and c == 'z' )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "2D models have no z coordinates");

     if ( dim == 1U and c != 'x' )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::AssignNodeCoordinatesTo",
                                    "1D models only have x coordinates");

     for ( auto& nit : node_vec_ )
       {
         if      ( c == 'x' || c == 'X' ) sc = nit->x();
         else if ( c == 'y' || c == 'Y' ) sc = nit->y();
         else if ( c == 'z' || c == 'Z' ) sc = nit->z();
         nit->Store( prop_key, sc );
       }

 } // end AssignNodeCoordinatesTo










// MIN-MAX PROPERTIES


/** Returns min/max of property values inside a region into its arguments.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::MinMaxOf( const char* property, double& gmin, double& gmax ) const
 {
    ErrorHandler&  csmp_error(ErrorHandler::Instance());

    if ( !pref_.IsDefined(property) ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::MinMaxOf", property, "is undefined; nothing could be done" );
          return;
      }
    csmp::Index  gprop_key(pref_.StorageKey(property));
    MinMaxOf( gprop_key, gmin, gmax );

 } // end MinMaxOf




/**
    Recovering and sorting the minimum and maximum Eigen values of the tensor variable.

    @todo SKM this wants to be a lambda function in the next method.
*/
template<uint32_t dim>
void minMaxEigenValues( const TensorVariable<dim>& ts, double& tmin, double& tmax )
 {
    VectorVariable<dim>  evals;
    ts.EigenValues( evals );
    std::set<double> min_max;
    for ( auto i{0U}; i<dim; i++ ) min_max.insert( evals[i] );
    tmin = (*min_max.begin());
    tmax = (*min_max.rbegin());
 }




/**
    Recovers ranges of variable values from regions, boundaries or split boundaries.
    
    @attention for vector variables the length range is returned.
    @attention for tensor variables the maximum Eigenvalue is recovered.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::MinMaxOf( const csmp::Index& prop_key, double& vmin, double& vmax ) const
 {
     if( prop_key.place == REGION || prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY || prop_key.place == MODEL )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::MinMaxOf",
                             "Model/Region/SplitBoundary/Boundary properties have constant values within individual model subdomains; read directly!");

   // the 'Model' region was already dealt with by Model::MinMaxOf

   // node properties of any kind
   if ( prop_key.place == NODE ) {
       typename vector<csmp::Node<dim>*>::const_iterator  nit(node_vec_.begin());
       assert( (*nit) != nullptr );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*nit)->Read(prop_key); nit++;
            while ( nit!=node_vec_.end() ) {
                 vmin = std::min( vmin, (*nit)->Read(prop_key) );
                 vmax = std::max( vmax, (*nit)->Read(prop_key) );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*nit)->Read(prop_key,vc); nit++;
            vmin = vmax = vc.Length();
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, vc );
                 const double  vlength(vc.Length());
                 vmin = std::min( vmin, vlength );
                 vmax = std::max( vmax, vlength );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*nit)->Read(prop_key,ts); nit++;
            minMaxEigenValues( ts, vmin, vmax );
            double tmin, tmax;
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, ts );
                 minMaxEigenValues( ts, tmin, tmax );
                 vmin = std::min( vmin, tmin );
                 vmax = std::max( vmax, tmax );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*nit)->Read(prop_key,a); nit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*nit)->Read(prop_key,a); nit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 nit++;
              }
            return;
         }
    } // end node properties


   // cell-based properties of any kind
   if ( prop_key.place == ELEMENT || prop_key.place == FACE || prop_key.place == INTER_FACE ) {
       typename vector<CELL<dim>*>::const_iterator  eit(cell_vec_.begin());
       assert( (*eit) != nullptr );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(prop_key); eit++;
            while ( eit!=cell_vec_.end() ) {
                 vmin = std::min( vmin, (*eit)->Read(prop_key) );
                 vmax = std::max( vmax, (*eit)->Read(prop_key) );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=cell_vec_.end() ) {
                 (*eit)->Read( prop_key, vc );
                 const double  vlength(vc.Length());
                 vmin = std::min( vmin, vlength );
                 vmax = std::max( vmax, vlength );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double tmin, tmax;
            while ( eit!=cell_vec_.end() ) {
                 (*eit)->Read( prop_key, ts );
                 minMaxEigenValues( ts, tmin, tmax );
                 vmin = std::min( vmin, tmin );
                 vmax = std::max( vmax, tmax );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                 (*eit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                 (*eit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 eit++;
              }
            return;
         }
    } // end cell, face, interface properties


   // cell-based integration-point properties of any kind
   if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
       typename vector<CELL<dim>*>::const_iterator  eit(cell_vec_.begin());
       assert( (*eit) != nullptr );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(0U,prop_key); eit++;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     vmin = std::min( vmin, (*eit)->Read(i,prop_key) );
                     vmax = std::max( vmax, (*eit)->Read(i,prop_key) );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(0U,prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=cell_vec_.end() ) {
                for ( auto i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, vc );
                     const double  vlength(vc.Length());
                     vmin = std::min( vmin, vlength );
                     vmax = std::max( vmax, vlength );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(0U,prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double tmin, tmax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, ts );
                     minMaxEigenValues( ts, tmin, tmax );
                     vmin = std::min( vmin, tmin );
                     vmax = std::max( vmax, tmax );
                 }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
    } // end cell integration point properties on cell, face or interface


   // FV facet integration-point properties of any kind placed on cell, face, or interface
   if ( prop_key.place == FACET_INTEGRATION_POINT ||
        prop_key.place == FACE_FACET_INTEGRATION_POINT ||
        prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT )
     {
       typename vector<CELL<dim>*>::const_iterator  eit(cell_vec_.begin());
       assert( (*eit) != nullptr );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(0U,0U,prop_key); eit++;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Facets(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     vmin = std::min( vmin, (*eit)->Read(i,j,prop_key) );
                     vmax = std::max( vmax, (*eit)->Read(i,j,prop_key) );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(0U,0U,prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Facets(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, vc );
                     const double  vlength(vc.Length());
                     vmin = std::min( vmin, vlength );
                     vmax = std::max( vmax, vlength );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(0U,0U,prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double tmin, tmax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Facets(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, ts );
                     minMaxEigenValues( ts, tmin, tmax );
                     vmin = std::min( vmin, tmin );
                     vmax = std::max( vmax, tmax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Facets(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Facets(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
    } // end FV facet integration point properties on cell, face or interface


   // FV sector integration-point properties of any kind placed on cell, face, or interface
   if ( prop_key.place == SECTOR_INTEGRATION_POINT ||
        prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
        prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT  )
     {
       typename vector<CELL<dim>*>::const_iterator  eit(cell_vec_.begin());
       assert( (*eit) != nullptr );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(0U,0U,prop_key); eit++;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Sectors(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     vmin = std::min( vmin, (*eit)->Read(i,j,prop_key) );
                     vmax = std::max( vmax, (*eit)->Read(i,j,prop_key) );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(0U,0U,prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Sectors(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, vc );
                     const double  vlength(vc.Length());
                     vmin = std::min( vmin, vlength );
                     vmax = std::max( vmax, vlength );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(0U,0U,prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double tmin, tmax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Sectors(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, ts );
                     minMaxEigenValues( ts, tmin, tmax );
                     vmin = std::min( vmin, tmin );
                     vmax = std::max( vmax, tmax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Sectors(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double  amin, amax;
            while ( eit!=cell_vec_.end() ) {
                for ( auto i{0U}; i<(*eit)->Sectors(); i++ )
                  for ( auto j{0U}; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
     } // end sector integration points for any cell, face or interface

    vmin = vmax = std::numeric_limits<double>::signaling_NaN();

    ErrorHandler&  csmp_error(ErrorHandler::Instance());
    csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::MinMaxOf(index)",
                      "placement of property coould not be indentified");

 } // end MinMaxOf(index)













// MANIPULATIONS WITH PROPERTIES




/// @todo (3-D) Refactor!! This is a potential bug if wrong integration properties are used!
template<uint32_t dim, template<uint32_t> class CELL>
template<typename Var>
void ModelSubDomain<dim,CELL>::InputPropertyValue( const char* input_prop,
                                                   const Var& var,
                                                   SUBDOMAIN_PART sdp )
 {
     const csmp::Index prop_key = pref_.StorageKey(input_prop);

     if ( sdp == PERIMETER  and (prop_key.place == REGION  or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY) )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                            input_prop, "placed on Region cannot be assigned just on perimeter");

     if ( prop_key.place == REGION || prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                            input_prop, "use Store() to assign Region/SplitBoundary/Boundary values");

     if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::InputPropertyValue",
                            input_prop, "for InterFace nodes, node property needs to be assigned with SplitBoundary::InputPropertyValue");

     if ( sdp == COMPLETE ) {
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( auto& it : cell_vec_ ) it->Store( prop_key, var );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( auto& it : cell_vec_ )
                  for ( auto i{0U}; i<it->IntegrationPoints(); i++ )
                    it->Store( i, prop_key, var );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( auto& it : cell_vec_ )
               for ( auto i{0U}; i<it->Facets(); i++ )
                 for ( auto j{0U}; j<it->IntegrationPointsPerFacet(); j++ )
                   it->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( auto& it : cell_vec_ )
               for ( auto i{0U}; i<it->Sectors(); i++ )
                 for ( auto j{0U}; j<it->IntegrationPointsPerSector(); j++ )
                   it->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( auto& nit : node_vec_ )
                  nit->Store( prop_key, var );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                                                "Property placement not recognized");
       }
     else if ( sdp == PERIMETER ) { // property is assigned only to perimeter of boundary
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
                  (*it)->Store( prop_key, var );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
                  for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ )
                    (*it)->Store( i, prop_key, var );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
               for ( auto i{0U}; i<(*it)->Facets(); i++ )
                 for ( auto j{0U}; j<(*it)->IntegrationPointsPerFacet(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
               for ( auto i{0U}; i<(*it)->Sectors(); i++ )
                 for ( auto j{0U}; j<(*it)->IntegrationPointsPerSector(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( auto nit=PerimeterNodesBegin(); nit!=node_vec_.end(); nit++ )
                  (*nit)->Store( prop_key, var );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                                                "Property placement not recognized");
       }
     else { // property is assigned only to interior of subdomain
           if ( prop_key.place == ELEMENT ) {
                for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
                     (*it)->Store( prop_key, var );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
                  for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ )
                    (*it)->Store( i, prop_key, var );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
               for ( auto i{0U}; i<(*it)->Facets(); i++ )
                 for ( auto j{0U}; j<(*it)->IntegrationPointsPerFacet(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for (  auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
               for ( auto i{0U}; i<(*it)->Sectors(); i++ )
                 for ( auto j{0U}; j<(*it)->IntegrationPointsPerSector(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                  (*nit)->Store( prop_key, var );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                                                    "Property placement not supported");
       }

 } // end InputPropertyValue

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART );



/**
    Helper function that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
*/


/**
    Assigns variable values to model subdomains, but only if the components of the variable do not have
    the flag identified as 'do_not_overwrite'.

    @attention Since tensor variables only have one flag per row, the overwrite protection is applied
    on a row-by-row basis.

    @attention This special treatment is not applied to FV integration points.

*/
template<uint32_t dim, template<uint32_t> class CELL>
template<typename Var>
void ModelSubDomain<dim,CELL>::InputPropertyValue( const char* input_prop,
                                                   const Var& var,
                                                   VARIABLE_FLAG do_not_overwrite,
                                                   SUBDOMAIN_PART sdp )
 {
     const csmp::Index prop_key = pref_.StorageKey(input_prop);

     if ( sdp == PERIMETER  and (prop_key.place == REGION  or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY) ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                            input_prop, "placed on Region cannot be assigned just on perimeter");
          return;
       }

     if( prop_key.place == REGION || prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                            input_prop, "use Store() to assign Region/SplitBoundary/Boundary valuese");

     if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::InputPropertyValue",
                            input_prop, "for InterFace nodes, node property needs to be assigned with SplitBoundary::InputPropertyValue");

     if ( sdp == COMPLETE ) {
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( typename vector<CELL<dim>*>::iterator
                      it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
                  writeVariableIf( (*it), prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<CELL<dim>*>::iterator
                      it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
                  for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ )
                    writeVariableIf( (*it), i, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<CELL<dim>*>::iterator
               it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
               for ( auto i{0U}; i<(*it)->Facets(); i++ )
                 for ( uint32_t j{0U}; j<(*it)->IntegrationPointsPerFacet(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<CELL<dim>*>::iterator
               it=cell_vec_.begin(); it!=cell_vec_.end(); it++ )
               for ( auto i{0U}; i<(*it)->Sectors(); i++ )
                 for ( uint32_t j{0U}; j<(*it)->IntegrationPointsPerSector(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=node_vec_.begin(); nit!=node_vec_.end(); nit++ )
                  writeVariableIf( (*nit), prop_key, var, DIRICH );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                                                    "Property placement not recognized");
       }
     else if ( sdp == PERIMETER ) { // property is assigned only to perimeter of boundary
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
                  writeVariableIf( (*it), prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
                  for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ )
                    writeVariableIf( (*it), i, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
               for ( auto i{0U}; i<(*it)->Facets(); i++ )
                 for ( uint32_t j{0U}; j<(*it)->IntegrationPointsPerFacet(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( auto it=PerimeterCellsBegin(); it!=cell_vec_.end(); it++ )
               for ( auto i{0U}; i<(*it)->Sectors(); i++ )
                 for ( uint32_t j{0U}; j<(*it)->IntegrationPointsPerSector(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( auto nit=PerimeterNodesBegin(); nit!=node_vec_.end(); nit++ )
                  writeVariableIf( (*nit), prop_key, var, do_not_overwrite );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                                                "Property placement not recognized");
       }
     else { // property is assigned only to interior of subdomain
           if ( prop_key.place == ELEMENT ) {
                for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
                  writeVariableIf( (*it), prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
                  for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ )
                    writeVariableIf( (*it), i, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
               for ( auto i{0U}; i<(*it)->Facets(); i++ )
                 for ( uint32_t j{0U}; j<(*it)->IntegrationPointsPerFacet(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( auto it=CellsBegin(); it!=PerimeterCellsBegin(); it++ )
               for ( auto i{0U}; i<(*it)->Sectors(); i++ )
                 for ( uint32_t j{0U}; j<(*it)->IntegrationPointsPerSector(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                  writeVariableIf( (*nit), prop_key, var, do_not_overwrite );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,CELL>::InputPropertyValue",
                                                    "Property placement not supported");
       }

 } // end InputPropertyValue

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const VectorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const VectorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const VectorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const VectorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const VectorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const VectorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const VectorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const VectorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const VectorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const TensorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const TensorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const TensorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const TensorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const TensorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const TensorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const TensorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const TensorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const TensorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );



/**
    Returns status of whole or parts of model subdomain if there is a single flag value; else ANY is returned.
    For scalars, a only a single variable flag can be returned, while for vectors and tensors the flag from component 0
    is provided by default, unless any other component number is asked for.
    Note: no matter which component is asked for in a vector or tensor, if any component has flag ANY,
    then the ANY flag is returned by default. (this is checked for in loops from 0 < dim for vectors and 0 <dim*dim for tensors)
*/
template<uint32_t dim, template<uint32_t> class CELL>
VARIABLE_FLAG  ModelSubDomain<dim,CELL>::PropertyStatus( const char* property, SUBDOMAIN_PART group_flag, uint32_t component ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( !pref_.IsDefined(property) ) {
        csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::PropertyStatus", property, "is undefined; nothing could be done" );
        return ANY;
    }
    const csmp::Index  prop_key = pref_.StorageKey(property);
    string src("ModelSubDomain<");
    string cache(to_string(dim));
    src += cache;
    src += ",";
    src += typeid(CELL<dim>).name();
    src +=">::PropertyStatus";

    if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
      throw csmp::Exception( ERROR, src.c_str(), "Use Status() to change flags of REGION/BOUNDARY/SPLIT_BUNDARY variables.");

     if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::PropertyStatus",
                            property, "for InterFace nodes, node property needs to be assigned with SplitBoundary::PropertyStatus");

    if ( prop_key.place == FACET_INTEGRATION_POINT ||
         prop_key.place == SECTOR_INTEGRATION_POINT ||
         prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
         prop_key.place == FACE_FACET_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT )
    {
        csmp_error.Note( ERROR, src.c_str(), "Method is not implemented for facet and sector integration points");
        return ANY;
    }

    if ( cell_vec_.empty() ) {
        csmp_error.Note( ERROR, src.c_str(), "ModelSubDomain is empty");
        return ANY;
    }

    if ( prop_key.type == SCALAR || prop_key.type == ARRAY )
    {
        // if all the flags in the region shall be changed
        if ( group_flag == COMPLETE ) {
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key));
                for ( const auto& eit : cell_vec_ )
                    if ( status0 != eit->Status(prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key));
                for ( const auto& eit : cell_vec_ )
                    for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                        if ( status0 != eit->Status(i,prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key));
                for ( const auto& nit : node_vec_  )
                    if ( status0 != nit->Status(prop_key) ) return ANY;
                return status0;
            }
        }
        // if the property shall only be changed on the Subdomain boundary
        else if ( group_flag == PERIMETER ) {
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(InteriorCells())->Status(prop_key));
                for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                    if ( status0 != (*eit)->Status(prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(InteriorCells())->Status(0,prop_key));
                for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                    for ( auto i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        if ( status0 != (*eit)->Status(i,prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(InteriorNodes())->Status(prop_key));
                for ( auto nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                    if ( status0 != (*nit)->Status(prop_key) ) return ANY;
                return status0;
            }
        }
        else if ( group_flag == INTERIOR ) {
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key));
                for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                    if ( status0 != (*eit)->Status(prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key));
                for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                    for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        if ( status0 != (*eit)->Status(i,prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key));
                for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                    if ( status0 != (*nit)->Status(prop_key) ) return ANY;
                return status0;
            }
        }

    } // END SCALARS, ARRAYS


    if ( prop_key.type == VECTOR || prop_key.type == TENSOR || prop_key.type == FLAGGEDARRAY )
    {
        size_t length( prop_key.flagDepth );

        // if all the flags in the region shall be changed
        if ( group_flag == COMPLETE ) {
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key,component));
                for ( const auto& eit : cell_vec_ )
                    for ( auto n=0U; n<length; n++ )
                        if ( status0 != eit->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key,component));
                for ( const auto& eit : cell_vec_ )
                    for ( auto i=0U; i<eit->IntegrationPoints(); i++ )
                        for ( auto n=0U; n<length; n++ )
                            if ( status0 != eit->Status(i,prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key,component));
                for ( const auto& nit : node_vec_ )
                    for ( auto n=0U; n<length; n++ )
                        if ( status0 != nit->Status(prop_key,n) ) return ANY;
                return status0;
            }
        }
        // if the property shall only be changed on the Subdomain boundary
        else if ( group_flag == PERIMETER ) {
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(InteriorCells())->Status(prop_key,component));
                for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                    for ( auto n=0U; n<length; n++ )
                        if ( status0 != (*eit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(InteriorCells())->Status(0,prop_key,component));
                for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                    for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        for ( auto n=0U; n<length; n++ )
                            if ( status0 != (*eit)->Status(i,prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(InteriorNodes())->Status(prop_key,component));
                for ( auto nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                    for ( auto n=0U; n<length; n++ )
                        if ( status0 != (*nit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
        }
        else if ( group_flag == INTERIOR ) {
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key,component));
                for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                    for ( auto n=0U; n<length; n++ )
                        if ( status0 != (*eit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key,component));
                for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                    for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        for ( auto n=0U; n<length; n++ )
                            if ( status0 != (*eit)->Status(i,prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key,component));
                for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                    for ( auto n=0U; n<length; n++ )
                        if ( status0 != (*nit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
        }

    } // END VECTORS, TENSORS, FLAGGEDARRAYS

    return ANY;

} // end PropertyStatus




template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::ChangePropertyStatus( const char* property,
                                                     VARIABLE_FLAG new_status_of_scalar,
                                                     SUBDOMAIN_PART sdpart )
 {
    vector<VARIABLE_FLAG>  status(1U,new_status_of_scalar);
    ChangePropertyStatus( property, status, sdpart );
 }



template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::ChangePropertyStatusWhere( const char* property,
                                                          VARIABLE_FLAG new_status_of_scalar,
                                                          double min_value_to_change,
                                                          double max_value_to_change )
 {
    vector<VARIABLE_FLAG>  status(1U,new_status_of_scalar);
    ChangePropertyStatusWhere( property, status, min_value_to_change, max_value_to_change );
 }



/// changes the variable flag to status for those group members which carry the group_flag.
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::ChangePropertyStatus( const char* property,
                                                     const vector<VARIABLE_FLAG>& status,
                                                     SUBDOMAIN_PART group_flag )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const csmp::Index  prop_key = pref_.StorageKey(property);
    string src("ModelSubDomain<");
    src += to_string(dim);
    src += ",";
    src += "Element/Face/Interface";
    src +=">::ChangePropertyStatus:";

    if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
      throw csmp::Exception( ERROR, src.c_str(), "Use Status() to change flags of REGION/BOUNDARY/SPLIT_BUNDARY variables.");

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, src.c_str(), "Method not implemented for tensor properties yet");

    if ( status.empty() )
      throw csmp::Exception( ERROR, src.c_str(), "Status vector has not been initialized");
    else if ( status.size() != 1U && status.size() != dim ) {
        if ( prop_key.type == SCALAR && prop_key.place == NODE && status.size() != node_vec_.size() )
          throw csmp::Exception( ERROR, src.c_str(), "Status vector (SCALAR,NODE) has the wrong size");

        if ( prop_key.type == VECTOR && prop_key.place == NODE && status.size() != node_vec_.size()*dim )
          throw csmp::Exception( ERROR, src.c_str(), "Status vector (VECTOR,NODE) has the wrong size");
      }
      
    if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::ChangePropertyStatus",
                            property, "for InterFace nodes, node property needs to be assigned with SplitBoundary::ChangePropertyStatus");

    if ( cell_vec_.empty() ) {
         csmp_error.Note( ERROR, src.c_str(), "Region is empty.");
         return;
      }

    if ( prop_key.type == SCALAR || prop_key.type == ARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto& eit : cell_vec_ )
                      eit->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto& eit : cell_vec_  )
                      for ( auto i=0U; i<eit->IntegrationPoints(); i++ )
                        eit->Status( i, prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto& nit : node_vec_ )
                      nit->Status( prop_key, status[0U] );
                    return;
                 }
            }
          // if the flags shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); ++eit )
                      (*eit)->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); ++eit )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=PerimeterNodesBegin(); nit!=NodesEnd(); ++nit )
                      (*nit)->Status( prop_key, status[0U] );
                    return;
                 }
            }
          // if the flags for all entities in the interior of the region shall be changed
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); ++eit )
                      (*eit)->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); ++eit )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); ++nit )
                      (*nit)->Status( prop_key, status[0U] );
                    return;
                 }
            }

      } // END SCALARS, ARRAYS

    if ( prop_key.type == VECTOR || prop_key.type == TENSOR || prop_key.type == FLAGGEDARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto& eit : cell_vec_ )
                      for ( auto n=0U; n<status.size(); n++ )
                        eit->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto& eit : cell_vec_ )
                      for ( auto i=0U; i<eit->IntegrationPoints(); i++ )
                        for ( auto n=0U; n<status.size(); n++ )
                          eit->Status( i, prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto& nit : node_vec_ )
                      for ( auto n=0U; n<status.size(); n++ )
                        nit->Status( prop_key, n, status[n] );
                    return;
                 }
            }
          // if the property shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                      for ( auto n=0U; n<status.size(); n++ )
                        (*eit)->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        for ( auto n=0U; n<status.size(); n++ )
                          (*eit)->Status( i, prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                      for ( auto n=0U; n<status.size(); n++ )
                        (*nit)->Status( prop_key, n, status[n] );
                    return;
                 }
            }
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                      for ( auto n=0U; n<status.size(); n++ )
                        (*eit)->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        for ( auto n=0U; n<status.size(); n++ )
                          (*eit)->Status( i, prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                      for ( auto n=0U; n<status.size(); n++ )
                        (*nit)->Status( prop_key, n, status[n] );
                    return;
                 }
            }

      } // END VECTORS,TENSORS,FLAGGEDARRAYS

    cout <<"\n'"<< property <<"' ";
    throw csmp::Exception( ERROR, src.c_str(), "Property placement not recognized");

 } // end ChangePropertyStatus






/// changes the variable flag to status for those group members which carry the group_flag.
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::ChangePropertyStatus( const char* property,
                                                     uint32_t position,
                                                     VARIABLE_FLAG status,
                                                     SUBDOMAIN_PART group_flag )
 {
    const csmp::Index  prop_key = pref_.StorageKey(property);
    string src("ModelSubDomain<");
    string cache(to_string(dim));
    src += cache;
    src += ",";
    src += typeid(CELL<dim>).name();
    src +=">::ChangePropertyStatus";

    if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
      throw csmp::Exception( ERROR, src.c_str(), "Use Status() to change flags of REGION/BOUNDARY/SPLIT_BUNDARY variables.");

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, src.c_str(), "Method not implemented for tensor properties yet");

    if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::ChangePropertyStatus",
                            property, "for InterFace nodes, node property needs to be assigned with SplitBoundary::ChangePropertyStatus");

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( cell_vec_.empty() ) {
         csmp_error.Note( ERROR, src.c_str(), "Region is empty");
         return;
      }

    if ( prop_key.type == SCALAR || prop_key.type == ARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto& eit : cell_vec_  )
                      eit->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto& eit : cell_vec_ )
                      for ( auto i=0U; i<eit->IntegrationPoints(); i++ )
                        eit->Status( i, prop_key, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto& nit : node_vec_ )
                      nit->Status( prop_key, status );
                    return;
                 }
            }
          // if the property shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                      (*eit)->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                      (*nit)->Status( prop_key, status );
                    return;
                 }
            }
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                      (*eit)->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto& nit : node_vec_ )
                      nit->Status( prop_key, status );
                    return;
                 }
            }

      } // END SCALARS, ARRAYS

    if ( prop_key.type == VECTOR || prop_key.type == TENSOR || prop_key.type == FLAGGEDARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto& eit : cell_vec_ )
                        eit->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto& eit : cell_vec_ )
                      for ( auto i=0U; i<eit->IntegrationPoints(); i++ )
                          eit->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto& nit : node_vec_ )
                        nit->Status( prop_key, position, status );
                    return;
                 }
            }
          // if the property shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                        (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=PerimeterCellsBegin(); eit!=CellsEnd(); eit++ )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                          (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                        (*nit)->Status( prop_key, position, status );
                    return;
                 }
            }
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                        (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( auto eit=CellsBegin(); eit!=PerimeterCellsBegin(); eit++ )
                      for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                          (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( auto nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                        (*nit)->Status( prop_key, position, status );
                    return;
                 }
            }

      } // END VECTORS,TENSORS, FLAGGEDARRAYS

    cout <<"\n'"<< property <<"' ";
    throw csmp::Exception( ERROR, src.c_str(), "Property placement not recognized");

 } // end ChangePropertyStatus






template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::ChangePropertyStatusWhere( const char* property,
                                                          const std::vector<VARIABLE_FLAG>& status,
                                                          double pmin, double pmax )
 {
    const csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere",
                            "Method not implemented for tensor properties yet");
    if ( status.empty() )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere",
                            "Status vector has not been initialized");

    if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::ChangePropertyStatusWhere",
                            property, "for InterFace nodes, node property needs to be assigned with SplitBoundary::ChangePropertyStatusWhere");

    if ( cell_vec_.empty() )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere", "Region is empty");

    if ( prop_key.type == SCALAR )
      {
        ScalarVariable  sc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     eit->Status( prop_key, status[0U] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i{0U}; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, sc );
                     if ( sc.IsWithinRange( pmin, pmax ) )
                       eit->Status( i, prop_key, status[0U] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     nit->Status( prop_key, status[0U] );
                }
           }
       } // END SCALARS

    if ( prop_key.type == ARRAY )
      {
         ArrayVariable  av;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     eit->Status( prop_key, status[0U] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, av );
                     if ( av.IsWithinRange( pmin, pmax ) )
                       eit->Status( i, prop_key, status[0U] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     nit->Status( prop_key, status[0U] );
                }
           }
       } // END ARRAYS

    if ( prop_key.type == VECTOR )
      {
         VectorVariable<dim>  vc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                     for ( auto n=0U; n<status.size(); n++ )
                       eit->Status( prop_key, n, status[n] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, vc );
                     if ( vc.IsWithinRange( pmin, pmax ) )
                       for ( auto n=0U; n<status.size(); n++ )
                         eit->Status( i, prop_key, n, status[n] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                     for ( auto n=0U; n<status.size(); n++ )
                       nit->Status( prop_key, n, status[n] );
                }
           }
       } // END VECTORS

    if ( prop_key.type == TENSOR )
      {
         TensorVariable<dim>  ts;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                     for ( auto n=0U; n<status.size(); n++ )
                       eit->Status( prop_key, n, status[n] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i{0U}; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, ts );
                     if ( ts.IsWithinRange( pmin, pmax ) )
                       for ( auto n=0U; n<status.size(); n++ )
                         eit->Status( i, prop_key, n, status[n] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                     for ( auto n=0U; n<status.size(); n++ )
                       nit->Status( prop_key, n, status[n] );
                }
           }
       } // END TENSORS

    if ( prop_key.type == FLAGGEDARRAY )
      {
         FlaggedArrayVariable fv;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                     for ( auto n=0U; n<status.size(); n++ )
                       eit->Status( prop_key, n, status[n] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i{0U}; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, fv );
                     if ( fv.IsWithinRange( pmin, pmax ) )
                       for ( auto n=0U; n<status.size(); n++ )
                         eit->Status( i, prop_key, n, status[n] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                     for ( auto n=0U; n<status.size(); n++ )
                       nit->Status( prop_key, n, status[n] );
                }
           }
       } // END FLAGGEDARRAYS

 } // end ChangePropertyStatusWhere (if in value range)









template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::ChangePropertyStatusWhere( const char* property,
                                                          uint32_t position,
                                                          VARIABLE_FLAG status,
                                                          double pmin, double pmax )
 {
    const csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere",
                            "Method not implemented for tensor properties yet");

    if ( prop_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::ChangePropertyStatusWhere",
                            property, "for InterFace nodes, node property needs to be assigned with SplitBoundary::ChangePropertyStatusWhere");

    if ( cell_vec_.empty() )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere", "Region is empty");

    if ( prop_key.type == SCALAR )
      {
        ScalarVariable  sc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     eit->Status( prop_key, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, sc );
                     if ( sc.IsWithinRange( pmin, pmax ) )
                       eit->Status( i, prop_key, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     nit->Status( prop_key, status );
                }
           }
       } // END SCALARS

    if ( prop_key.type == ARRAY )
      {
         ArrayVariable  av;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     eit->Status( prop_key, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, av );
                     if ( av.IsWithinRange( pmin, pmax ) )
                       eit->Status( i, prop_key, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     nit->Status( prop_key, status );
                }
           }
       } // END ARRAYS

    if ( prop_key.type == VECTOR )
      {
         VectorVariable<dim>  vc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                       eit->Status( prop_key, position, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, vc );
                     if ( vc.IsWithinRange( pmin, pmax ) )
                         eit->Status( prop_key, position, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                       nit->Status( prop_key, position, status );
                }
           }
       } // END VECTORS

    if ( prop_key.type == TENSOR )
      {
         TensorVariable<dim>  ts;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                       eit->Status( prop_key, position, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, ts );
                     if ( ts.IsWithinRange( pmin, pmax ) )
                         eit->Status( prop_key, position, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                       nit->Status( prop_key, position, status );
                }
           }
       } // END TENSORS

    if ( prop_key.type == FLAGGEDARRAY )
      {
         FlaggedArrayVariable fv;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( auto& eit : cell_vec_ ) {
                   eit->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                       eit->Status( prop_key, position, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( auto& eit : cell_vec_ )
                for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) {
                     eit->Read( i, prop_key, fv );
                     if ( fv.IsWithinRange( pmin, pmax ) )
                         eit->Status( prop_key, position, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( auto& nit : node_vec_ ) {
                   nit->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                       nit->Status( prop_key, position, status );
                }
           }
       } // END VECTORS

 } // end ChangePropertyStatusWhere (if in value range)




















// INTERPOLATIONS AND EXTRAPOLATIONS


/**

InterpolateNodePropertyToCellProperty() interpolates node property to the
'barycenter' of the triangle. The resulting
value is different from the result obtained by applying the interrelation
subclass NodeToCellProperty. The Calculate() method of the latter assigns
the average value of the 3 cell nodes to the cell property 'eprop'.

@section arguments Input Arguments

The two string arguments define the name of the node property which is
interpolated and the name of the cell property to which the
interpolated value is written.

@section application Application

This method was developed for two-dimensional convection calculations
where the fluid density as cell property must be exactly the same for in
the two adjacent triangles which make up a square in regular-gridded meshes.

@section messages Messages

If the variable placement or type of the specified properties fails to
match the specifications outlined above,
InterpolateNodePropertyToCellProperty() will report an error and
return without completing its task.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::InterpolateNodeToCellProperty( const char* nprop, const char* eprop )
 {
     const csmp::Index  e_key = pref_.StorageKey(eprop),
                        n_key = pref_.StorageKey(nprop);

     // 1. check whether conditions for operation are O.K.
     if ( e_key.place != ELEMENT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToCellProperty",
                                 "Property arg2 is not an cell property, nothing was done...");
          return;
       }
     if ( n_key.place != NODE ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToCellProperty",
                                 "Property arg1 is not a node property, nothing was done...");
          return;
       }
     if ( e_key.type != n_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToCellProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

    if ( n_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::InterpolateNodeToCellProperty",
                              nprop, "for InterFace nodes, node property needs to be assigned with SplitBoundary::InterpolateNodeToCellProperty");

     switch( e_key.type )
       {
           case SCALAR: {
                ScalarVariable  sce;
                for ( auto& eit : cell_vec_ )
                  {
                     eit->PropertyValueAtBaryCenter( n_key, sce );
                     eit->Store( e_key, sce );
                  }
                }
             break;
           case VECTOR: {
                VectorVariable<dim>  vce;
                for ( auto& eit : cell_vec_ )
                  {
                     eit->PropertyValueAtBaryCenter( n_key, vce );
                     eit->Store( e_key, vce );
                  }
               }
            break;
          case TENSOR: {
               TensorVariable<dim>  tse;
               for ( auto& eit : cell_vec_ )
                 {
                     eit->PropertyValueAtBaryCenter( n_key, tse );
                     eit->Store( e_key, tse );
                 }
              }
            break;
          default:
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToCellProperty",
                                   "Property type could not be identified, nothing was done...");

     } // end switch

 } // end InterpolateNodeToCellProperty





template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* cprop )
 {
     const csmp::Index  c_key = pref_.StorageKey(cprop),
                        n_key = pref_.StorageKey(nprop);

     // 1. check whether conditions for operation are O.K.
     if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToIntegrationPointProperty",
                                 "Property arg2 is not an constraint point property, nothing was done...");
          return;
       }
     if ( n_key.place != NODE ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToIntegrationPointProperty",
                                 "Property arg1 is not a node property, nothing was done...");
          return;
       }
     if ( c_key.type != n_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToIntegrationPointProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

    if ( n_key.place == NODE && is_same<CELL<dim>,InterFace<dim>>::value )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,InterFace>::InterpolateNodeToIntegrationPointProperty",
                              nprop, "for InterFace nodes, node property needs to be assigned with SplitBoundary::InterpolateNodeToIntegrationPointProperty");

     vector<double>  IPOL;

     switch( c_key.type )
       {
           case SCALAR: {
                for ( auto& eit : cell_vec_ )
                  for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                    {
                       eit->N_AtIntegrationPoint( i, IPOL );
                       double sc = IPOL[0] * eit->N(0)->Read( n_key );
                       for ( auto j{1U}; j<eit->Nodes(); j++ ) sc += IPOL[j] * eit->N(j)->Read( n_key );
                       eit->Store( i, c_key, makeScalar(eit->Status(i,c_key),sc) );
                    }
                }
             break;
           case VECTOR: {
                VectorVariable<dim>  vce, vce2;
                for ( auto& eit : cell_vec_ )
                  for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                    {
                       eit->N_AtIntegrationPoint( i, IPOL );
                       vce=0.;
                       for ( auto j{0U}; j<eit->Nodes(); j++ )
                         {
                            eit->N(j)->Read( n_key, vce2 );
                            vce += (vce2 * IPOL[j]);
                         }
                       eit->Store( i, c_key, vce );
                    }
                }
            break;
          case TENSOR: {
               TensorVariable<dim>  tse, tse2;
                for ( auto& eit : cell_vec_ )
                  for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                    {
                       eit->N_AtIntegrationPoint( i, IPOL );
                       tse=0.;
                       for ( auto j{0U}; j<eit->Nodes(); j++ )
                         {
                            eit->N(j)->Read( n_key, tse2 );
                            tse += (tse2 * IPOL[j]);
                         }
                       eit->Store( i, c_key, tse );
                    }
                }
            break;
          default:
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateNodeToIntegrationPointProperty",
                                   "Property type could not be identified, nothing was done...");

     } // end switch

 } // end InterpolateNodeToIntegrationPointProperty





/**

This method interpolates the specified constraint point property
to the target cell property.

@param cprop The name of the constraint point property which shall be interpolated to the
@param eprop cell property (arg2).

The results are returned to the Model.

@section implementation Implementation

Currently the method averages the IntegrationPoint variable values
to find the cell property value.

@section application Application

To visualise constraint point properties as CELL_CENTERED variables in
VTK they have to have a unique value in the cell = CELL.

Alternatively, one could integrate the property over the cell and
divide the result value by the cell area. However, this would work
only for scalar properties.

@section messages Messages

Consistency checks are performed on the placement and type of the
input variables.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::InterpolateIntegrationPointToCellProperty( const char* cprop, const char* eprop )
 {
     const csmp::Index  e_key = pref_.StorageKey(eprop);
     const csmp::Index  c_key = pref_.StorageKey(cprop);

     if ( e_key.place != ELEMENT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateIntegrationPointToCellProperty",
                                 "Property arg2 is not an cell property, nothing was done...");
          return;
       }
     if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateIntegrationPointToCellProperty",
                                 "Property arg1 is not a constraint point property, nothing was done...");
          return;
       }
     if ( e_key.type != c_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::InterpolateIntegrationPointToCellProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

    if ( e_key.type == SCALAR ) {
        ScalarVariable           sc;
          vector<ScalarVariable >  sc_vec;
          for ( auto& it : cell_vec_ ) {
              it->IntegrationPointPropertyVector( c_key, sc_vec );
              average( sc_vec, sc );
              it->Store( e_key, sc );
            }
      }
    else if ( e_key.type == VECTOR ) {
        VectorVariable<dim>           vc;
        vector<VectorVariable<dim> >  vc_vec;
          for ( auto& it : cell_vec_ ) {
              it->IntegrationPointPropertyVector( c_key, vc_vec );
              average( vc_vec, vc );
              it->Store( e_key, vc );
            }
      }
    if ( e_key.type == TENSOR ) {
        TensorVariable<dim>           ts;
        vector<TensorVariable<dim> >  ts_vec;
          for ( auto& it : cell_vec_ ) {
              it->IntegrationPointPropertyVector( c_key, ts_vec );
              average( ts_vec, ts );
              it->Store( e_key, ts );
            }
      }

 } // end InterpolateIntegrationPointToCellProperty

template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::ExtrapolateCellToIntegrationPointProperty( const char* eprop,
                                                                           const char* cprop )
{
    assert( !cell_vec_.empty() );

    const csmp::Index  e_key = pref_.StorageKey(eprop);
    const csmp::Index  c_key = pref_.StorageKey(cprop);

    // 1. check whether conditions for operation are O.K.
    if ( e_key.place != ELEMENT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToIntegrationPointProperty",
                                "Property arg1 is not an cell property, nothing was done...");
         return;
      }
    if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToIntegrationPointProperty",
                                "Property arg2 is not an integration point property, nothing was done...");
         return;
      }
    if ( e_key.type != c_key.type ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToIntegrationPointProperty",
                                "Properties are not of the same type, nothing was done...");
         return;
      }

    // 2. For each cell the property value is assigned for each integration point inside that cell

    // SCALAR VARIABLES
    // ----------------
    if ( e_key.type == SCALAR )
    {
         ScalarVariable sc;
         for ( auto& eit : cell_vec_ )
             {
                 // reading cell property
                 sc = eit->Read( e_key );

                 for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                     eit->Store( i, c_key, sc );
              }
    }

    // VECTOR VARIABLES
    // ----------------
    if ( e_key.type == VECTOR )
    {
         VectorVariable<dim>  vc;

         for ( auto& eit : cell_vec_ )
             {
                 // reading cell property
                 eit->Read( e_key, vc );

                 for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                     eit->Store( i, c_key, vc );
             }
    }

    // TENSOR VARIABLES
    // ----------------
    if ( e_key.type == TENSOR )
    {
         TensorVariable<dim>  ts;

         for ( auto& eit : cell_vec_ )
             {
                 eit->Read( e_key, ts );

                 for ( auto i{0U}; i<eit->IntegrationPoints(); i++ )
                     eit->Store( i, c_key, ts );
             }
    }

} // end ExtrapolateCellToIntegrationPointProperty







template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::ExtrapolateCellToFacetIntegrationPointProperty( const char* eprop,
                                                                                const char* fipprop )
{
    assert( !cell_vec_.empty() );

    const csmp::Index  e_key = pref_.StorageKey(eprop);
    const csmp::Index  fip_key = pref_.StorageKey(fipprop);

    // 1. check whether conditions for operation are O.K.
    if ( e_key.place != ELEMENT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToFacetIntegrationPointProperty",
                                "Property arg1 is not an cell property, nothing was done...");
         return;
      }
    if ( fip_key.place != FACET_INTEGRATION_POINT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToFacetIntegrationPointProperty",
                                "Property arg2 is not an integration point property, nothing was done...");
         return;
      }
    if ( e_key.type != fip_key.type ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToFacetIntegrationPointProperty",
                                "Properties are not of the same type, nothing was done...");
         return;
      }

    // 2. For each cell the property value is assigned for each facet integration point inside that cell

    // SCALAR VARIABLES
    // ----------------
    if ( e_key.type == SCALAR )
    {
         ScalarVariable sc;
         for ( auto& eit : cell_vec_ )
             {
                 // reading cell property
                 sc = eit->Read( e_key );

                 for ( auto i{0U}; i<eit->FV()->Facets(); i++ )
                   for ( auto j{0U}; j<eit->FV()->IntegrationPointsPerFacet(); j++ )
                     eit->Store( i, j, fip_key, sc );
              }
    }

    // VECTOR VARIABLES
    // ----------------
    if ( e_key.type == VECTOR )
    {
         VectorVariable<dim>  vc;

         for ( auto& eit : cell_vec_ )
             {
                 // reading cell property
                 eit->Read( e_key, vc );

                 for ( auto i{0U}; i<eit->FV()->Facets(); i++ )
                   for ( auto j{0U}; j<eit->IntegrationPointsPerFacet(); j++ )
                     eit->Store( i, j, fip_key, vc );
             }
    }

    // TENSOR VARIABLES
    // ----------------
    if ( e_key.type == TENSOR )
    {
         TensorVariable<dim>  ts;

         for ( auto& eit : cell_vec_ )
             {
                 eit->Read( e_key, ts );

                 for ( auto i{0U}; i<eit->FV()->Facets(); i++ )
                    for ( auto j{0U}; j<eit->FV()->IntegrationPointsPerFacet(); j++ )
                      eit->Store( i, j, fip_key, ts );
             }
    }

} // end ExtrapolateCellToFacetIntegrationPointProperty




/**
    Interpolates piecewise constant cell properties to the nodes.
    Since each node belongs to multiple cells their contributions to this node have to be weighted.
    This method offers two ways to do this.
    
    'by-distance' is the default. In this approach, contributions are weighted by the inverse of the 
    distance of each cell's barycentre from the node. Thus, small cells have a bigger role than
    big ones.
    
    'default' the contributions are weighted by cell length/area/volume. This weighs bigger cells
    more than small ones. In many cases this is not what one wants, but it is still offered for
    some specific applications.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::ExtrapolateCellToNodeProperty( const char* eprop,
                                                               const char* nprop,
                                                               bool  by_distance )
 {
     assert( !node_vec_.empty() );
     assert( !cell_vec_.empty() );

     const csmp::Index  e_key = pref_.StorageKey(eprop),
                        n_key = pref_.StorageKey(nprop);

     // 1. check whether conditions for operation are O.K.
     if ( e_key.place != ELEMENT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToNodeProperty",
                                 "Property arg1 is not an cell property, nothing was done...");
          return;
       }
     if ( n_key.place != NODE ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToNodeProperty",
                                 "Property arg2 is not a node property, nothing was done...");
          return;
       }
     if ( e_key.type != n_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateCellToNodeProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

     // 2. The cell property values are weighted by the area of the cell and then
     //    added for each node to a property array vector, the number of additions to each
     //    node is counted
     vector<double>  sum( Nodes(), 0. );

     RenumberNodes();

     // SCALAR VARIABLES
     // ----------------
     if ( e_key.type == SCALAR ) {
          vector<double>  sc_data( Nodes(), 0. );
          if ( !by_distance )
            {
               for ( auto& eit : cell_vec_ )
                 {
                    // getting the properties
                    double sc = eit->Read( e_key );

                    // weighting property value by cell volume
                    double volume = eit->Volume();
                    sc *= volume;

                    for ( auto i{0U}; i<eit->Nodes(); i++ ) {
                         sc_data[ eit->N(i)->Idx() ] += sc;
                         sum[ eit->N(i)->Idx() ]     += volume;
                      }
                 }
            }
          else
            {
                for ( auto& eit : cell_vec_ )
                  {
                     // reading cell property
                     double sc = eit->Read( e_key );
                     // measuring distances from cell to nodes to
                     // weight property values
                     Point<dim>  ctr(eit->BaryCenter());

                     for ( auto i=0U; i<eit->Nodes(); i++ ) {
                          double dist((ctr - eit->N(i)->Coordinate()).Length());
                          // memorizing distances and weighted values
                          sc_data[ eit->N(i)->Idx() ] += sc / dist;
                          sum[ eit->N(i)->Idx() ]     += 1. / dist;
                       }
                  }
            }

          // calculating the nodal averages and mapping them back to the nodes
          for ( auto& nit : node_vec_ )
            {
               ScalarVariable res( nit->Status( n_key ), sc_data[ nit->Idx() ] / sum[ nit->Idx() ] );
               nit->Store( n_key, res );
            }
       }



     // VECTOR VARIABLES
     // ----------------
     if ( e_key.type == VECTOR ) {
          VectorVariable<dim>  vc; vc = 0.;
          vector<VectorVariable<dim> >  vc_data(Nodes(), vc );
          if ( !by_distance )
            {
               for ( auto& eit : cell_vec_ )
                 {
                    eit->Read( e_key, vc );
                    double  volume = eit->Volume();
                    vc *= volume;

                    for ( auto i{0U}; i<eit->Nodes(); i++ )
                      {
                         vc_data[ eit->N(i)->Idx() ] += vc;
                         sum[ eit->N(i)->Idx() ] += volume;
                      }
                 }
            }
          else
            {
                for ( auto& eit : cell_vec_ )
                  {
                     // reading cell property
                     eit->Read( e_key, vc );
                     // measuring distances from cell to nodes to
                     // weight property values
                     Point<dim>  ctr(eit->BaryCenter());

                     for ( auto i=0U; i<eit->Nodes(); i++ ) {
                          double dist((ctr - eit->N(i)->Coordinate()).Length());
                          // memorizing distances and weighted values
                          vc_data[ eit->N(i)->Idx() ] += vc / dist;
                          sum[ eit->N(i)->Idx() ]     += 1. / dist;
                       }
                  }
            }
          for ( auto& nit : node_vec_ )
            {
               vc = vc_data[ nit->Idx() ] / sum[ nit->Idx() ];
               nit->Store( n_key, vc );
            }
       }



     // TENSOR VARIABLES
     // ----------------
     if ( e_key.type == TENSOR ) {
          TensorVariable<dim>  ts; ts = 0.;
          vector<TensorVariable<dim> >  ts_data(Nodes(), ts );
          if ( !by_distance )
            {
               for ( auto& eit : cell_vec_ )
                 {
                    eit->Read( e_key, ts );
                    double volume = eit->Volume();
                    ts *= volume;

                    for ( auto i{0U}; i<eit->Nodes(); i++ )
                      {
                         ts_data[ eit->N(i)->Idx() ] += ts;
                         sum[ eit->N(i)->Idx() ] += volume;
                      }
                 }
            }
          else
            {
                for ( auto& eit : cell_vec_ )
                  {
                     // reading cell property
                     eit->Read( e_key, ts );
                     // measuring distances from cell to nodes to
                     // weight property values
                     Point<dim>  ctr(eit->BaryCenter());

                     for ( auto i=0U; i<eit->Nodes(); i++ ) {
                          double dist((ctr - eit->N(i)->Coordinate()).Length());
                          // memorizing distances and weighted values
                          ts_data[ eit->N(i)->Idx() ] += ts / dist;
                          sum[ eit->N(i)->Idx() ]     += 1. / dist;
                       }
                  }
            }
          for ( auto& nit : node_vec_ )
            {
               ts = ts_data[ nit->Idx() ] / sum[ nit->Idx() ];
               nit->Store( n_key, ts );
            }
       }
    if ( verbose_ ) {
        cout <<"\nModel<"<<dim<<">::ExtrapolateCellToNodeProperty: ";
        cout <<"'" << eprop <<"' has been successfully extrapolated to '"<< nprop <<"'." << endl;
    }

 } // end ExtrapolateCellToNodeProperty



 

/**

Using the cell interpolation functions, this method extrapolates the
desired integration point property to the nodes. Constributions of adjacent
cells are averaged but no weighting by cell size or proximity of
barycentre to the node is applied.

@section arguments Input Arguments

The names of the targeted integration point and node variables.

The result of the extrapolation is returned into the Model property
storage.

@section implementation Implementation

The method depends on a corresponding function of the FiniteElement
which is for all numerically integrated cell types.
For cells with quadratic inpterpolation functions, the method
uses a linear extrapolation which is justified for solution variable derivatives
as they define a trilinear field on the IntegrationPoints.

@section application Application

Use this for numerically integrated finite cells, in particular ones with
an interpolation order greater than 1.

@section messages Messages

A consistency check on variable type and placement is performed.
 */
template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop )
 {
     const csmp::Index  n_key = pref_.StorageKey(nprop), c_key = pref_.StorageKey(cprop);

     ErrorHandler& csmp_error( ErrorHandler::Instance() );

     if ( n_key.place != NODE ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateIntegrationPointToNodeProperty",
                                 "Property arg2 is not a node property, nothing was done...");
          return;
       }
     if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateIntegrationPointToNodeProperty",
                                 "Property arg1 is not a constraint point property, nothing was done...");
          return;
       }
     if ( n_key.type != c_key.type ) {
          csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::ExtrapolateIntegrationPointToNodeProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

    RenumberNodes();

    // required vectors
    vector<double>  cp_var, n_var;

    // SCALAR PROPERTIES
    if ( n_key.type == SCALAR ) {
         vector<double> temp( Nodes(), 0. ), sum( Nodes(), 0. );
         const size_t  variable_components(1U);
         for ( auto& eit : cell_vec_ ) {
              // getting the integration-point property
              cp_var.resize( eit->IntegrationPoints() );
              for ( auto i=0U; i<eit->IntegrationPoints(); i++ ) cp_var[i] = eit->Read( i, c_key );
              // extrapolating it to nodes
              n_var.resize( eit->Nodes() );                 // components = 1
              eit->ExtrapolateIntegrationPointVariableToNodes( variable_components, cp_var, n_var );
              // weighting and accumulating it into vector<double> for later averaging
              const Point<dim>  barycenter(eit->BaryCenter());
              for ( auto i=0U; i<eit->Nodes(); i++ )
                {
                   // finding the distance of the node from the barycentre
                   const double distance = barycenter.DistanceTo(eit->N(i)->Coordinate());
                   temp[ eit->N(i)->Idx() ] += n_var[i] / distance;
                   sum[  eit->N(i)->Idx() ] += 1. / distance;
                }
           }
          // distance weighting and averaging the extrapolated node values
         for ( auto& nit : node_vec_ )
           nit->Store( n_key, makeScalar( nit->Status(n_key), temp[ nit->Idx() ] / sum[ nit->Idx() ]) );
      }

    // VECTOR PROPERTIES
    else if ( n_key.type == VECTOR ) {
         // creating a zero-initialized temporary vector<double>
         VectorVariable<dim>  zero_vec; zero_vec=0.;
         vector<VectorVariable<dim> > cpvec, temp( Nodes(), zero_vec ), sum( Nodes(), zero_vec );
         const auto  vcomponents{dim};
         for ( auto& it : cell_vec_ ) {
              // getting the constraint point property
              cpvec.resize( it->IntegrationPoints() );
              for ( auto i{0U}; i<it->IntegrationPoints(); i++ )
                it->Read( i, c_key, cpvec[i] );
              // rolling the vector<double> variables out into linear vector<double> cp_var
              cp_var.resize( it->IntegrationPoints() * vcomponents );
              for ( auto i{0U}; i<it->IntegrationPoints(); i++ )
                for ( auto j{0U}; j<vcomponents; j++ ) cp_var[ i * vcomponents + j ] = cpvec[i][j];
              // extrapolating constraint point property to nodes
              n_var.resize( it->Nodes() * vcomponents ); // components = 1
              it->ExtrapolateIntegrationPointVariableToNodes( vcomponents, cp_var, n_var );

              // weighting and accumulating it into vector<double> for later averaging
              const Point<dim>  barycenter(it->BaryCenter());
              // accumulating result into vector<double> for later averaging
              for ( auto i=0U; i<it->Nodes(); i++ )
                {
                   // finding the distance of the node from the barycentre
                   const double distance = barycenter.DistanceTo(it->N(i)->Coordinate());
                   for ( auto j{0U}; j<vcomponents; j++ ) {
                         temp[ it->N(i)->Idx() ](j) += n_var[ i * vcomponents + j ] / distance;
                         sum[  it->N(i)->Idx() ](j) += 1. / distance;
                     }
                }
           }
          // distance weighting and averaging the extrapolated node values
         for ( auto& it : node_vec_ )
           {
              // averaging the vector<double> variable
              for ( auto j{0U}; j<vcomponents; j++ ) {
                   temp[ it->Idx() ](j) /= sum[ it->Idx() ](j);
                   temp[ it->Idx() ].Flag(j) = it->Status( n_key, j );
                }
              // storing it
              it->Store( n_key, temp[ it->Idx() ] );
           }
      }

    // TENSOR PROPERTIES
    if ( n_key.type == TENSOR ) {
          csmp_error.Note( WARNING, "ModelSubDomain<dim,CELL>::ExtrapolateIntegrationPointToNodeProperty",
                            "distance weighting is not applied; tensor values are simply averaged at the nodes.");
         // creating a zero-initialized temporary vector<double>
         TensorVariable<dim>  zero_ts; zero_ts=0.;
         vector<TensorVariable<dim> > cpts, temp( Nodes(), zero_ts );
         const size_t  tcomponents(dim * dim);
         for ( auto& it : cell_vec_ ) {
              // getting the constraint point property
              cpts.resize( it->IntegrationPoints() );
              for ( auto i=0U; i<it->IntegrationPoints(); i++ )
                it->Read( i, c_key, cpts[i] );
              // rolling the tensor variable rows out sequentially into the linear vector<double> cp_var
              cp_var.resize( it->IntegrationPoints() * tcomponents );
              for ( auto i{0U}; i<it->IntegrationPoints(); i++ )
                for ( auto j{0U}; j<dim; j++ )
                  for ( auto k{0U}; k<dim; k++ ) cp_var[ i * tcomponents + j * dim + k ] = cpts[i](j,k);
              // extrapolating constraint point property to nodes
              n_var.resize( it->Nodes() * tcomponents ); // components = 1
              it->ExtrapolateIntegrationPointVariableToNodes( tcomponents, cp_var, n_var );
              // accumulating result into tensor variable vector<double> for later averaging
              for ( auto i=0U; i<it->Nodes(); i++ )
                for ( auto j{0U}; j<dim; j++ )
                  for ( auto k{0U}; k<dim; k++ )
                    temp[ it->N(i)->Idx() ](j,k) += n_var[ i * tcomponents + j * dim + k ];
           }
         // averaging the resulting node property and storing it
         for ( auto& it : node_vec_ ) {
              // averaging the tensor variable
              for ( auto j{0U}; j<dim; j++ )
                for ( auto k{0U}; k<dim; k++ )
                  temp[ it->Idx() ](j,k) /= static_cast<double>( it->Parents() );
              // storing it
              it->Store( n_key, temp[ it->Idx() ] );
           }
      }

 } // end ExtrapolateIntegrationPointToNodeProperty








// CALCULATIONS


/**

Takes the arithmetic (not cell size weighted) average of all property values of nodes, integration points,
or cells that belong to the region, depending on where the target property is placed.

@param prop The name of the property which shall be averaged.
*/
template<uint32_t dim, template<uint32_t> class CELL>
double  ModelSubDomain<dim,CELL>::Average( const char* prop ) const
 {
    const csmp::Index  idx = pref_.StorageKey(prop);
    size_t             counter(0U);

     if ( (cell_vec_.empty()) ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::RegionPropertyAverage",
                                        "Region is empty");
       }

    if ( idx.place == REGION or idx.place == BOUNDARY or idx.place == SPLIT_BOUNDARY )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::Average", "REGION/BOUNDARY/SPLIT_BOUNDARY variables have only a single value in a model subdomain.");

    switch( idx.place )
      {
        case ELEMENT:
        case FACE:
        case INTER_FACE:
            if ( idx.type == SCALAR ) {
                 ScalarVariable  sc;
                 double  avg(0.);
                 for ( typename vector<CELL<dim>*>::const_iterator
                       eit=cell_vec_.begin(); eit!=cell_vec_.end(); eit++ ) {
                      (*eit)->Read( idx, sc );
                      avg += sc();
                   }
                 return avg / static_cast<double>(cell_vec_.size());
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 double  avg(0.);
                 for ( typename vector<CELL<dim>*>::const_iterator
                       eit=cell_vec_.begin(); eit!=cell_vec_.end(); eit++ ) {
                      (*eit)->Read( idx, vc );
                      avg += vc.Length();
                   }
                 return avg /= static_cast<double>(cell_vec_.size());
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 double  avg(0.);
                 for ( typename vector<CELL<dim>*>::const_iterator
                       eit=cell_vec_.begin(); eit!=cell_vec_.end(); eit++ ) {
                      (*eit)->Read( idx, ts );
                      ts.EigenValues( evals );
                      double ts_avg(0.);
                      for ( auto j{0U}; j<dim; j++ ) ts_avg += evals[j];
                      avg += ts_avg / static_cast<double>(dim);
                   }
                 return avg / static_cast<double>(cell_vec_.size());
              }
          break;
        case ELEMENT_INTEGRATION_POINT:
            assert( pref_.VariableCount(ELEMENT_INTEGRATION_POINT) > 0U );
            counter = 0U;
            if ( idx.type == SCALAR ) {
                 ScalarVariable  sc;
                 double  avg(0.);
                 for ( typename vector<CELL<dim>*>::const_iterator
                       eit=cell_vec_.begin(); eit!=cell_vec_.end(); eit++ )
                   for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ ) {
                        (*eit)->Read( i, idx, sc );
                        counter++;
                        avg += sc();
                     }
                 return avg / static_cast<double>(counter);
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 double  avg(0.);
                 for ( typename vector<CELL<dim>*>::const_iterator
                       eit=cell_vec_.begin(); eit!=cell_vec_.end(); eit++ )
                   for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ ) {
                        (*eit)->Read( i, idx, vc );
                        counter++;
                        avg += vc.Length();
                     }
                 return avg /= static_cast<double>(counter);
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 double  avg(0.);
                 for ( typename vector<CELL<dim>*>::const_iterator
                       eit=cell_vec_.begin(); eit!=cell_vec_.end(); eit++ )
                   for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ ) {
                        (*eit)->Read( i, idx, ts );
                        counter++;
                        ts.EigenValues( evals );
                        double ts_avg(0.);
                        for ( auto j{0U}; j<dim; j++ ) ts_avg += evals[j];
                        avg += ts_avg / static_cast<double>(dim);
                     }
                 return avg / static_cast<double>(counter);
              }
          break;
        case NODE:
            if ( idx.type == SCALAR ) {
                 ScalarVariable  sc;
                 double  avg(0.);
                 for ( typename vector<csmp::Node<dim>*>::const_iterator
                       it=node_vec_.begin(); it!=node_vec_.end(); it++ ) {
                      (*it)->Read( idx, sc );
                      avg += sc();
                   }
                 return avg / static_cast<double>(node_vec_.size());
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 double  avg(0.);
                 for ( typename vector<csmp::Node<dim>*>::const_iterator
                       it=node_vec_.begin(); it!=node_vec_.end(); it++ ) {
                      (*it)->Read( idx, vc );
                      avg += vc.Length();
                   }
                 return avg /= static_cast<double>(node_vec_.size());
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 double  avg(0.);
                 for ( typename vector<csmp::Node<dim>*>::const_iterator
                       it=node_vec_.begin(); it!=node_vec_.end(); it++ ) {
                      (*it)->Read( idx, ts );
                      ts.EigenValues( evals );
                      double ts_avg(0.);
                      for ( auto j{0U}; j<dim; j++ ) ts_avg += evals[j];
                      avg += ts_avg / static_cast<double>(dim);
                   }
                 return avg / static_cast<double>(node_vec_.size());
              }
         break;
         default:
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::RegionAverage:", "property placement not handled yet.");
       }

    return std::numeric_limits<double>::signaling_NaN();

 } // end Average




/**

Provided that property 'a' is a scalar node variable and property 'b' is a
vector<double> variable placed on the cell, CopyGradientOfProperty_A_To_B() will
calculate the gradient of property 'a' for each cell and assign the
result to the cell variable 'b'.

@section arguments Input Arguments

The method takes two name strings as arguments. The first string specifies
the scalar node variable of which the gradient will be calculated and the
second string identifies the vector<double> variable placed on the cell which
will store the calculated gradient of property 'a'.

@section application Application

The definition of many geophysical or geochemical interrelations requires
a knowledge of property gradients which can be calculated with this method,
using the cell interpolation functions of type of finite cell which
was used to build the mesh. The method can also be used for post-processing
following the application of algorithms.

@section messages Messages

CopyGradientOfProperty_A_To_B() will return without completing any
calculations, if the input variables do not comply with the specifications
outlined above.

*/
template<uint32_t dim, template<uint32_t> class CELL>
bool  ModelSubDomain<dim,CELL>::CopyGradientOfProperty_A_To_B( const char* a, const char* b )
 {
    const csmp::Index  a_key = pref_.StorageKey(a), b_key = pref_.StorageKey(b);

    // 1. testing variable A for suitability
    // -------------------------------------
    if ( !pref_.IsDefined(a) ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property A does not exist. Now exciting..." << endl;
         return false;
      }
    if ( a_key.place != NODE ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property A is not a NODE variable. Gradient can't be calculated..." << endl;
         return false;
      }
    if ( a_key.type == TENSOR ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property A is a TENSOR variable. Nothing is done..." << endl;
         return false;
      }

    // 2. testing variable B for suitability
    // -------------------------------------
    if ( !pref_.IsDefined(b) ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property B does not exist. Now exciting..." << endl;
         return false;
      }
    if ( b_key.place != ELEMENT and b_key.place != ELEMENT_INTEGRATION_POINT ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property B is not an ELEMENT or ELEMENT_INTEGRATION_POINT variable. Gradient can't be calculated..." << endl;
         return false;
      }

    DenseMatrix<DM_MIN>  DN;

    // 3. SCALAR to VECTOR calculations
    // --------------------------------
    if ( a_key.type == SCALAR )
      {
         if ( b_key.type != VECTOR ) {
               cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
               cerr << "Property B is not a VECTOR variable. Nothing is done..." << endl;
               return false;
            }

         vector<ScalarVariable >  SC;
         VectorVariable<dim>      vc;
 
        if ( b_key.place == ELEMENT )
          {
             for ( auto& eit : cell_vec_ )
               {
                  if ( eit->FE()->UsesLocalCoordinates() ) eit->dN_AtBaryCenter( DN );
                  else eit->dN( DN );
                  eit->NodePropertyVector( a_key, SC );

                  vc = 0.;
                  const auto nodes(eit->Nodes());
                  for ( auto i=0U; i<nodes; i++ )
                    for ( auto j{0U}; j<dim; j++ )
                      vc(j) += DN(j,i) * SC[i]();
     
                  eit->Store( b_key, vc );
               }
             return true;
           }
        if ( b_key.place == ELEMENT_INTEGRATION_POINT )
          {
             for ( auto& eit : cell_vec_ )
               {
                  eit->NodePropertyVector( a_key, SC );

                  const auto ipoints(eit->IntegrationPoints());
                  const auto nodes(eit->Nodes());
                  for ( auto i=0U; i<ipoints; ++i )
                    {
                       eit->dN_AtIntegrationPoint( DN, i );
                       vc = 0.;
                       for ( auto l{0}; l<nodes; l++ )
                         for ( auto j{0U}; j<dim; j++ )
                           vc(j) += DN(j,l) * SC[i]();

                       eit->Store( i, b_key, vc );
                    }
               }
             return true;
           }
       } // end a.type=scalar
   
   
    // 4. VECTOR to TENSOR calculations
    // --------------------------------
    if ( a_key.type == VECTOR )
      {
         if ( b_key.type != TENSOR ) {
              cerr << "ModelSubDomain<"<<dim<<">::CopyGradientOfProperty_A_To_B: "<< endl;
              cerr << "Property A is a VECTOR so Property B should be a TENSOR variable. Nothing is done..." << endl;
              return false;
           }

         vector<VectorVariable<dim> > VC;
         TensorVariable<dim>          ts;

         if ( b_key.place == ELEMENT )
           {
             for ( auto& eit : cell_vec_ )
               {
                  if ( eit->FE()->UsesLocalCoordinates() ) eit->dN_AtBaryCenter( DN );
                  else eit->dN( DN );
                  eit->NodePropertyVector( a_key, VC );

                  ts = 0.;

                  // the gradients become rows of the tensor
                  const auto nodes(eit->Nodes());
                  for ( auto n=0U; n<nodes; n++ )
                    for ( auto i=0U; i<dim; i++ )
                      for ( auto j{0U}; j<dim; j++ )
                        ts(i,j) += DN(i,n) * VC[n][j];

                  // saving the resulting vector<double>
                  eit->Store( b_key, ts );
               }
             return true;
           }
         if ( b_key.place == ELEMENT_INTEGRATION_POINT )
           {
             for ( auto& eit : cell_vec_ )
               {
                  eit->NodePropertyVector( a_key, VC );
                  // the gradients become rows of the tensor
                  const size_t ipoints(eit->IntegrationPoints());
                  const size_t nodes(eit->Nodes());
                  for ( auto i{0U}; i<ipoints; ++i )
                    {
                      eit->dN_AtIntegrationPoint( DN, i );
                      ts = 0.;
                      for ( auto n=0U; n<nodes; n++ )
                        for ( auto l{0}; l<dim; l++ )
                          for ( auto j{0U}; j<dim; j++ )
                            ts(i,j) += DN(l,n) * VC[n][j];

                      // saving the resulting vector<double>
                      eit->Store( i, b_key, ts );
                    }
               }
             return true;
           }


      }

   return true;

 }  // end CopyGradientOfProperty_A_To_B
 

 
 

template<uint32_t dim, template<uint32_t> class CELL>
void  ModelSubDomain<dim,CELL>::CopyReplace( const char* from, const char* to )
 {
     if ( pref_.Type(from) == SCALAR ) {
          CopyReplaceVisitor<ScalarVariable,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == VECTOR ) {
          CopyReplaceVisitor<VectorVariable<dim>,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == TENSOR ) {
          CopyReplaceVisitor<TensorVariable<dim>,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == ARRAY ) {
          CopyReplaceVisitor<ArrayVariable,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == FLAGGEDARRAY ) {
          CopyReplaceVisitor<FlaggedArrayVariable,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else throw csmp::Exception( ERROR,
                                 "ModelSubDomain<dim,CELL>::CopyReplace",
                                 "Property type not supported yet");

 } // end CopyReplace




















// SCREEN OUTPUT



/**

Prints the values and flags of the distributed variable of interest on
stdout.

@param prop The name of the physical variable of interest.

@section application Application

To check variable values in small test problems.

@section messages Messages

OutputVariableToScreen() will report the variable type, its placement,
and property index. Then it will print the host object IDs followed
by the variable flags and values.

*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::OutputVariableToScreen( const char* prop ) const
 {
     const csmp::Index  prop_key = pref_.StorageKey(prop);

     if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
       throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::OutputVariableToScreen",
                             "Use direct variable access to display Region/Boundary/SplitBoundary properties.");

     string place_name = parsePlacement(prop_key.place);

     cout <<"\n\nModelSubDomain<"<< dim << ">::OutputVariableToScreen: output of ";
     switch(prop_key.type)
      {
         case SCALAR: cout <<"SCALAR variable: '";
           break;
         case VECTOR: cout <<"VECTOR variable: '";
           break;
         case TENSOR: cout <<"TENSOR variable: '";
           break;
         case ARRAY:  cout <<"ARRAY variable: '";
           break;
         case FLAGGEDARRAY:  cout <<"FLAGGEDARRAY variable: '";
           break;
         default:
           cout <<"\nModelSubDomain<"<< dim <<">::OutputVariableToScreen: type of '"<< prop <<"' not recognized."<< endl;
           return;
      }
     cout << prop <<"' placed on the: "<< place_name << endl;
     
     VectorVariable<dim> vc;
     TensorVariable<dim> ts;

     switch( prop_key.place )
       {
          case NODE:
              for ( auto nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                {
                  cout <<"\nIdx: " << (*nit)->Idx() <<"\t\t";
                  cout << string(parseBoundary((*nit)->AtBoundary())) <<" ";
                  switch (prop_key.type)
                    {
                       case SCALAR:
                           cout << (*nit)->Read( prop_key );
                         break;
                       case VECTOR:
                           (*nit)->Read( prop_key, vc );
                           cout << vc;
                         break;
                       case TENSOR:
                           (*nit)->Read( prop_key, ts );
                           cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            (*nit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            (*nit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::OutputVariableToScreen:",
                                               "type of Node variable not recognized." );
                    }
                }
              break;
          case ELEMENT_INTEGRATION_POINT:
              for ( auto eit=CellsBegin(); eit!=CellsEnd(); eit++ )
                {
                  cout <<"\nParent cell Idx: " << (*eit)->Idx() <<"\t\t";
                  for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                    switch (prop_key.type)
                      {
                         case SCALAR:
                             cout << i <<": "<< (*eit)->Read( i, prop_key );
                           break;
                         case VECTOR:
                             (*eit)->Read( i, prop_key, vc ); cout << i <<":"<< vc;
                           break;
                         case TENSOR:
                             (*eit)->Read( i, prop_key, ts ); cout << i <<":"<< ts;
                           break;
                       case ARRAY: {
                            ArrayVariable ary;
                            (*eit)->Read( i, prop_key, ary );
                            ary.Out();
                         }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            (*eit)->Read( i, prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::OutputVariableToScreen:",
                                               "type of Element integration point variable not recognized." );
                      }
                }
              break;
           case ELEMENT:
              for ( auto eit=CellsBegin(); eit!=CellsEnd(); eit++ )
                {
                  cout <<"\nIdx: " << (*eit)->Idx() <<"\t\t";
                  /// Roman, 2014 (Face&InterFace): Should Face contain AtBoundary flag?
                  //cout << string(parseBoundary((*eit)->AtBoundary())) <<" ";
                  switch (prop_key.type)
                    {
                       case SCALAR: cout << (*eit)->Read( prop_key );
                         break;
                       case VECTOR: (*eit)->Read( prop_key, vc ); cout << vc;
                         break;
                       case TENSOR: (*eit)->Read( prop_key, ts ); cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            (*eit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            (*eit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::OutputVariableToScreen:",
                                               "type of Element variable not recognized." );
                    }
                }
              break;
          default:
               cerr << "\nModelSubDomain<"<< dim <<">::OutputVariableToScreen: Property ";
               cerr << prop <<"  "<< place_name <<" cannot be retrieved from Region/Boundary/SplitBoundary objects; use Read().";
               break;
       }
    cout << endl;
    cout.flush();

 } // end OutputVariableToScreen







/**
    prints state of the object manifest in:
    
    subdomain_name_ - passed down when region is created so that it can be referred to

    cell_vec_ - doubly sorted, interior cells first
    
    bd_face_vec_ - face-IDs of cells at the domain boundary as in second segment of cell_vec_

    node_vec_ - doubly sorted, interior nodes first
    
    first_bd_node_ - first node in the boundary range

    verbose_ - or not

*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::Out() const
 {
    cout <<"\nModelSubDomain<dim,CELL>::Out(): name: '"<< subdomain_name_ <<"', unique subdomain index: "<< DomainIndex();
    cout <<"\n\tmember cells("<< cell_vec_.size() <<"): interior="<< InteriorCells();
    cout <<", perimeter="<< cell_vec_.size()-InteriorCells();
    cout <<"\n\tmember nodes ("<< node_vec_.size() <<"): interior nodes="<< node_vec_.size() - PerimeterNodes();
    cout <<", perimeter nodes="<< PerimeterNodes();
    size_t perimeter_faces(0U);
    for ( size_t i{0U}; i<bd_face_vec_.size(); ++i ) perimeter_faces += bd_face_vec_[i].size();
    cout <<"\n\tperimeter faces="<< perimeter_faces;
    cout <<"\n\tdetailed listing of cells and nodes:";

    for ( const auto& it : cell_vec_ ) {
         if ( it == nullptr )
           throw csmp::Exception( ERROR, "ModelSubDomain<dim,CELL>::Out",
                                 "member cell pointer not initialised");
//         else (*it)->Out();
      }

    cout <<"\n\n"<<"perimeter cells and their perimeter faces (current local numbering): "<< endl;
    auto  bit{ bd_face_vec_.begin() };
    for ( size_t i=InteriorCells(); i<cell_vec_.size(); i++, bit++ ) {
         cout <<"\n\t\t"<<"cell "<< i <<" ("<< parseFiniteElementType( cell_vec_[i]->FE_Type() ) <<"): edge face numbers: ";
         for ( auto ft=(*bit).begin(); ft!=(*bit).end(); ft++ ) cout << (*ft) <<" ";
      }

    cout <<"\n\n"<<"perimeter nodes: "<< node_vec_.size() - first_bd_node_ <<" (current local numbering):"<< endl;
    for ( size_t i=first_bd_node_; i<node_vec_.size(); i++ ) {
         if ( node_vec_[i] == nullptr )
           throw csmp::Exception( ERROR, "ModelSubDomain<dim>::Out", "member node pointer not initialised.");
         else cout <<"\n\t\t"<<"node "<< node_vec_[i]->Idx() <<" ("<< parseBoundary( node_vec_[i]->AtBoundary() ) <<")";
      }

    cout << endl;
   
 } // Out




/**
    Writes data block containing all information that is needed to reconstruct a ModelSubDomain
    
    @note this is not a function because template template parameters are not allowed for for functions.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::WriteDomainIndexesToBinaryFile( fstream& fp ) const
 {
    // 1. writing name of the region
    binaryFileWrite( fp, Name().c_str() );
   
    // 2. writing the interior cell records of the region
    std::vector<uint32_t> IDs( distance(CellsBegin(), PerimeterCellsBegin() ) );
    transform( CellsBegin(), PerimeterCellsBegin(),
               IDs.begin(), []( const CELL<dim>* const ptr ){ return ptr->Idx(); } ); // tested: OK
    binaryFileWrite( fp, IDs );

    // 3. writing the perimeter cell records of the region
    IDs.resize( distance(PerimeterCellsBegin(), CellsEnd()) );
    transform( PerimeterCellsBegin(), CellsEnd(),
               IDs.begin(), []( const CELL<dim>* const ptr ){ return ptr->Idx(); } );
    binaryFileWrite( fp, IDs );
   
    // 4. writing the boundary faces
    // (this not stored because pointer locations will change in reconstruction eliminating storage benefit)
    // -----------------------------------------------------------------------------------------------------
    /* 
        since this is vector of vectors predominated by single value entries,
        it is collapsed into a flat vector in which all entries that refer to 
        multiple values per cell are prefaced by a negative numer that indicates
        how many multiple faces per cell follow, for example
        1 5  5 3 -2 6 2 3 3 5 6
                    ^^^          marking the 2 local face indices that relate to an cell that has
        2 faces on the model boundary.
    */

    // 5. writing the interior nodes
    IDs.resize( InteriorNodes() );
    transform( NodesBegin(), PerimeterNodesBegin(),
               IDs.begin(), []( const Node<dim>* const ptr ){ return ptr->Idx(); } );
    binaryFileWrite( fp, IDs );
   
    // 6. writing the perimeter nodes
    IDs.resize( PerimeterNodes() );
    transform( PerimeterNodesBegin(), NodesEnd(),
               IDs.begin(), []( const Node<dim>* const ptr ){ return ptr->Idx(); } );
    binaryFileWrite( fp, IDs );
    
    // NB: the connectivity between the cells is not stored because it is handled by MeshManager
   
 } // end WriteDomainIndexesToBinaryFile


/* TESTING (code snippet)
cerr <<"\n\n\nregion: "<< Name() <<"\n";
cerr <<"interior nodes";
out( IDs );
cerr <<"\n perimeter nodes";
out( IDs );
*/




/**
    Reads all the data required to fully reconstruct a ModelSubDomain (without search operations)
 
@note SKM: refactored 23/8/2018: no longer uses boundary face vector because the
    the sorting of the cell vectors during the model reconstruction invalidates
    this vector. It is therefore cheaper to rebuild the vector from scratch
    during the reconstruction.
*/
void readDomainIndexesFromBinaryFile( uint32_t dim, fstream& fp, SubDomainInfo& info )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    // 1. reading name of the subdomain
    char name[INFO_STRING];
    binaryFileRead( fp, name );
    info.name = name;
    assert( !info.name.empty() );
   
    // 2. reading the interior cell records of the region
    binaryFileRead( fp, info.interior_elmts );
    if (dim > 2 && info.interior_elmts.empty() ) {
        csmp_error.Note( WARNING, "readDomainIndexesFromBinaryFile:",
                          "Model appears to have a region with no interior cells: ", name );
    }

    // 3. reading the perimeter cell records of the region
    binaryFileRead( fp, info.perimeter_elmts );
    assert( !info.perimeter_elmts.empty() );
   
    // 4. reading the interior nodes
    binaryFileRead( fp, info.interior_nodes );
  
    // 5. reading the perimeter nodes
    binaryFileRead( fp, info.perimeter_nodes );
    assert( !info.perimeter_nodes.empty() );
   
 } // end readDomainIndexesFromBinaryFile




/**
       Writes Node information to a text file for reading with Paraview (TableToPointData) or spreadsheet tools.
       The interior nodes are distinguished from the perimeter nodes by the keywords "INTERIOR" and "PERIMETER."
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::NodeAttributesToCSV()
  {
     ofstream ofs( Name() + "_node_attribute.csv" );

     //ofs <<"\n"<< Name();
     ofs <<"\nx,y,z,node_id,bflag,rflag";
     for ( auto n=NodesBegin(); n!=PerimeterNodesBegin(); ++n ) {
          ofs <<"\n"<< (*n)->x() <<","<< (*n)->y() <<","<< (*n)->z() <<","<< (*n)->Idx() <<",";
          ofs << parseBoundary((*n)->AtBoundary()) <<","<<"INTERIOR";
       }
     for ( auto n=PerimeterNodesBegin(); n!=NodesEnd(); ++n ) {
          ofs <<"\n"<< (*n)->x() <<","<< (*n)->y() <<","<< (*n)->z() <<","<< (*n)->Idx() <<",";
          ofs << parseBoundary((*n)->AtBoundary()) <<","<<"PERIMETER";
       }
     ofs << endl;

 } // end NodeAttributesToCSV



/**
    Removes cells and nodes and rebuilds the    bd_face_vec_  and   node_vec_ if necessary.
    @return size_t  the number of cells removed.
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::RemoveNullPointerCells()
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t n_cells{ cell_vec_.size() };
    
    // determining how the new number of interior cells
    const size_t n_interior_cells_new = InteriorCells() - count( CellsBegin(), PerimeterCellsBegin(), nullptr );
    
    // erasing cell vector without changing the relative number of its cells
    cell_vec_.erase( remove( cell_vec_.begin(), cell_vec_.end(), nullptr ), cell_vec_.end() );
    const size_t n_cells_new = cell_vec_.size();
    
    if ( n_cells_new < n_cells ) BuildPerimeterFaceVector( n_interior_cells_new );
    
    // removing node pointers if any
    size_t nodes_removed = node_vec_.size();
    node_vec_.erase( remove( node_vec_.begin(), node_vec_.end(), nullptr ), node_vec_.end() );
    nodes_removed -= node_vec_.size();
    
    if ( n_cells_new < n_cells == 0 && nodes_removed > 0 )
      csmp_error.Note( ERROR, "ModelSubDomain<dim,CELL>::RemoveNullPointerCells",
                        "removed nodes but not cells? - subdomain may be corrupt now.");
    
    return n_cells - n_cells_new;
 }






/**
    Deletes nullptr cells, rebuilds node vector, sorts everything and re-establishes the perimeter face vectors after modifications of cells.
    
    @attention method also deals with the case where the new cells have been added to the subdomain.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::RebuildSubDomainAfterChangeOfCellVector()
 {
    // erasing cell vector without changing the relative number of its cells
    cell_vec_.erase( remove( cell_vec_.begin(), cell_vec_.end(), nullptr ), cell_vec_.end() );
    cell_vec_.shrink_to_fit();
    
    // rebuild the node vector
    CreateNodePointerVector();
    // includes shrink to fit
    
    // sorting vectors and identifying perimeter cells and nodes
    IdentifyPerimeter();
    // the following happens inside of IdentifyPerimeter()->PartitionVectors()
    // BuildPerimeterFaceVector( InteriorCells() );
    
    rebuilt_needed_ = false;
    
 } // end RebuildSubDomainAfterChangeOfCellVector




    /// rebuilds subdomain on the basis of the cells that will be selected according to the supplied property constraints
template<uint32_t dim, template<uint32_t> class CELL>
void ModelSubDomain<dim,CELL>::UpdateCellMembershipApplyingConstraints( typename vector<CELL<dim>*>::const_iterator start,
                                                                        typename vector<CELL<dim>*>::const_iterator end,
                                                                        const PropertyConstraints& constraints )
 {
    cell_vec_.clear();
    
    // selecting the cells on the basis of the criteria specified in PropertyContraints
    while( start != end ) {
         if ( constraints.CheckConstraints( (*start) ) )
           cell_vec_.push_back( (*start) );
         ++start;
      }
    
    // rebuild the node vector
    CreateNodePointerVector();
    // includes shrink to fit
    
    // sorting vectors and identifying perimeter cells and nodes
    IdentifyPerimeter();
    // the following happens inside of IdentifyPerimeter()->PartitionVectors()
    // BuildPerimeterFaceVector( InteriorCells() );
    
    rebuilt_needed_ = false;

 } // end UpdateCellMembershipApplyingConstraints





/**
     @param start iterator to beginning of range that will be checked for region membership
     @param end iterator behind last element in the range
     @return returns the number of nodes in the supplied range, which also form part of the perimeter of this model subdomain
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t ModelSubDomain<dim,CELL>::SharedPerimeterNodes( typename vector<csmp::Node<dim>*>::const_iterator start,
                                                       typename vector<csmp::Node<dim>*>::const_iterator end ) const
 {
    if ( start == end ) return 0U;
 
    auto   first1(PerimeterNodesBegin());
    size_t shared_nodes(0U);

    // comparing the boundary nodes
    while ( first1 != NodesEnd() and start != end )
      {
        if ( *first1 < *start ) ++first1;
        else if ( *start < *first1 ) ++start;
        else {
             shared_nodes++;
             first1++;
             start++;
          }
      }

    return shared_nodes;
 }


// NON-MEMBER FUNCTIONS

/// distinguishes between Region, Boundary and SplitBoundary on the basis of the name string
PLACEMENT modelSubdomainType( const std::string& subdomain_name )
  {
     // 1. method isDiag.. needs upper case
     if ( isDiagnosticBoxBoundaryClassifier(subdomain_name) ) return BOUNDARY;
     
     // making the search case insensitive search
     string lowercase_name;
     for ( auto& it : subdomain_name )
       lowercase_name += tolower( it );

     // 2. split boundary first because boundary is a substring of splitboundary
     // (before boundaries to avoid ambiguity due to shared string 'boundary'
     if ( lowercase_name.find("splitboundary") != string::npos ||
          lowercase_name.find("split_boundary") != string::npos ||
          lowercase_name.find("split boundary") != string::npos ||
          lowercase_name.find("split-boundary") != string::npos )
       return SPLIT_BOUNDARY;

     // 3. boundaries as identified by name string
     if ( lowercase_name.find("boundary") != string::npos )
       return BOUNDARY;
     
     return REGION;

  } // end modelSubdomainType




/**
     Counts and returns the number of nodes shared between the two regions.
     Special attention is paid to the fact the nodes are sorted in two separate ranges for the interior and exterior.
     All potential combinations are considered.
     
     @attention function assumes that both subdomains are valid, containing no multiple nodes.
 */
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedNodes( const ModelSubDomain<dim,CELL>& g1, const ModelSubDomain<dim,CELL>& g2 )
 {
    // creating some copies of the sorted node vectors
    vector<Node<dim>*>  g1_nodes( g1.NodeVector() );
    vector<Node<dim>*>  g2_nodes( g2.NodeVector() );
    // not much extra sorting is required because the separate ranges were already sorted
    sort( g1_nodes.begin(), g1_nodes.end() );
    sort( g2_nodes.begin(), g2_nodes.end() );
    
    vector<Node<dim>*> shared_nodes;
    set_intersection( g1_nodes.begin(), g1_nodes.end(), g2_nodes.begin(), g2_nodes.end(),
                      back_inserter(shared_nodes) );

    return shared_nodes.size();

 } // end sharedNodes

template size_t sharedNodes( const ModelSubDomain<1U,Element>&, const ModelSubDomain<1U,Element>& );
template size_t sharedNodes( const ModelSubDomain<2U,Element>&, const ModelSubDomain<2U,Element>& );
template size_t sharedNodes( const ModelSubDomain<3U,Element>&, const ModelSubDomain<3U,Element>& );

template size_t sharedNodes( const ModelSubDomain<1U,Face>&, const ModelSubDomain<1U,Face>& );
template size_t sharedNodes( const ModelSubDomain<2U,Face>&, const ModelSubDomain<2U,Face>& );
template size_t sharedNodes( const ModelSubDomain<3U,Face>&, const ModelSubDomain<3U,Face>& );




/** 
    As for sharedNodes() - but application is restricted to nodes that sit
    on the perimeter / surface of the region.
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>& g1, const ModelSubDomain<dim,CELL>& g2 )
 {
    // creating some copies of the sorted node vectors
    vector<Node<dim>*>  g1_nodes( g1.PerimeterNodesBegin(), g1.NodesEnd() );
    vector<Node<dim>*>  g2_nodes( g2.PerimeterNodesBegin(), g2.NodesEnd() );
    
    // not sorting is required because the node ranges are already sorted
    vector<Node<dim>*> shared_nodes;
    set_intersection( g1_nodes.begin(), g1_nodes.end(), g2_nodes.begin(), g2_nodes.end(),
                      back_inserter(shared_nodes) );

#ifdef MODEL_SUBDOMAIN_DEBUG
if ( shared_nodes.size() == 2U ) {
    cout <<"\n"<<"Nodes shared between '"<< g1.Name() <<"' and '"<< g2.Name() <<"': ";
    cout << shared_nodes[0]->Idx() <<": "<< shared_nodes[0]->Coordinate() <<",  ";
    cout << shared_nodes[1]->Idx() <<": "<< shared_nodes[1]->Coordinate();
  }
#endif
    return shared_nodes.size();

 } // end sharedPerimeterNodes

template size_t sharedPerimeterNodes( const ModelSubDomain<1U,Element>&, const ModelSubDomain<1U,Element>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<2U,Element>&, const ModelSubDomain<2U,Element>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<3U,Element>&, const ModelSubDomain<3U,Element>& );

template size_t sharedPerimeterNodes( const ModelSubDomain<1U,Face>&, const ModelSubDomain<1U,Face>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<2U,Face>&, const ModelSubDomain<2U,Face>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<3U,Face>&, const ModelSubDomain<3U,Face>& );

template size_t sharedPerimeterNodes( const ModelSubDomain<1U,InterFace>&, const ModelSubDomain<1U,InterFace>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<2U,InterFace>&, const ModelSubDomain<2U,InterFace>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<3U,InterFace>&, const ModelSubDomain<3U,InterFace>& );




/**
    As method above, but returning pointers to the Nodes found into the argument vector.
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>& g1, const ModelSubDomain<dim,CELL>& g2, vector<Node<dim>*>& shared_nodes )
 {
    if ( !shared_nodes.empty() ) shared_nodes.clear();
 
    // creating some copies of the sorted node vectors
    vector<Node<dim>*>  g1_nodes( g1.PerimeterNodesBegin(), g1.NodesEnd() );
    vector<Node<dim>*>  g2_nodes( g2.PerimeterNodesBegin(), g2.NodesEnd() );
    
    // not sorting is required because the node ranges are already sorted
    set_intersection( g1_nodes.begin(), g1_nodes.end(), g2_nodes.begin(), g2_nodes.end(),
                      back_inserter(shared_nodes) );

    return shared_nodes.size();

 } // end sharedPerimeterNodes

template size_t sharedPerimeterNodes( const ModelSubDomain<1U,Element>&, const ModelSubDomain<1U,Element>&, vector<Node<1U>*>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<2U,Element>&, const ModelSubDomain<2U,Element>&, vector<Node<2U>*>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<3U,Element>&, const ModelSubDomain<3U,Element>&, vector<Node<3U>*>& );

template size_t sharedPerimeterNodes( const ModelSubDomain<1U,Face>&, const ModelSubDomain<1U,Face>&, vector<Node<1U>*>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<2U,Face>&, const ModelSubDomain<2U,Face>&, vector<Node<2U>*>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<3U,Face>&, const ModelSubDomain<3U,Face>&, vector<Node<3U>*>& );

template size_t sharedPerimeterNodes( const ModelSubDomain<1U,InterFace>&, const ModelSubDomain<1U,InterFace>&, vector<Node<1U>*>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<2U,InterFace>&, const ModelSubDomain<2U,InterFace>&, vector<Node<2U>*>& );
template size_t sharedPerimeterNodes( const ModelSubDomain<3U,InterFace>&, const ModelSubDomain<3U,InterFace>&, vector<Node<3U>*>& );








/**
       Finds the higher-dimensional perimeter cells that share their perimeter faces on the touching two subdomains.
       
       @param subdomain1 any unique (non-overlapping) Region, Boundary or Splitboundary
       @param subdomain2 subdomain that is expected to share outside faces with the first subdomain
       @param matching_cells  map of higher-dimensional cells that share a face at the contact between the two subdomains
       @return the number of shared perimeter faces that were found.
       
       @todo check for interpenetrating regions
       
       @author SKM
       @date 16.4.22
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedPerimeterCells( const ModelSubDomain<dim,CELL>& subdomain1, const ModelSubDomain<dim,CELL>& subdomain2,
                              vector<pair<pair<CELL<dim>*,uint32_t>,pair<CELL<dim>*,uint32_t> > >& matching_cells )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 0. some initial checks
    // ----------------------
    if ( subdomain1.Name() == subdomain2.Name() ) {
         csmp_error.Note( ERROR, "sharedPerimeterCells:", "argument regions are the same; nothing was done.");
         return 0U;
      }
    if ( subdomain1.Empty() ) {
         csmp_error.Note( ERROR, "sharedPerimeterCells:", subdomain1.Name(), "is empty; nothing was done.");
         return 0U;
      }
    if ( subdomain2.Empty() ) {
         csmp_error.Note( ERROR, "sharedPerimeterCells:", subdomain2.Name(), "is empty; nothing was done.");
         return 0U;
      }
    if ( subdomain2.Empty() ) {
         csmp_error.Note( ERROR, "sharedPerimeterCells:", subdomain2.Name(), "is empty; nothing was done.");
         return 0U;
      }
    if ( subdomain1.NeedsRebuilt() || subdomain2.NeedsRebuilt() ) {
         csmp_error.Note( ERROR, "sharedPerimeterCells:", "input subdomains have been flagged for rebuilt; nothing was done.");
         return 0U;
      }
    const auto cell_dim1 = subdomain1.SpatialDimensions();
    const auto cell_dim2 = subdomain2.SpatialDimensions();
    if ( cell_dim1.second != cell_dim2.second ) {
         cerr <<"\n\t'"<< subdomain1.Name() <<"': spatial dimension: "<< cell_dim1.second <<" vs. '";
         cerr << subdomain2.Name() <<"': spatial dimension: "<< cell_dim2.second;
         csmp_error.Note( WARNING, "sharedPerimeterCells:", "for successful processing input domains must have the same spatial dimension.");
         return 0U;
      }
    if ( cell_dim1.first > 1 ||  cell_dim2.first > 1 ) {
         cerr <<"\n\t'"<< subdomain1.Name() <<"': number of spatial dimensions: "<< cell_dim1.first <<" vs. '";
         cerr << subdomain2.Name() <<"': number of spatial dimensions: "<< cell_dim2.first;
         csmp_error.Note( WARNING, "sharedPerimeterCells:", "function can only process domains with elements of a single spatial dimension.");
         return 0U;
      }


    // 1. is there a shared interface? - looping over the perimeter faces of the adjacent regions
    // ------------------------------------------------------------------------------------------
    map<set<Node<dim>*>,pair<pair<CELL<dim>*,uint32_t>,pair<CELL<dim>*,uint32_t> > >  shared_perimeter_cells;
    // (determining contacts via shared corner-node-ptrs of cells adjacent to interface)
    // (not only these elements but also their face numbers will be recorded (inside elements will be first in pair)

    // starting with region1 - assuming that no face-matching can occur at this stage
    for ( size_t i{ subdomain1.InteriorCells() }; i<subdomain1.Cells(); i++ )
      for ( auto j{0U}; j<subdomain1.PerimeterFaces(i); j++ ) {
           auto pface = subdomain1.PerimeterFace(i,j);
           // making a search key from the corner nodes of the face and recording the perimeter element and its face number
           assert( !subdomain1.E(i)->IsLine() );
           shared_perimeter_cells.insert( make_pair( subdomain1.E(i)->CornerNodesOfFace(pface),
                                          make_pair( make_pair( subdomain1.E(i), pface ), make_pair( nullptr,NULL_IDX) ) ) );
        }
    // for the second region2, do the same, but when matching faces are found corresponding elements and face ids are assigned to second element-face pair
    size_t n_matching_faces{0U};
    for ( size_t i{ subdomain2.InteriorCells() }; i<subdomain2.Cells(); i++ )
      for ( auto j{0U}; j<subdomain2.PerimeterFaces(i); j++ ) {
           auto pface = subdomain2.PerimeterFace(i,j);
           // making a search key from the corner nodes of the face and recording the perimeter element and its face number
           assert( !subdomain2.E(i)->IsLine() );
           auto it = shared_perimeter_cells.insert( make_pair( subdomain2.E(i)->CornerNodesOfFace(pface),
                                                    make_pair( make_pair( subdomain2.E(i), pface ), make_pair( nullptr,NULL_IDX) ) ) );
           // if a matching face is found
           if ( it.second == false ) { // no new insertion could be made into map with unique keys
                (*it.first).second.second = make_pair( subdomain2.E(i), pface );
                n_matching_faces++;
             }
        }

    if ( n_matching_faces == 0U ) {
         string message( string(" input regions '") + subdomain1.Name() + "' and '" + subdomain2.Name() +"'");
         csmp_error.Note( WARNING, "sharedPerimeterCells:",
                            message, "do not share any faces, nothing could be done");
         return 0U;
      }

    // 2. Eliminating the single element entries from the 'shared_perimeter_faces' map
    // -------------------------------------------------------------------------------
    for ( auto it = shared_perimeter_cells.begin(); it != shared_perimeter_cells.end() /* not hoisted */; /* no increment */ ) {
         // there is only a single element in the Element-pointer pair, the map entry will be deleted
         if ( (*it).second.second.first == nullptr ) it = shared_perimeter_cells.erase(it); // since C++11
         else ++it;
      }

    // 3. collecting the required cell pairs from the 'shared_perimeter_faces' map
    // ---------------------------------------------------------------------------
    //vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > > matching_cells;
    matching_cells.reserve( shared_perimeter_cells.size() );
    for ( const auto& pit : shared_perimeter_cells )
      if ( pit.second.second.first != nullptr )
        matching_cells.emplace_back( pit.second );

    return matching_cells.size();
  
 } // end sharedPerimeterCells


template size_t sharedPerimeterCells( const ModelSubDomain<1U,Element>&, const ModelSubDomain<1U,Element>&, vector<pair<pair<Element<1>*,uint32_t>,pair<Element<1>*,uint32_t> > >& );
template size_t sharedPerimeterCells( const ModelSubDomain<2U,Element>&, const ModelSubDomain<2U,Element>&, vector<pair<pair<Element<2>*,uint32_t>,pair<Element<2>*,uint32_t> > >& );
template size_t sharedPerimeterCells( const ModelSubDomain<3U,Element>&, const ModelSubDomain<3U,Element>&, vector<pair<pair<Element<3>*,uint32_t>,pair<Element<3>*,uint32_t> > >& );

template size_t sharedPerimeterCells( const ModelSubDomain<1U,Face>&, const ModelSubDomain<1U,Face>&, vector<pair<pair<Face<1>*,uint32_t>,pair<Face<1>*,uint32_t> > >& );
template size_t sharedPerimeterCells( const ModelSubDomain<2U,Face>&, const ModelSubDomain<2U,Face>&, vector<pair<pair<Face<2>*,uint32_t>,pair<Face<2>*,uint32_t> > >& );
template size_t sharedPerimeterCells( const ModelSubDomain<3U,Face>&, const ModelSubDomain<3U,Face>&, vector<pair<pair<Face<3>*,uint32_t>,pair<Face<3>*,uint32_t> > >& );

/* MORE COMPLICATED BECAUSE face nodes might not be shared
template size_t sharedPerimeterCells( const ModelSubDomain<1U,InterFace>&, const ModelSubDomain<1U,InterFace>&, vector<pair<pair<Element<1>*,uint32_t>,pair<Element<1>*,uint32_t> > >& );
template size_t sharedPerimeterCells( const ModelSubDomain<2U,InterFace>&, const ModelSubDomain<2U,InterFace>&, vector<pair<pair<Element<2>*,uint32_t>,pair<Element<2>*,uint32_t> > >& );
template size_t sharedPerimeterCells( const ModelSubDomain<3U,InterFace>&, const ModelSubDomain<3U,InterFace>&, vector<pair<pair<Element<3>*,uint32_t>,pair<Element<3>*,uint32_t> > >& );
*/



 template class ModelSubDomain<1U,Element>;
 template class ModelSubDomain<2U,Element>;
 template class ModelSubDomain<3U,Element>;

 template class ModelSubDomain<1U,Face>;
 template class ModelSubDomain<2U,Face>;
 template class ModelSubDomain<3U,Face>;

 template class ModelSubDomain<1U,InterFace>;
 template class ModelSubDomain<2U,InterFace>;
 template class ModelSubDomain<3U,InterFace>;

} // end namespace

