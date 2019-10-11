#include "Boundary.h"
#include "Box.h"
#include "Region.h"

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "FEM_Data.h"

#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"

#include "Exception.h"
#include "ErrorHandler.h"

#include "binaryReadWrite.h"
#include "CSMP_highLevelUtilities.h"
#include "PL_Utilities.h"
#include "variableOperations.h"

#include "Visitor.h"
#include "FaceConstructionData.h"

//#define BOUNDARY_DEBUG

using namespace std;

namespace csmp {

/**
Constructor for an empty boundary (no faces yet).
*/
template<size_t dim>
Boundary<dim>::Boundary( std::string boundaryname, const PropertyDatabase<dim>& pref, BOX_BOUNDARY flag )
  : ModelSubDomain<dim, Face>( boundaryname, pref ),
    boundaryFlag_( flag )
{
  if ( dim == 1 )
    ErrorHandler::Instance().notice( ERROR, "Boundary<dim>::(custom constructor):", "boundary objects are only supported in 2 & 3D models." );
  this->ResizePropertyStorage( this->pref_.LocalVariablesAt( Placement() ) );
}


/// copy constructor
template<size_t dim>
Boundary<dim>::Boundary( const Boundary& ed )
  : ModelSubDomain<dim, Face>( ed ),
    boundaryFlag_( ed.boundaryFlag_ )
{
}


/// move constructor
template<size_t dim>
Boundary<dim>::Boundary( Boundary&& ed )
  : ModelSubDomain<dim, Face>( ed ),
    boundaryFlag_( ed.boundaryFlag_ )
{
}



/**
(re)-constructor of boundary from the mesh manager in the model

@attention expects that the face numbering is unique and that all faces
are numbered consecutivly.

@author SKM 4/5/2016

@test nodes are assigned correctly, Face objects are not
*/
template<size_t dim>
Boundary<dim>::Boundary( const PropertyDatabase<dim>& pref,
                         MeshManager<dim>& mesh,
                         const SubDomainInfo& info,
                         BOX_BOUNDARY bflag )
  : ModelSubDomain<dim, Face>( info.name, pref ),
    boundaryFlag_( bflag )
{
  // building the face vector
  // ------------------------
  const size_t elements( mesh.Elements() );
  this->elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( &mesh, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // assigning pointers to the interior faces
  for ( size_t i : info.interior_elmts )
    this->elmt_vec_.push_back( faces[i - elements] );

  // assigning pointers to the perimeter faces
  for ( size_t i : info.perimeter_elmts )
    this->elmt_vec_.push_back( faces[i - elements] );

  // building the node vector
  // ------------------------
  // assigning pointers to the interior and perimeter nodes
  this->first_bd_node_ = info.interior_nodes.size();
  this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );

  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // assigning pointers to the interior faces
  for ( size_t i : info.interior_nodes )
    this->node_vec_.push_back( nodes[i] );

  // assigning pointers to the perimeter faces
  for ( size_t i : info.perimeter_nodes )
    this->node_vec_.push_back( nodes[i] );

  this->SortVectors( info.interior_elmts.size(), info.interior_nodes.size() );

  // building the vector of vectors of those faces (edges) of the faces that lie on the subdomain perimeter
  this->BuildBoundaryFaceVector( info.interior_elmts.size() );

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( BOUNDARY ) );

} // end re-constructor







/**
    new constructor: re-constructor of boundary from index data stored in SubDomainInfo and faces from the MeshManager
*/
template<size_t dim>
Boundary<dim>::Boundary( const PropertyDatabase<dim>& pref,
                         const size_t& elements,
                         const std::deque<Node<dim>*>& nodes,
                         const std::deque<Face<dim>*>& faces,
                         const SubDomainInfo& info,
                         BOX_BOUNDARY bflag )
  : ModelSubDomain<dim, Face>( info.name, pref ),
    boundaryFlag_( bflag )
{
  // building the face vector
  // ------------------------	 
  this->elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );

  // assigning pointers to the interior faces
  for ( size_t i : info.interior_elmts )
    if( (i - elements) < faces.size() )
      this->elmt_vec_.push_back( faces[i - elements] );

  // assigning pointers to the perimeter faces
  for ( size_t i : info.perimeter_elmts )
    if ( (i - elements) < faces.size() )
      this->elmt_vec_.push_back( faces[i - elements] );

  // building the node vector
  // ------------------------
  // assigning pointers to the interior and perimeter nodes
  this->first_bd_node_ = info.interior_nodes.size();
  this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );

  // assigning pointers to the interior faces
  for ( size_t i : info.interior_nodes )
    this->node_vec_.push_back( nodes[i] );

  // assigning pointers to the perimeter faces
  for ( size_t i : info.perimeter_nodes )
    this->node_vec_.push_back( nodes[i] );

  this->SortVectors( info.interior_elmts.size(), info.interior_nodes.size() );

  // building the vector of vectors of those faces (edges) of the faces that lie on the subdomain perimeter
  this->BuildBoundaryFaceVector( info.interior_elmts.size() );

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( BOUNDARY ) );
}

template<size_t dim>
IntegrationPointVariables Boundary<dim>::FaceIntegrationPointVariables() const
{ return this->pref_.IntegrationPointVariablesAt( FACE ); }

template<size_t dim>
LocalVariables Boundary<dim>::FaceVariables() const
{ return this->pref_.LocalVariablesAt( FACE ); }



/**
Creates boundary from supplied vector of faces that must already
know their higher-dimensional parent elements and their equidimensional
neighbor faces.

@attention the neighborhood relations are not re-established.

@note this constructor is mainly used for edges of line-element faces.
*/
template<size_t dim>
Boundary<dim>::Boundary( const string& boundary_name,
                         const PropertyDatabase<dim>& dbase_ref,
                         typename std::vector<Face<dim>*>::iterator facesBegin,
                         typename std::vector<Face<dim>*>::iterator facesEnd,
                         BOX_BOUNDARY flag )
  : ModelSubDomain<dim, Face>( boundary_name, dbase_ref ),
    boundaryFlag_( flag )
{
  // moving the supplied Face pointers into the element storage
  this->elmt_vec_.assign( facesBegin, facesEnd );
  // initialising node pointer vector, sorting nodes and elements, and creating boundary face vector
  const bool updateFaceConnectivity( false ); // has been done before
  const bool updateIndexes( false );          // not necessarily needed
  Initialize( flag, updateFaceConnectivity, updateIndexes );
}


/**
Default dtor does not delete faces, even though constructed by Boundary. Need to call DeleteFaces explicitly!
@todo FIX THIS CONCEPTUAL MISUNDERSTANDING ALL THE WAY THROUGH:  the faces are managed by the MeshManager
*/
template<size_t dim>
Boundary<dim>::~Boundary()
{
}

template<size_t dim>
Boundary<dim>& Boundary<dim>::operator=( const Boundary<dim>& ed )
{
  if ( &ed != this )
  {
    *this = ed;
    boundaryFlag_ = ed.boundaryFlag_;
  }
  return *this;
}

// LOCAL VARIABLE STORAGE INTERFACE

// TODO: replace with non-member function
template<size_t dim>
bool Boundary<dim>::ValidVariable( const char* variableName ) const
{
  const PLACEMENT p( this->pref_.Placement( variableName ) );
  if ( p == NODE || p == FACE || p == BOUNDARY )
    return true;
  return false;
}


// VISITORS INTERFACE

/**
Visitation of boundary, like for any visitation,
the application level determines where the visitation starts
and where it ends.
*/
template<size_t dim>
void Boundary<dim>::Accept( Visitor<dim>& v )
{
  if ( v.ApplicationLevel() == MODEL || v.ApplicationLevel() == BOUNDARY )
    v.Visit( this );

  switch ( v.ApplicationTarget() ) {
    case MODEL:
      throw csmp::Exception( ERROR, "Region<dim>::Accept",
                             "ApplicationTarget MODEL; Visitor should have never arrived at this boundary" );
      break;
    case BOUNDARY:
      return;
    case FACE:
      for ( typename vector<Face<dim>*>::iterator
            it = this->ElementsBegin(); it != this->ElementsEnd(); it++ )
        (*it)->Accept( v );
      return;
    case NODE:
      for ( typename vector<csmp::Node<dim>*>::iterator
            nd_it = this->NodesBegin(); nd_it != this->NodesEnd(); nd_it++ )
        (*nd_it)->Accept( v );
      return;
    default:
      throw csmp::Exception( ERROR, "Boundary<dim>::Accept",
                             "ApplicationTarget was not resolved; nothing was done" );
  }

} // end Accept



  // -----------------------------------------------
  // Binary input/output
  // -----------------------------------------------




/**
Outs FEM_Data of provided variable type and placement to binary

@author  P. Lang
@date  9/30/2012

@param [in,out]  fp  If non-null, the binary file pointer
@param place         The placement of variables
@param vtype         The vtype of variables

@return  true if it succeeds, false if it fails.
*/
template<size_t dim>
template<class Var>
bool Boundary<dim>::Out( fstream& fp, PLACEMENT place, VARIABLE_TYPE vtype ) const
{
  set<string> propList;
  size_t vCount( 0 ), bytes( sizeof( size_t ) );

  this->pref_.ListProperties( place, vtype, propList );
  vCount = propList.size();
  fp.write( (char*)&vCount, bytes );
  if ( vCount != 0 )
  {
    FEM_Data<Var> femData;
    for ( set<string>::const_iterator it( propList.begin() ); it != propList.end(); ++it )
    {
      if ( !skm_C_fwrite( fp, it->c_str() ) )
        return false;
      OutputVariableTo( it->c_str(), femData );
      femData.OutBinary( fp );
    }
  }
  return true;
}





/**
Outputs specific property data to a FEM_Data container.
This function is a nested template:
the outer template provides double64 = data type and dim = dimension,
and Var the data type of the property
(ScalarVariable, VectorVariable, or TensorVariable).

@param data The FEM_Data container with the data that corresponds to the given property
is returned into the second method argument.

@section implementation Implementation

Only data pertinent to this local group is output.
A new data container containing only group relevant data is returned.

Thus, the data extracted only corresponds to -this- group and is organized
by local id numbers.

@section application Application

The idea is to extract the property data that corresponds to this group only.
*/
template<size_t dim>
template<class Var>
void Boundary<dim>::OutputVariableTo( const char* property, FEM_Data<Var>& data ) const
{
  csmp::Index  idx = this->pref_.StorageKey( property );
  Var          var;
  femDataOutputDispatch::initVariable( idx, var );

  switch ( idx.place ) {
    case FACE: {
      data.Reset( idx, this->elmt_vec_.size(), var );
      for ( typename vector<csmp::Face<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ ) {
        (*eit)->Read( idx, var );
        data[(*eit)->Idx()] = var;
      }
    }
               break;
    case NODE: {
      data.Reset( idx, this->node_vec_.size(), var );
      for ( typename vector<csmp::Node<dim>*>::const_iterator
            nit = this->node_vec_.begin(); nit != this->node_vec_.end(); nit++ ) {
        (*nit)->Read( idx, var );
        data[(*nit)->Idx()] = var;
      }
    }
               break;
    case ELEMENT_INTEGRATION_POINT: {
      data.Reset( idx, this->IntegrationPoints(), var );
      size_t counter( 0U );
      for ( typename vector<csmp::Face<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
      {
        for ( size_t i = 0U; i<(*eit)->IntegrationPoints(); i++ ) {
          (*eit)->Read( i, idx, var );
          data[counter + i] = var;
        }
        counter += (*eit)->IntegrationPoints();
      }
    }
                                    break;
    case BOUNDARY: {
      data.Reset( idx, 1U, var );
      this->Read( idx, var );
      data[0U] = var;
    }
                   break;
    default:break; {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      csmp_error.notice( WARNING, "Boundary::OutputVariableTo",
                         property, "placement not identified" );
    }
  }

} // end OutputVariableTo




template void Boundary<1>::OutputVariableTo<ScalarVariable >( const char*, FEM_Data<ScalarVariable >& ) const;
template void Boundary<2>::OutputVariableTo<ScalarVariable >( const char*, FEM_Data<ScalarVariable >& ) const;
template void Boundary<3>::OutputVariableTo<ScalarVariable >( const char*, FEM_Data<ScalarVariable >& ) const;

template void Boundary<1>::OutputVariableTo<ArrayVariable >( const char*, FEM_Data<ArrayVariable >& ) const;
template void Boundary<2>::OutputVariableTo<ArrayVariable >( const char*, FEM_Data<ArrayVariable >& ) const;
template void Boundary<3>::OutputVariableTo<ArrayVariable >( const char*, FEM_Data<ArrayVariable >& ) const;

template void Boundary<1>::OutputVariableTo<FlaggedArrayVariable >( const char*, FEM_Data<FlaggedArrayVariable >& ) const;
template void Boundary<2>::OutputVariableTo<FlaggedArrayVariable >( const char*, FEM_Data<FlaggedArrayVariable >& ) const;
template void Boundary<3>::OutputVariableTo<FlaggedArrayVariable >( const char*, FEM_Data<FlaggedArrayVariable >& ) const;

template void Boundary<1>::OutputVariableTo<VectorVariable<1U> >( const char*, FEM_Data<VectorVariable<1U> >& ) const;
template void Boundary<2>::OutputVariableTo<VectorVariable<2U> >( const char*, FEM_Data<VectorVariable<2U> >& ) const;
template void Boundary<3>::OutputVariableTo<VectorVariable<3U> >( const char*, FEM_Data<VectorVariable<3U> >& ) const;

template void Boundary<1>::OutputVariableTo<TensorVariable<1U> >( const char*, FEM_Data<TensorVariable<1U> >& ) const;
template void Boundary<2>::OutputVariableTo<TensorVariable<2U> >( const char*, FEM_Data<TensorVariable<2U> >& ) const;
template void Boundary<3>::OutputVariableTo<TensorVariable<3U> >( const char*, FEM_Data<TensorVariable<3U> >& ) const;




/**
Ins FEM_Data of provided variable type and placement from binary

Complimentary to Out(...) function.

@author  P. Lang
@date  9/30/2012

@param [in,out]  fp  If non-null, the binary file pointer

@return  true if it succeeds, false if it fails.
*/
template<size_t dim>
template<class Var>
bool Boundary<dim>::In( fstream& fp, PLACEMENT, VARIABLE_TYPE )
{
  size_t vCount( -1 );
  fp.read( (char*)&vCount, sizeof( size_t ) );
  if ( vCount == 0 )
    return true;
  FEM_Data<Var> femData;
  for ( size_t i( 0 ); i < vCount; ++i )
  {
    char propertyName[200];
    if ( !skm_C_fread( fp, propertyName ) )
      return false;
    femData.InBinary( fp );
    InputVariableFrom( propertyName, femData );
  }
  return true;
}


/**
Inputs geometry and property data from an FEM_Data object into this group.

The expectation is that nodes and element numbers in the VData set
correspond to these in the current Region and that elements / nodes
in the model have the same ordering as in the input dataset!

@section arguments Input Arguments

The name of the property which shall be input to VSet.

@section implementation Implementation

The current element and node indices are used to retrieve the data
from from the data container.

@section application Application

Enable to store data to disk to faciliate persistance.
*/
template<size_t dim>
template<class Var>
void Boundary<dim>::InputVariableFrom( const char* property,
                                       const FEM_Data<Var>& vdata )
{
  const csmp::Index  idx = this->pref_.StorageKey( property );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  switch ( idx.place ) {
    case FACE:
      //assert( this->Elements() == vdata.Size() );
      for ( typename vector<csmp::Face<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
        (*eit)->Store( idx, vdata[(*eit)->Idx()] );
      break;
    case NODE:
      //assert( this->Nodes() == vdata.Size() );
      for ( typename vector<csmp::Node<dim>*>::const_iterator
            nit = this->node_vec_.begin(); nit != this->node_vec_.end(); nit++ ) {
        (*nit)->Store( idx, vdata[(*nit)->Idx()] );
      }
      break;
    case ELEMENT_INTEGRATION_POINT: {
      // this will work only if the element ordering in the Region
      // from which the constraint point data were written out
      // with OutputDataTo() is exactly the same as in this Region.
      size_t  counter( 0U );
      //assert( this->IntegrationPoints() == vdata.Size() );
      for ( typename vector<csmp::Face<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
      {
        for ( size_t i = 0U; i<(*eit)->IntegrationPoints(); i++ )
          (*eit)->Store( i, idx, vdata[counter + i] );
        counter += (*eit)->IntegrationPoints();
      }
    }
    case BOUNDARY:
      if ( vdata.Size() != 1U )
        csmp_error.notice( WARNING, "Boundary<dim>::InputVariableFrom:",
                           "don't know which location in FEM_Data I should write the boundary variable to, using [0]." );

      this->Store( idx, vdata[0U] );
      break;
    default:
      break;/*
            throw csmp::Exception( ERROR, "Boundary<dim>::InputVariableFrom:",
            property, "Property placement could not be identified" );
            */
  }

} // end InputVariableFrom


template void Boundary<1>::InputVariableFrom<ScalarVariable>( const char*, const FEM_Data<ScalarVariable>& );
template void Boundary<2>::InputVariableFrom<ScalarVariable>( const char*, const FEM_Data<ScalarVariable>& );
template void Boundary<3>::InputVariableFrom<ScalarVariable>( const char*, const FEM_Data<ScalarVariable>& );

template void Boundary<1>::InputVariableFrom<ArrayVariable>( const char*, const FEM_Data<ArrayVariable>& );
template void Boundary<2>::InputVariableFrom<ArrayVariable>( const char*, const FEM_Data<ArrayVariable>& );
template void Boundary<3>::InputVariableFrom<ArrayVariable>( const char*, const FEM_Data<ArrayVariable>& );

template void Boundary<1>::InputVariableFrom<FlaggedArrayVariable>( const char*, const FEM_Data<FlaggedArrayVariable>& );
template void Boundary<2>::InputVariableFrom<FlaggedArrayVariable>( const char*, const FEM_Data<FlaggedArrayVariable>& );
template void Boundary<3>::InputVariableFrom<FlaggedArrayVariable>( const char*, const FEM_Data<FlaggedArrayVariable>& );

template void Boundary<1>::InputVariableFrom<VectorVariable<1U> >( const char*, const FEM_Data<VectorVariable<1U> >& );
template void Boundary<2>::InputVariableFrom<VectorVariable<2U> >( const char*, const FEM_Data<VectorVariable<2U> >& );
template void Boundary<3>::InputVariableFrom<VectorVariable<3U> >( const char*, const FEM_Data<VectorVariable<3U> >& );

template void Boundary<1>::InputVariableFrom<TensorVariable<1U> >( const char*, const FEM_Data<TensorVariable<1U> >& );
template void Boundary<2>::InputVariableFrom<TensorVariable<2U> >( const char*, const FEM_Data<TensorVariable<2U> >& );
template void Boundary<3>::InputVariableFrom<TensorVariable<3U> >( const char*, const FEM_Data<TensorVariable<3U> >& );






template<size_t dim>
bool Boundary<dim>::Divide( const Region<dim>& subsetRegion, Boundary<dim>& subsetBoundaryToForm )
{
  // make a map with faces of THIS boundary as value and their nodes as key
  map<set<Node<dim>*>, Face<dim>*>  thisFaces;
  set<Node<dim>*>                  faceNodes;
  for ( typename vector<Face<dim>*>::const_iterator it( this->ElementsBegin() ); it != this->ElementsEnd(); ++it )
  {
    faceNodes.clear();
    for ( size_t i = 0U; i< (*it)->Nodes(); i++ )
      faceNodes.insert( (*it)->N( i ) );
    thisFaces.insert( make_pair( faceNodes, (*it) ) );
  }
  faceNodes.clear();

  // make a set with sets of nodes of elements of subsetRegion
  set<set<Node<dim>*> > regionElementNodes;
  for ( typename vector<Element<dim>*>::const_iterator it( subsetRegion.ElementsBegin() ); it != subsetRegion.ElementsEnd(); ++it )
    regionElementNodes.insert( set<Node<dim>*>( (*it)->NodesBegin(), (*it)->NodesEnd() ) );

  // create a vector that contains those faces of THIS which have an equivalent key in the
  // region element nodes set --> simplexVector
  vector<Face<dim>*> subsetBoundaryToFormFaces;
  subsetBoundaryToFormFaces.reserve( subsetRegion.Elements() );
  for ( typename map<set<Node<dim>*>, Face<dim>*>::const_iterator it( thisFaces.begin() ); it != thisFaces.end(); ++it )
    if ( regionElementNodes.find( it->first ) != regionElementNodes.end() )
      subsetBoundaryToFormFaces.push_back( it->second );

  // create a method in Boundary that creates a boundary using a vector of faces
  // if empty(simplexVector) exit
  if ( subsetBoundaryToFormFaces.empty() )
    return false;

  return Divide( subsetBoundaryToFormFaces.begin(), subsetBoundaryToFormFaces.end(), subsetBoundaryToForm );
}




template<size_t dim>
bool Boundary<dim>::Divide( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                            const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                            Boundary<dim>& subsetBoundaryToForm )
{
  // use subsetBoundaryToForm.CreateFrom( simplexVector )
  subsetBoundaryToForm.CreateFrom( facesBegin, facesEnd, false /* do not update neighbor connectivity */, true /* update indexes */ );

  // remove faces in simplexVector from THIS boundary
  vector<Face<dim>*> subsetBoundaryToFormFaces( facesBegin, facesEnd );
  csmp::removeVectorElements( this->elmt_vec_, subsetBoundaryToFormFaces );

  if ( this->Elements() != 0 )
    Initialize( AtBoundary() /* box boundary flag */, false /* do not update neighbor connectivity*/, false /* do not update member indexes */ );

  return true;
}












// ------------------------------------------------------------------
// Building blocks
// -------------------------------------------------------------------




/**
Dispatched initialize method to establish node vector, perimeter entities,
bflags and entity sorting.

@attention this method will not change the BOX boundary flags of the nodes.
*/
template<size_t dim>
void Boundary<dim>::Initialize( BOX_BOUNDARY boxBoundary, bool updateNeighborConnectivity, bool updateIndexes )
{
  // establishing boundary node container
  this->CreateNodePointerVector();

  // identify entities on boundary perimeter
  if ( updateNeighborConnectivity )
    this->EstablishNeighborConnectivity( false );

  // sorts node and cell vectors into interior and exterior ranges; initialises boundary face vector bd_face_vec_
  this->IdentifyPerimeter();

  // initialize indices
  if ( updateIndexes ) this->UpdateMemberIndexes();

  // assigns boundary flags
  AtBoundary( boxBoundary );

} // end Initialize



  /*
  // TESTING - createLineFaceConnectivity()
  this->RenumberElements();
  this->RenumberNodes();
  cerr <<"\nBoundary<"<< dim <<">::Initialize: '"<< this->Name() <<"': printing element id(nodes): and neighbor ids";
  for ( auto it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); ++it )
  {
  cerr <<"\n\t"<< (*it)->Idx() <<" (";
  for ( size_t i=0U; i<(*it)->Nodes(); ++i ) cerr << (*it)->N(i)->Idx() <<",";
  cerr <<"): ";
  for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
  if ( (*it)->Neighbor(i) == nullptr ) cerr <<"nullptr ";
  else cerr << (*it)->Neighbor(i)->Idx() <<" ";
  }
  cerr << endl;
  */




template<size_t dim>
bool csmp::Boundary<dim>::IsExternal() const
{
    for ( typename vector<Face<dim>*>::const_iterator
          it=this->ElementsBegin(); it!=this->ElementsEnd(); ++it )
      // perhaps create method inside of Face to check whether it lies on the outside of the model
      if ( (*it)->OuterParent() == nullptr ) return false;
    return true;
}




/**
adds boundary flags to boundary nodes and elements if not already flagged as box boundary entity

@note changes BOX boundary flags for INTERIOR nodes along the boundary
*/
template<size_t dim>
void csmp::Boundary<dim>::AtBoundary( BOX_BOUNDARY boxBoundary )
{
  boundaryFlag_ = boxBoundary;
}



template<size_t dim>
void Boundary<dim>::CreateNodePointerVector()
{
  if ( this->elmt_vec_.empty() )
    throw csmp::Exception( ERROR, "Boundary<dim>::CreateNodePointerVector:",
                           this->Name(), "face vector is empty; nothing could be done." );

  if ( !this->node_vec_.empty() )
    this->node_vec_.clear();

  // creating the node index vector
  set<csmp::Node<dim>*>  nodes_set;
  for ( typename vector<Face<dim>*>::const_iterator it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); ++it )
    for ( size_t i = 0U; i<(*it)->Nodes(); ++i ) {
      assert( (*it)->N( i ) != nullptr );
      nodes_set.insert( (*it)->N( i ) );
    }

  this->node_vec_.assign( nodes_set.begin(), nodes_set.end() );
}





/**
Detects of how many spatial dimensions element types are contained in model.
It returns a pair: first value gives number of different spatial dimensions contained,
second value returns the highest spatial dimension contained.

@author SKM 1/11/2013
*/
template<size_t dim>
pair<int32, int32>  Boundary<dim>::FaceSpatialDimensions() const
{
  return this->SpatialDimensions();

} // end FaceSpatialDimensions








/**
Creates faces out of elements of supplied (lower dimensional) region

@author P. Lang - SKM refactored 2016
@date 30/8/2011

This method requires the supplied region to be of lower dimensional representation, throws otherwise.
Storage for both, the faces and boundary, is established. If a box boundary is supplied, member flags will be set accordingly;
whether it actually represents a valid box boundary is not checked for.

When establishing face connectivity(parent elements) the convention is as outlined:

@attention convention 1: If the element is on the model boundary(attached to a single higher dimensional parent element only),
the created face unit normal points outward.

@attention This method will work as well for a split boundary.
*/
template<size_t dim>
bool Boundary<dim>::CreateFrom( MeshManager<dim>& meshManager,
                                const Region<dim>& region,
                                const csmp::Index& mtrl_key, ///< to distinguish parent regions
                                BOX_BOUNDARY boxBoundary )
{
  // assert lower dimensional representation
  if ( !isOfLowerDimensionalRepresentation( region ) )
    throw csmp::Exception( ERROR, "Boundary<dim>::CreateFrom", "Region is not of lower dimensional representation." );

  // Thus the possible cases are:
  // 3D ( only surface elements ),
  // 3D ( surface + line elements),
  // 3D ( only line elements),
  // 2D ( only line elements )

  // for the case of a 3D model containing surface elements
  const bool regionContainsSurfaceElements( containsSurfaceElements( region ) );

  // LVS
  const LocalVariables lvsFaces( FaceVariables() );
  const IntegrationPointVariables lvsIntegrationPoints( FaceIntegrationPointVariables() );

  // prepping container for a max of total region element count
  this->elmt_vec_.reserve( region.Elements() );

  // looping over regions elements, assuring that it's an eligible face type, creating new face with variable storage,
  // establishing connectivity and inserting into boundary element container
  Face<dim>* root_face( NULL );
  const typename vector<Element<dim>*>::const_iterator regionElementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator it = region.ElementsBegin(); it != regionElementsEnd; ++it )
  {
    // in 3D, there still could be line elements in the region which are not eligible as face,
    // unless the Region consist only of line elements
    if ( dim == 3 && (*it)->IsLineElement() && regionContainsSurfaceElements )
      continue;

    // finding the higher-dimensional element that sits adjacent to the lower-dimensional one
    Face<dim>* faceObj( NULL );
    std::string region_name = region.Name();
    std::vector<csmp::Element<dim>*> inner_outter_elements;

    if( !higherDimensionalNeighbors( *(*it), inner_outter_elements ) ) continue;
    if( inner_outter_elements.size() == 2){
      // created a new face
      Face<dim> new_face( *(*it), inner_outter_elements[0], inner_outter_elements[1], lvsFaces, lvsIntegrationPoints );
      faceObj = meshManager.Add( new_face );
      // push back into face container
      this->elmt_vec_.emplace_back( faceObj );
    }
    else {
      // created a new face
      Face<dim> new_face( *(*it), const_cast<csmp::Element<dim>*>(inner_outter_elements[0]), nullptr, lvsFaces, lvsIntegrationPoints );
      faceObj = meshManager.Add( new_face );
      // push back into face container
      this->elmt_vec_.emplace_back( faceObj );
    }

    // the first face is assigned into the root face of this face group in the mesh
    if ( root_face == NULL ) {
      root_face = faceObj;
      meshManager.SetRootFace( faceObj );
    }

  } // region elements

    // free
  vector<Face<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  // initialize boundary essentials
  const bool update_nbor_connectivity( true ), update_member_indexes( true );
  Initialize( boxBoundary, update_nbor_connectivity, update_member_indexes );

  // done
  return true;

} // CreateFrom





/**
for post-processing the results of Divide
*/
template<size_t dim>
bool Boundary<dim>::CreateFrom( const typename vector<Face<dim>*>::const_iterator facesBegin,
                                const typename vector<Face<dim>*>::const_iterator facesEnd,
                                bool updateFaceConnectivity,
                                bool updateIndexes )
{

  this->elmt_vec_.assign( facesBegin, facesEnd );

  // initialize boundary essentials
  Initialize( AtBoundary(), updateFaceConnectivity, updateIndexes );

  return true;
}






/**
Creates faces around the region
Only for volume regions in 3D and surface regions in 2D.

@note any dim-2 elements are ignored: this means line elements
in a 3D model are ignored.

@date refactored by SKM 1/6/2016 (new face constructor etc.)

@author P. Lang, provisionally refactored SKM 7/6/2016
@date Aug 2011
*/
template<size_t dim>
bool Boundary<dim>::CreateAround( MeshManager<dim>& meshManager,
                                  const FiniteElementManager& finiteElementManager,
                                  const Region<dim>& region,
                                  BOX_BOUNDARY boxBoundary )
{
  // LVS
  const LocalVariables lvsFaces( FaceVariables() );
  const IntegrationPointVariables lvsIntegrationPoints( FaceIntegrationPointVariables() );

  // reserving storage for boundary elements
  // logic: the number of faces created cannot be larger than the number of boundary elements
  this->elmt_vec_.reserve( region.PerimeterElements() );

  // container for face nodes
  vector<size_t> faceNodes;

  // for the faces of the perimeter elements of the region
  Face<dim>* root_face( NULL );
  for ( size_t i = region.InteriorElements(); i<region.Elements(); ++i ) {
    // ignoring dim-2 elements
    if ( (dim == 3 and !region.E( i )->IsVolumeElement()) or
         (dim == 2 and !region.E( i )->IsSurfaceElement()) )
      continue;

    // for each perimter face
    const size_t perimeter_faces( region.PerimeterFaces( i ) );
    for ( size_t j = 0U; j<perimeter_faces; ++j )
    {
      // create the new face using the variables prepared above ( FV Stencil = NULL)
      Face<dim> new_face( *region.E( i ), finiteElementManager.E( region.E( i )->FE()->ElementTypeOfFace( region.PerimeterFace( i, j ) ) ), region.PerimeterFace( i, j ), lvsFaces, lvsIntegrationPoints );
      Face<dim>* faceObj = meshManager.Add( new_face );
      // push back into face container
      this->elmt_vec_.push_back( faceObj );

      // the first face is assigned into the root face of this face group in the mesh
      if ( root_face == NULL ) {
        root_face = faceObj;
        meshManager.SetRootFace( faceObj );
      }
    }
  }

  // free up excessive storage
  vector<Face<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  // initialize boundary essentials
  //  box boundary flag=true, update neighbor connectivity=true, update member indexes=true
  Initialize( boxBoundary, true, true );

  //done
  return true;

} // CreateAround


/**
@author P. Lang
@date 25/8/2011

This method assumes that elements are uniquely and throughgoingly numbered.
The Faces are build so that their normals point from region1 to region2.
Variable storage is assigned for both, the faces and boundary itself.

TODO: Fix! - this method does exactly (Face creation) what should be done by the MeshManager!

*/
template<size_t dim>
bool Boundary<dim>::CreateBetween( MeshManager<dim>& meshManager,
                                   const FiniteElementManager& finiteElementManager,
                                   const Region<dim>& region1,
                                   const Region<dim>& region2 )
{
  if ( string(region1.Name()) == region2.Name() )
    throw csmp::Exception( ERROR, "Boundary<dim>::CreateBetween:", "input region1 = input region2; nothing was done." );
    
  // LVS
  const LocalVariables lvsFaces( FaceVariables() );
  const IntegrationPointVariables lvsIntegrationPoints( FaceIntegrationPointVariables() );

  // reserving storage for boundary elements (logic: the number of faces created cannot be
  // larger than the minimum number boundary elements of the two neighboring groups)
  this->elmt_vec_.reserve( min( region1.PerimeterElements(), region2.PerimeterElements() ) );

  // searching for elements of region1 that are neighbors of ones in region2.
  // If so, there is a shared boundary and faces or interfaces are constructed.
  const size_t   n_elements( region1.Elements() );
  vector<size_t> fnids;

  Face<dim>* root_face( nullptr );
  // for the perimeter elements of the region
  for ( size_t i = region1.InteriorElements(); i < n_elements; ++i )
  {
    Element<dim>*  ePtr = region1.E( i );
    const size_t perimeter_faces( region1.PerimeterFaces(i) );
    for ( size_t j = 0U; j < perimeter_faces; ++j )
    {
      const size_t face = region1.PerimeterFace( i, j );
      Element<dim>*  ePtrNeighbor = ePtr->Neighbor( face );

// TODO: SKM DEBUGGING - perimeter faces are still not correctly identified by Region re-constructed from CSMP Binary
//cerr <<"\n"<< ePtr->Idx() <<": "<< parseBoundary( ePtr->AtBoundary() );
//if ( ePtr->AtBoundary() == NOT ) {
//      cout <<".";
//   }

      // checking whether neighbor element is part of the boundary of region2
      if ( ePtrNeighbor != nullptr )
        if ( region2.IsPerimeterElement( ePtrNeighbor ) )
        {
          // if the neighbor is in the boundary, the new Face is build
          FiniteElement* femPtr = finiteElementManager.E( ePtr->FE()->ElementTypeOfFace( face ) );
          Face<dim> new_face( femPtr, nullptr, lvsFaces, lvsIntegrationPoints );
          Face<dim>* faceObj = meshManager.AddIfUnique( new_face );

          // nodes are assigned to the new face
          ePtr->FE()->NodesOfFace( face, fnids );
          
          for ( size_t node = 0U; node < fnids.size(); ++node )
            faceObj->Assign( node, ePtr->N( fnids[node] ) );

          // the new face is connected to the elements it is sandwiched between
          // this assignment also includes connecting the face to its nodes
          //               inner        outer  element w.r.t. to normal of face
          faceObj->Assign( ePtr, ePtr->Neighbor( face ) );

          // added to boundary
          this->elmt_vec_.emplace_back( faceObj );

          // the first face is assigned into the root face of this face group in the mesh
          if ( root_face == nullptr ) {
            root_face = faceObj;
            meshManager.SetRootFace( faceObj );
          }
        } // neighboring elements

    } // perimeter faces

  } // perimeter elements

  // free
  vector<Face<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  // initialize boundary essentials
  Initialize( IRREGULAR, true, true );

  //done
  return true;

} // CreateBetween



  // CALCULATIONS

/**
    Returns the length of the perimeter line of the boundary.
    
    @return returns either the length of the perimeter of the boundary or NaN if the boundary is an edge.
    
    @attention method can only be applied in 3D and only on boundaries which are not lines themselves
    else the perimeter is not defined.
*/
template<size_t dim>
double64  Boundary<dim>::Perimeter() const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  if ( dim != 3U )
    csmp_error.notice( ERROR, "Boundary<dim>::Perimeter:",
                      "the perimeter of a Boundary is only defined when the boundary is a surface." );

  double64        perimeter_length( 0. );
  vector<size_t>  fnids;
  size_t          n( this->InteriorElements() );
  for ( typename vector<Face<dim>*>::const_iterator
        it = this->PerimeterElementsBegin(); it != this->ElementsEnd(); it++, n++ ) {
      if ( (*it)->IsLineElement() ) {
           return std::numeric_limits<double64>::quiet_NaN();
        }
      for ( size_t i = 0U; i<this->PerimeterFaces( n ); i++ ) {
        (*it)->FE()->NodesOfFace( this->PerimeterFace( n, i ), fnids );
        perimeter_length += ((*it)->N( fnids[1] )->Coordinate() -
                             (*it)->N( fnids[0] )->Coordinate()).Length();
      }
    }
    
  return perimeter_length;
}



template<size_t dim>
double64  Boundary<dim>::Area() const
{
  double64  integrated_area( 0. );

  for ( typename vector<Face<dim>*>::const_iterator
        it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
    integrated_area += (*it)->Area();

  return integrated_area;
}




template<size_t dim>
double64 Boundary<dim>::SurfaceIntegral( const PropertyDatabase<dim>& p, const char* property ) const
{
  csmp::Index prop_key = p.StorageKey( property );

  if ( prop_key.place == ELEMENT_INTEGRATION_POINT or prop_key.place == REGION ) {
    throw csmp::Exception( ERROR, "Boundary<dim>::SurfaceIntegral",
                           property, "placed on IntegrationPoint or Region cannot be assigned on boundary" );
    return std::numeric_limits<double64>::quiet_NaN();
  }
  if ( prop_key.type == TENSOR ) {
    throw csmp::Exception( ERROR, "Boundary<dim>::SurfaceIntegral",
                           property, "is a tensor property; this method does not know how to integrate it" );
    return std::numeric_limits<double64>::quiet_NaN();
  }

  double64  property_integral( 0. );

  // 1. if the property is a scalar
  if ( prop_key.type == SCALAR ) {
    if ( prop_key.place == ELEMENT ) {
      throw csmp::Exception( ERROR, "Boundary<dim>::SurfaceIntegral",
                             property, "is an Element property; this method does not know how to integrate it" );
    }
    else if ( prop_key.place == FACE ) {
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
        property_integral += (*it)->Volume() * (*it)->Read( prop_key );
    }
    else if ( prop_key.place == NODE ) { // for nodes on first side of interface
      ScalarVariable  sc;
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ ) {
        (*it)->PropertyValueAtBaryCenter( prop_key, sc );
        property_integral += (*it)->Volume() * sc();
      }
    }
    else {
      throw csmp::Exception( FATAL_ERROR, "Boundary<dim>::SurfaceIntegral",
                             "Property placement not recognized" );
    }
  }
  // 1. if the property is a scalar
  if ( prop_key.type == VECTOR ) {
    // the average of the values projected onto the normal are being used.
    if ( prop_key.place == FACE or prop_key.place == INTER_FACE ) {
      VectorVariable<dim>  unrml, vc;
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ ) {
        (*it)->UnitNormal( unrml );
        (*it)->Read( prop_key, vc );
        property_integral += dotProduct( unrml, vc );
      }
    }
    else if ( prop_key.place == NODE ) { // for nodes on first side of interface
      VectorVariable<dim>  unrml, vc;
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ ) {
        (*it)->UnitNormal( unrml );
        (*it)->PropertyValueAtBaryCenter( prop_key, vc );
        property_integral += dotProduct( unrml, vc );
      }
    }
    else throw csmp::Exception( FATAL_ERROR, "Boundary<dim>::SurfaceIntegral",
                                "Property placement not recognized" );
  }

  return property_integral;

} // end SurfaceIntegral







  // SCREEN OUTPUT




template<size_t dim>
void Boundary<dim>::Out() const
{
  // high-level output
  cout << "\n" << "Boundary<" << dim << ">::Out: (" << parseBoundary( boundaryFlag_ ) << ") '" << this->Name();
  cout << "', Face objects interior: " << this->InteriorElements() << ", perimeter: " << this->PerimeterElements() << endl;
  cout << "   Node objects interior: " << this->InteriorNodes() << ", perimeter: " << this->PerimeterNodes() << endl;

  // member faces
  cout << "\n\tFace objects, their area, nodes (,), and lower-(:) and higher-dimensional neighbor elements:\n";
  for ( const auto it : this->elmt_vec_ ) {
    if ( it == NULL ) throw csmp::Exception( ERROR, "Boundary<dim>::Out:", "member element pointer not initialized." );
    cout << "\t\t" << it->Idx() << ": " << it->Area() << ", ";
    // Nodes
    for ( size_t i = 0U; i<it->Nodes(); ++i ) {
      const string str = (it->N( i ) == nullptr) ? "none" : to_string( it->N( i )->Idx() );
      cout << str << ",";
    }
    cout << "\t\t ";
    // Face neighbors
    for ( size_t i = 0U; i<it->Neighbors(); ++i ) {
      const string str = (it->Neighbor( i ) == nullptr) ? "none" : to_string( it->Neighbor( i )->Idx() );
      if ( i<it->Neighbors() - 1 ) cout << str << ":";
      else cout << str;
    }
    cout << ",\t\t";
    // Element neighbors
    const string inner = (it->InnerParent() == nullptr) ? "none" : to_string( it->InnerParent()->Idx() );
    const string outer = (it->OuterParent() == nullptr) ? "none" : to_string( it->OuterParent()->Idx() );
    cout << inner << ":" << outer << "\n";
  }

  // printing the Faces
  //for ( auto it=this->ElementsBegin(); it!=this->ElementsEnd(); ++it )  (*it)->Out();

  cout << "\n\tperimeter Faces and edge numbers (current local numbering):\n";
  vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit( this->bd_face_vec_.begin() );
  for ( size_t i = this->InteriorElements(); i<this->elmt_vec_.size(); i++, bit++ ) {
    cout << i << ":";
    for ( vector<ONE_BYTE_NUMBER>::const_iterator
          ft = (*bit).begin(); ft != (*bit).end(); ft++ ) cout << (*ft) << " ";
  }
  cout << endl;
}


template class Boundary<1>;
template class Boundary<2>;
template class Boundary<3>;

} // end csmp
