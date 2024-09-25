#include "Boundary.h"
#include "Region.h"
#include "Box.h"

#include "writeVariableIf.h"
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
#include "MeshManagementUtilities.h"
#include "PL_Utilities.h"
#include "variableOperations.h"

#include "Visitor.h"
#include "FaceConstructionData.h"

using namespace std;

namespace csmp {

/**
Constructor for an empty boundary (no faces yet).
*/
template<uint32_t dim>
Boundary<dim>::Boundary( std::string boundaryname, const PropertyDatabase<dim>& pref, BOX_BOUNDARY flag )
  : ModelSubDomain<dim, Face>( boundaryname, pref ),
    boundaryFlag_( flag )
{
  if constexpr ( dim == 1 )
    ErrorHandler::Instance().Note( ERROR, "Boundary<dim>::(custom constructor):", "boundary objects are only supported in 2 & 3D models." );
  this->ResizePropertyStorage( this->pref_.LocalVariablesAt( Placement() ) );
}


/// copy constructor
template<uint32_t dim>
Boundary<dim>::Boundary( const Boundary& ed )
  : ModelSubDomain<dim, Face>( ed ),
    LocalVariableStorage<dim,Boundary>(ed),
    boundaryFlag_( ed.boundaryFlag_ )
{
}


/// move constructor
template<uint32_t dim>
Boundary<dim>::Boundary( Boundary&& ed )
  : ModelSubDomain<dim, Face>( std::move(ed) ),
    LocalVariableStorage<dim,Boundary>( std::move(ed) ),
    boundaryFlag_( std::move(ed.boundaryFlag_) )
{
}



/**

RECONSTRUCTOR of boundary from the mesh manager in the model

@attention expects that the face numbering is unique and that all faces
are numbered consecutivly.

@author SKM 4/5/2016

@attention assumes that the index for the first Face is equivalent to the number of elements in the model.

@test nodes are assigned correctly, Face objects are not
*/
template<uint32_t dim>
Boundary<dim>::Boundary( const PropertyDatabase<dim>& pref,
                         MeshManager<dim>& mesh,
                         const SubDomainInfo& info,
                         BOX_BOUNDARY bflag )
  : ModelSubDomain<dim, Face>( info.name, pref ),
    boundaryFlag_( bflag )
{
  // building the face vector
  // ------------------------
  this->cell_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );

  // assigning pointers to the interior faces
  const size_t  n_elements{ mesh.Elements() }; // needs to be subtracted accessing the face container
  for ( size_t i : info.interior_elmts )
    this->cell_vec_.push_back( &(*next(mesh.FacesBegin(),i-n_elements)) );

  // assigning pointers to the perimeter faces
  for ( size_t i : info.perimeter_elmts )
    this->cell_vec_.push_back( &(*next(mesh.FacesBegin(),i-n_elements)) );

  // building the node vector
  // ------------------------
  // assigning pointers to the interior and perimeter nodes
  this->first_bd_node_ = info.interior_nodes.size();
  this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );

  // assigning pointers to the interior faces
  for ( size_t i : info.interior_nodes )
    this->node_vec_.push_back( &(*next(mesh.NodesBegin(),i)) );

  // assigning pointers to the perimeter faces
  for ( size_t i : info.perimeter_nodes )
    this->node_vec_.push_back( &(*next(mesh.NodesBegin(),i)) );

  this->SortVectors( info.interior_elmts.size(), info.interior_nodes.size() );

  // building the vector of vectors of those faces (edges) of the faces that lie on the subdomain perimeter
  this->BuildPerimeterFaceVector( info.interior_elmts.size() );

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( BOUNDARY ) );

} // end re-constructor (MeshManager)







template<uint32_t dim>
IntegrationPointVariables Boundary<dim>::FaceIntegrationPointVariables() const
{ return this->pref_.IntegrationPointVariablesAt( FACE ); }

template<uint32_t dim>
LocalVariables Boundary<dim>::FaceVariables() const
{ return this->pref_.LocalVariablesAt( FACE ); }





/**
Creates boundary from supplied vector of faces that must already
know their higher-dimensional parent elements and their equidimensional
neighbor faces.

@attention the neighborhood relations are not re-established.

@note this constructor is mainly used for edges of line-element faces.
*/
template<uint32_t dim>
Boundary<dim>::Boundary( const string& boundary_name,
                         const PropertyDatabase<dim>& dbase_ref,
                         typename std::vector<Face<dim>*>::iterator facesBegin,
                         typename std::vector<Face<dim>*>::iterator facesEnd,
                         BOX_BOUNDARY boxBoundary )
  : ModelSubDomain<dim, Face>( boundary_name, dbase_ref ),
    boundaryFlag_( boxBoundary )
{
  // moving the supplied Face pointers into the element storage
  this->cell_vec_.assign( facesBegin, facesEnd );
  
  // initialize BOX_BOUNDARY of nodes
  InitializeBoundaryFlags( boxBoundary );
  
  // distinguishing interior from perimeter cells
  // (sorts node and cell vectors into interior and exterior ranges;
  //  initialises boundary face vector bd_face_vec_)
  this->IdentifyPerimeter();
  
  // trimming off extra capacity
  this->cell_vec_.shrink_to_fit();
  this->node_vec_.shrink_to_fit();
}




/**
Default dtor does not delete faces, even though constructed by Boundary. Need to call DeleteFaces explicitly!
@todo FIX THIS CONCEPTUAL MISUNDERSTANDING ALL THE WAY THROUGH:  the faces are managed by the MeshManager
*/
template<uint32_t dim>
Boundary<dim>::~Boundary()
{
}

template<uint32_t dim>
Boundary<dim>& Boundary<dim>::operator=( const Boundary<dim>& ed )
{
  if ( &ed != this ) {
      ModelSubDomain<dim,Face>::operator=( ed ); 
      this->LVS( ed.LVS() );
      boundaryFlag_ = ed.boundaryFlag_;
    }
  return *this;
}

// LOCAL VARIABLE STORAGE INTERFACE

template<uint32_t dim>
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
template<uint32_t dim>
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
      for ( auto it = this->CellsBegin(); it != this->CellsEnd(); it++ )
        (*it)->Accept( v );
      return;
    case NODE:
      for ( auto nd_it = this->NodesBegin(); nd_it != this->NodesEnd(); nd_it++ )
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
   for the assignment of properties that are unique to the instance of this subclass
   
      @author SKM
      @date 7/6/2020
*/
template<uint32_t dim>
template<typename Var>
void Boundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index prop_key(this->pref_.StorageKey(input_prop));
      if ( prop_key.place == BOUNDARY ) {
           if ( sd != COMPLETE )
             csmp_error.Note( WARNING, "Boundary<dim>::InputPropertyValue",
                              input_prop, "is a BOUNDARY property and no distinction between INTERIOR and PERIMETER can be made" );
           this->Store( prop_key, new_value );
           return;
        }
      
      // incorrect applications of method  
      if ( prop_key.place == REGION || prop_key.place == SPLIT_BOUNDARY || prop_key.place == MODEL || 
           prop_key.place == ELEMENT || prop_key.place == INTER_FACE )
        csmp_error.Note( ERROR, "Boundary<dim>::InputPropertyValue",
                           input_prop, "must be a BOUNDARY, FACE/IP or NODE property for this method to work" );
    
     // for any different property placement, the method of the base-class is called
     ModelSubDomain<dim,Face>::InputPropertyValue( input_prop, new_value, sd );
      
  } // end InputPropertyValue

// explicit instantiations
template void Boundary<1U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const VectorVariable<1U>&, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const VectorVariable<2U>&, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const VectorVariable<3U>&, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const TensorVariable<1U>&, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const TensorVariable<2U>&, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const TensorVariable<3U>&, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );




template<uint32_t dim>
template<typename Var>
void Boundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index key(this->pref_.StorageKey(input_prop));
      
      if ( key.place == BOUNDARY ) {
           if ( sd != COMPLETE )
             csmp_error.Note( WARNING, "Boundary<dim>::InputPropertyValue",
                                input_prop, "is a BOUNDARY property and no distinction between INTERIOR and PERIMETER can be made" );
                                
           // only overwriting those variable components / rows that are not flagged 'do_not_overwrite' 
           writeVariableIf( this, key, new_value, do_not_overwrite );
           return;
        }
      
      // incorrect applications of method  
      if ( key.place == REGION  || key.place == SPLIT_BOUNDARY || key.place == MODEL || 
           key.place == ELEMENT || key.place == INTER_FACE )
        csmp_error.Note( ERROR, "Boundary<dim>::InputPropertyValue",
                           input_prop, "must be a BOUNDARY, FACE/IP or NODE property for this method to work" );        
    
     // for any different property placement, the method of the base-class is called
     ModelSubDomain<dim,Face>::InputPropertyValue( input_prop, new_value, do_not_overwrite, sd );
      
  } // end InputPropertyValue

// explicit instantiations
template void Boundary<1U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const VectorVariable<1U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const VectorVariable<2U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const VectorVariable<3U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const TensorVariable<1U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const TensorVariable<2U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const TensorVariable<3U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<1U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<2U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Boundary<3U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );





/**
Outs FEM_Data of provided variable type and placement to binary

@author  P. Lang
@date  9/30/2012

@param [in,out]  fp  If non-null, the binary file pointer
@param place         The placement of variables
@param vtype         The vtype of variables

@return  true if it succeeds, false if it fails.
*/
template<uint32_t dim>
template<class Var>
bool Boundary<dim>::Out( fstream& fp, PLACEMENT place, VARIABLE_TYPE vtype ) const
{
  set<string> propList;
  size_t vCount( 0 ), bytes( sizeof( size_t ) );

  this->pref_.ListProperties( place, vtype, propList );
  vCount = propList.size();
  fp.write( reinterpret_cast<const char*>(&vCount), bytes );
  if ( vCount != 0 )
  {
    FEM_Data<Var> femData;
    for ( set<string>::const_iterator it( propList.begin() ); it != propList.end(); ++it )
    {
      if ( !binaryFileWrite( fp, it->c_str() ) )
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
the outer template provides double = data type and dim = dimension,
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
template<uint32_t dim>
template<class Var>
void Boundary<dim>::OutputVariableTo( const char* property, FEM_Data<Var>& data ) const
{
  csmp::Index  idx = this->pref_.StorageKey( property );
  Var          var;
  femDataOutputDispatch::initVariable( idx, var );

  switch ( idx.place ) {
    case FACE: {
      data.Reset( idx, this->cell_vec_.size(), var );
      for ( typename vector<csmp::Face<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ ) {
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
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
      {
        for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ ) {
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
      csmp_error.Note( WARNING, "Boundary::OutputVariableTo",
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
template<uint32_t dim>
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
    if ( !binaryFileRead( fp, propertyName ) )
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
template<uint32_t dim>
template<class Var>
void Boundary<dim>::InputVariableFrom( const char* property,
                                       const FEM_Data<Var>& vdata )
{
  const csmp::Index  idx = this->pref_.StorageKey( property );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  switch ( idx.place ) {
    case FACE:
      //assert( this->Cells() == vdata.Size() );
      for ( typename vector<csmp::Face<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
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
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
      {
        for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
          (*eit)->Store( i, idx, vdata[counter + i] );
        counter += (*eit)->IntegrationPoints();
      }
    }
    case BOUNDARY:
      if ( vdata.Size() != 1U )
        csmp_error.Note( WARNING, "Boundary<dim>::InputVariableFrom:",
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
















// ------------------------------------------------------------------
// Building blocks
// -------------------------------------------------------------------




/**
    Dispatched initialize method to establish node vector, perimeter entities,
    bflags and entity sorting.
*/
template<uint32_t dim>
void Boundary<dim>::InitializeBoundaryFlags( BOX_BOUNDARY boxBoundary )
{
  // assigns BOX_BOUNDARY flag to Boundary
  AtBoundary( boxBoundary );

  // establishing boundary node container
  this->CreateNodePointerVector();
  
  // assigning the box boundary flag to nodes where they do not already have another flag (like EDGE etc)
  for ( auto& nit : this->node_vec_ )
    if ( nit->AtBoundary() == NOT )
      nit->AtBoundary( boxBoundary );

} // end InitializeBoundaryFlags





/**
      If one of the Face objects in the boundary has no outer higher-dim neighbor,
      the boundary is located on the model boundary.
      
      @return returns true when part or all of the boundary is located on the model boundary.
*/
template<uint32_t dim>
bool csmp::Boundary<dim>::IsExternal() const
{
    for ( const auto& it : this->cell_vec_ )
      // perhaps create method inside of Face to check whether it lies on the outside of the model
      if ( it->OuterParent() == nullptr ) return true;
    return false;
}




/**
adds boundary flags to boundary nodes and elements if not already flagged as box boundary entity

@note changes BOX boundary flags for INTERIOR nodes along the boundary
*/
template<uint32_t dim>
void csmp::Boundary<dim>::AtBoundary( BOX_BOUNDARY boxBoundary )
{
  boundaryFlag_ = boxBoundary;
}




/**
Detects of how many spatial dimensions element types are contained in model.
It returns a pair: first value gives number of different spatial dimensions contained,
second value returns the highest spatial dimension contained.

@author SKM 1/11/2013
*/
template<uint32_t dim>
pair<int32_t, int32_t>  Boundary<dim>::FaceSpatialDimensions() const
{
  return this->SpatialDimensions();

} // end FaceSpatialDimensions








/**
Creates faces matching the elements of the supplied (lower dimensional) region.

@author P. Lang - SKM refactored 2016, 2021
@date 30/8/2011

This method requires the supplied region to be of lower dimensional representation, throws otherwise.
Storage for both, the faces and boundary, is established. If a box boundary is supplied, the flag of the Boundary will be set accordingly;
whether it actually represents a valid box boundary is not checked for.

When establishing face connectivity(parent elements) the convention is as outlined:

@attention If the element is on the model boundary(attached to a single higher dimensional parent element only),
the unit normal of the created face will point outward.

*/
template<uint32_t dim>
bool Boundary<dim>::CreateFrom( Region<dim>& region,
                                MeshManager<dim>& meshManager,
                                BOX_BOUNDARY boxBoundary )
{
throw csmp::Exception( ERROR, "Boundary<dim>::CreateFrom", "BROKEN: fix before using this method" );

  // assert lower dimensional representation
  if ( !hasLowerDimensionalRepresentation( region ) )
    throw csmp::Exception( ERROR, "Boundary<dim>::CreateFrom", "Region is not of lower dimensional representation." );

  // Thus the possible cases are:
  // 3D ( only surface elements ),
  // 3D ( surface + line elements),
  // 3D ( only line elements),
  // 2D ( only line elements )

  // LVS
  const LocalVariables lvsFaces( FaceVariables() );
  const IntegrationPointVariables lvsIntegrationPoints( FaceIntegrationPointVariables() );

  // prepping container for a max of total region element count
  this->cell_vec_.reserve( region.Cells() );

  // looping over regions elements, assuring that it's an eligible face type, creating new face with variable storage,
  // establishing connectivity and inserting into boundary element container
  const auto regionElementsEnd( region.CellsEnd() );
  size_t     cell_number{0};
  
  for ( auto it = region.CellsBegin(); it != regionElementsEnd; ++it )
    {
      // renumbering the elements so that they can be used to find neighbors
      (*it)->Idx( cell_number );
      
      // in 3D, there still could be line elements in the region which are not eligible to become faces,
      // unless the Region consist only of line elements
      if constexpr ( dim == 3u ) if ( (*it)->IsLine() ) continue;

      // finding the higher-dimensional element(s) that sit(s) adjacent to the lower-dimensional one
      auto new_nbors = parentElements<dim>( (*it)->NodesBegin(), (*it)->NodesEnd() );
      
      if ( new_nbors.first.first != nullptr && new_nbors.second.first != nullptr )
        {
          // create new internal face
          this->cell_vec_.push_back( meshManager.ReplaceElementByFace( (*it),
                                                                       new_nbors.first.first, new_nbors.second.first,
                                                                       new_nbors.first.second, new_nbors.second.second,
                                                                       lvsFaces, lvsIntegrationPoints ) );
        }
      else { // if this Face will be located at the model boundary
          // create new boundary face
          this->cell_vec_.push_back( meshManager.AddBoundaryFace( new_nbors.first.first,
                                                                  new_nbors.first.second,
                                                                  lvsFaces, lvsIntegrationPoints ) );
        }
      
       // giving the new Face the same number as that of the element
       this->cell_vec_.back()->Idx( cell_number++ );
       
    } // region elements


  // connecting the new faces to their neighbors
  size_t elmt{0};
  // processing the interior faces whose neighbors are all on the inside of the domain first
  for ( auto it=region.CellsBegin(); it!=region.PerimeterCellsBegin(); ++it ) {
       const size_t n_neighbors{ (*it)->Neighbors() };
       for ( auto i{0U}; i<n_neighbors; ++i ) {
            // since the mapping between elements and faces only exists in this subdomain
            // outside elements cannot be considered
            assert ( (*it)->Neighbor(i) != nullptr );
            this->cell_vec_[elmt]->Assign( i, this->cell_vec_[ (*it)->Neighbor(i)->Idx()] );
         }
       elmt++;
    }
  // processing perimeter faces
  elmt = region.InteriorCells();
  for ( auto it=region.PerimeterCellsBegin(); it!=region.CellsEnd(); ++it ) {
       const size_t n_neighbors{ (*it)->Neighbors() };
       for ( auto i{0U}; i<n_neighbors; ++i )
         // if the perimeter element neighbor is contained in the interior elements of region an assignment is made
         if ( (*it)->Neighbor(i) != nullptr && region.Contains( (*it)->Neighbor(i) ) )
           this->cell_vec_[elmt]->Assign( i, this->cell_vec_[ (*it)->Neighbor(i)->Idx()] );
       elmt++;
    }

  // initialize BOX_BOUNDARY of nodes
  InitializeBoundaryFlags( boxBoundary );

  // distinguishing interior from perimeter cells
  // (sorts node and cell vectors into interior and exterior ranges;
  //  initialises boundary face vector bd_face_vec_)
  this->IdentifyPerimeter();

  // done
  return true;

} // CreateFrom


// DEBUG - can all elements of the region be reached by a floodfill? (NOPE)
//set<Face<dim>*>  cells_contiguous_subset;
//findContiguousMeshPatch( *(this->cell_vec_[0]), cells_contiguous_subset );
//assert( cells_contiguous_subset.size() == this->cell_vec_.size() );
// DEBUG - neighbors after assignment
/*
cerr <<"\n\nRegion "<< region.Name() <<"\n";
for ( auto& it : region.CellVector() ) {
    cerr <<"\n\tElement "<< it->Idx() <<": nbors: ";
    for ( auto i{0U}; i<it->Neighbors(); ++i )
      if ( it->Neighbor(i) == nullptr ) cerr <<" null";
      else cerr <<" "<< it->Neighbor(i)->Idx();
  }
cerr << endl;

cerr <<"\n\nBoundary "<< region.Name() <<"\n";
for ( auto& it : this->cell_vec_ ) {
    cerr <<"\n\tFace "<< it->Idx() <<": nbors: ";
    for ( auto i{0U}; i<it->Neighbors(); ++i )
      if ( it->Neighbor(i) == nullptr ) cerr <<" null";
      else cerr <<" "<< it->Neighbor(i)->Idx();
  }
cerr << endl;
*/



/**
      Accumulates existing Face objects, identified by their ID number into a new Boundary object.
      The Face ids are expected to be in the the range between n_elements and n_interfaces-1.
      
      Method assumes that face numbers follow consecutively on the element numbers and before the interface numbers in the mesh.
      Thus, the first Face is expected to have the index  'n_elements'
*/
template<uint32_t dim>
size_t Boundary<dim>::AccumulateByNumber( MeshManager<dim>& mesh,
                                          vector<size_t>& cell_ids )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( cell_ids.empty() )
    csmp_error.Note( ERROR, "Boundary<dim>::AccumulateByNumber",
                      "user-supplied face-number vector is empty. Nothing is done." );

  if ( !this->cell_vec_.empty() ) {
      csmp_error.Note( WARNING, "Boundary<dim>::AccumulateByNumber",
                         "Boundary is not empty", "erasing all members..." );
      this->cell_vec_.clear();
    }

  // eliminating potential duplicates from element index vector
#ifdef DEBUG
  const size_t n_cells{cell_ids.size()};
  sort( cell_ids.begin(), cell_ids.end() );
  cell_ids.erase( unique( cell_ids.begin(), cell_ids.end() ), cell_ids.end() );
  if ( cell_ids.size() < n_cells )
    csmp_error.Note( WARNING, "Boundary<dim>::AccumulateByNumber",
                      "user-supplied face ID set contained duplicates which were removed." );
#endif
  // checking that the required faces already exist
  if ( mesh.Faces() == 0 )
    csmp_error.Note( FATAL_ERROR, "Boundary<dim>::AccumulateByNumber",
                       "n_faces=0! - This method expects that the Faces that shall be accumulated have already been created; please use other face creation method." );

  // checking that there are not more ids than there are faces in the model
  if ( cell_ids.size() > mesh.Faces() )
    csmp_error.Note( ERROR, "Boundary<dim>::AccumulateByNumber",
                       "user-supplied face-number vector is larger than range of index-to-element-pointer mapping." );

  // creating the face vector for the boundary
  // NB: assumes that the Faces are numbered consecutively from 0..n-1, while the supplied IDs start at the number of elements
  const auto offset = mesh.Elements();
  this->cell_vec_.reserve( cell_ids.size() );
  for ( auto& idx : cell_ids ) {
       const auto face = idx - offset;
       // TODO: this method makes no Faces, but depends on MeshManager for this task
       assert( mesh.Faces() != 0 );
       assert( face < mesh.Faces() );
       Face<dim>* fptr = &(*next(mesh.FacesBegin(),face));
       assert( fptr != nullptr );
       assert( fptr->Idx() == idx );
       this->cell_vec_.push_back( fptr );
    }
    
  // creating node vector
  if ( this->node_vec_.empty() ) this->node_vec_.clear();
  this->node_vec_.reserve( cell_ids.size() ); // just a loose measure, asuming that there will always be more elements than nodes
  // filling the vector
  for ( auto& it : this->cell_vec_ ) {
       const auto n_nodes{it->Nodes()};
       for ( auto i{0U}; i<n_nodes; ++i ) {
            assert( it->N(i) != nullptr );
            this->node_vec_.push_back( it->N(i) );
         }
     }
  // removing duplicates and trimming excess memory from node vector
  sort( this->node_vec_.begin(), this->node_vec_.end() );
  this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );

  this->IdentifyPerimeter();
  
  return this->cell_vec_.size();

} // end AccumulateByNumber





/**
    To create Boundary from an existing range of Face objecgts. It is assumed
    that the faces are already interconnected.
*/
template<uint32_t dim>
bool Boundary<dim>::CreateFrom( const typename vector<Face<dim>*>::const_iterator facesBegin,
                                const typename vector<Face<dim>*>::const_iterator facesEnd )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
  if ( distance(facesBegin,facesEnd) == 0 ) {
       csmp_error.Note( ERROR, "Boundary<dim>::CreateFrom", "supplied Face range is empty; nothing was done");
       return false;
    }
  if ( (*facesBegin)->ConnectedNeighbors() == 0 ) {
       csmp_error.Note( ERROR, "Boundary<dim>::CreateFrom", "some supplied Face objects do not have neighbors; nothing was done");
       return false;
    }
    
  this->cell_vec_.assign( facesBegin, facesEnd );

  // initialize BOX_BOUNDARY of nodes
  InitializeBoundaryFlags( boundaryFlag_ );
  
  // distinguishing interior from perimeter cells
  // (sorts node and cell vectors into interior and exterior ranges;
  //  initialises boundary face vector bd_face_vec_)
  this->IdentifyPerimeter();

  return true;
}







// ==========================================================================
//
// CALCULATIONS
//
// ==========================================================================



/**
    Returns the length of the perimeter line of the boundary.
    
    @return returns either the length of the perimeter of the boundary or NaN if the boundary is an edge.
    
    @attention method can only be applied in 3D and only on boundaries which are not lines themselves
    else the perimeter is not defined.
*/
template<uint32_t dim>
double  Boundary<dim>::Perimeter() const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  if ( dim != 3U )
    csmp_error.Note( ERROR, "Boundary<dim>::Perimeter:",
                      "the perimeter of a Boundary is only defined when the boundary is a surface." );

  double perimeter_length( 0. );
  size_t n( this->InteriorCells() );
  
  for ( auto it = this->PerimeterCellsBegin(); it != this->CellsEnd(); it++, n++ ) {
      if ( (*it)->IsLine() ) {
           return std::numeric_limits<double>::quiet_NaN();
        }
      for ( auto i{0U}; i<this->PerimeterFaces( n ); i++ ) {
        auto fnids = (*it)->FE()->NodesOfFace( this->PerimeterFace( n, i ) );
        perimeter_length += ((*it)->N( fnids[1] )->Coordinate() -
                             (*it)->N( fnids[0] )->Coordinate()).Length();
      }
    }
    
  return perimeter_length;
}


/**
    Sums the length of line elements in the boundary, if any.
*/
template<uint32_t dim>
double  Boundary<dim>::Length() const
{
  double  integrated_area( 0. );

  for ( const auto& it  : this->cell_vec_ )
    if ( it->IsLine() )
      integrated_area += it->Area();

  return integrated_area;
}


/**
    Sums the area of surface elements in the boundary, if any.
*/
template<uint32_t dim>
double  Boundary<dim>::Area() const
{
  double  integrated_area( 0. );

  for ( const auto& it  : this->cell_vec_ )
    if ( it->IsSurface() )
      integrated_area += it->Area();

  return integrated_area;
}




template<uint32_t dim>
double Boundary<dim>::SurfaceIntegral( const PropertyDatabase<dim>& p, const char* property ) const
{
  csmp::Index prop_key = p.StorageKey( property );

  if ( prop_key.place == ELEMENT_INTEGRATION_POINT or prop_key.place == REGION ) {
    throw csmp::Exception( ERROR, "Boundary<dim>::SurfaceIntegral",
                           property, "placed on IntegrationPoint or Region cannot be assigned on boundary" );
    return std::numeric_limits<double>::signaling_NaN();
  }
  if ( prop_key.type == TENSOR ) {
    throw csmp::Exception( ERROR, "Boundary<dim>::SurfaceIntegral",
                           property, "is a tensor property; this method does not know how to integrate it" );
    return std::numeric_limits<double>::signaling_NaN();
  }

  double  property_integral( 0. );

  // 1. if the property is a scalar
  if ( prop_key.type == SCALAR ) {
    if ( prop_key.place == ELEMENT ) {
      throw csmp::Exception( ERROR, "Boundary<dim>::SurfaceIntegral",
                             property, "is an Element property; this method does not know how to integrate it" );
    }
    else if ( prop_key.place == FACE ) {
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->cell_vec_.begin(); it != this->cell_vec_.end(); it++ )
        property_integral += (*it)->Volume() * (*it)->Read( prop_key );
    }
    else if ( prop_key.place == NODE ) { // for nodes on first side of interface
      ScalarVariable  sc;
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->cell_vec_.begin(); it != this->cell_vec_.end(); it++ ) {
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
            it = this->cell_vec_.begin(); it != this->cell_vec_.end(); it++ ) {
        (*it)->UnitNormal( unrml );
        (*it)->Read( prop_key, vc );
        property_integral += dotProduct( unrml, vc );
      }
    }
    else if ( prop_key.place == NODE ) { // for nodes on first side of interface
      VectorVariable<dim>  unrml, vc;
      for ( typename vector<Face<dim>*>::const_iterator
            it = this->cell_vec_.begin(); it != this->cell_vec_.end(); it++ ) {
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




template<uint32_t dim>
void Boundary<dim>::Out() const
{
  // high-level output
  cout << "\n"<<"Boundary<" << dim << ">::Out: (" << parseBoundary( boundaryFlag_ ) << ") '" << this->Name();
  if constexpr ( dim == 2U ) cout <<"', length "<< Area() <<" m";
  if constexpr ( dim == 3U ) cout <<"', area "<< Area() <<" m2";
  cout <<"\n\t"<<"Face objects interior: " << this->InteriorCells() << ", perimeter: " << this->PerimeterCells() << endl;
  cout <<"\n\t"<<"Node objects interior: " << this->InteriorNodes() << ", perimeter: " << this->PerimeterNodes() << endl;

  // member faces
  cout << "\n\tFace type followed by number: nodes (,), (inner:outer) parent elements, and (::) neighbor faces:\n";
  for ( const auto& it : this->cell_vec_ ) {
      if ( it == nullptr ) throw csmp::Exception( ERROR, "Boundary<dim>::Out:", "nullptr in face vector." );
      cout <<"\t\t"<< parseFiniteElementType( it->FE_Type() ) <<" "<< it->Idx() << ": ";
      // Nodes
      for ( auto i{0U}; i<it->Nodes(); ++i ) {
          const string str = (it->N(i) == nullptr) ? "none" : to_string( it->N(i)->Idx() );
          cout << str << ",";
        }
      cout << "\t\t";
      // Element neighbors
      const string inner = (it->InnerParent() == nullptr) ? "none" : to_string( it->InnerParent()->Idx() );
      const string outer = (it->OuterParent() == nullptr) ? "none" : to_string( it->OuterParent()->Idx() );
      cout << inner << ":" << outer;
      cout << ",\t\t";
      // Face neighbors
      for ( auto i{0U}; i<it->Neighbors(); ++i ) {
        const string str = (it->Neighbor(i) == nullptr) ? "none" : to_string( it->Neighbor(i)->Idx() );
        if ( i<it->Neighbors() - 1U ) cout << str << ":";
        else cout << str;
      }
     cout << endl;
  }

  // printing the Faces
  //for ( auto it=this->CellsBegin(); it!=this->CellsEnd(); ++it )  (*it)->Out();

  cout << "\n\tperimeter Faces and perimeter-edge numbers (local numbering):\n";
  auto  bit( this->bd_face_vec_.begin() );
  for ( auto i = this->InteriorCells(); i<this->cell_vec_.size(); i++, bit++ ) {
      cout << i << ":";
      for ( auto ft = (*bit).begin(); ft != (*bit).end(); ft++ ) cout << (*ft) << " ";
    }
  cout << endl;
  
} // end


template class Boundary<1>;
template class Boundary<2>;
template class Boundary<3>;

} // end csmp
