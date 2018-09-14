#include "FiniteElementPolicy.h"
#include "FiniteElement.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "variableOperations.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"

namespace csmp {

template<size_t dim, template<size_t> class SIMPLEX>
FiniteElementPolicy<dim,SIMPLEX>::FiniteElementPolicy( csmp::FiniteElement* eptr )
: fptr_(eptr)
{
}


template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::Assign( FiniteElement* fe_ptr )
 {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( fe_ptr != nullptr );
    assert( eptr->Nodes() == fe_ptr->Nodes() );
    fptr_ = fe_ptr;
 }



template<size_t dim, template<size_t> class SIMPLEX>
CSMP_FEM_TYPE FiniteElementPolicy<dim,SIMPLEX>::FE_Type() const
  {
    assert( fptr_ != nullptr );
    return fptr_->ElementType();
  }

template<size_t dim, template<size_t> class SIMPLEX>
FiniteElement* FiniteElementPolicy<dim,SIMPLEX>::FE() const
  {
    assert( fptr_ != nullptr );
    return fptr_;
  }


template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementPolicy<dim,SIMPLEX>::IsLineElement() const
  {
    assert( fptr_ != nullptr );
    return fptr_->IsLineElement();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementPolicy<dim,SIMPLEX>::IsSurfaceElement() const
  {
    assert( fptr_ != nullptr );
    return fptr_->IsSurfaceElement();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementPolicy<dim,SIMPLEX>::IsVolumeElement() const
  {
    assert( fptr_ != nullptr );
    return fptr_->IsVolumeElement();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t  FiniteElementPolicy<dim,SIMPLEX>::Interpolation() const
  {
    assert( fptr_ != nullptr );
    return fptr_->Interpolation();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool   FiniteElementPolicy<dim,SIMPLEX>::UsesLocalCoordinates() const
  {
    assert( fptr_ != nullptr );
    return fptr_->UsesLocalCoordinates();
  }


template<size_t dim, template<size_t> class SIMPLEX>
size_t   FiniteElementPolicy<dim,SIMPLEX>::Segments() const
 {
    assert( fptr_ != nullptr );
    return fptr_->Segments();
 }



/// Returns number of integration points per element, 0 if no FiniteElement assigned
template<size_t dim, template<size_t> class SIMPLEX>
size_t FiniteElementPolicy<dim,SIMPLEX>::IntegrationPoints() const
  {
    if( fptr_ != nullptr ) return fptr_->IntegrationPoints();
    return 0U;
  }

/**
     Returns the position of the element integration point i in physical space.
     The coordinate value is transformed into global coordinates by the
     finite element. Therefore the coordinate matrix must be initialized.
*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementPolicy<dim,SIMPLEX>::IntegrationPoint( size_t ip ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    // DNR matrix is used here as a temporary
    // so that no extra matrix has to be created (method uses NRST internally)
    fptr_->IntegrationPoint( ip, fptr_->DNR );
    return Point<dim>( fptr_->DNR );
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::WeightAtIntegrationPoint( size_t i ) const
  {
    assert( fptr_ != nullptr );
    return fptr_->WeightAtIntegrationPoint(i);
  }


// SHAPE FUNCTIONS AT DIFFERENT POINTS

/**
    Returns the value of the interpolation functions at a point in physical space.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::N_AtGlobalPoint( std::vector<double64>& Nn, const std::vector<double64>& xyz ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->N( Nn, xyz );
  }

template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::N_AtIntegrationPoint( size_t ipoint, std::vector<double64>& Nn ) const
  {
    assert( fptr_ != nullptr );
    fptr_->N_AtIntegrationPoint( ipoint, Nn );
  }


template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::N_AtBaryCenter( std::vector<double64>& Nn ) const
  {
    assert( fptr_ != nullptr );
    if ( !fptr_->Isoparametric() ) CoordinateMatrix();
    fptr_->N_AtBaryCenter( Nn );
  }


// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS
template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim, SIMPLEX>::Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K) const
{
	assert(fptr_ != nullptr);
	CoordinateMatrix();
	fptr_->Integral_dNT_K_dN(M, K);
}

template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::dN( DenseMatrix<DM_MIN>& M ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->dN( M );
  }
  

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::dN_AtNode( DenseMatrix<DM_MIN>& M, size_t nd, size_t dof )  const
  { 
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    if ( dof == 1U )
        return fptr_->dN_AtNode( M, nd );
    double64 detJ = fptr_->dN_AtNode( M, nd );
    dof == 2 ? dN_To2DOF( fptr_->Nodes(), M ) : dN_To3DOF( fptr_->Nodes(), M );
    return detJ;
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::dN_AtBaryCenter( DenseMatrix<DM_MIN>& M, size_t dof ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    if ( dof == 1U )
        return fptr_->dN_AtBarycenter( M );
    double64 detJ = fptr_->dN_AtBarycenter( M );
    dof == 2U ? dN_To2DOF( fptr_->Nodes(), M ) : dN_To3DOF( fptr_->Nodes(), M );
    return detJ;
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M,
                                                                  size_t gp, size_t dof )  const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    if ( dof == 1U )
        return fptr_->dN_AtIntegrationPoint( M, gp );
    double64 detJ = fptr_->dN_AtIntegrationPoint( M, gp );
    dof == 2U ? dN_To2DOF( fptr_->Nodes(), M ) : dN_To3DOF( fptr_->Nodes(), M );
    return detJ;
  }



// INTEGRATION

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::det_JINV_AtIntegrationPoint( size_t ipoint ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    // SKM_FIX
    fptr_->JacobianAtIntegrationPoint( ipoint );
    return fptr_->JacobianDeterminant();
  }

template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::IntegralNN( DenseMatrix<DM_MIN>& M ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->IntegralNN( M );
  }



/**

Writes the Elements node coordinates into a matrix of rows=nodes,
columns=spatial dimensions. The coordinate matrix is internally initialized
for connected FiniteElement subclass.

@section implementation Implementation

During the visitation of an element the coordinate matrix is initialized
only once, if the method is called without argument.

@attention  for efficiency this method checks an element id in the 
connected finite element against that of the current element that 
it is called by. If they are the same, nothing is done. 
Thus, YOU MUST SET THE ELEMENT ID for this method to have an effect.

*/
template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::CoordinateMatrix() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    assert( fptr_ != nullptr );
    // if this method is called the second time on the same element it will do nothing
    if ( fptr_->CurrentID() == eptr->Idx() ) return;

    eptr->NodeCoordinateMatrix( fptr_->XY );
    
    // resetting the element id on connected finite element
    fptr_->CurrentID( eptr->Idx() );

  } // end CoordinateMatrix




// PROPERTY INTERPOLATION, EXTRAPOLATION AND INTEGRATIONS

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::PropertyIntegral( const csmp::Index& prop_key ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( fptr_ != nullptr );

    if ( prop_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "FiniteElementPolicy<dim,SIMPLEX>::PropertyIntegral",
                                            "This method only integrates scalar properties" );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,SIMPLEX>::PropertyIntegral:",
                                 "encountered element with invalid element pointer.");
       }

    // if the element uses analytical integration, a linear element is assumed and the 
    // integration is performed for a property value interpolated to the element's barycentre
    if ( !fptr_->UsesLocalCoordinates() ||
         (prop_key.place == ELEMENT || prop_key.place == FACE || prop_key.place == INTER_FACE) )
      {
        ScalarVariable  sc;
        PropertyValueAtBaryCenter( prop_key, sc );
        return Volume() * sc();
      }

    if ( prop_key.place == ELEMENT_INTEGRATION_POINT )
      {
         double64 sum(0.);
         const size_t n_integration_points(fptr_->IntegrationPoints());
         CoordinateMatrix();
         for ( size_t i=0U; i <n_integration_points; i++ )
         {
             fptr_->JacobianAtIntegrationPoint( i );
             sum += fptr_->JacobianInverse() * fptr_->WeightAtIntegrationPoint(i) * eptr->Read( i, prop_key );
         }
         return sum;
      }

    // NODE properties
    // ---------------  
    // if the element is isoparametric we loop over the integration points interpolating
    // the property values to these positions and then integrating via multiplication with 
    // weights
    double64  sumN, sumI(0.);
    const size_t n_integration_points(fptr_->IntegrationPoints());

    CoordinateMatrix();
    for ( size_t i=0U; i <n_integration_points; i++ )
      {
        sumN = 0.;
        fptr_->JacobianAtIntegrationPoint(i);
        fptr_->N_AtIntegrationPoint( i, fptr_->NRST );
        for ( size_t j=0; j<fptr_->Nodes(); j++ )
            sumN += fptr_->NRST[j] * eptr->N(j)->Read( prop_key );
        sumI += sumN * fptr_->JacobianInverse() * fptr_->WeightAtIntegrationPoint(i);
      }  

    return sumI;  

  } // end PropertyIntegral





/**

Interpolates property value to a point within the element as specified
using physical coordinates.

@section implementation Implementation

If this method is used on an element with a local coordinate system
the point in local coordinates has to be found iteratively. This is
a much slower procedure than direct interpolation.

@attention If variable is an element property, method will
simply return this value.

 */
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAt( const csmp::Index& idx,
                                                        const std::vector<double64>& xyz,
                                                        Var& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE ) {
       std::cerr <<"\nFiniteElementPolicy<"<< dim;
       std::cerr <<">::PropertyValueAt: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       throw std::domain_error("FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAt");
    }

    Var  temp;
    temp.Resize( idx.dataDepth );
    var.Resize( idx.dataDepth, 0. );
    // initialisation for accumulation
    var = 0.;

    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->N( fptr_->NRST, xyz );

    const size_t  n_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
    {
       eptr->N(i)->Read( idx, temp );
       var += temp * fptr_->NRST[i];
    }
  }



// scalar version of previous method
template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAt( const csmp::Index& idx,
                                                            const std::vector<double64>& xyz ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         return eptr->Read( idx );
    }
    if ( idx.place != NODE ) {
       std::cerr <<"\nFiniteElementPolicy<"<< dim;
       std::cerr <<">::PropertyValueAt: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       throw std::domain_error("FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAt");
    }

    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->N( fptr_->NRST, xyz );

    double64 var(0.);
    const size_t  n_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
       var += eptr->N(i)->Read( idx ) * fptr_->NRST[i];

    return var;
}


/**

Node variables are interpolated to barycenter.
Integration point variables (ELEMENT_INTEGRATION_POINT) are averaged.

@attention If variable is an element property, method will
simply return this value.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtBaryCenter( const csmp::Index& idx,
                                                                  Var& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( fptr_ != nullptr );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE  and idx.place != ELEMENT_INTEGRATION_POINT and idx.place != SECTOR_INTEGRATION_POINT) {
       std::cerr <<"\nFiniteElementPolicy<"<< dim;
       std::cerr <<">::PropertyValueAtBaryCenter: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       std::cerr <<"values of ELEMENT_INTEGRATION_POINT properties are averaged."<< std::endl;
       throw std::domain_error("FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtBaryCenter");
    }

    Var  temp;
    temp.Resize( idx.dataDepth );
    var.Resize( idx.dataDepth, 0. );
    // initialisation for accumulation
    var = 0.;

    // simple averaging of integration point properties
    if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
       const size_t n_integration_points(IntegrationPoints());
       for ( size_t i=0U; i < n_integration_points; i++ ) {
            eptr->Read( i, idx, temp );
            var += temp;
         }
       var /= static_cast<double64>(fptr_->IntegrationPoints());
       return;
    }
    if ( idx.place == SECTOR_INTEGRATION_POINT ) {
       const size_t n_sector_integration_points(eptr->FV()->Sectors());
       for ( size_t i=0U; i < n_sector_integration_points; i++ ) {
            eptr->Read( i, 0U, idx, temp );
            var += temp;
         }
       var /= static_cast<double64>(n_sector_integration_points);
       return;
    }

    if ( !fptr_->UsesLocalCoordinates() ) CoordinateMatrix();
    fptr_->N_AtBaryCenter( fptr_->NRST );

    const size_t  n_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<n_nodes; i++ ) {
       eptr->N(i)->Read( idx, temp );
       var += temp * fptr_->NRST[i];
    }
  }




template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtBaryCenter( const csmp::Index& idx ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( fptr_ != nullptr );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE )  return eptr->Read( idx );

    if ( idx.place != NODE  and idx.place != ELEMENT_INTEGRATION_POINT and idx.place != SECTOR_INTEGRATION_POINT) {
       std::cerr <<"\nFiniteElementPolicy<"<< dim;
       std::cerr <<">::PropertyValueAtBaryCenter: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       std::cerr <<"values of ELEMENT_INTEGRATION_POINT properties are averaged."<< std::endl;
       throw std::domain_error("FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtBaryCenter");
    }

    double64 var(0.);

    // simple averaging of integration point properties
    if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
       const size_t n_integration_points(IntegrationPoints());
       for ( size_t i=0U; i < n_integration_points; i++ ) {
            var += eptr->Read( i, idx );
         }
       var /= static_cast<double64>(fptr_->IntegrationPoints());
       return var;
    }
    if ( idx.place == SECTOR_INTEGRATION_POINT ) {
       const size_t n_sector_integration_points(eptr->FV()->Sectors());
       for ( size_t i=0U; i < n_sector_integration_points; i++ ) {
            var += eptr->Read( i, 0U, idx );
         }
       var /= static_cast<double64>(n_sector_integration_points);
       return var;
    }

    if ( !fptr_->UsesLocalCoordinates() ) CoordinateMatrix();
    fptr_->N_AtBaryCenter( fptr_->NRST );

    const size_t  n_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<n_nodes; i++ ) {
       var += eptr->N(i)->Read( idx ) * fptr_->NRST[i];
    }
    
    return var;
    
  } // end scalar version





/**
Interpolates node properties to the integration points.

@attention If variable is an element property, method will
simply return this value.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtIntegrationPoint( const csmp::Index& idx,
                                                                        size_t ip,
                                                                        Var& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE ) {
         std::cerr <<"\nFiniteElementPolicy<"<< dim;
         std::cerr <<">::PropertyValueAtIntegrationPoint: This method only interpolates the ";
         std::cerr <<"values of NODE properties."<< std::endl;
         throw std::domain_error("FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtIntegrationPoint");
      }

    // potential size and value adjustments for array variables
    Var  temp;
    temp.Resize( idx.dataDepth );
    var.Resize( idx.dataDepth, 0. );
    // initialisation for accumulation
    var = 0.;

    assert( fptr_ != nullptr );
    fptr_->N_AtIntegrationPoint( ip, fptr_->NRST );

    const size_t  n_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
    {
       eptr->N(i)->Read( idx, temp );
       var += temp * fptr_->NRST[i];
    }
  }



/**
Interpolates SCALAR node properties to the integration points.

@attention If variable is an element property, method will
simply return this value.
*/
template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtIntegrationPoint( const csmp::Index& idx,
                                                                            size_t ip ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( (idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE) && idx.type == SCALAR ) return eptr->Read( idx );
    
    if ( idx.place != NODE or idx.type != SCALAR ) {
        std::cerr <<"\nFiniteElementPolicy<"<< dim;
        std::cerr <<">::PropertyValueAtIntegrationPoint: This method only interpolates ";
        std::cerr <<"scalar NODE properties."<< std::endl;
        throw std::domain_error("FiniteElementPolicy<dim,SIMPLEX>::PropertyValueAtIntegrationPoint");
     }

    assert( fptr_ != nullptr );
    fptr_->N_AtIntegrationPoint( ip, fptr_->NRST );

    double64  var(0.);
    const size_t  n_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
      var += fptr_->NRST[i] * eptr->N(i)->Read( idx );

   return var;
  }

/**

Returns the ELEMENT_INTEGRATION_POINT variable values into the parameter
vector.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void  FiniteElementPolicy<dim,SIMPLEX>::IntegrationPointPropertyVector( const csmp::Index& idx,
                                                                        std::vector<Var>& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,SIMPLEX>::IntegrationPointPropertyVector:",
                                 "encountered element with invalid element pointer.");
       }

    assert( fptr_->UsesLocalCoordinates() );
    assert( fptr_->IntegrationPoints() > 0U );

    // resizing V if necessary
    const size_t n_integration_points(fptr_->IntegrationPoints());
    var.resize( n_integration_points );

    for ( size_t i=0U; i<n_integration_points; i++ )
    {
        var[i].Resize( idx.dataDepth );
        var[i] = 0.;
        eptr->Read( i, idx, var[i] );
    }
  }


/**

    Linear extrapolation of a variable placed on the element's integration points
    to its nodes.

*/
template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>
::ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                              const std::vector<double64>& IVAR,
                                              std::vector<double64>& NVAR ) const
  {
    if( fptr_ == nullptr )
        throw csmp::Exception( ERROR,
                               "FiniteElementPolicy<dim,SIMPLEX>::ExtrapolateIntegrationPointVariableToNodes:",
                               "encountered element with invalid element pointer.");
    if( !fptr_->UsesLocalCoordinates() )
        throw csmp::Exception( ERROR,
                               "FiniteElementPolicy<dim,SIMPLEX>::ExtrapolateIntegrationPointVariableToNodes:",
                               "Extrapolation is not supported; only the isoparametric elements have integration points.");

    // NB: the following method checks that the supplied arrays have the right size
    // NB: NO coordinate matrix is required here
    fptr_->ExtrapolateIntegrationPointVariableToNodes( nvars, IVAR, NVAR );
  }



// ELEMENT GEOMETRY
template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::AspectRatio() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->AspectRatio();
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::InnerRadius() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->InnerRadius();
  }


template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::Volume() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->Volume();
  }





// SEGMENTS


/**
    Computes the length of the finite element edges = segments.
    
    @attention this method assumes that the first two nodes returned by
    NodesOfSegment() are the corner nodes, i.e. do not include midside nodes.
*/
template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementPolicy<dim,SIMPLEX>::SegmentLength( size_t segm ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,SIMPLEX>::SegmentLenghts:",
                                 "encountered element with invalid element pointer.");
       }

    std::vector<size_t>  snids;
    fptr_->NodesOfSegment( segm, snids );

    assert( snids.size() == 2U );
    return eptr->N(snids[1])->Coordinate().DistanceTo( eptr->N(snids[0])->Coordinate() );

  } // end BoundarySegmentLength



/// returns the length of all element segments=edges
template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementPolicy<dim,SIMPLEX>::SegmentLengths( std::vector<double64>& lengths ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->EdgeLengths( lengths );
    
  } // end SegmentLengths


/**
    returns the mid-point (in physical space) of the finite element segment = edge.
*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementPolicy<dim,SIMPLEX>::SegmentMidPoint( size_t segm ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( fptr_ != nullptr );
    fptr_->NodesOfSegment( segm, fptr_->IDX );
    assert( fptr_->IDX.size() == 2U );
    return midPoint(eptr->N(fptr_->IDX[0])->Coordinate(),eptr->N(fptr_->IDX[1])->Coordinate());
  }






// FACES


/**

Method calculates the area of the finite element face
using its corner nodes and returns it.

@attention works only for straight-sided elements.

@attention for lower-dimensional elements extra attributes are
needed to scale the returned parameter to arrive at the
correct spatial dimension.

@attention this method fails to detect element type for cubic or barycentric elements

@author SKM 4/6/2014

*/
template<size_t dim, template<size_t> class SIMPLEX>
double64  FiniteElementPolicy<dim,SIMPLEX>::FaceArea( size_t n ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,SIMPLEX>::FaceArea:",
                                 "encountered element with invalid element pointer.");
       }

    if ( fptr_->IsVolumeElement() ) {
         std::vector<size_t>  ids;
         fptr_->NodesOfFace( n, ids );

         // if the face is a triangle, function computes unit normal
         // tested: O.K.
         if ( ids.size() == 3U or ids.size() == 6U )
              return triangleArea( eptr->N(ids[0])->Coordinate(),
                                   eptr->N(ids[1])->Coordinate(),
                                   eptr->N(ids[2])->Coordinate() );

         // if the face is a quadrilateral unit normal is obtained this way
         // tested: O.K.
         if ( ids.size() == 4U or ids.size() == 8U or ids.size() == 9U )
              return facetArea4( eptr->N(ids[0])->Coordinate(),
                                 eptr->N(ids[1])->Coordinate(),
                                 eptr->N(ids[2])->Coordinate(),
                                 eptr->N(ids[3])->Coordinate() );
      }

    // if the element is a surface, the face is a line element
    // and its length will be returned
    if ( fptr_->IsSurfaceElement() ) {
         std::vector<size_t>  ids;
         fptr_->NodesOfFace( n, ids );
         assert( ids.size() == 2 );
         Point<dim>  edge = eptr->N(ids[1])->Coordinate() - eptr->N(ids[0])->Coordinate();
         return edge.Length();
      }

    // for a line element the area of a face is 1x1 m
    return 1.;

  } // end FaceArea





/**
    SKM new implementation
*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementPolicy<dim,SIMPLEX>::FaceBaryCenter( size_t face ) const
  {
    assert( fptr_ != nullptr );
    const SIMPLEX<dim>* e( static_cast<const SIMPLEX<dim>*>(this) );
    if ( fptr_ == nullptr ) {
          e->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,SIMPLEX>::FaceBaryCenter:",
                                "encountered element with invalid element pointer.");
       }
    fptr_->NodesOfFace( face, fptr_->IDX );
    const size_t face_nodes(fptr_->IDX.size());
    Point<dim>  barycenter(e->N(fptr_->IDX[0])->Coordinate());
    for ( size_t i=1U; i<face_nodes; ++i )
      barycenter += e->N(fptr_->IDX[i])->Coordinate();
    barycenter /= static_cast<double64>(face_nodes);

    return barycenter;
  }






// UNIT NORMALS



/**
    Computes outward pointing normal and returns the result into supplied vector variable.
    
    @note not as efficient as version that uses standard vector.
 
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementPolicy<dim,SIMPLEX>::UnitNormal( VectorVariable<dim>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    std::vector<double64> un( dim, 0. );
    fptr_->UnitNormal(un);
    for ( size_t d(0); d < dim; ++d )
        nrml(d) = un[d];
  }


template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementPolicy<dim,SIMPLEX>::UnitNormal( std::vector<double64>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->UnitNormal(nrml);
  }


template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementPolicy<dim,SIMPLEX>::UnitNormal() const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    std::vector<double64> un(dim,0.);
    fptr_->UnitNormal(un);
    return Point<dim>(un);
  }


/**

Method calculates the unit normal that is outward pointing to to the queried face.

@test SKM 22/9/2014: fixed 3D normals so that they are outward pointing.

*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementPolicy<dim,SIMPLEX>::UnitNormalToFace( size_t face ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    std::vector<double64> un(dim,0.);
    fptr_->UnitNormalToFace( face, un );
    return Point<dim>(un);

  } // end UnitNormalToFace



/** 
    Fastest version to compute outward-pointing normal to target face.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementPolicy<dim,SIMPLEX>::UnitNormalToFace( size_t face, std::vector<double64>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->UnitNormalToFace( face, nrml );
  }



template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementPolicy<dim,SIMPLEX>::UnitNormalToFace( size_t face, VectorVariable<dim>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    std::vector<double64> un(dim,0.);
    fptr_->UnitNormalToFace( face, un );
    for ( size_t d(0); d < dim; ++d )
      nrml(d) = un[d];
  }



// ELEMENT

// vector of properties from integration points
// scalar
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, size_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, size_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, size_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, size_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, size_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, size_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;


// FACE

// vector of properties from integration points
// scalar
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, size_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, size_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, size_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, size_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, size_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, size_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;

// INTERFACE

// vector of properties from integration points
// scalar
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, size_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, size_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, size_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, size_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, size_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, size_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;





template class FiniteElementPolicy<1U,Element>;
template class FiniteElementPolicy<1U,Face>;
template class FiniteElementPolicy<1U,InterFace>;

template class FiniteElementPolicy<2U,Element>;
template class FiniteElementPolicy<2U,Face>;
template class FiniteElementPolicy<2U,InterFace>;

template class FiniteElementPolicy<3U,Element>;
template class FiniteElementPolicy<3U,Face>;
template class FiniteElementPolicy<3U,InterFace>;

} // end csmp
