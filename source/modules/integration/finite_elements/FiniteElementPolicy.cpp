#include "FiniteElementPolicy.h"
#include "FiniteElement.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Edge.h"
#include "variableOperations.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include "Point.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
FiniteElementPolicy<dim,CELL>::FiniteElementPolicy( csmp::FiniteElement* eptr )
: fptr_(eptr)
{
}


template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::Assign( FiniteElement* fe_ptr )
 {
    assert( fe_ptr != nullptr );
    fptr_ = fe_ptr;
 }



template<uint32_t dim, template<uint32_t> class CELL>
CSMP_FEM_TYPE FiniteElementPolicy<dim,CELL>::FE_Type() const
  {
    assert( fptr_ != nullptr );
    return fptr_->ElementType();
  }

template<uint32_t dim, template<uint32_t> class CELL>
FiniteElement* FiniteElementPolicy<dim,CELL>::FE() const
  {
    return fptr_;
  }


template<uint32_t dim, template<uint32_t> class CELL>
bool FiniteElementPolicy<dim,CELL>::IsEquidimensional() const
 {
    if constexpr ( dim == 3U ) return IsVolume();
    if constexpr ( dim == 2U ) return IsSurface();
    if constexpr ( dim == 1U ) return IsLine();
    return false;
 }




template<uint32_t dim, template<uint32_t> class CELL>
bool FiniteElementPolicy<dim,CELL>::IsLine() const
  {
    if constexpr ( dim == 1U ) return true;
    assert( fptr_ != nullptr );
    return fptr_->IsLine();
  }

template<uint32_t dim, template<uint32_t> class CELL>
bool FiniteElementPolicy<dim,CELL>::IsSurface() const
  {
    assert( fptr_ != nullptr );
    return fptr_->IsSurface();
  }

template<uint32_t dim, template<uint32_t> class CELL>
bool FiniteElementPolicy<dim,CELL>::IsVolume() const
  {
    if constexpr ( dim != 3U ) return false;
    assert( fptr_ != nullptr );
    return fptr_->IsVolume();
  }

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  FiniteElementPolicy<dim,CELL>::Interpolation() const
  {
    assert( fptr_ != nullptr );
    return fptr_->Interpolation();
  }

template<uint32_t dim, template<uint32_t> class CELL>
bool   FiniteElementPolicy<dim,CELL>::UsesLocalCoordinates() const
  {
    assert( fptr_ != nullptr );
    return fptr_->UsesLocalCoordinates();
  }


template<uint32_t dim, template<uint32_t> class CELL>
uint32_t   FiniteElementPolicy<dim,CELL>::Segments() const
 {
    assert( fptr_ != nullptr );
    return fptr_->Segments();
 }



/// Returns number of integration points per element, 0 if no FiniteElement assigned
template<uint32_t dim, template<uint32_t> class CELL>
uint32_t FiniteElementPolicy<dim,CELL>::IntegrationPoints() const
  {
    if( fptr_ != nullptr ) return fptr_->IntegrationPoints();
    return 0U;
  }

/**
     Returns the position of the element integration point i in physical space.
     The coordinate value is transformed into global coordinates by the
     finite element. Therefore the coordinate matrix must be initialized.
*/
template<uint32_t dim, template<uint32_t> class CELL>
Point<dim>  FiniteElementPolicy<dim,CELL>::IntegrationPoint( uint32_t ip ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    // DNR matrix is used here as a temporary
    // so that no extra matrix has to be created (method uses NRST internally)
    fptr_->IntegrationPoint( ip, fptr_->DNR );
    return Point<dim>( fptr_->DNR );
  }


template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::WeightAtIntegrationPoint( uint32_t i ) const
  {
    assert( fptr_ != nullptr );
    return fptr_->WeightAtIntegrationPoint(i);
  }


// MAPPING BETWEEN LOCAL AND GLOBAL COORDINATES

/**
    Transforms point location from local to global coordinates.
 
    @note passing the point by value is intentional because the copy is reused by the function.
*/
template<uint32_t dim, template<uint32_t> class CELL>
Point<dim>  FiniteElementPolicy<dim,CELL>::RstToXYZ( Point<dim> rst ) const
{
   const CELL<dim>* e( static_cast<const CELL<dim>*>(this) );

   // initialises NRST
   N_At( rst );
  
   assert( e  != nullptr );
   rst = 0.;
   const auto nodes(e->Nodes());
   for ( auto i{0U}; i<nodes; i++ )
     for ( auto j{0U}; j<dim; ++j )
       // transformation of the coordinates
       rst[j] += e->FE()->NRST[i] * (*e->N(i))[j];
  
   // standard RVO
   return rst;
}




// SHAPE FUNCTIONS AT DIFFERENT POINTS



// SHAPE FUNCTIONS AT DIFFERENT POINTS


/**
     Calculates local interpolation function values at the supplied parametric
     coordinate->
*/
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::N_At( const Point<dim>& rst ) const
 {
    const CELL<dim>* e( static_cast<const CELL<dim>*>(this) );
    assert( e != nullptr );
 
    if ( e->IsVolume() ) {
         e->FE()->Nrst( rst[0], rst[1], rst[2], e->FE()->NRST );
         return;
      }
    if ( e->IsSurface() ) {
         e->FE()->Nrs( rst[0], rst[1], e->FE()->NRST );
         return;
      }
    // line element
    e->FE()->Nr( rst[0], e->FE()->NRST );
 }



/**
    Returns element interpolation functions at point rst in local coordinate system
    into argument vector IPOL.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::N_At( const Point<dim>& rst,
                                          vector<double>& IPOL ) const
 {

    const CELL<dim>* e( static_cast<const CELL<dim>*>(this) );
    assert( e != nullptr );
    if ( e->IsVolume() ) {
         e->FE()->Nrst( rst[0], rst[1], rst[2], IPOL );
         return;
      }
    if ( e->FE()->IsSurface() ) {
         e->FE()->Nrs( rst[0], rst[1], IPOL );
         return;
      }
    // line element
    e->FE()->Nr( rst[0], IPOL );
 }



/**
    Returns the value of the interpolation functions at a point in physical space.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::N_AtGlobalPoint( vector<double>& Nn, const vector<double>& xyz ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->N( Nn, xyz );
  }


template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::N_AtIntegrationPoint( uint32_t ipoint, vector<double>& Nn ) const
  {
    assert( fptr_ != nullptr );
    fptr_->N_AtIntegrationPoint( ipoint, Nn );
  }


template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::N_AtBaryCenter( vector<double>& Nn ) const
  {
    assert( fptr_ != nullptr );
    if ( !fptr_->Isoparametric() ) CoordinateMatrix();
    fptr_->N_AtBaryCenter( Nn );
  }


// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim, CELL>::Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K) const
{
	assert(fptr_ != nullptr);
	CoordinateMatrix();
	fptr_->Integral_dNT_K_dN(M, K);
}


template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::dN( DenseMatrix<DM_MIN>& M ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->dN( M );
  }
  

template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t nd, uint32_t dof )  const
  { 
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    if ( dof == 1U )
        return fptr_->dN_AtNode( M, nd );
    double detJ = fptr_->dN_AtNode( M, nd );
    dof == 2 ? dN_To2DOF( fptr_->Nodes(), M ) : dN_To3DOF( fptr_->Nodes(), M );
    return detJ;
  }

template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::dN_AtBaryCenter( DenseMatrix<DM_MIN>& M, uint32_t dof ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    if ( dof == 1U )
        return fptr_->dN_AtBarycenter( M );
    double detJ = fptr_->dN_AtBarycenter( M );
    dof == 2U ? dN_To2DOF( fptr_->Nodes(), M ) : dN_To3DOF( fptr_->Nodes(), M );
    return detJ;
  }

template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M,
                                                                  uint32_t gp, uint32_t dof )  const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    if ( dof == 1U )
        return fptr_->dN_AtIntegrationPoint( M, gp );
    double detJ = fptr_->dN_AtIntegrationPoint( M, gp );
    dof == 2U ? dN_To2DOF( fptr_->Nodes(), M ) : dN_To3DOF( fptr_->Nodes(), M );
    return detJ;
  }



// INTEGRATION

template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::det_JINV_AtIntegrationPoint( uint32_t ipoint ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    // SKM_FIX
    fptr_->JacobianAtIntegrationPoint( ipoint );
    return fptr_->JacobianDeterminant();
  }

template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::IntegralNN( DenseMatrix<DM_MIN>& M ) const
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
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::CoordinateMatrix() const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    assert( fptr_ != nullptr );
    // if this method is called the second time on the same element it will do nothing
    if ( fptr_->CurrentID() == eptr->Idx() ) return;

    eptr->NodeCoordinateMatrix( fptr_->XY );
    
    // resetting the element id on connected finite element
    fptr_->CurrentID( eptr->Idx() );

  } // end CoordinateMatrix




// PROPERTY INTERPOLATION, EXTRAPOLATION AND INTEGRATIONS

template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::PropertyIntegral( const csmp::Index& prop_key ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
    assert( fptr_ != nullptr );

    if ( prop_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "FiniteElementPolicy<dim,CELL>::PropertyIntegral",
                                            "This method only integrates scalar properties" );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,CELL>::PropertyIntegral:",
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
         double sum(0.);
         const auto n_integration_points(fptr_->IntegrationPoints());
         CoordinateMatrix();
         for ( auto i{0U}; i <n_integration_points; i++ )
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
    double  sumN, sumI(0.);
    const auto n_integration_points(fptr_->IntegrationPoints());

    CoordinateMatrix();
    for ( auto i{0U}; i <n_integration_points; i++ )
      {
        sumN = 0.;
        fptr_->JacobianAtIntegrationPoint(i);
        fptr_->N_AtIntegrationPoint( i, fptr_->NRST );
        for ( auto j{0U}; j<fptr_->Nodes(); j++ )
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
template<uint32_t dim, template<uint32_t> class CELL>
template<class Var>
void FiniteElementPolicy<dim,CELL>::PropertyValueAt( const csmp::Index& idx,
                                                        const vector<double>& xyz,
                                                        Var& var ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE ) {
       cerr <<"\nFiniteElementPolicy<"<< dim;
       cerr <<">::PropertyValueAt: This method only interpolates the ";
       cerr <<"values of NODE properties."<< endl;
       throw domain_error("FiniteElementPolicy<dim,CELL>::PropertyValueAt");
    }

    Var  temp;
    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( idx.dataDepth );
        var.Resize( idx.dataDepth, 0. );
      }
    // initialisation for accumulation
    var = 0.;

    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->N( fptr_->NRST, xyz );

    const auto n_nodes(fptr_->Nodes());
    for ( auto i{0U}; i<n_nodes; i++ )
    {
       eptr->N(i)->Read( idx, temp );
       var += temp * fptr_->NRST[i];
    }
  }



// scalar version of previous method
template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::PropertyValueAt( const csmp::Index& idx,
                                                            const vector<double>& xyz ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         return eptr->Read( idx );
    }
    if ( idx.place != NODE ) {
       cerr <<"\nFiniteElementPolicy<"<< dim;
       cerr <<">::PropertyValueAt: This method only interpolates the ";
       cerr <<"values of NODE properties."<< endl;
       throw domain_error("FiniteElementPolicy<dim,CELL>::PropertyValueAt");
    }

    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->N( fptr_->NRST, xyz );

    double var(0.);
    const auto  n_nodes(fptr_->Nodes());
    for ( auto i{0U}; i<n_nodes; i++ )
       var += eptr->N(i)->Read( idx ) * fptr_->NRST[i];

    return var;
}


/**

Node variables are interpolated to barycenter.
Integration point variables (ELEMENT_INTEGRATION_POINT) are averaged.

@attention If variable is an element property, method will
simply return this value.

*/
template<uint32_t dim, template<uint32_t> class CELL>
template<class Var>
void FiniteElementPolicy<dim,CELL>::PropertyValueAtBaryCenter( const csmp::Index& idx,
                                                               Var& var ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
    assert( fptr_ != nullptr );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE  and idx.place != ELEMENT_INTEGRATION_POINT and idx.place != SECTOR_INTEGRATION_POINT) {
       cerr <<"\nFiniteElementPolicy<"<< dim;
       cerr <<">::PropertyValueAtBaryCenter: This method only interpolates the ";
       cerr <<"values of NODE properties."<< endl;
       cerr <<"values of ELEMENT_INTEGRATION_POINT properties are averaged."<< endl;
       throw domain_error("FiniteElementPolicy<dim,CELL>::PropertyValueAtBaryCenter");
    }

    Var  temp;
    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( idx.dataDepth );
        var.Resize( idx.dataDepth, 0. );
      }
    // initialisation for accumulation
    var = 0.;

    // simple averaging of integration point properties
    if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
       const auto n_integration_points(IntegrationPoints());
       for ( auto i{0U}; i < n_integration_points; i++ ) {
            eptr->Read( i, idx, temp );
            var += temp;
         }
       var /= static_cast<double>(fptr_->IntegrationPoints());
       return;
    }
    if ( idx.place == SECTOR_INTEGRATION_POINT ) {
       const auto n_sector_integration_points(eptr->FV()->Sectors());
       for ( auto i{0U}; i < n_sector_integration_points; i++ ) {
            eptr->Read( i, 0U, idx, temp );
            var += temp;
         }
       var /= static_cast<double>(n_sector_integration_points);
       return;
    }

    if ( !fptr_->UsesLocalCoordinates() ) CoordinateMatrix();
    fptr_->N_AtBaryCenter( fptr_->NRST );

    const auto  n_nodes(fptr_->Nodes());
    for ( auto i{0U}; i<n_nodes; i++ ) {
       eptr->N(i)->Read( idx, temp );
       var += temp * fptr_->NRST[i];
    }
  }




template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::PropertyValueAtBaryCenter( const csmp::Index& idx ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
    assert( fptr_ != nullptr );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE )  return eptr->Read( idx );

    if ( idx.place != NODE  and idx.place != ELEMENT_INTEGRATION_POINT and idx.place != SECTOR_INTEGRATION_POINT) {
       cerr <<"\nFiniteElementPolicy<"<< dim;
       cerr <<">::PropertyValueAtBaryCenter: This method only interpolates the ";
       cerr <<"values of NODE properties."<< endl;
       cerr <<"values of ELEMENT_INTEGRATION_POINT properties are averaged."<< endl;
       throw domain_error("FiniteElementPolicy<dim,CELL>::PropertyValueAtBaryCenter");
    }

    double var(0.);

    // simple averaging of integration point properties
    if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
       const auto n_integration_points(IntegrationPoints());
       for ( auto i{0U}; i < n_integration_points; i++ ) {
            var += eptr->Read( i, idx );
         }
       var /= static_cast<double>(fptr_->IntegrationPoints());
       return var;
    }
    if ( idx.place == SECTOR_INTEGRATION_POINT ) {
       const auto n_sector_integration_points(eptr->FV()->Sectors());
       for ( auto i{0U}; i < n_sector_integration_points; i++ ) {
            var += eptr->Read( i, 0U, idx );
         }
       var /= static_cast<double>(n_sector_integration_points);
       return var;
    }

    if ( !fptr_->UsesLocalCoordinates() ) CoordinateMatrix();
    fptr_->N_AtBaryCenter( fptr_->NRST );

    const auto  n_nodes(fptr_->Nodes());
    for ( auto i{0U}; i<n_nodes; i++ ) {
       var += eptr->N(i)->Read( idx ) * fptr_->NRST[i];
    }
    
    return var;
    
  } // end scalar version





/**
Interpolates node properties to the integration points.

@attention If variable is an element property, method will
simply return this value.

*/
template<uint32_t dim, template<uint32_t> class CELL>
template<class Var>
void FiniteElementPolicy<dim,CELL>::PropertyValueAtIntegrationPoint( const csmp::Index& idx,
                                                                     uint32_t ip,
                                                                     Var& var ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    if ( idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE ) {
         cerr <<"\nFiniteElementPolicy<"<< dim;
         cerr <<">::PropertyValueAtIntegrationPoint: This method only interpolates the ";
         cerr <<"values of NODE properties."<< endl;
         throw domain_error("FiniteElementPolicy<dim,CELL>::PropertyValueAtIntegrationPoint");
      }

    // potential size and value adjustments for array variables
    // TODO: temp variable is needed only if the interpolated variable is not a scalar; fix!
    Var  temp;
    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( idx.dataDepth );
        var.Resize( idx.dataDepth, 0. );
      }
    // initialisation for accumulation
    var = 0.;

    assert( fptr_ != nullptr );
    fptr_->N_AtIntegrationPoint( ip, fptr_->NRST );

    const auto  n_nodes(fptr_->Nodes());
    for ( auto i{0U}; i<n_nodes; i++ )
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
template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::PropertyValueAtIntegrationPoint( const csmp::Index& idx,
                                                                       uint32_t ip ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    if ( (idx.place == ELEMENT or idx.place == FACE or idx.place == INTER_FACE) && idx.type == SCALAR ) return eptr->Read( idx );
    
    if ( idx.place != NODE or idx.type != SCALAR ) {
        cerr <<"\nFiniteElementPolicy<"<< dim;
        cerr <<">::PropertyValueAtIntegrationPoint: This method only interpolates ";
        cerr <<"scalar NODE properties."<< endl;
        throw domain_error("FiniteElementPolicy<dim,CELL>::PropertyValueAtIntegrationPoint");
     }

    assert( fptr_ != nullptr );
    fptr_->N_AtIntegrationPoint( ip, fptr_->NRST );

    double  var(0.);
    const auto  n_nodes(fptr_->Nodes());
    for ( auto i{0U}; i<n_nodes; i++ )
      var += fptr_->NRST[i] * eptr->N(i)->Read( idx );

   return var;
  }

/**

Returns the ELEMENT_INTEGRATION_POINT variable values into the parameter
vector.

*/
template<uint32_t dim, template<uint32_t> class CELL>
template<class Var>
void  FiniteElementPolicy<dim,CELL>::IntegrationPointPropertyVector( const csmp::Index& idx,
                                                                     vector<Var>& var ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,CELL>::IntegrationPointPropertyVector:",
                                 "encountered element with invalid element pointer.");
       }

    assert( fptr_->UsesLocalCoordinates() );
    assert( fptr_->IntegrationPoints() > 0U );

    // resizing V if necessary
    const auto n_integration_points(fptr_->IntegrationPoints());
    var.resize( n_integration_points );

    for ( auto i{0U}; i<n_integration_points; i++ )
      {
          if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
              var[i].Resize( idx.dataDepth );
            }
          var[i] = 0.;
          eptr->Read( i, idx, var[i] );
      }
  }


/**

    Linear extrapolation of a variable placed on the element's integration points
    to its nodes.

*/
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>
::ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                              const vector<double>& IVAR,
                                              vector<double>& NVAR ) const
  {
    if( fptr_ == nullptr )
        throw csmp::Exception( ERROR,
                               "FiniteElementPolicy<dim,CELL>::ExtrapolateIntegrationPointVariableToNodes:",
                               "encountered element with invalid element pointer.");
    if( !fptr_->UsesLocalCoordinates() )
        throw csmp::Exception( ERROR,
                               "FiniteElementPolicy<dim,CELL>::ExtrapolateIntegrationPointVariableToNodes:",
                               "Extrapolation is not supported; only the isoparametric elements have integration points.");

    // NB: the following method checks that the supplied arrays have the right size
    // NB: NO coordinate matrix is required here
    fptr_->ExtrapolateIntegrationPointVariableToNodes( nvars, IVAR, NVAR );
  }



    /// returns search key to match element faces
template<uint32_t dim, template<uint32_t> class CELL>
set<Node<dim>*>  FiniteElementPolicy<dim,CELL>::CornerNodesOfFace( uint32_t face_id ) const
 {
    assert( fptr_ != nullptr );
    assert( face_id < fptr_->Faces() );
 
    const CELL<dim>* const eptr( static_cast<const CELL<dim>* const>(this) );

    set<Node<dim>*>  temp;
    for ( const auto& nit : fptr_->CornerNodesOfFace(face_id) ) {
         assert( eptr->N(nit) != nullptr );
         temp.insert( eptr->N(nit) );
      }
      
    return temp;
 }




/** returns search key to match element faces; pointers in order so that they can be searched
 */
template<uint32_t dim, template<uint32_t> class CELL>
set<Node<dim>*> FiniteElementPolicy<dim,CELL>::CornerNodesConnectedTo( uint32_t node_id ) const
 {
    assert( fptr_ != nullptr );
    assert( node_id < fptr_->CornerNodes() );
 
    const CELL<dim>* const eptr( static_cast<const CELL<dim>* const>(this) );

    set<Node<dim>*>  temp;
    for ( const auto& nit : fptr_->NodesConnectedTo(node_id) ) {
         assert( eptr->N(nit) != nullptr );
         temp.insert( eptr->N(nit) );
      }
      
    return temp;
 }



// ELEMENT GEOMETRY
template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::AspectRatio() const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->AspectRatio();
  }

template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::InnerRadius() const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->InnerRadius();
  }


template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::Volume() const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
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
template<uint32_t dim, template<uint32_t> class CELL>
double FiniteElementPolicy<dim,CELL>::SegmentLength( uint32_t segm ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );

    if ( fptr_ == nullptr ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,CELL>::SegmentLenghts:",
                                 "encountered element with invalid element pointer.");
       }

    vector<uint32_t>  snids;
    fptr_->NodesOfSegment( segm, snids );

    assert( snids.size() == 2U );
    return eptr->N(snids[1])->Coordinate().DistanceTo( eptr->N(snids[0])->Coordinate() );

  } // end BoundarySegmentLength



/// returns the length of all element segments=edges
template<uint32_t dim, template<uint32_t> class CELL>
void FiniteElementPolicy<dim,CELL>::SegmentLengths( vector<double>& lengths ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->EdgeLengths( lengths );
    
  } // end SegmentLengths


/**
    returns the mid-point (in physical space) of the finite element segment = edge.
*/
template<uint32_t dim, template<uint32_t> class CELL>
Point<dim>  FiniteElementPolicy<dim,CELL>::SegmentMidPoint( uint32_t segm ) const
  {
    const CELL<dim>* eptr( static_cast<const CELL<dim>*>(this) );
    assert( fptr_ != nullptr );
    fptr_->NodesOfSegment( segm, fptr_->IDX );
    assert( fptr_->IDX.size() == 2U );
    return midPoint(eptr->N(fptr_->IDX[0])->Coordinate(),eptr->N(fptr_->IDX[1])->Coordinate());
  }






// FACES

// STUB FOR 1D normalOfTriangle
inline double triangleArea( const Point<1U>&, const Point<1U>&, const Point<1U>& ) {
     return 1.;
  }





/**
    Finds the centre of gravity of the face, using the coordinates of all of its nodes.
*/
template<uint32_t dim, template<uint32_t> class CELL>
Point<dim>  FiniteElementPolicy<dim,CELL>::FaceBaryCenter( uint32_t face ) const
  {
    assert( fptr_ != nullptr );
    const CELL<dim>* e( static_cast<const CELL<dim>*>(this) );
    if ( fptr_ == nullptr ) {
          e->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementPolicy<dim,CELL>::FaceBaryCenter:",
                                "encountered element with invalid element pointer.");
       }
    // not needed as the nodes are accessed directly: CoordinateMatrix();
    Point<dim>  barycenter; // initialises to zero
    assert( face < fptr_->Faces() );
    for ( const auto& i : fptr_->NodesOfFace(face) )
      barycenter += e->N(i)->Coordinate();
    barycenter /= static_cast<double>(fptr_->NodesPerFace(face));

    return barycenter;
  }






// UNIT NORMALS



/**
    Computes outward pointing normal and returns the result into supplied vector variable.
    
    @note not as efficient as version that uses standard vector.
 
*/
template<uint32_t dim, template<uint32_t> class CELL>
void  FiniteElementPolicy<dim,CELL>::UnitNormal( VectorVariable<dim>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    uint32_t count{0U};
    for ( const auto& d : fptr_->UnitNormal() )
      nrml( count++ ) = d;
  }


template<uint32_t dim, template<uint32_t> class CELL>
void  FiniteElementPolicy<dim,CELL>::UnitNormal( vector<double>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    nrml = fptr_->UnitNormal();
  }


template<uint32_t dim, template<uint32_t> class CELL>
Point<dim>  FiniteElementPolicy<dim,CELL>::UnitNormal() const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    return Point<dim>( fptr_->UnitNormal() );
  }


/**

Method calculates the unit normal that is outward pointing to to the queried face.

@test SKM 22/9/2014: fixed 3D normals so that they are outward pointing.

*/
template<uint32_t dim, template<uint32_t> class CELL>
Point<dim>  FiniteElementPolicy<dim,CELL>::UnitNormalToFace( uint32_t face ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    vector<double> un(dim,0.);
    fptr_->UnitNormalToFace( face, un );
    return Point<dim>(un);

  } // end UnitNormalToFace



/** 
    Fastest version to compute outward-pointing normal to target face.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void  FiniteElementPolicy<dim,CELL>::UnitNormalToFace( uint32_t face, vector<double>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    fptr_->UnitNormalToFace( face, nrml );
  }



template<uint32_t dim, template<uint32_t> class CELL>
void  FiniteElementPolicy<dim,CELL>::UnitNormalToFace( uint32_t face, VectorVariable<dim>& nrml ) const
  {
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    vector<double> un(dim,0.);
    fptr_->UnitNormalToFace( face, un );
    for ( uint32_t d(0u); d < dim; ++d )
      nrml(d) = un[d];
  }



    /// outputs finite element and discretised variable to VTK file
template<uint32_t dim, template<uint32_t> class CELL>
void  FiniteElementPolicy<dim,CELL>::OutputPropertyToVTK( const csmp::Index& key,
                                                          const char* file_name, const char* var_name ) const
  {
    const CELL<dim>* cell( static_cast<const CELL<dim>*>(this) );
    
    // to trigger creation of coordinate matrix
    auto idx = cell->Idx();
    cell->Idx( idx+1 );
    assert( fptr_ != nullptr );
    CoordinateMatrix();
    // resetting
    cell->Idx( idx-1 );

    // put property value into dense matrix (variable dimensions x nodes)
    assert( key.place != REGION );
    assert( key.place != ELEMENT_INTEGRATION_POINT ); // TODO: nice to have, perhaps extrapolate to nodes?
    assert( key.place != SECTOR_INTEGRATION_POINT ); // TODO: nice to have, perhaps extrapolate to nodes?
    assert( key.place != FACET_INTEGRATION_POINT ); // TODO: nice to have, perhaps extrapolate to nodes?
    
    DenseMatrix<DM_MIN> PROPMAT;
    
		if ( key.place == ELEMENT || key.place == FACE || key.place == INTER_FACE )
      {
        if ( key.type == SCALAR ) {
             PROPMAT.Resize(1,1);
             double val = cell->Read( key );
             for ( uint32_t i{0U}; i<cell->Nodes(); ++i )
               PROPMAT(0,i) = val;
          }
        else if ( key.type == VECTOR ) {
             PROPMAT.Resize(dim,1);
             VectorVariable<dim>  vc;
             cell->Read( key, vc );
             for ( uint32_t j{0U}; j<dim; ++j )
               PROPMAT(j,0) = vc[j];
          }
        else if ( key.type == TENSOR ) {
             PROPMAT.Resize(dim*dim,1);
             TensorVariable<dim>  ts;
             cell->Read( key, ts );
             uint32_t count{ 0U };
             for ( uint32_t j{0U}; j<dim; ++j )
               for ( uint32_t k{0U}; k<dim; ++k )
                 PROPMAT(count++,0) = ts(j,k);
          }
        else if ( key.type == ARRAY ) {
             ArrayVariable ar( key.dataDepth );
             cell->Read( key, ar );
             PROPMAT.Resize(ar.Size(),1);
             for ( uint32_t j{0U}; j<ar.Size(); ++j )
               PROPMAT(j,0) = ar(j);
          }
        else if ( key.type == FLAGGEDARRAY ) {
             FlaggedArrayVariable far( key.dataDepth );
             cell->Read( key, far );
             PROPMAT.Resize(far.Size(),1);
            for ( uint32_t j{0U}; j<far.Size(); ++j )
              PROPMAT(j,0) = far(j);
          }
      }
		else if ( key.place == NODE ) // if the operand is placed on the constraint-points
      {
        if ( key.type == SCALAR ) {
             PROPMAT.Resize(1,cell->Nodes());
             for ( uint32_t i{0U}; i<cell->Nodes(); ++i )
               PROPMAT(0,i) = cell->N(i)->Read( key );
          }
        else if ( key.type == VECTOR ) {
             PROPMAT.Resize(dim,cell->Nodes());
             VectorVariable<dim>  vc;
             for ( uint32_t i{0U}; i<cell->Nodes(); ++i ) {
                   cell->N(i)->Read( key, vc );
                   for ( uint32_t j{0U}; j<dim; ++j )
                     PROPMAT(j,i) = vc[j];
               }
          }
        else if ( key.type == TENSOR ) {
             PROPMAT.Resize(dim*dim,cell->Nodes());
             TensorVariable<dim>  ts;
             for ( uint32_t i{0U}; i<cell->Nodes(); ++i ) {
                  cell->N(i)->Read( key, ts );
                  uint32_t count{ 0U };
                  for ( uint32_t j{0U}; j<dim; ++j )
                    for ( uint32_t k{0U}; k<dim; ++k )
                      PROPMAT(count++,i) = ts(j,k);
               }
          }
        else if ( key.type == ARRAY ) {
             ArrayVariable ar( key.dataDepth );
             for ( uint32_t i{0U}; i<cell->Nodes(); ++i ) {
                  cell->N(i)->Read( key, ar );
                  PROPMAT.Resize(ar.Size(),cell->Nodes());
                  for ( uint32_t j{0U}; j<ar.Size(); ++j )
                    PROPMAT(j,i) = ar(j);
               }
          }
        else if ( key.type == FLAGGEDARRAY ) {
             FlaggedArrayVariable far( key.dataDepth );
             PROPMAT.Resize(far.Size(),cell->Nodes());
             for ( uint32_t i{0U}; i<cell->Nodes(); ++i ) {
                  cell->N(i)->Read( key, far );
                  for ( uint32_t j{0U}; j<far.Size(); ++j )
                    PROPMAT(j,i) = far(j);
               }
          }
		  }

    // replacing whitespace in string with underscores
    string variable_name{ var_name };
    for ( auto& it : variable_name )
      if ( it ==' ' || it == '\t' || it == '\n' ) it = '_';
      
    fptr_->OutputNodeDataToVTK( file_name, variable_name.c_str(), PROPMAT );

} //end OutDataToVTK



// ELEMENT

// vector of properties from integration points
// scalar
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Element>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Element>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, uint32_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, uint32_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, uint32_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, uint32_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, uint32_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, uint32_t, TensorVariable<3U>&) const;

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
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Face>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Face>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, uint32_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, uint32_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, uint32_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, uint32_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, uint32_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, uint32_t, TensorVariable<3U>&) const;

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
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, uint32_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, uint32_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, uint32_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, uint32_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, uint32_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, uint32_t, TensorVariable<3U>&) const;

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



// EDGE

// vector of properties from integration points
// scalar
template void FiniteElementPolicy<1U,Edge>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<2U,Edge>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
template void FiniteElementPolicy<3U,Edge>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, vector<ScalarVariable>&) const;
// array
template void FiniteElementPolicy<1U,Edge>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<2U,Edge>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
template void FiniteElementPolicy<3U,Edge>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementPolicy<1U,Edge>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<2U,Edge>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
template void FiniteElementPolicy<3U,Edge>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementPolicy<1U,Edge>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, vector<VectorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Edge>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, vector<VectorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Edge>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementPolicy<1U,Edge>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, vector<TensorVariable<1U> >&) const;
template void FiniteElementPolicy<2U,Edge>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, vector<TensorVariable<2U> >&) const;
template void FiniteElementPolicy<3U,Edge>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const vector<double>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const vector<double>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const vector<double>& xyz, TensorVariable<3U>& var) const;

// property value at integration point
// scalar
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, uint32_t, ScalarVariable&) const;
// array
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, uint32_t, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, uint32_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, uint32_t, VectorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, uint32_t, VectorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, uint32_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, uint32_t, TensorVariable<1U>&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, uint32_t, TensorVariable<2U>&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, uint32_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementPolicy<1U,Edge>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementPolicy<2U,Edge>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementPolicy<3U,Edge>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;


template class FiniteElementPolicy<1U,Element>;
template class FiniteElementPolicy<1U,Face>;
template class FiniteElementPolicy<1U,InterFace>;
template class FiniteElementPolicy<1U,Edge>;

template class FiniteElementPolicy<2U,Element>;
template class FiniteElementPolicy<2U,Face>;
template class FiniteElementPolicy<2U,InterFace>;
template class FiniteElementPolicy<2U,Edge>;

template class FiniteElementPolicy<3U,Element>;
template class FiniteElementPolicy<3U,Face>;
template class FiniteElementPolicy<3U,InterFace>;
template class FiniteElementPolicy<3U,Edge>;

} // end csmp
