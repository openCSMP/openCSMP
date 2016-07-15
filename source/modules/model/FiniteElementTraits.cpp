#include "FiniteElementTraits.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

namespace csmp {

template<size_t dim, template<size_t> class SIMPLEX>
FiniteElementTraits<dim,SIMPLEX>::FiniteElementTraits()
{
}

template<size_t dim, template<size_t> class SIMPLEX>
CSMP_FEM_TYPE FiniteElementTraits<dim,SIMPLEX>::FE_Type() const
  {
    assert( static_cast<const SIMPLEX<dim>*>(this)->FE() != NULL );
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->ElementType();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::IsLineElement() const
  {
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->IsLineElement();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::IsSurfaceElement() const
  {
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->IsSurfaceElement();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::IsVolumeElement() const
  {
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->IsVolumeElement();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t  FiniteElementTraits<dim,SIMPLEX>::Interpolation() const
  {
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->Interpolation();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool   FiniteElementTraits<dim,SIMPLEX>::UsesLocalCoordinates() const
  {
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->UsesLocalCoordinates();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t  FiniteElementTraits<dim,SIMPLEX>::Segments() const
  {
    assert( static_cast<const SIMPLEX<dim>*>(this)->FE() != NULL );
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->Segments();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t  FiniteElementTraits<dim,SIMPLEX>::Faces() const
  {
    assert( static_cast<const SIMPLEX<dim>*>(this)->FE() != NULL );
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->Faces();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t  FiniteElementTraits<dim,SIMPLEX>::Neighbors() const
  {
    assert( static_cast<const SIMPLEX<dim>*>(this)->FE() != NULL );
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->Neighbors();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t  FiniteElementTraits<dim,SIMPLEX>::Nodes() const
  {
    assert( static_cast<const SIMPLEX<dim>*>(this)->FE() != NULL );
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->Nodes();
  }







/**

Returns a matrix with 'nodes'-rows and 'coordinate-directions' columns.
This Meschach++ matrix defines the positions of the elements nodes for
the finite-element matrix assembly. Since the number of element nodes
may vary among different elements types, the number of rows in XY may
also vary from element to element.

@param XY A DenseMatrix<DM_MIN> class object (value type fT). This matrix is dynamically
resized if necessary but must have been constructed with a finite size
before passing it to CoordinateMatrix().

@return The node coordinates are returned into the supplied matrix.

@section application Application

Finite-element forms of differential equations require the global node
coordinates of the element to calculate the element constribution to the
global solution matrix. If the element uses local coordinates, the global
node coordinates will still be required to compute Jacobian (coordinate-
transformation) matrix.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementTraits<dim,SIMPLEX>::CoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    const size_t n_nodes( Nodes());
    XY.Resize( n_nodes, dim );
    for ( size_t i=0U; i<n_nodes; ++i )
        XY.AssignRow( i, eptr->N(i)->Coordinate() );

  } // end CoordinateMatrix



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
void  FiniteElementTraits<dim,SIMPLEX>::CoordinateMatrix() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    // if this method is called the second time on the same element it will do nothing
    if ( eptr->FE()->CurrentID() == eptr->Idx() ) return;

    const size_t  n_nodes( Nodes() );
    eptr->FE()->XY.Resize( n_nodes,dim );
    for ( size_t i=0U; i< n_nodes; ++i )
        eptr->FE()->XY.AssignRow( i, eptr->N(i)->Coordinate() );
    
    // resetting the element id on connected finite element
    eptr->FE()->CurrentID( eptr->Idx() );

  } // end CoordinateMatrix



/// Returns number of integration points per element, 0 if no FiniteElement assigned
template<size_t dim, template<size_t> class SIMPLEX>
size_t FiniteElementTraits<dim,SIMPLEX>::IntegrationPoints() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    if( eptr->FE() != NULL )
        return eptr->FE()->IntegrationPoints();
    return 0;
  }

/**
     Returns the position of the element integration point i in physical space.
     The coordinate value is transformed into global coordinates by the
     finite element. Therefore the coordinate matrix must be initialized.
*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementTraits<dim,SIMPLEX>::IntegrationPoint( size_t ip ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    eptr->CoordinateMatrix();
    // DNR matrix is used here as a temporary
    // so that no extra matrix has to be created (method uses NRST internally)
    eptr->FE()->IntegrationPoint( ip, eptr->FE()->DNR );
    return Point<dim>( eptr->FE()->DNR );
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::WeightAtIntegrationPoint( size_t i ) const
  {
    return static_cast<const SIMPLEX<dim>*>(this)->FE()->WeightAtIntegrationPoint(i);
  }


// SHAPE FUNCTIONS AT DIFFERENT POINTS

/**
    Returns the value of the interpolation functions at a point in physical space.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementTraits<dim,SIMPLEX>::N_AtGlobalPoint( std::vector<double64>& Nn, const std::vector<double64>& xyz ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    eptr->FE()->N( Nn, xyz );
  }

template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementTraits<dim,SIMPLEX>::N_AtIntegrationPoint( size_t ipoint, std::vector<double64>& Nn ) const
  {
    const SIMPLEX<dim>* eptr( static_cast< const SIMPLEX<dim>* const>(this) );
    eptr->FE()->N_AtIntegrationPoint( ipoint, Nn );
  }


template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementTraits<dim,SIMPLEX>::N_AtBaryCenter( std::vector<double64>& Nn ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    if ( !eptr->FE()->Isoparametric() )
        CoordinateMatrix();
    eptr->FE()->N_AtBaryCenter( Nn );
  }


// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS

template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementTraits<dim,SIMPLEX>::dN( DenseMatrix<DM_MIN>& M ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    eptr->FE()->dN( M );
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::dN_AtNode( DenseMatrix<DM_MIN>& M, size_t nd, size_t dof )  const
  { 
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    if ( dof == 1U )
        return eptr->FE()->dN_AtNode( M, nd );
    double64 detJ = eptr->FE()->dN_AtNode( M, nd );
    dof == 2 ? dN_To2DOF( Nodes(), M ) : dN_To3DOF( Nodes(), M );
    return detJ;
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::dN_AtBaryCenter( DenseMatrix<DM_MIN>& M, size_t dof ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    if ( dof == 1U )
        return eptr->FE()->dN_AtBarycenter( M );
    double64 detJ = eptr->FE()->dN_AtBarycenter( M );
    dof == 2U ? dN_To2DOF( Nodes(), M ) : dN_To3DOF( Nodes(), M );
    return detJ;
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M,
                                                                  size_t gp, size_t dof )  const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    if ( dof == 1U )
        return eptr->FE()->dN_AtIntegrationPoint( M, gp );
    double64 detJ = eptr->FE()->dN_AtIntegrationPoint( M, gp );
    dof == 2U ? dN_To2DOF( Nodes(), M ) : dN_To3DOF( Nodes(), M );
    return detJ;
  }

// INTEGRATION

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::det_JINV_AtIntegrationPoint( size_t ipoint ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    eptr->FE()->JacobianAtIntegrationPoint( ipoint );
    return eptr->FE()->JacobianInverse();
  }

template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementTraits<dim,SIMPLEX>::IntegralNN( DenseMatrix<DM_MIN>& M ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    eptr->FE()->IntegralNN( M );
  }




// PROPERTY INTERPOLATION, EXTRAPOLATION AND INTEGRATIONS

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::PropertyIntegral( const csmp::Index& prop_key ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( prop_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "FiniteElementTraits<dim,SIMPLEX>::PropertyIntegral",
                                            "This method only integrates scalar properties" );

    if ( eptr->FE() == NULL ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementTraits<dim,SIMPLEX>::PropertyIntegral:",
                                 "encountered element with invalid element pointer.");
       }

    // if the element uses analytical integration, a linear element is assumed and the 
    // integration is performed for a property value interpolated to the element's barycentre
    if ( !eptr->FE()->UsesLocalCoordinates() ||
         (prop_key.place == ELEMENT || prop_key.place == FACE || prop_key.place == INTER_FACE) )
      {
        ScalarVariable  sc;
        eptr->PropertyValueAtBaryCenter( prop_key, sc );
        return eptr->Volume() * sc.Value();
      }

    if ( prop_key.place == ELEMENT_INTEGRATION_POINT )
      {
         double64 sum(0.);
         const size_t n_integration_points(eptr->FE()->IntegrationPoints());
         CoordinateMatrix();
         for ( size_t i=0U; i <n_integration_points; i++ )
         {
             eptr->FE()->JacobianAtIntegrationPoint( i );
             sum += eptr->FE()->JacobianInverse() * eptr->FE()->WeightAtIntegrationPoint(i) * eptr->Read( i, prop_key );
         }
         return sum;
      }

    // NODE properties
    // ---------------  
    // if the element is isoparametric we loop over the integration points interpolating
    // the property values to these positions and then integrating via multiplication with 
    // weights
    double64  sumN, sumI(0.);
    const size_t n_integration_points(eptr->FE()->IntegrationPoints());

    CoordinateMatrix();
    for ( size_t i=0U; i <n_integration_points; i++ )
      {
        sumN = 0.;
        eptr->FE()->JacobianAtIntegrationPoint(i);
        eptr->FE()->N_AtIntegrationPoint( i, eptr->FE()->NRST );
        for ( size_t j=0; j<Nodes(); j++ )
            sumN += eptr->FE()->NRST[j] * eptr->N(j)->Read( prop_key );
        sumI += sumN * eptr->FE()->JacobianInverse() * eptr->FE()->WeightAtIntegrationPoint(i);
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
void FiniteElementTraits<dim,SIMPLEX>::PropertyValueAt( const csmp::Index& idx,
                                                        const std::vector<double64>& xyz,
                                                        Var& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    Var  temp;
    if ( idx.place == ELEMENT ) {
         eptr->Read( idx, temp );
         return;
    }
    if ( idx.place != NODE ) {
       std::cerr <<"\nFiniteElementTraits<"<< dim;
       std::cerr <<">::PropertyValueAt: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       throw std::domain_error("FiniteElementTraits<dim,SIMPLEX>::PropertyValueAt");
    }

    temp.Resize( idx.dataDepth );
    var.Resize( idx.dataDepth, 0. );
    // initialisation for accumulation
    var = 0.;

    CoordinateMatrix();
    eptr->FE()->N( eptr->FE()->NRST, xyz );

    const size_t  n_nodes(Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
    {
       eptr->N(i)->Read( idx, temp );
       var += temp * eptr->FE()->NRST[i];
    }
  }


/**

Node variables are interpolated to barycenter.
Integration point variables (ELEMENT_INTEGRATION_POINT) are averaged.

@attention If variable is an element property, method will
simply return this value.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void FiniteElementTraits<dim,SIMPLEX>::PropertyValueAtBaryCenter( const csmp::Index& idx,
                                                                  Var& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    Var  temp;
    if ( idx.place == ELEMENT ) {
         eptr->Read( idx, var );
         return;
    }
    if ( idx.place != NODE  and idx.place != ELEMENT_INTEGRATION_POINT and idx.place != SECTOR_INTEGRATION_POINT) {
       std::cerr <<"\nFiniteElementTraits<"<< dim;
       std::cerr <<">::PropertyValueAtBaryCenter: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       std::cerr <<"values of ELEMENT_INTEGRATION_POINT properties are averaged."<< std::endl;
       throw std::domain_error("FiniteElementTraits<dim,SIMPLEX>::PropertyValueAtBaryCenter");
    }

    temp.Resize( idx.dataDepth );
    var.Resize( idx.dataDepth, 0. );
    // initialisation for accumulation
    var = 0.;

    if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
       const size_t n_integration_points(IntegrationPoints());
       for ( size_t i=0U; i < n_integration_points; i++ ) {
            eptr->Read( i, idx, temp );
            var += temp;
         }
       var /= static_cast<double64>(eptr->FE()->IntegrationPoints());
       return;
    }
    if ( idx.place == SECTOR_INTEGRATION_POINT ) {
       const size_t n_sector_integration_points(eptr->FV_Stencil()->Sectors());
       for ( size_t i=0U; i < n_sector_integration_points; i++ ) {
            eptr->Read( i, 0U, idx, temp );
            var += temp;
         }
       var /= static_cast<double64>(n_sector_integration_points);
       return;
    }

    if ( !eptr->FE()->UsesLocalCoordinates() )
        CoordinateMatrix();
    eptr->FE()->N_AtBaryCenter( eptr->FE()->NRST );

    const size_t  n_nodes(Nodes());
    for ( size_t i=0U; i<n_nodes; i++ ) {
       eptr->N(i)->Read( idx, temp );
       var += temp * eptr->FE()->NRST[i];
    }
  }

/**
Interpolates node properties to the integration points.

@attention If variable is an element property, method will
simply return this value.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void FiniteElementTraits<dim,SIMPLEX>::PropertyValueAtIntegrationPoint( const csmp::Index& idx,
                                                                        size_t ip,
                                                                        Var& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    Var  temp;
    if ( idx.place == ELEMENT ) {
         eptr->Read( idx, temp );
         return;
    }
    if ( idx.place != NODE )
    {
       std::cerr <<"\nFiniteElementTraits<"<< dim;
       std::cerr <<">::PropertyValueAtIntegrationPoint: This method only interpolates the ";
       std::cerr <<"values of NODE properties."<< std::endl;
       throw std::domain_error("FiniteElementTraits<dim,SIMPLEX>::PropertyValueAtIntegrationPoint");
    }

    // potential size and value adjustments for array variables
    temp.Resize( idx.dataDepth );
    var.Resize( idx.dataDepth, 0. );
    // initialisation for accumulation
    var = 0.;

    eptr->FE()->N_AtIntegrationPoint( ip, eptr->FE()->NRST );

    const size_t  n_nodes(Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
    {
       eptr->N(i)->Read( idx, temp );
       var += temp * eptr->FE()->NRST[i];
    }
  }

/**
Interpolates SCALAR node properties to the integration points.

@attention If variable is an element property, method will
simply return this value.
*/
template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::PropertyValueAtIntegrationPoint( const csmp::Index& idx,
                                                                            size_t ip ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    double64  var(0.);
    if ( idx.place == ELEMENT && idx.type == SCALAR )
        return eptr->Read( idx );
    if ( idx.place != NODE or idx.type != SCALAR ) {
        std::cerr <<"\nFiniteElementTraits<"<< dim;
        std::cerr <<">::PropertyValueAtIntegrationPoint: This method only interpolates ";
        std::cerr <<"scalar NODE properties."<< std::endl;
        throw std::domain_error("FiniteElementTraits<dim,SIMPLEX>::PropertyValueAtIntegrationPoint");
     }

    eptr->FE()->N_AtIntegrationPoint( ip, eptr->FE()->NRST );

    const size_t  n_nodes(Nodes());
    for ( size_t i=0U; i<n_nodes; i++ )
      var += eptr->FE()->NRST[i] * eptr->N(i)->Read( idx );

   return var;
  }

/**

Returns the ELEMENT_INTEGRATION_POINT variable values into the parameter
vector.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void  FiniteElementTraits<dim,SIMPLEX>::IntegrationPointPropertyVector( const csmp::Index& idx,
                                                                        std::vector<Var>& var ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( eptr->FE() == NULL ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementTraits<dim,SIMPLEX>::IntegrationPointPropertyVector:",
                                 "encountered element with invalid element pointer.");
       }

    assert( eptr->FE()->UsesLocalCoordinates() );
    assert( eptr->FE()->IntegrationPoints() > 0U );

    // resizing V if necessary
    const size_t n_integration_points(eptr->FE()->IntegrationPoints());
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
to its nodes

*/


template<size_t dim, template<size_t> class SIMPLEX>
void FiniteElementTraits<dim,SIMPLEX>
::ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                              const std::vector<double64>& IVAR,
                                              std::vector<double64>& NVAR ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if( eptr->FE() == NULL )
        throw csmp::Exception( ERROR,
                               "FiniteElementTraits<dim,SIMPLEX>::ExtrapolateIntegrationPointVariableToNodes:",
                               "encountered element with invalid element pointer.");
    if( !eptr->FE()->UsesLocalCoordinates() )
        throw csmp::Exception( ERROR,
                               "FiniteElementTraits<dim,SIMPLEX>::ExtrapolateIntegrationPointVariableToNodes:",
                               "Extrapolation is not supported; only the isoparametric elements have integration points.");

    CoordinateMatrix();

    // NB: the following method checks that the supplied arrays have the right size
    eptr->FE()->ExtrapolateIntegrationPointVariableToNodes( nvars, IVAR, NVAR );

  }



/**

Outputs to vector which must be of a size which depends on where the
variable lives (size=nodes, size=element=1 etc.)

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<class Var>
void  FiniteElementTraits<dim,SIMPLEX>::NodePropertyVector( const csmp::Index& idx,
                                                            std::vector<Var>& V ) const
   {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if( eptr->FE() == NULL )
        throw csmp::Exception( ERROR,
                               "FiniteElementTraits<dim,SIMPLEX>::NodePropertyVector:",
                               "encountered element with invalid element pointer.");

    if ( idx.place != NODE ) {
         std::cerr <<"\nFiniteElementTraits<"<< dim;
         std::cerr <<">::NodePropertyVector: Requested property ";
         std::cerr <<"is not placed on the nodes; property Index: "<< std::endl;
         idx.Out();
         return;
      }

    // resizing V if necessary
    const size_t  n_nodes(Nodes());
    V.resize(n_nodes);

    for ( size_t i=0U; i<n_nodes; i++ )
      eptr->N(i)->Read( idx, V[i] );

   } // end













// ELEMENT GEOMETRY
template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::AspectRatio() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->AspectRatio();
  }

template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::InnerRadius() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->InnerRadius();
  }


template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::Volume() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    return eptr->FE()->Volume();
  }





/**

Measures the length of the element (in physical space) in a certain direction.

@param vecDirection direction in which the element is to be measured

@return Length of the element.

@sectin implementation Implementation

The main idea is to measure the height of the oriented bounding box
(oriented by the direction of vecDirection) which surrounds the element.
This is done by projecting each node onto the direction vector and calculating
the difference between the largest and smallest magnitude. The distance between both
projections will be given by the absolute value of this difference.

@section application Application

The length of the element in a certain dimension is used to weigh the error of the element,
when evaluating the quality of a certain mesh.
*/
template<size_t dim, template<size_t> class SIMPLEX>
double64  FiniteElementTraits<dim,SIMPLEX>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    //refresh coordinate matrix
    CoordinateMatrix();

    if ( dim == 1U ) return Volume();

    // find the projection magnitude (SKM??? - this is just the vecDirection.Length() ?)
    const double64 fMagnitudeOfDirection(vecDirection.Length()); // sqrt(vecDirection & vecDirection);

    double64 fMinTemp( static_cast<double64>( DBL_MAX) );
    double64 fMaxTemp( static_cast<double64>(-DBL_MAX) );

    //project each node onto the direction vector (calculate only magnitude)
    for ( size_t i=0; i<Nodes(); i++ )
    {
        // fTemp is the projection of the vector (0,0,0)-node(i) on the vector direction
        double64 fTemp(vecDirection.DotProduct( eptr->N(i)->Coordinate() ));
        fTemp /= fMagnitudeOfDirection;

        // update minimum value
        fMinTemp = std::min(fMinTemp, fTemp);
        // update maximum value
        fMaxTemp = std::max(fMaxTemp, fTemp);
    }

    //substract magnitudes
    return fMaxTemp - fMinTemp;
  }


/**

BaryCenter() calculates the node coordinate average for the element. This
coordinate value is equivalent to the center of gravity of the Element
type.

@section input Input Arguments

The barycentre is returned into a CSMP vector variable a reference to which
is supplied as single argument.

@section application Application

The element barycentre could be used for instance to output element material
properties from the Model as point data. This is done if you use the
Model OutputToTextFile() methods.

@return The VARIABLE_FLAG of the returned vector variable will not be changed by
BaryCentre().

@test O.K. SKM25/8/14 after refactoring loop

*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementTraits<dim,SIMPLEX>::BaryCenter() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( eptr->FE() == NULL ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementTraits<dim,SIMPLEX>::BaryCenter:",
                                 "encountered element with invalid element pointer.");
       }

    // for the local coordinate system
    if ( eptr->FE()->Isoparametric() ) {
         eptr->FE()->N_AtBaryCenter( eptr->FE()->NRST );
         Point<dim> pt(eptr->N(0U)->Coordinate() * eptr->FE()->NRST[0U]);
         const size_t n_nodes(eptr->FE()->NRST.size());
         for ( size_t i=1U; i<n_nodes; ++i )
           pt += eptr->N(i)->Coordinate() * eptr->FE()->NRST[i];

         return pt;
      }

    // global coordinate system
    Point<dim>    pt(eptr->N(0U)->Coordinate());
    const size_t  n_nodes(Nodes());
    for ( size_t i=1U; i<n_nodes; ++i )
      pt += eptr->N(i)->Coordinate();

    return pt / static_cast<double64>(Nodes());
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
double64  FiniteElementTraits<dim,SIMPLEX>::FaceArea( size_t n ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( eptr->FE() == NULL ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementTraits<dim,SIMPLEX>::FaceArea:",
                                 "encountered element with invalid element pointer.");
       }

    if ( eptr->FE()->IsVolumeElement() ) {
         std::vector<size_t>  ids;
         eptr->FE()->NodesOfFace( n, ids );

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
    if ( eptr->FE()->IsSurfaceElement() ) {
         std::vector<size_t>  ids;
         eptr->FE()->NodesOfFace( n, ids );
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
Point<dim>  FiniteElementTraits<dim,SIMPLEX>::FaceBaryCenter( size_t face ) const
  {
    const SIMPLEX<dim>* e( static_cast<const SIMPLEX<dim>*>(this) );
    if ( e->FE() == NULL ) {
          e->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementTraits<dim,SIMPLEX>::FaceBaryCenter:",
                                "encountered element with invalid element pointer.");
       }
    e->FE()->NodesOfFace( face, e->FE()->IDX );
    const size_t face_nodes(e->FE()->IDX.size());
    Point<dim>  barycenter(e->N(e->FE()->IDX[0])->Coordinate());
    for ( size_t i=1U; i<face_nodes; ++i )
      barycenter += e->N(e->FE()->IDX[i])->Coordinate();
    barycenter /= static_cast<double64>(face_nodes);

    return barycenter;
  }




// SEGMENTS


template<size_t dim, template<size_t> class SIMPLEX>
double64 FiniteElementTraits<dim,SIMPLEX>::SegmentLength( size_t segm ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    if ( eptr->FE() == NULL ) {
          eptr->Out();
          throw csmp::Exception( ERROR,
                                 "FiniteElementTraits<dim,SIMPLEX>::SegmentLenghts:",
                                 "encountered element with invalid element pointer.");
       }

    std::vector<size_t>  snids;
    eptr->FE()->NodesOfSegment( segm, snids );

    assert( snids.size() == 2U );
    return eptr->N(snids[1])->Coordinate().DistanceTo( eptr->N(snids[0])->Coordinate() );

  } // end BoundarySegmentLength



    /// returns the length of all element segments=edges
template<size_t dim, template<size_t> class SIMPLEX>
inline void FiniteElementTraits<dim,SIMPLEX>::SegmentLengths( std::vector<double64>& lengths ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    eptr->FE()->EdgeLengths( lengths );
  } // end SegmentLengths


template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementTraits<dim,SIMPLEX>::SegmentMidPoint( size_t segm ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>* >(this) );
    eptr->FE()->NodesOfSegment( segm, eptr->FE()->IDX );
    assert( eptr->FE()->IDX.size() == 2U );
    return midPoint(eptr->N(eptr->FE()->IDX[0])->Coordinate(),eptr->N(eptr->FE()->IDX[1])->Coordinate());
  }







// UNIT NORMALS



/// computes outward pointing normal
template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementTraits<dim,SIMPLEX>::UnitNormal( VectorVariable<dim>& nrml ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    std::vector<double64> un( dim, 0. );
    eptr->FE()->UnitNormal(un);
    for ( size_t d(0); d < dim; ++d )
        nrml(d) = un[d];
  }

template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementTraits<dim,SIMPLEX>::UnitNormal( std::vector<double64>& nrml ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    nrml.resize( dim, 0. );
    eptr->FE()->UnitNormal(nrml);
  }

template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementTraits<dim,SIMPLEX>::UnitNormal() const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    CoordinateMatrix();
    std::vector<double64> un( dim, 0. );
    eptr->FE()->UnitNormal(un);
    return Point<dim>(un);
  }


/**

Method calculates the unit normal to the queried face and checks whether it is inward or
outward pointing. The normal will be modified before it is returned as a point so that
it is outward pointing.

@test SKM 22/9/2014: fixed 3D normals so that they are outward pointing.

*/
template<size_t dim, template<size_t> class SIMPLEX>
Point<dim>  FiniteElementTraits<dim,SIMPLEX>::UnitNormalToFace( size_t face ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    std::vector<double64> unrml;
    eptr->FE()->UnitNormalToFace( face, unrml );
    
    return Point<dim>(unrml);

  } // end UnitNormalToFace



/// computes outward pointing normal to target face
template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementTraits<dim,SIMPLEX>::UnitNormalToFace( size_t face, std::vector<double64>& nrml ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    Point<dim>  p(UnitNormalToFace( face ));
    nrml = p.Coordinates();
  }

template<size_t dim, template<size_t> class SIMPLEX>
void  FiniteElementTraits<dim,SIMPLEX>::UnitNormalToFace( size_t face, VectorVariable<dim>& nrml ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );
    assert( eptr->FE() != NULL );
    CoordinateMatrix();
    Point<dim>  p(UnitNormalToFace( face ));
    for ( size_t d(0); d < dim; ++d )
        nrml(d) = p[d];
  }



// ELEMENT

// vector of properties from nodes
// scalar
template void FiniteElementTraits<1U,Element>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementTraits<2U,Element>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementTraits<3U,Element>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
// array
template void FiniteElementTraits<1U,Element>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementTraits<2U,Element>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementTraits<3U,Element>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementTraits<1U,Element>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementTraits<2U,Element>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementTraits<3U,Element>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementTraits<1U,Element>
::NodePropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Element>
::NodePropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Element>
::NodePropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementTraits<1U,Element>
::NodePropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Element>
::NodePropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Element>
::NodePropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// vector of properties from integration points
// scalar
template void FiniteElementTraits<1U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementTraits<2U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
template void FiniteElementTraits<3U,Element>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable>&) const;
// array
template void FiniteElementTraits<1U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementTraits<2U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
template void FiniteElementTraits<3U,Element>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable>&) const;
// flagged array
template void FiniteElementTraits<1U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementTraits<2U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
template void FiniteElementTraits<3U,Element>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable>&) const;
// vector
template void FiniteElementTraits<1U,Element>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Element>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Element>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementTraits<1U,Element>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Element>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Element>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementTraits<1U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementTraits<1U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementTraits<1U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementTraits<1U,Element>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementTraits<1U,Element>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<3U>& var) const;


// property value at integration point
// scalar
template void FiniteElementTraits<1U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
// array
template void FiniteElementTraits<1U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
// flagged array
template void FiniteElementTraits<1U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementTraits<1U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, size_t, VectorVariable<1U>&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, size_t, VectorVariable<2U>&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, size_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementTraits<1U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, size_t, TensorVariable<1U>&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, size_t, TensorVariable<2U>&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, size_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementTraits<1U,Element>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementTraits<1U,Element>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementTraits<1U,Element>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementTraits<1U,Element>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementTraits<1U,Element>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Element>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Element>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;






// FACE

// vector of properties from nodes
// scalar
template void FiniteElementTraits<1U,Face>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<2U,Face>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<3U,Face>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
// array
template void FiniteElementTraits<1U,Face>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<2U,Face>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<3U,Face>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
// flagged array
template void FiniteElementTraits<1U,Face>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<2U,Face>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<3U,Face>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
// vector
template void FiniteElementTraits<1U,Face>
::NodePropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Face>
::NodePropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Face>
::NodePropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementTraits<1U,Face>
::NodePropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Face>
::NodePropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Face>
::NodePropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// vector of properties from integration points
// scalar
template void FiniteElementTraits<1U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<2U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<3U,Face>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
// array
template void FiniteElementTraits<1U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<2U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<3U,Face>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
// flagged array
template void FiniteElementTraits<1U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<2U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<3U,Face>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
// vector
template void FiniteElementTraits<1U,Face>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Face>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Face>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementTraits<1U,Face>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementTraits<2U,Face>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementTraits<3U,Face>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementTraits<1U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementTraits<1U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementTraits<1U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementTraits<1U,Face>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementTraits<1U,Face>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<3U>& var) const;


// property value at integration point
// scalar
template void FiniteElementTraits<1U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
// array
template void FiniteElementTraits<1U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
// flagged array
template void FiniteElementTraits<1U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementTraits<1U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, size_t, VectorVariable<1U>&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, size_t, VectorVariable<2U>&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, size_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementTraits<1U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, size_t, TensorVariable<1U>&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, size_t, TensorVariable<2U>&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, size_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementTraits<1U,Face>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementTraits<1U,Face>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementTraits<1U,Face>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementTraits<1U,Face>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementTraits<1U,Face>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementTraits<2U,Face>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementTraits<3U,Face>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;



// INTERFACE


// vector of properties from nodes
// scalar
template void FiniteElementTraits<1U,InterFace>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<2U,InterFace>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<3U,InterFace>
::NodePropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
// array
template void FiniteElementTraits<1U,InterFace>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<2U,InterFace>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<3U,InterFace>
::NodePropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
// flagged array
template void FiniteElementTraits<1U,InterFace>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<2U,InterFace>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<3U,InterFace>
::NodePropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
// vector
template void FiniteElementTraits<1U,InterFace>
::NodePropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementTraits<2U,InterFace>
::NodePropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementTraits<3U,InterFace>
::NodePropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementTraits<1U,InterFace>
::NodePropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementTraits<2U,InterFace>
::NodePropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementTraits<3U,InterFace>
::NodePropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// vector of properties from integration points
// scalar
template void FiniteElementTraits<1U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<2U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
template void FiniteElementTraits<3U,InterFace>
::IntegrationPointPropertyVector<ScalarVariable >(const csmp::Index&, std::vector<ScalarVariable >&) const;
// array
template void FiniteElementTraits<1U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<2U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
template void FiniteElementTraits<3U,InterFace>
::IntegrationPointPropertyVector<ArrayVariable >(const csmp::Index&, std::vector<ArrayVariable >&) const;
// flagged array
template void FiniteElementTraits<1U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<2U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
template void FiniteElementTraits<3U,InterFace>
::IntegrationPointPropertyVector<FlaggedArrayVariable >(const csmp::Index&, std::vector<FlaggedArrayVariable >&) const;
// vector
template void FiniteElementTraits<1U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<1U> >(const csmp::Index&, std::vector<VectorVariable<1U> >&) const;
template void FiniteElementTraits<2U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<2U> >(const csmp::Index&, std::vector<VectorVariable<2U> >&) const;
template void FiniteElementTraits<3U,InterFace>
::IntegrationPointPropertyVector<VectorVariable<3U> >(const csmp::Index&, std::vector<VectorVariable<3U> >&) const;
// tensor
template void FiniteElementTraits<1U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<1U> >(const csmp::Index&, std::vector<TensorVariable<1U> >&) const;
template void FiniteElementTraits<2U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<2U> >(const csmp::Index&, std::vector<TensorVariable<2U> >&) const;
template void FiniteElementTraits<3U,InterFace>
::IntegrationPointPropertyVector<TensorVariable<3U> >(const csmp::Index&, std::vector<TensorVariable<3U> >&) const;

// property value at point
// scalar
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAt<ScalarVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ScalarVariable& var) const;
// array
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAt<ArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, ArrayVariable& var) const;
// flagged array
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAt<FlaggedArrayVariable >(const csmp::Index& idx, const std::vector<double64>& xyz, FlaggedArrayVariable& var) const;
// vector
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAt<VectorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<1U>& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAt<VectorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<2U>& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAt<VectorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAt<TensorVariable<1U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<1U>& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAt<TensorVariable<2U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<2U>& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAt<TensorVariable<3U> >(const csmp::Index& idx, const std::vector<double64>& xyz, TensorVariable<3U>& var) const;


// property value at integration point
// scalar
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtIntegrationPoint<ScalarVariable >(const csmp::Index&, size_t, ScalarVariable&) const;
// array
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtIntegrationPoint<ArrayVariable >(const csmp::Index&, size_t, ArrayVariable&) const;
// flagged array
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtIntegrationPoint<FlaggedArrayVariable >(const csmp::Index&, size_t, FlaggedArrayVariable&) const;
// vector
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<1U> >(const csmp::Index&, size_t, VectorVariable<1U>&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<2U> >(const csmp::Index&, size_t, VectorVariable<2U>&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtIntegrationPoint<VectorVariable<3U> >(const csmp::Index&, size_t, VectorVariable<3U>&) const;
// tensor
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<1U> >(const csmp::Index&, size_t, TensorVariable<1U>&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<2U> >(const csmp::Index&, size_t, TensorVariable<2U>&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtIntegrationPoint<TensorVariable<3U> >(const csmp::Index&, size_t, TensorVariable<3U>&) const;

// property value at bary center
// scalar
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx,ScalarVariable& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtBaryCenter<ScalarVariable >(const csmp::Index& idx, ScalarVariable& var) const;
// array
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtBaryCenter<ArrayVariable >(const csmp::Index&, ArrayVariable&) const;
// flagged array
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtBaryCenter<FlaggedArrayVariable >(const csmp::Index&, FlaggedArrayVariable&) const;
// vector
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtBaryCenter<VectorVariable<1U> >(const csmp::Index& idx, VectorVariable<1U>& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtBaryCenter<VectorVariable<2U> >(const csmp::Index& idx, VectorVariable<2U>& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtBaryCenter<VectorVariable<3U> >(const csmp::Index& idx, VectorVariable<3U>& var) const;
// tensor
template void FiniteElementTraits<1U,InterFace>
::PropertyValueAtBaryCenter<TensorVariable<1U> >(const csmp::Index& idx, TensorVariable<1U>& var) const;
template void FiniteElementTraits<2U,InterFace>
::PropertyValueAtBaryCenter<TensorVariable<2U> >(const csmp::Index& idx, TensorVariable<2U>& var) const;
template void FiniteElementTraits<3U,InterFace>
::PropertyValueAtBaryCenter<TensorVariable<3U> >(const csmp::Index& idx, TensorVariable<3U>& var) const;



// SKM: TODO: REMOVE these methods

/// Returns whether unit normal points outward (wrt inner/outer parent elements)
template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::isCorrectUnitNormalOrientation( Element<dim>* innerParent ) const
  {
//    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    VectorVariable<dim> vc;
    UnitNormal(vc);

    return FromOutside( innerParent, vc );
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::isCorrectUnitNormalOrientation( Element<dim>* innerParent, size_t face_id ) const
  {
//    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    VectorVariable<dim> vc;
    innerParent->UnitNormalToFace( face_id, vc );

    return FromInside( innerParent, vc );
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::FromInside( Element<dim>* parentElement, VectorVariable<dim> face_unit_normal ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    const Point<dim> bcFace( eptr->BaryCenter() ), bcParent( parentElement->BaryCenter() );
    VectorVariable<dim> faceToParent( ANY, 0. );
    for( size_t d(0); d < dim; ++d )
        faceToParent(d) = bcParent[d] - bcFace[d];
    if( (faceToParent&face_unit_normal) < 0. )
        return true;
    return false;
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool FiniteElementTraits<dim,SIMPLEX>::FromOutside( Element<dim>* parentElement, VectorVariable<dim> face_unit_normal ) const
  {
    const SIMPLEX<dim>* eptr( static_cast<const SIMPLEX<dim>*>(this) );

    const Point<dim> bcFace( eptr->BaryCenter() ), bcParent( parentElement->BaryCenter() );
    VectorVariable<dim> faceToParent( ANY, 0. );
    for( size_t d(0); d < dim; ++d )
        faceToParent(d) = bcParent[d] - bcFace[d];
    if( (faceToParent&face_unit_normal) < 0. )
        return false;
    return true;
  }







template class FiniteElementTraits<1U,Element>;
template class FiniteElementTraits<1U,Face>;
template class FiniteElementTraits<1U,InterFace>;

template class FiniteElementTraits<2U,Element>;
template class FiniteElementTraits<2U,Face>;
template class FiniteElementTraits<2U,InterFace>;

template class FiniteElementTraits<3U,Element>;
template class FiniteElementTraits<3U,Face>;
template class FiniteElementTraits<3U,InterFace>;

} // end csmp
