#include "BoundaryFluxVisitor.h"
#include "Boundary.h"
#include "PropertyDatabase.h"
#include "VectorVariable.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp{

template <uint32_t dim>
BoundaryFluxVisitor<dim>::BoundaryFluxVisitor( const PropertyDatabase<dim>& pdb, const char* property_name, const char* result_property_name )
: Visitor<dim>( BOUNDARY, BOUNDARY ), property_key_( pdb.StorageKey(property_name) ), result_property_key_( pdb.StorageKey(result_property_name) )
  {
    if( property_key_.place != ELEMENT || property_key_.type != VECTOR )
      throw csmp::Exception( ERROR, "BoundaryFluxVisitor", "Property has to be a vector placed on the element" );
    if( result_property_key_.place != BOUNDARY || result_property_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "BoundaryFluxVisitor", "Result property has to be a scalar placed on the boundary" );
  }

template <uint32_t dim>
void BoundaryFluxVisitor<dim>::Visit(Boundary<dim>* bd)
{
  double total_flux(0.);
  VectorVariable<dim> v(PLAIN, 0.), face_unit_normal(PLAIN, 0.);
  typename vector<Face<dim>*>::const_iterator facesEnd( bd->ElementsEnd() );

  for( typename vector<Face<dim>*>::const_iterator it( bd->ElementsBegin() ); it != facesEnd; ++it )
    {
      (*it)->Parent(INSIDE)->Read(property_key_, v);
      (*it)->UnitNormal(face_unit_normal);
      Point<dim> p( face_unit_normal.P() );
      total_flux += v.DotProduct(p)*(*it)->Volume();
    }
  bd->Store( result_property_key_, makeScalar(PLAIN, total_flux) );
}

template class BoundaryFluxVisitor<2>;
template class BoundaryFluxVisitor<3>;

}// csmp
