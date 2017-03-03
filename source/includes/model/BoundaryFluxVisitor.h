#ifndef BOUNDARY_FLUX_VISITOR_H
#define BOUNDARY_FLUX_VISITOR_H

#include "Visitor.h"

namespace csmp
{
  class Index;
  template<size_t> class PropertyDatabase;

  template <size_t dim>
  class BoundaryFluxVisitor:public Visitor<dim>
  {
  public:
    explicit BoundaryFluxVisitor( const PropertyDatabase<dim>&,
                                  const char* property_name,
                                  const char* result_property_name );
    
    virtual void Visit( Boundary<dim>* );   

  private:
    const csmp::Index property_key_;
    csmp::Index result_property_key_;
      
  };
}

#endif

