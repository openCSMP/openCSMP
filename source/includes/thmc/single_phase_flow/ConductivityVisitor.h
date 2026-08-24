#ifndef CONDUCTIVITY_VISITOR_H
#define CONDUCTIVITY_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#include "Index.h"

namespace csmp
{

  struct Index;
  template<uint32_t> class Model;

  /// PL 2012, Single-phase steady conductivity & transient diffusivity
  /// JEM 2015, modifications and updates
///  SKM 18/1/24 templatised and fixed open MP issue
template<uint32_t dim, template<uint32_t> class CELL=Element>
  class ConductivityVisitor final : public Visitor<dim> {
  public:
      ConductivityVisitor( Model<dim>&,
                           const char* specific_saturated_hydraulic_conductivity, // this is without mult. by density
                           const char* permeability,
                           const char* viscosity,
                           const char* density=NULL,
                           const char* saturated_hydraulic_conductivity=NULL,
                           const char* compressibility = NULL,
                           const char* porosity = NULL,
                           const char* specific_saturated_hydraulic_diffusivity=NULL); // this is without mult. by density

    /// for visiting elements, faces and interfaces
    virtual void Visit( CELL<dim>* ) override final;
    virtual void Visit( Model<dim>* ) override final {}

  protected:
    ConductivityVisitor();
    ConductivityVisitor( const ConductivityVisitor<dim>& );

  private:
    Model<dim>& model_;
    Index kKey_, muKey_, ctKey_, phiKey_, sshcKey_,shcKey_,rhoKey_,sshdKey_;
    ScalarVariable result,resultd, temp;
  };

} //csmp


#endif // CONDUCTIVITY_VISITOR_H
