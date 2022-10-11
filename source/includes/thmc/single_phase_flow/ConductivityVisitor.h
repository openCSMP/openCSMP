#ifndef CONDUCTIVITY_VISITOR_H
#define CONDUCTIVITY_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#if defined(_OPENMP )
#include "FiniteElementManager.h"
#endif

namespace csmp
{

  struct Index;
  template<uint32_t> class Model;

  /// PL 2012, Single-phase steady conductivity & transient diffusivity
  /// JEM 2015, modifications and updates
  template<uint32_t dim>
  class ConductivityVisitor : public Visitor<dim>
  {
  public:

      ConductivityVisitor(Model<dim>& model,
                           const char* specific_saturated_hydraulic_conductivity, // this is without mult. by density
                           const char* permeability,
                           const char* viscosity,
                           const char* density=NULL,
                           const char* saturated_hydraulic_conductivity=NULL,
                           const char* compressibility = NULL,
                           const char* porosity = NULL,
                           const char* specific_saturated_hydraulic_diffusivity=NULL); // this is without mult. by density

    virtual ~ConductivityVisitor() {}

    virtual void Visit(Element<dim>* element);
    void ComputeContribution(Element<dim>* element);
    virtual void Visit(Model<dim>* model);
    virtual void Visit(Region<dim>* region);

  protected:
    ConductivityVisitor();
    ConductivityVisitor( const ConductivityVisitor<dim>& );

  private:
    Model<dim>& model_;
    Index kKey_, muKey_, ctKey_, phiKey_, sshcKey_,shcKey_,rhoKey_,sshdKey_;
    ScalarVariable result,resultd, temp;
#if defined(_OPENMP )
    std::vector<ScalarVariable> thread_result_;
    std::vector<FiniteElementManager> femgrs_; // one manager per thread
#endif
  };

} //csmp


#endif // CONDUCTIVITY_VISITOR_H
