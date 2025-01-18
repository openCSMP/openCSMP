#ifndef STABILIZATION_PARAMETER_VISITOR_H
#define STABILIZATION_PARAMETER_VISITOR_H

#include "Visitor.h"
#include "PropertyDatabase.h"
#include "ScalarVariable.h"

namespace csmp {

/// for Stokes-lubrication equation where the same basis functions are used for pressure and velocity
template<uint32_t dim>
class StabilizationParameterVisitor : public Visitor<dim> {
  public:
    /// constant viscosity version
    StabilizationParameterVisitor( Model<dim>&,
                                   double mu,    ///< constant viscosity
                                   double coeff, ///< stabilization coefficient
                                   const char* stab_param );

    /// for 'viscosity' as element variable varying from element to element
    StabilizationParameterVisitor( Model<dim>&,
                                   double coeff,
                                   const char* stab_param );
  
    ~StabilizationParameterVisitor();
    
    virtual void Visit( Model<dim>* m );

    virtual void Visit( Element<dim>* e );
    
  private:
    const PropertyDatabase<dim>& pref;
    csmp::Index stparam_key;
    csmp::Index visc_key_;
    
    double viscosity;
    double coefficient;
    std::vector<double> segments;
    double min_segment;
    ScalarVariable sp;
};

} // namespace csmp

#endif
