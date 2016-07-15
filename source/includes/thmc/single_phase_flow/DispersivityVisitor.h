#ifndef SG_DISPERSIVITY_VISITOR_H
#define SG_DISPERSIVITY_VISITOR_H

#include "Visitor.h"
#include "Model.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Model;

/**
@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 2001 */
template<size_t dim>
class DispersivityVisitor : public Visitor< dim> {
  public:
    DispersivityVisitor( Model< dim>&, const char* dispersivity, const char* pore_velocity, 
                         const char* diffusivity,  const char* dispersion_long, const char* dispersion_trans );
                                                  
    DispersivityVisitor( Model< dim>&, const char* dispersivity, const char* pore_velocity, 
                          double64 diffusivity, double64 dispersion_long, double64 dispersion_trans );
                 
    virtual ~DispersivityVisitor();
    
    virtual void Visit(Element< dim>* n);   
    
  private:
    
    const PropertyDatabase<dim>&   pref;
    csmp::Index               dp_key, disp_key, at_key, al_key, vel_key;
    ScalarVariable            dp, alpha_t, alpha_l;
    double64                  v_abs, vx, vy, vz, off_diag;
    VectorVariable<dim>       vel;
    TensorVariable<dim>       disp;
    bool                      read_values;

};

} // csmp



#endif
