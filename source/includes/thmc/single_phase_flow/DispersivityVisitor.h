// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_DISPERSIVITY_VISITOR_H
#define CSMP_DISPERSIVITY_VISITOR_H

#include "Visitor.h"
#include "Model.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/**
@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 2001 */
template<uint32_t dim>
class DispersivityVisitor final : public Visitor< dim> {
  public:
    DispersivityVisitor( Model< dim>&, const char* dispersivity, const char* pore_velocity, 
                         const char* diffusivity,  const char* dispersion_long, const char* dispersion_trans );
                                                  
    DispersivityVisitor( Model< dim>&, const char* dispersivity, const char* pore_velocity, 
                          double diffusivity, double dispersion_long, double dispersion_trans );
    
    void Visit( Element<dim>* )  override final;
    void Visit( Model<dim>* )  override final {}
    
  private:
    
    const PropertyDatabase<dim>&   pref;
    csmp::Index               dp_key, disp_key, at_key, al_key, vel_key;
    ScalarVariable            dp, alpha_t, alpha_l;
    double                  v_abs, vx, vy, vz, off_diag;
    VectorVariable<dim>       vel;
    TensorVariable<dim>       disp;
    bool                      read_values;

};

} // csmp



#endif // CSMP_DISPERSIVITY_VISITOR_H
