#ifndef CSMP_HYDRO_FRACTURE_VISITOR_H
#define CSMP_HYDRO_FRACTURE_VISITOR_H

#include "Visitor.h"
#include "Model.h"
#include "AP_BoolVector.h"

namespace csmp {

/**

Adapts local permeability of rock-sequence such that fluid overpressure gets dissipated.

@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 1999
*/
template<uint32_t dim>
class HydroFractureVisitor : public Visitor<dim> {
  public:
    explicit HydroFractureVisitor( Model<dim>& );

    void Visit( Element<dim>* ) override final;
    void Visit( Model<dim>* ) override final;
    
    void HydroFracturedElements( std::vector<uint32_t>& ) const;

  private:
     BoolVector               fractured;
     csmp::Index              Pe_key, S_key, K_key, V_key; 
     double                   q, dp;
     ScalarVariable           K, Kf, Pe, Pf, S;
     VectorVariable<dim>      V, dPdxy;
     bool                     over_pressured;
     const double             GRAD_LIMIT, K_LIMIT;
     std::vector<double>      pres;
};

} // csmp

#endif

/*

copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


















