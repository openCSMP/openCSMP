#ifndef CSMP_HYDRO_FRACTURE_VISITOR_H
#define CSMP_HYDRO_FRACTURE_VISITOR_H

#include "Visitor.h"
#include "Model.h"
#include "AP_BoolVector.h"

namespace csmp {

/**
@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 1999
*/

template<uint32_t dim>
class HydroFractureVisitor : public Visitor<dim> {
  public:
    explicit HydroFractureVisitor( Model<dim>& sg );
    virtual ~HydroFractureVisitor();

    virtual void Visit( Element<dim>* );
    virtual void Visit( Model<dim>* );
    
    void    HydroFracturedElements( std::vector<uint32_t>& ) const;

  private:
     const PropertyDatabase<dim>&  pref;
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


















