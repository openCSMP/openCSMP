#ifndef QUARTZ_SOLUBILITY_NACL_H
#define QUARTZ_SOLUBILITY_NACL_H

#include "Interrelation.h"

namespace csmp {

/**

- aqueous solubility according to empirical formula of Fournier & Potter, GCA 46, p. 1969-73.

- correction for the presence sodium chloride according to Fournier 1983, GCA 47, p. 579-586.
*/

template<size_t dim>
class QuartzSolubilityNaCl : public Interrelation<dim> {
    Operand<dim>&    T; /// < temperature (oC)
    Operand<dim>&    R; /// < nodal density of saline fluid (kg m-3)
    Operand<dim>&    X; /// < weight fraction NaCl
    Operand<dim>&    S; /// < quartz solubility (kg silica / kg fluid)
    ScalarVariable  tC, rho, x_salt;
    double                  K, rho_e, A, B, C, 
                        m,    /// < molality of SiO2 in saline solution
                        logV, /// < decadic log of volume (g cm-3)
                        h,    /// < cation hydration number (denoting fraction of bound solvent)
                        F;    /// < weight fraction water

  public:
    QuartzSolubilityNaCl( const PropertyDatabase<dim>& p );
    ~QuartzSolubilityNaCl() {};
    void Calculate();
};


} // csmp

#endif
















