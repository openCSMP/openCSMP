#ifndef STOKES_EINSTEIN_DIFFUSIVITY_H
#define STOKES_EINSTEIN_DIFFUSIVITY_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class StokesEinsteinDiffusivity : public Interrelation<dim> {
    Operand<dim>&  T;         /// < temperature       [oC]
    Operand<dim>&  E;         /// < dynamic viscosity [Pa s-1]
    Operand<dim>&  D;         /// < diffusivity       [m2 s-1]
    const double64    Bo,     /// < Boltzman's constant (1.380662e-23 J K-1)
                      diffus, /// < basic diffusivity
                      sixPi,
                      omega;  /// < porous medium reduction factor
    double64          R0;     /// < solute radius
    std::string       solute; /// < diffusing species
    
    ScalarVariable  eta, Tc;
    
  public:
    StokesEinsteinDiffusivity( const PropertyDatabase<dim>& p, const char* species, 
                               double64 diffusivity, double64 Omega );
    ~StokesEinsteinDiffusivity() {};
    void Calculate();
};

} // csmp

#endif

