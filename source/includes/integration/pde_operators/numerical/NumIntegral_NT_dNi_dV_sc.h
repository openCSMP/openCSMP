#ifndef NUM_INTEGRAL_NT_DNI_DV_SC_H
#define NUM_INTEGRAL_NT_DNI_DV_SC_H

#include "MathOperatorLHS.h"
#include "CSMP_definitions.h"

namespace csmp {

enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_dNi_dV_sc : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_NT_dNi_dV_sc( const PropertyDatabase<dim>& pref,
                              const char* oper, 
                              const char* basic, 
                              const char* test );
    
    void ComputeContribution( const CELL<dim>& ) override final;
  
    /// to chose the spatial derivate direction of interest; default is Y-axis
    void SpatialDerivative( SPATIAL_DERIVATIVE );
  
    /// to transposed the element matrix that will get accumulated; default is false
    void Transposed();
  
    NumIntegral_NT_dNi_dV_sc<dim,CELL>* clone() const override final { return new NumIntegral_NT_dNi_dV_sc<dim,CELL> (*this); }
  
  private:
    NumIntegral_NT_dNi_dV_sc();
    
    SPATIAL_DERIVATIVE     xyz_;     ///< direction of partial derivative
    bool                   transp_;  ///< transpose element matrix or not
    std::vector<double>    IPOL;     ///< interpolaton function vector
    DenseMatrix<DM_MIN>    DN, TEMP; ///< temporary matrices
};

} // csmp

#endif
