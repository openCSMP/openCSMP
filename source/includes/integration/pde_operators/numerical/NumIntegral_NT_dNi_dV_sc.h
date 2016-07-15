#ifndef NUM_INTEGRAL_NT_DNI_DV_SC_H
#define NUM_INTEGRAL_NT_DNI_DV_SC_H

#include "MathOperatorLHS.h"
#include "CSMP_definitions.h"

namespace csmp {

enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

template<size_t dim, class SIMPLEX=Element<dim> >
class NumIntegral_NT_dNi_dV_sc : public MathOperatorLHS<dim> {
  public:
    NumIntegral_NT_dNi_dV_sc( const PropertyDatabase<dim>& pref,
                        const char* oper, 
                        const char* basic, 
                        const char* test );
    
    virtual ~NumIntegral_NT_dNi_dV_sc();
    
    virtual void ComputeContribution( SIMPLEX& e );
  
    /// to chose the spatial derivate direction of interest; default is Y-axis
    void SpatialDerivative( SPATIAL_DERIVATIVE );
  
    /// to transposed the element matrix that will get accumulated; default is false
    void Transposed();
  
    virtual NumIntegral_NT_dNi_dV_sc<dim,SIMPLEX>* clone() const { return new NumIntegral_NT_dNi_dV_sc<dim,SIMPLEX> (*this); }
  
  private:
    NumIntegral_NT_dNi_dV_sc();
    
    SPATIAL_DERIVATIVE     xyz_;     ///< direction of partial derivative
    bool                   transp_;  ///< transpose element matrix or not
    std::vector<double64>  IPOL;     ///< interpolaton function vector
    DenseMatrix<DM_MIN>    DN, TEMP; ///< temporary matrices
};

} // csmp

#endif
