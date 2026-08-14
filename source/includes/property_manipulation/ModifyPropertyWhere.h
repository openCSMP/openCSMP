#ifndef __ModifyPropertyWhere_h__
#define __ModifyPropertyWhere_h__

#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class ModifyPropertyWhere : public Interrelation<dim> {
    Operand<dim>&  k;
    double                if_value, to_value;
    
  public:
    ModifyPropertyWhere( const PropertyDatabase<dim>& p, 
                         const char* prop, double where, double to );

    void Calculate() override final;
};


template<uint32_t dim>
inline void ModifyPropertyWhere<dim>::Calculate()
 {
    if ( k == if_value ) k = to_value;
 } 

} // csmp

#endif

