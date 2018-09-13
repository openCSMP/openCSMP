#ifndef CSMP_TWOPHASE_TWOCOMPONENT_DES_TRANSPORT_H
#define CSMP_TWOPHASE_TWOCOMPONENT_DES_TRANSPORT_H

#include "FibonacciHeap.h"
#include "TwoPhaseDESTransport.h"
#include "CO2H2O_FunctionsModule1.h"

namespace csmp {

template<size_t dim>
class TwoPhaseTwoComponentDESTransport : public TwoPhaseDESTransport<dim,CO2H2O_FunctionsModule1> {

  public:
    TwoPhaseTwoComponentDESTransport ( Model<dim>& m, 
                                       const char* target_region,
                                       bool with_capillary_spreading,
                                       bool with_gravity_forces,
                                       double64 PEP_multiplier, 
                                       double64 cfl_multiplier);
  
    TwoPhaseTwoComponentDESTransport ( Model<dim>& m, 
                                       const char* target_region,
                                       bool with_capillary_spreading, 
                                       bool with_gravity_forces, 
                                       double64 PEP_multiplier,
                                       double64 cfl_multiplier, 
                                       double64 relaxing_factor);    
  
    
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;

    virtual void InitializeVariablesAndKeys(Model<dim>& m);
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t); 

  private:    
    // key_components - an ArrayVariable key for the variation rates of mass (kg/s) of different transport components:
    // [0] dissolved CO2 (CO2aq)
    // [1] evaporated water (H2Og)
    csmp::INDEX<ARRAY,NODE> key_components;

    double64 upper_CO2aq_, lower_CO2aq_, upper_H2Og_, lower_H2Og_; ///< range in which the result is allowed to vary     
};
  

}//end csmp

#endif 
