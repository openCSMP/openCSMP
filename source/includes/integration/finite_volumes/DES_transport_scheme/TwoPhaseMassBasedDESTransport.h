#ifndef CSMP_TWOPHASE_MASS_BASED_DES_TRANSPORT_H
#define CSMP_TWOPHASE_MASS_BASED_DES_TRANSPORT_H

#include "FibonacciHeap.h"
#include "TwoPhaseDESTransport.h"

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
class TwoPhaseMassBasedDESTransport : public TwoPhaseDESTransport<dim, FLOW_FUNCTIONS> {

  public:
    TwoPhaseMassBasedDESTransport ( Model<dim>& m, 
                           const char* target_region, 
                           bool with_capillary_spreading, 
                           bool with_gravity_forces,
                           bool tensor_k,
                           double64 PEP_multiplier, 
                           double64 cfl_multiplier);
                           
    TwoPhaseMassBasedDESTransport ( Model<dim>& m, 
                           const char* target_region, 
                           bool with_capillary_spreading, 
                           bool with_gravity_forces,
                           bool tensor_k, 
                           double64 PEP_multiplier,
                           double64 cfl_multiplier, 
                           double64 relaxing_factor);    
  
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;

    virtual void InitializeVariablesAndKeys(Model<dim>& m);
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual bool Schedule(Event<dim>* nd, double64 t_end);
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t); 

  private:    

    csmp::INDEX<SCALAR,NODE> key_dmCO2;
        
};
  

}//end csmp


#endif 
