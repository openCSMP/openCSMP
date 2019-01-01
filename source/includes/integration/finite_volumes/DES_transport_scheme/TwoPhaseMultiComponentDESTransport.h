#ifndef CSMP_TWOPHASE_MULTI_COMPONENT_DES_TRANSPORT_H
#define CSMP_TWOPHASE_MULTI_COMPONENT_DES_TRANSPORT_H

#include "FibonacciHeap.h"
#include "TwoPhaseDESTransport.h"

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
class TwoPhaseMultiComponentDESTransport : public TwoPhaseDESTransport<dim, FLOW_FUNCTIONS> {

  public:
    TwoPhaseMultiComponentDESTransport( Model<dim>& m,
                                         const char* target_region,
                                         FLOW_FUNCTIONS<dim>&,
                                         bool with_capillary_spreading,
                                         bool with_gravity_forces,
                                         bool tensor_k,
                                         double64 PEP_multiplier,
                                         double64 cfl_multiplier );
                           
    TwoPhaseMultiComponentDESTransport( Model<dim>& m,
                                         const char* target_region,
                                         FLOW_FUNCTIONS<dim>&,
                                         bool with_capillary_spreading,
                                         bool with_gravity_forces,
                                         bool tensor_k,
                                         double64 PEP_multiplier,
                                         double64 cfl_multiplier,
                                         double64 relaxing_factor );
  
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;

    virtual void InitializeVariablesAndKeys(Model<dim>& m);
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual bool Schedule(Event<dim>* nd, double64 t_end);
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t); 

  private:    
    //key for array variable storing variation rates of compositional mass of carbonic phase
    //[0] CO2 carbonic phase (supercritical/gasous CO2)
    //[1] H2O carbonic phase (evaporated water)
    csmp::INDEX<ARRAY,NODE> key_dcmCO2; 

    //key for array variable storing variation rates of compositional mass of aqueous phase
    //[0] H2O aqueous phase (liquid water)
    //[1] CO2 aqueous phase (dissolved CO2)
    //[2] NaCl aqeous phase (dissolved salt)     
    csmp::INDEX<ARRAY,NODE> key_dcmH2O;        
};
  

}//end csmp


#endif 
