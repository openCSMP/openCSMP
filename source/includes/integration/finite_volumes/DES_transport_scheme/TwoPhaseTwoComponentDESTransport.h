#ifndef CSMP_TWOPHASE_TWOCOMPONENT_DES_TRANSPORT_H
#define CSMP_TWOPHASE_TWOCOMPONENT_DES_TRANSPORT_H

#include "FibonacciHeap.h"
#include "TwoPhaseDESTransport.h"


namespace csmp {

template<size_t dim>
class TwoPhaseTwoComponentDESTransport : public TwoPhaseDESTransport<dim> {

  public:
    TwoPhaseTwoComponentDESTransport ( Model<dim>& m, 
                           const char* target_region, 
                           FlowFunctions<dim>& flowfunctions, 
                           bool with_capillary_spreading, 
                           bool with_gravity_forces,
                           double64 PEP_multiplier, 
                           double64 cfl_multiplier);
                           
    TwoPhaseTwoComponentDESTransport ( Model<dim>& m, 
                           const char* target_region, 
                           FlowFunctions<dim>& flowfunctions, 
                           bool with_capillary_spreading, 
                           bool with_gravity_forces, 
                           double64 PEP_multiplier,
                           double64 cfl_multiplier, 
                           double64 relaxing_factor);    
                           
    
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;

    virtual void initializeVariablsAndKeys(Model<dim>& m);
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t); 

  private:    
    // key_components - an ArrayVariable key for the variation rates of mass (kg/s) of different transport components:
    // [0] CO2 in carbonic phase (YCO2)
    // [1] H2O in carbonic phase (YH2O)
    // [2] CO2 in aqueous phase  (XCO2)
    // [3] H2O in aqueous phase  (XH2O)
    csmp::INDEX<ARRAY,NODE> key_components;
     
};
  

}//end csmp

#endif 
