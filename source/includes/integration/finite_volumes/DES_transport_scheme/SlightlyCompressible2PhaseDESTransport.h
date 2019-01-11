#ifndef CSMP_SLIGHTLY_COMPRESSIBLE_TWOPHASE_DES_TRANSPORT_H
#define CSMP_SLIGHTLY_COMPRESSIBLE_TWOPHASE_DES_TRANSPORT_H

#include "FibonacciHeap.h"
#include "TwoPhaseDESTransport.h"


namespace csmp {


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
class SlightlyCompressible2PhaseDESTransport : public TwoPhaseDESTransport<dim, FLOW_FUNCTIONS> {

  public:
    SlightlyCompressible2PhaseDESTransport ( Model<dim>& m, 
                                             const char* target_region,
                                             FLOW_FUNCTIONS<dim>&,
                                             bool with_capillary_spreading,
                                             bool with_gravity_forces,
                                             bool tensor_k,
                                             double64 PEP_multiplier, 
                                             double64 cfl_multiplier);
                           
    SlightlyCompressible2PhaseDESTransport ( Model<dim>& m, 
                                             const char* target_region,
                                             FLOW_FUNCTIONS<dim>&,
                                             bool with_capillary_spreading, 
                                             bool with_gravity_forces,
                                             bool tensor_k, 
                                             double64 PEP_multiplier,
                                             double64 cfl_multiplier, 
                                             double64 relaxing_factor);
  
    virtual ~SlightlyCompressible2PhaseDESTransport() {}
                           
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;


    virtual void ComputeRateofChange( Event<dim>* event );
    virtual bool Schedule(Event<dim>* nd, double64 t_end);
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t);
   
};
  

}//end csmp

#endif 
