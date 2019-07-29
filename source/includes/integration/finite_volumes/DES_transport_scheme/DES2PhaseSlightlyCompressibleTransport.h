#ifndef CSMP_DES_2PHASE_SLIGHTLY_COMPRESSIBLE_TRANSPORT_H
#define CSMP_DES_2PHASE_SLIGHTLY_COMPRESSIBLE_TRANSPORT_H

#include "FibonacciHeap.h"
#include "DES2PhaseTransport.h"


namespace csmp {


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
class DES2PhaseSlightlyCompressibleTransport : public DES2PhaseTransport<dim, FLOW_FUNCTIONS> {

  public:                           
    DES2PhaseSlightlyCompressibleTransport ( Model<dim>& m, 
                                             const char* target_region,
                                             bool with_gravity_forces,
                                             bool with_capillary_spreading, 
                                             double64 cfl_multiplier, 
                                             double64 PEP_multiplier,
                                             double64 relaxing_factor,
                                             bool tensor_k,
                                             FLOW_FUNCTIONS<dim>& ff);
  
    virtual ~DES2PhaseSlightlyCompressibleTransport() {}
    
    virtual void AdvectVariable_DES( double64 model_time, size_t num_threads=1 );
    virtual void AdvectVariable_TDS( double64 time_interval, size_t num_threads=1 );   
    
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;



 private:
    void InitializeVariablsAndKeys();
    virtual void InitializeEvents();
    virtual void ComputeGradients (Event<dim>* event );
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual bool Schedule(Event<dim>* nd, double64 t_end);
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t);
    virtual void Synchronize(Event<dim>* nd,double64 t_clock,double64& t_remove); 
    
    void AdvectVariable_DES_serial( double64 model_time );
    void AdvectVariable_DES_parallel ( double64 model_time, size_t num_threads ); 
    void AdvectVariable_TDS_serial( double64 time_interval );
    void AdvectVariable_TDS_parallel ( double64 time_interval, size_t num_threads );     
    
    csmp::INDEX<SCALAR,NODE> key_dsnw,key_NQV, key_sCO2_0;
    csmp::INDEX<TENSOR,ELEMENT> key_kk;
    
};
  

}//end csmp

#endif 
