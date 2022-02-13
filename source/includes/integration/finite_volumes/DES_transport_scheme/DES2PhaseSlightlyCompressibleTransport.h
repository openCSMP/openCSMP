#ifndef CSMP_DES_2PHASE_SLIGHTLY_COMPRESSIBLE_TRANSPORT_H
#define CSMP_DES_2PHASE_SLIGHTLY_COMPRESSIBLE_TRANSPORT_H

#include "FibonacciHeap.h"
#include "DES2PhaseTransport.h"


namespace csmp {


template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
class DES2PhaseSlightlyCompressibleTransport : public DES2PhaseTransport<dim, FLOW_FUNCTIONS> {

  public:                           
    DES2PhaseSlightlyCompressibleTransport ( Model<dim>& m, 
                                             const char* target_region,
                                             bool with_gravity_forces,
                                             bool with_capillary_spreading, 
                                             double cfl_multiplier, 
                                             double PEP_multiplier,
                                             double relaxing_factor,
                                             bool tensor_k,
                                             FLOW_FUNCTIONS<dim>& ff);
  
    virtual ~DES2PhaseSlightlyCompressibleTransport() {}
    
    virtual void AdvectVariable_DES( double model_time, size_t num_threads=1 );
    virtual void AdvectVariable_TDS( double time_interval, size_t num_threads=1 );   
    
    typedef ajb::detail::FibonacciHeap_Node<double,size_t> Heap_Node;



 private:
    void InitializeVariablsAndKeys();
    virtual void InitializeEvents();
    virtual void ComputeGradients (Event<dim>* event );
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual bool Schedule(Event<dim>* nd, double t_end);
    virtual void Update_DES(Event<dim>* nd, double t_clock);
    virtual void Update_TDS(Event<dim>* nd, double delta_t);
    virtual void Synchronize(Event<dim>* nd,double t_clock,double& t_remove); 
    
    void AdvectVariable_DES_serial( double model_time );
    void AdvectVariable_DES_parallel ( double model_time, size_t num_threads ); 
    void AdvectVariable_TDS_serial( double time_interval );
    void AdvectVariable_TDS_parallel ( double time_interval, size_t num_threads );     
    
    csmp::INDEX<SCALAR,NODE> key_dsnw,key_NQV, key_sCO2_0;
    csmp::INDEX<TENSOR,ELEMENT> key_kk;
    
};
  

}//end csmp

#endif 
