#ifndef CSMP_DES_ADVECTION_DIFFUSION_H
#define CSMP_DES_ADVECTION_DIFFUSION_H

#include "Node.h"
#include "Index.h"
#include "Event.h"
#include "FibonacciHeap.h"


namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;

template<size_t dim>
class DESAdvectionDiffusion {

  public:
    DESAdvectionDiffusion ( Model<dim>& m, const char* target_region, double cfl_multiplier, double PEP_multiplier, bool tensor_k);

    void AdvectVariable_DES( double model_time, size_t num_threads=1 );
    void AdvectVariable_TDS( double time_interval, size_t num_threads=1 );
    
    void AdvectVariable_DES_serial( double model_time );
    void AdvectVariable_DES_parallel ( double model_time, size_t num_threads ); 
    void AdvectVariable_TDS_serial( double time_interval );
    void AdvectVariable_TDS_parallel ( double time_interval, size_t num_threads );     
    
    virtual ~DESAdvectionDiffusion() {}
    
    typedef ajb::detail::FibonacciHeap_Node<double,size_t> Heap_Node;

  private:
    void InitializeVariablsAndKeys();
    void InitializeFiniteVolumeProperties();
    void InitializeEvents();
    void ComputeFluxBalanceAndCFL( Event<dim>* event);    
    void ComputeRateofChange( Event<dim>* event );
    bool Schedule(Event<dim>* event, double t_end);
    void Update_DES(Event<dim>* event, double t_clock);
    void Update_TDS(Event<dim>* event, double delta_t);
    void Synchronize(Event<dim>* event,double t_clock,double& t_remove);
    
    Model<dim>& sg_;
    Region<dim>& gref_;
    PropertyDatabase<dim>& db_;
    double upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    double CFL_multiplier_, PEP_multiplier_;
    bool tensor_k_= false;
    std::vector<Event<dim>*> PEPList, FullList;
    std::vector<Heap_Node*> HeapNodeFullList; 
    ajb::FibonacciHeap<double,size_t> EventHeap;
    size_t	rate_count_;
    size_t	update_count_;
    double T_RateOfChange_, T_Schedule_, T_InsertToHeap_, T_Update_, T_Synchronize_, T_RemoveFromHeap_, T_AdvectVariable_; //time recordings
    
    csmp::INDEX<SCALAR,NODE> key_EventIndex, key_update, key_rate, key_schedule, key_synchronize, key_fvPV, key_C0, key_FB, key_C1, key_NQC;
    csmp::INDEX<SCALAR,ELEMENT> key_phi, key_thi;
    csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA;
    csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; 
    csmp::INDEX<VECTOR,ELEMENT> key_V;
    
    // time_key - an ArrayVariable key for DES releated variables:
    // [0] current time stamp
    // [1] scheduled time stamp
    // [2] CFL time increment
    // [3] target time increment
    // [4] cumulative change of solution
    // [5] target change of solution
    csmp::INDEX<ARRAY,NODE> key_time;
    
};


}//end csmp

#endif 
