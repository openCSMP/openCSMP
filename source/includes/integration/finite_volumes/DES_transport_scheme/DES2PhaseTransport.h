#ifndef CSMP_DES_2PHASE_TRANSPORT_H
#define CSMP_DES_2PHASE_TRANSPORT_H

#include "Node.h"
#include "Index.h"
#include "Event.h"
#include "FibonacciHeap.h"


namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Region;
template<uint32_t> class Model;


template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
class DES2PhaseTransport {

  public:
    DES2PhaseTransport(  Model<dim>& m,
                         const char* target_region,
                         bool with_gravity_forces,                    
                         bool with_capillary_spreading,
                         double cfl_multiplier,
                         double PEP_multiplier,                          
                         double relaxing_factor,
                         bool tensor_k,
                         FLOW_FUNCTIONS<dim>& ff );
    
    virtual ~DES2PhaseTransport() {}

    virtual void AdvectVariable_DES( double time_increment, double model_time, size_t num_threads=1 ) = 0;
    virtual void AdvectVariable_TDS( double time_interval, size_t num_threads=1 ) = 0;
    
    void SetNoFlowBoundaryCondition(bool no_flow_boundary) {no_flow_boundary_ = no_flow_boundary;}; 
   
    typedef ajb::detail::FibonacciHeap_Node<double,size_t> Heap_Node;


  protected:
    void InitializeBasicVariablsAndKeys();
    void InitializeFiniteVolumeProperties();
    virtual void InitializeEvents() = 0;
    void ResetCFLMultiplier();
    virtual void ComputeRateofChange( Event<dim>* event ) = 0;
    virtual bool Schedule(Event<dim>* event, double t_end) = 0;
    virtual void Update_DES(Event<dim>* event, double t_clock) = 0;
    virtual void Update_TDS(Event<dim>* event, double delta_t) = 0;    
    virtual void Synchronize(Event<dim>* event,double t_clock,double& t_remove) = 0;
    

    Model<dim>& sg_;
    Region<dim>& gref_;
    PropertyDatabase<dim>& db_;
    double upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    double CFL_multiplier_, PEP_multiplier_;
    bool tensor_k_= false;
    std::vector<Event<dim>> FullList;
    std::vector<Event<dim>*> PEPList;
    //std::vector<Heap_Node*> HeapNodeFullList;
    ajb::FibonacciHeap<double,size_t> EventHeap;
    size_t	rate_count_;
    size_t	update_count_;
    double T_RateOfChange_, T_Schedule_, T_InsertToHeap_, T_Update_, T_Synchronize_, T_RemoveFromHeap_, T_AdvectVariable_; //time recordings
    
    FLOW_FUNCTIONS<dim>& flowfunctions_;
    bool with_capillary_spreading_, with_gravity_forces_;
    double relaxing_factor_; //relaxing factor for cfl multiplier at areas other than saturation front
    bool no_flow_boundary_ = true;     
    
    csmp::INDEX<SCALAR,NODE> key_EventIndex, key_update, key_rate, key_schedule, key_synchronize, key_fvPV, key_sCO2, key_sH2O, key_cut, key_CFL, key_pf, key_fv;
    csmp::INDEX<SCALAR,ELEMENT> key_phi, key_thi, key_ssH2O, key_k;
    csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA;
    csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; 
    csmp::INDEX<VECTOR,ELEMENT> key_gradSn, key_gradP;
    csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_sPV;
    
    std::set<csmp::Element<dim>*>  halo_stencils_;    
    
    
    // key_time - an ArrayVariable key for DES releated variables:
    // [0] current time stamp
    // [1] scheduled time stamp
    // [2] CFL time increment
    // [3] target time increment
    // [4] cumulative change of solution
    // [5] target change of solution
    // [6] CFL multiplier
    // [7] previous time stamp
    csmp::INDEX<ARRAY,NODE> key_time;
};


}//end csmp

#endif 
