#ifndef CSMP_TWOPHASE_DES_TRANSPORT_H
#define CSMP_TWOPHASE_DES_TRANSPORT_H

#include "CSMP_number_types.h"
#include "Node.h"
#include "Index.h"
#include "Event.h"
#include "Variables_TwoPhaseFlow.h"
#include "FlowFunctions.h"
#include "FibonacciHeap.h"


namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;



template<size_t dim>
class TwoPhaseDESTransport : public variables::Variables_TwoPhaseFlow {

  public:
    TwoPhaseDESTransport ( Model<dim>& m, const char* target_region, FlowFunctions<dim>& flowfunctions, bool with_capillary_spreading, bool with_gravity_forces);
    void AdvectVariable_DES_serial( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter );
    void AdvectVariable_TDS( double64 time_interval, double64 cfl_multiplication_factor, double64 PEP_parameter );
    
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;

  private:
    void initializeVariablsAndKeys(Model<dim>& m);
    void initializeFiniteVolumeProperties(Node<dim>* node);
    void ComputeRateofChange( Event<dim>* event );
    bool Schedule(Event<dim>* nd, double64 t_end, double64 cfl_multiplier);
    void Update_DES(Event<dim>* nd, double64 t_clock);
    void Update_TDS(Event<dim>* nd, double64 delta_t);
    void Synchronize(Event<dim>* nd,double64 t_clock,double64& t_remove);

    Region<dim>& gref_;
    FlowFunctions<dim>& flowfunctions_;
    bool with_capillary_spreading_, with_gravity_forces_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    std::vector<Event<dim>*> PEPList, FullList;
    std::vector<Heap_Node*> HeapNodeFullList; 
    ajb::FibonacciHeap<double64,size_t> EventHeap;
    size_t	rate_count_;
    size_t	update_count_;
    double64 T_RateOfChange_, T_Schedule_, T_InsertToHeap_, T_Update_, T_Synchronize_, T_RemoveFromHeap_, T_AdvectVariable_; //time recordings
    bool 	first_step_;
    
    csmp::INDEX<SCALAR,NODE> key_dsnw, key_EventIndex, key_update, key_rate, key_schedule, key_synchronize;
        
    // key_time - an ArrayVariable key for DES releated variables:
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
