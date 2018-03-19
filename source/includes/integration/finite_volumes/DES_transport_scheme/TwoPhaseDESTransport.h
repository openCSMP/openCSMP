#ifndef CSMP_TWOPHASE_DES_TRANSPORT_H
#define CSMP_TWOPHASE_DES_TRANSPORT_H

#include "CSMP_number_types.h"
#include "Node.h"
#include "Index.h"
#include "Event.h"
#include "Variables_TwoPhaseFlow.h"
#include "FlowFunctions.h"


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

  private:
    void initializeVariablsAndKeys(Model<dim>& m);
    void initializeFiniteVolumeProperties(Node<dim>* node);
    void ComputeRateofChange( Event<dim>* event );
    bool Schedule(Event<dim>* nd, double64 t_end, double64 cfl_multiplier);
    void Update_DES(Event<dim>* nd, double64 t_clock);
    void Update_TDS(Event<dim>* nd, double64 delta_t);
    void Synchronize(Event<dim>* nd,double64 t_clock);

    Region<dim>& gref_;
    FlowFunctions<dim>& flowfunctions_;
    bool with_capillary_spreading_, with_gravity_forces_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    std::vector<Event<dim>*>	PEPStack, EntireQueue, Queue;
    size_t	rate_count_;
    size_t	update_count_;
    double64 T_RateOfChange_, T_Schedule_, T_SortQueue_, T_Update_, T_Synchronize_, T_RemoveFromQueue_, T_AdvectVariable_; //time recordings
    
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


//class function for sorting event queue based on scheduled time stamps
template<size_t dim>
class sort_queue
{
    public:
    sort_queue() {}
    bool operator()(const Event<dim>* lhs, const Event<dim>* rhs)
    {
      return lhs->t_schedule() < rhs->t_schedule(); //scheduled time stamps
    }
};    

}//end csmp

#endif 
