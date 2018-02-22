#ifndef CSMP_TWOPHASE_DES_TRANSPORT_H
#define CSMP_TWOPHASE_DES_TRANSPORT_H

#include "CSMP_number_types.h"
#include "Node.h"
#include "Index.h"
#include "Event.h"


namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;


template<size_t dim>
class TwoPhaseDESTransport{

  public:
    TwoPhaseDESTransport ( Model<dim>& m, const char* target_region);
    void AdvectVariable_DES_serial( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter );
    void AdvectVariable_TDS( double64 time_interval, double64 cfl_multiplication_factor, double64 PEP_parameter );

  private:
    void initializeVariablsAndKeys(Model<dim>& m);
    void initializeFiniteVolumeProperties();  
    void ComputeRateofChange( Event<dim>* event );
    bool Schedule(Event<dim>* nd, double64 t_end, double64 cfl_multiplier);
    void Update_DES(Event<dim>* nd, double64 t_clock);
    void Update_TDS(Event<dim>* nd, double64 delta_t);
    void Synchronize(Event<dim>* nd,double64 t_clock);
    double64 ComputeNonWettingFractionalFlow( double64 sw, double64 swr, double64 snr, double64 muw, double64 mun, double64 lambda);

    Region<dim>& gref_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    std::vector<Event<dim>*>	PEPStack, EntireQueue, Queue;
    size_t	rate_count_;
    size_t	update_count_;
    double64 T_RateOfChange_, T_Schedule_, T_SortQueue_, T_Update_, T_Synchronize_, T_RemoveFromQueue_, T_AdvectVariable_; //time recordings
    
    csmp::Index phi_key, vD_key, fv_key, PV_key, sv_key, spv_key, fA_key, ff_key, fn_key, fb_key, nsrc_key, EventIndex_key, muw_key, mun_key, swr_key, snr_key, sw_key, sn_key, dsn_key, bcp_key;
    csmp::Index update_key, rate_key, schedule_key, synchronize_key; //keys for counting
    
    // time_key - an ArrayVariable key for DES releated variables:
    // [0] current time stamp
    // [1] scheduled time stamp
    // [2] CFL time increment
    // [3] target time increment
    // [4] cumulative change of solution
    // [5] target change of solution
    csmp::Index time_key;
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
