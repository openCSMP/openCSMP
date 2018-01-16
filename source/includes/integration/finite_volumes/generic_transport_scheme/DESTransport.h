#ifndef CSMP_DES_TRANSPORT_H
#define CSMP_DES_TRANSPORT_H

#include "CSMP_number_types.h"
#include "Node.h"
#include "Index.h"


namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;


template<size_t dim>
class DESTransport{

  public:
    DESTransport ( Model<dim>& m, const char* target_region );
    void AdvectVariable_DES( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter, size_t num_threads );
    void AdvectVariable_DES_openmp( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter, size_t num_threads );
    void AdvectVariable_DES_serial( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter );
    void AdvectVariable_TDS( double64 time_interval, double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter );

  private:
    void initializeKeys(Model<dim>& m);
    void initializeFiniteVolumeProperties();
    void ComputeFluxBalanceAndCFL( Node<dim>* nd );
    void ComputeRateofChange( Node<dim>* nd );
    bool Schedule(Node<dim>* nd, double64 t_end, double64 cfl_multiplier);
    void Update_DES(Node<dim>* nd, double64 t_clock);
    void Update_TDS(Node<dim>* nd, double64 delta_t);
    void Synchronize(Node<dim>* nd,double64 t_clock);
    void Synchronize_openmp(Node<dim>* nd,double64 t_clock, size_t num_threads);

    Region<dim>& gref_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    std::vector<Node<dim>*>	PEPStack;
    std::vector<Node<dim>*>	Queue;
    size_t	rate_count_;
    size_t	update_count_;
    double64 T_RateOfChange_, T_Schedule_, T_SortQueue_, T_Update_, T_Synchronize_, T_RemoveFromQueue_, T_AdvectVariable_; //time recordings
    
    csmp::Index phi_key, vD_key, fv_key, PV_key, sv_key, spv_key, fA_key, ff_key, fn_key, fb_key, nsrc_key, C0_key, C1_key;
    csmp::Index update_key, rate_key, schedule_key, synchronize_key; //keys for counting
};

}//end csmp

#endif 
