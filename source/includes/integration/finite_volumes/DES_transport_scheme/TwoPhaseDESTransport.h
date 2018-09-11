#ifndef CSMP_TWOPHASE_DES_TRANSPORT_H
#define CSMP_TWOPHASE_DES_TRANSPORT_H

#include "CSMP_number_types.h"
#include "Node.h"
#include "Index.h"
#include "Event.h"
#include "VariableSet_CO2GeoSequestration.h"
#include "FibonacciHeap.h"


namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;



template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
class TwoPhaseDESTransport : public variables::VariableSet_CO2GeoSequestration {

  public:
    TwoPhaseDESTransport ( Model<dim>& m, 
                           const char* target_region, 
                           bool with_capillary_spreading,
                           bool with_gravity_forces,
                           double64 PEP_multiplier, 
                           double64 cfl_multiplier);
                           
    TwoPhaseDESTransport ( Model<dim>& m, 
                           const char* target_region, 
                           bool with_capillary_spreading, 
                           bool with_gravity_forces, 
                           double64 PEP_multiplier,
                           double64 cfl_multiplier, 
                           double64 relaxing_factor);
  
    virtual ~TwoPhaseDESTransport() { /* nothing to do here? */ }
                           
    double64 VolumeIntegrateScalarFiniteVolumeVariable( const Model<dim>& sg, const char* property, bool take_porosity_into_account ) const;  
    void AdvectVariable_DES( double64 model_time, size_t num_threads );
    void AdvectVariable_DES_serial( double64 model_time );
    void AdvectVariable_DES_openmp ( double64 model_time, size_t num_threads );
    void AdvectVariable_TDS( double64 time_interval );
    
    typedef ajb::detail::FibonacciHeap_Node<double64,size_t> Heap_Node;

  //private:

    virtual void initializeVariablsAndKeys(Model<dim>& m);
 // TODO: if needed   void calculatePermeabilityProjections( Region<dim>& gref );
    void initializeFiniteVolumeProperties(Event<dim>* event);
    void ResetCFLMultiplier();
    void ComputeSaturationGradient (Event<dim>* event );
    virtual void ComputeRateofChange( Event<dim>* event );
    bool Schedule(Event<dim>* nd, double64 t_end);
    virtual void Update_DES(Event<dim>* nd, double64 t_clock);
    virtual void Update_TDS(Event<dim>* nd, double64 delta_t);
    void Synchronize(Event<dim>* nd,double64 t_clock,double64& t_remove);
    

  protected:
    Region<dim>& gref_;
    bool with_capillary_spreading_, with_gravity_forces_;
    FLOW_FUNCTIONS<dim>  flowfunctions_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    std::vector<Event<dim>*> PEPList, FullList;
    std::vector<Heap_Node*> HeapNodeFullList; 
    ajb::FibonacciHeap<double64,size_t> EventHeap;
    size_t	rate_count_;
    size_t	update_count_;
    double64 T_RateOfChange_, T_Schedule_, T_InsertToHeap_, T_Update_, T_Synchronize_, T_RemoveFromHeap_, T_AdvectVariable_; //time recordings
    bool 	first_step_;
    double64 PEP_multiplier_;
    double64 CFL_multiplier_; //cfl multiplier (applied on saturation front)
    double64 relaxing_factor_; //relaxing factor for cfl multiplier at areas other than saturation front
    bool	tensor_permeability_=false;
    
    csmp::INDEX<SCALAR,NODE> key_dsnw, key_EventIndex, key_update, key_rate, key_schedule, key_synchronize, key_CFL, key_ssn, key_ssw;
    csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_grad;
        
    // key_time - an ArrayVariable key for DES releated variables:
    // [0] current time stamp
    // [1] scheduled time stamp
    // [2] CFL time increment
    // [3] target time increment
    // [4] cumulative change of solution
    // [5] target change of solution
    // [6] CFL multiplier
    csmp::INDEX<ARRAY,NODE> key_time;
     
};
  

}//end csmp

#endif 
