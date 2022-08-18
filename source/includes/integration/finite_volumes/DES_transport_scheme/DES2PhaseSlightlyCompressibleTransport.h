#ifndef CSMP_DES_2PHASE_SLIGHTLY_COMPRESSIBLE_TRANSPORT_H
#define CSMP_DES_2PHASE_SLIGHTLY_COMPRESSIBLE_TRANSPORT_H

#include "FibonacciHeap.h"
#include "DES2PhaseTransport.h"
#include "NimbleRegion.h"


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
                                             bool second_order_in_space,
                                             FLOW_FUNCTIONS<dim>& ff);
  
    virtual ~DES2PhaseSlightlyCompressibleTransport() = default;

    virtual void AdvectVariable_TDS( double time_interval, size_t num_threads );
    virtual void AdvectVariable_DES( double time_increment, double model_time, size_t num_threads);
    void SetSecondOrderInSpace(bool second_order ) {second_order_in_space_ = second_order;};
    void ReinitializeEvents();
    
    typedef ajb::detail::FibonacciHeap_Node<double,size_t> Heap_Node;



 private:
    void InitializeVariablsAndKeys();
    virtual void InitializeEvents();
    virtual void ComputeRateofChange( Event<dim>* event );
    virtual bool Schedule(Event<dim>* nd, double t_end);
    virtual void Update_DES(Event<dim>* nd, double t_clock);
    virtual void Update_TDS(Event<dim>* nd, double delta_t);
    virtual void Synchronize(Event<dim>* nd,double t_clock,double& t_remove);
    void ComputePressureGradientAndFlowVelocities (Event<dim>* event );
    void ComputeSaturationGradient (Event<dim>* event );

    void AdvectVariable_TDS_serial( double time_interval );
    void AdvectVariable_TDS_parallel ( double time_interval, size_t num_threads );
    void AdvectVariable_DES_serial( double time_increment, double model_time );
    void AdvectVariable_DES_parallel ( double time_increment, double model_time, size_t num_threads );

    //splitboundary/interface transfer functions
    bool UpdateManifold(NodeManifold<dim>* md, double t_clock, bool use_DES);
    bool UpdateManifold2(NodeManifold<dim>* md, double t_clock, double dt, bool use_DES);
    void SynchronizeManifold(NodeManifold<dim>* md, double t_clock);
    bool UpdateContactStatus(Node<dim>* masterNode, Node<dim>* slaveNode);
    double ComputeFlowPotential(Node<dim>* node);
    Element<dim>* ParentElementOfManifoldNode(Node<dim>* node, const Index& index);
    double ComputeFluxBalance( Node<dim>* nd);
    void ComputeRateofChange( std::set<Node<dim>*> nodes, bool divergence_free_correction);
    void Update_TDS( std::set<Node<dim>*> nodes, double delta_t);
    double ComputeRateofChange( Node<dim>* nd, double input_sn );

    //pressure solving functions
    void ComputeFlowPropertiesAtBaryCenter(Element<dim>& e);
    void ComputeLocalFluidPressure( double time, double dt, bool transient);
    void SolveNimbleRegionPressureFullyImplicit( NimbleRegion<dim>& computation_domain, double time_increment, double time );
    void ComputeSteadyStatePressure();
    void SolveNimbleRegionSteadyStatePressure(NimbleRegion<dim>& computation_domain);
    double TimeLevel();


    //second order in space functions
    void ComputeRateofChange_2nd_order( Event<dim>* event );
    void LimitProperty_LSMGRAD( const Element<dim>& e, const char* prop, uint32_t inside_node, uint32_t outside_node,
                                uint32_t iFacet, const double prop_inside_node, const double prop_outside_node,
                                double& limited_prop_inside_node, double& limited_prop_outside_node );
    void CalculateCenterOfMass();
    void CalculateDistanceFacetFVBary();
    double CalculateSlopeLimiter(Node<dim>* nd, const char* prop, VectorVariable<dim> grad);
    void CalculateGenericNodalGradient (Node<dim>* nd, const char* prop, VectorVariable<dim>& grad);
    /// [EL_Id] [local_facet_id] [local_node_id], distance between facet with ID local_facet_id and node with ID local_node_id
    std::vector<std::vector<std::vector<VectorVariable<dim> > > > distance_facet_FVBarycenter_;
    //second order with MINMOD limiter
    double LimitProperty( double psi_hat_c, double psi_hat_d, double psi_facet, const std::pair<double,double>& SMINMAX_upstr, double xi=2. );
    void CalculateMinMax (Node<dim>* nd, const char* prop, std::pair<double,double>& minmax);
    double NVD_Function( double xi, double U_f, double U_c );

    //access keys
    csmp::INDEX<SCALAR,NODE> key_dsnw,key_NQV, key_sCO2_0, key_sCO2_initial;
    csmp::INDEX<TENSOR,ELEMENT> key_kk;
    csmp::INDEX<SCALAR,ELEMENT> key_srCO2, key_srH2O, key_pd;
    csmp::INDEX<SCALAR,NODE> key_breakthrough, key_continuous_p, key_rhoCO2, key_rhoH2O, key_compensate, key_muH2O, key_muCO2;
    csmp::INDEX<SCALAR,MODEL> key_g;
    csmp::INDEX<VECTOR,ELEMENT> key_dip, key_vt, key_vn, key_vw;
    csmp::INDEX<SCALAR,NODE> key_status;
    csmp::INDEX<SCALAR,NODE> key_p_count;
    csmp::INDEX<VECTOR,NODE> key_mc;
    csmp::INDEX<SCALAR,NODE> key_outrange;

    //node list that stores active nodes
    std::vector<Node<dim>*> NodeList_;

    bool second_order_in_space_ = false;

    //helper class for computing residual flux at contact pairs when interface is breakthrough
    class Compensator {
	  private:
		    DES2PhaseSlightlyCompressibleTransport<dim, FLOW_FUNCTIONS>* outer_ = nullptr;
        NodeManifold<dim>* manifold_ = nullptr;
        double model_time_;
	  public:
		    explicit Compensator(){}
        Compensator(DES2PhaseSlightlyCompressibleTransport<dim, FLOW_FUNCTIONS>* outer, NodeManifold<dim>* md, double time);
        void SetManifold(NodeManifold<dim>* md) {manifold_ = md;}
        void SetModelTime(const double& time) {model_time_ = time;}
		    double operator()(const double& saturation) const;
	  }; // end inner class Compensator   
    
    
};


template <class Function>
double brent_solve(Function& func, const double x1, const double x2, const double tol, const int ITMAX);
  

}//end csmp

#endif 
