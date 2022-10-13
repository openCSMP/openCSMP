#include "DES2PhaseSlightlyCompressibleTransport.h"
#include "NimbleRegion.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#include "compareFloats.h"
#include "FlowFunctionsModule.h"
#include "TransientDiffusor.h"
#include "NumIntegral_NT_op_N_dS.h"
#include "NumIntegral_NT_lhsop_N_dV.h"

#include "VTU_Interface.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#endif

#if defined(OPENMP)
#include "omp.h"
#endif

using namespace std;

namespace csmp {


template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::DES2PhaseSlightlyCompressibleTransport(  Model<dim>& m,
                                                                                                     const char* target_region,
                                                                                                     bool with_gravity_forces,
                                                                                                     bool with_capillary_spreading,
                                                                                                     double cfl_multiplier,
                                                                                                     double PEP_multiplier,
                                                                                                     double relaxing_factor,
                                                                                                     bool tensor_k,
                                                                                                     bool second_order_in_space,
                                                                                                     FLOW_FUNCTIONS<dim>& ff )
    : DES2PhaseTransport<dim,FLOW_FUNCTIONS>(m,target_region,with_gravity_forces,with_capillary_spreading,cfl_multiplier,PEP_multiplier,relaxing_factor,tensor_k,ff),
      second_order_in_space_(second_order_in_space)
{
    InitializeVariablsAndKeys();
    this->InitializeEvents();
    cout<<"\nDES2PhaseSlightlyCompressibleTransport constructed"<<endl;
    cout<<"with gravity forces = "<<this->with_gravity_forces_<<" (0=false, 1=true)"<<endl; 
    cout<<"with capillary spreading = "<<this->with_capillary_spreading_<<" (0=false, 1=true)"<<endl; 
    cout<<"CFL multiplier = "<<this->CFL_multiplier_<<endl;
    cout<<"PEP multiplier = "<<this->PEP_multiplier_<<endl;
    cout<<"relaxing factor = "<<this->relaxing_factor_<<endl;
    cout<<"tensor permeability = "<<tensor_k<<" (0=false, 1=true)\n"<<endl;    
    cout<<"second order in space = "<<second_order_in_space_<<" (0=false, 1=true)\n"<<endl;
} // end constructor 




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::InitializeVariablsAndKeys()
{
    //creating new variables if not defined yet
    if(!this->db_.IsDefined("variation rate nonwetting phase")) this->sg_.CreateProperty( "variation rate nonwetting phase", "vrnwp", "m3/(m3.s)", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);
    this->sg_.Region("Model").InputPropertyValue( "variation rate nonwetting phase", makeScalar(PLAIN,0), COMPLETE);
    if(!this->db_.IsDefined("compensation flux rate nonwetting phase")) this->sg_.CreateProperty( "compensation flux rate nonwetting phase", "cfrnwp", "m3/(m3.s)", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);
    this->sg_.Region("Model").InputPropertyValue( "compensation flux rate nonwetting phase", makeScalar(PLAIN,0), COMPLETE);
    if(!this->db_.IsDefined("old saturation carbonic phase")) this->sg_.CreateProperty( "old saturation carbonic phase", "sCO2_0", "m3/m3", SCALAR, NODE, 1, 0 ,1);
    this->sg_.CopyReplace( "saturation carbonic phase", "old saturation carbonic phase" );
    if(!this->db_.IsDefined("tensor permeability")) this->sg_.CreateProperty( "tensor permeability", "kk", "m2", TENSOR, ELEMENT, 3, 1E-21, 1.0e-5);
    if(!this->db_.IsDefined("nodal fluid volume source")) {
        this->sg_.CreateProperty( "nodal fluid volume source", "nfvq", "m3", SCALAR, NODE, 1, -1.00E+01, 1.00E+01);
        this->sg_.Region("Model").InputPropertyValue( "nodal fluid volume source", makeScalar(PLAIN,0), COMPLETE);
    }
    if(!this->db_.IsDefined("residual saturation carbonic phase")) this->sg_.CreateProperty( "residual saturation carbonic phase", "sCO2r", "m3/m3", SCALAR, ELEMENT, 1, 0., 1.);
    if(!this->db_.IsDefined("residual saturation aqueous phase")) this->sg_.CreateProperty( "residual saturation aqueous phase", "swr", "m3/m3", SCALAR, ELEMENT, 1, 0., 1.);
    if(!this->db_.IsDefined("breakthrough status")) this->sg_.CreateProperty( "breakthrough status", "bs", "none", SCALAR, NODE, 1, 0, 1);
    this->sg_.Region("Model").InputPropertyValue( "breakthrough status", makeScalar(PLAIN,0), COMPLETE);
    if(!this->db_.IsDefined("pressure continuity status")) this->sg_.CreateProperty( "pressure continuity status", "pfcs", "none", SCALAR, NODE, 1, 0, 1);
    this->sg_.Region("Model").InputPropertyValue( "pressure continuity status", makeScalar(PLAIN,0), COMPLETE);
    if(!this->db_.IsDefined("entry pressure")) this->sg_.CreateProperty( "entry pressure", "pd", "Pa", SCALAR, ELEMENT, 1, 0., 50000000.);
    if(!this->db_.IsDefined("acceleration gravity")) {
        this->sg_.CreateProperty( "acceleration gravity", "acc", "m/s2", SCALAR, MODEL, 1, 9.76, 9.83);
        this->sg_.Region("Model").InputPropertyValue( "acceleration gravity", makeScalar(PLAIN,9.8061), COMPLETE);
    }    
    if(!this->db_.IsDefined("density carbonic phase")) this->sg_.CreateProperty( "density carbonic phase", "rhoCO2", "kg/m3", SCALAR, NODE, 1, 50., 1000.);
    if(!this->db_.IsDefined("density aqueous phase")) this->sg_.CreateProperty( "density aqueous phase", "rhow", "kg/m3", SCALAR, NODE, 1, 500., 1500.);
    if(!this->db_.IsDefined("dip vector")) this->sg_.CreateProperty( "dip vector", "dip", "none", VECTOR, ELEMENT, 3, -1., 1.);
    if(!this->db_.IsDefined("total velocity")) this->sg_.CreateProperty( "total velocity", "vt", "m/s", VECTOR, ELEMENT, 3, -1.0e8, 1.0e8);
    if(!this->db_.IsDefined("viscosity carbonic phase")) this->sg_.CreateProperty( "viscosity carbonic phase", "muCO2", "Pa.s", SCALAR, NODE, 1, 1E-05, 0.1);
    if(!this->db_.IsDefined("viscosity aqueous phase")) this->sg_.CreateProperty( "viscosity aqueous phase", "muw", "Pa.s", SCALAR, NODE, 1, 1E-05, 0.0016);
    if(!this->db_.IsDefined("velocity carbonic phase")) this->sg_.CreateProperty( "velocity carbonic phase", "vCO2", "m/s", VECTOR, ELEMENT, 3, -1.0e8, 1.0e8);
    if(!this->db_.IsDefined("velocity aqueous phase")) this->sg_.CreateProperty( "velocity aqueous phase", "vw", "m/s", VECTOR, ELEMENT, 3, -1.0e8, 1.0e8);
    if(!this->db_.IsDefined("fluid pressure status")) this->sg_.CreateProperty( "fluid pressure status", "pfs", "none", SCALAR, NODE, 1, 0 ,8);
    if(!this->db_.IsDefined("local pressure solving count")) this->sg_.CreateProperty( "local pressure solving count", "lpfsc", "uint", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    this->sg_.Region("Model").InputPropertyValue( "local pressure solving count", makeScalar(PLAIN,0), COMPLETE);
    if(!this->db_.IsDefined("mass center")) this->sg_.CreateProperty( "mass center", "msc", "none", VECTOR, NODE, 3, -1.00E+10 ,1.00E+10);
    if(!this->db_.IsDefined("out range value count")) this->sg_.CreateProperty( "out range value count", "orvc", "uint", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    this->sg_.Region("Model").InputPropertyValue( "out range value count", makeScalar(PLAIN,0), COMPLETE);
    if(!this->db_.IsDefined("initial saturation carbonic phase")) this->sg_.CreateProperty( "initial saturation carbonic phase", "sCO2i", "m3/m3", SCALAR, NODE, 1, 0 ,1);
    this->sg_.CopyReplace( "saturation carbonic phase", "initial saturation carbonic phase" );

    //assigning keys 
    key_dsnw = INDEX<SCALAR,NODE> ( this->db_.StorageKey("variation rate nonwetting phase") );
    key_sCO2_0 = INDEX<SCALAR,NODE>( this->db_.StorageKey("old saturation carbonic phase") );
    key_NQV = INDEX<SCALAR,NODE>( this->db_.StorageKey("nodal fluid volume source") );
    key_kk = INDEX<TENSOR,ELEMENT>( this->db_.StorageKey("tensor permeability") );
    key_srCO2 = INDEX<SCALAR,ELEMENT>( this->db_.StorageKey("residual saturation carbonic phase") ); 
    key_srH2O = INDEX<SCALAR,ELEMENT>( this->db_.StorageKey("residual saturation aqueous phase") ); 
    key_breakthrough = INDEX<SCALAR,NODE>( this->db_.StorageKey("breakthrough status") );
    key_pd = INDEX<SCALAR,ELEMENT>( this->db_.StorageKey("entry pressure") ); 
    key_g = INDEX<SCALAR,MODEL>( this->db_.StorageKey("acceleration gravity") ); 
    key_rhoH2O = INDEX<SCALAR,NODE>( this->db_.StorageKey("density aqueous phase") );
    key_rhoCO2 = INDEX<SCALAR,NODE>( this->db_.StorageKey("density carbonic phase") );
    key_compensate = INDEX<SCALAR,NODE> ( this->db_.StorageKey("compensation flux rate nonwetting phase") );
    key_dip = INDEX<VECTOR,ELEMENT>( this->db_.StorageKey("dip vector") );
    key_vt = INDEX<VECTOR,ELEMENT>( this->db_.StorageKey("total velocity") );
    key_muH2O = INDEX<SCALAR,NODE>( this->db_.StorageKey("viscosity aqueous phase") );
    key_muCO2 = INDEX<SCALAR,NODE>( this->db_.StorageKey("viscosity carbonic phase") ); 	
    key_vn = INDEX<VECTOR,ELEMENT>( this->db_.StorageKey("velocity carbonic phase") );
    key_vw = INDEX<VECTOR,ELEMENT>( this->db_.StorageKey("velocity aqueous phase") );
    key_status = INDEX<SCALAR,NODE> ( this->db_.StorageKey("fluid pressure status") );
    key_p_count = INDEX<SCALAR,NODE>( this->db_.StorageKey("local pressure solving count") );
    key_mc = INDEX<VECTOR,NODE>( this->db_.StorageKey("mass center") );
    key_outrange = INDEX<SCALAR,NODE>( this->db_.StorageKey("out range value count") );
    key_sCO2_initial = INDEX<SCALAR,NODE>( this->db_.StorageKey("initial saturation carbonic phase") );
    
    //checking keys 
    if ( key_dsnw.place != NODE || key_dsnw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'variation rate nonwetting phase' variable must be SCALAR and placed on NODE"  );    
    if ( key_sCO2_0.place != NODE || key_sCO2_0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'old saturation carbonic phase' variable must be SCALAR and placed on NODE"  );    
    if ( key_NQV.place != NODE || key_NQV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );      
    if ( key_kk.place != ELEMENT || key_kk.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'tensor permeability' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_srCO2.place != ELEMENT || key_srCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'residual saturation carbonic phase' variable must be SCALAR  and placed on ELEMENT"  );  
    if ( key_srH2O.place != ELEMENT || key_srH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'residual saturation aqueous phase' variable must be SCALAR  and placed on ELEMENT"  );          
    if ( key_breakthrough.place != NODE || key_breakthrough.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'breakthrough status' variable must be SCALAR and placed on NODE"  );
    if ( key_pd.place != ELEMENT || key_pd.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'entry pressure' variable must be SCALAR  and placed on ELEMENT"  ); 
    if ( key_g.place != MODEL || key_g.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'acceleration gravity' variable must be SCALAR and placed on MODEL"  );   
    if ( key_rhoH2O.place != NODE || key_rhoH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'density aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoCO2.place != NODE || key_rhoCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'density carbonic phase' variable must be SCALAR and placed on NODE"  );     
    if ( key_compensate.place != NODE || key_compensate.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'compensation flux rate nonwetting phase' variable must be SCALAR and placed on NODE"  );
    if ( key_dip.place != ELEMENT || key_dip.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'dip vector' variable must be VECTOR and placed on ELEMENT"  );                                
    if ( key_vt.place != ELEMENT || key_vt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'total velocity' variable must be VECTOR and placed on ELEMENT"  ); 
    if ( key_muH2O.place != NODE || key_muH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'viscosity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_muCO2.place != NODE || key_muCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'viscosity carbonic phase' variable must be SCALAR and placed on NODE"  );		
    if ( key_vn.place != ELEMENT || key_vn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'velocity carbonic phase' variable must be VECTOR and placed on ELEMENT"  );    
    if ( key_vw.place != ELEMENT || key_vw.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'velocity aqueous phase' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_status.place != NODE || key_status.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'fluid pressure status' variable must be SCALAR and placed on NODE"  );   
    if ( key_p_count.place != NODE || key_p_count.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'local pressure solving count' variable must be SCALAR and placed on NODE"  );     
    if ( key_mc.place != NODE || key_mc.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'mass center' variable must be VECTOR and placed on NODE"  );     
    if ( key_outrange.place != NODE || key_outrange.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
        "The 'out range value count' variable must be SCALAR and placed on NODE"  );
    if ( key_sCO2_initial.place != NODE || key_sCO2_initial.type != SCALAR )
    throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport::InitializeVariablesAndKeys:",
                           "The 'initial saturation carbonic phase' variable must be SCALAR and placed on NODE"  );

}




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::InitializeEvents()
{   
    if(second_order_in_space_) {
        CalculateCenterOfMass();
        CalculateDistanceFacetFVBary();
    }

    //create events for all nodes and add them to event lists
    size_t index = 0;
    size_t dirich_count = 0;

    const auto nodes_end(this->gref_.NodesEnd());
    for ( auto nit=this->gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        if((*nit)->Status(  this->key_sCO2 ) != DIRICH) {
            (*nit)->Store( this->key_EventIndex, makeScalar( (*nit)->Status(this->key_EventIndex), index) );//event index 
            auto* event = new Event<dim>(*nit);
            this->PEPList.push_back(event);
            event->inPEPStack(true);
            ComputePressureGradientAndFlowVelocities(event);
            if(this->with_capillary_spreading_) ComputeSaturationGradient (event);
            if(second_order_in_space_) ComputeRateofChange_2nd_order(event);
            else ComputeRateofChange(event);
            Schedule(event, 0.);
            event->valid(false);
            auto* heap_node = new Heap_Node(event->t_schedule(),index);
            this->HeapNodeFullList.push_back(heap_node);
            this->FullList.push_back(event);
            event->inQueue(false);
            index++;
        } else {
            dirich_count++;
        }
    }

    cout<<this->FullList.size()<<" events created for all nodes, excluding "<<dirich_count<<" DIRICH nodes"<<endl;

    //update all contact status
    if(this->sg_.Mesh().NodeManifolds() > 0) {
        for(auto mit = this->sg_.Mesh().NodeManifoldsBegin();mit!=this->sg_.Mesh().NodeManifoldsEnd();mit++) {
            auto md = (*mit);
            auto master_node = md.N(0);
            for(size_t n=1;n<md.Branches();n++){
                auto slave_node = md.N(n);
                UpdateContactStatus(master_node, slave_node);
            }
        }
    }

}  



template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ReinitializeEvents()
{   
    //bring up all nodes to current state?

    //clean up all containers
    this->PEPList.clear();
    this->HeapNodeFullList.clear();
    this->FullList.clear();
    this->EventHeap.clear();

    //create events for all nodes and add them to event lists
    size_t index = 0;
    size_t dirich_count = 0;

    const auto nodes_end(this->gref_.NodesEnd());
    for ( auto nit=this->gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        if((*nit)->Status(  this->key_sCO2 ) != DIRICH) {
            (*nit)->Store( this->key_EventIndex, makeScalar( (*nit)->Status(this->key_EventIndex), index) );//event index 
            auto* event = new Event<dim>(*nit);
            this->PEPList.push_back(event);
            event->inPEPStack(true);
            ComputePressureGradientAndFlowVelocities(event);
            if(this->with_capillary_spreading_) ComputeSaturationGradient (event);
            if(second_order_in_space_) ComputeRateofChange_2nd_order(event);
            else ComputeRateofChange(event);
            Schedule(event, 0.);
            event->valid(false);
            auto* heap_node = new Heap_Node(event->t_schedule(),index);
            this->HeapNodeFullList.push_back(heap_node);
            this->FullList.push_back(event);
            event->inQueue(false);
            index++;
        } else {
            dirich_count++;
        }
    }

    cout<<this->FullList.size()<<" events created for all nodes, excluding "<<dirich_count<<" DIRICH nodes"<<endl;

    //update all contact status
    if(this->sg_.Mesh().NodeManifolds() > 0){
        for(auto mit = this->sg_.Mesh().NodeManifoldsBegin();mit!=this->sg_.Mesh().NodeManifoldsEnd();mit++) {
            auto md = (*mit);
            auto master_node = md.N(0);
            for(size_t n=1;n<md.Branches();n++){
                auto slave_node = md.N(n);
                UpdateContactStatus(master_node, slave_node);
            }
        }
    }
}  



// TODO: why is the saturation gradient not computed using the neighboring nodes?
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeSaturationGradient (Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != nullptr );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  
    //check if node is truncated by domain boundary
    int64_t truncated_node = static_cast<int64_t>(nd->Read(this->key_cut));
  
    const auto parent_elements(nd->Parents());
    for ( auto i{0U}; i<parent_elements; ++i )
    {
        Element<dim>* const eptr = nd->Parent(i);
        if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) { //ignore if parent element located outside domain
            continue;
        } else {
            DenseMatrix<DM_MIN> DN;
            eptr->dN_AtBaryCenter(DN);
            VectorVariable<dim> snw_gradient;
            snw_gradient = 0.;
            for ( auto j{0U}; j<eptr->Nodes(); j++ )
            {
                const double sn = eptr->N(j)->Read(this->key_sCO2);
                for ( auto k{0U}; k<dim; k++ ) snw_gradient(k) += DN(k,j) * sn;
            }
            eptr->Store(this->key_gradSn, snw_gradient);   
        }
    }
}




//this function shall be used when "total mobility permeability product" and "gravity term" have been computed and stored previously
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputePressureGradientAndFlowVelocities (Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != nullptr );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH);
// NOT USED    constexpr uint32_t  v( (dim==1u) ? 0u : 1u );
  
    //check if node is truncated by domain boundary
    int truncated_node = static_cast<int>(nd->Read(this->key_cut));
  
    const auto parent_elements(nd->Parents());
    for ( auto i{0U}; i<parent_elements; ++i )
    {
        Element<dim>* const eptr = nd->Parent(i);
        
        if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) { //ignore if parent element located outside domain
            continue;
        } else {
            DenseMatrix<DM_MIN> DN;
            eptr->dN_AtBaryCenter(DN);
            VectorVariable<dim> p_gradient;
            p_gradient = 0.;
            for ( auto j{0U}; j<eptr->Nodes(); j++ )
            {
                const double p = eptr->N(j)->Read(this->key_pf);
                for ( auto k{0U}; k<dim; k++ ) p_gradient(k) += -DN(k,j) * p;
            }
            eptr->Store(this->key_gradP, p_gradient);   

            VectorVariable<dim> e_vt(ANY, 0.);
            if(!this->tensor_k_) { //scalar permeability
                static Index  mobt_key(this->db_.StorageKey("total mobility permeability product"));
                double e_lt = eptr->Read(mobt_key); //total mobility permeability product
                if(isnan(e_lt)) 
                    throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputePressureGradientAndFlowVelocities:",
                    "total mobility permeability product has not been computed yet" );

                e_vt(0) = e_lt * p_gradient(0); //(krn/mun+krw/muw)*k*thi * gradP
                if ( dim != 1U ) e_vt(1) = e_lt * p_gradient(1);
                if ( dim == 3U ) e_vt(2) = e_lt * p_gradient(2); 

                if( this->with_gravity_forces_ ) {
                    static Index  gt_key(this->db_.StorageKey("gravity term"));
                    VectorVariable<dim> gravity;
                    eptr->Read( gt_key, gravity ); //(krn/mun*rhon+krw/muw*rhow)*k*thi*g
                    if(isnan(gravity(0))) {
                        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputePressureGradientAndFlowVelocities:", "gravity term has not been computed yet" );
                    }
                    e_vt += gravity;
                }  
                eptr->Store( this->key_vt, e_vt );
          
            } else { //tensor k
                static Index  LT_key(this->db_.StorageKey("tensor total mobility permeability product")); 
                TensorVariable<dim> LT;
                eptr->Read( LT_key, LT );
                if(isnan(LT(0,0)))
                    throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputePressureGradientAndFlowVelocities:",
                    "tensor total mobility permeability product has not been computed yet" );

                e_vt = LT * p_gradient;

                if( this->with_gravity_forces_ ) {
                    static Index  gt_key(this->db_.StorageKey("gravity term"));
                    VectorVariable<dim> gravity;
                    eptr->Read( gt_key, gravity ); //(krn/mun*rhon+krw/muw*rhow)*k*thi*g
                    if(isnan(gravity(0))) {
                        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputePressureGradientAndFlowVelocities:", "gravity term has not been computed yet" );
                    }
                    e_vt += gravity;
                }

                eptr->Store( this->key_vt, e_vt );
            }
        } 
    }
}




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != nullptr );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH){
      double PV = nd->Read(this->key_fvPV); //pore volume (m3)
      if(PV < numeric_limits<double>::epsilon())
          throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange:",
                                 "FV pore volume seems to be zero, has it been initialised yet?" );
      else
      {
          this->rate_count_++;//recording
          nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );

          //some variables to use
          double flux_balance(0.), outflow(0.), tot_inflow(0.), tot_outflow(0.), carb_accumulation(0.), aq_accumulation(0.);
          const size_t v( (dim==1u) ? 0u : 1u );
          VectorVariable<dim> facetNrml, gravity, gradP;
          vector<double> IPOL, NRST;

          const auto node_parent_elements(nd->Parents());
          double cfl_multiplier = this->CFL_multiplier_*this->relaxing_factor_; //default value
          long truncated_node = static_cast<long>(nd->Read(this->key_cut));//check if node is truncated by domain boundary

          for ( auto t{0U}; t<node_parent_elements; t++ )
          {
              Element<dim>* const eptr(nd->Parent(t));
              assert( eptr != nullptr );
      
              if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) //ignore if parent element located outside domain
                  continue;
        
              //element state and properties
              const auto nodes(eptr->Nodes());
              eptr->N_AtBaryCenter( IPOL );
              double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.);
              for ( auto i{0U}; i<nodes; ++i ) {
                  e_sw += IPOL[i] * eptr->N(i)->Read( this->key_sH2O );
		              e_muw += IPOL[i] * eptr->N(i)->Read( this->key_muH2O );
		              e_mun += IPOL[i] * eptr->N(i)->Read( this->key_muCO2 );
		              e_rhow += IPOL[i] * eptr->N(i)->Read( this->key_rhoH2O );
		              e_rhon += IPOL[i] * eptr->N(i)->Read( this->key_rhoCO2 );
                  ipol_sum += IPOL[i];
              }
              if(ipol_sum > 0.) {
		              e_sw *= 1. / ipol_sum;
		              e_muw *= 1. / ipol_sum;
		              e_mun *= 1. / ipol_sum;
		              e_rhow *= 1. / ipol_sum;
		              e_rhon *= 1. / ipol_sum;
              }
              e_sn = 1.0 - e_sw;

              const double e_swr = eptr->Read(key_srH2O);
              const double e_snr = eptr->Read(key_srCO2);

              eptr->Read(this->key_gradP, gradP); //pressure gradient
              double thickness = eptr->Read(this->key_thi); //thickness
              double e_k = eptr->Read(this->key_k); //permeability

              //read in total velocity computed previously
              VectorVariable<dim> e_vt(ANY, 0.);
              eptr->Read(this->key_vt, e_vt);

              double carb_inflow(0.);
              const auto pnid(nd->ParentNodeNumber(t));
              const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
              for ( auto i{0U}; i<sector_facets; i++ )
              {
                  const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                  const auto inside_node(eptr->FV()->InsideNode(iFacet));
                  const auto outside_node(eptr->FV()->OutsideNode(iFacet));
            
                  eptr->Read( iFacet, 0U,  this->key_fn, facetNrml );
                  const double facetArea = eptr->Read( iFacet, 0U,  this->key_fA );
                  const double sign = ( pnid == inside_node ) ? 1. : -1.;

                  //compute facet fluid flux
                  const double  vD_n = e_vt.DotProduct(facetNrml); //m/s
                  //update flux balance
                  double f_total = sign * vD_n * facetArea;
                  if(f_total != 0.) {
                      flux_balance += f_total;
                      //if(f_total < 0.) outflow -= f_total;
                      if(f_total > 0.) outflow += f_total;
                      //inflow += f_total;
                  }

                  //compute facet variables
                  eptr->N_AtFacetIntegrationPoint( iFacet, 0U, NRST );
                  double f_sw(0.), f_muw (0.), f_mun(0.), f_rhow(0.), f_rhon(0.);
                  for ( auto x{0U}; x<eptr->Nodes(); x++ ) {
                      f_sw += NRST[x] * eptr->N(x)->Read( this->key_sH2O );
                      f_muw += NRST[x] * eptr->N(x)->Read( this->key_muH2O );
                      f_mun += NRST[x] * eptr->N(x)->Read( this->key_muCO2 );
                      f_rhow += NRST[x] * eptr->N(x)->Read( this->key_rhoH2O );
                      f_rhon += NRST[x] * eptr->N(x)->Read( this->key_rhoCO2 );
                  }

                  //compute upstream directions at facet integration point
                  double vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
                  double vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
                  if( this->with_gravity_forces_ ) {
                      eptr->Read( key_dip, gravity );
                      if(isnan(gravity(v))) { //dip vector has not been initialised
                          if constexpr (dim==1u) {gravity(0u) = -1.;}
                          else if constexpr (dim==2u) {gravity(0u) = 0.; gravity(1u) = -1.;}
                          else {gravity(0u) = 0.; gravity(1u) = -1.; gravity(2u) = 0.;}
                          eptr->Store( key_dip, gravity );
                      }
               
                      if(!this->tensor_k_) { //scalar k
                          gravity *=  e_k * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
                      } else { //tensor k
                          TensorVariable<dim> kk;
                          eptr->Read( key_kk, kk );
                          VectorVariable<dim> kV;
                          kV = kk * gravity;
                          double kV_magnitude = kV.Length();
                          gravity *= kV_magnitude * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
                      }
                      double gravity_nrml = gravity.DotProduct(facetNrml);
                      vn_gravity_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * gravity_nrml;
                      vw_gravity_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * gravity_nrml;
                  }
            
                  if(this->with_capillary_spreading_) {
                      VectorVariable<dim> grad;
                      eptr->Read(this->key_gradSn, grad);
                      if(!this->tensor_k_) { //scalar k
                          grad *= ( -1.0* e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                      } else { //tensor k
                          TensorVariable<dim> kk;
                          eptr->Read( key_kk, kk );
                          grad = kk * grad;
                          grad *= ( -1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                      }
                      double grad_nrml = grad.DotProduct(facetNrml);
                      vn_capillary_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * grad_nrml;
                      vw_capillary_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * grad_nrml;
                  }
                       
                  double vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
                  double vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity - vw_capillary_component_of_velocity;

                  double total_f(0.), lw(0.), ln(0.);
                  //aqueous phase determination
                  if(vw_at_facet_int_point != 0.0) {
                      const auto upstream_node = (vw_at_facet_int_point  > 0.) ? inside_node : outside_node;
                      double sw = eptr->N(upstream_node)->Read(this->key_sH2O);
                      if(sw > e_swr) lw = this->flowfunctions_.krw_at(eptr, sw) / f_muw;
                  } else { //vw_at_facet_int_point == 0.0
                      double sw_inside = eptr->N(inside_node)->Read(this->key_sH2O);
                      double lw_inside (0.);
                      if(sw_inside > e_swr) lw_inside = this->flowfunctions_.krw_at(eptr, sw_inside) / f_muw;

                      double sw_outside = eptr->N(outside_node)->Read(this->key_sH2O);
                      double lw_outside (0.);
                      if(sw_outside > e_swr) lw_outside = this->flowfunctions_.krw_at(eptr, sw_outside) / f_muw;

                      lw = (lw_inside + lw_outside) * 0.5;
                  }
            
                  //carbonic phase determination
                  if(vn_at_facet_int_point != 0.0) {
                      const auto upstream_node = (vn_at_facet_int_point  > 0.) ? inside_node : outside_node;
                      double sw = eptr->N(upstream_node)->Read(this->key_sH2O);
                      double sn = eptr->N(upstream_node)->Read(this->key_sCO2);
                      if(sn > e_snr) ln = this->flowfunctions_.krn_at(eptr, sw) / f_mun;
                  } else { //vn_at_facet_int_point == 0.0
                      double sw_inside = eptr->N(inside_node)->Read(this->key_sH2O);
                      double sn_inside = eptr->N(inside_node)->Read(this->key_sCO2);
                      double ln_inside (0.);
                      if(sn_inside > e_snr) ln_inside = this->flowfunctions_.krn_at(eptr, sw_inside) / f_mun;

                      double sw_outside = eptr->N(outside_node)->Read(this->key_sH2O);
                      double sn_outside = eptr->N(outside_node)->Read(this->key_sCO2);
                      double ln_outside (0.);
                      if(sn_outside > e_snr) ln_outside = this->flowfunctions_.krn_at(eptr, sw_outside) / f_mun;

                      ln = (ln_inside + ln_outside) * 0.5;
                  }

                  //fractional flows determination
                  double total_mobility = lw + ln;
                  double upstream_fn=(total_mobility!=0.0 ? ln/total_mobility : 0.0);
                  double upstream_lambda_overbar=(total_mobility!=0.0 ? (ln*lw)/total_mobility : 0.0);
                  //viscous flux
                  double viscous_velocity_component = vD_n * upstream_fn;
                  //gravitational flux
                  double gravity_velocity_component(0.);
                  if( this->with_gravity_forces_ ) {
                      double gravity_nrml = gravity.DotProduct(facetNrml);
                      gravity_velocity_component = upstream_lambda_overbar * gravity_nrml;
                  }
                  //capillary flux
                  double capillary_velocity_component(0.);
                  if( this->with_capillary_spreading_ ) {
                      VectorVariable<dim> grad;
                      eptr->Read(this->key_gradSn, grad);
                      if(!this->tensor_k_) { //scalar k
                          grad *= ( -1.0* e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                      } else { //tensor k
                          TensorVariable<dim> kk;
                          eptr->Read( key_kk, kk );
                          grad = kk * grad;
                          grad *= ( -1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                      }
                      double grad_nrml = grad.DotProduct(facetNrml);
                      capillary_velocity_component = upstream_lambda_overbar * grad_nrml;
                  }

                  //total facet non-wetting phase flux
                  double f_n = sign * (viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;

                  //total facet wetting phase flux
                  double upstream_fw=(total_mobility!=0.0? lw/total_mobility : 0.0);
                  double viscous_velocity_component_w = vD_n * upstream_fw;
                  double f_w = sign * (viscous_velocity_component_w + gravity_velocity_component - capillary_velocity_component) * facetArea;

                  //update total outflow and inflow
                  total_f = f_n + f_w;
                  if(total_f < 0.) tot_inflow -= total_f;
                  if(total_f > 0.) tot_outflow += total_f;

                  //update phase accumulation
                  carb_accumulation += f_n;
                  carb_inflow += f_n;
                  aq_accumulation += f_w;

                  //determine cfl_multiplier based on non-wetting phase shock saturation
                  if (cfl_multiplier != this->CFL_multiplier_) {
                      if(vn_at_facet_int_point< 0.0) { //flowing in from outside node (upstream node)
                          double sn_shock = 1.0-eptr->Read(this->key_ssH2O); //sn at shock for outside node
                          double sn_outside_node = eptr->N(outside_node)->Read( this->key_sCO2 );
                          if (sn_outside_node >= sn_shock) { //upstream node passed shock saturation
                              double sn_inside_node = eptr->N(inside_node)->Read( this->key_sCO2 );
                              if (sn_inside_node < sn_shock) {//current node not yet reach shock saturation
                                  cfl_multiplier = this->CFL_multiplier_;
                              }
                          }
                      }
                  }
        
                  //determine cfl_multiplier based on wetting phase shock saturation
                  if (cfl_multiplier != this->CFL_multiplier_){
                      if(vw_at_facet_int_point< 0.0) { //flowing from outside node (upstream node)
                          double sw_shock = eptr->Read(this->key_ssH2O); //sw at shock for outside node
                          double sw_outside_node = eptr->N(outside_node)->Read( this->key_sH2O );
                          if (sw_outside_node >= sw_shock) { //upstream node passed shock saturation
                              double sw_inside_node = eptr->N(inside_node)->Read( this->key_sH2O );
                              if (sw_inside_node < sw_shock) {//current node not yet reach shock saturation
                                  cfl_multiplier = this->CFL_multiplier_;
                              }
                          }
                      }
                  }

              } //end sector_facets loop

              /*
              if (!this->no_flow_boundary_) {
                  //inflow/outflow compensation for truncated boundary node
                  if (truncated_node == 1) {
                      //if (inflow > 0.) {flux_balance += inflow; outflow += inflow;} //inflow compensation
                      //else if (inflow < 0.) {flux_balance -= inflow;}; //outflow compensation
                      if (inflow > 0.) {flux_balance += inflow;} //inflow compensation
                      else if (inflow < 0.) {flux_balance -= inflow; outflow -= inflow;}; //outflow compensation
                
                      if ( carb_inflow > 0. ) carb_accumulation += carb_inflow; //inflow compensation
                      else if ( carb_inflow < 0. ) carb_accumulation -= carb_inflow; //outflow compensation
                  }
              }
              */

          } //end parent element loop

          nd->Store( key_compensate, makeScalar( nd->Status( key_compensate ), carb_accumulation/nd->Read(this->key_fvPV) ) );

          /*
          //if(fabs(flux_balance - 0.) > 1.0e-12 )
              //cout<<"  local volume not conserved, event = "<<nd->Read(this->key_EventIndex)<<" flux_balance = "<<flux_balance<<endl;

          // divergence free correction (only when node is located inside domain and flux balance not equal to zero)
          if (nd->Manifold() != nullptr && truncated_node !=1 && fabs(flux_balance) > numeric_limits<double>::epsilon()) {
               // compute average fractional flow for the current finite volume
              double fw_avg(0), fn_avg(0.);
              double n_sn = nd->Read(this->key_sCO2);
              double n_sw = nd->Read(this->key_sH2O); //saturation aqueous phase at current node
              for ( size_t t=0U; t<node_parent_elements; t++ )
              {
                  Element<dim>* const eptr(nd->Parent(t));
                  double ep_swr = eptr->Read(this->key_srH2O);
                  double ep_snr = eptr->Read(this->key_srCO2);

                  const size_t nodes(eptr->Nodes());
                  vector<double> IPOL;
                  eptr->N_AtBaryCenter( IPOL );
                  double ipol_sum(0.), ep_sw(0.), ep_sn(0.), ep_muw (0.), ep_mun(0.), ep_rhow(0.), ep_rhon(0.);
                  for ( size_t i=0U; i<nodes; ++i ) {
                      ep_sw += IPOL[i] * eptr->N(i)->Read( this->key_sH2O );
		                  ep_muw += IPOL[i] * eptr->N(i)->Read( this->key_muH2O );
		                  ep_mun += IPOL[i] * eptr->N(i)->Read( this->key_muCO2 );
		                  ep_rhow += IPOL[i] * eptr->N(i)->Read( this->key_rhoH2O );
		                  ep_rhon += IPOL[i] * eptr->N(i)->Read( this->key_rhoCO2 );
                      ipol_sum += IPOL[i];
                  }
                  if(ipol_sum > 0.) {
		                  ep_sw *= 1. / ipol_sum;
		                  ep_muw *= 1. / ipol_sum;
		                  ep_mun *= 1. / ipol_sum;
		                  ep_rhow *= 1. / ipol_sum;
		                  ep_rhon *= 1. / ipol_sum;
                  }
                  ep_sn = 1.0 - ep_sw;

                  double lambda_w(0.), lambda_n(0.);
                  if(n_sw > ep_swr) lambda_w = this->flowfunctions_.krw_at(eptr, n_sw) / ep_muw;
                  if(n_sn > ep_snr) lambda_n = this->flowfunctions_.krn_at(eptr, n_sw) / ep_mun;
                  double lambda_t = lambda_w + lambda_n;
                  double ep_fw = (lambda_t!=0.0 ? lambda_w/lambda_t : 0.0);
                  double ep_fn = (lambda_t!=0.0 ? lambda_n/lambda_t : 0.0);
                  //fw_avg += ep_fw;
                  fn_avg += ep_fn;
              }
              //fw_avg /= static_cast<double>(node_parent_elements);
              fn_avg /= static_cast<double>(node_parent_elements);
        
              carb_accumulation -= (fn_avg*tot_accumulation);
              //aq_accumulation -= (fw_avg*tot_accumulation);

              //carb_accumulation -= (fn_avg*flux_balance);
              //aq_accumulation -= (fw_avg*flux_balance);

              //double new_tot_accumulation = carb_accumulation + aq_accumulation;
              //if(fabs(new_tot_accumulation - 0.) > numeric_limits<double>::epsilon())
                  //cout<<"  local volume still not conserved after correction, event = "<<nd->Read(this->key_EventIndex)<<" original divergence  = "<<tot_accumulation<<" corrected divergence = "<<new_tot_accumulation<<endl;
          }
          */

          //store rate of changes
          nd->Store( key_dsnw, makeScalar( nd->Status( key_dsnw ), carb_accumulation/nd->Read(  this->key_fvPV ) ) );

          //compute CFL time increment
          ArrayVariable array2;
          nd->Read(this->key_time, array2);

          //compute CFL time increment
          //tot_inflow += nd->Read(key_NQV) * nd->Read(this->key_fvPV);
          //if (tot_inflow < numeric_limits<double>::epsilon())
          if (tot_outflow < numeric_limits<double>::epsilon())
          //if (outflow < numeric_limits<double>::epsilon())
              array2.Component(2, numeric_limits<double>::max());
          else
              array2.Component(2, nd->Read(  this->key_fvPV ) / tot_outflow);
              //array2.Component(2, nd->Read(  this->key_fvPV ) / tot_inflow);
              //array2.Component(2, nd->Read(  this->key_fvPV ) / outflow);

          nd->Store( this->key_CFL, makeScalar( nd->Status( this->key_CFL ),cfl_multiplier ) );
          array2.Component(6, cfl_multiplier);
        
          nd->Store(this->key_time, array2);
    
      }  //end if(PV > numeric_limits<double>::epsilon())
  } //end if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH)
  
}



  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( std::set<Node<dim>*> input_nodes, bool divergence_free_correction)
  {

    for (auto nd : input_nodes) {
      if(nd == nullptr)
        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange(std::set<Node<dim>*> nodes):",
                               "target node is a nullptr" );
      if(nd->Status(this->key_sCO2) == DIRICH)
        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange(std::set<Node<dim>*> nodes):",
                               "this method does not work for DIRICH nodes" );
      if(!nd->IsManifold() )
        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange(std::set<Node<dim>*> nodes):",
                               "this method only works for manifold nodes or periodic nodes" );
    }

    //some variables to use
    double flux_balance(0.), outflow(0.), tot_inflow(0.), tot_outflow(0.), carb_accumulation(0.), aq_accumulation(0.), tot_PV(0.);
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> facetNrml, gravity, gradP;
    vector<double> IPOL, NRST;

    //loop over for nodes in the set
    for(auto nd : input_nodes) {
      double PV = nd->Read(this->key_fvPV); //pore volume (m3)
      if (PV < numeric_limits<double>::epsilon())
        throw csmp::Exception(FATAL_ERROR,
                              "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange:",
                              "FV pore volume seems to be zero, has it been initialised yet?");
      else {
        this->rate_count_++;//recording
        nd->Store(this->key_rate, makeScalar(nd->Status(this->key_rate), nd->Read(this->key_rate) + 1));

        tot_PV += PV;

        size_t node_parent_elements(nd->Parents());
        //double cfl_multiplier = this->CFL_multiplier_ * this->relaxing_factor_; //default value
        long truncated_node = static_cast<long>(nd->Read(this->key_cut));//check if node is truncated by domain boundary

        for (auto t{0U}; t < node_parent_elements; t++) {
          Element<dim> *const eptr(nd->Parent(t));
          assert(eptr != nullptr);

          if (truncated_node == 1 && (this->halo_stencils_.find(eptr) !=
                                      this->halo_stencils_.end())) //ignore if parent element located outside domain
            continue;

          //element state and properties
          const auto nodes(eptr->Nodes());
          eptr->N_AtBaryCenter(IPOL);
          double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw(0.), e_mun(0.), e_rhow(0.), e_rhon(0.);
          for (auto i = 0U; i < nodes; ++i) {
            e_sw += IPOL[i] * eptr->N(i)->Read(this->key_sH2O);
            e_muw += IPOL[i] * eptr->N(i)->Read(this->key_muH2O);
            e_mun += IPOL[i] * eptr->N(i)->Read(this->key_muCO2);
            e_rhow += IPOL[i] * eptr->N(i)->Read(this->key_rhoH2O);
            e_rhon += IPOL[i] * eptr->N(i)->Read(this->key_rhoCO2);
            ipol_sum += IPOL[i];
          }
          if (ipol_sum > 0.) {
            e_sw *= 1. / ipol_sum;
            e_muw *= 1. / ipol_sum;
            e_mun *= 1. / ipol_sum;
            e_rhow *= 1. / ipol_sum;
            e_rhon *= 1. / ipol_sum;
          }
          e_sn = 1.0 - e_sw;

          const double e_swr = eptr->Read(key_srH2O);
          const double e_snr = eptr->Read(key_srCO2);

          eptr->Read(this->key_gradP, gradP); //pressure gradient
          double thickness = eptr->Read(this->key_thi); //thickness
          double e_k = eptr->Read(this->key_k); //permeability

          //read in total velocity computed previously
          VectorVariable<dim> e_vt(ANY, 0.);
          eptr->Read(this->key_vt, e_vt);

          double carb_inflow(0.);
          const auto pnid(nd->ParentNodeNumber(t));
          const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
          for (auto i = 0U; i < sector_facets; i++) {
            const auto iFacet(eptr->FV()->FacetSurroundingSector(pnid, i));
            const auto inside_node(eptr->FV()->InsideNode(iFacet));
            const auto outside_node(eptr->FV()->OutsideNode(iFacet));

            eptr->Read(iFacet, 0U, this->key_fn, facetNrml);
            const double facetArea = eptr->Read(iFacet, 0U, this->key_fA);
            const double sign = (pnid == inside_node) ? 1. : -1.;

            //compute facet fluid flux
            const double vD_n = e_vt.DotProduct(facetNrml); //m/s
            //update flux balance
            double f_total = sign * vD_n * facetArea;
            if (f_total != 0.) {
              flux_balance += f_total;
              //if(f_total < 0.) outflow -= f_total;
              if (f_total > 0.) outflow += f_total;
              //inflow += f_total;
            }

            //compute facet variables
            eptr->N_AtFacetIntegrationPoint(iFacet, 0U, NRST);
            double f_sw(0.), f_muw(0.), f_mun(0.), f_rhow(0.), f_rhon(0.);
            for (auto x{0U}; x < eptr->Nodes(); x++) {
                f_sw += NRST[x] * eptr->N(x)->Read(this->key_sH2O);
                f_muw += NRST[x] * eptr->N(x)->Read(this->key_muH2O);
                f_mun += NRST[x] * eptr->N(x)->Read(this->key_muCO2);
                f_rhow += NRST[x] * eptr->N(x)->Read(this->key_rhoH2O);
                f_rhon += NRST[x] * eptr->N(x)->Read(this->key_rhoCO2);
              }

            //compute upstream directions at facet integration point
            double vn_gravity_component_of_velocity(0.0), vw_gravity_component_of_velocity(0.0);
            double vn_capillary_component_of_velocity(0.0), vw_capillary_component_of_velocity(0.0);

            if (this->with_gravity_forces_) {
              eptr->Read(key_dip, gravity);
              if (isnan(gravity(v))) { //dip vector has not been initialised
                if (dim == 1u) { gravity(0u) = -1.; }
                else if (dim == 2u) {
                  gravity(0u) = 0.;
                  gravity(1u) = -1.;
                }
                else {
                  gravity(0u) = 0.;
                  gravity(1u) = -1.;
                  gravity(2u) = 0.;
                }
                eptr->Store(key_dip, gravity);
              }

              if (!this->tensor_k_) { //scalar k
                gravity *= e_k * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
              } else { //tensor k
                TensorVariable<dim> kk;
                eptr->Read(key_kk, kk);
                VectorVariable<dim> kV;
                kV = kk * gravity;
                double kV_magnitude = kV.Length();
                gravity *= kV_magnitude * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
              }
              double gravity_nrml = gravity.DotProduct(facetNrml);
              vn_gravity_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * gravity_nrml;
              vw_gravity_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * gravity_nrml;
            }

            if (this->with_capillary_spreading_) {
              VectorVariable<dim> grad;
              eptr->Read(this->key_gradSn, grad);
              if (!this->tensor_k_) { //scalar k
                grad *= (-1.0 * e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw));
              } else { //tensor k
                TensorVariable<dim> kk;
                eptr->Read(key_kk, kk);
                grad = kk * grad;
                grad *= (-1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw));
              }
              double grad_nrml = grad.DotProduct(facetNrml);
              vn_capillary_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * grad_nrml;
              vw_capillary_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * grad_nrml;
            }

            double vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            double vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity - vw_capillary_component_of_velocity;

            double total_f(0.), lw(0.), ln(0.);
            //aqueous phase determination
            if (vw_at_facet_int_point != 0.0) {
              const auto upstream_node = (vw_at_facet_int_point > 0.) ? inside_node : outside_node;
              double sw = eptr->N(upstream_node)->Read(this->key_sH2O);
              if (sw > e_swr) lw = this->flowfunctions_.krw_at(eptr, sw) / f_muw;
            } else { //vw_at_facet_int_point == 0.0
              double sw_inside = eptr->N(inside_node)->Read(this->key_sH2O);
              double lw_inside(0.);
              if (sw_inside > e_swr) lw_inside = this->flowfunctions_.krw_at(eptr, sw_inside) / f_muw;

              double sw_outside = eptr->N(outside_node)->Read(this->key_sH2O);
              double lw_outside(0.);
              if (sw_outside > e_swr) lw_outside = this->flowfunctions_.krw_at(eptr, sw_outside) / f_muw;

              lw = (lw_inside + lw_outside) * 0.5;
            }

            //carbonic phase determination
            if (vn_at_facet_int_point != 0.0) {
              const auto upstream_node = (vn_at_facet_int_point > 0.) ? inside_node : outside_node;
              double sw = eptr->N(upstream_node)->Read(this->key_sH2O);
              double sn = eptr->N(upstream_node)->Read(this->key_sCO2);
              if (sn > e_snr) ln = this->flowfunctions_.krn_at(eptr, sw) / f_mun;
            } else { //vn_at_facet_int_point == 0.0
              double sw_inside = eptr->N(inside_node)->Read(this->key_sH2O);
              double sn_inside = eptr->N(inside_node)->Read(this->key_sCO2);
              double ln_inside(0.);
              if (sn_inside > e_snr) ln_inside = this->flowfunctions_.krn_at(eptr, sw_inside) / f_mun;

              double sw_outside = eptr->N(outside_node)->Read(this->key_sH2O);
              double sn_outside = eptr->N(outside_node)->Read(this->key_sCO2);
              double ln_outside(0.);
              if (sn_outside > e_snr) ln_outside = this->flowfunctions_.krn_at(eptr, sw_outside) / f_mun;

              ln = (ln_inside + ln_outside) * 0.5;
            }

            //fractional flows determination
            double total_mobility = lw + ln;
            double upstream_fn = (total_mobility != 0.0 ? ln / total_mobility : 0.0);
            double upstream_lambda_overbar = (total_mobility != 0.0 ? (ln * lw) / total_mobility : 0.0);
            //viscous flux
            double viscous_velocity_component = vD_n * upstream_fn;
            //gravitational flux
            double gravity_velocity_component(0.);
            if (this->with_gravity_forces_) {
              double gravity_nrml = gravity.DotProduct(facetNrml);
              gravity_velocity_component = upstream_lambda_overbar * gravity_nrml;
            }
            //capillary flux
            double capillary_velocity_component(0.);
            if (this->with_capillary_spreading_) {
              VectorVariable<dim> grad;
              eptr->Read(this->key_gradSn, grad);
              if (!this->tensor_k_) { //scalar k
                grad *= (-1.0 * e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw));
              } else { //tensor k
                TensorVariable<dim> kk;
                eptr->Read(key_kk, kk);
                grad = kk * grad;
                grad *= (-1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw));
              }
              double grad_nrml = grad.DotProduct(facetNrml);
              capillary_velocity_component = upstream_lambda_overbar * grad_nrml;
            }

            //total facet non-wetting phase flux
            double f_n =
            sign * (viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;

            //total facet wetting phase flux
            double upstream_fw = (total_mobility != 0.0 ? lw / total_mobility : 0.0);
            double viscous_velocity_component_w = vD_n * upstream_fw;
            double f_w =
            sign * (viscous_velocity_component_w + gravity_velocity_component - capillary_velocity_component) *
            facetArea;

            //update total outflow and inflow
            total_f = f_n + f_w;
            if (total_f < 0.) tot_inflow -= total_f;
            if (total_f > 0.) tot_outflow += total_f;

            //update phase accumulation
            carb_accumulation += f_n;
            carb_inflow += f_n;
            aq_accumulation += f_w;

          } //end sector_facets loop

        } //end parent element loop

      } //end parent element loop

    }

    for(auto nd : input_nodes)
      nd->Store(key_compensate, makeScalar(nd->Status(key_compensate), carb_accumulation / tot_PV));


    if (divergence_free_correction && fabs(flux_balance) > numeric_limits<double>::epsilon()) {
      // compute average fractional flow for the current finite volume
      double fn_avg(0.);
      uint32_t count(0U);
      for(auto nd : input_nodes) {
        double n_sn = nd->Read(this->key_sCO2);
        double n_sw = nd->Read(this->key_sH2O); //saturation aqueous phase at current node
        uint32_t node_parent_elements = nd->Parents();
        for (uint32_t t = 0U; t < node_parent_elements; t++) {
          Element<dim> *const eptr(nd->Parent(t));
          double ep_swr = eptr->Read(this->key_srH2O);
          double ep_snr = eptr->Read(this->key_srCO2);

          const auto nodes(eptr->Nodes());
          eptr->N_AtBaryCenter(IPOL);
          double ipol_sum(0.), ep_sw(0.), ep_sn(0.), ep_muw(0.), ep_mun(0.), ep_rhow(0.), ep_rhon(0.);
          for (auto i{0U}; i < nodes; ++i) {
            ep_sw += IPOL[i] * eptr->N(i)->Read(this->key_sH2O);
            ep_muw += IPOL[i] * eptr->N(i)->Read(this->key_muH2O);
            ep_mun += IPOL[i] * eptr->N(i)->Read(this->key_muCO2);
            ep_rhow += IPOL[i] * eptr->N(i)->Read(this->key_rhoH2O);
            ep_rhon += IPOL[i] * eptr->N(i)->Read(this->key_rhoCO2);
            ipol_sum += IPOL[i];
          }
          if (ipol_sum > 0.) {
            ep_sw *= 1. / ipol_sum;
            ep_muw *= 1. / ipol_sum;
            ep_mun *= 1. / ipol_sum;
            ep_rhow *= 1. / ipol_sum;
            ep_rhon *= 1. / ipol_sum;
          }
          ep_sn = 1.0 - ep_sw;

          double lambda_w(0.), lambda_n(0.);
          if (n_sw > ep_swr) lambda_w = this->flowfunctions_.krw_at(eptr, n_sw) / ep_muw;
          if (n_sn > ep_snr) lambda_n = this->flowfunctions_.krn_at(eptr, n_sw) / ep_mun;
          double lambda_t = lambda_w + lambda_n;
 // not used:         double ep_fw = (lambda_t != 0.0 ? lambda_w / lambda_t : 0.0);
          double ep_fn = (lambda_t != 0.0 ? lambda_n / lambda_t : 0.0);
          //fw_avg += ep_fw;
          fn_avg += ep_fn;

          count++;
        }
      }

      //fw_avg /= static_cast<double>(node_parent_elements);
      fn_avg /= static_cast<double>(count);

      carb_accumulation -= (fn_avg*flux_balance);
      //aq_accumulation -= (fw_avg*flux_balance);
    }


    for(auto nd : input_nodes) {
      //store rate of changes
      nd->Store( key_dsnw, makeScalar( nd->Status( key_dsnw ), carb_accumulation/tot_PV ) );
      //nd->Store( key_dsnw, makeScalar( nd->Status( key_dsnw ), aq_accumulation/tot_PV ) );

      //compute CFL time increment
      ArrayVariable array;
      nd->Read(this->key_time, array);

      //compute CFL time increment
      //tot_inflow += nd->Read(key_NQV) * nd->Read(this->key_fvPV);
      //if (tot_inflow < numeric_limits<double>::epsilon())
      if (tot_outflow < numeric_limits<double>::epsilon()) {
        //if (outflow < numeric_limits<double>::epsilon())
        array.Component(2, numeric_limits<double>::max());
      } else {
        array.Component(2, tot_PV / tot_outflow);
        //array.Component(2, nd1->Read(  this->key_fvPV )  / tot_outflow);
        //array.Component(2, nd->Read(  this->key_fvPV ) / tot_inflow);
      }

      nd->Store( this->key_CFL, makeScalar( nd->Status( this->key_CFL ),this->CFL_multiplier_ ) );
      array.Component(6, this->CFL_multiplier_);

      nd->Store(this->key_time, array);

    }

  }



//schedule an event associated with a node/FV
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
bool DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Schedule(Event<dim>* event, double t_end)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != nullptr );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH) {
      ArrayVariable array;
      nd->Read(this->key_time, array);

      nd->Store( this->key_schedule, makeScalar( nd->Status( this->key_schedule), nd->Read( this->key_schedule) + 1 ) );
      event->valid(true);
      //compute target change
      double dt_CFL = array[2];//CFL time increment
      double ChangeRate = nd->Read( key_dsnw);//rate of change
      double source = nd->Read(key_NQV);
      double dC_CFL = dt_CFL*this->CFL_multiplier_*(-ChangeRate+source);//targe change

      if (fabs(dC_CFL) < numeric_limits<double>::epsilon()){//idle node/FV
          array.Component(5, numeric_limits<double>::epsilon());//target change of solution
          array.Component(3, numeric_limits<double>::max());//target time increment
      } else {
          array.Component(5, dC_CFL);//target change of solution
          array.Component(3, dt_CFL*array[6]);//target time increment
      };

      double t_current = array[0];//current time stamp
      double dt_target = array[3];//target time increment
      if ((dt_target + t_current) >= t_end) {
          nd->Store(this->key_time, array);
          return false;
      } else {
          event->t_schedule(t_current + dt_target);
          array.Component(1, t_current + dt_target);//schedule time stamp
          nd->Store(this->key_time,array);
          return true;
      }
  }
  return false;
}




//update solution and check it against the specified range (with TDS)
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double delta_t)
{
    Node<dim>* nd = event->getNode();
    assert( nd  != nullptr );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH)
    {
      this->update_count_++;//recording
      const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));

      double ChangeRate = nd->Read(key_dsnw);//variation rate
      double solution = nd->Read(this->key_sCO2);//old solution
      nd->Store(key_sCO2_0, makeScalar( status, solution ));//store old solution

      double new_solution = solution - delta_t * ChangeRate;//compute new solution
      const double source(nd->Read( key_NQV));
      new_solution += source * delta_t;//add source to new solution.

      //check new solution value against range and stored it to key_sCO2
      if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
      } else {
        if ( new_solution > this->upper_limit_+1.0e-8 or new_solution < this->lower_limit_-1.0e-8 ) {
          cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_;
          if(nd->IsManifold()) cerr<<", is a manifold node"<<endl;
          else cerr<<endl;
        }
        if ( new_solution > this->upper_limit_) new_solution = this->upper_limit_;
        else if ( new_solution < this->lower_limit_) new_solution = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
        nd->Store(this->key_outrange, makeScalar( nd->Status(this->key_outrange), nd->Read(this->key_outrange) + 1 ) );
      }

      nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
    }
}




//update solution and check it against the specified range (with TDS)
  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_TDS(std::set<Node<dim>*> input_nodes, double delta_t)
  {
    for (auto nd : input_nodes) {
      if(nd == nullptr)
        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_TDS(std::set<Node<dim>*> input_nodes, double delta_t):",
                               "target node is a nullptr" );
      if(nd->Status(this->key_sCO2) == DIRICH)
        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_TDS(std::set<Node<dim>*> input_nodes, double delta_t):",
                               "this method does not work for DIRICH nodes" );
      if(!nd->IsManifold())
        throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange(std::set<Node<dim>*> nodes):",
                               "this method only works for manifold nodes or periodic nodes" );
    }

    double tot_PV(0.), new_tot_vCO2(0.);
    for (auto nd : input_nodes) {
      this->update_count_++;//recording

      double ChangeRate = nd->Read(key_dsnw);//variation rate
      double ori_solution = nd->Read(this->key_sCO2);//old solution
      nd->Store(key_sCO2_0, makeScalar(nd->Status(this->key_sCO2), ori_solution));//store old solution

      double PV = nd->Read(this->key_fvPV);
      tot_PV += PV;

      double ori_vCO2 = PV * ori_solution; //m3
      new_tot_vCO2 += (ori_vCO2 - delta_t * ChangeRate * PV); //m3

      const double source(nd->Read(key_NQV));
      new_tot_vCO2 += PV * source * delta_t;//add source to new solution, m3
    }

    double new_solution = new_tot_vCO2 / tot_PV;

    //check new solution value against range and stored it to key_sCO2
    if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
      for (auto nd: input_nodes) {
        nd->Store(this->key_sCO2, makeScalar(nd->Status(this->key_sCO2), new_solution));
        nd->Store(this->key_sH2O, makeScalar(nd->Status(this->key_sH2O), 1. - new_solution));
        nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
      }
    } else {
      if ( new_solution > this->upper_limit_+1.0e-8 or new_solution < this->lower_limit_-1.0e-8 ) {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_;
        auto nit = input_nodes.begin();
        if((*nit)->IsManifold()) cerr<<", is a manifold node" << endl;
        else cerr<<endl;
      }
      if ( new_solution > this->upper_limit_) new_solution = this->upper_limit_;
      else if ( new_solution < this->lower_limit_) new_solution = this->lower_limit_;
      for (auto nd: input_nodes) {
        nd->Store(this->key_sCO2, makeScalar(nd->Status(this->key_sCO2), new_solution));
        nd->Store(this->key_sH2O, makeScalar(nd->Status(this->key_sH2O), 1. - new_solution));
        nd->Store(this->key_outrange, makeScalar( nd->Status(this->key_outrange), nd->Read(this->key_outrange) + 1 ) );
        nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
      }
    }
  }




//update solution and check it against the specified range (with DES)
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != nullptr );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH)
  {
      this->update_count_++;//recording
      const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
      ArrayVariable array;
      nd->Read(this->key_time, array);

      double ChangeRate = nd->Read( key_dsnw);//variaition rate
      double solution = nd->Read( this->key_sCO2);//old solution
      nd->Store(key_sCO2_0, makeScalar( status, solution ));//store old solution
      double t_current = array[0]; //current time stamp
      double new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
      const double source(nd->Read(key_NQV));
      new_solution += source * (t_clock - t_current);//add source to new solution

      //check new solution value against range and stored it to key_sCO2
      if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
          nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
          nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
      } else {
          if ( new_solution > this->upper_limit_+1.0e-8 or new_solution < this->lower_limit_-1.0e-8 ) {
              cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_;
              if(nd->IsManifold()) cerr<<", is a manifold node"<<endl;
              else cerr<<endl;
          }
          if ( new_solution > this->upper_limit_ ) new_solution = this->upper_limit_;
          else if ( new_solution < this->lower_limit_ ) new_solution = this->lower_limit_;
          nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
          nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
          nd->Store(this->key_outrange, makeScalar( nd->Status(this->key_outrange), nd->Read(this->key_outrange) + 1 ) );
      }
      new_solution = nd->Read(this->key_sCO2);//stored new solution
      double dsn_cumulative = array[4];
      array.Component(4, dsn_cumulative + (new_solution-solution));//update cumulative change
        
      array.Component(0, t_clock); //current time stamp
      array.Component(7, t_current); //previous time stamp
    
      nd->Store(this->key_time, array);

      nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
  }
}




//Synchronize neighbor nodes/FVs
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Synchronize(Event<dim>* event, double t_clock, double& t_remove )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != nullptr );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH)
  {
      nd->Store( this->key_synchronize, makeScalar( nd->Status(this->key_synchronize), nd->Read(this->key_synchronize) + 1 ) );
      event->valid(false);
      ArrayVariable array;
      nd->Read(this->key_time, array);
      array.Component(4, 0.); //reset cumulative change of solution
      nd->Store(this->key_time, array);
      for ( auto n{0U}; n<nd->Neighbors(); ++n ) {
          Node<dim>* neighbor_node = nd->Neighbor(n);
          if( neighbor_node != nullptr && neighbor_node->Status(  this->key_sCO2 ) != DIRICH) {
              size_t index = neighbor_node->Read(this->key_EventIndex);
              if(index >= 0 && index < this->FullList.size()){
                  Event<dim>* neighbor_event = this->FullList[index];
                  assert( neighbor_event  != nullptr );
                  if (neighbor_event != nullptr && neighbor_event->inPEPStack() == false) {
                      this->PEPList.push_back(neighbor_event);
                      neighbor_event->inPEPStack(true);
                      Update_DES(neighbor_event,t_clock);
                      ArrayVariable neighbor_array;
                      neighbor_node->Read(this->key_time, neighbor_array);
                      double dC_cumulative = neighbor_array[4];//cumulative change of solution
                      double dC_target = neighbor_array[5];//target change of solution
                      if (fabs(dC_cumulative) >= fabs(dC_target)) {
                          #if defined(OPENMP)
                          double t_begin = omp_get_wtime();
                          #else
                          clock_t t_begin = clock();
                          #endif
                          if (neighbor_event->inQueue()){
                              Heap_Node* neighbor_heap_node = this->HeapNodeFullList[index];
                              this->EventHeap.remove(neighbor_heap_node);
                              neighbor_event->inQueue(false);
                          }
                          #if defined(OPENMP)
                          t_remove += omp_get_wtime() - t_begin;
                          #else
                          t_remove += clock() - t_begin;
                          #endif
                          Synchronize (neighbor_event, t_clock,t_remove);
                      }
                  }
              }
          }
      }
  }
}




//advect variable with TDS (time driven simulation)
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS( double time_interval, size_t num_threads )
{
    const auto stack_end(this->PEPList.end());
    for ( auto it=this->PEPList.begin(); it!=stack_end; ++it )
    {   
        ComputePressureGradientAndFlowVelocities ((*it));    
        if(this->with_capillary_spreading_) ComputeSaturationGradient ((*it));
    }

#if defined(OPENMP)
    if (num_threads <= 0) {
        cerr <<"WARNING: input number of threads is less than 1, reset to 1"<<endl;
        num_threads = 1;
        AdvectVariable_TDS_parallel( time_interval, num_threads);
    } 
    else if (num_threads > omp_get_max_threads()) {
        cerr <<"WARNING: input number of threads is larger than maximum available threads (" << omp_get_max_threads() << "), reset to "<< omp_get_max_threads() <<endl;
        num_threads = omp_get_max_threads();
        AdvectVariable_TDS_parallel( time_interval, num_threads); 
    }  
    else if (num_threads == 1) {
        AdvectVariable_TDS_serial( time_interval );
    } 
    else {
    AdvectVariable_TDS_parallel( time_interval, num_threads); 
    }            
#else
    if (num_threads > 1)
        cerr <<"WARNING: OPENMP is not available, using serial mode"<<endl;
    AdvectVariable_TDS_serial( time_interval );
#endif
}




//advect variable with DES (discrete event simulation)
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES( double time_increment, double model_time, size_t num_threads )
{

    const auto stack_end(this->FullList.end());
    for ( auto it=this->FullList.begin(); it!=stack_end; ++it )
    {            
        ComputePressureGradientAndFlowVelocities ((*it));    
        if(this->with_capillary_spreading_) ComputeSaturationGradient ((*it));
    }

#if defined(OPENMP)
    if (num_threads <= 0) {
        cerr <<"WARNING: input number of threads is less than 1, reset to 1"<<endl;
        num_threads = 1;
        AdvectVariable_DES_parallel( time_increment, model_time, num_threads);
    } 
    else if (num_threads > omp_get_max_threads()) {
        cerr <<"WARNING: input number of threads is larger than maximum available threads (" << omp_get_max_threads() << "), reset to "<< omp_get_max_threads() <<endl;
        num_threads = omp_get_max_threads();
        AdvectVariable_DES_parallel( time_increment, model_time, num_threads);
    }  
    else if (num_threads == 1) {
        AdvectVariable_DES_serial( time_increment, model_time );
    } 
    else {
        AdvectVariable_DES_parallel( time_increment, model_time, num_threads);
    }            
#else
    if (num_threads > 1) cerr <<"WARNING: OPENMP is not available, using serial mode"<<endl;
    AdvectVariable_DES_serial( time_increment, model_time );
#endif
}





//advect variable with TDS (time-driven simulation), serial mode
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_serial( double time_interval)
{
    cout<<"\nStart DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_serial "<<endl;

    bool verbose(false);
    bool update_pressure(false);

    clock_t begin=clock();
    double time_increment(time_interval); 

    //determine time increment
    clock_t T_begin= clock();
    const auto stack_end(this->PEPList.end());
    for ( auto it=this->PEPList.begin(); it!=stack_end; ++it )
    {
        if(this->with_capillary_spreading_) ComputeSaturationGradient ((*it));
        if(second_order_in_space_) ComputeRateofChange_2nd_order((*it));
        else ComputeRateofChange((*it));
        ArrayVariable array;
        (*it)->getNode()->Read(this->key_time, array);
        double dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
    }

    //if(with_periodic_boundary_) AdjustFlowsAtPeriodicBoundaries(false);
    
    this->T_RateOfChange_ += clock() - T_begin; 

    const double one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double time(0.);
    
    while (time < time_interval)
    {
        T_begin= clock();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        cout <<"\n\tadvection (sub)step: "<< substep <<" of total steps "<<max(floor(time_interval/time_increment),one)<< endl;
        time += time_increment;
        //update saturation
        for ( auto it=this->PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment); //update to time level: time + time_increment
            
            ArrayVariable array;
            (*it)->getNode()->Read(this->key_time, array);
            array.Component(7, time - time_increment); //previous time stamp
            (*it)->getNode()->Store(this->key_time, array);
            
        };

        this->T_Update_ += clock() - T_begin;

        //update node manifolds
        size_t pressure_changed(0);
        if(this->sg_.Mesh().HasNodeManifolds()) {
            for (auto mit = this->sg_.Mesh().NodeManifoldsBegin(); mit != this->sg_.Mesh().NodeManifoldsEnd(); mit++) {
              auto md = (*mit);
              auto master_node = md.N(0);
              for (size_t n = 1; n < md.Branches(); n++) {
                auto slave_node = md.N(n);
                if (UpdateContactStatus(master_node, slave_node)) pressure_changed++;
              }
            }
        }

        size_t updated_manifolds(0);
        if(this->sg_.Mesh().HasNodeManifolds()) {
            for (auto mit = this->sg_.Mesh().NodeManifoldsBegin(); mit != this->sg_.Mesh().NodeManifoldsEnd(); mit++) {
                auto md = (*mit);
                if (UpdateManifold(&md, time, false)) updated_manifolds++;
                //if (UpdateManifold2(&md, time, time_increment, false)) updated_manifolds++; //added
            }
        }
        if(updated_manifolds > 0) cout<<"    updated "<<updated_manifolds<<" manifolds"<<endl;

        //update pressure field
        if(update_pressure && pressure_changed > 0) {
            for ( auto eit = this->sg_.Region("Model").CellsBegin(); eit != this->sg_.Region("Model").CellsEnd(); ++eit )
                ComputeFlowPropertiesAtBaryCenter(*(*eit));
            ComputeSteadyStatePressure();
            for ( auto it=this->PEPList.begin(); it!=this->PEPList.end(); ++it )
                ComputePressureGradientAndFlowVelocities ((*it));
            cout<<"    updated fluid pressure"<<endl;
        }

        //determine new time increment
        T_begin= clock();
        for ( auto it=this->PEPList.begin(); it!=stack_end; ++it )
        {
            if(this->with_capillary_spreading_) ComputeSaturationGradient ((*it));
            if(second_order_in_space_) ComputeRateofChange_2nd_order((*it));
            else ComputeRateofChange((*it));
            ArrayVariable array2;
            (*it)->getNode()->Read(this->key_time, array2);
            double dt_CFL = array2[2];//CFL time increment
            time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
        };

        this->T_RateOfChange_ += clock() - T_begin; 

        substep++;
    };
    
    this->T_AdvectVariable_+= clock() - begin;

    cout<<"\nFinish DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_serial "<<endl;
    if(verbose){
        cout <<"rate_count_ = "<<this->rate_count_<<endl;
        cout <<"update_count_ = "<<this->update_count_<<endl;
        cout <<"T_Schedule_ = "<< this->T_Schedule_ /double(CLOCKS_PER_SEC) << endl;
        cout <<"T_Update_  = "<< this->T_Update_  /double(CLOCKS_PER_SEC) << endl;
        cout <<"T_Synchronize_ = "<< this->T_Synchronize_/double(CLOCKS_PER_SEC) << endl;
        cout <<"T_RateOfChange_ = "<< this->T_RateOfChange_ /double(CLOCKS_PER_SEC) << endl;
        cout <<"T_InsertToHeap_ = "<< this->T_InsertToHeap_ /double(CLOCKS_PER_SEC) << endl;
        cout <<"T_RemoveFromHeap_ = "<< this->T_RemoveFromHeap_ /double(CLOCKS_PER_SEC) << endl;
        cout <<"T_AdvectVariable_ = "<< this->T_AdvectVariable_ /double(CLOCKS_PER_SEC) << endl;
    }
}

  


#if defined(OPENMP)
//advect variable with TDS (time-driven simulation), parallel mode
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_parallel( double time_interval, size_t num_threads )
{
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim>::AdvectVariable_TDS_parallel "<<endl;
    cout<<"Using threads = "<<num_threads<<" Maximum available threads ="<< omp_get_max_threads() << endl;

    bool verbose(false);
    bool update_pressure(false); 

    double begin=omp_get_wtime();
    double time_increment(time_interval); 

    //determine time increment
    double T_begin;
    T_begin = omp_get_wtime();
    
    size_t PEPList_size = this->PEPList.size();

    if(this->with_capillary_spreading_) {
        for(size_t i = 0U; i < PEPList_size; ++i) {
            auto it = this->PEPList.begin()+i;
            ComputeSaturationGradient ((*it));
        }
    }
        
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(dynamic)
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            if(second_order_in_space_) ComputeRateofChange_2nd_order((*it));
            else ComputeRateofChange((*it));
        };
    }

    for(size_t i = 0U; i < PEPList_size; ++i)
    {
        auto it = this->PEPList.begin()+i;
        Event<dim>* event = *it;     
        ArrayVariable array;
        (*it)->getNode()->Read(this->key_time, array);
        double dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
    }
    
    this->T_RateOfChange_ += omp_get_wtime() - T_begin;

    const double one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double time(0.);
    
    while (time < time_interval)
    {
        T_begin = omp_get_wtime();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        cout <<"\n\tadvection (sub)step: "<< substep <<" of total steps "<<max(floor(time_interval/time_increment),one)<< endl;
        time += time_increment;

        //update saturation
        const auto stack_end(this->PEPList.end());
        for ( auto it=this->PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment); //update to time level: time + time_increment
            
            ArrayVariable array;
            (*it)->getNode()->Read(this->key_time, array);
            array.Component(7, time - time_increment); //previous time stamp
            (*it)->getNode()->Store(this->key_time, array);
        };        
        
        this->T_Update_ += omp_get_wtime() - T_begin;

        //update node manifolds
        size_t pressure_changed(0);
        size_t num_manifolds = this->sg_.Mesh().NodeManifolds();
        if(num_manifolds > 0) {
            auto manifold_begin = this->sg_.Mesh().NodeManifoldsBegin();

            vector<size_t> num_pressure_changed(num_threads, 0);
            #pragma omp parallel num_threads(num_threads)
            {
                size_t id = omp_get_thread_num();
                size_t local_pressure_changed(0);
                #pragma omp for schedule(dynamic)
                for(size_t i = 0U; i < num_manifolds; ++i)
                {
                    auto mit = manifold_begin + i;
                    auto md = (*mit);
                    auto master_node = md->N(0);
                    for(size_t n=1;n<md->Nodes();n++) {
                        auto slave_node = md->N(n);
                        if(UpdateContactStatus(master_node, slave_node)) local_pressure_changed++;
                    }
                }
                num_pressure_changed[id] = local_pressure_changed;
            }
            for(size_t n=0; n<num_threads; n++) pressure_changed += num_pressure_changed[n];

            size_t updated_manifolds(0);
            vector<size_t> num_updated_manifolds(num_threads, 0);
            #pragma omp parallel num_threads(num_threads)
            {
                size_t id = omp_get_thread_num();
                size_t local_updated_manifolds(0);
                #pragma omp for schedule(dynamic)
                for(size_t i = 0U; i < num_manifolds; ++i)
                {
                    auto mit = manifold_begin + i;
                    auto md = (*mit);
                    if(UpdateManifold(md, time, false)) local_updated_manifolds++;
                }
                num_updated_manifolds[id] = local_updated_manifolds;
            }
            for(size_t n=0; n<num_threads; n++) updated_manifolds += num_updated_manifolds[n];

            if(updated_manifolds > 0) cout<<"    updated "<<updated_manifolds<<" manifolds"<<endl;
        }


        //update pressure field
        if(update_pressure && pressure_changed > 0) {
            for ( auto eit = this->sg_.Region("Model").CellsBegin(); eit != this->sg_.Region("Model").CellsEnd(); ++eit )
                ComputeFlowPropertiesAtBaryCenter(*(*eit));
            ComputeSteadyStatePressure();
            for ( auto it=this->PEPList.begin(); it!=this->PEPList.end(); ++it )
                ComputePressureGradientAndFlowVelocities ((*it));
            cout<<"    updated fluid pressure"<<endl;
        }

        //compute new time increment
        T_begin = omp_get_wtime();
        time_increment = time_interval; 
        if(this->with_capillary_spreading_) {
            for(size_t i = 0U; i < PEPList_size; ++i) {
                auto it = this->PEPList.begin()+i;
                ComputeSaturationGradient ((*it));
            }
        }
                
        #pragma omp parallel num_threads(num_threads)
        {        
            #pragma omp for schedule(dynamic)
            for(size_t i = 0U; i < PEPList_size; ++i)
            {
                auto it = this->PEPList.begin()+i;
                if(second_order_in_space_) ComputeRateofChange_2nd_order((*it));
                else ComputeRateofChange((*it));
            }
        }        
        
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            ArrayVariable array2;
            (*it)->getNode()->Read(this->key_time, array2);
            double dt_CFL = array2[2];//CFL time increment
            time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
        };

        this->T_RateOfChange_ += omp_get_wtime() - T_begin; 


        substep++;
    };
    
    this->T_AdvectVariable_+= omp_get_wtime() - begin;

    if(verbose) {
        cout << "Finish DES2PhaseSlightlyCompressibleTransport<dim>::AdvectVariable_TDS_parallel " << endl;
        cout << "rate_count_ = " << this->rate_count_ << endl;
        cout << "update_count_ = " << this->update_count_ << endl;
        cout << "T_Schedule_ = " << this->T_Schedule_ << endl;
        cout << "T_Update_  = " << this->T_Update_ << endl;
        cout << "T_Synchronize_ = " << this->T_Synchronize_ << endl;
        cout << "T_RateOfChange_ = " << this->T_RateOfChange_ << endl;
        cout << "T_InsertToHeap_ = " << this->T_InsertToHeap_ << endl;
        cout << "T_RemoveFromHeap_ = " << this->T_RemoveFromHeap_ << endl;
        cout << "T_AdvectVariable_ = " << this->T_AdvectVariable_ << endl;
    }
}
#endif




//advect variable with DES (discrete event simulation), serial mode
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial( double time_increment, double model_time)
{
    clock_t begin=clock();
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial "<<endl;

    /*
    vector<Event<dim>*> source_events;
    for ( auto it=this->FullList.begin(); it!=this->FullList.end(); ++it )       
    {    
        auto event = *it;
        auto nd = event->getNode();
        if(nd->Read(this->key_NQV)>0.) source_events.push_back(event);
    }
    */

    bool verbose(false);
    bool update_pressure(false);

    double time(0.);
    bool Finished = false;

    //uncomment for recording events at each time interval
    /*
    const auto nodes_end(gref_.NodesEnd());
    for ( auto nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
        (*nit)->Store(  this->count_key, makeScalar( (*nit)->Status( this->count_key), 0 ) );
    */

    this->ResetCFLMultiplier();
    
    //used for local pressure solves
    time = model_time - time_increment;
    double time_level(model_time); //set initial time level to advection end time
    double tolerance = 1.; //1 sec tolerance
    size_t pressure_changed(0);

    while (!Finished)
    {
        //performing Iterative local transient pressure computations
        if(update_pressure && pressure_changed > 0) {
            double new_time_level (0.);
            if(!this->EventHeap.empty()) time_level = this->EventHeap.minimum()->getK(); //set time level to first event in the current queue
            double ori_time_level = time_level;
            size_t iter(0);
            //remember starting pressure and saturation
            this->sg_.CopyReplace("fluid pressure","previous fluid pressure");
            do {
                iter++;
                this->sg_.CopyReplace("previous fluid pressure","fluid pressure"); //reset pressure before solving
                double dt = time_level - time;
                new_time_level = time_level;
                if(dt > numeric_limits<double>::epsilon()) ComputeLocalFluidPressure( new_time_level, dt, false );
                time_level = TimeLevel();
                time_level = max(time_level, time); //larger than current time level
                time_level = min(time_level, model_time); //smaller than advection end time level
            } while(fabs(time_level-new_time_level)>tolerance);
            cout<<endl<<"  Iterative local transient pressure computations finished in "<<iter<<" steps, time level adjusted from "<<ori_time_level<<" to "<<time_level<<" dt = "<<time_level - time<<endl;
        }
        
        clock_t T_begin{ clock() };
        const auto stack_end(this->PEPList.end());
        for ( auto it=this->PEPList.begin(); it!=stack_end; ++it )
        {   
            Event<dim>* event = *it;
            T_begin= clock();
            if(this->with_capillary_spreading_) ComputeSaturationGradient (event);
            if(second_order_in_space_) ComputeRateofChange_2nd_order((*it));
            else ComputeRateofChange((*it));
            this->T_RateOfChange_ += clock() - T_begin;         
            if ((*it)->valid() == false) {
                T_begin= clock();
                bool isactive = Schedule(event, model_time);   
                this->T_Schedule_ += clock() - T_begin;              
                if (isactive) {
                    T_begin= clock();
                    double scheduled_time = event->t_schedule();
                    int index = static_cast<int>(event->getNode()->Read(this->key_EventIndex));
                    auto* heap_node = new Heap_Node(scheduled_time,index);
                    this->EventHeap.insert(heap_node);
                    this->HeapNodeFullList[index] = heap_node;
                    event->inQueue(true);
                    this->T_InsertToHeap_ += clock() - T_begin; 
                }  
            } 
            event->inPEPStack(false);         
        };    
        this->T_RateOfChange_ += clock() - T_begin; 

        //if(with_periodic_boundary_) AdjustFlowsAtPeriodicBoundaries(true);
            
        if (this->EventHeap.empty()) time=model_time;
        else time = this->EventHeap.minimum()->getK();
        
        cout<<"  time = "<<time<<" model_time = "<<model_time<<" PEPList size = "<< this->PEPList.size() <<" Queue size = "<< this->EventHeap.size()<<endl;

        if (time == model_time) {
            Finished = true;
            const auto End(this->PEPList.end());
            for ( auto e=this->PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        this->PEPList.clear();
    
        double dt_PEP=numeric_limits<double>::max();
        size_t count = 0U;
        while (!this->EventHeap.empty())
        {           
            Heap_Node* root_node = this->EventHeap.minimum();
            size_t top_index = root_node->getV();
            Event<dim>* top_event = this->FullList[top_index];

            if(top_event->valid() == false) {
                T_begin= clock();
                this->EventHeap.remove(root_node);
                top_event->inQueue(false);
                this->T_RemoveFromHeap_ += clock() - T_begin;
                continue;
            }
            
            count++;            
            ArrayVariable array;
            top_event->getNode()->Read(this->key_time, array);
            double dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, this->PEP_multiplier_*dt_target);
            double t_schedule = array[1];//scheduled time stamp
            if (t_schedule > (time+dt_PEP)) break;            
            if (top_event->inPEPStack() == false) {
                this->PEPList.push_back(top_event);
                top_event->inPEPStack(true);
                T_begin= clock();
                Update_DES(top_event,time);
                this->T_Update_ += clock() - T_begin;                
            };  
            
            T_begin= clock();
            this->EventHeap.remove(root_node); 
            top_event->inQueue(false);
            this->T_RemoveFromHeap_ += clock() - T_begin;
            
            T_begin= clock();
            double t_remove(0.);    
            Synchronize(top_event,time,t_remove);

            //if(this->with_periodic_boundary_) SynchronizePeriodicNode(top_event, time);
            this->T_RemoveFromHeap_ += t_remove;
            this->T_Synchronize_ += clock() - T_begin - t_remove;
        };


        /*
        //update nodal source nodes if not updated yet
        size_t source_update_count(0U);
        for(auto event : source_events)   
        {
            auto nd = event->getNode();

            if (event->inPEPStack() == false) {
                this->PEPList.push_back(event);
                event->inPEPStack(true); 
                Update_DES(event,time);
                if (event->inQueue()) {
                    int index = static_cast<int>(event->getNode()->Read(this->key_EventIndex));
                    if(index >= 0 && index < this->FullList.size()) {
                        Heap_Node* heap_node = this->HeapNodeFullList[index];
                        this->EventHeap.remove(heap_node);
                        event->inQueue(false);
                    }
                }
                double t_remove(0.);
                Synchronize(event,time,t_remove);

                source_update_count++;
             }
         }  
         cout<<"    updated and synchronised "<<source_update_count<<" source nodes"<<endl;
         */
        
        size_t n_manifolds = this->sg_.Mesh().NodeManifolds();
        pressure_changed = 0;
        if(n_manifolds > 0) {
            std::set<NodeManifold<dim>*> manifolds;
            for(auto event : this->PEPList) {
                auto md = event->getNode()->Manifold();
                if(md != nullptr) manifolds.insert(md);
            }
              
            size_t updated_manifolds(0);
            if(manifolds.size() > 0) {
                for(auto md : manifolds) {
                    auto master_node = md->N(0);
                    for(size_t n=1;n<md->Branches();n++){
                        auto slave_node = md->N(n);
                        if(UpdateContactStatus(master_node, slave_node)) pressure_changed++;
                    }                   
                    if(UpdateManifold(md, time, true)) {updated_manifolds++;}
                }

                for(auto md : manifolds) SynchronizeManifold(md, time);
                if(updated_manifolds > 0) cout<<"    updated "<<updated_manifolds<<" manifolds"<<endl; 
            }
        }   

        //if(with_periodic_boundary_) AveragePeriodicBoundaryNodeSaturations();

    };
    /*
    //reset all events and add them to PEPList (for advection at next integration step)
    this->EventHeap.clear();
    this->PEPList.clear();
    const auto stack_end(this->FullList.end());
    for ( auto it=this->FullList.begin(); it!=stack_end; ++it )
    {  
        Event<dim>* event = *it;
        PEPList.push_back(event);
        event->inPEPStack(true);          
        event->valid(false);
        event->inQueue(false);
    }    
    */
    

    //update all contact status before next global pressure solving
    if(this->sg_.Mesh().NodeManifolds() > 0) {
        for(auto mit = this->sg_.Mesh().NodeManifoldsBegin();mit!=this->sg_.Mesh().NodeManifoldsEnd();mit++) {
            auto md = (*mit);     
            auto master_node = md.N(0);
            for(size_t n=1;n<md.Branches();n++){
                auto slave_node = md.N(n);
                UpdateContactStatus(master_node, slave_node);
            }
        }
    }
    
    this->T_AdvectVariable_+= clock() - begin;

    if(verbose) {
        cout << "Finish DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial " << endl;
        cout << "rate_count_ = " << this->rate_count_ << endl;
        cout << "update_count_ = " << this->update_count_ << endl;
        cout << "T_Schedule_ = " << this->T_Schedule_ / double(CLOCKS_PER_SEC) << endl;
        cout << "T_Update_  = " << this->T_Update_ / double(CLOCKS_PER_SEC) << endl;
        cout << "T_Synchronize_ = " << this->T_Synchronize_ / double(CLOCKS_PER_SEC) << endl;
        cout << "T_RateOfChange_ = " << this->T_RateOfChange_ / double(CLOCKS_PER_SEC) << endl;
        cout << "T_InsertToHeap_ = " << this->T_InsertToHeap_ / double(CLOCKS_PER_SEC) << endl;
        cout << "T_RemoveFromHeap_ = " << this->T_RemoveFromHeap_ / double(CLOCKS_PER_SEC) << endl;
        cout << "T_AdvectVariable_ = " << this->T_AdvectVariable_ / double(CLOCKS_PER_SEC) << endl;
    }
}   



#if defined(OPENMP)
//advect variable with DES (discrete event simulation), parallel mode
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_parallel( double time_increment,double model_time, size_t num_threads)
{
    double begin=omp_get_wtime();
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_parallel "<<endl;
    cout << "Using threads = "<<num_threads<<" Maximum available threads ="<< omp_get_max_threads() << endl;

    bool verbose(false);
    bool update_pressure(false); 

    double time(0.);
    bool Finished = false;

    //uncomment for recording events at each time interval
    /*
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    {     
        (*nit)->Store(  this->count_key, makeScalar( (*nit)->Status( this->count_key), 0 ) );
    }
    */
    
    this->ResetCFLMultiplier();
    
    //variables used for iterative local pressure solving
    time = model_time - time_increment;
    double previous_time = time;
    double time_level(model_time); //set initial time level to advection end time
    double tolerance = 1.; //1 sec tolerance
    size_t pressure_changed(0);

    while (!Finished)
    {
        //aterative local transient pressure computations
        if(update_pressure && pressure_changed > 0) {
            double new_time_level (0.);
            if(!this->EventHeap.empty()) time_level = this->EventHeap.minimum()->getK(); //set time level to first event in the current queue
            double ori_time_level = time_level;
            size_t iter(0);
            //remember starting pressure and saturation
            this->sg_.CopyReplace("fluid pressure","previous fluid pressure");
            do {
                iter++;
                this->sg_.CopyReplace("previous fluid pressure","fluid pressure"); //reset pressure before solving
                double dt = time_level - time;
                new_time_level = time_level;
                if(dt > numeric_limits<double>::epsilon()) ComputeLocalFluidPressure( new_time_level, dt, false );
                time_level = TimeLevel();
                time_level = max(time_level, time); //larger than current time level
                time_level = min(time_level, model_time); //smaller than advection end time level
            } while(fabs(time_level-new_time_level)>tolerance);
            cout<<endl<<"  Iterative local transient pressure computations finished in "<<iter<<" steps, time level adjusted from "<<ori_time_level<<" to "<<time_level<<" dt = "<<time_level - time<<endl;
        }

        double T_begin;
        
        T_begin = omp_get_wtime();
                
        size_t PEPList_size = this->PEPList.size();
        
        if(this->with_capillary_spreading_) {
            for(size_t i = 0U; i < PEPList_size; ++i) {
                auto it = this->PEPList.begin()+i;
                Event<dim>* event = *it;
                ComputeSaturationGradient (event);
            }
        }

        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp for schedule(dynamic)
            //#pragma omp for schedule(static)
            for(size_t i = 0U; i < PEPList_size; ++i)
            {
                auto it = this->PEPList.begin()+i;
                Event<dim>* event = *it;
                if(second_order_in_space_) ComputeRateofChange_2nd_order((*it));
                else ComputeRateofChange((*it));
            };
        } 
        
        this->T_RateOfChange_ += omp_get_wtime() - T_begin; 
        
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            Event<dim>* event = *it;        

            if ((*it)->valid() == false) {
                T_begin= omp_get_wtime();
                bool isactive = Schedule(event, model_time);   
                this->T_Schedule_ += omp_get_wtime() - T_begin;              
                if (isactive) {
                    T_begin= omp_get_wtime();
                    double scheduled_time = event->t_schedule();
                    int index = static_cast<int>(event->getNode()->Read(this->key_EventIndex));
                    auto* heap_node = new Heap_Node(scheduled_time,index);
                    this->EventHeap.insert(heap_node);
                    this->HeapNodeFullList[index] = heap_node;
                    event->inQueue(true);
                    this->T_InsertToHeap_ += omp_get_wtime() - T_begin; 
                }  
            } 
            event->inPEPStack(false);  
        }          
        
        //if(with_periodic_boundary_) AdjustFlowsAtPeriodicBoundaries(true);
            
        if (this->EventHeap.empty()) time=model_time;
        else time = this->EventHeap.minimum()->getK();
        
        cout<<"  time = "<<time<<" model_time = "<<model_time<<" PEPList size = "<< this->PEPList.size() <<" Queue size = "<< this->EventHeap.size()<<endl;

        if (time == model_time) {
            Finished = true;
            const auto End(this->PEPList.end());
            for ( auto e=this->PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        this->PEPList.clear();
    
        double dt_PEP=numeric_limits<double>::max();
        size_t count = 0U;
        while (!this->EventHeap.empty())
        {           
            Heap_Node* root_node = this->EventHeap.minimum();
            size_t top_index = root_node->getV();
            Event<dim>* top_event = this->FullList[top_index];
            
            if(top_event->valid() == false) {
                T_begin= omp_get_wtime();
                this->EventHeap.remove(root_node);
                top_event->inQueue(false);
                this->T_RemoveFromHeap_ += omp_get_wtime() - T_begin;
                continue;
            }
            
            count++;            
            ArrayVariable array;
            top_event->getNode()->Read(this->key_time, array);
            double dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, this->PEP_multiplier_*dt_target);
            double t_schedule = array[1];//scheduled time stamp
            if (t_schedule > (time+dt_PEP)) break;            
            if (top_event->inPEPStack() == false) {
                this->PEPList.push_back(top_event);
                top_event->inPEPStack(true);
                T_begin= omp_get_wtime();
                Update_DES(top_event,time);
                this->T_Update_ += omp_get_wtime() - T_begin;                
            };  
            
            T_begin= omp_get_wtime();
            this->EventHeap.remove(root_node); 
            top_event->inQueue(false);
            this->T_RemoveFromHeap_ += omp_get_wtime() - T_begin;
            
            T_begin= omp_get_wtime();
            double t_remove(0.);    
            Synchronize(top_event,time,t_remove);
            //if(this->with_periodic_boundary_) SynchronizePeriodicNode(top_event, time);
            this->T_RemoveFromHeap_ += t_remove;
            this->T_Synchronize_ += omp_get_wtime() - T_begin - t_remove;

        };

        size_t n_manifolds = this->sg_.Mesh().NodeManifolds();
        if(n_manifolds > 0) {
            std::set<NodeManifold<dim>*> manifolds;
            for(auto event : this->PEPList) {
                auto md = event->getNode()->ParentManifold();
                if(md != nullptr) manifolds.insert(md);
            }

            vector<NodeManifold<dim>*> manifolds_vector;
            manifolds_vector.assign(manifolds.begin(), manifolds.end()); 
            size_t num_manifolds = manifolds_vector.size();

            size_t pressure_changed(0);
            auto manifold_begin = manifolds_vector.begin();

            vector<size_t> num_pressure_changed(num_threads, 0);
            #pragma omp parallel num_threads(num_threads)
            {        
                size_t id = omp_get_thread_num();
                size_t local_pressure_changed(0);
                #pragma omp for schedule(dynamic)    
                //#pragma omp for schedule(static)  
                for(size_t i = 0U; i < num_manifolds; ++i)
                {
                    auto mit = manifold_begin + i;
                    auto md = (*mit); 
                    auto master_node = md->N(0);
                    for(size_t n=1;n<md->Branches();n++){
                        auto slave_node = md->N(n);
                        if(UpdateContactStatus(master_node, slave_node)) local_pressure_changed++;
                    }
                }
                num_pressure_changed[id] = local_pressure_changed;
            }  
            for(size_t n=0; n<num_threads; n++) pressure_changed += num_pressure_changed[n]; 

            size_t total_updated_manifolds(0);
            vector<size_t> num_updated_manifolds(num_threads, 0);
            #pragma omp parallel num_threads(num_threads)
            {        
                size_t id = omp_get_thread_num();
                size_t local_updated_manifolds(0);
                #pragma omp for schedule(dynamic)    
                //#pragma omp for schedule(static)  
                for(size_t i = 0U; i < num_manifolds; ++i)
                {
                    auto mdit = manifold_begin  + i;
                    auto md = *mdit; 
                    if(UpdateManifold(md, time, true)) local_updated_manifolds++; 
                }
                num_updated_manifolds[id] = local_updated_manifolds;
            }               
            for(size_t n=0; n<num_threads; n++) total_updated_manifolds += num_updated_manifolds[n]; 
            for(auto md : manifolds) SynchronizeManifold(md, time);
            if(total_updated_manifolds > 0) cout<<"    updated "<<total_updated_manifolds<<" manifolds"<<endl;
        }      

        //if(with_periodic_boundary_) AveragePeriodicBoundaryNodeSaturations();

    };
    /*
    //reset all events and add them to PEPList (for advection at next integration step)
    this->EventHeap.clear();
    this->PEPList.clear();
    const typename vector<Event<dim>*>::iterator stack_end(this->FullList.end());
    for ( typename vector<Event<dim>*>::iterator it=this->FullList.begin(); it!=stack_end; ++it )       
    {  
        Event<dim>* event = *it;
        this->PEPList.push_back(event);
        event->inPEPStack(true);          
        event->valid(false);
        event->inQueue(false);
    } 
    */

    //update all contact status before next global pressure solving
    if(this->sg_.Mesh().NodeManifolds() > 0){
        for(auto mit = this->sg_.Mesh().NodeManifoldsBegin();mit!=this->sg_.Mesh().NodeManifoldsEnd();mit++) {
            auto md = (*mit);     
            auto master_node = md.N(0);
            for(size_t n=1;n<md.Branches();n++){
                auto slave_node = md.N(n);
                UpdateContactStatus(master_node, slave_node);
            }
        }
    }   
        
    this->T_AdvectVariable_+= omp_get_wtime() - begin;

    if(verbose) {
        cout << "Finish DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_parallel " << endl;
        cout << "rate_count_ = " << this->rate_count_ << endl;
        cout << "update_count_ = " << this->update_count_ << endl;
        cout << "T_Schedule_ = " << this->T_Schedule_ << endl;
        cout << "T_Update_  = " << this->T_Update_ << endl;
        cout << "T_Synchronize_ = " << this->T_Synchronize_ << endl;
        cout << "T_RateOfChange_(including_T_Schedule_) = " << this->T_RateOfChange_ << endl;
        cout << "T_InsertToHeap_ = " << this->T_InsertToHeap_ << endl;
        cout << "T_RemoveFromHeap_ = " << this->T_RemoveFromHeap_ << endl;
        cout << "T_AdvectVariable_ = " << this->T_AdvectVariable_ << endl;
    }
}

#endif




//determine node manifold status
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
bool DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::UpdateContactStatus( Node<dim>* masterNode,Node<dim>* slaveNode )
{

    double master_sn = masterNode->Read(this->key_sCO2);
    double slave_sn = slaveNode->Read(this->key_sCO2);

//    Element<dim>* master_e = ParentElementOfManifoldNode(masterNode, this->key_pd);
    Element<dim>* slave_e = ParentElementOfManifoldNode(slaveNode, this->key_pd);

    bool pressureChanged(false);
    double numEpsilon(1.0e-14);

    // 1. wetting phase at slave nodes has increased, auto breakthrough
    double slave_initial_sn = slaveNode->Read(key_sCO2_initial);
    if(slave_sn > slave_initial_sn + numEpsilon) {
        slaveNode->Store( this->key_breakthrough, makeScalar( slaveNode->Status( this->key_breakthrough), 1 ) );
        VARIABLE_FLAG old_p_statues = slaveNode->Status(this->key_pf);
        if(slaveNode->Status(this->key_pf) != DIRICH) slaveNode->Status(this->key_pf, ROBIN);
        VARIABLE_FLAG new_p_statues = slaveNode->Status(this->key_pf);
        if(old_p_statues != new_p_statues) pressureChanged = true;
        //cerr<<"breakthrough occurs at condition 1"<<endl; //added
        return pressureChanged;
    }

    // 2. only water at both
    if (slave_sn < numEpsilon && master_sn < numEpsilon) {
        slaveNode->Store( this->key_breakthrough, makeScalar( slaveNode->Status( this->key_breakthrough), 0 ) );
        VARIABLE_FLAG old_p_statues = slaveNode->Status(this->key_pf);
        if(slaveNode->Status(this->key_pf) != DIRICH) slaveNode->Status(this->key_pf, ROBIN);
        VARIABLE_FLAG new_p_statues = slaveNode->Status(this->key_pf);
        if(old_p_statues != new_p_statues) pressureChanged = true;
        return pressureChanged;
    }

    // 3. oil dams up at slave side then it breakthroughs
    double master_sw = masterNode->Read(this->key_sH2O);
    Element<dim>* master_e = ParentElementOfManifoldNode(masterNode, this->key_pd);
    double master_pc = this->flowfunctions_.pc_at(master_e, master_sw);
    double slave_pe = slave_e->Read(this->key_pd);
    //double master_pf = masterNode->Read(this->key_pf);
    //double slave_pf = slaveNode->Read(this->key_pf);
    if (master_pc > slave_pe + numEpsilon) {
    //if (master_pf - slave_pf > slave_pe + numEpsilon) {
        slaveNode->Store( this->key_breakthrough, makeScalar( slaveNode->Status( this->key_breakthrough), 1 ) );
        VARIABLE_FLAG old_p_statues = slaveNode->Status(this->key_pf);
        if(slaveNode->Status(this->key_pf) != DIRICH) slaveNode->Status(this->key_pf, ROBIN);
        VARIABLE_FLAG new_p_statues = slaveNode->Status(this->key_pf);
        if(old_p_statues != new_p_statues) pressureChanged = true;
        //cerr<<"breakthrough occurs at condition 3, diff_P = "<<master_pf - slave_pf<<", pd = "<<slave_pe<<endl; //added
        return pressureChanged;
    }
/*
    // 4. flow from lower permeability to higher permeability which is partially saturated by non-wetting phase
    if (ComputeFlowPotential(masterNode) < ComputeFlowPotential(slaveNode)) {
        slaveNode->Store( this->key_breakthrough, makeScalar( slaveNode->Status( this->key_breakthrough), 0 ) );
        VARIABLE_FLAG old_p_statues = slaveNode->Status(this->key_pf);
        if(slaveNode->Status(this->key_pf) != DIRICH) slaveNode->Status(this->key_pf, ROBIN);
        VARIABLE_FLAG new_p_statues = slaveNode->Status(this->key_pf);
        if(old_p_statues != new_p_statues) pressureChanged = true;
        return pressureChanged;
    }
*/
    // 5. last case, no exchange at all
    slaveNode->Store( this->key_breakthrough, makeScalar( slaveNode->Status( this->key_breakthrough), 0 ) );
    VARIABLE_FLAG old_p_statues = slaveNode->Status(this->key_pf);
    if(slaveNode->Status(this->key_pf) != DIRICH) slaveNode->Status(this->key_pf, PLAIN);
    VARIABLE_FLAG new_p_statues = slaveNode->Status(this->key_pf);
    if(old_p_statues != new_p_statues) pressureChanged = true;
    return pressureChanged;

};




//find parent element of a manifold node based on largest property value
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
Element<dim>* DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ParentElementOfManifoldNode(Node<dim>* node, const Index& index)
{
    Element<dim>* parent_elmt = node->Parent(0);
    double largest_value = node->Parent(0)->Read(index);
    for( auto e{1u}; e < node->Parents(); e++) {
        double value = node->Parent(e)->Read(index);
        if(value > largest_value) {parent_elmt = node->Parent(e); largest_value = value; };
    }

    return parent_elmt;
}



//compute average potential that close to a node
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeFlowPotential(Node<dim>* node)
{
    double accumulatedPressure = 0;
    double accumulatedGravityTerm = 0;
    size_t neighborNodes = 0;

    for ( auto i{0U}; i < node->Parents(); ++i) {
        auto const& elem = node->Parent(i);
        for ( auto j{0U}; j < elem->Nodes(); ++j) {
            // accumulate pressure, node can be accumulated multiple times
            auto const& nd = elem->N(j);
            accumulatedPressure += nd->Read(this->key_pf);
            neighborNodes ++;
        }
        // accumulate gravity term g*h_center*relative_density
        double relativeDensity = this->flowfunctions_.f(elem, 0) * elem->PropertyValueAtBaryCenter(this->key_rhoH2O)
                               + this->flowfunctions_.f(elem, 1) * elem->PropertyValueAtBaryCenter(this->key_rhoCO2);
        accumulatedGravityTerm += this->sg_.Read(key_g) * relativeDensity * elem->BaryCenter()[1];
    }

    return accumulatedPressure / double(neighborNodes) + accumulatedGravityTerm / node->Parents();
}



//update node manifold
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
bool DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::UpdateManifold(NodeManifold<dim>* md, double t_clock, bool use_DES)
{
    assert(md != nullptr);

    //breakthrough status should have been determined previously by UpdateContactStatus function
    size_t breakthrough_count(0);
    for(size_t n=1;n<md->Branches();n++){
        auto slave_node = md->N(n);
        if(slave_node->Read(this->key_breakthrough) == 1) breakthrough_count++;
    }
    
    if(breakthrough_count <= 0) { //if there is no breakthrough
        return false;
      
    } else { //if there is breakthrough occurs
        //master node
        auto master_node = md->N(0);
	      Element<dim>* master_e = ParentElementOfManifoldNode(master_node, this->key_pd);
        //double master_swr = master_e->Read(this->key_srH2O);
        //double upperSaturation = 1 - master_swr;
        //double lowerSaturation = master_e->Read(this->key_srCO2);
        double upperSaturation(1.0), lowerSaturation(0.0);
        Compensator compensator(this, md, t_clock);
        double absoluteConvergenCriteria = 1e-14;
        int maxBrentSearch = 50;
        double sn = brent_solve<Compensator>(compensator, lowerSaturation, upperSaturation,  absoluteConvergenCriteria, maxBrentSearch);
        //if(sn>1.-master_e->Read(this->key_srH2O) or sn<master_e->Read(this->key_srCO2)) {
        if(sn > 1. or sn < 0.) {
            double sn_ori = sn;
            //sn = min(sn, 1.-master_e->Read(this->key_srH2O));
            //sn = max(sn, master_e->Read(this->key_srCO2));
            sn = min(sn, 1.0);
            sn = max(sn, 0.0);
            cerr<<"    master_sn = "<<sn_ori<<", adjusted to "<<sn<<endl;
        }
        double old_sn = master_node->Read(this->key_sCO2);
        master_node->Store(this->key_sCO2, makeScalar(master_node->Status(this->key_sCO2), sn));
        master_node->Store(this->key_sH2O, makeScalar(master_node->Status(this->key_sH2O), 1. - sn));

        if(use_DES) {
            ArrayVariable array;
            master_node->Read(this->key_time, array);
            array.Component(0, t_clock); //current time stamp
            array.Component(7, t_clock); //previous time stamp
            double dsn_cumulative = array[4];
            array.Component(4, dsn_cumulative + (sn-old_sn));//update cumulative change
            master_node->Store(this->key_time, array);
        }

        //slave nodes
        for(size_t n=1;n<md->Branches();n++){
            auto slave_node = md->N(n); 
            int isbreakthrough = static_cast<int>(slave_node->Read(this->key_breakthrough));
            if(isbreakthrough == 1) {
                Element<dim>* slave_e = ParentElementOfManifoldNode(slave_node, this->key_pd);
                double master_sn = master_node->Read(this->key_sCO2);
                double slave_sn = slave_node->Read(this->key_sCO2);
	              double master_pc = this->flowfunctions_.pc_at(master_e, 1.0 - master_sn);
                double slave_pd = slave_e->Read(this->key_pd);
                double slave_sw_at_master_pc = this->flowfunctions_.seff_to_sw(slave_e, 1.);
                if(master_pc >= slave_pd) {
                    slave_sw_at_master_pc = this->flowfunctions_.sw_from_pc_at(slave_e, master_pc, 1. - slave_sn);
	                  //if(slave_sw_at_master_pc > 1 - slave_e->Read(this->key_srCO2) or slave_sw_at_master_pc < slave_e->Read(this->key_srH2O)) {
                    if(slave_sw_at_master_pc > 1.0 or slave_sw_at_master_pc < 0.0) {
                        double ori_sw = slave_sw_at_master_pc;
	                      //slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1 - slave_e->Read(this->key_srCO2));
                        //slave_sw_at_master_pc = max(slave_sw_at_master_pc, slave_e->Read(this->key_srH2O));
                        slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1.0);
                        slave_sw_at_master_pc = max(slave_sw_at_master_pc, 0.0);
                        cerr<<"    slave_sw = "<<ori_sw<<", adjusted to "<<slave_sw_at_master_pc<<endl;
	                  }
                }
                slave_node->Store(this->key_sCO2, makeScalar(slave_node->Status(this->key_sCO2), 1. - slave_sw_at_master_pc));
                slave_node->Store(this->key_sH2O, makeScalar(slave_node->Status(this->key_sH2O), slave_sw_at_master_pc));   
               
                if(use_DES) {
                    ArrayVariable array;
                    slave_node->Read(this->key_time, array);
                    array.Component(0, t_clock); //current time stamp
                    array.Component(7, t_clock); //previous time stamp
                    double dsn_cumulative = array[4];
                    double slave_new_sn = slave_node->Read(this->key_sCO2);
                    array.Component(4, dsn_cumulative + (slave_new_sn-slave_sn));//update cumulative change
                    slave_node->Store(this->key_time, array);
                }                                
            }
        }
        
        return true;
    }
}


//update node manifold
  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  bool DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::UpdateManifold2(NodeManifold<dim>* md, double t_clock, double dt, bool use_DES) {
    assert(md != nullptr);
    std::set<Node<dim> *> connected_nodes;
    std::set<Node<dim> *> nodes_to_update;

    //breakthrough status should have been determined previously by UpdateContactStatus function
    size_t breakthrough_count(0);
    for (size_t n = 1; n < md->Branches(); n++) {
      auto slave_node = md->N(n);
      if (slave_node->Read(this->key_breakthrough) == 1) {
        breakthrough_count++;
        connected_nodes.insert(slave_node);
      } else
        nodes_to_update.insert(slave_node);
    }

    if (connected_nodes.size() > 0)
      connected_nodes.insert(md->N(0));
    else
      nodes_to_update.insert(md->N(0));

    if (connected_nodes.size() > 1) {
      ComputeRateofChange(connected_nodes, false);
      Update_TDS(connected_nodes, dt);
    }

    if (nodes_to_update.size() > 0) {
      for (auto nd: nodes_to_update) {
        size_t index = nd->Read(this->key_EventIndex);
        auto event = this->FullList[index];
        Update_TDS(event, dt);
      }
    }

    if(breakthrough_count <= 0)  //if there is no breakthrough
      return false;

    if(use_DES) {
      ArrayVariable array;
      for (uint32_t n = 0; n < md->Branches(); n++) {
        auto nd = md->N(n);
        nd->Read(this->key_time, array);
        array.Component(0, t_clock); //current time stamp
        array.Component(7, t_clock); //previous time stamp
        double dsn_cumulative = array[4];
        double old_sn = nd->Read(this->key_sCO2_0);
        double sn = nd->Read(this->key_sCO2);
        array.Component(4, dsn_cumulative + (sn - old_sn));//update cumulative change
        nd->Store(this->key_time, array);
      }
    }

    return true;
  }




//synchronize manifold
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::SynchronizeManifold(NodeManifold<dim>* md, double t_clock)
{
    assert(md != nullptr);
    //breakthrough status should have been determined previously by UpdateContactStatus function
    size_t breakthrough_count(0);
    for(size_t n=1;n<md->Branches();n++){
        auto slave_node = md->N(n);
        if(slave_node->Read(this->key_breakthrough) == 1) breakthrough_count++;
    }
    
    if(breakthrough_count > 0) {
        for(size_t n=0;n<md->Branches();n++){
            auto nd = md->N(n);
            if(n!=0 && nd->Read(this->key_breakthrough) == 0) continue;    
            size_t nd_index = nd->Read(this->key_EventIndex);
            assert((nd_index >= 0 && nd_index < this->FullList.size()));
            if(nd_index >= 0 && nd_index < this->FullList.size()){
                Event<dim>* event = this->FullList[nd_index];
                assert( event  != nullptr );
                if (event != nullptr && event->inPEPStack() == false) {
                    this->PEPList.push_back(event);
                    event->inPEPStack(true);
                    ArrayVariable array;
                    nd->Read(this->key_time, array);
                    double dC_cumulative = array[4];//cumulative change of solution
                    double dC_target = array[5];//target change of solution
                    if (fabs(dC_cumulative) >= fabs(dC_target)) {
                        if (event->inQueue()){
                            Heap_Node* heap_node = this->HeapNodeFullList[nd_index];
                            this->EventHeap.remove(heap_node);
                            event->inQueue(false);
                        }                
                        clock_t T_begin= clock();
                        double t_remove(0.);
                        Synchronize(event,t_clock,t_remove);
                        this->T_RemoveFromHeap_ += t_remove;
                        this->T_Synchronize_ += clock() - T_begin - t_remove;
                    }
                }
            }
        }
    } 
}
 



template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Compensator::Compensator(DES2PhaseSlightlyCompressibleTransport<dim, FLOW_FUNCTIONS>* outer, NodeManifold<dim>* md, double time)
: outer_(outer),
  manifold_(md),
  model_time_(time)
{ 

}


/*
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Compensator::operator()(const double& sn) const
{
    // masterNode compenastion fluxes
    double sum(0);
    auto master_node = manifold_->N(0);
    assert(master_node != nullptr);
    double master_sn = master_node->Read(outer_->key_sCO2);
    double variation_rate = master_node->Read( outer_->key_dsnw );//m3/(m3 s)
    double PV = master_node->Read( outer_->key_fvPV ); //m3
    sum += variation_rate * PV; //m3/s

    ArrayVariable array;
    master_node->Read(outer_->key_time, array);
    double time_increment = model_time_ - array[7];
    assert(time_increment != 0);
    sum -= (master_sn - sn) * PV / time_increment; // compensate value (m3/s)

    // slaveNode compensation fluxes
    for (size_t n(1); n < manifold_->Branches(); n++) {
        auto slave_node = manifold_->N(n);
        int is_breakthough = static_cast<int>(slave_node->Read(outer_->key_breakthrough));
        if(is_breakthough == 1) {
	          Element<dim>* master_parent = outer_->ParentElementOfManifoldNode(master_node, outer_->key_pd);
	          double master_pc = outer_->flowfunctions_.pc_at(master_parent, 1.0 - sn);
            Element<dim>* slave_parent = outer_->ParentElementOfManifoldNode(slave_node, outer_->key_pd);
            double slave_pd = slave_parent->Read(outer_->key_pd);
            double slave_sw_at_master_pc = outer_->flowfunctions_.seff_to_sw(slave_parent, 1.);
            if(master_pc >= slave_pd) {
                slave_sw_at_master_pc = outer_->flowfunctions_.sw_from_pc_at(slave_parent, master_pc, 1. - slave_node->Read(outer_->key_sCO2));
                //slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1 - slave_parent->Read(outer_->key_srCO2));
	              //slave_sw_at_master_pc = max(slave_sw_at_master_pc, slave_parent->Read(outer_->key_srH2O));
                slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1.0);
                slave_sw_at_master_pc = max(slave_sw_at_master_pc, 0.0);
            }
            double slave_sn = 1. - slave_sw_at_master_pc;
            variation_rate = slave_node->Read( outer_->key_dsnw );//m3/(m3 s)
            PV = slave_node->Read( outer_->key_fvPV ); //m3
            sum += variation_rate * PV; //m3/s
            //slave_node->Read(outer_->key_time, array);
            //time_increment = model_time_ - array[7];
	          sum += (slave_sn - slave_node->Read(outer_->key_sCO2)) * PV / time_increment;
        }
    }

    return sum;    
}
*/


  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Compensator::operator()(const double& sn) const
  {
    // masterNode compenastion fluxes
    double sum(0);
    auto master_node = manifold_->N(0);
    assert(master_node != nullptr);
    //double master_sn = master_node->Read(outer_->key_sCO2);
    double variation_rate = outer_->ComputeRateofChange(master_node, sn);
    double PV = master_node->Read( outer_->key_fvPV ); //m3
    sum += variation_rate * PV; //m3/s

    // slaveNode compensation fluxes
    for (size_t n(1); n < manifold_->Branches(); n++) {
      auto slave_node = manifold_->N(n);
      int is_breakthough = static_cast<int>(slave_node->Read(outer_->key_breakthrough));
      if(is_breakthough == 1) {
        Element<dim>* master_parent = outer_->ParentElementOfManifoldNode(master_node, outer_->key_pd);
        double master_pc = outer_->flowfunctions_.pc_at(master_parent, 1.0 - sn);
        Element<dim>* slave_parent = outer_->ParentElementOfManifoldNode(slave_node, outer_->key_pd);
        double slave_pd = slave_parent->Read(outer_->key_pd);
        double slave_sw_at_master_pc = outer_->flowfunctions_.seff_to_sw(slave_parent, 1.);
        if(master_pc >= slave_pd) {
          slave_sw_at_master_pc = outer_->flowfunctions_.sw_from_pc_at(slave_parent, master_pc, 1. - slave_node->Read(outer_->key_sCO2));
          //slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1 - slave_parent->Read(outer_->key_srCO2));
          //slave_sw_at_master_pc = max(slave_sw_at_master_pc, slave_parent->Read(outer_->key_srH2O));

          //slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1.0);
          //slave_sw_at_master_pc = max(slave_sw_at_master_pc, 0.0);

          if(slave_sw_at_master_pc > 1.0 or slave_sw_at_master_pc < 0.0) {
            double ori_sw = slave_sw_at_master_pc;
            //slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1 - slave_e->Read(this->key_srCO2));
            //slave_sw_at_master_pc = max(slave_sw_at_master_pc, slave_e->Read(this->key_srH2O));
            slave_sw_at_master_pc = min(slave_sw_at_master_pc, 1.0);
            slave_sw_at_master_pc = max(slave_sw_at_master_pc, 0.0);
            cerr<<"    slave_sw = "<<ori_sw<<", adjusted to "<<slave_sw_at_master_pc<<endl;
          }
        }
        double slave_sn = 1. - slave_sw_at_master_pc;
        variation_rate = outer_->ComputeRateofChange(slave_node, slave_sn);
        PV = slave_node->Read( outer_->key_fvPV ); //m3
        sum += variation_rate * PV; //m3/s
      }
    }

    return sum;
  }


  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeFluxBalance( Node<dim>* nd )
  {
    //some variables to use
    VectorVariable<dim> facetNrml, gravity, gradP;
    vector<double> IPOL, NRST;

    const size_t node_parent_elements(nd->Parents());
    long truncated_node = static_cast<long>(nd->Read(this->key_cut));//check if node is truncated by domain boundary

    double flux_balance(0.);

    for ( auto t{0U}; t<node_parent_elements; t++ )
    {
      Element<dim>* const eptr(nd->Parent(t));
      assert( eptr != nullptr );

      if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) //ignore if parent element located outside domain
        continue;

      //read in total velocity computed previously
      VectorVariable<dim> e_vt(ANY, 0.);
      eptr->Read(this->key_vt, e_vt);

      const auto pnid(nd->ParentNodeNumber(t));
      const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
      for ( auto i{0U}; i<sector_facets; i++ )
      {
        const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
        const auto inside_node(eptr->FV()->InsideNode(iFacet));
//        const auto outside_node(eptr->FV()->OutsideNode(iFacet));
        eptr->Read( iFacet, 0U,  this->key_fn, facetNrml );
        const double facetArea = eptr->Read( iFacet, 0U,  this->key_fA );
        const double sign = ( pnid == inside_node ) ? 1. : -1.;

        //compute facet fluid flux
        const double  vD_n = e_vt.DotProduct(facetNrml); //m/s
        //update flux balance
        double f_total = sign * vD_n * facetArea;
        flux_balance += f_total;

      } //end sector_facets loop

    } //end parent element loop

    return flux_balance;

  }




  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Node<dim>* nd, double input_sn )
  {
// not used:    double PV = nd->Read(this->key_fvPV); //pore volume (m3)
    this->rate_count_++;//recording
    nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );

    //store input sn value to the node
    double ori_sn = nd->Read(this->key_sCO2);
    nd->Store( this->key_sCO2, makeScalar( nd->Status(this->key_sCO2), input_sn ) );
    nd->Store( this->key_sH2O, makeScalar( nd->Status(this->key_sH2O), 1.0-input_sn ) );

    //some variables to use
    double carb_accumulation(0.), aq_accumulation(0.);
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> facetNrml, gravity, gradP;
    vector<double> IPOL, NRST;

    const auto node_parent_elements(nd->Parents());
    long truncated_node = static_cast<long>(nd->Read(this->key_cut));//check if node is truncated by domain boundary

    for ( auto t{0U}; t<node_parent_elements; t++ )
    {
      Element<dim>* const eptr(nd->Parent(t));
      assert( eptr != nullptr );

      if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) //ignore if parent element located outside domain
        continue;

      //element state and properties
      const size_t nodes(eptr->Nodes());
      eptr->N_AtBaryCenter( IPOL );
      double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.);
      for ( auto i{0U}; i<nodes; ++i ) {
        e_sw += IPOL[i] * eptr->N(i)->Read( this->key_sH2O );
        e_muw += IPOL[i] * eptr->N(i)->Read( this->key_muH2O );
        e_mun += IPOL[i] * eptr->N(i)->Read( this->key_muCO2 );
        e_rhow += IPOL[i] * eptr->N(i)->Read( this->key_rhoH2O );
        e_rhon += IPOL[i] * eptr->N(i)->Read( this->key_rhoCO2 );
        ipol_sum += IPOL[i];
      }
      if(ipol_sum > 0.) {
        e_sw *= 1. / ipol_sum;
        e_muw *= 1. / ipol_sum;
        e_mun *= 1. / ipol_sum;
        e_rhow *= 1. / ipol_sum;
        e_rhon *= 1. / ipol_sum;
      }
      e_sn = 1.0 - e_sw;

      const double e_swr = eptr->Read(key_srH2O);
      const double e_snr = eptr->Read(key_srCO2);

      eptr->Read(this->key_gradP, gradP); //pressure gradient
      double thickness = eptr->Read(this->key_thi); //thickness
      double e_k = eptr->Read(this->key_k); //permeability
//      double e_phi = eptr->Read(this->key_phi); //porosity

      //read in total velocity computed previously
      VectorVariable<dim> e_vt(ANY, 0.);
      eptr->Read(this->key_vt, e_vt);

      double carb_inflow(0.);
      const auto pnid(nd->ParentNodeNumber(t));
      const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
      for ( auto i{0U}; i<sector_facets; i++ )
      {
        const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
        const auto inside_node(eptr->FV()->InsideNode(iFacet));
        const auto outside_node(eptr->FV()->OutsideNode(iFacet));

        eptr->Read( iFacet, 0U,  this->key_fn, facetNrml );
        const double facetArea = eptr->Read( iFacet, 0U,  this->key_fA );
        const double sign = ( pnid == inside_node ) ? 1. : -1.;

        //compute facet fluid flux
        const double  vD_n = e_vt.DotProduct(facetNrml); //m/s
        //update flux balance
// not used:        double f_total = sign * vD_n * facetArea;

        //compute facet variables
        eptr->N_AtFacetIntegrationPoint( iFacet, 0U, NRST );
        double f_sw(0.), f_muw (0.), f_mun(0.), f_rhow(0.), f_rhon(0.);
        for ( auto x{0U}; x<eptr->Nodes(); x++ ) {
            f_sw += NRST[x] * eptr->N(x)->Read( this->key_sH2O );
            f_muw += NRST[x] * eptr->N(x)->Read( this->key_muH2O );
            f_mun += NRST[x] * eptr->N(x)->Read( this->key_muCO2 );
            f_rhow += NRST[x] * eptr->N(x)->Read( this->key_rhoH2O );
            f_rhon += NRST[x] * eptr->N(x)->Read( this->key_rhoCO2 );
          }

        //compute upstream directions at facet integration point
        double vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
        double vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );

        if( this->with_gravity_forces_ ) {
          eptr->Read( key_dip, gravity );
          if(isnan(gravity(v))) { //dip vector has not been initialised
            if(dim==1u) {gravity(0u) = -1.;}
            else if(dim==2u) {gravity(0u) = 0.; gravity(1u) = -1.;}
            else {gravity(0u) = 0.; gravity(1u) = -1.; gravity(2u) = 0.;}
            eptr->Store( key_dip, gravity );
          }

          if(!this->tensor_k_) { //scalar k
            gravity *=  e_k * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
          } else { //tensor k
            TensorVariable<dim> kk;
            eptr->Read( key_kk, kk );
            VectorVariable<dim> kV;
            kV = kk * gravity;
            double kV_magnitude = kV.Length();
            gravity *= kV_magnitude * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
          }
          double gravity_nrml = gravity.DotProduct(facetNrml);
          vn_gravity_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * gravity_nrml;
          vw_gravity_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * gravity_nrml;
        }

        VectorVariable<dim> sn_gradient;
        sn_gradient = 0.;
        if(this->with_capillary_spreading_) {
          DenseMatrix<DM_MIN> DN;
          eptr->dN_AtBaryCenter(DN);
          for ( auto j{0U}; j < eptr->Nodes(); j++) {
            const double sn = eptr->N(j)->Read(this->key_sCO2);
            for ( auto k{0U}; k < dim; k++) sn_gradient(k) += DN(k, j) * sn;
          }
        }

        if(this->with_capillary_spreading_) {
          //eptr->Read(this->key_gradSn, grad);
          VectorVariable<dim> grad;
          grad = sn_gradient;

          if(!this->tensor_k_) { //scalar k
            grad *= ( -1.0* e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
          } else { //tensor k
            TensorVariable<dim> kk;
            eptr->Read( key_kk, kk );
            grad = kk * grad;
            grad *= ( -1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
          }
          double grad_nrml = grad.DotProduct(facetNrml);
          vn_capillary_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * grad_nrml;
          vw_capillary_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * grad_nrml;
        }

        double vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
        double vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity - vw_capillary_component_of_velocity;

        double lw(0.), ln(0.);
        //aqueous phase determination
        if(vw_at_facet_int_point != 0.0) {
          const auto upstream_node = (vw_at_facet_int_point  > 0.) ? inside_node : outside_node;
          double sw = eptr->N(upstream_node)->Read(this->key_sH2O);
          if(sw > e_swr) lw = this->flowfunctions_.krw_at(eptr, sw) / f_muw;
        } else { //vw_at_facet_int_point == 0.0
          double sw_inside = eptr->N(inside_node)->Read(this->key_sH2O);
          double lw_inside (0.);
          if(sw_inside > e_swr) lw_inside = this->flowfunctions_.krw_at(eptr, sw_inside) / f_muw;

          double sw_outside = eptr->N(outside_node)->Read(this->key_sH2O);
          double lw_outside (0.);
          if(sw_outside > e_swr) lw_outside = this->flowfunctions_.krw_at(eptr, sw_outside) / f_muw;

          lw = (lw_inside + lw_outside) * 0.5;
        }

        //carbonic phase determination
        if(vn_at_facet_int_point != 0.0) {
          const auto upstream_node = (vn_at_facet_int_point  > 0.) ? inside_node : outside_node;
          double sw = eptr->N(upstream_node)->Read(this->key_sH2O);
          double sn = eptr->N(upstream_node)->Read(this->key_sCO2);
          if(sn > e_snr) ln = this->flowfunctions_.krn_at(eptr, sw) / f_mun;
        } else { //vn_at_facet_int_point == 0.0
          double sw_inside = eptr->N(inside_node)->Read(this->key_sH2O);
          double sn_inside = eptr->N(inside_node)->Read(this->key_sCO2);
          double ln_inside (0.);
          if(sn_inside > e_snr) ln_inside = this->flowfunctions_.krn_at(eptr, sw_inside) / f_mun;

          double sw_outside = eptr->N(outside_node)->Read(this->key_sH2O);
          double sn_outside = eptr->N(outside_node)->Read(this->key_sCO2);
          double ln_outside (0.);
          if(sn_outside > e_snr) ln_outside = this->flowfunctions_.krn_at(eptr, sw_outside) / f_mun;

          ln = (ln_inside + ln_outside) * 0.5;
        }

        //fractional flows determination
        double total_mobility = lw + ln;
        double upstream_fn=(total_mobility!=0.0 ? ln/total_mobility : 0.0);
        double upstream_lambda_overbar=(total_mobility!=0.0 ? (ln*lw)/total_mobility : 0.0);
        //viscous flux
        double viscous_velocity_component = vD_n * upstream_fn;
        //gravitational flux
        double gravity_velocity_component(0.);
        if( this->with_gravity_forces_ ) {
          double gravity_nrml = gravity.DotProduct(facetNrml);
          gravity_velocity_component = upstream_lambda_overbar * gravity_nrml;
        }
        //capillary flux
        double capillary_velocity_component(0.);
        if( this->with_capillary_spreading_ ) {
          VectorVariable<dim> grad;
          //eptr->Read(this->key_gradSn, grad);
          grad = sn_gradient;
          if(!this->tensor_k_) { //scalar k
            grad *= ( -1.0* e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
          } else { //tensor k
            TensorVariable<dim> kk;
            eptr->Read( key_kk, kk );
            grad = kk * grad;
            grad *= ( -1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
          }
          double grad_nrml = grad.DotProduct(facetNrml);
          capillary_velocity_component = upstream_lambda_overbar * grad_nrml;
        }

        //total facet non-wetting phase flux
        double f_n = sign * (viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;

        //total facet wetting phase flux
        double upstream_fw=(total_mobility!=0.0? lw/total_mobility : 0.0);
        double viscous_velocity_component_w = vD_n * upstream_fw;
        double f_w = sign * (viscous_velocity_component_w + gravity_velocity_component - capillary_velocity_component) * facetArea;

        //update phase accumulation
        carb_accumulation += f_n;
        carb_inflow += f_n;
        aq_accumulation += f_w;

      } //end sector_facets loop

      /*
      if (!this->no_flow_boundary_) {
          //inflow/outflow compensation for truncated boundary node
          if (truncated_node == 1) {
              //if (inflow > 0.) {flux_balance += inflow; outflow += inflow;} //inflow compensation
              //else if (inflow < 0.) {flux_balance -= inflow;}; //outflow compensation
              if (inflow > 0.) {flux_balance += inflow;} //inflow compensation
              else if (inflow < 0.) {flux_balance -= inflow; outflow -= inflow;}; //outflow compensation

              if ( carb_inflow > 0. ) carb_accumulation += carb_inflow; //inflow compensation
              else if ( carb_inflow < 0. ) carb_accumulation -= carb_inflow; //outflow compensation
          }
      }
      */

    } //end parent element loop

    //restore original saturations to node
    nd->Store( this->key_sCO2, makeScalar( nd->Status(this->key_sCO2), ori_sn ) );
    nd->Store( this->key_sH2O, makeScalar( nd->Status(this->key_sH2O), 1.0-ori_sn ) );

    return carb_accumulation/nd->Read(  this->key_fvPV );

  }



  template <class Function>
double brent_solve(Function& func, const double x1, const double x2, const double tol, const int ITMAX) 
{
    bool verbose(false);
    size_t iteration(0);
	  constexpr double  EPS = std::numeric_limits<double>::epsilon();
	  double  a = x1, b = x2, c = x2, d{ std::numeric_limits<double>::quiet_NaN() },
            e{0.}, fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;
	  if (abs(fa) < EPS) return a;
	  if (abs(fb) < EPS) return b;

    /*			
	  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		    //a = 0.000 - 1e-6;
		    //b = 1.000 + 1e-6;
		    a -= 1e-6;
		    b += 1e-6;
	    	fa = func(a);
		    fb = func(b);
        if (abs(fa) < EPS) return a;
        if (abs(fb) < EPS) return b;
	  }
    */
        
	  fc = fb;
	  for (int iter = 0; iter < ITMAX; iter++) {
		    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			      c = a;
            fc = fa;
			      e = d = b - a;
		    }
	    	if (abs(fc) < abs(fb)) {
			      a = b;
			      b = c;
			      c = a;
			      fa = fb;
			      fb = fc;
			      fc = fa;
		    }
		    tol1 = 2.0*EPS*abs(b) + 0.5*tol;
		    xm = 0.5*(c - b);
		    if (abs(xm) <= tol1 || abs(fb) < tol) {
		      	return b;
		    }
	    	if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
		      	s = fb / fa;
		      	if (fabs(a - c) < tol) {
				        p = 2.0*xm*s;
				        q = 1.0 - s;
		      	}
			      else {
				        q = fa / fc;
				        r = fb / fc;
				        p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
				        q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			      }
			      if (p > 0.0) q = -q;
			      p = abs(p);
			      if (2 * p < std::min(3 * xm*q - fabs(tol1*q), fabs(e*q))) {
				        e = d;
				        d = p / q;
			      }
			      else {
				        d = xm;
				        e = d;
			      }
		    }
		    else {
			      d = xm;
			      e = d;
		    }
		    a = b;
		    fa = fb;
		    if (abs(d) > tol1)
			      b += d;
		    else
			      b += std::copysignf(tol1, xm);
		    fb = func(b);
		
		    iteration ++;

		    if(iter == ITMAX -1) cout<<"    max brent iterations ( "<< ITMAX <<" ) reached, solution may not be converged"<<endl;
    }

	  if(verbose) cout<<"   brent solve finished in "<<iteration<<" iterations"<<endl;

	  return b;

}



template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeLocalFluidPressure( double time, double dt, bool transient )
{
    bool verbose (false);

    NodeList_.clear();
    std::set<NodeManifold<dim>*> manifolds;
    std::set<Node<dim>*> manifolds_nodes;
    //insert all nodes in PEPList to NodeList_
    //insert all manifold nodes in PEPList to manifolds set
    for(size_t i = 0U; i < this->PEPList.size(); ++i)
    {
        auto it = this->PEPList.begin()+i;
        auto nd = (*it)->getNode();
        NodeList_.push_back(nd);
        auto md = nd->Manifold();
        if(md != nullptr) manifolds.insert(md);
    }

    //insert all nodes in manifolds set to manifolds_nodes set
    for(auto md : manifolds) {
        for(size_t n=0;n<md->Branches();n++) {
            auto node = md->N(n);
            manifolds_nodes.insert(node);
        }
    }

    //as nodes with DIRICH sCO2 conditions are excluded from DES, they need to be added back to NodeList_
    size_t missing_nodes(0);
    for(auto nd : manifolds_nodes) {
        if(std::find(NodeList_.begin(), NodeList_.end(), nd) != NodeList_.end()) {
            if(nd->Status(  this->key_sCO2 ) != DIRICH) {
                NodeList_.push_back(nd);
                missing_nodes++;
            }
        }
    }
    if(missing_nodes > 0) cout<<"    added "<<missing_nodes<<" missingg nodes to NimbleRegion"<<endl; //added
 
    //create an active nimble region from NodeList_
    auto active_region = new NimbleRegion<dim>( this->db_, NodeList_.begin(), NodeList_.end() );
    //compute flow properties in nimble region
    for ( auto eit = active_region->CellsBegin(); eit != active_region->CellsEnd(); ++eit )
        ComputeFlowPropertiesAtBaryCenter(*(*eit));

    //firstly store original pressure status for perimeter nodes in nimble region
    //then set their status to DIRICH
    for ( auto nit = active_region->PerimeterNodesBegin(); nit != active_region->NodesEnd(); ++nit ) {
        VARIABLE_FLAG flag = (*nit)->Status(this->key_pf);
        int status(0);
        if ( flag == PLAIN ) status = 0;
        else if ( flag == INIT_GUESS ) status = 1;
        else if ( flag == INIT_COND ) status = 2;
        else if ( flag == DIRICH ) status = 3;
        else if ( flag == NEUMANN ) status = 4;
        else if ( flag == ROBIN ) status = 5;
        else if ( flag == PERIODIC ) status = 6;
        else if ( flag == FIELD_DATA ) status = 7;
        else if ( flag == ANY ) status = 8;
        (*nit)->Store( key_status, makeScalar( (*nit)->Status( key_status), status ) );
        //set to Dirichlet
        (*nit)->Status( this->key_pf, DIRICH );
    }

    //solve transient or steady state pressure in nimble region
    if(transient) SolveNimbleRegionPressureFullyImplicit( *active_region, dt, time );
    else SolveNimbleRegionSteadyStatePressure(*active_region);

    //check whether solved pressure is in reasonable range
    //if so, update pressure gradents and flow velocities
    double min_pf_computed, max_pf_computed;
    double min_pf(0.), max_pf(6e9);
    this->sg_.MinMaxOf( "fluid pressure", min_pf_computed, max_pf_computed );
    if(min_pf_computed<min_pf || max_pf_computed>max_pf) {
        cout<<"warning: computed fluid pressures ("<<min_pf_computed<<" - "<<max_pf_computed<<") out of range, local pressure solving skipped"<<endl;
        this->sg_.CopyReplace("previous fluid pressure","fluid pressure"); //reset reduced pressure
    }
    else {
        for(size_t i = 0U; i < this->PEPList.size(); ++i) {
            auto it = this->PEPList.begin()+i;
            ComputePressureGradientAndFlowVelocities(*it);
        }
    }

    //restore original pressure status of perimeter nodes in nimble region
    for ( auto nit = active_region->PerimeterNodesBegin(); nit != active_region->NodesEnd(); ++nit ) {
        int status = (*nit)->Read(this->key_status);
        assert(status>=0);
        assert(status<=8);
        if ( status == 0 ) (*nit)->Status( this->key_pf, PLAIN );
        else if ( status == 1 ) (*nit)->Status( this->key_pf, INIT_GUESS );
        else if ( status == 2 ) (*nit)->Status( this->key_pf, INIT_COND );
        else if ( status == 3 ) (*nit)->Status( this->key_pf, DIRICH );
        else if ( status == 4 ) (*nit)->Status( this->key_pf, NEUMANN );
        else if ( status == 5 ) (*nit)->Status( this->key_pf, ROBIN );
        else if ( status == 6 ) (*nit)->Status( this->key_pf, PERIODIC );
        else if ( status == 7 ) (*nit)->Status( this->key_pf, FIELD_DATA );
        else if ( status == 8 ) (*nit)->Status( this->key_pf, ANY );
    }

    //update pressure solving count
    for ( auto nit = active_region->NodesBegin(); nit != active_region->NodesEnd(); ++nit )
        (*nit)->Store( key_p_count, makeScalar( (*nit)->Status( key_p_count), (*nit)->Read(key_p_count) + 1 ) );

    if(verbose) {
        cout << "  Local P solved, Nimble region input N = " << NodeList_.size()
             << " interior N = " << active_region->InteriorNodes()
             << " perimeter N = " << active_region->PerimeterNodes()
             << " total N = " << active_region->Nodes()
             << " total E = " << active_region->Cells() << endl;
    }
        
    delete active_region;
}




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeFlowPropertiesAtBaryCenter(Element<dim>& e)
{
    const size_t v( (dim==1u) ? 0u : 1u );

    static Index  cw_key(this->sg_.Database().StorageKey("compressibility aqueous phase"));
    static Index  cn_key(this->sg_.Database().StorageKey("compressibility carbonic phase"));
    static Index  cR_key(this->sg_.Database().StorageKey("compressibility rock"));
    static Index  ct_key(this->sg_.Database().StorageKey("total system compressibility"));
    static Index  gt_key(this->sg_.Database().StorageKey("gravity term"));

    // bulk interpolating all values for the flow property calculations
    const size_t nodes(e.Nodes());
    vector<double> IPOL;
    e.N_AtBaryCenter( IPOL );
    double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.), e_cw(0.), e_cn(0.);
    for ( auto i{0U}; i<nodes; ++i ) {
        e_sw += IPOL[i] * e.N(i)->Read( this->key_sH2O );
        e_muw += IPOL[i] * e.N(i)->Read( this->key_muH2O );
	      e_mun += IPOL[i] * e.N(i)->Read( this->key_muCO2 );
	      e_rhow += IPOL[i] * e.N(i)->Read( this->key_rhoH2O );
	      e_rhon += IPOL[i] * e.N(i)->Read( this->key_rhoCO2 );
        e_cw += IPOL[i] * e.N(i)->Read( cw_key );
        e_cn += IPOL[i] * e.N(i)->Read( cn_key );			
        ipol_sum += IPOL[i];
    }
    if(ipol_sum > 0.) {
	      e_sw *= 1. / ipol_sum;
	      e_muw *= 1. / ipol_sum;
	      e_mun *= 1. / ipol_sum;
	      e_rhow *= 1. / ipol_sum;
	      e_rhon *= 1. / ipol_sum;
	      e_cw *= 1. / ipol_sum;
	      e_cn *= 1. / ipol_sum;
    }
    e_sn = 1.0 - e_sw;

    const double e_swr = e.Read(this->key_srH2O);
    const double e_snr = e.Read(this->key_srCO2); 
    double thickness = e.Read(this->key_thi); //thickness
    if(isnan(thickness)) thickness = 1.;
    double e_phi = e.Read(this->key_phi); //porosity

    if(!this->tensor_k_) { //scalar k
        static Index  mobt_key(this->sg_.Database().StorageKey("total mobility permeability product"));
        double e_k = e.Read(this->key_k); //permeability
        //phase mobilities
        double e_lw(0.), e_ln(0.), e_lt(0.);
        if(e_sw > e_swr) e_lw = this->flowfunctions_.krw_at(&e, e_sw) / e_muw; //krw/muw
        if(e_sn > e_snr) e_ln = this->flowfunctions_.krn_at(&e, e_sw) / e_mun; //krn/mun
        //total mobility permeability product
        e_lt = (e_lw + e_ln) *  e_k * thickness;
        e.Store( mobt_key, makeScalar( e.Status( mobt_key ), e_lt ) );

        if( this->with_gravity_forces_ ) {
            VectorVariable<dim> gravity;
            e.Read( this->key_dip, gravity );
            if(isnan(gravity(1))) { //dip vector has not been initialised
                gravity(0u) = 0.; gravity(1u) = -1.;
                if(dim==3u) gravity(2u) = 0.;
                e.Store( key_dip, gravity );
            }
            double gravity_w(0.), gravity_n(0.);
            if(e_sw > e_swr) gravity_w = this->sg_.Read(key_g) * e_k * thickness * this->flowfunctions_.krw_at(&e, e_sw) / e_muw * (e_rhow);
            if(e_sn > e_snr) gravity_n = this->sg_.Read(key_g) * e_k * thickness * this->flowfunctions_.krn_at(&e, e_sw) / e_mun * (e_rhon);
            gravity_w *= fabs(gravity[v]);
            gravity_n *= fabs(gravity[v]); 
            gravity *= (gravity_w + gravity_n);
            e.Store( gt_key, gravity );        
        }     
         
        // weighted compressibility
        const double e_cR = e.Read( cR_key );
        double e_ct = (1. - e_phi) * e_cR + e_phi * (e_sw * e_cw + (1.-e_sw) * e_cn);
        e_ct *= thickness;    
        e.Store( ct_key, makeScalar( e.Status( ct_key ), e_ct ) );

    } else { //tensor_k
        static Index  kk_key(this->sg_.Database().StorageKey("tensor permeability"));
        static Index  LT_key(this->sg_.Database().StorageKey("tensor total mobility permeability product"));
        //compute phase mobilities
        double e_lw(0.), e_ln(0.), e_lt(0.);
        if(e_sw > e_swr) e_lw = thickness * this->flowfunctions_.krw_at(&e, e_sw) / e_muw; //krw/muw
        if(e_sn > e_snr) e_ln = thickness * this->flowfunctions_.krn_at(&e, e_sw) / e_mun; //krn/mun
        e_lt = e_lw + e_ln;
        TensorVariable<dim> kk;
        e.Read( kk_key, kk );
        kk *= e_lt;
        e.Store( LT_key, kk );

        if( this->with_gravity_forces_ ) {
            VectorVariable<dim> gravity;
            e.Read( key_dip, gravity );
            if(isnan(gravity(1))) { //dip vector has not been initialised
                gravity(0u) = 0.; gravity(1u) = -1.;
                if(dim==3u) gravity(2u) = 0.;
                e.Store( key_dip, gravity );
            }
            // projection of tensor on the dip vector
            e.Read( key_kk, kk );
            VectorVariable<dim> kV;
            kV = kk * gravity;
            double kV_magnitude = kV.Length();
            double gravity_w(0.), gravity_n(0.);
            if(e_sw > e_swr) gravity_w = this->sg_.Read(key_g) * kV_magnitude * thickness * this->flowfunctions_.krw_at(&e, e_sw) / e_muw * (e_rhow);
            if(e_sn > e_snr) gravity_n = this->sg_.Read(key_g) * kV_magnitude * thickness * this->flowfunctions_.krn_at(&e, e_sw) / e_mun * (e_rhon); 
            gravity_w *= fabs(gravity[v]);
            gravity_n *= fabs(gravity[v]);      
            gravity *= (gravity_w + gravity_n);      
            e.Store( gt_key, gravity );
        }  
            
        // weighted compressibility
        const double e_cR = e.Read( cR_key );
        double e_ct = (1. - e_phi) * e_cR + e_phi * (e_sw * e_cw + (1.-e_sw) * e_cn);
        e_ct *= thickness;    
        e.Store( ct_key, makeScalar( e.Status( ct_key ), e_ct ) );
    }

 
    // -----------------------------------------------------------------------------------------
    // 5. capillary diffusivity (capillary spreading can occur only where 2 phases are present)
    // -----------------------------------------------------------------------------------------
    // uses mobility product which is zero at single phase conditions
    // the capillary diffusion multiplier automatically uses the saturation at the element barycentre
    //double lambda_overbar_rho = l_rho_H2O * l_rho_CO2 / lt_rho;
    //double diffusivity = std::max( this->flowfunctions_.dpcds_at(&e,sw) * k * lambda_overbar_rho, min_val_for_diffusion );
    //assert( !isnan(diffusivity) );
    //e.Store( this->key_diffpc, makeScalar( e.Status( this->key_diffpc ), diffusivity ) );
   
 } // end ComputeFlowPropertiesAtBaryCenter




//solve pressure on entire model - used in TDS scheme
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeSteadyStatePressure()
{
    bool verbose(false);

    std::string	conductance_operator;
    if(!this->tensor_k_) conductance_operator = "total mobility permeability product";
    else conductance_operator = "tensor total mobility permeability product";

    if(verbose) {
      cout << "\nDES2PhaseSlightlyCompressibleTransport::ComputeSteadyStatePressure: Input parameters: " << endl;
      printRangeOfVariable( this->sg_, "fluid pressure" );
      printRangeOfVariable( this->sg_, conductance_operator.c_str() );
      printRangeOfVariable( this->sg_, "permeability" );
      printRangeOfVariable( this->sg_, "fluid volume source");
      //printRangeOfVariable( this->sg_, "nodal fluid volume source" );
      printRangeOfVariable( this->sg_, "boundary influx" );
      if(this->with_gravity_forces_ ) printRangeOfVariable( this->sg_, "gravity term" );
    }

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                  samg_solver( &settings );
    PDE_Integrator<dim,Element>  steady_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER linear_solver;
    PDE_Integrator<dim,Element>  steady_pressure(linear_solver);
#endif
    NumIntegral_dNT_op_dN_dV<dim>  conductance( this->sg_.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<dim>    elmt_volume_source( this->sg_.Database(), "fluid volume source", "fluid pressure" );
    //PointSource_rhsop<dim>          nodal_volume_source( this->sg_.Database(), "nodal fluid volume source", "fluid pressure" );
    //NumIntegral_NT_op_N_dS<dim,Face>     influx( this->sg_.Database(), "boundary influx", "fluid pressure" );
    //influx.LumpedFormulation(true);
    steady_pressure.Add( &conductance );
    steady_pressure.Add( &elmt_volume_source );
    //steady_pressure.Add( &nodal_mass_source );
    //steady_pressure.AddBoundaryIntegral( &influx );
    //nodal_mass_source.AddAccumulate();
    steady_pressure.Verbose(verbose);

    NumIntegral_dNT_op_dV<dim>* gravity(nullptr);
    if ( this->with_gravity_forces_ ) {
        gravity = new NumIntegral_dNT_op_dV<dim>( this->sg_.Database(), "gravity term", "fluid pressure" );
        steady_pressure.Add( gravity );
      }

    if(verbose) {
        cout << "\n\nDES2PhaseSlightlyCompressibleTransport:ComputeSteadyStatePressure: ";
        cout << " Computing '" << "steady state fluid pressure" << "'" << endl;
      }

    this->sg_.Apply( steady_pressure );

    delete gravity;

} // end ComputeSteadyStatePressure



//solve steady state pressure on nimble region
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::SolveNimbleRegionSteadyStatePressure( NimbleRegion<dim>& computation_domain)
{
    bool verbose(false);

    std::string	conductance_operator;
    if(!this->tensor_k_) conductance_operator = "total mobility permeability product";
    else conductance_operator = "tensor total mobility permeability product";

    if(verbose) {
      cout << "\nDES2PhaseSlightlyCompressibleTransport::SolveNimbleRegionSteadyStatePressure: Input parameters: " << endl;
      printRangeOfVariable( this->sg_, "fluid pressure" );
      printRangeOfVariable( this->sg_, conductance_operator.c_str() );
      printRangeOfVariable( this->sg_, "permeability" );
      printRangeOfVariable( this->sg_, "fluid volume source");
      //printRangeOfVariable( this->sg_, "nodal fluid volume source" );
      //printRangeOfVariable( this->sg_, "boundary influx" );
      if(this->with_gravity_forces_ ) printRangeOfVariable( this->sg_, "gravity term" );
    }

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                       samg_solver( &settings );
    // PDE_Integrator<dim,NimbleRegion>  steady_pressure(samg_solver);
    PDE_Integrator<dim,Element>  steady_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER        linear_solver;
    // PDE_Integrator<dim,NimbleRegion>  steady_pressure(linear_solver);
    PDE_Integrator<dim,Element>  steady_pressure(linear_solver);
#endif
    NumIntegral_dNT_op_dN_dV<dim>   conductance( this->sg_.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<dim>     elmt_volume_source( this->sg_.Database(), "fluid volume source", "fluid pressure" );
    //PointSource_rhsop<dim>          nodal_volume_source( this->sg_.Database(), "nodal fluid volume source", "fluid pressure" );
    //NumIntegral_NT_op_N_dS<dim,Face>     influx( this->sg_.Database(), "boundary influx", "fluid pressure" );
    //influx.LumpedFormulation(true);
    steady_pressure.Add( &conductance );
    steady_pressure.Add( &elmt_volume_source );
    //steady_pressure.Add( &nodal_volume_source );
    //steady_pressure.AddBoundaryIntegral( &influx );
    //nodal_volume_source.AddAccumulate();
    steady_pressure.Verbose(verbose);

    NumIntegral_dNT_op_dV<dim>* gravity(nullptr);
    if ( this->with_gravity_forces_ ) {
         gravity = new NumIntegral_dNT_op_dV<dim>( this->sg_.Database(), "gravity term", "fluid pressure" );
         steady_pressure.Add( gravity );
      }

    if (verbose) {
        cout << "\n\nDES2PhaseSlightlyCompressibleTransport:SolveNimbleRegionSteadyStatePressure: ";
        cout << " Computing '" << "steady state fluid pressure" << "'" << endl;
      }

    steady_pressure.IntegrateOver(computation_domain, false);

    delete gravity;

  } // end SolveNimbleRegionSteadyStatePressure




//solve transient pressure on nimble region
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::SolveNimbleRegionPressureFullyImplicit( NimbleRegion<dim>& computation_domain, double time_increment, double time )
{
    bool verbose(false);

    std::string	conductance_operator;
    if(!this->tensor_k_) conductance_operator = "total mobility permeability product";
    else conductance_operator = "tensor total mobility permeability product";

    if(verbose) {
        cout << "\nDES2PhaseSlightlyCompressibleTransport::SolveNimbleRegionPressureFullyImplicit: Input parameters: " << endl;
        printRangeOfVariable( this->sg_, "fluid pressure" );
        printRangeOfVariable( this->sg_, conductance_operator.c_str() );
        printRangeOfVariable( this->sg_, "permeability" );
        printRangeOfVariable( this->sg_, "total system compressibility" );
        printRangeOfVariable( this->sg_, "fluid volume source");
        //printRangeOfVariable( this->sg_, "nodal fluid volume source" );
        //printRangeOfVariable( this->sg_, "boundary influx" );
        if(this->with_gravity_forces_ ) printRangeOfVariable( this->sg_, "gravity term" );
    }

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                       samg_solver( &settings );
    PDE_Integrator<dim,Element>  transient_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER        linear_solver;
//    PDE_Integrator<dim,NimbleRegion>  transient_pressure(linear_solver);
    PDE_Integrator<dim,Element>  transient_pressure(linear_solver);
#endif
    NumIntegral_dNT_op_dN_dV<dim> conductance( this->sg_.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_lhsop_N_dV<dim> capacitance_lhs( this->sg_.Database(), "total system compressibility", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<dim> capacitance_rhs( this->sg_.Database(), "total system compressibility", "fluid pressure" );
    NumIntegral_NT_op_N_dV<dim> elmt_volume_source( this->sg_.Database(), "fluid volume source", "fluid pressure" );
    //PointSource_rhsop<dim> nodal_volume_source( this->sg_.Database(), "nodal fluid volume source", "fluid pressure" );

    capacitance_lhs.MultiplyWithTimeIncrement( true );
    capacitance_rhs.MultiplyWithTimeIncrement( true );
    capacitance_lhs.LumpedFormulation( true );
    //capacitance_rhs.LumpedFormulation( true );
    elmt_volume_source.AddAccumulateLater();
    //nodal_volume_source.AddAccumulateLater();

    transient_pressure.Verbose(verbose);
    transient_pressure.Add( &conductance );
    transient_pressure.Add( &capacitance_lhs );
    transient_pressure.Add( &capacitance_rhs );
    transient_pressure.Add( &elmt_volume_source );
    //reduced_pressure.Add( &nodal_volume_source );

    // optional inclusion of the gravity term
    NumIntegral_dNT_op_dV<dim>* gravity(nullptr);
    if ( this->with_gravity_forces_) {
          gravity = new NumIntegral_dNT_op_dV<dim>( this->sg_.Database(), "gravity term", "fluid pressure" );
          gravity->AddAccumulateLater();
          transient_pressure.Add( gravity );
      }

    transient_pressure.TimeIncrement( 1. / time_increment );

    if(verbose) {
        cout << "\n\nDES2PhaseSlightlyCompressibleTransport::SolveNimbleRegionPressureFullyImplicit: ";
        cout << " Computing '" << "transient fluid pressure" << "'" << endl;
      }

    transient_pressure.IntegrateOver(computation_domain, false);

    delete gravity;

} // SolveNimbleRegionPressureFullyImplicit





//Find predicted time level based on smallest scheduled time stamp - used in DES scheme
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::TimeLevel()
{
    this->ResetCFLMultiplier();
    //compute pressure gradients       
    size_t PEPList_size = this->PEPList.size();
    for(size_t i = 0U; i < PEPList_size; ++i)
    {
        auto it = this->PEPList.begin()+i;
        Event<dim>* event = *it;
        ComputePressureGradientAndFlowVelocities (event);
        if(this->with_capillary_spreading_) ComputeSaturationGradient (event);
    }
    
    //compute rate of change
#if defined(OPENMP)   
    double T_begin = omp_get_wtime();
    size_t num_threads = omp_get_max_threads();
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(dynamic)
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            Event<dim>* event = *it;
            if(second_order_in_space_) ComputeRateofChange_2nd_order(event);
            else ComputeRateofChange(event);
        };
    }
    this->T_RateOfChange_ += omp_get_wtime() - T_begin; 
#else
    clock_t T_begin = clock();
    for(size_t i = 0U; i < PEPList_size; ++i)
    {
        auto it = this->PEPList.begin()+i;
        Event<dim>* event = *it;
        if(second_order_in_space_) ComputeRateofChange_2nd_order(event);
        else ComputeRateofChange(event);
    };
    this->T_RateOfChange_ += clock() - T_begin; 
#endif    

    //find smallest time level
    double min_t (numeric_limits<double>::max());
    for(size_t i = 0U; i < PEPList_size; ++i)
    {
        auto it = this->PEPList.begin()+i;
        Event<dim>* event = *it;
        if ((*it)->valid() == false) {    
            Node<dim>* nd = event->getNode();
            if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH){
                ArrayVariable array;
                nd->Read(this->key_time, array);
                double dt_target = array[2]*array[6];//CFL time increment
                double t_current = array[0];//current time stamp
                double t_schedule = t_current + dt_target; //scheduled time stamp
                min_t = min(min_t, t_schedule);
            }
        }
    }
  
    return min_t;
}




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange_2nd_order( Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != nullptr );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    if(nd  != nullptr && nd->Status( this->key_sCO2 ) != DIRICH) {
        double PV = nd->Read(this->key_fvPV); //pore volume (m3)
        if(PV < numeric_limits<double>::epsilon())
            throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange_2nd_order:",
                                   "FV pore volume seems to be zero, has it been initialised yet?" );
        else
        {
            this->rate_count_++;//recording
            nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );

            //some variables to use
            double flux_balance(0.), outflow(0.), tot_inflow(0.), tot_outflow(0.), carb_accumulation(0.), aq_accumulation(0.);
            const uint32_t v( (dim==1u) ? 0u : 1u );
            VectorVariable<dim> facetNrml, gravity, gradP;
            vector<double> IPOL, NRST;

            const size_t node_parent_elements(nd->Parents());
            double cfl_multiplier = this->CFL_multiplier_*this->relaxing_factor_; //default value
            long truncated_node = static_cast<long>(nd->Read(this->key_cut));//check if node is truncated by domain boundary

            for ( auto t{0U}; t<node_parent_elements; t++ )
            {
                Element<dim>* const eptr(nd->Parent(t));
                assert( eptr != nullptr );
      
                if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) //ignore if parent element located outside domain
                    continue;
        
                //element state and properties
                const auto nodes{ eptr->Nodes() };
                eptr->N_AtBaryCenter( IPOL );
                double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.);
                for ( auto i{0U}; i<nodes; ++i ) {
                    e_sw += IPOL[i] * eptr->N(i)->Read( this->key_sH2O );
		                e_muw += IPOL[i] * eptr->N(i)->Read( this->key_muH2O );
		                e_mun += IPOL[i] * eptr->N(i)->Read( this->key_muCO2 );
		                e_rhow += IPOL[i] * eptr->N(i)->Read( this->key_rhoH2O );
		                e_rhon += IPOL[i] * eptr->N(i)->Read( this->key_rhoCO2 );
                    ipol_sum += IPOL[i];
                }
                if(ipol_sum > 0.) {
		                e_sw *= 1. / ipol_sum;
		                e_muw *= 1. / ipol_sum;
		                e_mun *= 1. / ipol_sum;
		                e_rhow *= 1. / ipol_sum;
		                e_rhon *= 1. / ipol_sum;
                }
                e_sn = 1.0 - e_sw;

                const double e_swr = eptr->Read(key_srH2O);
                const double e_snr = eptr->Read(key_srCO2);

        
                eptr->Read(this->key_gradP, gradP); //pressure gradient
                double thickness = eptr->Read(this->key_thi); //thickness
                double e_k = eptr->Read(this->key_k); //permeability
// not used:                double e_phi = eptr->Read(this->key_phi); //porosity

                //read in total velocity computed previously
                VectorVariable<dim> e_vt(ANY, 0.);
                eptr->Read(this->key_vt, e_vt);

                double carb_inflow(0.);
                const auto pnid(nd->ParentNodeNumber(t));
                const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
                for ( auto i{0U}; i<sector_facets; i++ )
                {
                    const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                    const auto inside_node(eptr->FV()->InsideNode(iFacet));
                    const auto outside_node(eptr->FV()->OutsideNode(iFacet));
            
                    eptr->Read( iFacet, 0U,  this->key_fn, facetNrml );
                    const double facetArea = eptr->Read( iFacet, 0U,  this->key_fA );
                    const double sign = ( pnid == inside_node ) ? 1. : -1.;

                    //compute facet fluid flux
                    const double  vD_n = e_vt.DotProduct(facetNrml); //m/s
                    //update flux balance
                    double f_total = sign * vD_n * facetArea;
                    if(f_total != 0.) {
                        flux_balance += f_total;
                        //if(f_total < 0.) outflow -= f_total;
                        if(f_total > 0.) outflow += f_total;
                        //inflow += f_total;
                    }

                    //facet properties
                    eptr->N_AtFacetIntegrationPoint( iFacet, 0U );
                    double f_sw(0.), f_sn(0.), f_muw (0.), f_mun(0.), f_rhow(0.), f_rhon(0.);
                    for ( auto x{0U}; x<eptr->Nodes(); x++ ) {
                        f_sw += eptr->FE()->NRST[x] * eptr->N(x)->Read( this->key_sH2O );
                        f_muw += eptr->FE()->NRST[x] * eptr->N(x)->Read( this->key_muH2O );
                        f_mun += eptr->FE()->NRST[x] * eptr->N(x)->Read( this->key_muCO2 );
                        f_rhow += eptr->FE()->NRST[x] * eptr->N(x)->Read( this->key_rhoH2O );
                        f_rhon += eptr->FE()->NRST[x] * eptr->N(x)->Read( this->key_rhoCO2 );
                      }
                    f_sn = 1.0 - f_sw;

                    //compute upstream directions at facet integration point
                    double vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
                    double vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
                    if( this->with_gravity_forces_ ) {
                          eptr->Read( key_dip, gravity );
                          if(isnan(gravity(v))) { //dip vector has not been initialised
                              if(dim==1u) {gravity(0u) = -1.;}
                              else if(dim==2u) {gravity(0u) = 0.; gravity(1u) = -1.;}
                              else {gravity(0u) = 0.; gravity(1u) = -1.; gravity(2u) = 0.;}
                              eptr->Store( key_dip, gravity );
                          }
                          gravity *=  e_k * thickness * this->sg_.Read(key_g) * (f_rhow - f_rhon); //positive
                          double gravity_nrml = gravity.DotProduct(facetNrml);
                          vn_gravity_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * gravity_nrml;
                          vw_gravity_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * gravity_nrml;
                    }
            
                    if(this->with_capillary_spreading_){
                        VectorVariable<dim> grad;
                        eptr->Read(this->key_gradSn, grad);

                        grad *= ( -1.0* e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                        double grad_nrml = grad.DotProduct(facetNrml);
                        vn_capillary_component_of_velocity = this->flowfunctions_.krw_at(eptr, f_sw) / f_muw * grad_nrml;
                        vw_capillary_component_of_velocity = this->flowfunctions_.krn_at(eptr, f_sw) / f_mun * grad_nrml;
                    }
            
                    double vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
                    double vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity - vw_capillary_component_of_velocity;

                    //compute limited values for saturations
                    double sn_inside = eptr->N(inside_node)->Read(this->key_sCO2);
                    double sn_outside = eptr->N(outside_node)->Read(this->key_sCO2);
                    double limited_sn_inside(0.0), limited_sn_outside(0.0);

                    //LSM method
                    //LimitProperty_LSMGRAD( *eptr, "saturation carbonic phase", inside_node, outside_node, iFacet, sn_inside, sn_outside, limited_sn_inside, limited_sn_outside );
            
                    //MINMOD method
                    std::pair<double,double> minmax_inside, minmax_outside;
                    CalculateMinMax (eptr->N(inside_node), "saturation carbonic phase", minmax_inside);
                    CalculateMinMax (eptr->N(outside_node), "saturation carbonic phase", minmax_outside);
                    limited_sn_inside = LimitProperty( sn_inside, sn_outside, f_sn, minmax_inside );
                    limited_sn_outside = LimitProperty( sn_outside, sn_inside, f_sn, minmax_outside );
            
                    //1st order for boundary nodes
                    //if(eptr->N(inside_node)->AtBoundary() != NOT) limited_sn_inside = sn_inside;
                    //if(eptr->N(outside_node)->AtBoundary() != NOT) limited_sn_outside = sn_outside;

                    double total_f(0.), lw(0.), ln(0.);
                    //aqueous phase determination
                    if(vw_at_facet_int_point != 0.0) {
                        double sw(0.);
                        if (vw_at_facet_int_point  > 0.) sw = 1.0 - limited_sn_inside;
                        else sw = 1.0 - limited_sn_outside;
                        if(sw > e_swr) lw = this->flowfunctions_.krw_at(eptr, sw) / f_muw;
                    } else { //vw_at_facet_int_point == 0.0
                        double sw_inside = 1.0 - limited_sn_inside;
                        double lw_inside (0.);
                        if(sw_inside > e_swr) lw_inside = this->flowfunctions_.krw_at(eptr, sw_inside) / f_muw;
                        double sw_outside = 1.0 - limited_sn_outside;
                        double lw_outside (0.);
                        if(sw_outside > e_swr) lw_outside = this->flowfunctions_.krw_at(eptr, sw_outside) / f_muw;
                        lw = (lw_inside + lw_outside) * 0.5;
                    }
            
                    //carbonic phase determination
                    if(vn_at_facet_int_point != 0.0) {
                        double sw(0.), sn(0.);
                        if (vn_at_facet_int_point  > 0.) sn = limited_sn_inside;
                        else sn = limited_sn_outside;
                        sw = 1.0 - sn;
                        if(sn > e_snr) ln = this->flowfunctions_.krn_at(eptr, sw) / f_mun;
                    } else { //vn_at_facet_int_point == 0.0
                        double inside_sn = limited_sn_inside;
                        double inside_sw = 1.0 - inside_sn;
                        double inside_ln (0.);
                        if(inside_sn > e_snr) inside_ln = this->flowfunctions_.krn_at(eptr, inside_sw) / f_mun;
                
                        double outside_sn = limited_sn_outside;
                        double outside_sw = 1.0 - outside_sn;
                        double outside_ln (0.);
                        if(outside_sn > e_snr) outside_ln = this->flowfunctions_.krn_at(eptr, outside_sw) / f_mun;
                                
                        ln = (inside_ln + outside_ln) * 0.5;
                    }

                    //fractional flows determination
                    double total_mobility = lw + ln;
                    double upstream_fn=(total_mobility!=0.0? ln/total_mobility : 0.0);
                    double upstream_lambda_overbar=(total_mobility!=0.0 ? (ln*lw)/total_mobility : 0.0);
                    //viscous flux
                    double viscous_velocity_component = vD_n * upstream_fn;
                    //gravitational flux
                    double gravity_velocity_component(0.);
                    if( this->with_gravity_forces_ ) {
                        //gravity already computed previously
                        double gravity_nrml = gravity.DotProduct(facetNrml);
                        gravity_velocity_component = upstream_lambda_overbar * gravity_nrml;
                    }
                    //capillary flux
                    double capillary_velocity_component(0.);
                    if( this->with_capillary_spreading_ ) {
                        VectorVariable<dim> grad;
                        eptr->Read(this->key_gradSn, grad);
                        if(!this->tensor_k_) { //scalar k
                            grad *= ( -1.0* e_k * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                        } else { //tensor k
                            TensorVariable<dim> kk;
                            eptr->Read( key_kk, kk );
                            grad = kk * grad;
                            grad *= ( -1.0 * thickness * this->flowfunctions_.dpcds_at(eptr, f_sw) );
                        }
                        double grad_nrml = grad.DotProduct(facetNrml);
                        capillary_velocity_component = upstream_lambda_overbar * grad_nrml;
                    }

                    //total facet non-wetting phase flux
                    double f_n = sign * (viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;

                    //total facet wetting phase flux
                    double upstream_fw=(total_mobility!=0.0? lw/total_mobility : 0.0);
                    double viscous_velocity_component_w = vD_n * upstream_fw;
                    double f_w = sign * (viscous_velocity_component_w + gravity_velocity_component - capillary_velocity_component) * facetArea;

                    //update total outflow and inflow
                    total_f = f_n + f_w;
                    if(total_f < 0.) tot_inflow -= total_f;
                    if(total_f > 0.) tot_outflow += total_f;

                    //update phase accumulation
                    carb_accumulation += f_n;
                    carb_inflow += f_n;
                    aq_accumulation += f_w;

                    //determine cfl_multiplier based on non-wetting phase shock saturation
                    if (cfl_multiplier != this->CFL_multiplier_){
                        if(vn_at_facet_int_point< 0.0) { //flowing in from outside node (upstream node)
                            double sn_shock = 1.0-eptr->Read(this->key_ssH2O); //sn at shock for outside node
                            double sn_outside_node = eptr->N(outside_node)->Read( this->key_sCO2 );
                            if (sn_outside_node >= sn_shock) { //upstream node passed shock saturation
                                double sn_inside_node = eptr->N(inside_node)->Read( this->key_sCO2 );
                                if (sn_inside_node < sn_shock) {//current node not yet reach shock saturation
                                    cfl_multiplier = this->CFL_multiplier_;
                                }
                            }
                        }
                    }
        
                    //determine cfl_multiplier based on wetting phase shock saturation
                    if (cfl_multiplier != this->CFL_multiplier_){
                        if(vw_at_facet_int_point< 0.0) { //flowing from outside node (upstream node)
                            double sw_shock = eptr->Read(this->key_ssH2O); //sw at shock for outside node
                            double sw_outside_node = eptr->N(outside_node)->Read( this->key_sH2O );
                            if (sw_outside_node >= sw_shock) { //upstream node passed shock saturation
                                double sw_inside_node = eptr->N(inside_node)->Read( this->key_sH2O );
                                if (sw_inside_node < sw_shock) {//current node not yet reach shock saturation
                                    cfl_multiplier = this->CFL_multiplier_;
                                }
                            }
                        }
                    }
                } //end sector_facets loop

                /*
                if (!this->no_flow_boundary_) {
                    //inflow/outflow compensation for truncated boundary node
                    if (truncated_node == 1) {
                        //if (inflow > 0.) {flux_balance += inflow; outflow += inflow;} //inflow compensation
                        //else if (inflow < 0.) {flux_balance -= inflow;}; //outflow compensation
                        if (inflow > 0.) {flux_balance += inflow;} //inflow compensation
                        else if (inflow < 0.) {flux_balance -= inflow; outflow -= inflow;}; //outflow compensation
                
                        if ( carb_inflow > 0. ) carb_accumulation += carb_inflow; //inflow compensation
                        else if ( carb_inflow < 0. ) carb_accumulation -= carb_inflow; //outflow compensation
                    }
                }
                */
        
            } //end parent element loop

            /*
            //compensate flow for boundary node
            if( nd->AtBoundary() != NOT && nd->Status(this->key_pf) == DIRICH ) { //added
               //double outflux = ComputeOutflowAtBoundary(event); //added
               double outflux = carb_accumulation; //added
               nd->Store( key_outflow, makeScalar( nd->Status( key_outflow ), outflux ) ); //added
               //carb_accumulation -= outflux; //added
            } //added
            */

            /*
            if(with_periodic_boundary_) {
                for (auto node_pair : this->periodic_nodes_) {
                    if(nd == node_pair.second) { //outlet boundary node
                        double outflux = carb_accumulation; //added
                        nd->Store( key_outflow, makeScalar( nd->Status( key_outflow ), outflux ) ); //added
                        carb_accumulation -= outflux; //added
                        break;
                    }
                }
            }
            */

            nd->Store( key_compensate, makeScalar( nd->Status( key_compensate ), carb_accumulation/nd->Read(this->key_fvPV) ) );

            /*
            if(fabs(flux_balance - 0.) > 1.0e-12 )
                cout<<"  local volume not conserved, event = "<<nd->Read(this->key_EventIndex)<<" flux_balance = "<<flux_balance<<endl;
        
            // divergence free correction (only when node is located inside domain and flux balance not equal to zero)
            if (nd->Manifold() != nullptr && truncated_node !=1 && fabs(flux_balance) > numeric_limits<double>::epsilon()) {
                // compute average fractional flow for the current finite volume
                double fw_avg(0), fn_avg(0.);
                double n_sn = nd->Read(this->key_sCO2);
                double n_sw = nd->Read(this->key_sH2O); //saturation aqueous phase at current node
                for ( size_t t=0U; t<node_parent_elements; t++ )
                {
                    Element<dim>* const eptr(nd->Parent(t));
                    double ep_swr = eptr->Read(this->key_srH2O);
                    double ep_snr = eptr->Read(this->key_srCO2);

                    const size_t nodes(eptr->Nodes());
                    eptr->N_AtBaryCenter( IPOL );
                    double ipol_sum(0.), ep_sw(0.), ep_sn(0.), ep_muw (0.), ep_mun(0.), ep_rhow(0.), ep_rhon(0.);
                    for ( size_t i=0U; i<nodes; ++i ) {
                        ep_sw += IPOL[i] * eptr->N(i)->Read( this->key_sH2O );
		                    ep_muw += IPOL[i] * eptr->N(i)->Read( this->key_muH2O );
		                    ep_mun += IPOL[i] * eptr->N(i)->Read( this->key_muCO2 );
		                    ep_rhow += IPOL[i] * eptr->N(i)->Read( this->key_rhoH2O );
		                    ep_rhon += IPOL[i] * eptr->N(i)->Read( this->key_rhoCO2 );
                        ipol_sum += IPOL[i];
                    }
                    if(ipol_sum > 0.) {
		                    ep_sw *= 1. / ipol_sum;
		                    ep_muw *= 1. / ipol_sum;
		                    ep_mun *= 1. / ipol_sum;
		                    ep_rhow *= 1. / ipol_sum;
		                    ep_rhon *= 1. / ipol_sum;
                    }
                    ep_sn = 1.0 - ep_sw;

                    double lambda_w(0.), lambda_n(0.);
                    if(n_sw > ep_swr) lambda_w = this->flowfunctions_.krw_at(eptr, n_sw) / ep_muw;
                    if(n_sn > ep_snr) lambda_n = this->flowfunctions_.krn_at(eptr, n_sw) / ep_mun;
                    double lambda_t = lambda_w + lambda_n;
                    double ep_fw = (lambda_t!=0.0 ? lambda_w/lambda_t : 0.0);
                    double ep_fn = (lambda_t!=0.0 ? lambda_n/lambda_t : 0.0);
                    //fw_avg += ep_fw;
                    fn_avg += ep_fn;
                }
                //fw_avg /= static_cast<double>(node_parent_elements);
                fn_avg /= static_cast<double>(node_parent_elements);
        
                //carb_accumulation -= (fn_avg*tot_accumulation); //added
                //aq_accumulation -= (fw_avg*tot_accumulation); //added

                carb_accumulation -= (fn_avg*flux_balance); //added
                //aq_accumulation -= (fw_avg*flux_balance); //added

                //double new_tot_accumulation = carb_accumulation + aq_accumulation;
                //if(fabs(new_tot_accumulation - 0.) > 1.0e-6)
                //if(fabs(new_tot_accumulation - 0.) > numeric_limits<double>::epsilon())
                  //cout<<"  local volume still not conserved after correction, event = "<<nd->Read(this->key_EventIndex)<<" original divergence  = "<<tot_accumulation<<" corrected divergence = "<<new_tot_accumulation<<endl;
            }
             */
    
            //store rate of changes
            nd->Store( key_dsnw, makeScalar( nd->Status( key_dsnw ), carb_accumulation/nd->Read(  this->key_fvPV ) ) );

            //compute CFL time increment
            ArrayVariable array2;
            nd->Read(this->key_time, array2);
            //if (tot_inflow < numeric_limits<double>::epsilon())
            if (tot_outflow < numeric_limits<double>::epsilon())
                array2.Component(2, numeric_limits<double>::max());
            else
                //array2.Component(2, nd->Read(  this->key_fvPV ) / outflow);
                //array2.Component(2, nd->Read(  this->key_fvPV ) / tot_inflow);
                array2.Component(2, nd->Read(  this->key_fvPV ) / tot_outflow);
    
            nd->Store( this->key_CFL, makeScalar( nd->Status( this->key_CFL ),cfl_multiplier ) );
            array2.Component(6, cfl_multiplier);
        
            nd->Store(this->key_time, array2);
    
        }  //end if(PV > numeric_limits<double>::epsilon())
    } //end if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH)
}





template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::LimitProperty_LSMGRAD(  const Element<dim>& e,
                                                                                         const char* prop,
                                                                                         uint32_t inside_node, uint32_t outside_node, uint32_t iFacet,
                                                                                         const double prop_inside_node, const double prop_outside_node,
                                                                                         double& limited_prop_inside_node, 
                                                                                         double& limited_prop_outside_node)
{
    // Calculate the coordinates of the Facet Integration Point
    const Point<dim> local_c_point(e.FV()->FacetIntegrationPoint( iFacet, 0U ));
    std::vector<double> temp(e.Nodes()); //has the local interp. function values
    std::vector<double> global_c(dim),local_c(local_c_point.Coordinates());

    if( e.IsLine() ){
        e.FE()->Nr( local_c[0], temp );
    }else if( e.IsSurface()){
        e.FE()->Nrs( local_c[0], local_c[1], temp );
    }else{
        e.FE()->Nrst( local_c[0], local_c[1], local_c[2], temp );
    }

    e.CoordinateMatrix();

    global_c.assign( dim, 0.0);

    // transform local c's to global c's
    for ( uint32_t m = 0; m<e.Nodes(); m++)
        for ( uint32_t n = 0; n<dim; n++){
          global_c[n] += e.FE()->XY(m,n)*temp[m];
    }
    
    VectorVariable<dim> mass_center_inside_node(PLAIN,0.0);
    VectorVariable<dim> mass_center_outside_node(PLAIN,0.0);
    e.N(inside_node)->Read( key_mc, mass_center_inside_node);
    e.N(outside_node)->Read( key_mc, mass_center_outside_node);
    
    VectorVariable<dim> distance_inside_node(PLAIN,0.0);
    VectorVariable<dim> distance_outside_node(PLAIN,0.0);
    for( uint32_t l=0U;l<dim;l++) {
        distance_inside_node.Component(l,global_c[l]-mass_center_inside_node[l]);
        distance_outside_node.Component(l,global_c[l]-mass_center_outside_node[l]);
    }    

    Node<dim>* in_node = e.N(inside_node);
    Node<dim>* out_node = e.N(outside_node);
    
    VectorVariable<dim> grad_prop_inside_node(PLAIN,0.0);
    VectorVariable<dim> grad_prop_outside_node(PLAIN,0.0);
         
    CalculateGenericNodalGradient (in_node, prop, grad_prop_inside_node);
    CalculateGenericNodalGradient (out_node, prop, grad_prop_outside_node);
    
    double limiter_prop_inside_node = CalculateSlopeLimiter(in_node, prop, grad_prop_inside_node);
    double limiter_prop_outside_node = CalculateSlopeLimiter(out_node, prop, grad_prop_outside_node);
    //double limiter_prop_inside_node = 1.0;
    //double limiter_prop_outside_node = 1.0;

    double prop_linear_increment_inside_node=0.0;
    double prop_linear_increment_outside_node=0.0;
    for( uint32_t l=0U;l<dim;l++) {
        prop_linear_increment_inside_node+=grad_prop_inside_node[l]*distance_inside_node[l];
        prop_linear_increment_outside_node+=grad_prop_outside_node[l]*distance_outside_node[l];  
    }

    limited_prop_inside_node = prop_inside_node + limiter_prop_inside_node * prop_linear_increment_inside_node;
    limited_prop_outside_node = prop_outside_node + limiter_prop_outside_node * prop_linear_increment_outside_node;

}




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::CalculateGenericNodalGradient (Node<dim>* nd, const char* prop, VectorVariable<dim>& grad)
{
    grad = 0.0;
    csmp::Index u_key = this->sg_.Database().StorageKey( prop );
    VARIABLE_FLAG status = nd->Status(u_key);
    
    // gradient is zero if nodes are flagged DIRICH or NEUMANN
    if ( status == PLAIN || status == ANY) {
        double tolerance(1.0e-12);
        VectorVariable<dim> xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0), dcxyz(PLAIN,0.0);
        double val1, val2, dc, det;
        
        nd->Read(key_mc, xyz1);
        val1 = nd->Read(u_key );
        
        if(dim == 1U) {
            double sum_x2(0.);
            for ( uint32_t i=0U; i<nd->Neighbors(); i++ ) {
                nd->Neighbor(i)->Read(key_mc, xyz2);
                dxyz.Component(0, xyz2[0]-xyz1[0]);
                sum_x2 += dxyz[0] * dxyz[0];
        
                val2 = nd->Neighbor(i)->Read(u_key );
                dc = (val2 -val1);    
                dcxyz.Component(0, dcxyz[0]+dc*dxyz[0]);    
            }
        
            // calculate determinante
            det = sum_x2;
            if(det!=0.0){
                grad(0) = dcxyz[0]/det;
                if(fabs(grad(0)) < tolerance) grad(0) = 0.0;    
            } else { 
                grad = 0.0;
            }
            
        } else if (dim == 2U) {
            double sum_x2(0.), sum_y2(0.), sum_xy(0.);
            for ( uint32_t i=0U; i<nd->Neighbors(); i++ ) {
                nd->Neighbor(i)->Read(key_mc, xyz2);
                dxyz.Component(0, xyz2[0]-xyz1[0]);
                dxyz.Component(1, xyz2[1]-xyz1[1]);
                sum_x2 += dxyz[0] * dxyz[0];
                sum_y2 += dxyz[1] * dxyz[1];
                sum_xy += dxyz[0] * dxyz[1];                
        
                val2 = nd->Neighbor(i)->Read(u_key );
                dc = (val2 -val1);    
                dcxyz.Component(0, dcxyz[0]+dc*dxyz[0]); 
                dcxyz.Component(1, dcxyz[1]+dc*dxyz[1]);   
            }
        
            // calculate determinante
            det = sum_x2*sum_y2 - sum_xy*sum_xy;
            
            size_t zero_grad_index = 6U;
            if(det == 0.0){
                if(sum_x2 == 0.0)
                    zero_grad_index = 0U;
                else if(sum_y2 == 0.0)
                    zero_grad_index = 1U;
                else
                    zero_grad_index = 6U;
            }            
            
            if(det != 0.0){
                grad(0)  = (dcxyz[0] * sum_y2 - dcxyz[1] * sum_xy) / det;
                grad(1)  = (dcxyz[1] * sum_x2 - dcxyz[0] * sum_xy) / det;
                if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
            }else if(zero_grad_index == 0U){
                grad(0) = 0.0;
                grad(1) = dcxyz[1]/ sum_y2;
                if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
            }else if(zero_grad_index == 1U){
                grad(0) = dcxyz[0]/ sum_x2;
                grad(1) = 0.0;
                if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
            }else {
                grad = 0.0;
            }  
            
        } else if (dim == 3U) {   
            double sum_x2(0.), sum_y2(0.), sum_z2(0.), sum_xy(0.), sum_xz(0.), sum_yz(0.);
            for ( uint32_t i=0U; i<nd->Neighbors(); i++ ) {
                nd->Neighbor(i)->Read(key_mc, xyz2);
                dxyz.Component(0, xyz2[0]-xyz1[0]);
                dxyz.Component(1, xyz2[1]-xyz1[1]);
                dxyz.Component(2, xyz2[2]-xyz1[2]);
                sum_x2 += dxyz[0] * dxyz[0];
                sum_y2 += dxyz[1] * dxyz[1];
                sum_z2 += dxyz[2] * dxyz[2];
                sum_xy += dxyz[0] * dxyz[1];    
                sum_xz += dxyz[0] * dxyz[2];
                sum_yz += dxyz[1] * dxyz[2];                            
        
                val2 = nd->Neighbor(i)->Read(u_key );
                dc = (val2 -val1);    
                dcxyz.Component(0, dcxyz[0]+dc*dxyz[0]); 
                dcxyz.Component(1, dcxyz[1]+dc*dxyz[1]);   
                dcxyz.Component(2, dcxyz[2]+dc*dxyz[2]); 
            }
        
            // calculate determinante
            det = (sum_x2*sum_y2 - sum_xy*sum_xy) * sum_z2;
            det -= (sum_x2*sum_yz - sum_xy*sum_xz) * sum_yz;
            det += (sum_xy*sum_yz - sum_y2*sum_xz) * sum_xz;
        
            // the indexes of zero gradient components
            //0 - 0, 1- 1, 2 - 2, 0&1 - 3, 0&2 - 4, 1&2 - 5, else 6
            size_t zero_grad_index = 6U;
            if(det == 0.0){
                if(sum_x2 == 0.0){
                    if(sum_y2 == 0.0)
                        zero_grad_index = 3U;
                    else if(sum_z2 == 0.0)
                        zero_grad_index = 4U;
                    else
                        zero_grad_index = 0U;
                } else if(sum_y2 == 0.0){
                    if(sum_x2 == 0.0)
                        zero_grad_index = 3U;
                    else if(sum_z2 == 0.0)
                        zero_grad_index = 5U;
                    else
                        zero_grad_index = 1U;
                } else if(sum_z2 == 0.0){
                    if(sum_x2 == 0.0)
                        zero_grad_index = 4U;
                    else if(sum_y2 == 0.0)
                        zero_grad_index = 5U;
                    else
                        zero_grad_index = 2U;
                } else {
                    zero_grad_index = 6U;
                }
            }          
            
            if(det != 0.0){
                grad(0) = ( dcxyz[0]*sum_y2 - dcxyz[1]*sum_xy ) * sum_z2;
                grad(0) -= ( dcxyz[0]*sum_yz - dcxyz[1]*sum_xz ) *sum_yz;
                grad(0) += ( sum_xy*sum_yz - sum_y2*sum_xz ) * dcxyz[2];
                grad(0) /= det;

                grad(1) = ( sum_x2*dcxyz[1] - sum_xy*dcxyz[0] ) * sum_z2;
                grad(1) -= ( sum_x2*sum_yz - sum_xy*sum_xz ) * dcxyz[2];
                grad(1) += ( dcxyz[0]*sum_yz - dcxyz[1]*sum_xz ) * sum_xz;
                grad(1) /= det;

                grad(2) = ( sum_x2*sum_y2 - sum_xy*sum_xy ) * dcxyz[2];
                grad(2) -= ( sum_x2*dcxyz[1] - sum_xy*dcxyz[0] ) * sum_yz;
                grad(2) += ( sum_xy*dcxyz[1] - sum_y2*dcxyz[0] ) * sum_xz;
                grad(2) /= det;

                if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
                if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;
                

            }else if(zero_grad_index == 0U){
                grad(0) = 0.0;
                
                grad(1) = ( dcxyz[1]*sum_z2 - dcxyz[2]*sum_yz );
                grad(1) /= ( sum_y2*sum_z2 - sum_yz*sum_yz );

                grad(2) = ( sum_y2*dcxyz[2] - sum_yz*dcxyz[1] );
                grad(2) /= ( sum_y2*sum_z2 - sum_yz*sum_yz );

                if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
                if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

            }else if(zero_grad_index == 1U){
                grad(0) = ( dcxyz[0]*sum_z2 - dcxyz[2]*sum_xz );
                grad(0) /= ( sum_x2*sum_z2 - sum_xz*sum_xz );

                grad(1) = 0.0;

                grad(2) = ( sum_x2*dcxyz[2] - sum_xz*dcxyz[0] );
                grad(2) /= ( sum_x2*sum_z2 - sum_xz*sum_xz );

                if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

            }else if(zero_grad_index == 2U){
                if(( sum_x2*sum_y2 - sum_xy*sum_xy ) != 0.){
                    grad(0) = ( dcxyz[0]*sum_y2 - dcxyz[1]*sum_xy );
                    grad(0) /= ( sum_x2*sum_y2 - sum_xy*sum_xy );
                    
                    grad(1) = ( sum_x2*dcxyz[1] - sum_xy*dcxyz[0] );
                    grad(1) /= ( sum_x2*sum_y2 - sum_xy*sum_xy );
                } else {
                    grad(0) = 0.0;
                    grad(1) = 0.0;
                }

                grad(2) = 0.0;

                if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
                

            }else if(zero_grad_index == 3U){
                grad(0) = 0.0;

                grad(1) = 0.0;

                grad(2) = dcxyz[2]/sum_z2;

                if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

            }else if(zero_grad_index == 4U){
                grad(0) = 0.0;

                grad(1) = dcxyz[1]/sum_y2;

                grad(2) = 0.0;

                if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;

            }else if(zero_grad_index == 5U){
                grad(0) = dcxyz[0]/sum_x2;

                grad(1) = 0.0;

                grad(2) = 0.0;

                if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;

            }else {
                grad = 0.0;
            }  
            
        } else {
            throw csmp::Exception( FATAL_ERROR, "DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::CalculateGenericNodalGradient:",
                "wrong dimension, dim must be 1, 2 or 3" );
        }
        
    } //end if ( status == PLAIN || status == ANY) 

}



template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::CalculateSlopeLimiter(Node<dim>* nd, const char* prop, VectorVariable<dim> grad)
{
    csmp::Index u_key = this->sg_.Database().StorageKey( prop );
    VARIABLE_FLAG status = nd->Status(u_key);

    // if finite volume is at a node with DIRICHLET or NEUMANN boundary conditions,
    // its phi value is zero (first order upwind scheme)
    double phi(0.);
    if ( status == PLAIN || status == ANY) {
        double tolerance(1.0e-12);
        VectorVariable<dim> dist;
        double val1 = nd->Read(u_key);
        phi = 1.0; // set phi to 1.0   
        
        //including current FV value  
        double min = val1;
        double max = val1;
        
        //excluding current FV value
        //double min = numeric_limits<double>::max();
        //double max = numeric_limits<double>::min();
        
        double val2;
        for ( uint32_t i=0U; i<nd->Neighbors(); i++ ) {
            val2 = nd->Neighbor(i)->Read(u_key );
            if ( val2 < min ) min = val2;
            if ( val2 > max ) max = val2;
        }
        
        // sitting at a node, loop over parent el's
        for(  uint32_t p{0U}; p< nd->Parents(); p++ ){
            Element<dim>* current_el = nd->Parent(p);  
            size_t   global_el_id = current_el->Idx(); // get the global parent id:
            uint32_t nloc_id = nd->ParentNodeNumber( p ); // get local node number
            // at that parent element, loop over all facets that belong to the current node/fv
            for ( uint32_t i=0U; i < current_el->FV()->FacetsPerSector(nloc_id); i++ ) {
                uint32_t local_facet_id = current_el->FV()->FacetSurroundingSector( nloc_id,i ); //get local facet_id
                // get distance barycenter - facetcenter for that facet:
                dist = distance_facet_FVBarycenter_[ global_el_id ][ local_facet_id ][ nloc_id ];

                // construct linear interpolant:
                double val_left  = val1;
                for ( uint32_t j=0U;j<dim;j++) val_left += grad[j]*dist[j];
                
                //phi_temp = limitProperty(val_left, val1, min, max);    
                double phi_temp(1.);
                if ( val_left > val1 ) {
                    phi_temp  = max;
                    phi_temp -= val1;
                    phi_temp /= ( val_left - val1 );
                } else if ( val_left < val1 ) {
                    phi_temp  = min;
                    phi_temp -= val1;
                    phi_temp /= ( val_left - val1 );
                }

                // if necessary, set new value to phi
                if ( phi_temp < phi ) phi=phi_temp;
                if ( fabs( phi ) < tolerance ) phi = 0.0;
                if ( essentiallyEqual( phi, 0.0 ) ) break;
            } //end facets
        } //end parent elements
                                                    
    } //end if
    
    return phi;
}



//compute center of mass for all nodes
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::CalculateCenterOfMass()
{
    size_t  glob_n_id;
    double volume_;
    std::vector<double>   current_bc; // stencil barycenter in global coord's
    std::vector<size_t> ids;
    VectorVariable<dim> temp_;
    Point<dim> tmp_p;
 
    std::vector<std::pair<VectorVariable<dim>, double> > center_of_mass;
    center_of_mass.resize(( this->sg_.Region("Model").Nodes() ));
    std::pair<VectorVariable<dim>, double>   init_pair;
    init_pair.first = VectorVariable<dim>( PLAIN, 0.0);
    init_pair.second = 0.0;
    fill( center_of_mass.begin(), center_of_mass.end(), init_pair ); 
 
    // loop over stencils
    // ------------------
    for ( auto eit=this->sg_.Region("Model").CellsBegin();
          eit!=this->sg_.Region("Model").CellsEnd(); eit++ )
    {
        auto e = *(*eit);
        ids.resize( e.Nodes() );
        // get vector of global id's
        for( uint32_t i=0U;i<e.Nodes();i++)
            ids[i]=e.N(i)->Idx();

        // loop over sectors
        // -------------------
        for( uint32_t i=0U; i < e.FV()->Sectors(); i++ ){

              //get the global node id for the current segment
              glob_n_id = ids[ i ];

              if( this->sg_.Region("Model").N( glob_n_id)->AtBoundary() ){
                  tmp_p = this->sg_.Region("Model").N( glob_n_id)->Coordinate();
                  temp_ = VectorVariable<dim>( PLAIN,  0.0);
                  for( uint32_t k=0;k<dim;k++)
                      temp_.Component(k,tmp_p.Coordinates()[k]);
                  center_of_mass[ glob_n_id ].first  = temp_;
                  center_of_mass[ glob_n_id ].second = 1.;

              }else{
                  //get the barycenter of current segment
                  //ConvertToGlobalCoordinates( e, e.FV()->SectorIntegrationPoint( i, 0U), current_bc );
                  const Point<dim>& local_c_point = e.FV()->SectorIntegrationPoint( i, 0U);
                  std::vector<double> temp(e.Nodes()); //has the local interp. function values
                  std::vector<double> local_c(local_c_point.Coordinates());

                  if( e.IsLine() )
                      e.FE()->Nr( local_c[0], temp );
                  else if( e.IsSurface())
                      e.FE()->Nrs( local_c[0], local_c[1], temp );
                  else
                      e.FE()->Nrst( local_c[0], local_c[1], local_c[2], temp );
  
                  e.CoordinateMatrix();
                  current_bc.assign( dim, 0.0);
                  // transform local c's to global c's
                  for ( uint32_t k{0U}; k<e.Nodes(); k++ )
                      for ( uint32_t j{0U}; j<dim; j++ )
                          current_bc[j] += e.FE()->XY(k,j) * temp[k];
 
                  //get volume of current sector
                  volume_ = e.SectorVolume( i );

                  temp_ = VectorVariable<dim>( PLAIN,  0.0);
                  for( uint32_t k=0U;k<dim;k++)
                      temp_.Component(k,current_bc[k] * volume_);

                  center_of_mass[ glob_n_id ].first += temp_;
                  center_of_mass[ glob_n_id ].second += volume_;
             }
         }
     }

     for( size_t i = 0U;  i< this->sg_.Region("Model").Nodes(); i++ )
     {
        //  center_of_mass_[ i ].first.Out();
        //  calculate center: sum_i(x_i * A_i) / sum_i(A_i)
        center_of_mass[ i ].first *=  1. / center_of_mass[ i ].second;

        // write back to node:
        this->sg_.Region("Model").N( i )->Store( key_mc, center_of_mass[ i ].first );
     }

}



template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::CalculateDistanceFacetFVBary()
{
    //initialization of distances facet centers-barycenters
    distance_facet_FVBarycenter_.resize( this->sg_.Region("Model").Cells() ); // resize outer vector
    for ( auto eit=this->sg_.Region("Model").CellsBegin();eit!=this->sg_.Region("Model").CellsEnd(); eit++ ) {
        auto e = *(*eit);
        //resize middle vector
        distance_facet_FVBarycenter_[ e.Idx() ].resize( e.FV()->Facets() );
        // loop over facets
        // -------------------
        for( size_t fi=0U; fi < e.FV()->Facets(); fi++ ){
            // resize inner vector
            distance_facet_FVBarycenter_[ e.Idx() ][ fi ].resize( e.Nodes() );
            fill( distance_facet_FVBarycenter_[ e.Idx() ][ fi ].begin(), distance_facet_FVBarycenter_[ e.Idx() ][ fi ].end(), VectorVariable<dim>( PLAIN, 0.0));
        }
    }    
    
    //calculation of distances facet centers-barycenters
    uint32_t inside_node_, outside_node_; // local node id's:
    std::vector<double>  facet_bc(dim); // global coordinates of facet barycenter
    VectorVariable<dim>  mass_center;
    // loop over elements
    for ( auto eit=this->sg_.Region("Model").CellsBegin(); eit!=this->sg_.Region("Model").CellsEnd(); eit++ ){
        auto e = *(*eit);
        // loop over facets
        for( uint32_t fi=0U; fi < e.FV()->Facets(); fi++ ){
            // get facet barycenter in global coordinates:
            //ConvertToGlobalCoordinates( e, e.FV()->FacetIntegrationPoint( fi, 0U ),  facet_bc);
            const Point<dim>& local_c_point = e.FV()->FacetIntegrationPoint( fi, 0U );
            std::vector<double> temp(e.Nodes()); //has the local interp. function values
            std::vector<double> local_c(local_c_point.Coordinates());

            if( e.IsLine() )
                e.FE()->Nr( local_c[0], temp );
            else if( e.IsSurface())
                e.FE()->Nrs( local_c[0], local_c[1], temp );
            else
                e.FE()->Nrst( local_c[0], local_c[1], local_c[2], temp );
  
            e.CoordinateMatrix();
            facet_bc.assign( dim, 0.0);
            // transform local c's to global c's
            for ( uint32_t i = 0U; i<e.Nodes(); i++)
                for (uint32_t j = 0U; j<dim; j++)
                    facet_bc[j] += e.FE()->XY(i,j) * temp[i];              


            // get local node id's
            e.FV()->FacetEdgeNodes( fi, inside_node_, outside_node_ );

            // get mass center for FV of inside_node_:
            e.N(inside_node_)->Read(key_mc, mass_center );

            VectorVariable<dim> face_dist_inside_node( PLAIN,  0.0);
            for( uint32_t k=0U;k<dim;k++)
                face_dist_inside_node.Component(k,facet_bc[k] - mass_center[k]);

            distance_facet_FVBarycenter_[ e.Idx() ][ fi ][ inside_node_ ] =face_dist_inside_node;

            // get mass center for FV of outside_node_:
            e.N(outside_node_)->Read(key_mc, mass_center );

            // calculate distance vector:
            VectorVariable<dim> face_dist_outside_node( PLAIN,  0.0);
            for( uint32_t k=0U;k<dim;k++)
                face_dist_outside_node.Component(k,facet_bc[k] - mass_center[k]);

            distance_facet_FVBarycenter_[ e.Idx() ][ fi ][ outside_node_ ] = face_dist_outside_node;
        }
    }
}




/*****************************/
//2nd order with MINMOD limiter
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::CalculateMinMax (Node<dim>* nd, const char* prop, std::pair<double,double>& minmax)
{
    csmp::Index u_key = this->sg_.Database().StorageKey( prop );
// not used:    VARIABLE_FLAG status = nd->Status(u_key);

    minmax.first = minmax.second = nd->Read( u_key );
    for ( uint32_t i=0U; i<nd->Neighbors(); i++ ) {
        const double adv_var(nd->Neighbor(i)->Read( u_key ));
        // if element value is smaller the current minimum is assigned etc.
        minmax.first  = std::min( minmax.first,  adv_var );
        minmax.second = std::max( minmax.second, adv_var );
    }

}



/**

Calculates the isotropic limiter value using the MINMOD scheme
for a slope limiter xi with default value 2.

@section arguments Input Arguments

SMINMAX stores min and max values of the advected property in the
neighborhood of the upstream node, including its own value.

@return U_tilde_f at segment.

@section application Application

In the second-order transport methods of the DES2PhaseSlightlyCompressibleTransport class

*/
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::LimitProperty( double psi_hat_c, double psi_hat_d, double psi_facet,
                                                                                  const std::pair<double,double>& SMINMAX_upstr, double xi )
 {
    // isotropic multi-dimensional limiting of facet saturations      min                   max
    const double psi_hat_u = (psi_hat_d > psi_hat_c) ? SMINMAX_upstr.first : SMINMAX_upstr.second;

    // finding min, max values of advected property in neighborhood of upstream node
    double psi_diff = psi_hat_d - psi_hat_u;
    // avoiding division by zero

    if ( std::fabs(psi_diff) < std::numeric_limits<double>::epsilon() ) psi_diff = 1.0e-20 * ( psi_diff < 0.0 ? -1.0 : 1.0);

    // computing U_f, U_c factors for slope limiting (Pain et al. paper)
    // 'psi_dash_f' = fsn (advected property values extrapolated to facet integration point)
    const double U_f = (psi_facet  - psi_hat_u) / psi_diff;
    const double U_c = (psi_hat_c  - psi_hat_u) / psi_diff;

    // computing normalized (always positive) slope-limiter value at face (eqn 38, Pain et al. 2001)
    const double U_tilde_f = NVD_Function( xi, U_f, U_c );

    // finding the slope limited value of the advected property at the segment according to eqn. 36
    // linearly interpolate between two values
    //const double psi_hat = lerp(U_tilde_f, psi_hat_u, psi_hat_d);
    const double psi_hat = (1.0 - U_tilde_f) * psi_hat_u + U_tilde_f * psi_hat_d;

    return psi_hat;
 }




/**

Computes the spatially limited variable value from Xi, U_f the value of the
property at the segment and the upstream value of the advected property
U_c.

@return U_tilde_f at segment.

@section application Application

In the method which calculates the isotropic limiter values.
*/
template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
double DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::NVD_Function( double xi, double U_f, double U_c )
 {
    constexpr double one(1.), zero(0.);
    if ( U_c > one || U_c < zero ) return U_c;
    double  U_tilde_f = fmax( zero, U_f );
    return fmin( fmin(xi * U_c, U_tilde_f), one );
 }




//BrooksCoreySaturationFunctions
template class DES2PhaseSlightlyCompressibleTransport<1U,FlowFunctionsModule1>;
template class DES2PhaseSlightlyCompressibleTransport<2U,FlowFunctionsModule1>;
template class DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule1>;
/*
//BrooksCoreySaturationFunctionsWithHysteresis
template class DES2PhaseSlightlyCompressibleTransport<1U,FlowFunctionsModule2>;
template class DES2PhaseSlightlyCompressibleTransport<2U,FlowFunctionsModule2>;
template class DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule2>;
*/
//ExperimentalSaturationFunctions
template class DES2PhaseSlightlyCompressibleTransport<1U,FlowFunctionsModule3>;
template class DES2PhaseSlightlyCompressibleTransport<2U,FlowFunctionsModule3>;
template class DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule3>;


} // end csmp 



                                                                       
