#include "DES2PhaseSlightlyCompressibleTransport.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#include "FlowFunctionsModule.h"
#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace std;

namespace csmp {


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::DES2PhaseSlightlyCompressibleTransport(  Model<dim>& m,
                                                                                                     const char* target_region,
                                                                                                     bool with_gravity_forces,
                                                                                                     bool with_capillary_spreading,
                                                                                                     double64 cfl_multiplier,
                                                                                                     double64 PEP_multiplier,
                                                                                                     double64 relaxing_factor,
                                                                                                     bool tensor_k,
                                                                                                     FLOW_FUNCTIONS<dim>& ff )
    : DES2PhaseTransport<dim,FLOW_FUNCTIONS>(m,target_region,with_gravity_forces,with_capillary_spreading,cfl_multiplier,PEP_multiplier,relaxing_factor,tensor_k,ff)
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
} // end constructor 




template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::InitializeVariablsAndKeys()
{
    //creating new variables if not defined yet
    if(!this->db_.IsDefined("variation rate nonwetting phase")) this->sg_.CreateProperty( "variation rate nonwetting phase", "m3/(m3.s)", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);
    if(!this->db_.IsDefined("old saturation carbonic phase")) this->sg_.CreateProperty( "old saturation carbonic phase", "m3/m3", SCALAR, NODE, 1, 0 ,1);
    if(!this->db_.IsDefined("tensor permeability")) this->sg_.CreateProperty( "tensor permeability", "m2", TENSOR, ELEMENT, 3, 1E-21, 1.0e-5);
    if(!this->db_.IsDefined("nodal fluid volume source")) {
        this->sg_.CreateProperty( "nodal fluid volume source", "m3", SCALAR, NODE, 1, -1.00E+01, 1.00E+01);
        this->sg_.Region("Model").InputPropertyValue( "nodal fluid volume source", makeScalar(PLAIN,0), COMPLETE);
    }    
    
    //assigning keys 
    key_dsnw = INDEX<SCALAR,NODE> ( this->db_.StorageKey("variation rate nonwetting phase") );
    key_sCO2_0 = INDEX<SCALAR,NODE>( this->db_.StorageKey("old saturation carbonic phase") );
    key_NQV = INDEX<SCALAR,NODE>( this->db_.StorageKey("nodal fluid volume source") );
    key_kk = INDEX<TENSOR,ELEMENT>( this->db_.StorageKey("tensor permeability") );
    
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
        
    // model-wide initialisation
    this->sg_.Region("Model").InputPropertyValue( "variation rate nonwetting phase", makeScalar(PLAIN,0), COMPLETE);
}




template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::InitializeEvents()
{ 
    //create events for all nodes and add them to event lists
    size_t index = 0;
    size_t dirich_count = 0;
    const typename vector<Node<dim>*>::const_iterator  nodes_end(this->gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=this->gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        if((*nit)->Status(  this->key_sCO2 ) != DIRICH) {
            (*nit)->Store( this->key_EventIndex, makeScalar( (*nit)->Status(this->key_EventIndex), index) );//event index 
            Event<dim>* event = new Event<dim>(*nit);
            this->PEPList.push_back(event);
            event->inPEPStack(true);
            ComputeGradients (event );
            ComputeRateofChange(event);
            Schedule(event, 0.);
            event->valid(false);
            Heap_Node* heap_node = new Heap_Node(event->t_schedule(),index);
            this->HeapNodeFullList.push_back(heap_node);
            this->FullList.push_back(event);
            event->inQueue(false);
                
            index++;
        
        } else {
            dirich_count++;
        }
    }
    cout<<this->FullList.size()<<" events created for all nodes, excluding "<<dirich_count<<" DIRICH nodes"<<endl;
}  




//Compute non-wetting phase saturaiton gradient at parement elements, for capillary component computation.
// TODO: super expensive approach - use values from neighboring nodes 
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeGradients (Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  
    //check if node is truncated by domain boundary
    int truncated_node = static_cast<int>(nd->Read(this->key_cut));
  
    const size_t parent_elements(nd->Parents());      
    for ( size_t i=0U; i<parent_elements; ++i ) {
        Element<dim>* const eptr = nd->Parent(i);
        
        if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) { //ignore if parent element located outside domain
            continue;
        
        } else {        
        
            DenseMatrix<DM_MIN> DN;
            eptr->dN_AtBaryCenter(DN);
            VectorVariable<dim> snw_gradient, p_gradient;
            snw_gradient = 0.;
            p_gradient = 0.;     
        
            for ( size_t j=0U; j<eptr->Nodes(); j++ ) {
                const double64 sn = eptr->N(j)->Read(this->key_sCO2);
                const double64 p = eptr->N(j)->Read(this->key_pf);
                //const double64 p = eptr->N(j)->Read(this->this->key_rpf);
                for ( size_t k=0U; k<dim; k++ ) {
                    if(this->with_capillary_spreading_) snw_gradient(k) += DN(k,j) * sn;
                    p_gradient(k) += -DN(k,j) * p;               
                }
            }
        
            if(this->with_capillary_spreading_) eptr->Store(this->key_gradSn, snw_gradient);
            eptr->Store(this->key_gradP, p_gradient);       
        }
    }
}




//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){
    this->rate_count_++;//recording
    nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );
    
    double64 accumulation(0.), flux_balance(0.), outflow(0.);
        
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD, facetNrml, gravity;
    const size_t node_parent_elements(nd->Parents());
    
    double64 cfl_multiplier = this->CFL_multiplier_*this->relaxing_factor_; //default value
    
    int truncated_node = static_cast<int>(nd->Read(this->key_cut));//check if node is truncated by domain boundary
    
    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
      Element<dim>* const eptr(nd->Parent(t));
      assert( eptr != NULL );
        
      if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) { //ignore if parent element located outside domain
        continue;
        
      } else {        
        
        const size_t pnid(nd->ParentNodeNumber(t));
        
        //eptr->Read( this->this->key_vt, vD);
        
        //compute total velocity (without gravity)
        VectorVariable<dim> gradP;
        eptr->Read(this->key_gradP, gradP); //pressure gradient
        double64 lambda_t = this->flowfunctions_.TotalMobility(eptr);
        double64 thickness = eptr->Read(this->key_thi); //thickness
        if(!this->tensor_k_) { //scalar permeability
            double64 k = eptr->Read( this->key_k ); //permeability
            k *= lambda_t;
            if (!isnan(thickness)) k *= thickness;
            vD(0) = k * gradP(0);
            if ( dim != 1U ) vD(1) = k * gradP(1);
            if ( dim == 3U ) vD(2) = k * gradP(2);
        } else { //tensor permeability
            TensorVariable<dim> kk;
            eptr->Read( key_kk, kk );
            kk *= lambda_t;
            if (!isnan(thickness)) kk *= thickness;
            vD= kk * gradP;
        }
        /*
        if( this->with_gravity_forces_ ){ //take into account gravity effect
            double64 delta_rho = rho_w - rho_n;
            double64 gravity_t = lambda_t * delta_rho * ACC_GRAVITY;
            vD(v) += gravity_t;
        }
        */         
        
        double64 inflow(0.), CO2_inflow (0.);
        const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
        for ( size_t i=0U; i<sector_facets; i++ )
        {
            const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
            const size_t inside_node(eptr->FV()->InsideNode(iFacet));
            const size_t outside_node(eptr->FV()->OutsideNode(iFacet));
            
            eptr->Read( iFacet, 0U,  this->key_fn, facetNrml );
            const double64  vD_n = vD.DotProduct(facetNrml);
            const double64  facetArea = eptr->Read( iFacet, 0U,  this->key_fA ); 
            
            const double64 sign = ( pnid == inside_node ) ? 1. : -1.;
            //compute facet fluid flux
            double64 facet_flux = sign * vD_n * facetArea;
            //update flux balance
            flux_balance += facet_flux;
            inflow += facet_flux;
            //update outflow
            if ( facet_flux > 0. ) outflow += facet_flux;
                         
            //compute inside and outside node mobilities, by using their saturations
            const double64 sn_inside_node = eptr->N(inside_node)->Read(  this->key_sCO2 );
            const double64 sw_inside_node = 1.-sn_inside_node;
            const double64 ln_inside_node = this->flowfunctions_.Mobility_at(eptr,1U,1.0-sn_inside_node);
            const double64 lw_inside_node = this->flowfunctions_.Mobility_at(eptr,0U,1.0-sn_inside_node);
        
            const double64 sn_outside_node = eptr->N(outside_node)->Read(  this->key_sCO2 );
            const double64 sw_outside_node = 1.-sn_outside_node;
            const double64 ln_outside_node = this->flowfunctions_.Mobility_at(eptr,1U,1.0-sn_outside_node);
            const double64 lw_outside_node = this->flowfunctions_.Mobility_at(eptr,0U,1.0-sn_outside_node);
            
            //compute phase velocities at facet integration point                      
            double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
            
            if( this->with_gravity_forces_ ){   
                this->flowfunctions_.GravityMultiplier(eptr, gravity);                    
                //vn_gravity_component_of_velocity = this->flowfunctions_.Mobility(eptr, 0U) * gravity[v] * facetNrml[v];
                //vw_gravity_component_of_velocity = this->flowfunctions_.Mobility(eptr, 1U) * gravity[v] * facetNrml[v];
                double64 gravity_nrml = gravity.DotProduct(facetNrml);
                vn_gravity_component_of_velocity = this->flowfunctions_.Mobility(eptr, 0U) * gravity_nrml;
                vw_gravity_component_of_velocity = this->flowfunctions_.Mobility(eptr, 1U) * gravity_nrml;                
            }         
            
            if(this->with_capillary_spreading_){  
                VectorVariable<dim> grad;
                eptr->Read(this->key_gradSn, grad);
                double64 dsdn = grad.DotProduct(facetNrml);
            
                if(!isnan(dsdn)){
                    vn_capillary_component_of_velocity = -dsdn*this->flowfunctions_.CapillaryDiffusionMultiplier_Phase(eptr,0U);
                    vw_capillary_component_of_velocity = -dsdn*this->flowfunctions_.CapillaryDiffusionMultiplier_Phase(eptr,1U);
                }   
            }  
            
            double64 vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            double64 vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity; 
        
            //determine upstream mobilities
            double64 upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);    
              
            if(vn_at_facet_int_point>0.0)
                upstream_mobility_n=ln_inside_node;
            else if (vn_at_facet_int_point<0.0)
                upstream_mobility_n=ln_outside_node;
            else 
                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
            
            if(vw_at_facet_int_point>0.0)
                upstream_mobility_w=lw_inside_node;
            else if (vw_at_facet_int_point<0.0)
                upstream_mobility_w=lw_outside_node;
            else
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);
                            
            total_mobility=upstream_mobility_n+upstream_mobility_w;
            double64 upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            double64 upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);        

            //compute each velocity component     
            double64 viscous_velocity_component(0.0), capillary_velocity_component(0.0), gravity_velocity_component(0.0);  
                     
            viscous_velocity_component = vD_n * upstream_fn;
                                       
            if( this->with_gravity_forces_ ) {
                this->flowfunctions_.GravityMultiplier(eptr, gravity);
                gravity_velocity_component = upstream_lambda_overbar * gravity[v] * facetNrml[v];
            }

            if( this->with_capillary_spreading_ )
                capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity; 
        
            //update non-wetting flux accumulation   
            double64 fn = sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;                                        
            accumulation += fn; 
            CO2_inflow += fn;
            

            //determine cfl_multiplier based on non-wetting phase shock saturation
            if (cfl_multiplier != this->CFL_multiplier_){
                if(vn_at_facet_int_point < 0.0) { //flowing in from outside node (upstream node)
                    double64 sn_shock = 1.0-eptr->Read(this->key_ssH2O); //sn at shock for outside node
                    if (sn_outside_node >= sn_shock) { //upstream node passed shock saturation
                        if (sn_inside_node < sn_shock) {//current node not yet reach shock saturation   
                            cfl_multiplier = this->CFL_multiplier_;
                        }
                    }
                }
            }
        
            //determine cfl_multiplier based on wetting phase shock saturation
            if (cfl_multiplier != this->CFL_multiplier_){
                if(vw_at_facet_int_point < 0.0) { //flowing from outside node (upstream node)
                    double64 sw_shock = eptr->Read(this->key_ssH2O); //sw at shock for outside node
                    if (sw_outside_node >= sw_shock) { //upstream node passed shock saturation
                        if (sw_inside_node < sw_shock) {//current node not yet reach shock saturation   
                            cfl_multiplier = this->CFL_multiplier_;
                        }
                    }
                }        
            }
            
        } //end sector_facets loop
        
        if (!this->no_flow_boundary_) {
            //inflow/outflow compensation for truncated boundary node 
            if (truncated_node == 1) {      
                if (inflow > 0.) {flux_balance += inflow; outflow += inflow;} //inflow compensation
                else if (inflow < 0.) {flux_balance -= inflow;}; //outflow compensation
            
                if ( CO2_inflow > 0. ) accumulation += CO2_inflow; //inflow compensation
                else if ( CO2_inflow < 0. ) accumulation -= CO2_inflow; //outflow compensation
            }
        }
        
      } //end else        
    } //end parent element loop
        
    ArrayVariable array2;
    nd->Read(this->key_time, array2);
   
    //compute CFL time increment 
    if (outflow < numeric_limits<double64>::epsilon())
        array2.Component(2, numeric_limits<double64>::max());
    else 
        array2.Component(2, nd->Read(  this->key_fvPV ) / outflow);
    
    nd->Store( this->key_CFL, makeScalar( nd->Status( this->key_CFL ),cfl_multiplier ) );
    array2.Component(6, cfl_multiplier);    
        
    nd->Store(this->key_time, array2);
    
    
    // divergence free correction (only when node is located inside domain and flux balance not equal to zero)
    if (truncated_node !=1 && fabs(flux_balance) > numeric_limits<double64>::epsilon()) {
         // compute average fractional flow for the current finite volume
        double64 fn_avg = 0.;
        double64 sw = 1. - nd->Read(this->key_sCO2); //saturation aqueous phase at current node
        for ( size_t t=0U; t<node_parent_elements; t++ )
        {
            Element<dim>* const eptr(nd->Parent(t));
            fn_avg += this->flowfunctions_.f_at(eptr,1U,sw);
        }
        fn_avg /= static_cast<double64>(node_parent_elements);
        accumulation -= fn_avg*flux_balance;
    }     
    
    //compute and store variation rate    
    double64 PV = nd->Read(  this->key_fvPV );
    nd->Store( key_dsnw, makeScalar( nd->Status( key_dsnw ), accumulation/PV ) );    
  }
    
}  



//schedule an event associated with a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
bool DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Schedule(Event<dim>* event, double64 t_end)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    ArrayVariable array;
    nd->Read(this->key_time, array);

    nd->Store( this->key_schedule, makeScalar( nd->Status( this->key_schedule), nd->Read( this->key_schedule) + 1 ) );
    event->valid(true);
    //compute target change
    double64 dt_CFL = array[2];//CFL time increment
    double64 ChangeRate = nd->Read( key_dsnw);//rate of change
    double64 source = nd->Read(key_NQV);
    double64 dC_CFL = dt_CFL*this->CFL_multiplier_*(-ChangeRate+source);//targe change

    if (fabs(dC_CFL) < numeric_limits<double64>::epsilon()){//idle node/FV
        array.Component(5, numeric_limits<double64>::epsilon());//target change of solution
        array.Component(3, numeric_limits<double64>::max());//target time increment          
    } else {
        array.Component(5, dC_CFL);//target change of solution
        array.Component(3, dt_CFL*array[6]);//target time increment  
    };

    double64 t_current = array[0];//current time stamp
    double64 dt_target = array[3];//target time increment
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



//update solution and check it against the specified range (with DES)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double64 t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    ArrayVariable array;
    nd->Read(this->key_time, array);

    double64 ChangeRate = nd->Read( key_dsnw);//variaition rate   
    double64 solution = nd->Read( this->key_sCO2);//old solution
    nd->Store(key_sCO2_0, makeScalar( status, solution ));//store old solution
    double64 t_current = array[0]; //current time stamp
    double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
    const double64 source(nd->Read(key_NQV));
    new_solution += source * (t_clock - t_current);//add source to new solution
        
    //check new solution value against range and stored it to key_sCO2
    if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    } else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_ << endl;
        if ( new_solution > this->upper_limit_ ) new_solution = this->upper_limit_;
        else if ( new_solution < this->lower_limit_ ) new_solution = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    }  
    new_solution = nd->Read(this->key_sCO2);//stored new solution
    double64 dsn_cumulative = array[4];
    array.Component(4, dsn_cumulative + (new_solution-solution));//update cumulative change
        
    array.Component(0, t_clock); //current time stamp
    
    nd->Store(this->key_time, array);

    nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
  }
}



//update solution and check it against the specified range (with TDS)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double64 delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );   
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording    
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    
    double64 ChangeRate = nd->Read(key_dsnw);//variation rate  
    double64 solution = nd->Read(this->key_sCO2);//old solution
    nd->Store(key_sCO2_0, makeScalar( status, solution ));//store old solution
    double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
    const double64 source(nd->Read( key_NQV));
    new_solution += source * delta_t;//add source to new solution.
                
    //check new solution value against range and stored it to key_sCO2
    if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    } else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_ << endl;
        if ( new_solution > this->upper_limit_ ) new_solution = this->upper_limit_;
        else if ( new_solution < this->lower_limit_ ) new_solution = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    }
        
    nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );   
  }
}




//Synchronize neighbor nodes/FVs
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::Synchronize(Event<dim>* event, double64 t_clock, double64& t_remove )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL ); 
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    nd->Store( this->key_synchronize, makeScalar( nd->Status(this->key_synchronize), nd->Read(this->key_synchronize) + 1 ) );
    event->valid(false);
    ArrayVariable array;
    nd->Read(this->key_time, array);
    array.Component(4, 0.); //reset cumulative change of solution
    nd->Store(this->key_time, array);
    for ( size_t n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* neighbor_node = nd->Neighbor(n);
        if( neighbor_node != NULL && neighbor_node->Status(  this->key_sCO2 ) != DIRICH){
            size_t index = neighbor_node->Read(this->key_EventIndex);
            if(index >= 0 && index < this->FullList.size()){
                Event<dim>* neighbor_event = this->FullList[index];  
                assert( neighbor_event  != NULL ); 
                if (neighbor_event != NULL && neighbor_event->inPEPStack() == false) {
                    this->PEPList.push_back(neighbor_event);
                    neighbor_event->inPEPStack(true);
                    Update_DES(neighbor_event,t_clock);
                    ArrayVariable neighbor_array;
                    neighbor_node->Read(this->key_time, neighbor_array); 
                    double64 dC_cumulative = neighbor_array[4];//cumulative change of solution
                    double64 dC_target = neighbor_array[5];//target change of solution
                    if (fabs(dC_cumulative) >= fabs(dC_target)) {
                        #if defined(_OPENMP)
                        double64 t_begin = omp_get_wtime();
                        #else
                        clock_t t_begin = clock();
                        #endif
                        if (neighbor_event->inQueue()){
                            Heap_Node* neighbor_heap_node = this->HeapNodeFullList[index];
                            this->EventHeap.remove(neighbor_heap_node);
                            neighbor_event->inQueue(false);
                        }
                        #if defined(_OPENMP)
                        t_remove += omp_get_wtime() - t_begin;
                        #else
                        t_remove += clock() - t_begin; 
                        #endif
                        Synchronize (neighbor_event, t_clock,t_remove); 
                    };
                };
            };
        };
    };
  }
}




//advect variable with TDS (time driven simulation)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS( double64 time_interval, size_t num_threads )
{
#if defined(_OPENMP)
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
        cerr <<"WARNING: OpenMP is not available, using serial mode"<<endl;
    AdvectVariable_TDS_serial( time_interval );
#endif
}



//advect variable with DES (discrete event simulation)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES( double64 model_time, size_t num_threads )
{
#if defined(_OPENMP)
    if (num_threads <= 0) {
        cerr <<"WARNING: input number of threads is less than 1, reset to 1"<<endl;
        num_threads = 1;
        AdvectVariable_DES_parallel( model_time, num_threads);
    } 
    else if (num_threads > omp_get_max_threads()) {
        cerr <<"WARNING: input number of threads is larger than maximum available threads (" << omp_get_max_threads() << "), reset to "<< omp_get_max_threads() <<endl;
        num_threads = omp_get_max_threads();
        AdvectVariable_DES_parallel( model_time, num_threads); 
    }  
    else if (num_threads == 1) {
        AdvectVariable_DES_serial( model_time );
    } 
    else {
    AdvectVariable_DES_parallel( model_time, num_threads); 
    }            
#else
    if (num_threads > 1)
        cerr <<"WARNING: OpenMP is not available, using serial mode"<<endl;
    AdvectVariable_DES_serial( model_time );
#endif
}





//advect variable with TDS (time-driven simulation), serial mode
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_serial( double64 time_interval)
{
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_serial "<<endl;
    
    double64 begin=clock();
    double64 time_increment(time_interval); 
    
    clock_t T_begin= clock();
    const typename vector<Event<dim>*>::iterator stack_end(this->PEPList.end());
    for ( typename vector<Event<dim>*>::iterator it=this->PEPList.begin(); it!=stack_end; ++it )     
    {   
        ComputeGradients ((*it));
        ComputeRateofChange((*it));  
        ArrayVariable array;
        (*it)->getNode()->Read(this->key_time, array);
        double64 dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
    }
    this->T_RateOfChange_ += clock() - T_begin; 

    const double64 one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double64 time(0.);
    
    while (time < time_interval)
    {
        cout <<"\n\tadvection (sub)step: "<< substep <<" of total steps "<<max(floor(time_interval/time_increment),one)<< endl;
        T_begin= clock();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        for ( typename vector<Event<dim>*>::iterator it=this->PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };
        this->T_Update_ += clock() - T_begin;
        
        T_begin= clock();
        for ( typename vector<Event<dim>*>::iterator it=this->PEPList.begin(); it!=stack_end; ++it )
        { 
            ComputeGradients ((*it));
            ComputeRateofChange((*it));
            ArrayVariable array2;
            (*it)->getNode()->Read(this->key_time, array2);
            double64 dt_CFL = array2[2];//CFL time increment
            time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
        };
        this->T_RateOfChange_ += clock() - T_begin; 

        time += time_increment;
        substep++;
    };
    
    this->T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS_serial "<<endl;
    cout <<"rate_count_ = "<<this->rate_count_<<endl;
    cout <<"update_count_ = "<<this->update_count_<<endl; 
    cout <<"T_Schedule_ = "<< this->T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< this->T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< this->T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RateOfChange_ = "<< this->T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_InsertToHeap_ = "<< this->T_InsertToHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RemoveFromHeap_ = "<< this->T_RemoveFromHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_AdvectVariable_ = "<< this->T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl;
}




#if defined(_OPENMP)
//advect variable with TDS (time-driven simulation), parallel mode
template<size_t dim>
void DES2PhaseSlightlyCompressibleTransport<dim>::AdvectVariable_TDS_parallel( double64 time_interval, size_t num_threads )
{
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim>::AdvectVariable_TDS_parallel "<<endl;
    cout<<"Using threads = "<<num_threads<<" Maximum available threads ="<< omp_get_max_threads() << endl;

    double64 begin=omp_get_wtime();
    double64 time_increment(time_interval); 
    
    double64 T_begin;
    T_begin = omp_get_wtime();
    
    size_t PEPList_size = this->PEPList.size();
    
    for(size_t i = 0U; i < PEPList_size; ++i)
    {
        auto it = this->PEPList.begin()+i;
        this->ComputeGradients ((*it));
    } 
        
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(dynamic)
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            ComputeRateofChange((*it));  
        };
    }

    for(size_t i = 0U; i < PEPList_size; ++i)
    {
        auto it = this->PEPList.begin()+i;
        Event<dim>* event = *it;     
        ArrayVariable array;
        (*it)->getNode()->Read(this->key_time, array);
        double64 dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
    }
    
    this->T_RateOfChange_ += omp_get_wtime() - T_begin;

    const double64 one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double64 time(0.);
    
    while (time < time_interval)
    {
        cout <<"\n\tadvection (sub)step: "<< substep <<" of total steps "<<max(floor(time_interval/time_increment),one)<< endl;
        T_begin = omp_get_wtime();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        for ( typename vector<Event<dim>*>::iterator it=this->PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };
        this->T_Update_ += omp_get_wtime() - T_begin;
        
        T_begin = omp_get_wtime();
        time_increment = time_interval; 
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            this->ComputeGradients ((*it));
        }  
                
        #pragma omp parallel num_threads(num_threads)
        {        
            #pragma omp for schedule(dynamic)     
            for(size_t i = 0U; i < PEPList_size; ++i)
            {
                auto it = this->PEPList.begin()+i;
                ComputeRateofChange((*it));;
            }
        }        
        
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            ArrayVariable array2;
            (*it)->getNode()->Read(this->key_time, array2);
            double64 dt_CFL = array2[2];//CFL time increment
            time_increment=min(time_increment, dt_CFL*this->CFL_multiplier_);
        };

        this->T_RateOfChange_ += omp_get_wtime() - T_begin; 

        time += time_increment;
        substep++;
    };
    
    this->T_AdvectVariable_+= omp_get_wtime() - begin;

    cout <<"Finish DES2PhaseSlightlyCompressibleTransport<dim>::AdvectVariable_TDS_parallel "<<endl;
    cout <<"rate_count_ = "<<this->rate_count_<<endl;
    cout <<"update_count_ = "<<this->update_count_<<endl; 
    cout <<"T_Schedule_ = "<< this->T_Schedule_ << endl;
    cout <<"T_Update_  = "<< this->T_Update_ << endl;
    cout <<"T_Synchronize_ = "<< this->T_Synchronize_ << endl;
    cout <<"T_RateOfChange_ = "<< this->T_RateOfChange_ << endl;
    cout <<"T_InsertToHeap_ = "<< this->T_InsertToHeap_ << endl; 
    cout <<"T_RemoveFromHeap_ = "<< this->T_RemoveFromHeap_ << endl; 
    cout <<"T_AdvectVariable_ = "<< this->T_AdvectVariable_ << endl;
}
#endif




//advect variable with DES (discrete event simulation), serial mode
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial( double64 model_time)
{
    double64 begin=clock();
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial "<<endl;
    double64 time(0.);
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
    
    while (!Finished)
    {   
        clock_t T_begin;
        const typename vector<Event<dim>*>::iterator stack_end(this->PEPList.end());
        for ( typename vector<Event<dim>*>::iterator it=this->PEPList.begin(); it!=stack_end; ++it )       
        {   
            Event<dim>* event = *it;
            T_begin= clock();
            ComputeGradients (event );
            ComputeRateofChange((*it));    
            this->T_RateOfChange_ += clock() - T_begin;         
            if ((*it)->valid() == false) {
                T_begin= clock();
                bool isactive = Schedule(event, model_time);   
                this->T_Schedule_ += clock() - T_begin;              
                if (isactive) {
                    T_begin= clock();
                    double64 scheduled_time = event->t_schedule();
                    int index = static_cast<int>(event->getNode()->Read(this->key_EventIndex));
                    Heap_Node* heap_node = new Heap_Node(scheduled_time,index);
                    this->EventHeap.insert(heap_node);
                    this->HeapNodeFullList[index] = heap_node;
                    event->inQueue(true);
                    this->T_InsertToHeap_ += clock() - T_begin; 
                }  
            } 
            event->inPEPStack(false);         
        };    
        this->T_RateOfChange_ += clock() - T_begin; 

            
        if (this->EventHeap.empty()) time=model_time;
        else time = this->EventHeap.minimum()->getK();
        
        cout<<"  time = "<<time<<" model_time = "<<model_time<<" PEPList size = "<< this->PEPList.size() <<" Queue size = "<< this->EventHeap.size()<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Event<dim>*>::iterator End(this->PEPList.end());
            for ( typename vector<Event<dim>*>::iterator e=this->PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        this->PEPList.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
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
            double64 dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, this->PEP_multiplier_*dt_target);
            double64 t_schedule = array[1];//scheduled time stamp
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
            double64 t_remove(0.);    
            Synchronize(top_event,time,t_remove);
            this->T_RemoveFromHeap_ += t_remove;
            this->T_Synchronize_ += clock() - T_begin - t_remove;
        };
            
        //cout<<"  iteration count = "<<count<<endl;
    };
    /*
    //reset all events and add them to PEPList (for advection at next integration step)
    this->EventHeap.clear();
    this->PEPList.clear();
    const typename vector<Event<dim>*>::iterator stack_end(this->FullList.end());
    for ( typename vector<Event<dim>*>::iterator it=this->FullList.begin(); it!=stack_end; ++it )       
    {  
        Event<dim>* event = *it;
        PEPList.push_back(event);
        event->inPEPStack(true);          
        event->valid(false);
        event->inQueue(false);
    }    
    */
    this->T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial "<<endl;
    cout <<"rate_count_ = "<<this->rate_count_<<endl;
    cout <<"update_count_ = "<<this->update_count_<<endl; 
    cout <<"T_Schedule_ = "<< this->T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< this->T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< this->T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RateOfChange_ = "<< this->T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_InsertToHeap_ = "<< this->T_InsertToHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RemoveFromHeap_ = "<< this->T_RemoveFromHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_AdvectVariable_ = "<< this->T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl;
}   



#if defined(_OPENMP)
//advect variable with DES (discrete event simulation), parallel mode
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_parallel( double64 model_time, size_t num_threads)
{
    double64 begin=omp_get_wtime();
    cout<<"Start DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_parallel "<<endl;
    cout << "Using threads = "<<num_threads<<" Maximum available threads ="<< omp_get_max_threads() << endl;
    double64 time(0.);
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
   
    while (!Finished)
    {   
        double64 T_begin;
        
        T_begin = omp_get_wtime();
                
        size_t PEPList_size = this->PEPList.size();
        
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = this->PEPList.begin()+i;
            Event<dim>* event = *it;   
            ComputeGradients (event );
        }

        std::vector<Event<dim>*> tempList;
        #pragma omp parallel num_threads(num_threads)
        {
            std::vector<Event<dim>*> privateList;
            
            #pragma omp for schedule(dynamic)
            for(size_t i = 0U; i < PEPList_size; ++i)
            {
                
                auto it = this->PEPList.begin()+i;
                Event<dim>* event = *it;                                
                ComputeRateofChange((*it)); 
                
                if (event->valid() == false) 
                    if (Schedule(event, model_time)) privateList.push_back(event);
                event->inPEPStack(false); 
            };
            
            #pragma omp critical
            tempList.insert(tempList.end(), privateList.begin(), privateList.end());            
        } 
        
        size_t tempList_size = tempList.size();
        for(size_t i = 0U; i < tempList_size; ++i)
        {
            auto it = tempList.begin()+i;
            Event<dim>* event = *it;                 
            double64 scheduled_time = event->t_schedule();
            size_t index = event->getNode()->Read(this->key_EventIndex);
            Heap_Node* heap_node = new Heap_Node(scheduled_time,index);                    
            this->EventHeap.insert(heap_node);
            this->HeapNodeFullList[index] = heap_node;
            event->inQueue(true);
        }
        tempList.clear();
        
        
        this->T_RateOfChange_ += omp_get_wtime() - T_begin; 
            
        if (this->EventHeap.empty()) time=model_time;
        else time = this->EventHeap.minimum()->getK();
        
        cout<<"  time = "<<time<<" model_time = "<<model_time<<" PEPList size = "<< this->PEPList.size() <<" Queue size = "<< this->EventHeap.size()<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Event<dim>*>::iterator End(this->PEPList.end());
            for ( typename vector<Event<dim>*>::iterator e=this->PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        this->PEPList.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
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
            double64 dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, this->PEP_multiplier_*dt_target);
            double64 t_schedule = array[1];//scheduled time stamp
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
            double64 t_remove(0.);    
            Synchronize(top_event,time,t_remove);
            this->T_RemoveFromHeap_ += t_remove;
            this->T_Synchronize_ += omp_get_wtime() - T_begin - t_remove;
        };
            
        //cout<<"  iteration count = "<<count<<endl;
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
    this->T_AdvectVariable_+= omp_get_wtime() - begin;

    cout<<"Finish DES2PhaseSlightlyCompressibleTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_parallel "<<endl;
    cout <<"rate_count_ = "<<this->rate_count_<<endl;
    cout <<"update_count_ = "<<this->update_count_<<endl; 
    cout <<"T_Schedule_ = "<< this->T_Schedule_  << endl;
    cout <<"T_Update_  = "<< this->T_Update_  << endl;
    cout <<"T_Synchronize_ = "<< this->T_Synchronize_<< endl;
    cout <<"T_RateOfChange_(including_T_Schedule_) = "<< this->T_RateOfChange_ << endl;
    cout <<"T_InsertToHeap_ = "<< this->T_InsertToHeap_ << endl; 
    cout <<"T_RemoveFromHeap_ = "<< this->T_RemoveFromHeap_ << endl; 
    cout <<"T_AdvectVariable_ = "<< this->T_AdvectVariable_ << endl;
} 
#endif





template class DES2PhaseSlightlyCompressibleTransport<1U,FlowFunctionsModule1>;
template class DES2PhaseSlightlyCompressibleTransport<2U,FlowFunctionsModule1>;
template class DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule1>;

template class DES2PhaseSlightlyCompressibleTransport<1U,FlowFunctionsModule2>;
template class DES2PhaseSlightlyCompressibleTransport<2U,FlowFunctionsModule2>;
template class DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule2>;

template class DES2PhaseSlightlyCompressibleTransport<1U,FlowFunctionsModule3>;
template class DES2PhaseSlightlyCompressibleTransport<2U,FlowFunctionsModule3>;
template class DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule3>;


} // end csmp 



                                                                       
