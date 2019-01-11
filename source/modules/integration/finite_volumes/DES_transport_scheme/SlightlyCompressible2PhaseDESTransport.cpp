#include "SlightlyCompressible2PhaseDESTransport.h"
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
SlightlyCompressible2PhaseDESTransport<dim,FLOW_FUNCTIONS>::SlightlyCompressible2PhaseDESTransport( Model<dim>& m,
                                                                                                    const char* target_region,
                                                                                                    FLOW_FUNCTIONS<dim>& ff,
                                                                                                    bool with_capillary_spreading,
                                                                                                    bool with_gravity_forces,
                                                                                                    bool tensor_k,
                                                                                                    double64 PEP_multiplier,
                                                                                                    double64 cfl_multiplier )
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,ff,with_capillary_spreading,with_gravity_forces,tensor_k,PEP_multiplier,cfl_multiplier)
{
    cout<<"SlightlyCompressible2PhaseDESTransport constructed"<<endl;
} // end constructor  


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
SlightlyCompressible2PhaseDESTransport<dim,FLOW_FUNCTIONS>::SlightlyCompressible2PhaseDESTransport(  Model<dim>& m,
                                                                                                     const char* target_region,
                                                                                                     FLOW_FUNCTIONS<dim>& ff,
                                                                                                     bool with_capillary_spreading,
                                                                                                     bool with_gravity_forces,
                                                                                                     bool tensor_k,
                                                                                                     double64 PEP_multiplier,
                                                                                                     double64 cfl_multiplier,
                                                                                                     double64 relaxing_factor )
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,ff,with_capillary_spreading,with_gravity_forces,tensor_k,PEP_multiplier,cfl_multiplier, relaxing_factor)
{
    cout<<"SlightlyCompressible2PhaseDESTransport constructed"<<endl;
} // end constructor 




//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void SlightlyCompressible2PhaseDESTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
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
        
        //eptr->Read( this->key_vt, vD);
        
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
            eptr->Read( this->key_kk, kk );
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
    nd->Store( this->key_dsnw, makeScalar( nd->Status( this->key_dsnw ), accumulation/PV ) );    
  }
    
}  



//schedule an event associated with a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
bool SlightlyCompressible2PhaseDESTransport<dim,FLOW_FUNCTIONS>::Schedule(Event<dim>* event, double64 t_end)
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
    double64 ChangeRate = nd->Read( this->key_dsnw);//rate of change
    double64 source = nd->Read(this->key_nQV);
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
void SlightlyCompressible2PhaseDESTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double64 t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    ArrayVariable array;
    nd->Read(this->key_time, array);

    double64 ChangeRate = nd->Read( this->key_dsnw);//variaition rate   
    double64 solution = nd->Read( this->key_sCO2);//old solution
    nd->Store(this->key_sCO2_0, makeScalar( status, solution ));//store old solution
    double64 t_current = array[0]; //current time stamp
    double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
    const double64 source(nd->Read(this->key_nQV));
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
void SlightlyCompressible2PhaseDESTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double64 delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );   
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording    
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    
    double64 ChangeRate = nd->Read(this->key_dsnw);//variation rate  
    double64 solution = nd->Read(this->key_sCO2);//old solution
    nd->Store(this->key_sCO2_0, makeScalar( status, solution ));//store old solution
    double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
    const double64 source(nd->Read( this->key_nQV));
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


template class SlightlyCompressible2PhaseDESTransport<1U,FlowFunctionsModule1>;
template class SlightlyCompressible2PhaseDESTransport<2U,FlowFunctionsModule1>;
template class SlightlyCompressible2PhaseDESTransport<3U,FlowFunctionsModule1>;

template class SlightlyCompressible2PhaseDESTransport<1U,FlowFunctionsModule2>;
template class SlightlyCompressible2PhaseDESTransport<2U,FlowFunctionsModule2>;
template class SlightlyCompressible2PhaseDESTransport<3U,FlowFunctionsModule2>;

template class SlightlyCompressible2PhaseDESTransport<1U,FlowFunctionsModule3>;
template class SlightlyCompressible2PhaseDESTransport<2U,FlowFunctionsModule3>;
template class SlightlyCompressible2PhaseDESTransport<3U,FlowFunctionsModule3>;


} // end csmp 



                                                                       
