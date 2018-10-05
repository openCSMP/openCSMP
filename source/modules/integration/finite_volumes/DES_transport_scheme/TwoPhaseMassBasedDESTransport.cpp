#include "TwoPhaseMassBasedDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "CO2H2O_FunctionsModule1.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseMassBasedDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier)
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,with_capillary_spreading,with_gravity_forces,PEP_multiplier,cfl_multiplier) 
{
    InitializeVariablesAndKeys(m);
    cout<<"TwoPhaseMassBasedDESTransport constructed."<<endl;
} // end constructor  






template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseMassBasedDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier,
                                                 double64 relaxing_factor)
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,with_capillary_spreading,with_gravity_forces,PEP_multiplier,cfl_multiplier, relaxing_factor)
{
    InitializeVariablesAndKeys(m);
    cout<<"TwoPhaseMassBasedDESTransport constructed."<<endl;
} // end constructor 

  



template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys(Model<dim>& m)
{
    //creating new variables if not defined yet from input file
    if(!m.Database().IsDefined("mass variation rate carbonic phase")) m.CreateProperty( "mass variation rate carbonic phase", "kg/s", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);
    
    key_dmCO2 = INDEX<SCALAR,NODE> ( m.Database().StorageKey("mass variation rate carbonic phase") );
    
    //checking keys 
    if ( key_dmCO2.place != NODE || key_dmCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys:",
        "The 'mass variation rate carbonic phase' variable must be SCALAR and placed on NODE"  );       

    // model-wide initialisation
    m.Region("Model").InputPropertyValue( "mass variation rate carbonic phase", makeScalar(PLAIN,0.), COMPLETE);  
}



//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){ 
    this->rate_count_++;//recording
    nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );
    
    FLOW_FUNCTIONS<dim> flowfunctions(this->db_);

    double64 flux_balance(0.), outflow(0.);
    double64 mCO2_accumulation(0.);
        
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD, facetNrml;
    const size_t node_parent_elements(nd->Parents());
    
    double64 cfl_multiplier = this->CFL_multiplier_*this->relaxing_factor_; //default value
    
    size_t truncated_node = nd->Read(this->key_cut);//check if node is truncated by domain boundary
           
    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
      Element<dim>* const eptr(nd->Parent(t));
      assert( eptr != NULL );
        
      if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) { //ignore if parent element located outside domain
        continue;
        
      } else {
        flowfunctions.InitialiseBrooksCoreyParameters(eptr); 
        
        const size_t pnid(nd->ParentNodeNumber(t));
        
        double64 rho_n = eptr->PropertyValueAtBaryCenter(this->key_rhoCO2);
        double64 rho_w = eptr->PropertyValueAtBaryCenter(this->key_rhoH2O);    

        //compute total velocity (without gravity)
        VectorVariable<dim> gradP;
        eptr->Read(this->key_gradP, gradP); //pressure gradient
        double64 k = eptr->Read( this->key_k ); //permeability
        double64 lambda_t = flowfunctions.TotalMobility(eptr);
        lambda_t *= k;
        double64 thickness = eptr->Read(this->key_thi); //thickness
        if (!isnan(thickness)) lambda_t *= thickness;
        vD(0) = lambda_t * gradP(0);
        if ( dim != 1U ) vD(1) = lambda_t * gradP(1);
        if ( dim == 3U ) vD(2) = lambda_t * gradP(2);   

        /*
        if( this->with_gravity_forces_ ){ //take into account gravity effect
            double64 delta_rho = rho_w - rho_n;
            double64 gravity_t = lambda_t * delta_rho * ACC_GRAVITY;
            vD(v) += gravity_t;
        }
        */     
              
        double64 inflow(0.), mCO2_inflow (0.);
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
            const double64 ln_inside_node = flowfunctions.Mobility_at(eptr,1U,1.0-sn_inside_node); 
            const double64 lw_inside_node = flowfunctions.Mobility_at(eptr,0U,1.0-sn_inside_node); 
        
            const double64 sn_outside_node = eptr->N(outside_node)->Read(  this->key_sCO2 );
            const double64 sw_outside_node = 1.-sn_outside_node;
            const double64 ln_outside_node = flowfunctions.Mobility_at(eptr,1U,1.0-sn_outside_node); 
            const double64 lw_outside_node = flowfunctions.Mobility_at(eptr,0U,1.0-sn_outside_node);   
            
            //compute phase velocities at facet integration point                      
            double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
            if( this->with_gravity_forces_ ){                       
                vn_gravity_component_of_velocity = flowfunctions.Mobility(eptr, 0U) * flowfunctions.GravityTerm(eptr) * facetNrml[v];
                vw_gravity_component_of_velocity = flowfunctions.Mobility(eptr, 1U) * flowfunctions.GravityTerm(eptr) * facetNrml[v];
            }         
            
            if(this->with_capillary_spreading_){  
                VectorVariable<dim> grad;
                eptr->Read(this->key_gradSn, grad);
                double64 dsdn = grad.DotProduct(facetNrml);
            
                if(!isnan(dsdn)){
                    vn_capillary_component_of_velocity = -dsdn*flowfunctions.CapillaryDiffusionMultiplier_Phase(eptr,0U);
                    vw_capillary_component_of_velocity = -dsdn*flowfunctions.CapillaryDiffusionMultiplier_Phase(eptr,1U);
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
            
            //double64 upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);
            double64 total_mobility_rho = upstream_mobility_n * rho_n + upstream_mobility_w * rho_w;
            double64 upstream_lambda_overbar_rho=(total_mobility_rho!=0.0? (upstream_mobility_n*rho_n*upstream_mobility_w*rho_w)/total_mobility_rho : 0.0);
            
            //compute each velocity component     
            double64 viscous_velocity_component(0.0), capillary_velocity_component(0.0), gravity_velocity_component(0.0); 
                     
            viscous_velocity_component = vD_n * upstream_fn * rho_n;
                                       
            if( this->with_gravity_forces_ ) {
                //gravity_velocity_component = rho_n * upstream_lambda_overbar * flowfunctions.GravityTerm(eptr) * facetNrml[v];
                gravity_velocity_component = upstream_lambda_overbar_rho * flowfunctions.GravityTerm(eptr) * facetNrml[v]; 
            }
            
            if( this->with_capillary_spreading_ ) {
                double64 mobility_product = flowfunctions.MobilityProduct(eptr);
                capillary_velocity_component=(mobility_product!=0.0 ? upstream_lambda_overbar_rho*flowfunctions.CapillaryDiffusionMultiplier(eptr)/mobility_product : 0.0);
            }
                
            //update non-wetting flux accumulation 
            double64 fn = sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea; // m3/s
            mCO2_accumulation += fn; // kg/s
            mCO2_inflow += fn; // kg/s           
            
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
            
                if ( mCO2_inflow > 0. ) mCO2_accumulation += mCO2_inflow; //inflow compensation
                else if ( mCO2_inflow < 0. ) mCO2_accumulation -= mCO2_inflow; //outflow compensation
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
            flowfunctions.InitialiseBrooksCoreyParameters(eptr);   
            fn_avg += flowfunctions.f_at(eptr,1U,sw);
        }
        fn_avg /= static_cast<double64>(node_parent_elements);
        double64 rhon = nd->Read(this->key_rhoCO2);
        mCO2_accumulation -= fn_avg*flux_balance*rhon; // kg/s
    } 
    
    //compute and store variation rate    
    nd->Store( this->key_dmCO2, makeScalar( nd->Status( this->key_dmCO2 ), mCO2_accumulation ) ); // kg/s 
  }
}  



//schedule an event associated with a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
bool TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::Schedule(Event<dim>* event, double64 t_end)
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
    
    double64 mCO2_rate = nd->Read( key_dmCO2); //variaition rate of mass carbonic phase (kg/s)
    double64 rhoCO2 = nd->Read(this->key_rhoCO2); //density carbonic phase (kg/m3)
    double64 PV = nd->Read(this->key_fvPV); //pore volume (m3)
    double64 sCO2_rate = mCO2_rate/rhoCO2/PV; //variaition rate of saturation carbonic phase (m3/(m3.s)) 
    
    double64 source = nd->Read(this->key_nQV);
    double64 dC_CFL = dt_CFL*this->CFL_multiplier_*(-sCO2_rate+source);//targe change

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



//update solution and check it against the specified range (with DES). Also update mass fraction of each
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double64 t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording   
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    ArrayVariable array;
    nd->Read(this->key_time, array);
    
    //update mass for carbonic phase
    double64 PV = nd->Read(this->key_fvPV); //pore volume (m3)
    double64 rhoCO2 = nd->Read(this->key_rhoCO2); //density carbonic phase (kg/m3)
    double64 sCO2_old = nd->Read( this->key_sCO2);//old saturation carbonic phase (m3/m3)
    nd->Store(this->key_sCO2_0, makeScalar( status, sCO2_old )); //store old saturation  
      
    double64 mCO2_old = PV * sCO2_old * rhoCO2; //old mass carbonic phase (kg)
    double64 mCO2_rate = nd->Read( key_dmCO2); //variaition rate of mass carbonic phase (kg/s)
    
    double64 t_current = array[0]; //current time stamp (s)
    double64 mCO2_new = mCO2_old - mCO2_rate * (t_clock - t_current); //new mass carbonic phase (kg)
        
    double64 source = nd->Read(this->key_nQV); // nodal fluid volume source (m3/(m3.s))
    mCO2_new += PV * source * rhoCO2 * (t_clock - t_current); //add source (kg)
    
    //compute new saturations from updated mass
    double64 sCO2_new = mCO2_new/rhoCO2/PV;
    
    //check new saturation value against range and stored it to this->key_sCO2
    if ( sCO2_new <= this->upper_limit_ && sCO2_new >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, sCO2_new ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - sCO2_new ));
    } else {
        cerr <<"saturation carbonic phase value: "<< sCO2_new <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_ << endl;
        if ( sCO2_new > this->upper_limit_ ) sCO2_new = this->upper_limit_;
        else if ( sCO2_new < this->lower_limit_ ) sCO2_new = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, sCO2_new ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - sCO2_new ));
    }
      
    sCO2_new = nd->Read(this->key_sCO2);//stored new saturation
    double64 dsn_cumulative = array[4];
    array.Component(4, dsn_cumulative + (sCO2_new-sCO2_old));//update cumulative change of saturation        
    
    array.Component(0, t_clock); //current time stamp
    nd->Store(this->key_time, array);

    nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
  }
}


//update solution and check it against the specified range (with TDS)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseMassBasedDESTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double64 delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );   
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording     
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    
    //update mass for carbonic phase
    double64 PV = nd->Read(this->key_fvPV); //pore volume (m3)
    double64 rhoCO2 = nd->Read(this->key_rhoCO2); //density carbonic phase (kg/m3)
    double64 sCO2_old = nd->Read( this->key_sCO2);//old saturation carbonic phase (m3/m3)
    nd->Store(this->key_sCO2_0, makeScalar( status, sCO2_old )); //store old saturation  
      
    double64 mCO2_old = PV * sCO2_old * rhoCO2; //old mass carbonic phase (kg)
    double64 mCO2_rate = nd->Read( key_dmCO2); //variaition rate of mass carbonic phase (kg/s)
    double64 mCO2_new = mCO2_old - mCO2_rate * delta_t; //new mass carbonic phase (kg)
        
    double64 source = nd->Read(this->key_nQV); // nodal fluid volume source (m3/(m3.s))
    mCO2_new += PV * source * rhoCO2 * delta_t; //add source (kg)
    
    //compute new saturations from updated mass
    double64 sCO2_new = mCO2_new/rhoCO2/PV;
    //check new saturation value against range and stored it to this->key_sCO2
    if ( sCO2_new <= this->upper_limit_ && sCO2_new >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, sCO2_new ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - sCO2_new ));
    } else {
        cerr <<"saturation carbonic phase value: "<< sCO2_new <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_ << endl;
        if ( sCO2_new > this->upper_limit_ ) sCO2_new = this->upper_limit_;
        else if ( sCO2_new < this->lower_limit_ ) sCO2_new = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, sCO2_new ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - sCO2_new ));
    }  

    nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
  }
}

template class TwoPhaseMassBasedDESTransport<1U,CO2H2O_FunctionsModule0>;
template class TwoPhaseMassBasedDESTransport<2U,CO2H2O_FunctionsModule0>;
template class TwoPhaseMassBasedDESTransport<3U,CO2H2O_FunctionsModule0>;

template class TwoPhaseMassBasedDESTransport<1U,CO2H2O_FunctionsModule1>;
template class TwoPhaseMassBasedDESTransport<2U,CO2H2O_FunctionsModule1>;
template class TwoPhaseMassBasedDESTransport<3U,CO2H2O_FunctionsModule1>;

template class TwoPhaseMassBasedDESTransport<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseMassBasedDESTransport<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseMassBasedDESTransport<3U,CO2H2O_FunctionsModule2>;

} // end csmp 



                                                                       
