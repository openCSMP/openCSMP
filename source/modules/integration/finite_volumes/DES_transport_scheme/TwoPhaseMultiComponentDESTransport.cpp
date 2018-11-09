#include "TwoPhaseMultiComponentDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "CO2H2O_FunctionsModule1.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseMultiComponentDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 bool tensor_k,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier)
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,with_capillary_spreading,with_gravity_forces,tensor_k,PEP_multiplier,cfl_multiplier) 
{
    this->SetEquilibration(true);
    InitializeVariablesAndKeys(m);
    cout<<"TwoPhaseMultiComponentDESTransport constructed."<<endl;
} // end constructor  






template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseMultiComponentDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 bool tensor_k,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier,
                                                 double64 relaxing_factor)
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,with_capillary_spreading,with_gravity_forces,tensor_k,PEP_multiplier,cfl_multiplier, relaxing_factor)
{
    this->SetEquilibration(true);
    InitializeVariablesAndKeys(m);
    cout<<"TwoPhaseMultiComponentDESTransport constructed."<<endl;
} // end constructor 

  



template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys(Model<dim>& m)
{
    //creating new variables if not defined yet from input file
    if(!m.Database().IsDefined("mass variation rates carbonic phase")) m.CreateProperty( "mass variation rates carbonic phase", "kg/s", ARRAY, NODE, 2, -1.00E+08 ,1.00E+08);
    if(!m.Database().IsDefined("mass variation rates aqueous phase")) m.CreateProperty( "mass variation rates aqueous phase", "kg/s", ARRAY, NODE, 3, -1.00E+08 ,1.00E+08);
    
    key_dcmCO2 = INDEX<ARRAY,NODE> ( m.Database().StorageKey("mass variation rates carbonic phase") );
    key_dcmH2O = INDEX<ARRAY,NODE> ( m.Database().StorageKey("mass variation rates aqueous phase") );
    
    //checking keys 
    if ( key_dcmCO2.place != NODE || key_dcmCO2.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys:",
        "The 'mass variation rates carbonic phase' variable must be ARRAY and placed on NODE"  );       
    if ( key_dcmH2O.place != NODE || key_dcmH2O.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys:",
        "The 'mass variation rates aqueous phase' variable must be ARRAY and placed on NODE"  );         
        
    // model-wide initialisation
    m.Region("Model").InputPropertyValue( "mass variation rates carbonic phase", ArrayVariable(2,0.,PLAIN), COMPLETE);
    m.Region("Model").InputPropertyValue( "mass variation rates aqueous phase", ArrayVariable(3,0.,PLAIN), COMPLETE);  
}



//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){ 
    this->rate_count_++;//recording
    nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );
    
    double64 flux_balance(0.), outflow(0.);
    double64 YCO2_accumulation(0.), YH2O_accumulation(0.), XH2O_accumulation(0.), XCO2_accumulation(0.), XNaCl_accumulation(0.);
    double64 dsdn;
        
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
        this->flowfunctions_->InitialiseBrooksCoreyParameters(eptr); 
        
        const size_t pnid(nd->ParentNodeNumber(t));
        
        double64 rho_n = eptr->PropertyValueAtBaryCenter(this->key_rhoCO2);
        double64 rho_w = eptr->PropertyValueAtBaryCenter(this->key_rhoH2O);    

        //compute total velocity (without gravity)
        VectorVariable<dim> gradP;
        eptr->Read(this->key_gradP, gradP); //pressure gradient
        double64 lambda_t = this->flowfunctions_->TotalMobility(eptr);
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
              
        double64 inflow(0.), YCO2_inflow (0.), YH2O_inflow (0.), XH2O_inflow (0.), XCO2_inflow (0.), XNaCl_inflow (0.);
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
            const double64 ln_inside_node = this->flowfunctions_->Mobility_at(eptr,1U,1.0-sn_inside_node); 
            const double64 lw_inside_node = this->flowfunctions_->Mobility_at(eptr,0U,1.0-sn_inside_node); 
            
            const double64 sn_outside_node = eptr->N(outside_node)->Read(  this->key_sCO2 );
            const double64 sw_outside_node = 1.-sn_outside_node;
            const double64 ln_outside_node = this->flowfunctions_->Mobility_at(eptr,1U,1.0-sn_outside_node); 
            const double64 lw_outside_node = this->flowfunctions_->Mobility_at(eptr,0U,1.0-sn_outside_node);   
            
            //compute phase velocities at facet integration point                      
            double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
            if( this->with_gravity_forces_ ){                       
                vn_gravity_component_of_velocity = this->flowfunctions_->Mobility(eptr, 0U) * this->flowfunctions_->GravityTerm(eptr) * facetNrml[v];
                vw_gravity_component_of_velocity = this->flowfunctions_->Mobility(eptr, 1U) * this->flowfunctions_->GravityTerm(eptr) * facetNrml[v];
            }         
            
            if(this->with_capillary_spreading_){  
                VectorVariable<dim> grad;
                eptr->Read(this->key_gradSn, grad);
                dsdn = grad.DotProduct(facetNrml);
            
                if(!isnan(dsdn)){
                    vn_capillary_component_of_velocity = -dsdn*this->flowfunctions_->CapillaryDiffusionMultiplier_Phase(eptr,0U);
                    vw_capillary_component_of_velocity = -dsdn*this->flowfunctions_->CapillaryDiffusionMultiplier_Phase(eptr,1U);
                }   
            }  
            
            double64 vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            double64 vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity; 
        
            //determine upstream mobilities
            double64 upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);   
            
            ArrayVariable upstream_CO2_comp, upstream_H2O_comp; 
            
            if(vn_at_facet_int_point>0.0) { 
                upstream_mobility_n=ln_inside_node;
                eptr->N(inside_node)->Read(this->key_CO2_comp, upstream_CO2_comp); 
            } else if (vn_at_facet_int_point<0.0) {
                upstream_mobility_n=ln_outside_node;
                eptr->N(outside_node)->Read(this->key_CO2_comp, upstream_CO2_comp);
            } else {
                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                ArrayVariable inside_CO2_comp, outside_CO2_comp;
                eptr->N(inside_node)->Read(this->key_CO2_comp, inside_CO2_comp);
                eptr->N(outside_node)->Read(this->key_CO2_comp, outside_CO2_comp); 
                upstream_CO2_comp = (inside_CO2_comp+outside_CO2_comp)/2.; 
            }
            
            if(vw_at_facet_int_point>0.0) {
                upstream_mobility_w=lw_inside_node;
                eptr->N(inside_node)->Read(this->key_H2O_comp, upstream_H2O_comp); 
            } else if (vw_at_facet_int_point<0.0) {
                upstream_mobility_w=lw_outside_node;
                eptr->N(outside_node)->Read(this->key_H2O_comp, upstream_H2O_comp);
            } else {
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);
                ArrayVariable inside_H2O_comp, outside_H2O_comp;
                eptr->N(inside_node)->Read(this->key_H2O_comp, inside_H2O_comp);
                eptr->N(outside_node)->Read(this->key_H2O_comp, outside_H2O_comp);
                upstream_H2O_comp = (inside_H2O_comp+outside_H2O_comp)/2.;               
            }
            
            total_mobility=upstream_mobility_n+upstream_mobility_w;
            double64 upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            double64 upstream_fw=(total_mobility!=0.0? upstream_mobility_w/total_mobility : 0.0);
            
            double64 total_mobility_rho = upstream_mobility_n * rho_n + upstream_mobility_w * rho_w;
            double64 upstream_lambda_overbar_rho=(total_mobility_rho!=0.0? (upstream_mobility_n*rho_n*upstream_mobility_w*rho_w)/total_mobility_rho : 0.0);
            
            //compute each velocity component     
            double64 viscous_velocity_component_n(0.0), capillary_velocity_component_n(0.0), gravity_velocity_component(0.0); 
            double64 viscous_velocity_component_w(0.0), capillary_velocity_component_w(0.0); 
  
            viscous_velocity_component_n = vD_n * upstream_fn * rho_n;
            viscous_velocity_component_w = vD_n * upstream_fw * rho_w;

            if( this->with_gravity_forces_ ) {
                //gravity_velocity_component = rho_n * upstream_lambda_overbar * this->flowfunctions_->GravityTerm(eptr) * facetNrml[v];
                gravity_velocity_component = upstream_lambda_overbar_rho * this->flowfunctions_->GravityTerm(eptr) * facetNrml[v];              
            }
            
            if( this->with_capillary_spreading_) {
                capillary_velocity_component_n=vn_capillary_component_of_velocity/this->flowfunctions_->Mobility(eptr, 0U)*upstream_lambda_overbar_rho;
                capillary_velocity_component_w=vw_capillary_component_of_velocity/this->flowfunctions_->Mobility(eptr, 1U)*upstream_lambda_overbar_rho;   
            }
                      
            double64 fn = sign*(viscous_velocity_component_n-gravity_velocity_component-capillary_velocity_component_n)*facetArea; // kg/s
            double64 fw = sign*(viscous_velocity_component_w+gravity_velocity_component+capillary_velocity_component_w)*facetArea; // kg/s
            double64 f_YCO2=fn*upstream_CO2_comp[0];
            double64 f_YH2O=fn*upstream_CO2_comp[1];
            double64 f_XH2O=fw*upstream_H2O_comp[0];
            double64 f_XCO2=fw*upstream_H2O_comp[1];
            double64 f_XNaCl=fw*upstream_H2O_comp[2];
            
            YCO2_accumulation += f_YCO2; // kg/s
            YH2O_accumulation += f_YH2O; // kg/s
            XH2O_accumulation += f_XH2O; // kg/s
            XCO2_accumulation += f_XCO2; // kg/s
            XNaCl_accumulation += f_XNaCl; // kg/s
            
            YCO2_inflow += f_YCO2; // kg/s
            YH2O_inflow += f_YH2O; // kg/s
            XH2O_inflow += f_XH2O; // kg/s
            XCO2_inflow += f_XCO2; // kg/s
            XNaCl_inflow += f_XNaCl; // kg/s            
          
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
   
                if ( YCO2_inflow > 0. ) YCO2_accumulation += YCO2_inflow; //inflow compensation
                else if ( YCO2_inflow < 0. ) YCO2_accumulation -= YCO2_inflow; //outflow compensation  
                
                if ( YH2O_inflow > 0. ) YH2O_accumulation += YH2O_inflow; //inflow compensation
                else if ( YH2O_inflow < 0. ) YH2O_accumulation -= YH2O_inflow; //outflow compensation 
                
                if ( XH2O_inflow > 0. ) XH2O_accumulation += XH2O_inflow; //inflow compensation
                else if ( XH2O_inflow < 0. ) XH2O_accumulation -= XH2O_inflow; //outflow compensation 
                
                if ( XCO2_inflow > 0. ) XCO2_accumulation += XCO2_inflow; //inflow compensation
                else if ( XCO2_inflow < 0. ) XCO2_accumulation -= XCO2_inflow; //outflow compensation      
                
                if ( XNaCl_inflow > 0. ) XNaCl_accumulation += XNaCl_inflow; //inflow compensation
                else if ( XNaCl_inflow < 0. ) XNaCl_accumulation -= XNaCl_inflow; //outflow compensation         
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
        double64 fw_avg = 0.;
        double64 sw = 1. - nd->Read(this->key_sCO2); //saturation aqueous phase at current node
        ArrayVariable CO2_comp, H2O_comp;
        nd->Read(this->key_CO2_comp, CO2_comp); 
        nd->Read(this->key_H2O_comp, H2O_comp);
        for ( size_t t=0U; t<node_parent_elements; t++ )
        {
            Element<dim>* const eptr(nd->Parent(t));
            this->flowfunctions_->InitialiseBrooksCoreyParameters(eptr);   
            fn_avg += this->flowfunctions_->f_at(eptr,1U,sw);
            fw_avg += this->flowfunctions_->f_at(eptr,0U,sw);
        }
        
        fn_avg /= static_cast<double64>(node_parent_elements);
        fw_avg /= static_cast<double64>(node_parent_elements);
        double64 rhon = nd->Read(this->key_rhoCO2);
        double64 rhow = nd->Read(this->key_rhoH2O);
        
        YCO2_accumulation -= fn_avg*flux_balance*rhon*CO2_comp[0]; // kg/s
        YH2O_accumulation -= fn_avg*flux_balance*rhon*CO2_comp[1]; // kg/s
        XH2O_accumulation -= fw_avg*flux_balance*rhow*H2O_comp[0]; // kg/s
        XCO2_accumulation -= fw_avg*flux_balance*rhow*H2O_comp[1]; // kg/s
        XNaCl_accumulation -= fw_avg*flux_balance*rhow*H2O_comp[2]; // kg/s
    } 
    
    
    //compute and store variation rate    
    ArrayVariable dcm_CO2, dcm_H2O;
    nd->Read(key_dcmCO2, dcm_CO2);
    nd->Read(key_dcmH2O, dcm_H2O);
    dcm_CO2.Component(0, YCO2_accumulation);
    dcm_CO2.Component(1, YH2O_accumulation);
    dcm_H2O.Component(0, XH2O_accumulation);
    dcm_H2O.Component(1, XCO2_accumulation);
    dcm_H2O.Component(2, XNaCl_accumulation);
    nd->Store(key_dcmCO2, dcm_CO2);
    nd->Store(key_dcmH2O, dcm_H2O);    
  }
}  



//schedule an event associated with a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
bool TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::Schedule(Event<dim>* event, double64 t_end)
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
    
    ArrayVariable dcm_CO2;
    nd->Read(key_dcmCO2, dcm_CO2);
    double64 YCO2_rate = dcm_CO2[0]; //mass variaition rate of CO2 carbonic phase (kg/s)
    double64 YH2O_rate = dcm_CO2[1]; //mass variaition rate of H2O carbonic phase (kg/s)
    
    double64 rhoCO2 = nd->Read(this->key_rhoCO2); //density carbonic phase (kg/m3)
    double64 PV = nd->Read(this->key_fvPV); //pore volume (m3)
    double64 sCO2_rate = (YCO2_rate+YH2O_rate)/rhoCO2/PV; //variaition rate of saturation carbonic phase (m3/(m3.s))
    
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
void TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double64 t_clock)
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
    double64 rhoH2O = nd->Read(this->key_rhoH2O); //density aqueous phase (kg/m3)
    double64 sCO2_old = nd->Read( this->key_sCO2);//old saturation carbonic phase (m3/m3)
    double64 sH2O_old = nd->Read( this->key_sH2O);//old saturation aqeous phase (m3/m3)
    nd->Store(this->key_sCO2_0, makeScalar( status, sCO2_old )); //store old saturation  
      
    double64 mCO2_old = PV * sCO2_old * rhoCO2; //old mass carbonic phase (kg)
    double64 mH2O_old = PV * sH2O_old * rhoH2O; //old mass aqueous phase (kg)

    ArrayVariable CO2_comp, H2O_comp, dcm_CO2, dcm_H2O;
    nd->Read(this->key_CO2_comp, CO2_comp); 
    nd->Read(this->key_H2O_comp, H2O_comp);
    
    double64 YCO2 = CO2_comp[0]; //mass fraction CO2 carbonic phase
    double64 YH2O = CO2_comp[1]; //mass fraction H2O carbonic phase
    double64 XH2O = H2O_comp[0]; //mass fraction H2O aqueous phase
    double64 XCO2 = H2O_comp[1]; //mass fraction CO2 aqueous phase
    double64 XNaCl = H2O_comp[2]; //mass fraction NaCl aqueous phase    
    
    double64 mYCO2 = mCO2_old * YCO2; //mass of CO2 carbonic phase (kg)
    double64 mYH2O = mCO2_old * YH2O; //mass of H2O carbonic phase (kg)
    double64 mXH2O = mH2O_old * XH2O; //mass of H2O aqueous phase (kg)
    double64 mXCO2 = mH2O_old * XCO2; //mass of CO2 aqueous phase (kg)
    double64 mXNaCl = mH2O_old * XNaCl; //mass of NaCl aqueous phase (kg)
        
    nd->Read(key_dcmCO2, dcm_CO2);
    nd->Read(key_dcmH2O, dcm_H2O);    
    double64 YCO2_rate = dcm_CO2[0]; //mass variaition rate of CO2 carbonic phase (kg/s)
    double64 YH2O_rate = dcm_CO2[1]; //mass variaition rate of H2O carbonic phase (kg/s)
    double64 XH2O_rate = dcm_H2O[0]; //mass variaition rate of H2O aqueous phase (kg/s)
    double64 XCO2_rate = dcm_H2O[1]; //mass variaition rate of CO2 aqueous phase (kg/s)
    double64 XNaCl_rate = dcm_H2O[2]; //mass variaition rate of NaCl aqueous phase (kg/s)        
    
    double64 t_current = array[0]; //current time stamp (s)

    mYCO2 -= YCO2_rate * (t_clock - t_current); // updated mass of CO2 carbonic phase (kg)
    mYH2O -= YH2O_rate * (t_clock - t_current); //updated mass of H2O carbonic phase (kg)
    mXH2O -= XH2O_rate * (t_clock - t_current); //updated mass of H2O aqueous phase (kg)
    mXCO2 -= XCO2_rate * (t_clock - t_current); //updated mass of CO2 aqueous phase (kg)
    mXNaCl -= XNaCl_rate * (t_clock - t_current); //updated mass of NaCl aqueous phase (kg)
    
    double64 mCO2_new = mYCO2 + mYH2O; //updated mass carbonic phase
    double64 mH2O_new = mXH2O + mXCO2 + mXNaCl; //updated mass aqueous phase
    
    if(mCO2_new != 0.) {
        YCO2 = mYCO2 / mCO2_new; //updated mass fraction CO2 carbonic phase
        YH2O = mYH2O / mCO2_new; //updated mass fraction H2O carbonic phase
        CO2_comp.Component(0, YCO2);
        CO2_comp.Component(1, YH2O);
    }
        
    if(mH2O_new != 0.) {    
        XH2O = mXH2O / mH2O_new; //updated mass fraction H2O aqueous phase
        XCO2 = mXCO2 / mH2O_new; //updated mass fraction CO2 aqueous phase
        XNaCl = mXNaCl / mH2O_new; //updated mass fraction NaCl aqueous phase
        H2O_comp.Component(0, XH2O);
        H2O_comp.Component(1, XCO2);
        H2O_comp.Component(2, XNaCl);
    }   
    
    //compute new saturations from updated mass
    double64 sCO2_new = mCO2_new/rhoCO2/PV;

    double64 source = nd->Read(this->key_nQV); // nodal fluid volume source (m3/(m3.s))
    sCO2_new += source * (t_clock - t_current); //add source (m3/m3)
        
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
void TwoPhaseMultiComponentDESTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double64 delta_t)
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
    double64 rhoH2O = nd->Read(this->key_rhoH2O); //density aqueous phase (kg/m3) 
    double64 sCO2_old = nd->Read( this->key_sCO2);//old saturation carbonic phase (m3/m3)
    double64 sH2O_old = nd->Read( this->key_sH2O);//old saturation aqeous phase (m3/m3) 
    nd->Store(this->key_sCO2_0, makeScalar( status, sCO2_old )); //store old saturation  
      
    double64 mCO2_old = PV * sCO2_old * rhoCO2; //old mass carbonic phase (kg)
    double64 mH2O_old = PV * sH2O_old * rhoH2O; //old mass aqueous phase (kg) 

    ArrayVariable CO2_comp, H2O_comp, dcm_CO2, dcm_H2O;
    nd->Read(this->key_CO2_comp, CO2_comp); 
    nd->Read(this->key_H2O_comp, H2O_comp);
    
    double64 YCO2 = CO2_comp[0]; //mass fraction CO2 carbonic phase
    double64 YH2O = CO2_comp[1]; //mass fraction H2O carbonic phase
    double64 XH2O = H2O_comp[0]; //mass fraction H2O aqueous phase
    double64 XCO2 = H2O_comp[1]; //mass fraction CO2 aqueous phase
    double64 XNaCl = H2O_comp[2]; //mass fraction NaCl aqueous phase    
    
    double64 mYCO2 = mCO2_old * YCO2; //mass of CO2 carbonic phase (kg)
    double64 mYH2O = mCO2_old * YH2O; //mass of H2O carbonic phase (kg)
    double64 mXH2O = mH2O_old * XH2O; //mass of H2O aqueous phase (kg)
    double64 mXCO2 = mH2O_old * XCO2; //mass of CO2 aqueous phase (kg)
    double64 mXNaCl = mH2O_old * XNaCl; //mass of NaCl aqueous phase (kg)
        
    nd->Read(key_dcmCO2, dcm_CO2);
    nd->Read(key_dcmH2O, dcm_H2O);    
    double64 YCO2_rate = dcm_CO2[0]; //mass variaition rate of CO2 carbonic phase (kg/s)
    double64 YH2O_rate = dcm_CO2[1]; //mass variaition rate of H2O carbonic phase (kg/s)
    double64 XH2O_rate = dcm_H2O[0]; //mass variaition rate of H2O aqueous phase (kg/s)
    double64 XCO2_rate = dcm_H2O[1]; //mass variaition rate of CO2 aqueous phase (kg/s)
    double64 XNaCl_rate = dcm_H2O[2]; //mass variaition rate of NaCl aqueous phase (kg/s)     
    
    mYCO2 -= YCO2_rate * delta_t; // updated mass of CO2 carbonic phase (kg)
    mYH2O -= YH2O_rate * delta_t; //updated mass of H2O carbonic phase (kg)
    mXH2O -= XH2O_rate * delta_t; //updated mass of H2O aqueous phase (kg)
    mXCO2 -= XCO2_rate * delta_t; //updated mass of CO2 aqueous phase (kg)
    mXNaCl -= XNaCl_rate * delta_t; //updated mass of NaCl aqueous phase (kg)
    
    double64 mCO2_new = mYCO2 + mYH2O; //updated mass carbonic phase
    double64 mH2O_new = mXH2O + mXCO2 + mXNaCl; //updated mass aqueous phase
    
    if(mCO2_new != 0.) {
        YCO2 = mYCO2 / mCO2_new; //updated mass fraction CO2 carbonic phase
        YH2O = mYH2O / mCO2_new; //updated mass fraction H2O carbonic phase
        CO2_comp.Component(0, YCO2);
        CO2_comp.Component(1, YH2O);
    }
        
    if(mH2O_new != 0.) {    
        XH2O = mXH2O / mH2O_new; //updated mass fraction H2O aqueous phase
        XCO2 = mXCO2 / mH2O_new; //updated mass fraction CO2 aqueous phase
        XNaCl = mXNaCl / mH2O_new; //updated mass fraction NaCl aqueous phase
        H2O_comp.Component(0, XH2O);
        H2O_comp.Component(1, XCO2);
        H2O_comp.Component(2, XNaCl);
    }    
    
    double64 sCO2_new = mCO2_new/rhoCO2/PV;
        
    double64 source = nd->Read(this->key_nQV); // nodal fluid volume source (m3/(m3.s))
    //mCO2_new += PV * source * rhoCO2 * delta_t; //add source (kg)
    sCO2_new += source * delta_t; //add source (m3/m3)
    
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

template class TwoPhaseMultiComponentDESTransport<1U,CO2H2O_FunctionsModule0>;
template class TwoPhaseMultiComponentDESTransport<2U,CO2H2O_FunctionsModule0>;
template class TwoPhaseMultiComponentDESTransport<3U,CO2H2O_FunctionsModule0>;

template class TwoPhaseMultiComponentDESTransport<1U,CO2H2O_FunctionsModule1>;
template class TwoPhaseMultiComponentDESTransport<2U,CO2H2O_FunctionsModule1>;
template class TwoPhaseMultiComponentDESTransport<3U,CO2H2O_FunctionsModule1>;

template class TwoPhaseMultiComponentDESTransport<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseMultiComponentDESTransport<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseMultiComponentDESTransport<3U,CO2H2O_FunctionsModule2>;

} // end csmp 



                                                                       
