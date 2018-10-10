#include "TwoPhaseTwoComponentDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "CO2H2O_FunctionsModule1.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseTwoComponentDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseTwoComponentDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 bool tensor_k,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier)
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,with_capillary_spreading,with_gravity_forces,tensor_k,PEP_multiplier,cfl_multiplier) 
{
    InitializeVariablesAndKeys(m);
    m.Database().RangeOf( m.Database().Name(this->key_CO2aq), lower_CO2aq_, upper_CO2aq_ );
    m.Database().RangeOf( m.Database().Name(this->key_H2Og), lower_H2Og_, upper_H2Og_ );    
    cout<<"TwoPhaseTwoComponentDESTransport constructed."<<endl;
} // end constructor  






template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseTwoComponentDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseTwoComponentDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 bool tensor_k,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier,
                                                 double64 relaxing_factor)
    : TwoPhaseDESTransport<dim,FLOW_FUNCTIONS> (m,target_region,with_capillary_spreading,with_gravity_forces,tensor_k,PEP_multiplier,cfl_multiplier, relaxing_factor)
{
    InitializeVariablesAndKeys(m);
    m.Database().RangeOf( m.Database().Name(this->key_CO2aq), lower_CO2aq_, upper_CO2aq_ );
    m.Database().RangeOf( m.Database().Name(this->key_H2Og), lower_H2Og_, upper_H2Og_ );     
    cout<<"TwoPhaseTwoComponentDESTransport constructed."<<endl;
} // end constructor 

  



template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseTwoComponentDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys(Model<dim>& m)
{
    if(!m.Database().IsDefined("component mass variation rates array")) m.CreateProperty( "component mass variation rates array", "kg/s", ARRAY, NODE, 2, -1.00E+10,1.00E+10);   
    this->key_components = INDEX<ARRAY,NODE>( m.Database().StorageKey("component mass variation rates array") );
    
    if ( this->key_components.place != NODE || this->key_components.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseTwoComponentDESTransport::initializeKeys:",
        "The 'component mass variation rates array' variable must be ARRAY and placed on NODE"  ); 
        
    const typename vector<Node<dim>*>::const_iterator  nodes_end(this->gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=this->gref_.NodesBegin(); nit!=nodes_end; ++nit )
    {         
        ArrayVariable arrayVariable( 2, 0., PLAIN );         
        (*nit)->Store( this->key_components, arrayVariable );
    }                                                 
}


//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseTwoComponentDESTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){ 
    this->rate_count_++;//recording
    nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );
    
    //FLOW_FUNCTIONS<dim> flowfunctions(this->db_);

    double64 accumulation(0.), flux_balance(0.), outflow(0.);
    double64 accumulation_CO2aq(0.), accumulation_H2Og(0.);
        
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD, facetNrml;
    const size_t node_parent_elements(nd->Parents());
    
    double64 cfl_multiplier = this->CFL_multiplier_*this->relaxing_factor_; //default value

    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        
        this->flowfunctions_->InitialiseBrooksCoreyParameters(eptr); 
        
        const size_t pnid(nd->ParentNodeNumber(t));
        eptr->Read( this->key_vt, vD);

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
            //update outflow
            if ( facet_flux > 0. ) outflow += facet_flux;
                         
            //compute inside and outside node mobilities, by using their saturations
            const double64 sn_inside_node = eptr->N(inside_node)->Read(  this->key_sCO2 );
            const double64 sw_inside_node = 1.-sn_inside_node;
            const double64 ln_inside_node = this->flowfunctions_->Mobility_at(eptr,1U,1.0-sn_inside_node); 
            const double64 lw_inside_node = this->flowfunctions_->Mobility_at(eptr,0U,1.0-sn_inside_node); 
            const double64 CO2aq_inside_node = eptr->N(inside_node)->Read( this->key_CO2aq ); //dissolved CO2
            const double64 H2Og_inside_node = eptr->N(inside_node)->Read( this->key_H2Og ); //evaporated water            
        
            const double64 sn_outside_node = eptr->N(outside_node)->Read(  this->key_sCO2 );
            const double64 sw_outside_node = 1.-sn_outside_node;
            const double64 ln_outside_node = this->flowfunctions_->Mobility_at(eptr,1U,1.0-sn_outside_node); 
            const double64 lw_outside_node = this->flowfunctions_->Mobility_at(eptr,0U,1.0-sn_outside_node);   
            const double64 CO2aq_outside_node = eptr->N(outside_node)->Read( this->key_CO2aq ); //dissolved CO2
            const double64 H2Og_outside_node = eptr->N(outside_node)->Read( this->key_H2Og ); //evaporated water            
            
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
                double64 dsdn = grad.DotProduct(facetNrml);
            
                if(!isnan(dsdn)){
                    vn_capillary_component_of_velocity = -dsdn*this->flowfunctions_->CapillaryDiffusionMultiplier_Phase(eptr,0U);
                    vw_capillary_component_of_velocity = -dsdn*this->flowfunctions_->CapillaryDiffusionMultiplier_Phase(eptr,1U);
                }   
            }  
            
            double64 vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            double64 vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity; 
        
            //determine upstream mobilities
            double64 upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);   
            double64 upstream_CO2aq(0.0), upstream_H2Og(0.0);  
            
            if(vn_at_facet_int_point>0.0) {
                upstream_mobility_n=ln_inside_node;
                upstream_H2Og=H2Og_inside_node; 
            } else if (vn_at_facet_int_point<0.0) {
                upstream_mobility_n=ln_outside_node;
                upstream_H2Og=H2Og_outside_node;
            } else { 
                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                upstream_H2Og=0.5*(H2Og_inside_node+H2Og_outside_node);
            }    
            
            if(vw_at_facet_int_point>0.0) {
                upstream_mobility_w=lw_inside_node;
                upstream_CO2aq=CO2aq_inside_node; 
            } else if (vw_at_facet_int_point<0.0) {
                upstream_mobility_w=lw_outside_node;
                upstream_CO2aq=CO2aq_outside_node;  
            } else {
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);
                upstream_CO2aq=0.5*(CO2aq_inside_node+CO2aq_outside_node);
            }
                            
            total_mobility=upstream_mobility_n+upstream_mobility_w;
            double64 upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            double64 upstream_fw=(total_mobility!=0.0? upstream_mobility_w/total_mobility : 0.0); 
            double64 upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);
            
            //compute each velocity component     
            double64 viscous_velocity_component(0.0), capillary_velocity_component(0.0), gravity_velocity_component(0.0); 
            double64 viscous_velocity_component_w(0.0), capillary_velocity_component_w(0.0), gravity_velocity_component_w(0.0);  
                     
            viscous_velocity_component = vD_n * upstream_fn;
            viscous_velocity_component_w = vD_n * upstream_fw;
                                       
            if( this->with_gravity_forces_ ) {
                gravity_velocity_component = upstream_lambda_overbar * this->flowfunctions_->GravityTerm(eptr) * facetNrml[v];
                gravity_velocity_component_w = upstream_lambda_overbar * this->flowfunctions_->GravityTerm(eptr) * facetNrml[v];
            }

            if( this->with_capillary_spreading_ ) {
                capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity; 
                capillary_velocity_component_w = upstream_fw * vw_capillary_component_of_velocity;
            }
                
            //update non-wetting flux accumulation 
            double64 fn = sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;
            accumulation += fn;
            accumulation_H2Og +=  fn * upstream_H2Og; 
        
            double64 fw = sign * ( viscous_velocity_component_w + gravity_velocity_component_w + capillary_velocity_component_w) * facetArea;
            accumulation_CO2aq +=  fw * upstream_CO2aq;

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
        }
    }
        
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
    
     /* 
    // divergence free correction
    // compute average fractional flow for the current finite volume
    if (fabs(flux_balance) > numeric_limits<double64>::epsilon()) {
        double64 fn_avg = 0.;
        double64 fw_avg = 0.; 
        size_t number = 0;
        for (auto sip : nd->AllSectorIntegrationPoints()) {
            fn_avg += this->flowfunctions_.f(sip, 1U);
            fw_avg += this->flowfunctions_.f(sip, 0U); 
            number ++;
        }
        fn_avg /= double64(number);
        fw_avg /= double64(number); 
        accumulation -= fn_avg*flux_balance;
        accumulation_H2Og -= fn_avg*flux_balance*nd->Read(this->key_H2Og); 
        accumulation_CO2aq -= fw_avg*flux_balance*nd->Read(this->key_CO2aq); 
        
    } 
    */   
    
    //compute and store variation rate    
    double64 PV = nd->Read(  this->key_fvPV );
    nd->Store( this->key_dsnw, makeScalar( nd->Status( this->key_dsnw ), accumulation/PV ) );    
    
    ArrayVariable component_array;
    nd->Read(this->key_components, component_array);
    component_array.Component(0, accumulation_CO2aq);
    component_array.Component(1, accumulation_H2Og);
    nd->Store(this->key_components, component_array);    
  }
}  




//update solution and check it against the specified range (with DES). Also update mass fraction of each
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseTwoComponentDESTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double64 t_clock)
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
    nd->Store(this->key_sCO2_0, makeScalar( status, solution )); //store old solution
    double64 t_current = array[0]; //current time stamp
    double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
    const double64 source(nd->Read(this->key_nQV));
    new_solution += source * (t_clock - t_current);//add source to new solution
    
    double64 mCO2aq = nd->Read(this->key_CO2aq)* nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O);
    double64 mH2Og = nd->Read(this->key_H2Og)* nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2);
    
    //update mass fraction for each component
    ArrayVariable component_array;
    nd->Read(this->key_components, component_array);    
    mCO2aq -= component_array[0] * (t_clock - t_current);
    mH2Og -= component_array[1] * (t_clock - t_current);
    
    //check new solution value against range and stored it to this->key_sCO2
    if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    } else {
        cerr <<"saturation carbonic phase value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_ << endl;
        if ( new_solution > this->upper_limit_ ) new_solution = this->upper_limit_;
        else if ( new_solution < this->lower_limit_ ) new_solution = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    }  
    new_solution = nd->Read(this->key_sCO2);//stored new solution
    double64 dsn_cumulative = array[4];
    array.Component(4, dsn_cumulative + (new_solution-solution));//update cumulative change
    
    double64 CO2aq = ((1.-new_solution) <= 0.) ? 0. : mCO2aq/(nd->Read(this->key_fvPV)*(1.-new_solution));
    if ( CO2aq <= upper_CO2aq_ && CO2aq >= lower_CO2aq_  ) {
        nd->Store(this->key_CO2aq, makeScalar( nd->Status(this->key_CO2aq), CO2aq));
    } else {
        cerr <<"dissolved CO2 value: "<< CO2aq <<" versus range from PropertyDatabase: "<< lower_CO2aq_  <<"-"<< upper_CO2aq_ << endl;
        if ( CO2aq > upper_CO2aq_ ) nd->Store(this->key_CO2aq, makeScalar( nd->Status(this->key_CO2aq), upper_CO2aq_));
        else if ( CO2aq < lower_CO2aq_ ) nd->Store(this->key_CO2aq, makeScalar( nd->Status(this->key_CO2aq), lower_CO2aq_));    
    }
    
    double64 H2Og = (new_solution <= 0.) ? 0. : mH2Og/(nd->Read(this->key_fvPV)*new_solution);
    if ( H2Og <= upper_H2Og_ && H2Og >= lower_H2Og_  ) {
        nd->Store(this->key_H2Og, makeScalar( nd->Status(this->key_H2Og), H2Og));
    } else {
        cerr <<"evaporated H2O value: "<< H2Og <<" versus range from PropertyDatabase: "<< lower_H2Og_  <<"-"<< upper_H2Og_ << endl;
        if ( H2Og > upper_H2Og_ ) nd->Store(this->key_H2Og, makeScalar( nd->Status(this->key_H2Og), upper_H2Og_));
        else if ( H2Og < lower_H2Og_ ) nd->Store(this->key_H2Og, makeScalar( nd->Status(this->key_H2Og), lower_H2Og_));    
    }
        
    array.Component(0, t_clock); //current time stamp
    
    nd->Store(this->key_time, array);

    nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );
  }
}



//update solution and check it against the specified range (with TDS)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseTwoComponentDESTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double64 delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );   
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    this->update_count_++;//recording     
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    
    double64 ChangeRate = nd->Read(this->key_dsnw);//variation rate  
    double64 solution = nd->Read(this->key_sCO2);//old solution
    nd->Store(this->key_sCO2_0, makeScalar( status, solution )); //store old solution
    double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
    const double64 source(nd->Read( this->key_nQV));
    new_solution += source * delta_t;//add source to new solution.
    
    double64 mCO2aq = nd->Read(this->key_CO2aq)* nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O);
    double64 mH2Og = nd->Read(this->key_H2Og)* nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2);
    
    ArrayVariable component_array;
    nd->Read(this->key_components, component_array);    
    mCO2aq -= component_array[0] * delta_t;
    mH2Og -= component_array[1] * delta_t;  
            
    //check new solution value against range and stored it to this->key_sCO2
    if ( new_solution <= this->upper_limit_ && new_solution >= this->lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    } else {
        cerr <<"saturation carbonic phase value: "<< new_solution <<" versus range from PropertyDatabase: "<< this->lower_limit_ <<"-"<< this->upper_limit_ << endl;
        if ( new_solution > this->upper_limit_ ) new_solution = this->upper_limit_;
        else if ( new_solution < this->lower_limit_ ) new_solution = this->lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    }
    
    new_solution = nd->Read(this->key_sCO2);//stored new solution
    double64 CO2aq = ((1.-new_solution) <= 0.) ? 0. : mCO2aq/(nd->Read(this->key_fvPV)*(1.-new_solution));
    if ( CO2aq <= upper_CO2aq_ && CO2aq >= lower_CO2aq_  ) {
        nd->Store(this->key_CO2aq, makeScalar( nd->Status(this->key_CO2aq), CO2aq));
    } else {
        cerr <<"dissolved CO2 value: "<< CO2aq <<" versus range from PropertyDatabase: "<< lower_CO2aq_  <<"-"<< upper_CO2aq_ << endl;
        if ( CO2aq > upper_CO2aq_ ) nd->Store(this->key_CO2aq, makeScalar( nd->Status(this->key_CO2aq), upper_CO2aq_));
        else if ( CO2aq < lower_CO2aq_ ) nd->Store(this->key_CO2aq, makeScalar( nd->Status(this->key_CO2aq), lower_CO2aq_));    
    }
    
    double64 H2Og = (new_solution <= 0.) ? 0. : mH2Og/(nd->Read(this->key_fvPV)*new_solution);
    if ( H2Og <= upper_H2Og_ && H2Og >= lower_H2Og_  ) {
        nd->Store(this->key_H2Og, makeScalar( nd->Status(this->key_H2Og), H2Og));
    } else {
        cerr <<"evaporated H2O value: "<< H2Og <<" versus range from PropertyDatabase: "<< lower_H2Og_  <<"-"<< upper_H2Og_ << endl;
        if ( H2Og > upper_H2Og_ ) nd->Store(this->key_H2Og, makeScalar( nd->Status(this->key_H2Og), upper_H2Og_));
        else if ( H2Og < lower_H2Og_ ) nd->Store(this->key_H2Og, makeScalar( nd->Status(this->key_H2Og), lower_H2Og_));    
    }      
        
    nd->Store( this->key_update, makeScalar( nd->Status(this->key_update), nd->Read(this->key_update) + 1 ) );   
  }
}

template class TwoPhaseTwoComponentDESTransport<1U,CO2H2O_FunctionsModule0>;
template class TwoPhaseTwoComponentDESTransport<2U,CO2H2O_FunctionsModule0>;
template class TwoPhaseTwoComponentDESTransport<3U,CO2H2O_FunctionsModule0>;

template class TwoPhaseTwoComponentDESTransport<1U,CO2H2O_FunctionsModule1>;
template class TwoPhaseTwoComponentDESTransport<2U,CO2H2O_FunctionsModule1>;
template class TwoPhaseTwoComponentDESTransport<3U,CO2H2O_FunctionsModule1>;

template class TwoPhaseTwoComponentDESTransport<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseTwoComponentDESTransport<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseTwoComponentDESTransport<3U,CO2H2O_FunctionsModule2>;

} // end csmp 



                                                                       
