#include "TwoPhaseTwoComponentDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "DenseMatrix.h"
#include "CSMP_mathUtilities.h"
#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim>
TwoPhaseTwoComponentDESTransport<dim>::TwoPhaseTwoComponentDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 FlowFunctions<dim>& flowfunctions, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier)
    : TwoPhaseDESTransport<dim> (m,target_region,flowfunctions,with_capillary_spreading,with_gravity_forces,PEP_multiplier,cfl_multiplier) 
{
    initializeVariablsAndKeys(m);
    cout<<"TwoPhaseTwoComponentDESTransport constructed"<<endl;
} // end constructor  


template<size_t dim>
TwoPhaseTwoComponentDESTransport<dim>::TwoPhaseTwoComponentDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 FlowFunctions<dim>& flowfunctions, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces,
                                                 double64 PEP_multiplier,
                                                 double64 cfl_multiplier,
                                                 double64 relaxing_factor)
    : TwoPhaseDESTransport<dim> (m,target_region,flowfunctions,with_capillary_spreading,with_gravity_forces,PEP_multiplier,cfl_multiplier, relaxing_factor)
{
    initializeVariablsAndKeys(m);
    cout<<"TwoPhaseTwoComponentDESTransport constructed"<<endl;
} // end constructor 

  

template<size_t dim>
void TwoPhaseTwoComponentDESTransport<dim>::initializeVariablsAndKeys(Model<dim>& m)
{
    if(!m.Database().IsDefined("component mass variation rates array")) m.CreateProperty( "component mass variation rates array", "kg/s", ARRAY, NODE, 4, -1.00E+10,1.00E+10);   
    this->key_components = INDEX<ARRAY,NODE>( m.Database().StorageKey("component mass variation rates array") );
    
    if ( this->key_components.place != NODE || this->key_components.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseTwoComponentDESTransport::initializeKeys:",
        "The 'component mass variation rates array' variable must be ARRAY and placed on NODE"  ); 
        
    const typename vector<Node<dim>*>::const_iterator  nodes_end(this->gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=this->gref_.NodesBegin(); nit!=nodes_end; ++nit )
    {         
        ArrayVariable arrayVariable( 7, 0., PLAIN );         
        (*nit)->Store( this->key_time, arrayVariable );
    }                                                 
}



//Compute the rate of change of non-wetting phase in a node/FV, as well as the mass variation rate for each of components
template<size_t dim>
void TwoPhaseTwoComponentDESTransport<dim>::ComputeRateofChange( Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH);

    this->rate_count_++;//recording
    nd->Store(  this->key_rate, makeScalar( nd->Status( this->key_rate), nd->Read( this->key_rate) + 1 ) );

    double64 accumulation(0.), flux_balance(0.), outflow(0.);
    double64 accumulation_XCO2(0.), accumulation_XH2O(0.), accumulation_YCO2(0.), accumulation_YH2O(0.);
        
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD;
    
    double64 cfl_multiplier = this->CFL_multiplier_*this->relaxing_factor_; //default value
    size_t index = 0;
    for (auto fip : nd->AllFacetIntegrationPoints()) {
        const double64 sign = fip.FromInside() ? 1. : -1.;
        //compute facet flux
        fip.Obtain( this->key_vt, vD );
        //double64 vD_n = fip.ProjectOntoFacetNormal(vD);  
        double64 facetArea = event->facetAreaCollection[index];
        Point<dim> facetNrml = event->facetNormalCollection[index];
        index++;
        double64 vD_n = dotProduct(vD.P(), facetNrml); 
              
        const double64 facet_flux = sign * vD_n * facetArea;
        flux_balance += facet_flux;       
        if ( facet_flux > 0. ) outflow += facet_flux ; 
        
        //compute mobilities on inside and outside nodes          
        auto inside_node = fip.InsideNode();
        auto outside_node = fip.OutsideNode(); 
       
        const double64 sn_inside_node = inside_node.Read( this->key_sCO2 );
        const double64 sw_inside_node = 1.-sn_inside_node;
        const double64 ln_inside_node = this->flowfunctions_.Mobility_at(fip,1U,1.0-sn_inside_node); 
        const double64 lw_inside_node = this->flowfunctions_.Mobility_at(fip,0U,1.0-sn_inside_node); 
        const double64 rhon_inside_node = inside_node.Read( this->key_rhoCO2 ); //density carbonic phase
        const double64 rhow_inside_node = inside_node.Read( this->key_rhoH2O ); //density aqueous phase
        const double64 XCO2_inside_node = inside_node.Read( this->key_XCO2 ); //mass fraction CO2 aqueous phase
        const double64 XH2O_inside_node = inside_node.Read( this->key_XH2O ); //mass fraction H2O aqueous phase
        const double64 YCO2_inside_node = inside_node.Read( this->key_YCO2 ); //mass fraction CO2 carbonic phase
        const double64 YH2O_inside_node = inside_node.Read( this->key_YH2O ); //mass fraction H2O carbonic phase
        
        const double64 sn_outside_node = outside_node.Read( this->key_sCO2 );
        const double64 sw_outside_node = 1.-sn_outside_node;
        const double64 ln_outside_node = this->flowfunctions_.Mobility_at(fip,1U,1.0-sn_outside_node); 
        const double64 lw_outside_node = this->flowfunctions_.Mobility_at(fip,0U,1.0-sn_outside_node);      
        const double64 rhon_outside_node = outside_node.Read( this->key_rhoCO2 ); //density carbonic phase
        const double64 rhow_outside_node = outside_node.Read( this->key_rhoH2O ); //density aqueous phase
        const double64 XCO2_outside_node = outside_node.Read( this->key_XCO2 ); //mass fraction CO2 aqueous phase
        const double64 XH2O_outside_node = outside_node.Read( this->key_XH2O ); //mass fraction H2O aqueous phase
        const double64 YCO2_outside_node = outside_node.Read( this->key_YCO2 ); //mass fraction CO2 carbonic phase
        const double64 YH2O_outside_node = outside_node.Read( this->key_YH2O ); //mass fraction H2O carbonic phase
         
        //compute velocities at facet integration point                      
        double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
        double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
        if( this->with_gravity_forces_ ){                       
            vn_gravity_component_of_velocity = this->flowfunctions_.Mobility(fip, 0U) * this->flowfunctions_.GravityTerm(fip) * facetNrml[v];
            vw_gravity_component_of_velocity = this->flowfunctions_.Mobility(fip, 1U) * this->flowfunctions_.GravityTerm(fip) * facetNrml[v];
        }         
            
        if(this->with_capillary_spreading_){  
            VectorVariable<dim> grad;
            fip.Read(this->key_grad, grad);
            Point<dim> snw_gradient=grad.P();       
            
            double64 dsdn = dotProduct(snw_gradient, facetNrml);
            
            if(!isnan(dsdn)){
                vn_capillary_component_of_velocity = -dsdn*this->flowfunctions_.CapillaryDiffusionMultiplier_Phase(fip,0U);
                vw_capillary_component_of_velocity = -dsdn*this->flowfunctions_.CapillaryDiffusionMultiplier_Phase(fip,1U);
            }   
        }  
            
        double64 vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
        double64 vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity;  
  
        //determine upstream mobilities
        double64 upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);  
        double64 upstream_density_n(0.0), upstream_density_w(0.0), upstream_XCO2(0.0), upstream_XH2O(0.0), upstream_YCO2(0.0), upstream_YH2O(0.0);  
              
        if(vn_at_facet_int_point>0.0) {
            upstream_mobility_n=ln_inside_node;
            upstream_density_n=rhon_inside_node; 
            upstream_YCO2=YCO2_inside_node;
            upstream_YH2O=YH2O_inside_node; 
        } else if (vn_at_facet_int_point<0.0) {
            upstream_mobility_n=ln_outside_node;
            upstream_density_n=rhon_outside_node; 
            upstream_YCO2=YCO2_outside_node; 
            upstream_YH2O=YH2O_outside_node;             
        } else { 
            upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
            upstream_density_n=0.5*(rhon_inside_node+rhon_outside_node);
            upstream_YCO2=0.5*(YCO2_inside_node+YCO2_outside_node);
            upstream_YH2O=0.5*(YH2O_inside_node+YH2O_outside_node);
        }    
            
        if(vw_at_facet_int_point>0.0) {
            upstream_mobility_w=lw_inside_node;
            upstream_density_w=rhow_inside_node; 
            upstream_XCO2=XCO2_inside_node; 
            upstream_XH2O=XH2O_inside_node;          
        } else if (vw_at_facet_int_point<0.0) {
            upstream_mobility_w=lw_outside_node;
            upstream_density_w=rhow_outside_node; 
            upstream_XCO2=XCO2_outside_node; 
            upstream_XH2O=XH2O_outside_node;              
        } else {
            upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);
            upstream_density_w=0.5*(rhow_inside_node+rhow_outside_node); 
            upstream_XCO2=0.5*(XCO2_inside_node+XCO2_outside_node); 
            upstream_XH2O=0.5*(XH2O_inside_node+XH2O_outside_node); 
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
            gravity_velocity_component = upstream_lambda_overbar * this->flowfunctions_.GravityTerm(fip) * facetNrml[v];
            gravity_velocity_component_w = upstream_lambda_overbar * this->flowfunctions_.GravityTerm(fip) * facetNrml[v];
        }

        if( this->with_capillary_spreading_ ) {
            capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity; 
            capillary_velocity_component_w = upstream_fw * vw_capillary_component_of_velocity;
        }
        
        //update non-wetting flux accumulation 
        double64 fn = sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;
        accumulation += fn;
        accumulation_YCO2 +=  fn * upstream_density_n * upstream_YCO2/100.; 
        accumulation_YH2O +=  fn * upstream_density_n * upstream_YH2O/100.; 
        
        double64 fw = sign * ( viscous_velocity_component_w + gravity_velocity_component_w + capillary_velocity_component_w) * facetArea;
        accumulation_XCO2 +=  fw * upstream_density_w * upstream_XCO2/100.; 
        accumulation_XH2O +=  fw * upstream_density_w * upstream_XH2O/100.;          
        
        //determine cfl_multiplier based on non-wetting phase shock saturation
        if (cfl_multiplier != this->CFL_multiplier_){
            if(vn_at_facet_int_point < 0.0) { //flowing in from outside node (upstream node)
                double64 sn_outside_shock = outside_node.Read(this->key_ssn); //sn at shock for outside node
                if (sn_outside_node >= sn_outside_shock) { //upstream node passed shock saturation
                    double64 sn_inside_shock = inside_node.Read(this->key_ssn); //sn at shock for inside node (current node)
                    if (sn_inside_node < sn_inside_shock) {//current node not yet reach shock saturation   
                        cfl_multiplier = this->CFL_multiplier_;
                    }
                }
            }
        }
        
        //determine cfl_multiplier based on wetting phase shock saturation
        if (cfl_multiplier != this->CFL_multiplier_){
            if(vw_at_facet_int_point < 0.0) { //flowing from outside node (upstream node)
            //if(vD_n < 0.0) { //flowing from outside node (upstream node)
                double64 sw_outside_shock = outside_node.Read(this->key_ssw); //sw at shock for outside node
                if (sw_outside_node >= sw_outside_shock) { //upstream node passed shock saturation
                    double64 sw_inside_shock = inside_node.Read(this->key_ssw); //sn at shock for inside node (current node)
                    if (sw_inside_node < sw_inside_shock) {//current node not yet reach shock saturation   
                        cfl_multiplier = this->CFL_multiplier_;
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
        accumulation_YCO2 -= fn_avg*flux_balance*nd->Read(this->key_rhoCO2)*nd->Read(this->key_YCO2)/100.; 
        accumulation_YH2O -= fn_avg*flux_balance*nd->Read(this->key_rhoCO2)*nd->Read(this->key_YH2O)/100.; 
        accumulation_XCO2 -= fw_avg*flux_balance*nd->Read(this->key_rhoH2O)*nd->Read(this->key_XCO2)/100.; 
        accumulation_XH2O -= fw_avg*flux_balance*nd->Read(this->key_rhoH2O)*nd->Read(this->key_XH2O)/100.;
        
    }     
    
    //compute and store variation rate    
    double64 PV = nd->Read(  this->key_fvPV );
    nd->Store( this->key_dsnw, makeScalar( nd->Status( this->key_dsnw ), accumulation/PV ) );    
    
    ArrayVariable component_array;
    nd->Read(this->key_components, component_array);
    component_array.Component(0, accumulation_YCO2);
    component_array.Component(1, accumulation_YH2O);
    component_array.Component(2, accumulation_XCO2);
    component_array.Component(3, accumulation_XH2O);    
    nd->Store(this->key_components, component_array);
}  



//update solution and check it against the specified range (with DES). Also update mass fraction of each
template<size_t dim>
void TwoPhaseTwoComponentDESTransport<dim>::Update_DES(Event<dim>* event, double64 t_clock)
{
    this->update_count_++;//recording
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );
    assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    ArrayVariable array;
    nd->Read(this->key_time, array);

    double64 ChangeRate = nd->Read( this->key_dsnw);//variaition rate   
    double64 solution = nd->Read( this->key_sCO2);//old solution
    double64 t_current = array[0]; //current time stamp
    double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
    const double64 source(nd->Read(this->key_nQV));
    new_solution += source * (t_clock - t_current);//add source to new solution
    
    double64 mYCO2 = nd->Read(this->key_rhoCO2) * nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2) * nd->Read(this->key_YCO2) / 100.;
    double64 mYH2O = nd->Read(this->key_rhoCO2) * nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2) * nd->Read(this->key_YH2O) / 100.;
    double64 mXCO2 = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O) * nd->Read(this->key_XCO2) / 100.;
    double64 mXH2O = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O) * nd->Read(this->key_XH2O) / 100.;
    //double64 mXSalt = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O) * (100. - nd->Read(this->key_XCO2) - nd->Read(this->key_XH2O)) / 100.;
    
    //update mass fraction for each component
    ArrayVariable component_array;
    nd->Read(this->key_components, component_array);    
    mYCO2 -= component_array[0] * (t_clock - t_current);
    mYH2O -= component_array[1] * (t_clock - t_current);
    mXCO2 -= component_array[2] * (t_clock - t_current);
    mXH2O -= component_array[3] * (t_clock - t_current);  
    //double64 mYtot = nd->Read(this->key_rhoCO2) * nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2);
    //double64 mXtot = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O);
    double64 mYtot = mYCO2 + mYH2O;
    double64 mXtot = mXCO2 + mXH2O;    
    nd->Store(this->key_YCO2, makeScalar( nd->Status(this->key_YCO2), mYCO2/mYtot*100.));
    nd->Store(this->key_YH2O, makeScalar( nd->Status(this->key_YH2O), mYH2O/mYtot*100.));
    nd->Store(this->key_XCO2, makeScalar( nd->Status(this->key_XCO2), mXCO2/mXtot*100.));
    nd->Store(this->key_XH2O, makeScalar( nd->Status(this->key_XH2O), mXH2O/mXtot*100.));
        
    //check new solution value against range and stored it to this->key_sCO2
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



//update solution and check it against the specified range (with TDS)
template<size_t dim>
void TwoPhaseTwoComponentDESTransport<dim>::Update_TDS(Event<dim>* event, double64 delta_t)
{
    this->update_count_++;//recording
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );   
    assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    
    double64 ChangeRate = nd->Read(this->key_dsnw);//variation rate  
    double64 solution = nd->Read(this->key_sCO2);//old solution
    double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
    const double64 source(nd->Read( this->key_nQV));
    new_solution += source * delta_t;//add source to new solution.
    
    //update mass fraction for each of the components
    double64 mYCO2 = nd->Read(this->key_rhoCO2) * nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2) * nd->Read(this->key_YCO2) / 100.;
    double64 mYH2O = nd->Read(this->key_rhoCO2) * nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2) * nd->Read(this->key_YH2O) / 100.;
    double64 mXCO2 = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O) * nd->Read(this->key_XCO2) / 100.;
    double64 mXH2O = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O) * nd->Read(this->key_XH2O) / 100.;
    //double64 mXSalt = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O) * (100. - nd->Read(this->key_XCO2) - nd->Read(this->key_XH2O)) / 100.;
    
    ArrayVariable component_array;
    nd->Read(this->key_components, component_array);    
    mYCO2 -= component_array[0] * delta_t;
    mYH2O -= component_array[1] * delta_t;
    mXCO2 -= component_array[2] * delta_t;
    mXH2O -= component_array[3] * delta_t;  
    //double64 mYtot = nd->Read(this->key_rhoCO2) * nd->Read(this->key_fvPV) * nd->Read(this->key_sCO2);
    //double64 mXtot = nd->Read(this->key_rhoH2O) * nd->Read(this->key_fvPV) * nd->Read(this->key_sH2O);
    double64 mYtot = mYCO2 + mYH2O;
    double64 mXtot = mXCO2 + mXH2O;    
    nd->Store(this->key_YCO2, makeScalar( nd->Status(this->key_YCO2), mYCO2/mYtot*100.));
    nd->Store(this->key_YH2O, makeScalar( nd->Status(this->key_YH2O), mYH2O/mYtot*100.));
    nd->Store(this->key_XCO2, makeScalar( nd->Status(this->key_XCO2), mXCO2/mXtot*100.));
    nd->Store(this->key_XH2O, makeScalar( nd->Status(this->key_XH2O), mXH2O/mXtot*100.));
                
    //check new solution value against range and stored it to this->key_sCO2
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


template class TwoPhaseTwoComponentDESTransport<1U>;
template class TwoPhaseTwoComponentDESTransport<2U>;
template class TwoPhaseTwoComponentDESTransport<3U>;

} // end csmp 



                                                                       
