#include "TwoPhaseDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "DenseMatrix.h"
#include "CSMP_mathUtilities.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim>
TwoPhaseDESTransport<dim>::TwoPhaseDESTransport( Model<dim>& m, 
                                                 const char* target_region, 
                                                 FlowFunctions<dim>& flowfunctions, 
                                                 bool with_capillary_spreading, 
                                                 bool with_gravity_forces)
    : variables::Variables_TwoPhaseFlow(m.Database()),
      gref_(m.Region(target_region)),
      flowfunctions_(flowfunctions),
      with_capillary_spreading_(with_capillary_spreading),
      with_gravity_forces_(with_gravity_forces),
      upper_limit_(1.), lower_limit_(0.), rate_count_(0U), update_count_(0U),
      T_RateOfChange_(0.), T_Schedule_(0.), T_SortQueue_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromQueue_(0.), T_AdvectVariable_(0.)
{
    m.InstantiateFiniteVolumes();
    initializeVariablsAndKeys(m);
    //create events for all nodes and add them to PEPStack and EntireQueue
    int index = 0;
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        (*nit)->Store( key_EventIndex, makeScalar( (*nit)->Status(key_EventIndex), index) );//event index        
        (*nit)->Store( key_update, makeScalar( (*nit)->Status( key_update), 0 ) ); //update count
        (*nit)->Store( key_rate, makeScalar( (*nit)->Status( key_rate), 0 ) ); //changerate count
        (*nit)->Store( key_schedule, makeScalar( (*nit)->Status( key_schedule), 0 ) ); //schedule count
        (*nit)->Store( key_synchronize, makeScalar( (*nit)->Status( key_synchronize), 0 ) ); //synchronize count   
        
        ArrayVariable arrayVariable( 6, 0., PLAIN );         
        (*nit)->Store( key_time, arrayVariable );  
        
        initializeFiniteVolumeProperties(*nit);           
        Event<dim>* event = new Event<dim>(*nit);
        event->inPEPStack(true);
        event->valid(false);
        PEPStack.push_back(event); //add event to PEPStack
        EntireQueue.push_back(event); //add event to EntireQueue
        index++;
    }
         
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_sCO2), lower_limit_, upper_limit_ );
    cout<<"events created for all nodes and added to PEPStack - size = "<<PEPStack.size()<<endl;
    cout<<"TwoPhaseDESTransport constructed"<<endl;
} // end constructor  



template<size_t dim>
void TwoPhaseDESTransport<dim>::initializeVariablsAndKeys(Model<dim>& m)
{
    //creating new variables if not defined yet from input file
    if(!m.Database().IsDefined("variation rate nonwetting phase")) m.CreateProperty( "variation rate nonwetting phase", "m3/(m3.s)", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);
    if(!m.Database().IsDefined("nodal fluid volume source")) m.CreateProperty( "nodal fluid volume source", "m3", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);    
    if(!m.Database().IsDefined("event index")) m.CreateProperty( "event index", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("update count")) m.CreateProperty( "update count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("rate count")) m.CreateProperty( "rate count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("schedule count")) m.CreateProperty( "schedule count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("synchronize count")) m.CreateProperty( "synchronize count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("nonwetting phase timing array")) m.CreateProperty( "nonwetting phase timing array", "none", ARRAY, NODE, 6, -1.00E+10 ,1.00E+10); 
    
    // model-wide initialisation
    m.Region("Model").InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );    
    
    //assigning keys 
    key_dsnw = INDEX<SCALAR,NODE> ( m.Database().StorageKey("variation rate nonwetting phase") );
    key_NQV = INDEX<SCALAR,NODE> ( m.Database().StorageKey("nodal fluid volume source") );
    key_EventIndex = INDEX<SCALAR,NODE>( m.Database().StorageKey("event index") );
    key_update = INDEX<SCALAR,NODE>( m.Database().StorageKey("update count") );
    key_rate = INDEX<SCALAR,NODE>( m.Database().StorageKey("rate count") );
    key_schedule = INDEX<SCALAR,NODE>( m.Database().StorageKey("schedule count") );
    key_synchronize = INDEX<SCALAR,NODE>( m.Database().StorageKey("synchronize count") );
    key_time = INDEX<ARRAY,NODE>( m.Database().StorageKey("nonwetting phase timing array") );
    
    //checking keys 
    if ( key_dsnw.place != NODE || key_dsnw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'variation rate nonwetting phase' variable must be SCALAR and placed on NODE"  );    
    if ( key_NQV.place != NODE || key_NQV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );                                    
    if ( key_EventIndex.place != NODE || key_EventIndex.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'event index' variable must be SCALAR and placed on NODE"  );       
    if ( key_update.place != NODE || key_update.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( key_rate.place != NODE || key_rate.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( key_schedule.place != NODE || key_schedule.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( key_synchronize.place != NODE || key_synchronize.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  ); 
    if ( key_time.place != NODE || key_time.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'nonwetting phase timing array' variable must be ARRAY and placed on NODE"  );                
}


template<size_t dim>
void TwoPhaseDESTransport<dim>::initializeFiniteVolumeProperties(Node<dim>* nd)
{
    double64 pore_volume  = 0.;
    for (auto sip : nd->AllSectorIntegrationPoints()) {
        const double64 sector_volume = sip.SectorVolume();
        double64 phi = sip.Interpolate( this->key_phi );
        double64 thickness = sip.Interpolate( this->key_thi );
        if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
        pore_volume += sector_volume * phi;
    }
    nd->Store( this->key_fvPV, makeScalar(PLAIN, pore_volume) );
          
} // end initializeFiniteVolumeProperties



//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim>
void TwoPhaseDESTransport<dim>::ComputeRateofChange( Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );

    //ignore DIRICH node
    if(nd->Status(  this->key_sCO2 ) == DIRICH) {
        nd->Store(  this->key_fb, makeScalar( nd->Status(  this->key_fb ), 0. ) );//flux balance
        ArrayVariable array;
        nd->Read(key_time, array);
        array.Component(2, numeric_limits<double64>::max()); //CFL time increment
        nd->Store(key_time, array);    
        return;
    };

    rate_count_++;//recording
    nd->Store(  key_rate, makeScalar( nd->Status( key_rate), nd->Read( key_rate) + 1 ) );

    double64 accumulation(0.), flux_balance(0.), outflow(0.);
    
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD;
    
    for (auto fip : nd->AllFacetIntegrationPoints()) {
        const double64 sign = fip.FromInside() ? 1. : -1.;
        fip.Interpolate( this->key_vt, vD );
        const double64 vD_nA = vD.DotProduct(fip.DirectedArea());
        const double64 facet_flux = sign * vD_nA;
        flux_balance += facet_flux;       
        if ( facet_flux > 0. ) outflow += facet_flux ;    
              
        auto inside_node = fip.InsideNode();
        auto outside_node = fip.OutsideNode();  
        double64 vD_n = fip.ProjectOntoFacetNormal(vD);  
        double64 facetArea = fip.FacetArea();  
        Point<dim> facetNrml = fip.FacetNormal();   
        double64 viscous_velocity_component(0.0), capillary_velocity_component(0.0), gravity_velocity_component(0.0);
        
        if ( vD_n != 0. ) {
            const double64 sn_inside_node  = inside_node.Read( this->key_sCO2 );
            const double64 sn_outside_node = outside_node.Read( this->key_sCO2 );
                          
            const double64  psi_hat_c( (vD_n < 0.) ? sn_outside_node : sn_inside_node );         
            double64 f_n = flowfunctions_.f(fip, 1U, 1.-psi_hat_c);                         
            viscous_velocity_component = vD_n* f_n;  
        }
            
        if(with_gravity_forces_)
            gravity_velocity_component = flowfunctions_.GravityMultiplier_G(fip) * facetNrml[v];
        
        if(with_capillary_spreading_) {
            Point<dim> snw_gradient = fip.Gradient (this->key_sCO2);
            double64 dsdn = dotProduct(snw_gradient, facetNrml);
            capillary_velocity_component = flowfunctions_.CapillaryDiffusionMultiplier(fip) * (-dsdn);
        }          
                                                   
        accumulation += sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;                                  
    } 
    
    //compute CFL time increment 
    ArrayVariable array2;
    nd->Read(key_time, array2);
    if (outflow < numeric_limits<double64>::epsilon())
        array2.Component(2, numeric_limits<double64>::max());
    else 
        array2.Component(2, nd->Read(  this->key_fvPV ) / outflow);
    nd->Store(key_time, array2);
    
     
    // divergence free correction
    // compute average fractional flow for the current finite volume
    double64 fn_avg = 0.;
    if (fabs(flux_balance) > numeric_limits<double64>::epsilon()) {
        for (auto sip : nd->AllSectorIntegrationPoints()) {//??????? not sure if it should be sector based
            fn_avg += flowfunctions_.f(sip, 1U);
        }
        fn_avg /= nd->Parents();  
        accumulation -= fn_avg*flux_balance;     
    }     
    
    double64 PV = nd->Read(  this->key_fvPV );
    //store variation rate
    nd->Store(  key_dsnw, makeScalar( nd->Status( key_dsnw ), accumulation/PV ) );         
        
}        


//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim>
void TwoPhaseDESTransport<dim>::ComputeRateofChange_upstream( Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );

    //ignore DIRICH node
    if(nd->Status(  this->key_sCO2 ) == DIRICH) {
        nd->Store(  this->key_fb, makeScalar( nd->Status(  this->key_fb ), 0. ) );//flux balance
        ArrayVariable array;
        nd->Read(key_time, array);
        array.Component(2, numeric_limits<double64>::max()); //CFL time increment
        nd->Store(key_time, array);    
        return;
    };

    rate_count_++;//recording
    nd->Store(  key_rate, makeScalar( nd->Status( key_rate), nd->Read( key_rate) + 1 ) );

    double64 accumulation(0.), flux_balance(0.), outflow(0.);
    
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD;
    
    for (auto fip : nd->AllFacetIntegrationPoints()) {
        const double64 sign = fip.FromInside() ? 1. : -1.;
        fip.Interpolate( this->key_vt, vD );
        const double64 vD_nA = vD.DotProduct(fip.DirectedArea());
        const double64 facet_flux = sign * vD_nA;
        flux_balance += facet_flux;       
        if ( facet_flux > 0. ) outflow += facet_flux ;    
              
        auto inside_node = fip.InsideNode();
        auto outside_node = fip.OutsideNode();  
        double64 vD_n = fip.ProjectOntoFacetNormal(vD);  
        double64 facetArea = fip.FacetArea();  
        Point<dim> facetNrml = fip.FacetNormal();   
        double64 viscous_velocity_component(0.0), capillary_velocity_component(0.0), gravity_velocity_component(0.0);
        if( !with_gravity_forces_ && !with_capillary_spreading_ ){ //viscous effect only                                           
            if ( vD_n != 0. ) {
                const double64 sn_inside_node  = inside_node.Read( this->key_sCO2 );
                const double64 sn_outside_node = outside_node.Read( this->key_sCO2 );
                          
                const double64  psi_hat_c( (vD_n < 0.) ? sn_outside_node : sn_inside_node );         
                double64 f_n = flowfunctions_.f(fip, 1U, 1.-psi_hat_c);                         
                viscous_velocity_component = vD_n* f_n;  
            }
            
        } else { // with gravity or capillary effects            
            double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
            if( with_gravity_forces_ ){                       
                vn_gravity_component_of_velocity = flowfunctions_.Mobility(fip, 0U) * flowfunctions_.GravityTerm(fip) * facetNrml[v];
                vw_gravity_component_of_velocity = flowfunctions_.Mobility(fip, 1U) * flowfunctions_.GravityTerm(fip) * facetNrml[v];
            }         
            
            if(with_capillary_spreading_){                
                Point<dim> snw_gradient = fip.Gradient (this->key_sCO2);
                                             
                double64 dsdn = dotProduct(snw_gradient, facetNrml);
                
                vn_capillary_component_of_velocity = -dsdn*flowfunctions_.CapillaryDiffusionMultiplier(fip)*flowfunctions_.TotalMobility(fip)/flowfunctions_.Mobility(fip, 1U);
                vw_capillary_component_of_velocity = -dsdn*flowfunctions_.CapillaryDiffusionMultiplier(fip)*flowfunctions_.TotalMobility(fip)/flowfunctions_.Mobility(fip, 0U);              
            }  
            
            double64 vn_at_facet_int_point = vD_n - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            double64 vw_at_facet_int_point = vD_n + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity;  
                    
            // mixture moving inside the CV        mixture moving outside the CV
            //
            //       |  vt                                 /|\ vt
            //       |                                      |
            //      \|/              /|\                    |
            //     -----              |  N_up             -----
            //   /       \                              /       \
            //   \       /                              \       /
            //     -----              |                   -----
            //      /|\              \|/  N_down            |
            //       |                                      |
            //       |  vt                                 \|/ vt

            double64 upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);

            //computes mobilities on inside and outside nodes
            const double64 ln_inside_node  = flowfunctions_.Mobility(inside_node, 1U);
            const double64 lw_inside_node  = flowfunctions_.Mobility(inside_node, 0U);
                    
            const double64 ln_outside_node  = flowfunctions_.Mobility(outside_node, 1U);
            const double64 lw_outside_node  = flowfunctions_.Mobility(outside_node, 0U);     
            
            double64 zero(0.0);
            if((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point>zero)){
                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_inside_node;
            }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero)){
                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_outside_node;
            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point>zero)){
                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_inside_node;
            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point<zero)){
                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_outside_node;
            }else{
                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);
            }                     
                    
            total_mobility=upstream_mobility_n+upstream_mobility_w;
            double64 upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            double64 upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);
                    
            viscous_velocity_component = vD_n * upstream_fn;
                                       
            if( with_gravity_forces_ )
                gravity_velocity_component = upstream_lambda_overbar * flowfunctions_.GravityTerm(fip) * facetNrml[v];

            if( with_capillary_spreading_ )
                capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity; 
        }               
                                                   
        accumulation += sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;                                  
    } 
    
    //compute CFL time increment 
    ArrayVariable array2;
    nd->Read(key_time, array2);
    if (outflow < numeric_limits<double64>::epsilon())
        array2.Component(2, numeric_limits<double64>::max());
    else 
        array2.Component(2, nd->Read(  this->key_fvPV ) / outflow);
    nd->Store(key_time, array2);
    
     
    // divergence free correction
    // compute average fractional flow for the current finite volume
    double64 fn_avg = 0.;
    if (fabs(flux_balance) > numeric_limits<double64>::epsilon()) {
        for (auto sip : nd->AllSectorIntegrationPoints()) {//??????? not sure if it should be sector based
            fn_avg += flowfunctions_.f(sip, 1U);
        }
        fn_avg /= nd->Parents();  
        accumulation -= fn_avg*flux_balance;     
    }     
    
    double64 PV = nd->Read(  this->key_fvPV );
    //store variation rate
    nd->Store(  key_dsnw, makeScalar( nd->Status( key_dsnw ), accumulation/PV ) );         
        
}  



//schedule an event associated with a node/FV
template<size_t dim>
bool TwoPhaseDESTransport<dim>::Schedule(Event<dim>* event, double64 t_end, double64 cfl_multiplier)
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );
    ArrayVariable array;
    nd->Read(key_time, array);
    //ignore DIRICH node
    if(nd->Status(  this->key_sCO2 ) == DIRICH) {
        return false;
    } else {
        nd->Store( key_schedule, makeScalar( nd->Status( key_schedule), nd->Read( key_schedule) + 1 ) );
        event->valid(true);
        //compute target change
        double64 CFL = array[2];//CFL number
        double64 ChangeRate = nd->Read( key_dsnw);//rate of change
  
        double64 dC_CFL = -CFL*cfl_multiplier*ChangeRate;//targe change

        if (fabs(dC_CFL) < numeric_limits<double64>::epsilon()){//idle node/FV
            array.Component(5, numeric_limits<double64>::epsilon());//target change of solution
            array.Component(3, numeric_limits<double64>::max());//target time increment          
        } else {
            array.Component(5, dC_CFL);//target change of solution
            array.Component(3, cfl_multiplier*CFL);//target time increment          
        };

        double64 t_current = array[0];//current time stamp
        double64 dt_target = array[3];//target time increment
        if ((dt_target + t_current) >= t_end) {
            nd->Store(key_time, array);
            return false;
        } else {
            event->t_schedule(t_current + dt_target);
            array.Component(1, t_current + dt_target);//schedule time stamp
            nd->Store(key_time,array);
            return true;
        };
    };
}



//update solution and check it against the specified range (with DES)
template<size_t dim>
void TwoPhaseDESTransport<dim>::Update_DES(Event<dim>* event, double64 t_clock)
{
    update_count_++;//recording
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );
    ArrayVariable array;
    nd->Read(key_time, array);
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    if ( status != DIRICH )
    {    
        double64 ChangeRate = nd->Read( key_dsnw);//variaition rate   
        double64 solution = nd->Read( this->key_sCO2);//old solution
        double64 t_current = array[0]; //current time stamp
        double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
        const double64 source(nd->Read(this->key_NQV));
        new_solution += source * (t_clock - t_current);//add source to new solution
        
        //check new solution value against range and stored it to key_sCO2
        if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        else {
            cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( new_solution > upper_limit_ ) nd->Store( this->key_sCO2, makeScalar( status, upper_limit_ ) );
            else if ( new_solution < lower_limit_ ) nd->Store( this->key_sCO2, makeScalar( status, lower_limit_ ) );
        }
        
        new_solution = nd->Read(this->key_sCO2);//stored new solution
        double64 dsn_cumulative = array[4];
        array.Component(4, dsn_cumulative + (new_solution-solution));//update cumulative change
        
        array.Component(0, t_clock); //current time stamp
        nd->Store(key_time, array);

        nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) );
    }
}



//update solution and check it against the specified range (with TDS)
template<size_t dim>
void TwoPhaseDESTransport<dim>::Update_TDS(Event<dim>* event, double64 delta_t)
{
    update_count_++;//recording
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );    
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    if ( status != DIRICH )
    {    
        double64 ChangeRate = nd->Read(key_dsnw);//variation rate  
        double64 solution = nd->Read(this->key_sCO2);//old solution
        double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
        const double64 source(nd->Read( this->key_NQV));
        new_solution += source * delta_t;//add source to new solution.
                
        //check new solution value against range and stored it to key_sCO2
        if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        else {
            cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( new_solution > upper_limit_ ) nd->Store( this->key_sCO2, makeScalar( status, upper_limit_ ) );
            else if ( new_solution < lower_limit_ ) nd->Store( this->key_sCO2, makeScalar( status, lower_limit_ ) );
        }
        
        nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) );   
    }
}



//Synchronize neighbor nodes/FVs
template<size_t dim>
void TwoPhaseDESTransport<dim>::Synchronize(Event<dim>* event,double64 t_clock)
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL ); 
    nd->Store( key_synchronize, makeScalar( nd->Status(key_synchronize), nd->Read(key_synchronize) + 1 ) );
    event->valid(false);
    ArrayVariable array;
    nd->Read(key_time, array);
    array.Component(4, 0.); //reset cumulative change of solution
    nd->Store(key_time, array);
    for ( size_t n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* neighbor_node = nd->Neighbor(n);
        int index = neighbor_node->Read(key_EventIndex);
        Event<dim>* neighbor_event = EntireQueue[index];  
        assert( neighbor_event  != NULL ); 
        if (neighbor_event->inPEPStack() == false) {
            PEPStack.push_back(neighbor_event);
            neighbor_event->inPEPStack(true);
            Update_DES(neighbor_event,t_clock);
            ArrayVariable neighbor_array;
            neighbor_node->Read(key_time, neighbor_array); 
            double64 dC_cumulative = neighbor_array[4];//cumulative change of solution
            double64 dC_target = neighbor_array[5];//target change of solution
            if (fabs(dC_cumulative) >= fabs(dC_target)) Synchronize (neighbor_event, t_clock); 
        };
    };
}


//advect variable with TDS (time-driven simulation)
template<size_t dim>
void TwoPhaseDESTransport<dim>::AdvectVariable_TDS( double64 time_interval, double64 cfl_multiplication_factor, double64 PEP_parameter )
{
    double64 begin=clock();
    double64 time_increment(time_interval); 
    
    clock_t T_begin= clock();
    const typename vector<Event<dim>*>::iterator stack_end(PEPStack.end());
    for ( typename vector<Event<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )   
    {    
        ComputeRateofChange((*it));  
        ArrayVariable array;
        (*it)->getNode()->Read(key_time, array);
        double64 dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*cfl_multiplication_factor);
    }
    T_RateOfChange_ += clock() - T_begin; 

    const double64 one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double64 time(0.);
    
    while (time < time_interval)
    {
        cout <<"\n\tadvection (sub)step: "<< substep << endl;
        T_begin= clock();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        for ( typename vector<Event<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };
        T_Update_ += clock() - T_begin;
        
        T_begin= clock();
        double64 new_time_increment(time_interval); 
        for ( typename vector<Event<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )
        { 
            ComputeRateofChange((*it));
            ArrayVariable array2;
            (*it)->getNode()->Read(key_time, array2);
            double64 dt_CFL = array2[2];//CFL time increment
            new_time_increment=min(new_time_increment, dt_CFL*cfl_multiplication_factor);
        };
        T_RateOfChange_ += clock() - T_begin; 

        time += time_increment;
        substep++;
    };
    
    T_AdvectVariable_+= clock() - begin;

    cout <<"Finish DESTransport<dim>::AdvectVariable_TDS "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl;  
    cout <<"T_Schedule_ = "<< T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_SortQueue_ = "<< T_SortQueue_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RemoveFromQueue_ = "<< T_RemoveFromQueue_/double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;  
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl; 
}



//advect variable with DES (discrete event simulation)
template<size_t dim>
void TwoPhaseDESTransport<dim>::AdvectVariable_DES_serial( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter )
{
    double64 begin=clock();
    cout<<"Start DESTransport<dim>::AdvectVariable_DES_serial "<<endl;
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
    
    while (!Finished)
    {    
        clock_t T_begin= clock();
        const typename vector<Event<dim>*>::iterator stack_end(PEPStack.end());
        for ( typename vector<Event<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )       
        {   
            ComputeRateofChange((*it));            
            if ((*it)->valid() == false)
                if (Schedule((*it), model_time, cfl_multiplication_factor))
                    Queue.push_back((*it));    
           (*it)->inPEPStack(false);          
        };    
        T_RateOfChange_ += clock() - T_begin; 
        cout <<" Queue size = "<< Queue.size()<<endl;   

        T_begin= clock();        
        if (Queue.empty()) {
            time=model_time;
        } else {
            sort(Queue.begin(),Queue.end(),sort_queue<dim>());                       
            const typename vector<Event<dim>*>::iterator begin(Queue.begin());
            ArrayVariable begin_array;
            (*begin)->getNode()->Read(key_time, begin_array);
            time = begin_array[1]; //scheduled time stamp
        };
        T_SortQueue_ += clock() - T_begin;
        cout<<"  time = "<<time<<" model_time = "<<model_time<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Event<dim>*>::iterator End(PEPStack.end());
            for ( typename vector<Event<dim>*>::iterator e=PEPStack.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        PEPStack.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
        size_t count = 0U;
        while (!Queue.empty())
        {           
            const typename vector<Event<dim>*>::iterator top(Queue.begin()); 
            if((*top)->valid() == false) {
                Queue.erase(top);
                continue;
            }
            count++;            
            ArrayVariable array;
            (*top)->getNode()->Read(key_time, array);
            double64 dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, PEP_parameter*dt_target);
            double64 t_schedule = array[1];//scheduled time stamp
            if (t_schedule > (time+dt_PEP)) break;            
            if ((*top)->inPEPStack() == false) {
                PEPStack.push_back((*top));
                (*top)->inPEPStack(true);
                T_begin= clock();
                Update_DES((*top),time);
                T_Update_ += clock() - T_begin;                
            }; 
            T_begin= clock();
            Synchronize((*top),time);
            T_Synchronize_ += clock() - T_begin;
            Queue.erase(top); 
        };
        //remove invalid events/nodes from Queue
        T_begin= clock();
        if(!Queue.empty()) {
            size_t queue_size = Queue.size();
            for ( size_t n = 0U; n<queue_size; n++ )
            {  
                if (Queue[n]->valid() == false) {
                    Queue.erase(Queue.begin()+n);
                    queue_size --;
                    n--;
                };
            };
        };
        T_RemoveFromQueue_ += clock() - T_begin;        
        cout<<"  count =  "<<count<<endl;
    };
    T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DESTransport<dim>::AdvectVariable_DES_serial "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_SortQueue_ = "<< T_SortQueue_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RemoveFromQueue_ = "<< T_RemoveFromQueue_/double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;  
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl;
}   


template class TwoPhaseDESTransport<1U>;
template class TwoPhaseDESTransport<2U>;
template class TwoPhaseDESTransport<3U>;

} // end csmp 



                                                                       
