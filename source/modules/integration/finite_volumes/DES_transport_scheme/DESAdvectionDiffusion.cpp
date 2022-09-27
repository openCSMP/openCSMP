#include "DESAdvectionDiffusion.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#if defined(_OPENMP )
#include "omp.h"
#endif
#include "Exception.h"
#include "ErrorHandler.h"


using namespace std;

namespace csmp {

template<uint32_t dim>
DESAdvectionDiffusion<dim>::DESAdvectionDiffusion( Model<dim>& m, const char* target_region, double cfl_multiplier, double PEP_multiplier, bool tensor_k )
  : sg_(m), 
    gref_(m.Region(target_region)),
    db_(m.Database()),
    CFL_multiplier_(cfl_multiplier), PEP_multiplier_(PEP_multiplier),
    tensor_k_(tensor_k),
    rate_count_(0U), update_count_(0U),
    T_RateOfChange_(0.), T_Schedule_(0.), T_InsertToHeap_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromHeap_(0.), T_AdvectVariable_(0.)
{
    InitializeVariablsAndKeys();
    InitializeFiniteVolumeProperties();
    InitializeEvents();
    cout<<"\nDESAdvectionDiffusion constructed"<<endl;
    cout<<"CFL multiplier = "<<CFL_multiplier_<<endl;
    cout<<"PEP multiplier = "<<PEP_multiplier_<<endl;
    cout<<"tensor permeability = "<<tensor_k<<" (0=false, 1=true)\n"<<endl;
}
  


template<uint32_t dim>
void DESAdvectionDiffusion<dim>::InitializeVariablsAndKeys()
{
    //creating new variables if not defined yet 
    if(!db_.IsDefined("event index")) sg_.CreateProperty( "event index", "EI", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    if(!db_.IsDefined("update count")) sg_.CreateProperty( "update count", "UC", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    if(!db_.IsDefined("rate count")) sg_.CreateProperty( "rate count", "RC", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    if(!db_.IsDefined("schedule count")) sg_.CreateProperty( "schedule count", "SDC", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    if(!db_.IsDefined("synchronize count")) sg_.CreateProperty( "synchronize count", "SC", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    if(!db_.IsDefined("DES array")) sg_.CreateProperty( "DES array", "DESa", "none", ARRAY, NODE, 6, -1.00E+10 ,1.00E+10);
    if(!db_.IsDefined("porosity")) sg_.CreateProperty( "porosity", "phi", "none", SCALAR, ELEMENT, 1, 1.00E-05, 1.00E+01);
    if(!db_.IsDefined("thickness")) sg_.CreateProperty( "thickness", "thi", "m", SCALAR, ELEMENT, 1, 0.0E+0, 1.00E+10);
    if(!db_.IsDefined("facet area")) sg_.CreateProperty( "facet area", "fA", "m2", SCALAR, FACET_INTEGRATION_POINT, 1, -1.00E+10, 1.00E+10);
    if(!db_.IsDefined("facet normal")) sg_.CreateProperty( "facet normal", "fN", "m2 s-1", VECTOR, FACET_INTEGRATION_POINT, 3, 0., 1.);
    if(!db_.IsDefined("concentration")) sg_.CreateProperty( "concentration", "C", "kg m-3", SCALAR, NODE, 1, -5.00E-01, 1.00E+03);
    if(!db_.IsDefined("new concentration")) sg_.CreateProperty( "new concentration", "Cn", "kg m-3", SCALAR, NODE, 1, -5.00E-01, 1.00E+03);
    if(!db_.IsDefined("flux balance")) sg_.CreateProperty( "flux balance", "fb", "m3 s-1", SCALAR, NODE, 1, -1.00E+10, 1.00E+10);
    if(!db_.IsDefined("nodal concentration source")) sg_.CreateProperty( "nodal concentration source", "nCq", "kg/m3 s", SCALAR, NODE, 1, -1.00E+04, 1.00E+04);
    if(!db_.IsDefined("velocity")) sg_.CreateProperty( "velocity", "v", "m s-1", VECTOR, ELEMENT, 3, -1.00E+02, 1.00E+02);
    if(!db_.IsDefined("FV pore volume")) sg_.CreateProperty( "FV pore volume", "fvpV", "m3", SCALAR, NODE, 1, 0.00E+00 ,1.00E+8);
    
    //assigning keys  
    key_EventIndex = INDEX<SCALAR,NODE>( db_.StorageKey("event index") );
    key_update = INDEX<SCALAR,NODE>( db_.StorageKey("update count") );
    key_rate = INDEX<SCALAR,NODE>( db_.StorageKey("rate count") );
    key_schedule = INDEX<SCALAR,NODE>( db_.StorageKey("schedule count") );
    key_synchronize = INDEX<SCALAR,NODE>( db_.StorageKey("synchronize count") );
    key_time = INDEX<ARRAY,NODE>( db_.StorageKey("DES array") ); 
    key_phi = INDEX<SCALAR,ELEMENT>( db_.StorageKey("porosity") );  
    key_thi = INDEX<SCALAR,ELEMENT>( db_.StorageKey("thickness") ); 
    key_fA = INDEX<SCALAR,FACET_INTEGRATION_POINT>( db_.StorageKey("facet area") );
    key_fn = INDEX<VECTOR,FACET_INTEGRATION_POINT>( db_.StorageKey("facet normal") );    
    key_C0 = INDEX<SCALAR,NODE>( db_.StorageKey("concentration") );
    key_FB = INDEX<SCALAR,NODE>( db_.StorageKey("flux balance") );
    key_C1 = INDEX<SCALAR,NODE>( db_.StorageKey("new concentration") );
    key_NQC = INDEX<SCALAR,NODE>( db_.StorageKey("nodal concentration source") );
    key_V = INDEX<VECTOR,ELEMENT>( db_.StorageKey("velocity") );
    key_fvPV = INDEX<SCALAR,NODE>( db_.StorageKey("FV pore volume") );
    
    //checking keys
    if ( key_EventIndex.place != NODE || key_EventIndex.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'event index' variable must be SCALAR and placed on NODE"  );       
    if ( key_update.place != NODE || key_update.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( key_rate.place != NODE || key_rate.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( key_schedule.place != NODE || key_schedule.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( key_synchronize.place != NODE || key_synchronize.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  );    
    if ( key_time.place != NODE || key_time.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'DES array' variable must be ARRAY and placed on NODE"  );   
    if ( key_phi.place != ELEMENT || key_phi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );    
    if ( key_thi.place != ELEMENT || key_thi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  ); 
    if ( key_fA.place != FACET_INTEGRATION_POINT || key_fA.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );       
    if ( key_fn.place != FACET_INTEGRATION_POINT || key_fn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );       
    if ( key_C0.place != NODE || key_C0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'concentration' variable must be SCALAR and placed on NODE"  );
    if ( key_FB.place != NODE || key_FB.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );   
    if ( key_C1.place != NODE || key_C1.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'new concentration' variable must be SCALAR and placed on NODE"  );     
    if ( key_NQC.place != NODE || key_NQC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'nodal concentration source' variable must be SCALAR and placed on NODE"  ); 
    if ( key_V.place != ELEMENT || key_V.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DESAdvectionDiffusion::InitializeVariablsAndKeys:",
        "The 'velocity' variable must be VECTOR and placed on ELEMENT"  ); 
    if ( key_fvPV.place != NODE || key_fvPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );                               
        
    // model-wide initialisation
    sg_.Region("Model").InputPropertyValue( "update count", makeScalar(PLAIN,0.), COMPLETE );
    sg_.Region("Model").InputPropertyValue( "rate count", makeScalar(PLAIN,0.), COMPLETE );   
    sg_.Region("Model").InputPropertyValue( "schedule count", makeScalar(PLAIN,0.), COMPLETE ); 
    sg_.Region("Model").InputPropertyValue( "synchronize count", makeScalar(PLAIN,0.), COMPLETE ); 
    sg_.Region("Model").InputPropertyValue( "DES array", ArrayVariable(6,0.,PLAIN), COMPLETE);     
    sg_.Region("Model").InputPropertyValue( "new concentration", makeScalar(PLAIN,0.), COMPLETE );
    sg_.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );    
    db_.RangeOf( db_.Name(key_C0), lower_limit_, upper_limit_ );
              
}



template<uint32_t dim>
void DESAdvectionDiffusion<dim>::InitializeFiniteVolumeProperties()
 {
    //zeroing FV pore volumes for accumulation in element loop
    gref_.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE ); 
    
    // For the interior elements of the region compute relevant variable values
    const auto it_end(gref_.CellsEnd());
    for ( auto it=gref_.CellsBegin(); it!=it_end; ++it )
    {
         const auto sectors((*it)->Sectors());
         const auto facets((*it)->Facets());

         // computing sector pore volumes
         double phi = (*it)->Read( key_phi);
         const double thickness = (*it)->Read( key_thi );
         if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
         
         for ( auto i{0U}; i<sectors; ++i ) {
              const double sector_volume = (*it)->SectorVolume(i);  
              double pore_volume   = (*it)->N(i)->Read( key_fvPV );
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( key_fvPV, makeScalar(PLAIN,pore_volume) );
         }

         // computing facet normals and areas
         for ( auto j{0U}; j<facets; ++j ) {
              const double facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, key_fA, makeScalar( PLAIN, facet_area ) );
              Point<dim> nrml = (*it)->FacetNormal(j);
              VectorVariable<dim>  fnrml;
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, key_fn, fnrml );
         }            
   }

   // initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // (here the pore volumes do not include the sectors outside the region)
   const auto nit_end(gref_.NodesEnd());
   for ( auto nit=gref_.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const auto parent_elements((*nit)->Parents());
        for ( auto i{0U}; i<parent_elements; ++i ) {
             Element<dim>* const eptr = (*nit)->Parent(i);
             // computing facet normals and areas
             const auto facets(eptr->Facets());
             for ( auto j{0U}; j<facets; ++j ) {
                  const double facet_area = eptr->FacetArea(j);
                  eptr->Store( j, 0U, key_fA, makeScalar( PLAIN, facet_area ) );
                  Point<dim> nrml = eptr->FacetNormal(j);
                  VectorVariable<dim>  fnrml;
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, key_fn, fnrml );
             }             
        }
   }
 } // end initializeFiniteVolumeProperties




//advect variable with TDS (time driven simulation)
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::AdvectVariable_TDS( double time_interval, size_t num_threads )
{
#if defined(_OPENMP )
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
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::AdvectVariable_DES( double model_time, size_t num_threads )
{
#if defined(_OPENMP )
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





template<uint32_t dim>
void DESAdvectionDiffusion<dim>::InitializeEvents()
{        
    //create events for all nodes and add them to event lists
    size_t index = 0;
    size_t dirich_count = 0;
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        if((*nit)->Status(  key_C0 ) != DIRICH) {
            (*nit)->Store( key_EventIndex, makeScalar( (*nit)->Status(key_EventIndex), index) );//event index 
            Event<dim>* event = new Event<dim>(*nit);
            PEPList.push_back(event);
            event->inPEPStack(true);
            ComputeFluxBalanceAndCFL(event);
            Schedule(event, 0.);
            event->valid(false);
            Heap_Node* heap_node = new Heap_Node(event->t_schedule(),index);
            HeapNodeFullList.push_back(heap_node);
            FullList.push_back(event);
            event->inQueue(false);
                
            index++;
        
        } else {
            dirich_count++;
        }
    }
    cout<<FullList.size()<<" events created for all nodes, excluding "<<dirich_count<<" DIRICH nodes"<<endl;
}  




//Compute the flux balance and CFL time increment for a node/FV
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::ComputeFluxBalanceAndCFL( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status( key_C0 ) != DIRICH);

  if(nd  != NULL && nd->Status( key_C0 ) != DIRICH) {
    double flux_balance(0.), outflow(0.);
    
    VectorVariable<dim> vD, facetNrml;
    const auto node_parent_elements(nd->Parents());
    
    for ( auto t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        const auto pnid(nd->ParentNodeNumber(t));
        eptr->Read( key_V, vD);

        const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
        for ( auto i{0U}; i<sector_facets; i++ )
        {
            const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
            const auto inside_node(eptr->FV()->InsideNode(iFacet));
                        
            eptr->Read( iFacet, 0U,  key_fn, facetNrml );
            const double  vD_n = vD.DotProduct(facetNrml);
            const double  facetArea = eptr->Read( iFacet, 0U, key_fA ); 
            
            const double sign = ( pnid == inside_node ) ? 1. : -1.;
            //compute facet fluid flux
            double facet_flux = sign * vD_n * facetArea;
            //update flux balance
            flux_balance += facet_flux;
            //update outflow
            if ( facet_flux > 0. ) outflow += facet_flux;
        }
    }
        
    nd->Store(  key_FB, makeScalar( nd->Status(  key_FB ), flux_balance ) );//flux balance
    
    //CFL time increment 
    ArrayVariable array;
    nd->Read(key_time, array);
    if (outflow < numeric_limits<double>::epsilon())
        array.Component(2, numeric_limits<double>::max());
    else 
        array.Component(2, nd->Read( key_fvPV ) / outflow);  
    nd->Store(key_time, array);
  }
}



//Compute the rate of change in a node/FV
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status(  key_C0 ) != DIRICH) {
    rate_count_++;//recording
    nd->Store( key_rate, makeScalar( nd->Status(key_rate), nd->Read(key_rate) + 1 ) );

    double accumulation(0.);
    VectorVariable<dim> vD, facetNrml;
    const auto node_parent_elements(nd->Parents());
    
    for ( auto t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        const auto pnid(nd->ParentNodeNumber(t));
        eptr->Read( key_V, vD);

        const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
        for ( auto i{0U}; i<sector_facets; i++ )
        {
            const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
            const auto inside_node(eptr->FV()->InsideNode(iFacet));
            const auto outside_node(eptr->FV()->OutsideNode(iFacet));
            
            eptr->Read( iFacet, 0U, key_fn, facetNrml );
            const double  vD_n = vD.DotProduct(facetNrml);
            const double  facetArea = eptr->Read( iFacet, 0U, key_fA ); 
            
            const double sign = ( pnid == inside_node ) ? 1. : -1.;
            //compute facet fluid flux
            double facet_flux = sign * vD_n * facetArea;
            // finding the upstream node
            const auto upstream_node = (vD_n < 0.) ? outside_node : inside_node;
            const double C_upstream = eptr->N(upstream_node)->Read( key_C0 );
            accumulation += facet_flux * C_upstream;  
        }
    }    

    nd->Store(  key_C1, makeScalar( nd->Status( key_C1 ), accumulation ) );    
  }
} 




//schedule an event associated with a node/FV
template<uint32_t dim>
bool DESAdvectionDiffusion<dim>::Schedule(Event<dim>* event, double t_end)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status( key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( key_C0 ) != DIRICH) {    
    ArrayVariable array;
    nd->Read(key_time, array);

    nd->Store( key_schedule, makeScalar( nd->Status(key_schedule), nd->Read(key_schedule) + 1 ) );
    event->valid(true);
    //compute target change
    double CFL = array[2];//CFL number
    double PV = nd->Read(key_fvPV);//Pore volume
    double ChangeRate = nd->Read( key_C1);//rate of change
    double C0 = nd->Read( key_C0);//concentration
    double flux_balance = nd->Read( key_FB);//flux balance

    double dC_CFL = -CFL*CFL_multiplier_/PV*(ChangeRate-C0*flux_balance);//targe change

    if (fabs(dC_CFL) < numeric_limits<double>::epsilon()){//idle node/FV
        array.Component(5, numeric_limits<double>::epsilon());//target change of solution
        array.Component(3, numeric_limits<double>::max());//target time increment          
    } else {
        array.Component(5, dC_CFL);//target change of solution
        array.Component(3, CFL_multiplier_*CFL);//target time increment      
    };

    double t_current = array[0];//current time stamp
    double dt_target = array[3];//target time increment
    if ((dt_target + t_current) >= t_end) {
        event->t_schedule(t_end);
        array.Component(1, t_end);//schedule time stamp
        nd->Store(key_time, array);
        return false;
    } else {
        event->t_schedule(t_current + dt_target);
        array.Component(1, t_current + dt_target);//schedule time stamp
        nd->Store(key_time,array);
        return true;
    }
  }
  return false;
}


//update solution and check it against the specified range (with DES)
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::Update_DES(Event<dim>* event, double t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status( key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( key_C0 ) != DIRICH) {  
    update_count_++;//recording  
    const VARIABLE_FLAG status(nd->Status( key_C0 ));
    ArrayVariable array;
    nd->Read(key_time, array);

    // 1. starting with the sum of facet flux-concentration products stored in C1
    double ChangeRate = nd->Read( key_C1); 
    double solution = nd->Read( key_C0);//concentration
    // 2. correcting this sum for div vD using 'flux balance'   
    ChangeRate -= solution * nd->Read( key_FB ); 
    // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
    double t_current = array[0]; //current time stamp
    double new_solution = solution - ((t_clock - t_current)/nd->Read(key_fvPV)) * ChangeRate;
    // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
    // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
    const double source(nd->Read(key_NQC)); //kg/m3/s
    // new concentration
    new_solution += source * (t_clock - t_current);

    double C_last = nd->Read(key_C0);//last concentration
        
    if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(key_C0, makeScalar( status, new_solution ));//store solution value to key_C0
    else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
        if ( new_solution > upper_limit_ ) nd->Store( key_C0, makeScalar( status, upper_limit_ ) );
        else if ( new_solution < lower_limit_ ) nd->Store( key_C0, makeScalar( status, lower_limit_ ) );
    }
        
    double C_current = nd->Read(key_C0);//current concentration
    double dC_cumulative = array[4];
    array.Component(4, dC_cumulative + (C_current-C_last));//update cumulative change
        
    array.Component(0, t_clock); //current time stamp
    nd->Store(key_time, array);

    nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) );
  }
}


//update solution and check it against the specified range (with TDS)
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::Update_TDS(Event<dim>* event, double delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );    
  assert( nd->Status( key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( key_C0 ) != DIRICH) {  
    update_count_++;//recording     
    const VARIABLE_FLAG status(nd->Status( key_C0 ));
    
    // 1. starting with the sum of facet flux-concentration products stored in C1
    double ChangeRate = nd->Read(key_C1);   
    double solution = nd->Read(key_C0);//concentration
    // 2. correcting this sum for div vD using 'flux balance'   
    ChangeRate -= solution * nd->Read( key_FB ); 
    // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
    double new_solution = solution - (delta_t/nd->Read(key_fvPV)) * ChangeRate;
    // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
    // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
    const double source(nd->Read( key_NQC)); //kg/m3/s
    // new concentration
    new_solution += source * delta_t;
 
    if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(key_C0, makeScalar( status, new_solution ));//store solution value to key_C0
    else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
        if ( new_solution > upper_limit_ ) nd->Store( key_C0, makeScalar( status, upper_limit_ ) );
        else if ( new_solution < lower_limit_ ) nd->Store( key_C0, makeScalar( status, lower_limit_ ) );
    }
    
    nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) ); 
  }
}



//Synchronize neighbor nodes/FVs
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::Synchronize(Event<dim>* event, double t_clock, double& t_remove )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL ); 
  assert( nd->Status( key_C0 ) != DIRICH);
  if(nd  != NULL && nd->Status( key_C0 ) != DIRICH){     
    nd->Store( key_synchronize, makeScalar( nd->Status(key_synchronize), nd->Read(key_synchronize) + 1 ) );
    event->valid(false);
    ArrayVariable array;
    nd->Read(key_time, array);
    array.Component(4, 0.); //reset cumulative change of solution
    nd->Store(key_time, array);
    for ( auto n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* neighbor_node = nd->Neighbor(n);
        if( neighbor_node != NULL && neighbor_node->Status( key_C0 ) != DIRICH){
            auto index = neighbor_node->Read(key_EventIndex);
            if(index >= 0 && index < FullList.size()){
                Event<dim>* neighbor_event = FullList[index];  // TODO: deal with implicit type conversion
                assert( neighbor_event  != NULL ); 
                if (neighbor_event != NULL && neighbor_event->inPEPStack() == false) {
                    PEPList.push_back(neighbor_event);
                    neighbor_event->inPEPStack(true);
                    Update_DES(neighbor_event,t_clock);
                    ArrayVariable neighbor_array;
                    neighbor_node->Read(key_time, neighbor_array); 
                    double dC_cumulative = neighbor_array[4];//cumulative change of solution
                    double dC_target = neighbor_array[5];//target change of solution
                    if (fabs(dC_cumulative) >= fabs(dC_target)) {
                        #if defined(_OPENMP)
                        double t_begin = omp_get_wtime();
                        #else
                        clock_t t_begin = clock();
                        #endif
                        if (neighbor_event->inQueue()){
                            Heap_Node* neighbor_heap_node = HeapNodeFullList[index]; // TODO: deal with implicit type conversion
                            EventHeap.remove(neighbor_heap_node);
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



//advect variable with TDS (time-driven simulation), serial mode
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::AdvectVariable_TDS_serial( double time_interval )
{
    cout<<"Start DESAdvectionDiffusion<dim>::AdvectVariable_TDS_serial "<<endl;

    double begin=clock();
    double time_increment(time_interval); 
    
    clock_t T_begin= clock();
    const typename vector<Event<dim>*>::iterator stack_end(PEPList.end());
    for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )   
    {     
        ComputeRateofChange((*it));  
        ArrayVariable array;
        (*it)->getNode()->Read(key_time, array);
        double dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*CFL_multiplier_);
    }
    T_RateOfChange_ += clock() - T_begin; 

    const double one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double time(0.);
    
    while (time < time_interval)
    {
        cout <<"\n\tadvection (sub)step: "<< substep <<" of total steps "<<max(floor(time_interval/time_increment),one)<< endl;
        T_begin= clock();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };
        T_Update_ += clock() - T_begin;
        
        T_begin= clock();
        double new_time_increment(time_interval); 
        for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )
        { 
            ComputeRateofChange((*it));
            ArrayVariable array2;
            (*it)->getNode()->Read(key_time, array2);
            double dt_CFL = array2[2];//CFL time increment
            new_time_increment=min(new_time_increment, dt_CFL*CFL_multiplier_);
        };
        T_RateOfChange_ += clock() - T_begin; 

        time += time_increment;
        substep++;
    };
    
    T_AdvectVariable_+= clock() - begin;

    cout <<"Finish DESAdvectionDiffusion<dim>::AdvectVariable_TDS_serial "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_ /double(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double(CLOCKS_PER_SEC) << endl;
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double(CLOCKS_PER_SEC) << endl;
    cout <<"T_InsertToHeap_ = "<< T_InsertToHeap_ /double(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RemoveFromHeap_ = "<< T_RemoveFromHeap_ /double(CLOCKS_PER_SEC) << endl; 
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double(CLOCKS_PER_SEC) << endl;
}



#if defined(_OPENMP)
//advect variable with TDS (time-driven simulation), parallel mode
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::AdvectVariable_TDS_parallel( double time_interval )
{
    cout<<"Start DESAdvectionDiffusion<dim>::AdvectVariable_TDS_parallel "<<endl;
    cout<<"Using threads = "<<num_threads<<" Maximum available threads ="<< omp_get_max_threads() << endl;

    double begin=omp_get_wtime();
    double time_increment(time_interval); 
    
    double T_begin;
    T_begin = omp_get_wtime();
    
    size_t PEPList_size = PEPList.size();
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(dynamic)
        for( size_t i{0U}; i < PEPList_size; ++i)
        {
            auto it = PEPList.begin()+i;
            ComputeRateofChange((*it));  
        };
    }

    for( size_t i{0U}; i < PEPList_size; ++i)
    {
        auto it = PEPList.begin()+i;
        Event<dim>* event = *it;     
        ArrayVariable array;
        (*it)->getNode()->Read(key_time, array);
        double dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*CFL_multiplier_);
    }
    
    T_RateOfChange_ += omp_get_wtime() - T_begin;

    const double one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t substep(1);
    double time(0.);
    
    while (time < time_interval)
    {
        cout <<"\n\tadvection (sub)step: "<< substep <<" of total steps "<<max(floor(time_interval/time_increment),one)<< endl;
        T_begin = omp_get_wtime();
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };
        T_Update_ += omp_get_wtime() - T_begin;
        
        T_begin = omp_get_wtime();
        time_increment = time_interval; 
        #pragma omp parallel num_threads(num_threads)
        {        
            #pragma omp for schedule(dynamic)     
            for( size_t i{0U}; i < PEPList_size; ++i)
            {
                auto it = PEPList.begin()+i;
                ComputeRateofChange((*it));;
            }
        }        
        
        for( size_t i{0U}; i < PEPList_size; ++i)
        {
            auto it = PEPList.begin()+i;
            ArrayVariable array2;
            (*it)->getNode()->Read(key_time, array2);
            double dt_CFL = array2[2];//CFL time increment
            time_increment=min(time_increment, dt_CFL*CFL_multiplier_);
        };

        T_RateOfChange_ += omp_get_wtime() - T_begin; 

        time += time_increment;
        substep++;
    };
    
    T_AdvectVariable_+= omp_get_wtime() - begin;

    cout <<"Finish DESAdvectionDiffusion<dim>::AdvectVariable_TDS_parallel "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_ << endl;
    cout <<"T_Update_  = "<< T_Update_ << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_ << endl;
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ << endl;
    cout <<"T_InsertToHeap_ = "<< T_InsertToHeap_ << endl; 
    cout <<"T_RemoveFromHeap_ = "<< T_RemoveFromHeap_ << endl; 
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ << endl;
}
#endif





//advect variable with DES (discrete event simulation), serial mode
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::AdvectVariable_DES_serial( double model_time )
{
    double begin=clock();
    cout<<"Start DESAdvectionDiffusion<dim>::AdvectVariable_DES_serial "<<endl;
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
    
    while (!Finished)
    { 
        clock_t T_begin;
        const typename vector<Event<dim>*>::iterator stack_end(PEPList.end());
        for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )          
        {   
            Event<dim>* event = *it;
            T_begin= clock();
            ComputeRateofChange(event);  
            T_RateOfChange_ += clock() - T_begin;          
            if (event->valid() == false) {
                T_begin= clock();
                bool isactive = Schedule(event, model_time );
                T_Schedule_ += clock() - T_begin; 
                if (isactive) {
                    T_begin= clock();
                    double scheduled_time = event->t_schedule();
                    size_t index = event->getNode()->Read(key_EventIndex); // TODO: deal with implicit type conversion
                    Heap_Node* heap_node = new Heap_Node(scheduled_time,index);
                    EventHeap.insert(heap_node);
                    HeapNodeFullList[index] = heap_node;
                    event->inQueue(true);
                    T_InsertToHeap_ += clock() - T_begin; 
                }
            }
            event->inPEPStack(false);          
        };   

        if (EventHeap.empty()) time=model_time;
        else time = EventHeap.minimum()->getK();
        
        cout<<"  time = "<<time<<" model_time = "<<model_time<<" PEPList size = "<< PEPList.size() <<" Queue size = "<< EventHeap.size()<<endl;

        if (time == model_time) {
            Finished = true;
            
            const typename vector<Event<dim>*>::iterator End(PEPList.end());
            for ( typename vector<Event<dim>*>::iterator e=PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);        
            break;
        };

        PEPList.clear();
    
        double dt_PEP=numeric_limits<double>::max();
        size_t count = 0U;
        while (!EventHeap.empty())
        {
            Heap_Node* root_node = EventHeap.minimum();
            size_t top_index = root_node->getV();
            Event<dim>* top_event = FullList[top_index];
            
            if(top_event->valid() == false) {
                T_begin= clock();
                EventHeap.remove(root_node);
                top_event->inQueue(false);
                T_RemoveFromHeap_ += clock() - T_begin;
                continue;
            }
            
            count++;            
            ArrayVariable array;
            top_event->getNode()->Read(key_time, array);
            double dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, PEP_multiplier_*dt_target);
            double t_schedule = array[1];//scheduled time stamp
            if (t_schedule > (time+dt_PEP)) break;            
            if (top_event->inPEPStack() == false) {
                PEPList.push_back(top_event);
                top_event->inPEPStack(true);
                T_begin= clock();
                Update_DES(top_event,time);
                T_Update_ += clock() - T_begin;                
            }; 

            T_begin= clock();
            EventHeap.remove(root_node); 
            top_event->inQueue(false);
            T_RemoveFromHeap_ += clock() - T_begin;
            
            T_begin= clock();
            double t_remove(0.);    
            Synchronize(top_event,time,t_remove);
            T_RemoveFromHeap_ += t_remove;
            T_Synchronize_ += clock() - T_begin - t_remove;
        };
        
        //cout<<"  iteration count = "<<count<<endl;
    }
    T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DESAdvectionDiffusion<dim>::AdvectVariable_DES_serial "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_ /double(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double(CLOCKS_PER_SEC) << endl;
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double(CLOCKS_PER_SEC) << endl;
    cout <<"T_InsertToHeap_ = "<< T_InsertToHeap_ /double(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RemoveFromHeap_ = "<< T_RemoveFromHeap_ /double(CLOCKS_PER_SEC) << endl; 
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double(CLOCKS_PER_SEC) << endl;
}   



#if defined(_OPENMP)
//advect variable with DES (discrete event simulation), parallel mode
template<uint32_t dim>
void DESAdvectionDiffusion<dim>::AdvectVariable_DES_parallel( double model_time, size_t num_threads)
{
    double begin=omp_get_wtime();
    cout<<"Start DESAdvectionDiffusion<dim>::AdvectVariable_DES_parallel "<<endl;
    cout << "Using threads = "<<num_threads<<" Maximum available threads ="<< omp_get_max_threads() << endl;
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

    while (!Finished)
    {   
        double T_begin;
        T_begin = omp_get_wtime();
                
        size_t PEPList_size = PEPList.size();
           
        std::vector<Event<dim>*> tempList;
        #pragma omp parallel num_threads(num_threads)
        {
            std::vector<Event<dim>*> privateList;
            
            #pragma omp for schedule(dynamic)
            for( size_t i{0U}; i < PEPList_size; ++i)
            {
                auto it = PEPList.begin()+i;
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
        for( size_t i{0U}; i < tempList_size; ++i)
        {
            auto it = tempList.begin()+i;
            Event<dim>* event = *it;                 
            double scheduled_time = event->t_schedule();
            size_t index = event->getNode()->Read(key_EventIndex);
            Heap_Node* heap_node = new Heap_Node(scheduled_time,index);                    
            EventHeap.insert(heap_node);
            HeapNodeFullList[index] = heap_node;
            event->inQueue(true);
        }
        tempList.clear();
        
        T_RateOfChange_ += omp_get_wtime() - T_begin; 
            
        if (EventHeap.empty()) time=model_time;
        else time = EventHeap.minimum()->getK();
        
        cout<<"  time = "<<time<<" model_time = "<<model_time<<" PEPList size = "<< PEPList.size() <<" Queue size = "<< EventHeap.size()<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Event<dim>*>::iterator End(PEPList.end());
            for ( typename vector<Event<dim>*>::iterator e=PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        PEPList.clear();
    
        double dt_PEP=numeric_limits<double>::max();
        size_t count = 0U;
        while (!EventHeap.empty())
        {           
            Heap_Node* root_node = EventHeap.minimum();
            size_t top_index = root_node->getV();
            Event<dim>* top_event = FullList[top_index];
            
            if(top_event->valid() == false) {
                T_begin= omp_get_wtime();
                EventHeap.remove(root_node);
                top_event->inQueue(false);
                T_RemoveFromHeap_ += omp_get_wtime() - T_begin;
                continue;
            }
            
            count++;            
            ArrayVariable array;
            top_event->getNode()->Read(key_time, array);
            double dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, PEP_multiplier_*dt_target);
            double t_schedule = array[1];//scheduled time stamp
            if (t_schedule > (time+dt_PEP)) break;            
            if (top_event->inPEPStack() == false) {
                PEPList.push_back(top_event);
                top_event->inPEPStack(true);
                T_begin= omp_get_wtime();
                Update_DES(top_event,time);
                T_Update_ += omp_get_wtime() - T_begin;                
            };  
            
            T_begin= omp_get_wtime();
            EventHeap.remove(root_node); 
            top_event->inQueue(false);
            T_RemoveFromHeap_ += omp_get_wtime() - T_begin;
            
            T_begin= omp_get_wtime();
            double t_remove(0.);    
            Synchronize(top_event,time,t_remove);
            T_RemoveFromHeap_ += t_remove;
            T_Synchronize_ += omp_get_wtime() - T_begin - t_remove;
        };
            
        cout<<"  iteration count = "<<count<<endl;
    };
    /*
    //reset all events and add them to PEPList (for advection at next integration step)
    EventHeap.clear();
    PEPList.clear();
    const typename vector<Event<dim>*>::iterator stack_end(FullList.end());
    for ( typename vector<Event<dim>*>::iterator it=FullList.begin(); it!=stack_end; ++it )       
    {  
        Event<dim>* event = *it;
        PEPList.push_back(event);
        event->inPEPStack(true);          
        event->valid(false);
        event->inQueue(false);
    } 
    */
    T_AdvectVariable_+= omp_get_wtime() - begin;

    cout<<"Finish DESAdvectionDiffusion<dim>::AdvectVariable_DES_parallel "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_  << endl;
    cout <<"T_Update_  = "<< T_Update_  << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_<< endl;
    cout <<"T_RateOfChange_(including_T_Schedule_) = "<< T_RateOfChange_ << endl;
    cout <<"T_InsertToHeap_ = "<< T_InsertToHeap_ << endl; 
    cout <<"T_RemoveFromHeap_ = "<< T_RemoveFromHeap_ << endl; 
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ << endl;
} 
#endif



template class DESAdvectionDiffusion<1U>;
template class DESAdvectionDiffusion<2U>;
template class DESAdvectionDiffusion<3U>;

} // end csmp  
    
    
