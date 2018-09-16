#include "DESTransport.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim>
DESTransport<dim>::DESTransport( Model<dim>& m, const char* target_region )
  : variables::VariableSet_TracerTransfer(m.Database()),
    gref_(m.Region(target_region)),
    upper_limit_(1.), lower_limit_(0.), rate_count_(0U), update_count_(0U),
    T_RateOfChange_(0.), T_Schedule_(0.), T_InsertToHeap_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromHeap_(0.), T_AdvectVariable_(0.),
    first_step_(true)
{
    m.InstantiateFiniteVolumes();
    initializeVariablsAndKeys(m);
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_C0), lower_limit_, upper_limit_ );
    cout<<"DESTransport constructed"<<endl;
}


template<size_t dim>
void DESTransport<dim>::initializeVariablsAndKeys(Model<dim>& m)
{
    //creating new variables if not defined yet from input file   
    if(!m.Database().IsDefined("event index")) m.CreateProperty( "event index", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("update count")) m.CreateProperty( "update count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("rate count")) m.CreateProperty( "rate count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("schedule count")) m.CreateProperty( "schedule count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("synchronize count")) m.CreateProperty( "synchronize count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("concentration timing array")) m.CreateProperty( "concentration timing array", "none", ARRAY, NODE, 6, -1.00E+10 ,1.00E+10);
    
    // model-wide initialisation
    m.Region("Model").InputPropertyValue( "new concentration", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );    
    
    //assigning keys 
    key_EventIndex = INDEX<SCALAR,NODE>( m.Database().StorageKey("event index") );
    key_update = INDEX<SCALAR,NODE>( m.Database().StorageKey("update count") );
    key_rate = INDEX<SCALAR,NODE>( m.Database().StorageKey("rate count") );
    key_schedule = INDEX<SCALAR,NODE>( m.Database().StorageKey("schedule count") );
    key_synchronize = INDEX<SCALAR,NODE>( m.Database().StorageKey("synchronize count") );
    key_time = INDEX<ARRAY,NODE>( m.Database().StorageKey("concentration timing array") );    
    
    //checking keys
    if ( key_EventIndex.place != NODE || key_EventIndex.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'event index' variable must be SCALAR and placed on NODE"  );       
    if ( key_update.place != NODE || key_update.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( key_rate.place != NODE || key_rate.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( key_schedule.place != NODE || key_schedule.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( key_synchronize.place != NODE || key_synchronize.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  );    
    if ( key_time.place != NODE || key_time.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'concentration timing array' variable must be ARRAY and placed on NODE"  );                
}


template<size_t dim>
void DESTransport<dim>::initializeFiniteVolumeProperties()
 {
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref_.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref_.ElementsBegin(); it!=it_end; ++it )
    {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // computing sector pore volumes
         double64 phi = (*it)->Read( this->key_THI);
         const double64 thickness = (*it)->Read( this->key_THI );
         if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
         
         for ( size_t i=0U; i<sectors; ++i ) {
              const double64 sector_volume = (*it)->SectorVolume(i);       
              double64 pore_volume   = (*it)->N(i)->Read( this->key_FVPV );
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( this->key_FVPV, makeScalar(PLAIN,pore_volume) );
         }

         // computing facet normals and areas
         for ( size_t j=0U; j<facets; ++j ) {
              const double64 facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, this->key_fA, makeScalar( PLAIN, facet_area ) );
              Point<dim> nrml = (*it)->FacetNormal(j);
              VectorVariable<dim>  fnrml;
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, this->key_fn, fnrml );
         }            
   }

   // initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // (here the pore volumes do not include the sectors outside the region)
   const typename vector<Node<dim>*>::iterator nit_end(gref_.NodesEnd());
   for ( typename vector<Node<dim>*>::iterator nit=gref_.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const size_t parent_elements((*nit)->Parents());      
        for ( size_t i=0U; i<parent_elements; ++i ) {
             Element<dim>* const eptr = (*nit)->Parent(i);
             // computing facet normals and areas
             const size_t facets(eptr->Facets());
             for ( size_t j=0U; j<facets; ++j ) {
                  const double64 facet_area = eptr->FacetArea(j);
                  eptr->Store( j, 0U, this->key_fA, makeScalar( PLAIN, facet_area ) );
                  Point<dim> nrml = eptr->FacetNormal(j);
                  VectorVariable<dim>  fnrml;
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, this->key_fn, fnrml );
             }             
        }
   }
 } // end initializeFiniteVolumeProperties




//Compute the flux balance and CFL number for a node/FV
template<size_t dim>
void DESTransport<dim>::ComputeFluxBalanceAndCFL( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_C0 ) != DIRICH);

  if(nd  != NULL && nd->Status(  this->key_C0 ) != DIRICH) {
    double64 flux_balance(0.), outflow(0.);
    
    VectorVariable<dim> vD, facetNrml;
    const size_t node_parent_elements(nd->Parents());
    
    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        const size_t pnid(nd->ParentNodeNumber(t));
        eptr->Read( this->key_V, vD);

        const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
        for ( size_t i=0U; i<sector_facets; i++ )
        {
            const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
            const size_t inside_node(eptr->FV()->InsideNode(iFacet));
                        
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
        }
    }
        
    nd->Store(  this->key_FB, makeScalar( nd->Status(  this->key_FB ), flux_balance ) );//flux balance
    
    //CFL time increment 
    ArrayVariable array2;
    nd->Read(key_time, array2);
    if (outflow < numeric_limits<double64>::epsilon())
        array2.Component(2, numeric_limits<double64>::max());
    else 
        array2.Component(2, nd->Read(  this->key_FVPV ) / outflow);  
    nd->Store(key_time, array2);
  }
}



//Compute the rate of change in a node/FV
template<size_t dim>
void DESTransport<dim>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status(  this->key_C0 ) != DIRICH) {
    rate_count_++;//recording
    nd->Store(  key_rate, makeScalar( nd->Status( key_rate), nd->Read( key_rate) + 1 ) );

    double64 accumulation(0.);
    VectorVariable<dim> vD, facetNrml;
    const size_t node_parent_elements(nd->Parents());
    
    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        const size_t pnid(nd->ParentNodeNumber(t));
        eptr->Read( this->key_V, vD);

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
            // finding the upstream node
            const size_t upstream_node = (vD_n < 0.) ? outside_node : inside_node;
            const double64 C_upstream = eptr->N(upstream_node)->Read( this->key_C0 );
            accumulation += facet_flux * C_upstream;  
        }
    }    

    nd->Store(  this->key_C1, makeScalar( nd->Status( this->key_C1 ), accumulation ) );    
  }
} 




//schedule an event associated with a node/FV
template<size_t dim>
bool DESTransport<dim>::Schedule(Event<dim>* event, double64 t_end, double64 cfl_multiplier)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status(  this->key_C0 ) != DIRICH) {    
    ArrayVariable array;
    nd->Read(key_time, array);

    nd->Store( key_schedule, makeScalar( nd->Status( key_schedule), nd->Read( key_schedule) + 1 ) );
    event->valid(true);
    //compute target change
    double64 CFL = array[2];//CFL number
    double64 PV = nd->Read( this->key_FVPV);//Pore volume
    double64 ChangeRate = nd->Read( this->key_C1);//rate of change
    double64 C0 = nd->Read(  this->key_C0);//concentration
    double64 flux_balance = nd->Read(  this->key_FB);//flux balance

    double64 dC_CFL = -CFL*cfl_multiplier/PV*(ChangeRate-C0*flux_balance);//targe change

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
        event->t_schedule(t_end);
        array.Component(1, t_end);//schedule time stamp
        nd->Store(key_time, array);
        return false;
    } else {
        event->t_schedule(t_current + dt_target);
        array.Component(1, t_current + dt_target);//schedule time stamp
        nd->Store(key_time,array);
        return true;
    };
  }
}


//update solution and check it against the specified range (with DES)
template<size_t dim>
void DESTransport<dim>::Update_DES(Event<dim>* event, double64 t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status(  this->key_C0 ) != DIRICH) {  
    update_count_++;//recording  
    const VARIABLE_FLAG status(nd->Status( this->key_C0 ));
    ArrayVariable array;
    nd->Read(key_time, array);

    // 1. starting with the sum of facet flux-concentration products stored in C1
    double64 ChangeRate = nd->Read( this->key_C1);   
    double64 solution = nd->Read( this->key_C0);//concentration
    // 2. correcting this sum for div vD using 'flux balance'   
    ChangeRate -= solution * nd->Read( this->key_FB ); 
    // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
    double64 t_current = array[0]; //current time stamp
    double64 new_solution = solution - ((t_clock - t_current)/nd->Read( this->key_FVPV)) * ChangeRate;
    // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
    // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
    const double64 source(nd->Read(this->key_NQV));
    // new concentration
    new_solution += source * (t_clock - t_current);

    double64 C_last = nd->Read(this->key_C0);//last concentration
        
    if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(this->key_C0, makeScalar( status, new_solution ));//store solution value to this->key_C0
    else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
        if ( new_solution > upper_limit_ ) nd->Store( this->key_C0, makeScalar( status, upper_limit_ ) );
        else if ( new_solution < lower_limit_ ) nd->Store( this->key_C0, makeScalar( status, lower_limit_ ) );
    }
        
    double64 C_current = nd->Read(this->key_C0);//current concentration
    double64 dC_cumulative = array[4];
    array.Component(4, dC_cumulative + (C_current-C_last));//update cumulative change
        
    array.Component(0, t_clock); //current time stamp
    nd->Store(key_time, array);

    nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) );
  }
}


//update solution and check it against the specified range (with TDS)
template<size_t dim>
void DESTransport<dim>::Update_TDS(Event<dim>* event, double64 delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );    
  assert( nd->Status(  this->key_C0 ) != DIRICH);
    
  if(nd  != NULL && nd->Status(  this->key_C0 ) != DIRICH) {  
    update_count_++;//recording     
    const VARIABLE_FLAG status(nd->Status( this->key_C0 ));
    
    // 1. starting with the sum of facet flux-concentration products stored in C1
    double64 ChangeRate = nd->Read(this->key_C1);   
    double64 solution = nd->Read(this->key_C0);//concentration
    // 2. correcting this sum for div vD using 'flux balance'   
    ChangeRate -= solution * nd->Read(  this->key_FB ); 
    // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
    double64 new_solution = solution - (delta_t/nd->Read( this->key_FVPV)) * ChangeRate;
    // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
    // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
    const double64 source(nd->Read( this->key_NQV));
    // new concentration
    new_solution += source * delta_t;
 
    if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(this->key_C0, makeScalar( status, new_solution ));//store solution value to this->key_C0
    else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
        if ( new_solution > upper_limit_ ) nd->Store( this->key_C0, makeScalar( status, upper_limit_ ) );
        else if ( new_solution < lower_limit_ ) nd->Store( this->key_C0, makeScalar( status, lower_limit_ ) );
    }
    
    nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) ); 
  }
}



//Synchronize neighbor nodes/FVs
template<size_t dim>
void DESTransport<dim>::Synchronize(Event<dim>* event,double64 t_clock,double64& t_remove )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL ); 
  assert( nd->Status(  this->key_C0 ) != DIRICH);
  if(nd  != NULL && nd->Status( this->key_C0 ) != DIRICH){     
    nd->Store( key_synchronize, makeScalar( nd->Status(key_synchronize), nd->Read(key_synchronize) + 1 ) );
    event->valid(false);
    ArrayVariable array;
    nd->Read(key_time, array);
    array.Component(4, 0.); //reset cumulative change of solution
    nd->Store(key_time, array);
    for ( size_t n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* neighbor_node = nd->Neighbor(n);
        if( neighbor_node != NULL && neighbor_node->Status(  this->key_C0 ) != DIRICH){
            size_t index = neighbor_node->Read(key_EventIndex);
            if(index >= 0 && index < FullList.size()){
                Event<dim>* neighbor_event = FullList[index];  
                assert( neighbor_event  != NULL ); 
                if (neighbor_event != NULL && neighbor_event->inPEPStack() == false) {
                    PEPList.push_back(neighbor_event);
                    neighbor_event->inPEPStack(true);
                    Update_DES(neighbor_event,t_clock);
                    ArrayVariable neighbor_array;
                    neighbor_node->Read(key_time, neighbor_array); 
                    double64 dC_cumulative = neighbor_array[4];//cumulative change of solution
                    double64 dC_target = neighbor_array[5];//target change of solution
                    if (fabs(dC_cumulative) >= fabs(dC_target)) {
                        #if defined(_OPENMP)
                        double64 t_begin = omp_get_wtime();
                        #else
                        clock_t t_begin = clock();
                        #endif
                        if (neighbor_event->inQueue()){
                            Heap_Node* neighbor_heap_node = HeapNodeFullList[index];
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



//advect variable with TDS (time-driven simulation)
template<size_t dim>
void DESTransport<dim>::AdvectVariable_TDS( double64 time_interval, double64 cfl_multiplication_factor, double64 PEP_parameter )
{
    if(first_step_){
        initializeFiniteVolumeProperties();
        //create events for all nodes and add them to PEPList
        size_t dirich_count = 0;
        const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
        for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
        { 
            (*nit)->Store( key_update, makeScalar( (*nit)->Status( key_update), 0 ) ); //update count
            (*nit)->Store( key_rate, makeScalar( (*nit)->Status( key_rate), 0 ) ); //changerate count
            (*nit)->Store( key_schedule, makeScalar( (*nit)->Status( key_schedule), 0 ) ); //schedule count
            (*nit)->Store( key_synchronize, makeScalar( (*nit)->Status( key_synchronize), 0 ) ); //synchronize count   
        
            ArrayVariable arrayVariable( 6, 0., PLAIN );         
            (*nit)->Store( key_time, arrayVariable );  
                   
            //initializeFiniteVolumeProperties(*nit);
            if((*nit)->Status(  this->key_C0 ) != DIRICH) {
                Event<dim>* event = new Event<dim>(*nit);
                event->valid(false);
                ComputeFluxBalanceAndCFL(event);
                PEPList.push_back(event);
                event->inPEPStack(true);
            } else {
                dirich_count++;
            }
        }
        cout<<FullList.size()<<" events created for all nodes, excluding "<<dirich_count<<" DIRICH nodes"<<endl;
        first_step_=false;
    } 


    double64 begin=clock();
    double64 time_increment(time_interval); 
    
    clock_t T_begin= clock();
    const typename vector<Event<dim>*>::iterator stack_end(PEPList.end());
    for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )   
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
        for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };
        T_Update_ += clock() - T_begin;
        
        T_begin= clock();
        double64 new_time_increment(time_interval); 
        for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )
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

    cout<<"Finish DESTransport<dim>::AdvectVariable_TDS "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_InsertToHeap_ = "<< T_InsertToHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RemoveFromHeap_ = "<< T_RemoveFromHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl;
}



//advect variable with DES (discrete event simulation)
template<size_t dim>
void DESTransport<dim>::AdvectVariable_DES( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter, size_t num_threads )
{
#if defined(_OPENMP )
    if (num_threads <= 0) {
        cerr <<"WARNING: input number of threads is less than 1, reset to 1"<<endl;
        num_threads = 1;
        AdvectVariable_DES_openmp( model_time, cfl_multiplication_factor, PEP_parameter, num_threads);
    } 
    else if (num_threads > omp_get_max_threads()) {
        cerr <<"WARNING: input number of threads is larger than maximum available threads (" << omp_get_max_threads() << "), reset to "<< omp_get_max_threads() <<endl;
        num_threads = omp_get_max_threads();
        AdvectVariable_DES_openmp( model_time, cfl_multiplication_factor, PEP_parameter, num_threads); 
    }  
    else if (num_threads == 1) {
        AdvectVariable_DES_serial( model_time, cfl_multiplication_factor, PEP_parameter );
    } 
    else {
    AdvectVariable_DES_openmp( model_time, cfl_multiplication_factor, PEP_parameter, num_threads); 
    }            
#else
    if (num_threads > 1)
        cerr <<"WARNING: OpenMP is not available, using serial mode"<<endl;
    AdvectVariable_DES_serial( model_time, cfl_multiplication_factor, PEP_parameter );
#endif
}




#if defined(_OPENMP)
//advect variable with DES (discrete event simulation), with openmp
template<size_t dim>
void DESTransport<dim>::AdvectVariable_DES_openmp( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter, size_t num_threads)
{
    double64 begin=omp_get_wtime();
    cout<<"Start DESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_openmp "<<endl;
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
    
    if(first_step_){
        initializeFiniteVolumeProperties();
        //create events for all nodes and add them to event lists
        size_t index = 0;
        size_t dirich_count = 0;
        const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
        for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
        { 
            (*nit)->Store( key_update, makeScalar( (*nit)->Status( key_update), 0 ) ); //update count
            (*nit)->Store( key_rate, makeScalar( (*nit)->Status( key_rate), 0 ) ); //changerate count
            (*nit)->Store( key_schedule, makeScalar( (*nit)->Status( key_schedule), 0 ) ); //schedule count
            (*nit)->Store( key_synchronize, makeScalar( (*nit)->Status( key_synchronize), 0 ) ); //synchronize count   
        
            ArrayVariable arrayVariable( 7, 0., PLAIN );         
            (*nit)->Store( key_time, arrayVariable );  
                   
            if((*nit)->Status(  this->key_C0 ) != DIRICH) {
                (*nit)->Store( key_EventIndex, makeScalar( (*nit)->Status(key_EventIndex), index) );//event index 
                Event<dim>* event = new Event<dim>(*nit);
                PEPList.push_back(event);
                event->inPEPStack(true);
                ComputeRateofChange(event);
                Schedule(event, model_time, cfl_multiplication_factor);
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
        first_step_=false;
    }  

    while (!Finished)
    {   
        double64 T_begin;
        
        T_begin = omp_get_wtime();
                
        size_t PEPList_size = PEPList.size();
           
        std::vector<Event<dim>*> tempList;
        #pragma omp parallel num_threads(num_threads)
        {
            std::vector<Event<dim>*> privateList;
            
            #pragma omp for schedule(dynamic)
            for(size_t i = 0U; i < PEPList_size; ++i)
            {
                auto it = PEPList.begin()+i;
                Event<dim>* event = *it;                                
                ComputeRateofChange((*it)); 
                
                if (event->valid() == false) 
                    if (Schedule(event, model_time, cfl_multiplication_factor)) privateList.push_back(event);
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
            size_t index = event->getNode()->Read(key_EventIndex);
            Heap_Node* heap_node = new Heap_Node(scheduled_time,index);                    
            EventHeap.insert(heap_node);
            HeapNodeFullList[index] = heap_node;
            event->inQueue(true);
        }
        tempList.clear();
        
        
        T_RateOfChange_ += omp_get_wtime() - T_begin; 
        cout <<"  PEPList size = " << PEPList.size() << "  Queue size = "<< EventHeap.size()<<endl;          

            
        if (EventHeap.empty()) time=model_time;
        else time = EventHeap.minimum()->getK();
        cout<<"  time = "<<time<<" model_time = "<<model_time<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Event<dim>*>::iterator End(PEPList.end());
            for ( typename vector<Event<dim>*>::iterator e=PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);
            break;
        };

        PEPList.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
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
            double64 dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, PEP_parameter*dt_target);
            double64 t_schedule = array[1];//scheduled time stamp
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
            double64 t_remove(0.);    
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

    cout<<"Finish DESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_openmp "<<endl;
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




//advect variable with DES (discrete event simulation)
template<size_t dim>
void DESTransport<dim>::AdvectVariable_DES_serial( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter )
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
    
    
    if(first_step_){
        initializeFiniteVolumeProperties();
        //create events for all nodes and add them to PEPStack and EntireQueue
        size_t index = 0;
        size_t dirich_count = 0;
        const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
        for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
        { 
            (*nit)->Store( key_update, makeScalar( (*nit)->Status( key_update), 0 ) ); //update count
            (*nit)->Store( key_rate, makeScalar( (*nit)->Status( key_rate), 0 ) ); //changerate count
            (*nit)->Store( key_schedule, makeScalar( (*nit)->Status( key_schedule), 0 ) ); //schedule count
            (*nit)->Store( key_synchronize, makeScalar( (*nit)->Status( key_synchronize), 0 ) ); //synchronize count   
        
            ArrayVariable arrayVariable( 6, 0., PLAIN );         
            (*nit)->Store( key_time, arrayVariable );  
                   
            //initializeFiniteVolumeProperties(*nit);
            if((*nit)->Status(  this->key_C0 ) != DIRICH) {
                (*nit)->Store( key_EventIndex, makeScalar( (*nit)->Status(key_EventIndex), index) );//event index 
                Event<dim>* event = new Event<dim>(*nit);
                ComputeFluxBalanceAndCFL(event);
                PEPList.push_back(event);
                event->inPEPStack(true);
                ComputeRateofChange(event);
                Schedule(event, model_time, cfl_multiplication_factor);
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
        first_step_=false;
    }  
        
    
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
                bool isactive = Schedule(event, model_time, cfl_multiplication_factor);
                T_Schedule_ += clock() - T_begin; 
                if (isactive) {
                    T_begin= clock();
                    double64 scheduled_time = event->t_schedule();
                    size_t index = event->getNode()->Read(key_EventIndex);
                    Heap_Node* heap_node = new Heap_Node(scheduled_time,index);
                    EventHeap.insert(heap_node);
                    HeapNodeFullList[index] = heap_node;
                    event->inQueue(true);
                    T_InsertToHeap_ += clock() - T_begin; 
                }
            }
            event->inPEPStack(false);          
        };   
        cout <<"  PEPList size = " << PEPList.size() << "  Queue size = "<< EventHeap.size()<<endl;   

        if (EventHeap.empty()) time=model_time;
        else time = EventHeap.minimum()->getK();
        cout<<"  time = "<<time<<" model_time = "<<model_time<<endl;

        if (time == model_time) {
            Finished = true;
            
            const typename vector<Event<dim>*>::iterator End(PEPList.end());
            for ( typename vector<Event<dim>*>::iterator e=PEPList.begin(); e!=End; ++e )
                (*e)->valid(false);        
            break;
        };

        PEPList.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
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
            double64 dt_target = array[3];//target time stamp
            dt_PEP = min(dt_PEP, PEP_parameter*dt_target);
            double64 t_schedule = array[1];//scheduled time stamp
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
            double64 t_remove(0.);    
            Synchronize(top_event,time,t_remove);
            T_RemoveFromHeap_ += t_remove;
            T_Synchronize_ += clock() - T_begin - t_remove;
        };
        
        cout<<"  iteration count = "<<count<<endl;
    };
    T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DESTransport<dim>::AdvectVariable_DES_serial "<<endl;
    cout <<"rate_count_ = "<<rate_count_<<endl;
    cout <<"update_count_ = "<<update_count_<<endl; 
    cout <<"T_Schedule_ = "<< T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_InsertToHeap_ = "<< T_InsertToHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RemoveFromHeap_ = "<< T_RemoveFromHeap_ /double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl;
}   


template class DESTransport<1U>;
template class DESTransport<2U>;
template class DESTransport<3U>;

} // end csmp  
    
    
