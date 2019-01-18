#include "TwoPhaseDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "DenseMatrix.h"
#include "CSMP_mathUtilities.h"
#include "FlowFunctionsModule.h"
#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseDESTransport( Model<dim>& m,
                                                                 const char* target_region,
                                                                 FLOW_FUNCTIONS<dim>& ff,
                                                                 bool with_capillary_spreading,
                                                                 bool with_gravity_forces,
                                                                 bool tensor_k,
                                                                 double64 PEP_multiplier,
                                                                 double64 cfl_multiplier )
    : variables::VariableSet_CO2GeoSequestration(m.Database()),
      gref_(m.Region(target_region)),
      db_(m.Database()),
      with_capillary_spreading_(with_capillary_spreading),
      with_gravity_forces_(with_gravity_forces),
      tensor_k_(tensor_k),
      upper_limit_(1.), lower_limit_(0.), rate_count_(0U), update_count_(0U),
      T_RateOfChange_(0.), T_Schedule_(0.), T_InsertToHeap_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromHeap_(0.), T_AdvectVariable_(0.),
      first_step_(true),
      PEP_multiplier_(PEP_multiplier),
      CFL_multiplier_(cfl_multiplier),
      relaxing_factor_(10.),
      flowfunctions_(ff)
{
    m.InstantiateFiniteVolumes();
    InitializeVariablesAndKeys(m);

// TODO: perhaps only where you have to    calculatePermeabilityProjections(m.Region(target_region));
  
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_sCO2), lower_limit_, upper_limit_ );
    cout<<"TwoPhaseDESTransport constructed"<<endl;
} // end constructor  


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseDESTransport(  Model<dim>& m,
                                                                 const char* target_region,
                                                                 FLOW_FUNCTIONS<dim>& ff,
                                                                 bool with_capillary_spreading,
                                                                 bool with_gravity_forces,
                                                                 bool tensor_k,
                                                                 double64 PEP_multiplier,
                                                                 double64 cfl_multiplier,
                                                                 double64 relaxing_factor )
    : variables::VariableSet_CO2GeoSequestration(m.Database()),
      gref_(m.Region(target_region)),
      db_(m.Database()),
      with_capillary_spreading_(with_capillary_spreading),
      with_gravity_forces_(with_gravity_forces),
      tensor_k_(tensor_k),
      upper_limit_(1.), lower_limit_(0.), rate_count_(0U), update_count_(0U),
      T_RateOfChange_(0.), T_Schedule_(0.), T_InsertToHeap_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromHeap_(0.), T_AdvectVariable_(0.),
      first_step_(true),
      PEP_multiplier_(PEP_multiplier),
      CFL_multiplier_(cfl_multiplier),
      relaxing_factor_(relaxing_factor),
      flowfunctions_(ff)
{
    m.InstantiateFiniteVolumes();
    InitializeVariablesAndKeys(m);

// TODO: perhaps only where you have to        calculatePermeabilityProjections(m.Region(target_region));
  
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_sCO2), lower_limit_, upper_limit_ );
    cout<<"TwoPhaseDESTransport constructed"<<endl;
} // end constructor 




  

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::InitializeVariablesAndKeys(Model<dim>& m)
{
    //creating new variables if not defined yet from input file
    if(!m.Database().IsDefined("variation rate nonwetting phase")) m.CreateProperty( "variation rate nonwetting phase", "m3/(m3.s)", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08);
    if(!m.Database().IsDefined("event index")) m.CreateProperty( "event index", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("update count")) m.CreateProperty( "update count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("rate count")) m.CreateProperty( "rate count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("schedule count")) m.CreateProperty( "schedule count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("synchronize count")) m.CreateProperty( "synchronize count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("nonwetting phase timing array")) m.CreateProperty( "nonwetting phase timing array", "none", ARRAY, NODE, 7, -1.00E+10 ,1.00E+10);
    if(!m.Database().IsDefined("cfl multiplier")) m.CreateProperty( "cfl multiplier", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("shock saturation aqueous phase")) m.CreateProperty( "shock saturation aqueous phase", "none", SCALAR, ELEMENT, 1, -5.00E-02 ,1.05E+00);     
    if(!m.Database().IsDefined("saturation gradient")) m.CreateProperty( "saturation gradient", "none", VECTOR, ELEMENT, 3, -1.00E+08 ,1.00E+08);
    if(!m.Database().IsDefined("pressure gradient")) m.CreateProperty( "pressure gradient", "none", VECTOR, ELEMENT, 3, -1.00E+10 ,1.00E+10);
    if(!m.Database().IsDefined("truncated FV")) m.CreateProperty( "truncated FV", "none", SCALAR, NODE, 1, 0 ,1);
    
    // model-wide initialisation
    m.Region("Model").InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );    
    m.Region("Model").InputPropertyValue( "truncated FV", makeScalar(PLAIN,0), COMPLETE); 
    
    //assigning keys 
    key_dsnw = INDEX<SCALAR,NODE> ( m.Database().StorageKey("variation rate nonwetting phase") );
    key_EventIndex = INDEX<SCALAR,NODE>( m.Database().StorageKey("event index") );
    key_update = INDEX<SCALAR,NODE>( m.Database().StorageKey("update count") );
    key_rate = INDEX<SCALAR,NODE>( m.Database().StorageKey("rate count") );
    key_schedule = INDEX<SCALAR,NODE>( m.Database().StorageKey("schedule count") );
    key_synchronize = INDEX<SCALAR,NODE>( m.Database().StorageKey("synchronize count") );
    key_time = INDEX<ARRAY,NODE>( m.Database().StorageKey("nonwetting phase timing array") );
    key_CFL = INDEX<SCALAR,NODE>( m.Database().StorageKey("cfl multiplier") );
    key_gradSn = INDEX<VECTOR,ELEMENT>( m.Database().StorageKey("saturation gradient") );
    key_gradP = INDEX<VECTOR,ELEMENT>( m.Database().StorageKey("pressure gradient") );
    key_cut = INDEX<SCALAR,NODE> ( m.Database().StorageKey("truncated FV") );
    
    //checking keys 
    if ( key_dsnw.place != NODE || key_dsnw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'variation rate nonwetting phase' variable must be SCALAR and placed on NODE"  );    
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
    if ( key_CFL.place != NODE || key_CFL.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'cfl multiplier' variable must be SCALAR and placed on NODE"  );
    if ( key_gradSn.place != ELEMENT || key_gradSn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'saturation gradient' variable must be VECTOR and placed on ELEMENT"  );  
    if ( key_gradP.place != ELEMENT || key_gradP.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'pressure gradient' variable must be VECTOR and placed on ELEMENT"  );          
    if ( key_cut.place != NODE || key_cut.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'truncated FV' variable must be SCALAR and placed on NODE"  );                                                                   
}



template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::InitializeFiniteVolumeProperties()
 {
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref_.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref_.ElementsBegin(); it!=it_end; ++it )
    {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // computing sector pore volumes
         double64 phi = (*it)->Read( this->key_phi);
         const double64 thickness = (*it)->Read( this->key_thi );
         if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
         
         for ( size_t i=0U; i<sectors; ++i ) {
              const double64 sector_volume = (*it)->SectorVolume(i);       
              double64 pore_volume   = (*it)->N(i)->Read( this->key_fvPV );
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( this->key_fvPV, makeScalar(PLAIN,pore_volume) );
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
         
         //compute wetting phase saturation at shock
         double64 sw_shock = flowfunctions_.ShockHeight(*it);
         (*it)->Store( this->key_ssH2O, makeScalar( (*it)->Status( this->key_ssH2O), sw_shock ) );                
   }

   // initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // (here the pore volumes do not include the sectors outside the region)
   const typename vector<Node<dim>*>::iterator nit_end(gref_.NodesEnd());
   for ( typename vector<Node<dim>*>::iterator nit=gref_.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const size_t parent_elements((*nit)->Parents()); 
        bool truncated_node = false;   
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
             //compute wetting phase saturation at shock
             double64 sw_shock = flowfunctions_.ShockHeight(eptr);
             eptr->Store( this->key_ssH2O, makeScalar( eptr->Status( this->key_ssH2O), sw_shock ) );    
             //determine whether FV node is truncated
             if(!gref_.Contains(eptr)) { //parent elment located outside domain
                 halo_stencils_.insert(eptr);
                 truncated_node = true;
             }     
        }
        if (truncated_node) (*nit)->Store( key_cut, makeScalar( (*nit)->Status( key_cut), 1 ) );
   }
   
 } // end initializeFiniteVolumeProperties



//resest cfl multipliers to default value = CFL_multiplier_*relaxing_factor_ for all nodes
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::ResetCFLMultiplier()
{
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        (*nit)->Store( key_CFL, makeScalar( (*nit)->Status( key_CFL), CFL_multiplier_*relaxing_factor_ ) ); //CFL multiplier 
        ArrayVariable array;
        (*nit)->Read(key_time, array);
        array.Component(6, CFL_multiplier_*relaxing_factor_);
        (*nit)->Store( key_time, array);
    }
} 



//Compute non-wetting phase saturaiton gradient at parement elements, for capillary component computation.
// TODO: super expensive approach - use values from neighboring nodes 
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::ComputeGradients (Event<dim>* event )
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
            for ( size_t k=0U; k<dim; k++ ) {
                if(with_capillary_spreading_) snw_gradient(k) += DN(k,j) * sn;
                p_gradient(k) += -DN(k,j) * p;               
            }
        }
        
        if(with_capillary_spreading_) eptr->Store(this->key_gradSn, snw_gradient);
        eptr->Store(this->key_gradP, p_gradient);       
        
        }
    }
}




//Synchronize neighbor nodes/FVs
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::Synchronize(Event<dim>* event,double64 t_clock,double64& t_remove)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL ); 
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    nd->Store( key_synchronize, makeScalar( nd->Status(key_synchronize), nd->Read(key_synchronize) + 1 ) );
    event->valid(false);
    ArrayVariable array;
    nd->Read(key_time, array);
    array.Component(4, 0.); //reset cumulative change of solution
    nd->Store(key_time, array);
    for ( size_t n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* neighbor_node = nd->Neighbor(n);
        if( neighbor_node != NULL && neighbor_node->Status(  this->key_sCO2 ) != DIRICH){
            int index = static_cast<int>(neighbor_node->Read(key_EventIndex));
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
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS( double64 time_interval)
{
    if(first_step_){
        InitializeFiniteVolumeProperties();
        //create events for all nodes and add them to PEPList
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
                   
            if((*nit)->Status(  this->key_sCO2 ) != DIRICH) {
                Event<dim>* event = new Event<dim>(*nit);
                event->valid(false);
                //initializeFiniteVolumeProperties(event);
                PEPList.push_back(event);
                event->inPEPStack(true);
            } else {
                dirich_count++;
            }
        }
        cout<<PEPList.size()<<" events created for all nodes, excluding "<<dirich_count<<" DIRICH nodes"<<endl;
        first_step_=false;
    } 

    double64 begin=clock();
    double64 time_increment(time_interval); 
    
    clock_t T_begin= clock();
    const typename vector<Event<dim>*>::iterator stack_end(PEPList.end());
    for ( typename vector<Event<dim>*>::iterator it=PEPList.begin(); it!=stack_end; ++it )     
    {   
        //if (with_capillary_spreading_) ComputeSaturationGradient ((*it)); 
        ComputeGradients ((*it));
        ComputeRateofChange((*it));  
        ArrayVariable array;
        (*it)->getNode()->Read(key_time, array);
        double64 dt_CFL = array[2];//CFL time increment
        time_increment=min(time_increment, dt_CFL*CFL_multiplier_);
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
            ComputeGradients ((*it));
            ComputeRateofChange((*it));
            ArrayVariable array2;
            (*it)->getNode()->Read(key_time, array2);
            double64 dt_CFL = array2[2];//CFL time increment
            new_time_increment=min(new_time_increment, dt_CFL*CFL_multiplier_);
        };
        T_RateOfChange_ += clock() - T_begin; 

        time += time_increment;
        substep++;
    };
    
    T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS "<<endl;
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
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES( double64 model_time, size_t num_threads )
{
#if defined(_OPENMP)
    if (num_threads <= 0) {
        cerr <<"WARNING: input number of threads is less than 1, reset to 1"<<endl;
        num_threads = 1;
        //AdvectVariable_DES_serial( model_time );
        AdvectVariable_DES_openmp( model_time, num_threads); 
    } 
    else if (num_threads > omp_get_max_threads()) {
        cerr <<"WARNING: input number of threads is larger than maximum available threads (" << omp_get_max_threads() << "), reset to "<< omp_get_max_threads() <<endl;
        num_threads = omp_get_max_threads();
        AdvectVariable_DES_openmp( model_time, num_threads); 
    }  
    else if (num_threads == 1) {
        AdvectVariable_DES_serial( model_time );
    } 
    else {
    AdvectVariable_DES_openmp( model_time, num_threads); 
    }            
#else
    if (num_threads > 1)
        cerr <<"WARNING: OpenMP is not available, using serial mode"<<endl;
    AdvectVariable_DES_serial( model_time );
#endif
}



//advect variable with DES (discrete event simulation), serial version
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial( double64 model_time)
{
    double64 begin=clock();
    cout<<"Start DESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial "<<endl;
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

    ResetCFLMultiplier();
    if(first_step_){
        InitializeFiniteVolumeProperties();
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
                   
            if((*nit)->Status(  this->key_sCO2 ) != DIRICH) {
                (*nit)->Store( key_EventIndex, makeScalar( (*nit)->Status(key_EventIndex), index) );//event index 
                Event<dim>* event = new Event<dim>(*nit);
                PEPList.push_back(event);
                event->inPEPStack(true);
                ComputeGradients (event );
                ComputeRateofChange(event);
                Schedule(event, model_time);
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
            ComputeGradients (event );
            ComputeRateofChange((*it));    
            T_RateOfChange_ += clock() - T_begin;         
            if ((*it)->valid() == false) {
                T_begin= clock();
                bool isactive = Schedule(event, model_time);   
                T_Schedule_ += clock() - T_begin;              
                if (isactive) {
                    T_begin= clock();
                    double64 scheduled_time = event->t_schedule();
                    int index = static_cast<int>(event->getNode()->Read(key_EventIndex));
                    Heap_Node* heap_node = new Heap_Node(scheduled_time,index);
                    EventHeap.insert(heap_node);
                    HeapNodeFullList[index] = heap_node;
                    event->inQueue(true);
                    T_InsertToHeap_ += clock() - T_begin; 
                }  
            } 
            event->inPEPStack(false);         
        };    
        T_RateOfChange_ += clock() - T_begin; 
        cout <<"  PEPList size = "<< PEPList.size() <<"  Queue size = "<< EventHeap.size() << endl;

            
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
            dt_PEP = min(dt_PEP, PEP_multiplier_*dt_target);
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
    T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_serial "<<endl;
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



#if defined(_OPENMP)
//advect variable with DES (discrete event simulation), parallel version
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_DES_openmp( double64 model_time, size_t num_threads)
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
    
    ResetCFLMultiplier();
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
                   
            if((*nit)->Status(  this->key_sCO2 ) != DIRICH) {
                (*nit)->Store( key_EventIndex, makeScalar( (*nit)->Status(key_EventIndex), index) );//event index 
                Event<dim>* event = new Event<dim>(*nit);
                PEPList.push_back(event);
                event->inPEPStack(true);
                ComputeGradients (event );
                ComputeRateofChange(event);
                Schedule(event, model_time);
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
        
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = PEPList.begin()+i;
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
            dt_PEP = min(dt_PEP, PEP_multiplier_*dt_target);
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

template class TwoPhaseDESTransport<1U,FlowFunctionsModule1>;
template class TwoPhaseDESTransport<2U,FlowFunctionsModule1>;
template class TwoPhaseDESTransport<3U,FlowFunctionsModule1>;

template class TwoPhaseDESTransport<1U,FlowFunctionsModule2>;
template class TwoPhaseDESTransport<2U,FlowFunctionsModule2>;
template class TwoPhaseDESTransport<3U,FlowFunctionsModule2>;

template class TwoPhaseDESTransport<1U,FlowFunctionsModule3>;
template class TwoPhaseDESTransport<2U,FlowFunctionsModule3>;
template class TwoPhaseDESTransport<3U,FlowFunctionsModule3>;

template class TwoPhaseDESTransport<1U,FlowFunctionsModule4>;
template class TwoPhaseDESTransport<2U,FlowFunctionsModule4>;
template class TwoPhaseDESTransport<3U,FlowFunctionsModule4>;

template class TwoPhaseDESTransport<1U,FlowFunctionsModule5>;
template class TwoPhaseDESTransport<2U,FlowFunctionsModule5>;
template class TwoPhaseDESTransport<3U,FlowFunctionsModule5>;

template class TwoPhaseDESTransport<1U,FlowFunctionsModule6>;
template class TwoPhaseDESTransport<2U,FlowFunctionsModule6>;
template class TwoPhaseDESTransport<3U,FlowFunctionsModule6>;

} // end csmp 



                                                                       
