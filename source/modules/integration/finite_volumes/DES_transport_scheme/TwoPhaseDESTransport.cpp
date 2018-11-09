#include "TwoPhaseDESTransport.h"
#include "Region.h"
#include "Model.h"
#include "DenseMatrix.h"
#include "CSMP_mathUtilities.h"
#include "CO2H2O_FunctionsModule1.h"
#include "equilibrateH2O_CO2_NaCl.h"
#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseDESTransport( Model<dim>& m,
                                                                 const char* target_region,
                                                                 bool with_capillary_spreading,
                                                                 bool with_gravity_forces,
                                                                 bool tensor_k,
                                                                 double64 PEP_multiplier,
                                                                 double64 cfl_multiplier)
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
      relaxing_factor_(10.)
{
    flowfunctions_ = new FLOW_FUNCTIONS<dim> (db_);
    m.InstantiateFiniteVolumes();
    InitializeVariablesAndKeys(m);

// TODO: perhaps only where you have to    calculatePermeabilityProjections(m.Region(target_region));
  
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_sCO2), lower_limit_, upper_limit_ );
    cout<<"TwoPhaseDESTransport constructed"<<endl;
} // end constructor  


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::TwoPhaseDESTransport( Model<dim>& m,
                                                                 const char* target_region,
                                                                 bool with_capillary_spreading,
                                                                 bool with_gravity_forces,
                                                                 bool tensor_k,
                                                                 double64 PEP_multiplier,
                                                                 double64 cfl_multiplier,
                                                                 double64 relaxing_factor)
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
      relaxing_factor_(relaxing_factor)
{
    flowfunctions_ = new FLOW_FUNCTIONS<dim> (db_);
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
    if(!m.Database().IsDefined("saturation gradient")) m.CreateProperty( "saturation gradient", "none", VECTOR, ELEMENT, -1.00E+08 ,1.00E+08);  
    if(!m.Database().IsDefined("pressure gradient")) m.CreateProperty( "pressure gradient", "none", VECTOR, ELEMENT, -1.00E+10 ,1.00E+10);
    if(!m.Database().IsDefined("truncated FV")) m.CreateProperty( "truncated FV", "none", SCALAR, NODE, 1, 0 ,1);
    if(!m.Database().IsDefined("equilibrate")) m.CreateProperty( "equilibrate", "none", SCALAR, NODE, 1, 0 ,1);
    if(!m.Database().IsDefined("update phi and k")) m.CreateProperty( "update phi and k", "none", SCALAR, ELEMENT, 1, 0 ,1);
    
    // model-wide initialisation
    m.Region("Model").InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );    
    m.Region("Model").InputPropertyValue( "truncated FV", makeScalar(PLAIN,0), COMPLETE); 
    m.Region("Model").InputPropertyValue( "equilibrate", makeScalar(PLAIN,0), COMPLETE); 
    m.Region("Model").InputPropertyValue( "update phi and k", makeScalar(PLAIN,0), COMPLETE); 
    
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
    key_equilibrate = INDEX<SCALAR,NODE> ( m.Database().StorageKey("equilibrate") );
    key_UpdatePhiK = INDEX<SCALAR,ELEMENT> ( m.Database().StorageKey("update phi and k") );
    
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
    if ( key_equilibrate.place != NODE || key_equilibrate.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'equilibrate' variable must be SCALAR and placed on NODE"  );  
    if ( key_UpdatePhiK.place != ELEMENT || key_UpdatePhiK.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'update phi and k' variable must be SCALAR and placed on ELEMENT"  );                                                                          
}



template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::initializeFiniteVolumeProperties()
 {
    //FLOW_FUNCTIONS<dim>  flowfunctions(db_);
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref_.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref_.ElementsBegin(); it!=it_end; ++it )
    {
         flowfunctions_->InitialiseBrooksCoreyParameters(*it);
            
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
         double64 sw_shock = flowfunctions_->ShockHeight(*it);
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
             flowfunctions_->InitialiseBrooksCoreyParameters(eptr);
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
             double64 sw_shock = flowfunctions_->ShockHeight(eptr);
             eptr->Store( this->key_ssH2O, makeScalar( eptr->Status( this->key_ssH2O), sw_shock ) );    
             //determine whether FV node is truncated
             if(!gref_.Contains(eptr)) { //parent elment located outside domain
                 halo_stencils_.insert(eptr);
                 truncated_node = true;
             }     
        }
        if(truncated_node) (*nit)->Store( key_cut, makeScalar( (*nit)->Status( key_cut), 1 ) ); 
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
    
    size_t truncated_node = nd->Read(this->key_cut);//check if node is truncated by domain boundary
    
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



//update hysteretic brooks-corey parameters
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::UpdateBCParameters (Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){ 
    //FLOW_FUNCTIONS<dim> flowfunctions(db_);
    const size_t node_parent_elements(nd->Parents());
    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        flowfunctions_->UpdateBrooksCoreyParameters(eptr);
    }
  }
}


//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::ComputeRateofChange( Event<dim>* event )
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
    
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){
    rate_count_++;//recording
    nd->Store(  key_rate, makeScalar( nd->Status( key_rate), nd->Read( key_rate) + 1 ) );
    
    //FLOW_FUNCTIONS<dim> flowfunctions(db_);

    double64 accumulation(0.), flux_balance(0.), outflow(0.);
        
    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> vD, facetNrml;
    const size_t node_parent_elements(nd->Parents());
    
    double64 cfl_multiplier = CFL_multiplier_*relaxing_factor_; //default value
    
    size_t truncated_node = nd->Read(this->key_cut);//check if node is truncated by domain boundary

    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
      Element<dim>* const eptr(nd->Parent(t));
      assert( eptr != NULL );
        
      if(truncated_node == 1 && (this->halo_stencils_.find(eptr) != this->halo_stencils_.end())) { //ignore if parent element located outside domain
        continue;
        
      } else {        
        
        flowfunctions_->InitialiseBrooksCoreyParameters(eptr); 
        
        const size_t pnid(nd->ParentNodeNumber(t));
        
        //eptr->Read( this->key_vt, vD);
        
        //compute total velocity (without gravity)
        VectorVariable<dim> gradP;
        eptr->Read(this->key_gradP, gradP); //pressure gradient
        double64 lambda_t = flowfunctions_->TotalMobility(eptr);
        double64 thickness = eptr->Read(this->key_thi); //thickness
        if(!tensor_k_) { //scalar permeability
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
            const double64 ln_inside_node = flowfunctions_->Mobility_at(eptr,1U,1.0-sn_inside_node); 
            const double64 lw_inside_node = flowfunctions_->Mobility_at(eptr,0U,1.0-sn_inside_node); 
        
            const double64 sn_outside_node = eptr->N(outside_node)->Read(  this->key_sCO2 );
            const double64 sw_outside_node = 1.-sn_outside_node;
            const double64 ln_outside_node = flowfunctions_->Mobility_at(eptr,1U,1.0-sn_outside_node); 
            const double64 lw_outside_node = flowfunctions_->Mobility_at(eptr,0U,1.0-sn_outside_node);   
            
            //compute phase velocities at facet integration point                      
            double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
            
            if( with_gravity_forces_ ){                       
                vn_gravity_component_of_velocity = flowfunctions_->Mobility(eptr, 0U) * flowfunctions_->GravityTerm(eptr) * facetNrml[v];
                vw_gravity_component_of_velocity = flowfunctions_->Mobility(eptr, 1U) * flowfunctions_->GravityTerm(eptr) * facetNrml[v];
            }         
            
            if(with_capillary_spreading_){  
                VectorVariable<dim> grad;
                eptr->Read(this->key_gradSn, grad);
                double64 dsdn = grad.DotProduct(facetNrml);
            
                if(!isnan(dsdn)){
                    vn_capillary_component_of_velocity = -dsdn*flowfunctions_->CapillaryDiffusionMultiplier_Phase(eptr,0U);
                    vw_capillary_component_of_velocity = -dsdn*flowfunctions_->CapillaryDiffusionMultiplier_Phase(eptr,1U);
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
                                       
            if( with_gravity_forces_ )
                gravity_velocity_component = upstream_lambda_overbar * flowfunctions_->GravityTerm(eptr) * facetNrml[v];

            if( with_capillary_spreading_ )
                capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity; 
        
            //update non-wetting flux accumulation   
            double64 fn = sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;                                        
            accumulation += fn; 
            CO2_inflow += fn;

            //determine cfl_multiplier based on non-wetting phase shock saturation
            if (cfl_multiplier != CFL_multiplier_){
                if(vn_at_facet_int_point < 0.0) { //flowing in from outside node (upstream node)
                    double64 sn_shock = 1.0-eptr->Read(this->key_ssH2O); //sn at shock for outside node
                    if (sn_outside_node >= sn_shock) { //upstream node passed shock saturation
                        if (sn_inside_node < sn_shock) {//current node not yet reach shock saturation   
                            cfl_multiplier = CFL_multiplier_;
                        }
                    }
                }
            }
        
            //determine cfl_multiplier based on wetting phase shock saturation
            if (cfl_multiplier != CFL_multiplier_){
                if(vw_at_facet_int_point < 0.0) { //flowing from outside node (upstream node)
                    double64 sw_shock = eptr->Read(this->key_ssH2O); //sw at shock for outside node
                    if (sw_outside_node >= sw_shock) { //upstream node passed shock saturation
                        if (sw_inside_node < sw_shock) {//current node not yet reach shock saturation   
                            cfl_multiplier = CFL_multiplier_;
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
    nd->Read(key_time, array2);
   
    //compute CFL time increment 
    if (outflow < numeric_limits<double64>::epsilon())
        array2.Component(2, numeric_limits<double64>::max());
    else 
        array2.Component(2, nd->Read(  this->key_fvPV ) / outflow);
    
    nd->Store( key_CFL, makeScalar( nd->Status( key_CFL ),cfl_multiplier ) );
    array2.Component(6, cfl_multiplier);    
        
    nd->Store(key_time, array2);
    
    
    // divergence free correction (only when node is located inside domain and flux balance not equal to zero)
    if (truncated_node !=1 && fabs(flux_balance) > numeric_limits<double64>::epsilon()) {
         // compute average fractional flow for the current finite volume
        double64 fn_avg = 0.;
        double64 sw = 1. - nd->Read(this->key_sCO2); //saturation aqueous phase at current node
        for ( size_t t=0U; t<node_parent_elements; t++ )
        {
            Element<dim>* const eptr(nd->Parent(t));
            flowfunctions_->InitialiseBrooksCoreyParameters(eptr);   
            fn_avg += flowfunctions_->f_at(eptr,1U,sw);
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
bool TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::Schedule(Event<dim>* event, double64 t_end)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH);
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    ArrayVariable array;
    nd->Read(key_time, array);

    nd->Store( key_schedule, makeScalar( nd->Status( key_schedule), nd->Read( key_schedule) + 1 ) );
    event->valid(true);
    //compute target change
    double64 dt_CFL = array[2];//CFL time increment
    double64 ChangeRate = nd->Read( key_dsnw);//rate of change
    double64 source = nd->Read(this->key_nQV);
    double64 dC_CFL = dt_CFL*CFL_multiplier_*(-ChangeRate+source);//targe change

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
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::Update_DES(Event<dim>* event, double64 t_clock)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    update_count_++;//recording
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    ArrayVariable array;
    nd->Read(key_time, array);

    double64 ChangeRate = nd->Read( key_dsnw);//variaition rate   
    double64 solution = nd->Read( this->key_sCO2);//old solution
    nd->Store(this->key_sCO2_0, makeScalar( status, solution ));//store old solution
    double64 t_current = array[0]; //current time stamp
    double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
    const double64 source(nd->Read(this->key_nQV));
    new_solution += source * (t_clock - t_current);//add source to new solution
        
    //check new solution value against range and stored it to key_sCO2
    if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    } else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
        if ( new_solution > upper_limit_ ) new_solution = upper_limit_;
        else if ( new_solution < lower_limit_ ) new_solution = lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    }  
    new_solution = nd->Read(this->key_sCO2);//stored new solution
    double64 dsn_cumulative = array[4];
    array.Component(4, dsn_cumulative + (new_solution-solution));//update cumulative change
        
    array.Component(0, t_clock); //current time stamp
    
    nd->Store(key_time, array);

    nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) );
    nd->Store( key_equilibrate, makeScalar( nd->Status(key_equilibrate), 1 ) );
  }
}



//update solution and check it against the specified range (with TDS)
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::Update_TDS(Event<dim>* event, double64 delta_t)
{
  Node<dim>* nd = event->getNode();
  assert( nd  != NULL );   
  assert( nd->Status(  this->key_sCO2 ) != DIRICH); 
  if(nd  != NULL && nd->Status( this->key_sCO2 ) != DIRICH){     
    update_count_++;//recording    
    const VARIABLE_FLAG status(nd->Status( this->key_sCO2 ));
    
    double64 ChangeRate = nd->Read(key_dsnw);//variation rate  
    double64 solution = nd->Read(this->key_sCO2);//old solution
    nd->Store(this->key_sCO2_0, makeScalar( status, solution ));//store old solution
    double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
    const double64 source(nd->Read( this->key_nQV));
    new_solution += source * delta_t;//add source to new solution.
                
    //check new solution value against range and stored it to key_sCO2
    if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) {
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    } else {
        cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
        if ( new_solution > upper_limit_ ) new_solution = upper_limit_;
        else if ( new_solution < lower_limit_ ) new_solution = lower_limit_;
        nd->Store(this->key_sCO2, makeScalar( status, new_solution ));
        nd->Store(this->key_sH2O, makeScalar( status, 1. - new_solution ));
    }
        
    nd->Store( key_update, makeScalar( nd->Status(key_update), nd->Read(key_update) + 1 ) );   
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
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::AdvectVariable_TDS( double64 time_interval)
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
        UpdateBCParameters ((*it));
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
            UpdateBCParameters ((*it));
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

    if (equilibration_) {
        equilibrateFluid();
        updatePorosityandPermeability();
        updatePoreVolume();
    }
}


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::equilibrateFluid()
{
    size_t equilibrate;
    variables::VariableSet_CO2GeoSequestration props(db_);
    for ( auto nit=gref_.NodesBegin(); nit!=gref_.NodesEnd(); ++nit )
    { 
        equilibrate = (*nit)->Read(key_equilibrate);
        if(equilibrate == 1) {
            equilibrateH2O_CO2_NaCl(props, *(*nit) );

            for ( size_t t=0U; t<(*nit)->Parents(); t++ )
            {
                Element<dim>* const eptr((*nit)->Parent(t));
                eptr->Store( key_UpdatePhiK, makeScalar( (*nit)->Status( key_UpdatePhiK), 1 ) );                 
            }
        }
    }  
} 


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::updatePorosityandPermeability()
{
    size_t update;
    variables::VariableSet_CO2GeoSequestration props(db_);
    for ( auto eit=gref_.ElementsBegin(); eit!=gref_.ElementsEnd(); ++eit )
    {
        update = (*eit)->Read(key_UpdatePhiK);
        if(update == 1) {
            double64 phi = (*eit)->Read(this->key_phi);
            double64 permeability = (*eit)->Read(this->key_k);
            double64 phi3(phi * phi * phi);
            double64 t1((1.-phi) * (1.-phi));
            double64 constant = permeability*t1/phi3 ;     
            
            porosityWithSalt( props, *(*eit) );
            
            phi = (*eit)->Read(this->key_phi);
            phi3 = phi * phi * phi;
            t1 = (1.-phi) * (1.-phi);
            permeability = constant*phi3/t1;
            
            (*eit)->Store( this->key_k, makeScalar( (*eit)->Status( this->key_k), permeability ) );
            (*eit)->Store( key_UpdatePhiK, makeScalar( (*eit)->Status( key_UpdatePhiK), 0 ) );
        }
    }
}


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void TwoPhaseDESTransport<dim,FLOW_FUNCTIONS>::updatePoreVolume()
{
    size_t equilibrate;
    for ( auto nit=gref_.NodesBegin(); nit!=gref_.NodesEnd(); ++nit )
    {    
        equilibrate = (*nit)->Read(key_equilibrate);
        if(equilibrate == 1) {   
            double64 pore_volume (0.); 
            for ( size_t t=0U; t<(*nit)->Parents(); t++ )
            {
                Element<dim>* const eptr((*nit)->Parent(t));
                double64 phi = eptr->Read( this->key_phi);
                const double64 thickness = eptr->Read( this->key_thi );
                if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
                const size_t pnid((*nit)->ParentNodeNumber(t));
                const double64 sector_volume = eptr->SectorVolume(pnid);
                pore_volume   += phi * sector_volume;
            }
            (*nit)->Store( this->key_fvPV, makeScalar( (*nit)->Status( this->key_fvPV), pore_volume ) ); 
            (*nit)->Store( key_equilibrate, makeScalar( (*nit)->Status( key_equilibrate), 0 ) ); 
        }         
    }
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
                UpdateBCParameters (event );
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
            UpdateBCParameters (event);
            ComputeRateofChange((*it));    
            T_RateOfChange_ += clock() - T_begin;         
            if ((*it)->valid() == false) {
                T_begin= clock();
                bool isactive = Schedule(event, model_time);   
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
        T_RateOfChange_ += clock() - T_begin; 
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
                UpdateBCParameters (event);
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
        
        for(size_t i = 0U; i < PEPList_size; ++i)
        {
            auto it = PEPList.begin()+i;
            Event<dim>* event = *it;   
            UpdateBCParameters (event );
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



template class TwoPhaseDESTransport<1U,CO2H2O_FunctionsModule0>;
template class TwoPhaseDESTransport<2U,CO2H2O_FunctionsModule0>;
template class TwoPhaseDESTransport<3U,CO2H2O_FunctionsModule0>;

template class TwoPhaseDESTransport<1U,CO2H2O_FunctionsModule1>;
template class TwoPhaseDESTransport<2U,CO2H2O_FunctionsModule1>;
template class TwoPhaseDESTransport<3U,CO2H2O_FunctionsModule1>;

template class TwoPhaseDESTransport<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseDESTransport<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseDESTransport<3U,CO2H2O_FunctionsModule2>;

} // end csmp 



                                                                       
