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
TwoPhaseDESTransport<dim>::TwoPhaseDESTransport( Model<dim>& m, const char* target_region, bool with_capillary_spreading, bool with_gravity_forces)
    : gref_(m.Region(target_region)),
      with_capillary_spreading_(with_capillary_spreading),
      with_gravity_forces_(with_gravity_forces),
      upper_limit_(1.), lower_limit_(0.), rate_count_(0U), update_count_(0U),
      T_RateOfChange_(0.), T_Schedule_(0.), T_SortQueue_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromQueue_(0.), T_AdvectVariable_(0.)
{
    m.InstantiateFiniteVolumes();
    initializeVariablsAndKeys(m);
    initializeFiniteVolumeProperties();
    //create events for all nodes and add them to PEPStack and EntireQueue
    int index = 0;
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        (*nit)->Store( EventIndex_key, makeScalar( (*nit)->Status(EventIndex_key), index) );//event index        
        (*nit)->Store( update_key, makeScalar( (*nit)->Status( update_key), 0 ) ); //update count
        (*nit)->Store( rate_key, makeScalar( (*nit)->Status( rate_key), 0 ) ); //changerate count
        (*nit)->Store( schedule_key, makeScalar( (*nit)->Status( schedule_key), 0 ) ); //schedule count
        (*nit)->Store( synchronize_key, makeScalar( (*nit)->Status( synchronize_key), 0 ) ); //synchronize count   
        
        ArrayVariable arrayVariable( 6, 0., PLAIN );         
        (*nit)->Store( time_key, arrayVariable );  
                   
        Event<dim>* event = new Event<dim>(*nit);
        event->inPEPStack(true);
        event->valid(false);
        PEPStack.push_back(event); //add event to PEPStack
        EntireQueue.push_back(event); //add event to EntireQueue
        index++;
    }
         
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(sn_key), lower_limit_, upper_limit_ );
    cout<<"events created for all nodes and added to PEPStack - size = "<<PEPStack.size()<<endl;
    cout<<"TwoPhaseDESTransport constructed"<<endl;
} // end constructor  



template<size_t dim>
void TwoPhaseDESTransport<dim>::initializeVariablsAndKeys(Model<dim>& m)
{
    //creating new variables if not defined yet from input file
    if(!m.Database().IsDefined("event index")) m.CreateProperty( "event index", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("update count")) m.CreateProperty( "update count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!m.Database().IsDefined("rate count")) m.CreateProperty( "rate count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("schedule count")) m.CreateProperty( "schedule count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("synchronize count")) m.CreateProperty( "synchronize count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!m.Database().IsDefined("variation rate")) m.CreateProperty( "variation rate", "s-1", SCALAR, NODE, 1, -1.00E+08 ,1.00E+08); 
    
    // model-wide initialisation
    //m.Region("Model").InputPropertyValue( "new concentration", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "finite volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );    
    
    //assigning keys 
    phi_key = m.Database().StorageKey("porosity");
    vD_key  = m.Database().StorageKey("velocity");
    fv_key  = m.Database().StorageKey("finite volume");
    PV_key  = m.Database().StorageKey("FV pore volume");
    fA_key  = m.Database().StorageKey("facet area");
    ff_key  = m.Database().StorageKey("facet flux");
    fn_key  = m.Database().StorageKey("facet normal");
    fb_key  = m.Database().StorageKey("flux balance");
    nsrc_key = m.Database().StorageKey("nodal fluid volume source");

    sw_key = m.Database().StorageKey("saturation water"); // water saturation (node)
    sn_key = m.Database().StorageKey("saturation oil"); // non-wetting phase saturation (node)  
    dsn_key = m.Database().StorageKey("variation rate"); // variation rate of non-wetting phase saturation (node)      
    muw_key = m.Database().StorageKey("viscosity water"); // fluid viscosity (node)
    mun_key = m.Database().StorageKey("viscosity oil"); // oil viscosity (node)
    rhw_key = m.Database().StorageKey("density water"); // water density (node)
    rhn_key = m.Database().StorageKey("density oil"); // oil density (node)
    swr_key = m.Database().StorageKey("residual saturation wetting phase"); // irreducible saturation wetting phase (element)
    snr_key = m.Database().StorageKey("residual saturation non-wetting phase"); // residual saturation non-wetting phase (element)
    bcp_key = m.Database().StorageKey("brooks corey parameter"); // Brooks-Corey parameter (element) 
    k_key = m.Database().StorageKey("permeability"); // permeability (element)    
    pd_key = m.Database().StorageKey("entry pressure"); // entry pressure (element)      
    
    EventIndex_key = m.Database().StorageKey("event index");     
    update_key = m.Database().StorageKey("update count");
    rate_key = m.Database().StorageKey("rate count");
    schedule_key = m.Database().StorageKey("schedule count");
    synchronize_key = m.Database().StorageKey("synchronize count");   
    time_key  = m.Database().StorageKey("non-wetting phase timing array");      
                                                    
    if ( phi_key.place != ELEMENT || phi_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( vD_key.place != ELEMENT || vD_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( fv_key.place != NODE || fv_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'finite volume' variable must be SCALAR and placed on NODE"  );
    if ( PV_key.place != NODE || PV_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( fA_key.place != FACET_INTEGRATION_POINT || fA_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( ff_key.place != FACET_INTEGRATION_POINT || ff_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'facet flux' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( fn_key.place != FACET_INTEGRATION_POINT || fn_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( fb_key.place != NODE || fb_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );
    if ( nsrc_key.place != NODE || nsrc_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( sw_key.place != NODE || sw_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'saturation water' variable must be SCALAR and placed on NODE"  );
    if ( sn_key.place != NODE || sn_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'saturation oil' variable must be SCALAR and placed on NODE"  );    
    if ( dsn_key.place != NODE || dsn_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'variation rate' variable must be SCALAR and placed on NODE"  );
    if ( muw_key.place != NODE || muw_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'viscosity water' variable must be SCALAR and placed on NODE"  );
    if ( mun_key.place != NODE|| mun_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'viscosity oil' variable must be SCALAR and placed on NODE"  );
    if ( swr_key.place != ELEMENT || swr_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'residual saturation wetting phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( snr_key.place != ELEMENT || snr_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'residual saturation non-wetting phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( bcp_key.place != ELEMENT || bcp_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'brooks corey parameter' variable must be SCALAR and placed on ELEMENT"  );          
    if ( k_key.place != ELEMENT || k_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( pd_key.place != ELEMENT || pd_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'entry pressure' variable must be SCALAR and placed on ELEMENT"  );                               
    if ( EventIndex_key.place != NODE || EventIndex_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'event index' variable must be SCALAR and placed on NODE"  );       
    if ( update_key.place != NODE || update_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( rate_key.place != NODE || rate_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( schedule_key.place != NODE || schedule_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( synchronize_key.place != NODE || synchronize_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  );              
    if ( sn_key.place != NODE || sn_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'saturation oil' variable must be SCALAR and placed on NODE"  );
    if ( dsn_key.place != NODE || dsn_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'variation rate' variable must be SCALAR and placed on NODE"  );     
    if ( time_key.place != NODE || time_key.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseDESTransport::initializeKeys:",
        "The 'non-wetting phase timing array' variable must be ARRAY and placed on NODE"  );                
}


template<size_t dim>
void TwoPhaseDESTransport<dim>::initializeFiniteVolumeProperties()
 {
    Point<dim>           nrml;
    VectorVariable<dim>  fnrml, vt;

    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref_.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref_.ElementsBegin(); it!=it_end; ++it )
      {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // element-based total velocity
         (*it)->Read( vD_key, vt );

         // 1. computing sector pore volumes
         // --------------------------------
         const double64 phi = (*it)->Read( phi_key );
         for ( size_t i=0U; i<sectors; ++i ) {
              // sector pore volume
              const double64 sector_volume = (*it)->SectorVolume(i);       
              // sector volume is added to  pore volume of FV's containing this sector
              double64 finite_volume = (*it)->N(i)->Read( fv_key );
              double64 pore_volume   = (*it)->N(i)->Read( PV_key );
              // sector volume from FV traits
              finite_volume += sector_volume;
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( fv_key, makeScalar(PLAIN,finite_volume) );
              (*it)->N(i)->Store( PV_key, makeScalar(PLAIN,pore_volume) );
           }

         // 2. computing facet normals and areas
         // ------------------------------------
         for ( size_t j=0U; j<facets; ++j ) {
              // computing facet areas
              const double64 facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, fA_key, makeScalar( PLAIN, facet_area ) );
              // computing facet normals
              nrml = (*it)->FacetNormal(j);
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, fn_key, fnrml );
           }
      }

   // 4. initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // -------------------------------------------------------------------------------------------------------------------
   // (here the pore volumes do not include the sectors outside the region)
   const typename vector<Node<dim>*>::iterator nit_end(gref_.NodesEnd());
   for ( typename vector<Node<dim>*>::iterator nit=gref_.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const size_t parent_elements((*nit)->Parents());      
        for ( size_t i=0U; i<parent_elements; ++i ) {
             Element<dim>* const eptr = (*nit)->Parent(i);
             // ---------------------------------
             // computing facet normals and areas
             // ---------------------------------
             const size_t facets(eptr->Facets());
             for ( size_t j=0U; j<facets; ++j ) {
                  // computing facet areas
                  const double64 facet_area = eptr->FacetArea(j);
                  eptr->Store( j, 0U, fA_key, makeScalar( PLAIN, facet_area ) );
                  // computing facet normals
                  nrml = eptr->FacetNormal(j);
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, fn_key, fnrml );
               }
          }
     }
 } // end initializeFiniteVolumeProperties



//Compute the rate of change of non-wetting phase in a node/FV
template<size_t dim>
void TwoPhaseDESTransport<dim>::ComputeRateofChange( Event<dim>* event )
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );

    //ignore DIRICH node
    if(nd->Status(  sn_key ) == DIRICH) {   
        nd->Store(  fb_key, makeScalar( nd->Status(  fb_key ), 0. ) );//flux balance
        ArrayVariable array;
        nd->Read(time_key, array);
        array.Component(2, numeric_limits<double64>::max()); //CFL time increment
        nd->Store(time_key, array);    
        return;
    };

    rate_count_++;//recording
    nd->Store(  rate_key, makeScalar( nd->Status( rate_key), nd->Read( rate_key) + 1 ) );

    double64 accumulation(0.), flux_balance(0.), outflow(0.);;
    const size_t node_parent_elements(nd->Parents());
    VectorVariable<dim>  nrml, vD;

    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        const size_t pnid(nd->ParentNodeNumber(t));
        eptr->Read(  vD_key, vD);
        
        //check if vD is nan
        bool isNAN(false);
        for (size_t n=0U;n<dim;++n) {
            if (isnan(vD(n))){
                isNAN=true;
                break;
            };
        };

        if (!isNAN){
        
            const size_t v( (dim==1u) ? 0u : 1u );
            //read in variables placed on element
            const double64 swr = eptr->Read(  swr_key );//irreducible saturation wetting phase  
            const double64 snr = eptr->Read(  snr_key );//irreducible saturation non-wetting phase    
            const double64 lambda = eptr->Read(  bcp_key );//Brooks-Corey parameter   
            const double64 k = eptr->Read(  k_key );//permeability
            const double64 pd = eptr->Read(  pd_key );//entry pressure
            
            std::vector<double64> dsdn_(dim);
            DenseMatrix<DM_MIN> DN_;
            if(with_capillary_spreading_){                
                fill( dsdn_.begin(), dsdn_.end(), 0. );
                eptr->dN_AtBaryCenter( DN_ );
                for ( size_t j=0U; j<eptr->Nodes(); j++ ) {
                    const double64 sn = eptr->N(j)->Read( sn_key);
                    for ( size_t k=0U; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
                }
            }
             
            const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
            for ( size_t i=0U; i<sector_facets; i++ )
            {
                const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                const size_t inside_node(eptr->FV()->InsideNode(iFacet));
                const size_t outside_node(eptr->FV()->OutsideNode(iFacet));

                eptr->Read( iFacet, 0U,  fn_key, nrml );
                const double64  vD_n = vD.DotProduct(nrml);
                const double64  facetArea = eptr->Read( iFacet, 0U,  fA_key );       
                
                const size_t upstream_node = (vD_n < 0.) ? outside_node : inside_node;// finding the upstream node
                //read in variables placed on node
                const double64 sn = eptr->N(upstream_node)->Read(  sn_key );//non-wetting phase saturation                                         
                const double64 sw = 1.-sn;//wetting phase saturation
                const double64 muw = eptr->N(upstream_node)->Read(  muw_key );//viscosity wetting phase
                const double64 mun = eptr->N(upstream_node)->Read(  mun_key );//viscosity non-wetting phase
                const double64 rhw = eptr->N(upstream_node)->Read(  rhw_key );//density wetting phase
                const double64 rhn = eptr->N(upstream_node)->Read(  rhn_key );//density non-wetting phase
                                                                                   
                const double64 sign = ( pnid == inside_node ) ? 1. : -1.;
                //compute facet fluid flux
                double64 flux = sign * vD_n * facetArea;
                //update flux balance
                flux_balance += flux;
                //update outflow
                if ( flux > 0. ) outflow += flux;  
               
                double64 viscous_velocity_component(0.0), capillary_velocity_component(0.0),gravity_velocity_component(0.0);
                if( !with_gravity_forces_ && !with_capillary_spreading_ ){ //viscous effect only               
                    //compute fractional flow for non-wetting phase
                    double64 f_n = ComputeNonWettingFractionalFlow( sw, swr, snr, muw, mun, lambda);               
                    //compute viscous effect
                    viscous_velocity_component = vD_n* f_n;
                    
                } else { // with gravity or capillary effects
                    //compute gravity and capillary effects
                    double64 vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
                    double64 vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );
                    if( with_gravity_forces_ ){                       
                        vn_gravity_component_of_velocity = ComputeMobilityPhase(2U,sw,swr,snr,muw,mun,lambda)*ComputeGravityTerm(k, rhw, rhn) * nrml(v);
                        vw_gravity_component_of_velocity = ComputeMobilityPhase(1U,sw,swr,snr,muw,mun,lambda)*ComputeGravityTerm(k, rhw, rhn) * nrml(v);
                    }                    
                    
                    if(with_capillary_spreading_){
                        //compute dsdn       
                        double64  dsdn = dsdn_[0]*nrml[0];
                        if (dim != 1) dsdn += dsdn_[1]*nrml[1];
                        if (dim == 3) dsdn += dsdn_[2]*nrml[2]; 
                        double64 dpcdn = -dsdn*Compute_dpcds(sw,swr,snr,muw,mun,lambda,pd); 
                        vn_capillary_component_of_velocity = ComputeMobilityPhase(2U,sw,swr,snr,muw,mun,lambda) * k * dpcdn;
                        vw_capillary_component_of_velocity = ComputeMobilityPhase(1U,sw,swr,snr,muw,mun,lambda) * k * dpcdn;
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

                    double64 upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0)/*, upstream_sn(0.0), upstream_sw(0.0)*/;
                    //read in variables placed on inside and outside nodes
                    const double64 sn_inside_node = eptr->N(inside_node)->Read(  sn_key );//non-wetting phase saturation                                         
                    const double64 sw_inside_node = 1.-sn_inside_node;//wetting phase saturation
                    const double64 muw_in = eptr->N(inside_node)->Read(  muw_key );//viscosity wetting phase
                    const double64 mun_in = eptr->N(inside_node)->Read(  mun_key );//viscosity non-wetting phase     
                    
                    const double64 sn_outside_node = eptr->N(outside_node)->Read(  sn_key );//non-wetting phase saturation                                         
                    const double64 sw_outside_node = 1.-sn_outside_node;//wetting phase saturation
                    const double64 muw_out = eptr->N(outside_node)->Read(  muw_key );//viscosity wetting phase
                    const double64 mun_out = eptr->N(outside_node)->Read(  mun_key );//viscosity non-wetting phase   
                    
                    //computes mobilities on inside and outside nodes
                    const double64 ln_inside_node  = ComputeMobilityPhase(2U,sw_inside_node,swr,snr,muw_in,mun_in,lambda);
                    const double64 lw_inside_node  = ComputeMobilityPhase(1U,sw_inside_node,swr,snr,muw_in,mun_in,lambda);
                    
                    const double64 ln_outside_node  = ComputeMobilityPhase(2U,sw_outside_node,swr,snr,muw_out,mun_out,lambda);
                    const double64 lw_outside_node  = ComputeMobilityPhase(1U,sw_outside_node,swr,snr,muw_out,mun_out,lambda);
                    
                    double64 zero(0.0);
                    if((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point>zero)){
                        upstream_mobility_n=ln_inside_node;
                        upstream_mobility_w=lw_inside_node;
                        //upstream_sn = sn_inside_node;
                        //upstream_sw = 1.0 - sn_inside_node;
                    }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero)){
                        upstream_mobility_n=ln_inside_node;
                        upstream_mobility_w=lw_outside_node;
                        //upstream_sn = sn_inside_node;
                        //upstream_sw = 1.0 - sn_outside_node;
                    }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point>zero)){
                        upstream_mobility_n=ln_outside_node;
                        upstream_mobility_w=lw_inside_node;
                        //upstream_sn = sn_outside_node;
                        //upstream_sw = 1.0 - sn_inside_node;
                    }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point<zero)){
                        upstream_mobility_n=ln_outside_node;
                        upstream_mobility_w=lw_outside_node;
                        //upstream_sn = sn_outside_node;
                        //upstream_sw = 1.0 - sn_outside_node;
                    }else{
                        upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                        upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);
                        //upstream_sn = 0.5*(sn_inside_node + sn_outside_node);
                        //upstream_sw = 1.0 - upstream_sn;
                    }                    
                    
                    total_mobility=upstream_mobility_n+upstream_mobility_w;
                    double64 upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
                    double64 upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);
                    
                    viscous_velocity_component = vD_n * upstream_fn;
                                       
                    if( with_gravity_forces_ )
                        gravity_velocity_component = upstream_lambda_overbar * ComputeGravityTerm(k, rhw, rhn) * nrml(v);

                    if( with_capillary_spreading_ )
                        capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity; 
                }                
                                                        
                accumulation += sign * ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component) * facetArea;                  
            }
        }
    }
    
    //compute CFL time increment 
    ArrayVariable array2;
    nd->Read(time_key, array2);
    if (outflow < numeric_limits<double64>::epsilon())
        array2.Component(2, numeric_limits<double64>::max());
    else 
        array2.Component(2, nd->Read(  PV_key ) / outflow);  
    nd->Store(time_key, array2);
    
    // divergence free correction
    // compute average fractional flow for the current finite volume
    double64 fn_avg = static_cast<double64>(0.);
    if (fabs(flux_balance) > numeric_limits<double64>::epsilon()) {
        //read in nodal variables
        const double64 sn = nd->Read(  sn_key );//non-wetting phase saturation                                         
        const double64 sw = 1.-sn;//wetting phase saturation
        const double64 muw = nd->Read(  muw_key );//viscosity wetting phase
        const double64 mun = nd->Read(  mun_key );//viscosity non-wetting phase                      
        // for each finite volume, f is evaluated on a sector by sector basis 
        for ( size_t t=0U; t<node_parent_elements; t++ ) {
            Element<dim>* const eptr(nd->Parent(t));
            assert( eptr != NULL );
            //read in element variables
            const double64 swr = eptr->Read(  swr_key );//irreducible saturation wetting phase  
            const double64 snr = eptr->Read(  snr_key );//irreducible saturation non-wetting phase    
            const double64 lambda = eptr->Read(  bcp_key );//Brooks-Corey parameter                  
            // weighted for each specific FV sector
            fn_avg += ComputeNonWettingFractionalFlow( sw, swr, snr, muw, mun, lambda); 
        }
        fn_avg /= node_parent_elements;  
        accumulation -= fn_avg*flux_balance;
    } 
       
    double64 PV = nd->Read(  PV_key );
    //store variation rate
    nd->Store(  dsn_key, makeScalar( nd->Status( dsn_key ), accumulation/PV ) );    
} 


//Computes the mobility of phase
template<size_t dim>
double64 TwoPhaseDESTransport<dim>::ComputeMobilityPhase( size_t phase, double64 sw, double64 swr, double64 snr, double64 muw, double64 mun, double64 lambda)
 {
    //computes effective saturation for wetting phase
    double64 seff = std::min( std::max( (sw - swr) / (1. - swr - snr), 0. ), 1. );
    
    //computes relative k for non-wetting and wetting phases (based on brooks corey model)
    double64 krn, krw;
    if ( lambda == static_cast<double64>(0.) ) { //  switch to linear relperm model if lambda = 0
        krw = seff;
        krn = 1. - seff;
    } else {
        krw = std::pow( seff, 2./lambda + 3.0 );
        double64  seffn = 1.-seff;
        krn = (seffn * seffn) * (1. - pow( seff, 2./lambda + 1.0) );
    }    
        
    if ( phase == 1U ) return krw / muw;
    return krn / mun; 
 } 
 
 
//Computes dpcds
template<size_t dim>
double64 TwoPhaseDESTransport<dim>::Compute_dpcds(double64 sw, double64 swr, double64 snr, double64 muw, double64 mun, double64 lambda, double64 pd)
{
    const double64 dSedSw( 1.0/ (1.0 - swr - snr) );
      
    //computes effective saturation for wetting phase
    double64 seff = std::min( std::max( (sw - swr) / (1. - swr - snr), 0. ), 1. );  
    
    double64 dpcds;
    double64 h = 0.00001;
    if( seff < 0.+h )
        dpcds = ( Compute_pc(lambda, seff + h, pd, swr, snr ) - Compute_pc(lambda, seff, pd, swr, snr ) ) / h * dSedSw;
    else if( seff > 1.-h )
        dpcds = ( Compute_pc(lambda, seff, pd, swr, snr ) - Compute_pc(lambda, seff - h, pd, swr, snr ) ) / h * dSedSw;
    else
        dpcds = ( Compute_pc(lambda, seff + h, pd, swr, snr ) - Compute_pc(lambda, seff - h, pd, swr, snr ) )/ (2.0*h) * dSedSw;

  return dpcds;

}


//Computes pc
template<size_t dim>
double64 TwoPhaseDESTransport<dim>::Compute_pc( double64 lambda, double64 se, double64 entry_pressure, double64 swr, double64 snr)
{
    // linear relperm model
    double64 MAX_CAPILLARY_PRESSURE_= 4e7;
    double64 MAX_CAPILLARY_PRESSURE_SLOPE_ = 1e7;
    if ( lambda == 0. ){
        if ( se <= 0. )
            return MAX_CAPILLARY_PRESSURE_;

        if ( se >= 1. )
            return entry_pressure;
		
        if (entry_pressure == MAX_CAPILLARY_PRESSURE_)
	    return entry_pressure;

        return entry_pressure + ( 1. - se ) * ( MAX_CAPILLARY_PRESSURE_ - entry_pressure );
   }

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure == 0. )
       return 0.;

   const double64 se_mult( 1.0/ (1.0 - swr - snr ) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = entry_pressure * pow( ( entry_pressure / ( lambda * MAX_CAPILLARY_PRESSURE_SLOPE_/se_mult ) ),
                                           ( -1. / ( 1. + lambda ) ) );

   if ( se <= std::pow( pcmax / entry_pressure, -lambda ) )
   {
      // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
       const double64 Se_min =  pow( ( entry_pressure / ( lambda * MAX_CAPILLARY_PRESSURE_SLOPE_/se_mult ) ),
                                ( lambda / ( 1. + lambda ) ) );
       // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
       return pcmax + ( Se_min - se ) * MAX_CAPILLARY_PRESSURE_SLOPE_/se_mult;
   }

   return entry_pressure * std::pow( se, -1. / lambda );
}



//Computes the fractional flow of the non-wetting phase (gravitational and capillary effects are not included)
template<size_t dim>
double64 TwoPhaseDESTransport<dim>::ComputeNonWettingFractionalFlow( double64 sw, double64 swr, double64 snr, double64 muw, double64 mun, double64 lambda)
 {
    //computes effective saturation for wetting phase
    double64 seff = std::min( std::max( (sw - swr) / (1. - swr - snr), 0. ), 1. );
    
    //computes relative k for non-wetting and wetting phases (based on brooks corey model)
    double64 krn, krw;
    if ( lambda == static_cast<double64>(0.) ) { //  switch to linear relperm model if lambda = 0
        krw = seff;
        krn = 1. - seff;
    } else {
        krw = std::pow( seff, 2./lambda + 3.0 );
        double64  seffn = 1.-seff;
        krn = (seffn * seffn) * (1. - pow( seff, 2./lambda + 1.0) );
    }    
    
    //computes mobility for non-wetting phase and the total mobility
    double64 mobility = krn / mun;
    double64 total_mobility = krn / mun + krw / muw;
    
    //computes fractional flow for non-wetting phase
    double64 f_n = mobility/total_mobility;
    return f_n;    
 } // ComputeNonWettingFractionalFlow


//Computes the gravity term
template<size_t dim>
double64 TwoPhaseDESTransport<dim>::ComputeGravityTerm( double64 k, double64 rhw, double64 rhn)
 {
  // note that the projected gravity acts opposite the y-axis
  const double64 k_g_drho = k * -9.8066 * (rhw - rhn);

    // economizing the calculation
    if ( std::fabs(k_g_drho) < 1.0e-17 ) return static_cast<double64>(0.);

    // else compute result using G saturation derivative
    return k_g_drho;
}



//schedule an event associated with a node/FV
template<size_t dim>
bool TwoPhaseDESTransport<dim>::Schedule(Event<dim>* event, double64 t_end, double64 cfl_multiplier)
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );
    ArrayVariable array;
    nd->Read(time_key, array);
    //ignore DIRICH node
    if(nd->Status(  sn_key ) == DIRICH) {
        return false;
    } else {
        nd->Store( schedule_key, makeScalar( nd->Status( schedule_key), nd->Read( schedule_key) + 1 ) );
        event->valid(true);
        //compute target change
        double64 CFL = array[2];//CFL number
        double64 ChangeRate = nd->Read( dsn_key);//rate of change
  
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
            nd->Store(time_key, array);
            return false;
        } else {
            event->t_schedule(t_current + dt_target);
            array.Component(1, t_current + dt_target);//schedule time stamp
            nd->Store(time_key,array);
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
    nd->Read(time_key, array);
    const VARIABLE_FLAG status(nd->Status( sn_key ));
    if ( status != DIRICH )
    {    
        double64 ChangeRate = nd->Read( dsn_key);//variaition rate   
        double64 solution = nd->Read( sn_key);//old solution
        double64 t_current = array[0]; //current time stamp
        double64 new_solution = solution - (t_clock - t_current) * ChangeRate;//compute new solution
        const double64 source(nd->Read(nsrc_key));
        new_solution += source * (t_clock - t_current);//add source to new solution
        
        //check new solution value against range and stored it to sn_key
        if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(sn_key, makeScalar( status, new_solution ));
        else {
            cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( new_solution > upper_limit_ ) nd->Store( sn_key, makeScalar( status, upper_limit_ ) );
            else if ( new_solution < lower_limit_ ) nd->Store( sn_key, makeScalar( status, lower_limit_ ) );
        }
        
        new_solution = nd->Read(sn_key);//stored new solution
        double64 dsn_cumulative = array[4];
        array.Component(4, dsn_cumulative + (new_solution-solution));//update cumulative change
        
        array.Component(0, t_clock); //current time stamp
        nd->Store(time_key, array);

        nd->Store( update_key, makeScalar( nd->Status(update_key), nd->Read(update_key) + 1 ) );
    }
}



//update solution and check it against the specified range (with TDS)
template<size_t dim>
void TwoPhaseDESTransport<dim>::Update_TDS(Event<dim>* event, double64 delta_t)
{
    update_count_++;//recording
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL );    
    const VARIABLE_FLAG status(nd->Status( sn_key ));
    if ( status != DIRICH )
    {    
        double64 ChangeRate = nd->Read(dsn_key);//variation rate  
        double64 solution = nd->Read(sn_key);//old solution
        double64 new_solution = solution - delta_t * ChangeRate;//compute new solution
        const double64 source(nd->Read( nsrc_key));
        new_solution += source * delta_t;//add source to new solution
        
        //check new solution value against range and stored it to sn_key
        if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(sn_key, makeScalar( status, new_solution ));
        else {
            cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( new_solution > upper_limit_ ) nd->Store( sn_key, makeScalar( status, upper_limit_ ) );
            else if ( new_solution < lower_limit_ ) nd->Store( sn_key, makeScalar( status, lower_limit_ ) );
        }
        
        nd->Store( update_key, makeScalar( nd->Status(update_key), nd->Read(update_key) + 1 ) );    
    }
}



//Synchronize neighbor nodes/FVs
template<size_t dim>
void TwoPhaseDESTransport<dim>::Synchronize(Event<dim>* event,double64 t_clock)
{
    Node<dim>* nd = event->getNode();
    assert( nd  != NULL ); 
    nd->Store( synchronize_key, makeScalar( nd->Status(synchronize_key), nd->Read(synchronize_key) + 1 ) );
    event->valid(false);
    ArrayVariable array;
    nd->Read(time_key, array);
    array.Component(4, 0.); //reset cumulative change of solution
    nd->Store(time_key, array);
    for ( size_t n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* neighbor_node = nd->Neighbor(n);
        int index = neighbor_node->Read(EventIndex_key);
        Event<dim>* neighbor_event = EntireQueue[index];  
        assert( neighbor_event  != NULL ); 
        if (neighbor_event->inPEPStack() == false) {
            PEPStack.push_back(neighbor_event);
            neighbor_event->inPEPStack(true);
            Update_DES(neighbor_event,t_clock);
            ArrayVariable neighbor_array;
            neighbor_node->Read(time_key, neighbor_array); 
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
        (*it)->getNode()->Read(time_key, array);
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
            (*it)->getNode()->Read(time_key, array2);
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
            (*begin)->getNode()->Read(time_key, begin_array);
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
            (*top)->getNode()->Read(time_key, array);
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



                                                                       
