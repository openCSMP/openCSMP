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
  : gref_(m.Region(target_region)),
    upper_limit_(1.), lower_limit_(0.), rate_count_(0U), update_count_(0U)
{
    m.InstantiateFiniteVolumes();
    initializeKeys(m);
    initializeFiniteVolumeProperties();

    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(C0_key), lower_limit_, upper_limit_ );
    
    //add all nodes to PEPStack
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    { 
        ComputeFluxBalanceAndCFL((*nit)); //compute flux balance and CFL number
        PEPStack.push_back((*nit)); //add node to PEPStack
        (*nit)->inPEPStack_ = true;
        (*nit)->valid_ = false;
        (*nit)->Store(  update_key, makeScalar( (*nit)->Status( update_key), 0 ) ); //update count
        (*nit)->Store(  rate_key, makeScalar( (*nit)->Status( rate_key), 0 ) ); //changerate count
        (*nit)->Store(  schedule_key, makeScalar( (*nit)->Status( schedule_key), 0 ) ); //schedule count
        (*nit)->Store(  synchronize_key, makeScalar( (*nit)->Status( synchronize_key), 0 ) ); //synchronize count
        
    }; 
    cout<<"all nodes added to PEPStack - size = "<<PEPStack.size()<<endl;
    cout<<"DESTransport constructed"<<endl;
}


template<size_t dim>
void DESTransport<dim>::initializeKeys(Model<dim>& m)
{
    // model-wide initialisation: results will be accumulated into this variable
    m.Region("Model").InputPropertyValue( "new concentration", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "finite volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    m.Region("Model").InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );
    
    update_key = m.Database().StorageKey("update count");
    rate_key = m.Database().StorageKey("rate count");
    schedule_key = m.Database().StorageKey("schedule count");
    synchronize_key = m.Database().StorageKey("synchronize count");   
    phi_key = m.Database().StorageKey("porosity");
    vD_key  = m.Database().StorageKey("velocity");
    fv_key  = m.Database().StorageKey("finite volume");
    PV_key  = m.Database().StorageKey("FV pore volume");
    sv_key  = m.Database().StorageKey("sector volume");
    spv_key = m.Database().StorageKey("sector pore volume");
    fA_key  = m.Database().StorageKey("facet area");
    ff_key  = m.Database().StorageKey("facet flux");
    fn_key  = m.Database().StorageKey("facet normal");
    fb_key  = m.Database().StorageKey("flux balance");
    nsrc_key = m.Database().StorageKey("nodal fluid volume source");
    C0_key = m.Database().StorageKey("concentration");
    C1_key = m.Database().StorageKey("new concentration");

    if ( update_key.place != NODE || update_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( rate_key.place != NODE || rate_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( schedule_key.place != NODE || schedule_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( synchronize_key.place != NODE || synchronize_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  );                                
    if ( phi_key.place != ELEMENT || phi_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( vD_key.place != ELEMENT || vD_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( fv_key.place != NODE || fv_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'finite volume' variable must be SCALAR and placed on NODE"  );
    if ( PV_key.place != NODE || PV_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( sv_key.place != SECTOR_INTEGRATION_POINT || sv_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'sector volume' variable must be SCALAR and placed on SECTOR_INTEGRATION_POINT"  );
    if ( spv_key.place != SECTOR_INTEGRATION_POINT || spv_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'sector pore volume' variable must be SCALAR and placed on SECTOR_INTEGRATION_POINT"  );
    if ( fA_key.place != FACET_INTEGRATION_POINT || fA_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( ff_key.place != FACET_INTEGRATION_POINT || ff_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'facet flux' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( fn_key.place != FACET_INTEGRATION_POINT || fn_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( fb_key.place != NODE || fb_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );
    if ( nsrc_key.place != NODE || nsrc_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( C0_key.place != NODE || C0_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'concentration' variable must be SCALAR and placed on NODE"  );
    if ( C1_key.place != NODE || C1_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DESTransport::initializeKeys:",
        "The 'new concentration' variable must be SCALAR and placed on NODE"  );
}


template<size_t dim>
void DESTransport<dim>::initializeFiniteVolumeProperties()
 {
    const bool initialize_flux(true);
    
    Point<dim>           nrml;
    VectorVariable<dim>  fnrml, vt;

    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref_.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref_.ElementsBegin(); it!=it_end; ++it )
      {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // element-based total velocity
         if ( initialize_flux ) (*it)->Read( vD_key, vt );

         // 1. computing sector pore volumes
         // --------------------------------
         const double64 phi = (*it)->Read( phi_key );
         for ( size_t i=0U; i<sectors; ++i ) {
              // sector pore volume
              const double64 sector_volume = (*it)->SectorVolume(i);       
              (*it)->Store( i, 0U, sv_key, makeScalar( PLAIN, sector_volume ) );
              (*it)->Store( i, 0U, spv_key, makeScalar( PLAIN, phi * sector_volume ) );
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

              // 3. computing total facet fluxes and flux balance
              // ------------------------------------------------
              if ( initialize_flux ) {
                   double64 facet_flux(nrml[0] * vt[0]);
                   if ( dim != 1U ) facet_flux += nrml[1] * vt[1];
                   if ( dim == 3U ) facet_flux += nrml[2] * vt[2];
                   facet_flux *= facet_area;
                   (*it)->Store( j, 0U, ff_key, makeScalar((*it)->Status( j, 0U, ff_key),facet_flux) );
                }
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
             // ---------------------------------------
             // computing sector volumes & pore volumes
             // ---------------------------------------
             const double64 porosity = eptr->Read( phi_key );
             const size_t sectors(eptr->Sectors());
             for ( size_t j=0U; j<sectors; ++j ) {
                  const double64 sector_volume = eptr->SectorVolume(j);
                  eptr->Store( j, 0U, sv_key, makeScalar( PLAIN, sector_volume ) );
                  eptr->Store( j, 0U, spv_key, makeScalar( PLAIN, sector_volume * porosity ) );
               }
          }
     }

   // 5. computing FV flux balances over the complete stencils
   // --------------------------------------------------------
   if ( initialize_flux ) {
        // loop over FV stencils, computing the relevant variable values

        const typename vector<Node<dim>*>::iterator nit_end(gref_.NodesEnd());
        //double64 bmin(1e30), bmax(-1e30);

        for ( typename vector<Node<dim>*>::iterator nit=gref_.NodesBegin(); nit!=nit_end; ++nit )
          if ( (*nit)->AtBoundary() != NOT )
            {
               const size_t parent_elements((*nit)->Parents());
               double64 flux_balance(0.);
               for ( size_t i=0U; i<parent_elements; ++i ) {
                    const Element<dim>* const eptr = (*nit)->Parent(i);
                    const size_t sector_node      = (*nit)->ParentNodeNumber(i);
                    for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
                         const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
                         const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
                         const double64 facet_flux = sign * eptr->Read( facet, 0U, ff_key );
                         flux_balance += facet_flux;
                      }
                 }
               (*nit)->Store( fb_key, makeScalar((*nit)->Status(fb_key),flux_balance) );

               //bmin = std::min( bmin, (*nit)->Read( fb_key ) );
               //bmax = std::max( bmax, (*nit)->Read( fb_key ) );            
            }
        //cout <<"\nDESTransport: initializeFiniteVolumeProperties: initial flux balance: "<< std::max(fabs(bmin), fabs(bmax)) << endl;
     }
 } // end initializeFiniteVolumeProperties


//Compute the flux balance and CFL number for a node/FV
template<size_t dim>
void DESTransport<dim>::ComputeFluxBalanceAndCFL( Node<dim>* nd )
{
    assert( nd  != NULL );

    double64 flux_balance(0.), outflow(0.);
    const size_t node_parent_elements(nd->Parents());

    //ignore DIRICH node
    if(nd->Status(  C0_key ) == DIRICH) {    
        nd->Store(  fb_key, makeScalar( nd->Status(  fb_key ), 0. ) );//flux balance
        nd->Store(  C1_key, makeScalar( nd->Status(  C1_key ), nd->Read(  C0_key ) ) ); //new concentration=concentration
        nd->dt_CFL_ = numeric_limits<double64>::max(); //CFL time increment
        return;
    };

    VectorVariable<dim>  nrml, vD;

    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
        Element<dim>* const eptr(nd->Parent(t));
        assert( eptr != NULL );
        const size_t pnid(nd->ParentNodeNumber(t));
        eptr->Read(  vD_key, vD);

        //check is vD is nan
        bool isNAN(false);
        for (size_t n=0U;n<dim;++n) {
            if (isnan(vD(n))){
                isNAN=true;
                break;
            };
        };

        if (!isNAN){
            const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
            for ( size_t i=0U; i<sector_facets; i++ )
            {
                const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                const size_t inside_node(eptr->FV()->InsideNode(iFacet));
                eptr->Read( iFacet, 0U,  fn_key, nrml );
                const double64  vD_n = vD.DotProduct(nrml);
                const double64  facetArea = eptr->Read( iFacet, 0U,  fA_key );       
                const double64  sign = ( pnid == inside_node ) ? 1. : -1.;
                flux_balance += sign * vD_n * facetArea;
                if ( sign * vD_n * facetArea > 0. ) outflow += sign * vD_n * facetArea;           
            }
        }
    }
    nd->Store(  fb_key, makeScalar( nd->Status(  fb_key ), flux_balance ) );//flux balance
    if (outflow < numeric_limits<double64>::epsilon()) nd->dt_CFL_ = numeric_limits<double64>::max();
    else nd->dt_CFL_ = nd->Read(  PV_key ) / outflow; //CFL time increment 
}



//Compute the rate of change in a node/FV
template<size_t dim>
void DESTransport<dim>::ComputeRateofChange( Node<dim>* nd )
{
    assert( nd  != NULL );

    //ignore DIRICH node
    if(nd->Status(  C0_key ) == DIRICH) {   
        return;
    };

    rate_count_++;//recording
    nd->Store(  rate_key, makeScalar( nd->Status( rate_key), nd->Read( rate_key) + 1 ) );

    double64 accumulation(0.);
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
            const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
            for ( size_t i=0U; i<sector_facets; i++ )
            {
                const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                const size_t inside_node(eptr->FV()->InsideNode(iFacet));
                const size_t outside_node(eptr->FV()->OutsideNode(iFacet));

                eptr->Read( iFacet, 0U,  fn_key, nrml );
                const double64  vD_n = vD.DotProduct(nrml);
                const double64  facetArea = eptr->Read( iFacet, 0U,  fA_key );       
                // finding the upstream concentration
                const double64 C_upstream = (vD_n < 0.) ? eptr->N(outside_node)->Read(  C0_key ) :
                                                          eptr->N(inside_node)->Read(  C0_key );
                  
                const double64 sign = ( pnid == inside_node ) ? 1. : -1.;
                accumulation += sign * vD_n * facetArea * C_upstream;           
            }
        }
    }
    nd->Store(  C1_key, makeScalar( nd->Status(  C1_key ), accumulation ) );     
} 



//schedule an event associated with a node/FV
template<size_t dim>
bool DESTransport<dim>::Schedule(Node<dim>* nd, double64 t_end, double64 cfl_multiplier)
{
    //ignore DIRICH node
    if(nd->Status(  C0_key ) == DIRICH) {
        return false;
    } else {
        nd->Store(  schedule_key, makeScalar( nd->Status( schedule_key), nd->Read( schedule_key) + 1 ) );
        nd->valid_ = true;
        //compute target change
        double64 CFL = nd->dt_CFL_;//CFL number
        double64 PV = nd->Read( PV_key);//Pore volume
        double64 ChangeRate = nd->Read(  C1_key);//rate of change
        double64 C0 = nd->Read(  C0_key);//concentration
        double64 flux_balance = nd->Read(  fb_key);//flux balance

        double64 dC_CFL = -CFL*cfl_multiplier/PV*(ChangeRate-C0*flux_balance);//targe change

        if (fabs(dC_CFL) < numeric_limits<double64>::epsilon()){//idle node/FV
            nd->dC_tr_ = numeric_limits<double64>::epsilon();
            nd->dt_tr_ = numeric_limits<double64>::max();
        } else {
            nd->dC_tr_ = dC_CFL;//targe concentration change
            nd->dt_tr_ = cfl_multiplier*nd->dt_CFL_;//targe timing
        };

        if ((nd->dt_tr_+nd->t_current_) >= t_end) {
            return false;
        } else {
            nd->t_next_ = nd->t_current_ + nd->dt_tr_; 
            return true;
        };
    };
}


//update solution and check it against the specified range (with DES)
template<size_t dim>
void DESTransport<dim>::Update_DES(Node<dim>* nd, double64 t_clock)
{
    update_count_++;//recording
    const VARIABLE_FLAG status(nd->Status(  C0_key ));
    if ( status != DIRICH )
    {    
        // 1. starting with the sum of facet flux-concentration products stored in C1
        double64 ChangeRate = nd->Read( C1_key);   
        double64 solution = nd->Read( C0_key);//concentration
        // 2. correcting this sum for div vD using 'flux balance'   
        ChangeRate -= solution * nd->Read(  fb_key ); 
        // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
        double64 new_solution = solution - ((t_clock - nd->t_current_)/nd->Read( PV_key)) * ChangeRate;
        // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
        // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
        const double64 source(nd->Read( nsrc_key));
        // new concentration
        new_solution += source * (t_clock - nd->t_current_);

        double64 C_last = nd->Read(  C0_key);//last concentration
        
        if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(  C0_key, makeScalar( status, new_solution ));//store solution value to C0_key
        else {
            cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( new_solution > upper_limit_ ) nd->Store(  C0_key, makeScalar( status, upper_limit_ ) );
            else if ( new_solution < lower_limit_ ) nd->Store(  C0_key, makeScalar( status, lower_limit_ ) );
        }
        
        double64 C_current = nd->Read(  C0_key);//current concentration
        nd->dC_cumulative_ += (C_current-C_last);//update cumulative change
        nd->t_current_ = t_clock;

        nd->Store(  update_key, makeScalar( nd->Status( update_key), nd->Read( update_key) + 1 ) );
    }
}


//update solution and check it against the specified range (with TDS)
template<size_t dim>
void DESTransport<dim>::Update_TDS(Node<dim>* nd, double64 delta_t)
{
    update_count_++;//recording
    const VARIABLE_FLAG status(nd->Status(  C0_key ));
    if ( status != DIRICH )
    {    
        // 1. starting with the sum of facet flux-concentration products stored in C1
        double64 ChangeRate = nd->Read( C1_key);   
        double64 solution = nd->Read( C0_key);//concentration
        // 2. correcting this sum for div vD using 'flux balance'   
        ChangeRate -= solution * nd->Read(  fb_key ); 
        // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
        double64 new_solution = solution - (delta_t/nd->Read( PV_key)) * ChangeRate;
        // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
        // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
        const double64 source(nd->Read( nsrc_key));
        // new concentration
        new_solution += source * delta_t;
 
        if ( new_solution <= upper_limit_ && new_solution >= lower_limit_ ) nd->Store(  C0_key, makeScalar( status, new_solution ));//store solution value to C0_key
        else {
            cerr <<"value: "<< new_solution <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( new_solution > upper_limit_ ) nd->Store(  C0_key, makeScalar( status, upper_limit_ ) );
            else if ( new_solution < lower_limit_ ) nd->Store(  C0_key, makeScalar( status, lower_limit_ ) );
        }
    }
}



//Synchronize neighbor nodes/FVs
template<size_t dim>
void DESTransport<dim>::Synchronize(Node<dim>* nd,double64 t_clock)
{
    nd->Store(  synchronize_key, makeScalar( nd->Status( synchronize_key), nd->Read( synchronize_key) + 1 ) );
    nd->valid_ = false;
    nd->dC_cumulative_ = 0.;
    for ( size_t n=0U; n<nd->Neighbors(); ++n ) {
        Node<dim>* e = nd->Neighbor(n);
        if (e->inPEPStack_ == false) {
            PEPStack.push_back(e);
            e->inPEPStack_ = true;
            Update_DES(e,t_clock);
            if (fabs(e->dC_cumulative_) >= fabs(e->dC_tr_)) Synchronize (e, t_clock); 
        };
    };
}


#if defined(_OPENMP )
//Synchronize neighbor nodes/FVs
template<size_t dim>
void DESTransport<dim>::Synchronize_openmp(Node<dim>* nd,double64 t_clock, size_t num_threads)
{
    std::vector<Node<dim>*> UpdateList, SynList;
    SynList.push_back(nd);
    
    while (SynList.size() > 0) {
        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp for schedule(dynamic)
            for(size_t i = 0U; i < SynList.size(); ++i)
            {
                auto it1 = SynList.begin()+i;
                (*it1)->valid_ = false;
                (*it1)->dC_cumulative_ = 0.;        
                for ( size_t n=0U; n<(*it1)->Neighbors(); ++n ) {
                    Node<dim>* e = (*it1)->Neighbor(n);
                    #pragma omp critical
                    {    
                        if (e->inPEPStack_ == false) {
                            PEPStack.push_back(e);
                            e->inPEPStack_ = true;
                            UpdateList.push_back(e);
                        }
                    }
                }
            }
            #pragma omp single    
            SynList.clear();    

            #pragma omp for schedule(dynamic)
            for(size_t i = 0U; i < UpdateList.size(); ++i)
            {
                auto it2 = UpdateList.begin()+i;
                Update_DES(*it2,t_clock);
                if (fabs((*it2)->dC_cumulative_) >= fabs((*it2)->dC_tr_)) {
                    #pragma omp critical
                    SynList.push_back(*it2);
                }
            }
            #pragma omp single 
            UpdateList.clear();
        }
    }
}
#endif


//advect variable with TDS (time-driven simulation)
template<size_t dim>
void DESTransport<dim>::AdvectVariable_TDS( double64 time_interval, double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter )
{
    double64 time_increment(time_interval); 
    const typename vector<Node<dim>*>::iterator stack_end(PEPStack.end());
    for ( typename vector<Node<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )
    {     
        ComputeRateofChange((*it));  
        time_increment=min(time_increment, (*it)->dt_CFL_*cfl_multiplication_factor);
    }

    const double64 one(1.);
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);

    size_t   substep(1);
    double64 time(0.);
    
    while (time < time_interval)
    {
        cout <<"\n\tadvection (sub)step: "<< substep << endl;
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
        for ( typename vector<Node<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )
        { 
            Update_TDS((*it),time_increment);
        };

        double64 new_time_increment(time_interval); 
        for ( typename vector<Node<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )
        { 
            ComputeRateofChange((*it));
            new_time_increment=min(new_time_increment, (*it)->dt_CFL_*cfl_multiplication_factor);
        };

        time += time_increment;
        substep++;
    };
    cout<<endl<<"rate_count_ = "<<rate_count_<<" update_count_ = "<<update_count_<<endl;
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
    cerr <<"WARNING: OpenMP is not available, using serial mode"<<endl;
    AdvectVariable_DES_serial( model_time, cfl_multiplication_factor, PEP_parameter );
#endif
}


#if defined(_OPENMP )
//advect variable with DES (discrete event simulation), with openmp
template<size_t dim>
void DESTransport<dim>::AdvectVariable_DES_openmp( double64 model_time, double64 cfl_multiplication_factor, double64 PEP_parameter, size_t num_threads)
{
    double64 begin=omp_get_wtime();
    cout<<"Start DESTransport<dim>::AdvectVariable_DES_openmp "<<endl;
    double64 time(0.);
    bool Finished = false;
    //uncomment for recording events at each time interval
    /*
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
    {     
        (*nit)->Store(  update_key, makeScalar( (*nit)->Status( update_key), 0 ) );
    }
    */

    double64 T_begin;

    while (!Finished)
    {
        T_begin = omp_get_wtime();
        size_t PEPStack_size = PEPStack.size();
        cout << "Number of threads = "<< omp_get_max_threads() << " PEPStack_size = " << PEPStack_size;
        #pragma omp parallel num_threads(num_threads)
        {
        #pragma omp for schedule(dynamic)
            for(size_t i = 0U; i < PEPStack_size; ++i)
            {
                auto it = PEPStack.begin()+i;
                ComputeRateofChange(*it);            
                if ((*it)->valid_ == false)
                    if (Schedule(*it, model_time, cfl_multiplication_factor)){
                    #pragma omp critical
                        Queue.push_back(*it);
                    }         
                (*it)->inPEPStack_ = false;             
            };
        }   
        T_RateOfChange_ += omp_get_wtime() - T_begin;
        cout <<" Queue size = "<< Queue.size()<<endl;   

        T_begin= omp_get_wtime();        
        if (Queue.empty()) {
            time=model_time;
        } else {
            sort(Queue.begin(),Queue.end());
            const typename vector<Node<dim>*>::iterator begin(Queue.begin());
            time = (*begin)->t_next_;
        };
        T_SortQueue_ += omp_get_wtime() - T_begin;
        cout<<"  time = "<<time<<" model_time = "<<model_time<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Node<dim>*>::iterator End(PEPStack.end());
            for ( typename vector<Node<dim>*>::iterator nd=PEPStack.begin(); nd!=End; ++nd )
                (*nd)->valid_ = false;
            break;
        };

        PEPStack.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
        size_t count = 0U;
        while (!Queue.empty())
        {
            count++;
            const typename vector<Node<dim>*>::iterator top(Queue.begin());
            dt_PEP = min(dt_PEP, PEP_parameter*(*top)->dt_tr_);
            if ((*top)->t_next_ > (time+dt_PEP)) break;
            if ((*top)->inPEPStack_ == false) {
                PEPStack.push_back((*top));
                (*top)->inPEPStack_ = true;
            }; 
            T_begin= omp_get_wtime();
            Update_DES((*top),time);
            T_Update_ += omp_get_wtime() - T_begin;
            T_begin= omp_get_wtime();
            Synchronize_openmp((*top),time,num_threads);
            T_Synchronize_ += omp_get_wtime() - T_begin;
            Queue.erase(top); 
            //remove invalid events/nodes from Queue
            T_begin= omp_get_wtime();
            if(!Queue.empty()) {
                size_t queue_size = Queue.size();
                for ( size_t n = 0U; n<queue_size; n++ )
                {  
                    if (Queue[n]->valid_ == false) {
                        Queue.erase(Queue.begin()+n);
                        queue_size --;
                    };
                };
            };
            T_RemoveFromQueue_ += omp_get_wtime() - T_begin;
        };
        cout<<"  count =  "<<count<<endl;
    };
    T_AdvectVariable_+=omp_get_wtime()-begin; 

    cout<<"Finish DESTransport<dim>::AdvectVariable "<<endl;
    cout<<"rate_count_ = "<<rate_count_<<" update_count_ = "<<update_count_<<endl;    
    cout <<"T_Schedule_ = "<< T_Schedule_ << endl;
    cout <<"T_SortQueue_ = "<< T_SortQueue_ << endl;
    cout <<"T_Update_  = "<< T_Update_ << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_ << endl;
    cout <<"T_RemoveFromQueue_ = "<< T_RemoveFromQueue_ << endl;   
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ << endl;   
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
        (*nit)->Store(  update_key, makeScalar( (*nit)->Status( update_key), 0 ) );
    }
    */
    
    while (!Finished)
    {    
        clock_t T_begin= clock();
        const typename vector<Node<dim>*>::iterator stack_end(PEPStack.end());
        for ( typename vector<Node<dim>*>::iterator it=PEPStack.begin(); it!=stack_end; ++it )
        {
            ComputeRateofChange((*it));            
            if ((*it)->valid_ == false)
                if (Schedule((*it), model_time, cfl_multiplication_factor))
                    Queue.push_back((*it));      
           (*it)->inPEPStack_ = false;
        };    
        T_RateOfChange_ += clock() - T_begin; 
        cout <<" Queue size = "<< Queue.size()<<endl;   

        T_begin= clock();        
        if (Queue.empty()) {
            time=model_time;
        } else {
            sort(Queue.begin(),Queue.end());
            const typename vector<Node<dim>*>::iterator begin(Queue.begin());
            time = (*begin)->t_next_;
        };
        T_SortQueue_ += clock() - T_begin;
        cout<<"  time = "<<time<<" model_time = "<<model_time<<endl;

        if (time == model_time) {
            Finished = true;
            const typename vector<Node<dim>*>::iterator End(PEPStack.end());
            for ( typename vector<Node<dim>*>::iterator nd=PEPStack.begin(); nd!=End; ++nd )
                (*nd)->valid_ = false;
            break;
        };

        PEPStack.clear();
    
        double64 dt_PEP=numeric_limits<double64>::max();
        size_t count = 0U;
        while (!Queue.empty())
        {
            count++;
            const typename vector<Node<dim>*>::iterator top(Queue.begin());
            dt_PEP = min(dt_PEP, PEP_parameter*(*top)->dt_tr_);
            if ((*top)->t_next_ > (time+dt_PEP)) break;
            if ((*top)->inPEPStack_ == false) {
                PEPStack.push_back((*top));
                (*top)->inPEPStack_ = true;
            }; 
            T_begin= clock();
            Update_DES((*top),time);
            T_Update_ += clock() - T_begin;
            T_begin= clock();
            Synchronize((*top),time);
            T_Synchronize_ += clock() - T_begin;
            Queue.erase(top); 
            //remove invalid events/nodes from Queue
            T_begin= clock();
            if(!Queue.empty()) {
                size_t queue_size = Queue.size();
                for ( size_t n = 0U; n<queue_size; n++ )
                {  
                    if (Queue[n]->valid_ == false) {
                        Queue.erase(Queue.begin()+n);
                        queue_size --;
                    };
                };
            };
            T_RemoveFromQueue_ += clock() - T_begin;
        };
        cout<<"  count =  "<<count<<endl;
    };
    T_AdvectVariable_+= clock() - begin;

    cout<<"Finish DESTransport<dim>::AdvectVariable "<<endl;
    cout<<"rate_count_ = "<<rate_count_<<" update_count_ = "<<update_count_<<endl;    
    cout <<"T_Schedule_ = "<< T_Schedule_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_SortQueue_ = "<< T_SortQueue_ /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Update_  = "<< T_Update_  /double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_Synchronize_ = "<< T_Synchronize_/double64(CLOCKS_PER_SEC) << endl;
    cout <<"T_RemoveFromQueue_ = "<< T_RemoveFromQueue_/double64(CLOCKS_PER_SEC) << endl; 
    cout <<"T_RateOfChange_ = "<< T_RateOfChange_ /double64(CLOCKS_PER_SEC) << endl;  
    cout <<"T_AdvectVariable_ = "<< T_AdvectVariable_ /double64(CLOCKS_PER_SEC) << endl;
}   


template class DESTransport<1U>;
template class DESTransport<2U>;
template class DESTransport<3U>;

} // end csmp  
    
    
