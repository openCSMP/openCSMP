#include "NodeCenteredFiniteVolumeTransport.h"
#include "NodeCenteredFiniteVolumeAlgorithm.h"
#include "Element.h"
#include "StencilProcessor.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include "meshManagementUtilities.h"
#include "Box.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "StencilProcessor.h"
#include "TwoPhaseModel.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using namespace std;

namespace csmp {

/** Passive (linear) advection of tracer in single-phase flow.

@section arguments Input Arguments

The model and basic transport-related variables.
*/
template<uint32_t dim>
NodeCenteredFiniteVolumeTransport<dim>::NodeCenteredFiniteVolumeTransport( const char* group_name,
                                                                           Model<dim>& sg,
                                                                           const char* porosity,
                                                                           const char* advected_prop,
                                                                           const char* transp_velocity,
                                                                           const char* nodal_source,
                                                                           bool second_order_in_space,
                                                                           bool second_order_in_time,
                                                                           const char* elmt_thickness_attribute,
                                                                           const char* velocity_multiplier)
    : pref_(sg.Database()),
      gref_(sg.Region(group_name)),
      mref_(sg),
      advected_variable_(advected_prop),
      phi_key_(sg.Database().StorageKey(porosity)),
      // dif_key_: place=MODEL, index=ULONG_MAX
      adv1_key_(sg.Database().StorageKey(advected_prop)),
      vel_key_(sg.Database().StorageKey(transp_velocity)),
      src_key_(sg.Database().StorageKey(nodal_source)),
      thi_key_( ((elmt_thickness_attribute==NULL) ? csmp::Index() : sg.Database().StorageKey(elmt_thickness_attribute)) ),
      velo_mult_key_( ((velocity_multiplier==NULL) ? csmp::Index() : sg.Database().StorageKey(velocity_multiplier)) ),
      stencil_( adv1_key_, vel_key_, src_key_, velo_mult_key_ ),
      cfl_multiplier_(1.),
      firstCall_(true),
      baseAdvector_(NULL),
      grad_advprop_limiter_(NULL),
      with_lsmgrad_limiter_(false),
      target_nonlinear_limiting_case_residual_(1.0e-3),
      max_nonlinear_limiting_case_iterations_(10U),verbose_(true)
{
    /// Thickness should not be multiplied if velocity was assumed to be de-scaled.
    /// for legacy purposes, the transport class assumes that the supplied velocity has already been scaled
    /// to adjust for facet areas coming from the thickness attribute.  This is perhaps a little unintuitive,
    /// since any new user needs to know beforehand that this has to happen.  The user should always provide the
    /// "real" velocity (not scaled) and the transport class should deal with thickness through internal scaling of its
    /// facet areas.
    /// -- Julian M. 23-09-2013

    CheckTransportVariables();

    gref_.UpdateMemberIndexes();

    bool multiply_with_thickness_attribute = (elmt_thickness_attribute==NULL) ? false : true;
    InitializeFiniteVolumeData( multiply_with_thickness_attribute );

    InitializeArraysForFirstOrderMethod();

    if ( second_order_in_space || second_order_in_time )
        InitializeArraysForSecondOrderMethod( second_order_in_space, second_order_in_time );

    size_t allocated_memory(MeasureAllocatedMemory());

    if(this->Verbose()){
        cout <<"\nNodeCenteredFiniteVolumeTransport(constructing module for single-phase passive solute transport): ";

        if ( allocated_memory >= 1e6 )
            cout <<"\nBasic storage allocated data by transport algorithm (MB): "<< allocated_memory/1e6 <<"."<< endl;
        else
            cout <<"\nBasic storage allocated data by transport algorithm (bytes): "<< allocated_memory <<"."<< endl;
        cout <<"\nNodeCenteredFiniteVolumeTransport(constructor): Constructed successfully."<< endl << endl;
    }

} // end constructor (solute advection only)



/** Passive (linear) advection of tracer in single-phase flow, 1st-order method.

@section arguments Input Arguments

The model and basic transport-related variables.

@section implementation Implementation

Maybe later this can be made more flexible by turning the vectors into
deque's that can grow on either side.
 */
template<uint32_t dim>
NodeCenteredFiniteVolumeTransport<dim>::NodeCenteredFiniteVolumeTransport( const char* group_name,
                                                                           Model<dim>& sg,
                                                                           const char* porosity,
                                                                           const char* diffusivity,
                                                                           const char* advected_prop,
                                                                           const char* transp_velocity,
                                                                           const char* nodal_source,
                                                                           bool second_order_in_space,
                                                                           bool second_order_in_time,
                                                                           const char* elmt_thickness_attribute,
                                                                           const char* velocity_multiplier)
    : pref_(sg.Database()),
      gref_(sg.Region(group_name)),
      mref_(sg),
      advected_variable_(advected_prop),
      phi_key_(sg.Database().StorageKey(porosity)),
      diff_key_(sg.Database().StorageKey(diffusivity)),
      adv1_key_(sg.Database().StorageKey(advected_prop)),
      vel_key_(sg.Database().StorageKey(transp_velocity)),
      src_key_(sg.Database().StorageKey(nodal_source)),
      thi_key_( ((elmt_thickness_attribute==NULL) ? csmp::Index() : sg.Database().StorageKey(elmt_thickness_attribute)) ),
      velo_mult_key_( ((velocity_multiplier==NULL) ? csmp::Index() : sg.Database().StorageKey(velocity_multiplier)) ),
      stencil_( adv1_key_, vel_key_, diff_key_, src_key_, velo_mult_key_ ),
      cfl_multiplier_(1.),
      firstCall_(true),
      baseAdvector_(NULL),
      grad_advprop_limiter_(NULL),
      with_lsmgrad_limiter_(false),
      target_nonlinear_limiting_case_residual_(1.0e-3),
      max_nonlinear_limiting_case_iterations_(10U),verbose_(true)
{
    /// Thickness should not be multiplied if velocity was assumed to be de-scaled.
    /// for legacy purposes, the transport class assumes that the supplied velocity has already been scaled
    /// to adjust for facet areas coming from the thickness attribute.  This is perhaps a little unintuitive,
    /// since any new user needs to know beforehand that this has to happen.  The user should always provide the
    /// "real" velocity (not scaled) and the transport class should deal with thickness through internal scaling of its
    /// facet areas.
    /// -- Julian M. 23-09-2013

    CheckTransportVariables();

    gref_.UpdateMemberIndexes();

    bool multiply_with_thickness_attribute = (elmt_thickness_attribute==NULL) ? false : true;
    InitializeFiniteVolumeData( multiply_with_thickness_attribute );

    InitializeArraysForFirstOrderMethod();

    if ( second_order_in_space || second_order_in_time )
        InitializeArraysForSecondOrderMethod( second_order_in_space, second_order_in_time );

    size_t allocated_memory(MeasureAllocatedMemory());

    if (this->Verbose()){
    cout <<"\nNodeCenteredFiniteVolumeTransport(constructing module for single-phase passive solute transport): ";
    if ( allocated_memory >= 1.E+6 )
        cout <<"\nBasic storage allocated data by transport algorithm (MB): "<< allocated_memory/1e6 <<"."<< endl;
    else
        cout <<"\nBasic storage allocated data by transport algorithm (bytes): "<< allocated_memory <<"."<< endl;
    cout <<"\nNodeCenteredFiniteVolumeTransport(constructor): Constructed successfully."<< endl << endl;
    }

} // end constructor (solute transport)





template<uint32_t dim>
NodeCenteredFiniteVolumeTransport<dim>::~NodeCenteredFiniteVolumeTransport()
{
    delete baseAdvector_;
    delete grad_advprop_limiter_;
}

template<uint32_t dim>
CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS& NodeCenteredFiniteVolumeTransport<dim>::GetSolverSettings()
{

    if (firstCall_){
        baseAdvector_= new NodeCenteredFiniteVolumeAlgorithm<dim>( gref_ );
        firstCall_ = false;
    }

    return baseAdvector_->GetSolverSettings();
}

template<uint32_t dim>
CSMP_DEFAULT_LINEAR_SOLVER* NodeCenteredFiniteVolumeTransport<dim>::GetSolver()
{

    if (firstCall_){
        baseAdvector_= new NodeCenteredFiniteVolumeAlgorithm<dim>( gref_ );
        firstCall_ = false;
    }

    return baseAdvector_->GetSolver();
}

template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AdjustSolverSettings()
{
#ifdef CSMP_WITH_SAMG_SOLVER
    std::cout <<"\n\n*** NodeCenteredFiniteVolumeTransport::AdjustSolverSettings ***\n\n";

    GetSolverSettings().Set_iout1( 0 );
    GetSolverSettings().Set_iout2( 0 );

    /// SAMG solution criterion
    GetSolverSettings().Set_eps(0.);
    GetSolverSettings().Set_rel_eps(1.E-10);

    /// SAMG output to file
    //GetSolverSettings().Set_idmp( 8 );        // define SAMG command and file output
    //GetSolverSettings().Set_ioform( "f" );    // define SAMG file output format for reduced file size, idmp > 1 is required
    //GetSolverSettings().Set_filnam_dump( "SAMG_Transport" ); // set filename for SAMG file output other than default "level", idmp > 1 is required
#else
    /// add extra functinality for alternative solver
#endif
} // end AdjustSolverSettings


/** Uses the initializations of internal variables to detect whether the
transport object is setup for second-order in time computations.
*/
template<uint32_t dim>
bool NodeCenteredFiniteVolumeTransport<dim>::SecondOrderInTime() const
{
    // all thats needed for the second order in time scheme
    if ( !FACETFLUXES0.empty() && !LTDSATS0.empty() && !SAT0.empty() ) return true;

    return false;
}



/** If a node connectivity graph is present in the algorithm it is suited
for second-order in space computations.  */
template<uint32_t dim>
bool NodeCenteredFiniteVolumeTransport<dim>::SecondOrderInSpace() const
{
    // all that is needed for the second order in time scheme
    if ( !SMINMAX.empty() ) return true;

    return false;
}


template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::MaxNonlinear2ndOrderIterations(size_t max_nonlinear_iterations)
{
    max_nonlinear_limiting_case_iterations_= max_nonlinear_iterations;
}

template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::TargetNonlinear2ndOrderResidual(double target_residual)
{
    target_nonlinear_limiting_case_residual_ = target_residual;
}




/**

Build the finite volume stencil data array and computes properties that
stay constant, like the facet areas, the sector volumes, and the normals
to the facets.

The finite volume stencil data are stored such that the stencils inside
the target group are numbered contiguously. If the group lies somewhere
inside the model, the halo stencils are appended as a further contiguous
set to the stencil data.

As option, sector pore volumes of lower dimensional elements can be multiplied with
the thickness attribute.

*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::InitializeFiniteVolumeData( bool multiply_pore_volumes_with_thickness )
{
    STENCIL_DATA.resize( gref_.Cells() );
    vector<FV_Parameter>(STENCIL_DATA).swap(STENCIL_DATA);

    // verifying that the thickness key has actually been initialized
    if ( multiply_pore_volumes_with_thickness )
        assert( thi_key_ != csmp::Index() );

    double poro;
    double thi(1.0);
    
    // initializing FV stencil data
    const bool store_normals(true);
    for ( size_t i{0U}; i<gref_.Cells(); i++ ) {
        // resizing the data vectors
        STENCIL_DATA[i].Resize( gref_.E(i)->FV()->Sectors(),
                                gref_.E(i)->FV()->Facets(), dim, store_normals );

        // computing the sector pore volumes
        if ( multiply_pore_volumes_with_thickness ) thi = gref_.E(i)->Read( thi_key_ );

        if (phi_key_.place == ELEMENT)
        {
            poro = gref_.E(i)->Read( phi_key_ );
            for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
                STENCIL_DATA[i].SectorVolume( j, poro *thi* ( *gref_.E(i) ).SectorVolume(j) );
        }
        else if (phi_key_.place == SECTOR_INTEGRATION_POINT)
        {
            for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
            {
                poro = gref_.E(i)->Read( j, 0U,  phi_key_ );
                STENCIL_DATA[i].SectorVolume( j, poro *thi* ( *gref_.E(i) ).SectorVolume(j) );
            }
        }
        else if (phi_key_.place == NODE)
        {
            for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
            {
                poro = gref_.E(i)->N(j)->Read( phi_key_ );
                STENCIL_DATA[i].SectorVolume( j, poro *thi* ( *gref_.E(i) ).SectorVolume(j) );
            }
        }
        
        // computing the facet areas
        for ( auto j{0U}; j<gref_.E(i)->FV()->Facets(); j++ )
            STENCIL_DATA[i].FacetArea( j, ( *gref_.E(i) ).FacetArea(j) );

        // computing the facet normals
        for ( auto j{0U}; j<gref_.E(i)->FV()->Facets(); j++ )
            STENCIL_DATA[i].FacetNormal( j, ( *gref_.E(i) ).FacetNormal(j).Coordinates() );
    }

} // end InitializeFiniteVolumeData






/** 
      As in InitializeFiniteVolumeData(), but all quantities are mapped from parametric to physical space.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::InitializeFiniteVolumeDataParametricToPhysical( bool store_normals )
{
    STENCIL_DATA.resize( gref_.Cells() );
    vector<FV_Parameter>(STENCIL_DATA).swap(STENCIL_DATA);

    // initializing FV stencil data
    for ( size_t i{0U}; i<gref_.Cells(); i++ ) {
        // resizing the data vectors
        STENCIL_DATA[i].Resize( gref_.E(i)->FV()->Sectors(),
                                gref_.E(i)->FV()->Facets(), dim, store_normals );

        // computing the sector pore volumes
        double phi = gref_.E(i)->Read( phi_key_ );

        for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
            STENCIL_DATA[i].SectorVolume( j, phi * ( *gref_.E(i) ).SectorVolume(j) );

        // computing the facet areas
        for ( auto j{0U}; j<gref_.E(i)->FV()->Facets(); j++ )
            STENCIL_DATA[i].FacetArea( j, ( *gref_.E(i) ).FacetAreaMapped(j) );

        // computing the facet normals
        if ( store_normals )
            for ( auto j{0U}; j<gref_.E(i)->FV()->Facets(); j++ )
                STENCIL_DATA[i].FacetNormal( j, ( *gref_.E(i) ).FacetNormalMapped(j).Coordinates() );
    }

} // end InitializeFiniteVolumeDataParametricToPhysical




/**
     Computes sector pore volumes, by multiplying them with the porosity and a thickness attribute
     that comes into play for lower dimensional elements.
     
     @attention thickness attribute should be 1 for elements of dim=dim, and the true thickness (or
     cross-sectional area fraction of 1) for lower dimensional elements.
     
     The results are stored in the FV stencil data array.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::InitializeSectorPoreVolumeData( bool multiply_pore_volumes_with_thickness )
{
    assert( STENCIL_DATA.size() == gref_.Cells() );
    if (multiply_pore_volumes_with_thickness){
        assert( thi_key_.type == SCALAR );
        assert( thi_key_.place == ELEMENT );
    }

    double poro;
    double thi(1.0);
    
    if (phi_key_.place == ELEMENT)
    {
        for ( size_t i{0U}; i<gref_.Cells(); i++ )
        {
            // computing the sector pore volumes
            poro = gref_.E(i)->Read( phi_key_ );
            if ( multiply_pore_volumes_with_thickness ) thi = gref_.E(i)->Read( thi_key_ );
            for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
                STENCIL_DATA[i].SectorVolume( j, poro * thi * ( *gref_.E(i) ).SectorVolume(j) );
        }
    }
    else if (phi_key_.place == SECTOR_INTEGRATION_POINT)
    {
        for ( size_t i{0U}; i<gref_.Cells(); i++ )
        {
            // computing the sector pore volumes
            if ( multiply_pore_volumes_with_thickness ) thi = gref_.E(i)->Read( thi_key_ );
            for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
            {
                poro = gref_.E(i)->Read( j, 0U, phi_key_ );
                STENCIL_DATA[i].SectorVolume( j, poro * thi * ( *gref_.E(i) ).SectorVolume(j) );
            }
        }
    }
    else if (phi_key_.place == NODE)
    {
        for ( size_t i{0U}; i<gref_.Cells(); i++ )
        {
            // computing the sector pore volumes
            if ( multiply_pore_volumes_with_thickness ) thi = gref_.E(i)->Read( thi_key_ );
            for ( auto j{0U}; j<gref_.E(i)->FV()->Sectors(); j++ )
            {
                poro = gref_.E(i)->N(j)->Read( phi_key_ );
                STENCIL_DATA[i].SectorVolume( j, poro * thi * ( *gref_.E(i) ).SectorVolume(j) );
            }
        }
    }

    // the halo stencils are not included into the STENCIL_DATA

} // end InitializeSectorPoreVolumeData






/**

Initializes flux balance and pore volume vectors for computation. The
method assumes that SectorVolume() refers to the pore volume of the
sectors. The stencils are initialized this way by the method
InitializeFiniteVolumeData(), see above.

*/
template<uint32_t dim>
bool NodeCenteredFiniteVolumeTransport<dim>::InitializeArraysForFirstOrderMethod()
{
    FVPOREVOL.resize( gref_.Nodes() );
    vector<double>(FVPOREVOL).swap(FVPOREVOL);
    fill( FVPOREVOL.begin(), FVPOREVOL.end(), 0. );
    FLUX_BALANCE.resize( gref_.Nodes() );
    vector<double>(FLUX_BALANCE).swap(FLUX_BALANCE);
    fill( FLUX_BALANCE.begin(), FLUX_BALANCE.end(), 0. );

    // computing pore volume of each finite volume from its sector volumes (these were already multiplied with phi)
    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
      for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
        FVPOREVOL[ (*eit)->N(i)->Idx() ] += STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);

    var_ncomponents_=1;

    if (adv1_key_.type != SCALAR) {
        if (adv1_key_.type == ARRAY) {
            ArrayVariable av;
            (*gref_.NodesBegin())->Read(adv1_key_,av);
            var_ncomponents_=av.Size();
         }
        else if (adv1_key_.type == FLAGGEDARRAY){
            FlaggedArrayVariable fav;
            (*gref_.NodesBegin())->Read(adv1_key_,fav);
            var_ncomponents_=fav.Size();
         }
        else
        throw csmp::Exception(FATAL_ERROR, "NodeCenteredFiniteVolumeTransport<dim>::InitializeArraysForFirstOrderMethod",
            "Variable type not accepted!");
    }

    return true;

} // end InitializeArraysForFirstOrderMethod





/** 
     Computes min/max of saturation at upstream nodes, and remembers transported variable values
     from current timestep.
 */
template<uint32_t dim>
bool NodeCenteredFiniteVolumeTransport<dim>::InitializeArraysForSecondOrderMethod(
        bool second_order_in_space,
        bool second_order_in_time )
{
    if ( second_order_in_space ) {
        SMINMAX.resize( gref_.Nodes() );  // advected quantity min/max in the neighborhood of each node
        vector<pair<double,double> >(SMINMAX).swap(SMINMAX);
        cout <<"\nNodeCenteredFiniteVolumeTransport::InitializeArraysForSecondOrderMethod: ";
    }

    if ( second_order_in_time ) {
        FACETFLUXES0.resize( gref_.Cells() );
        vector<vector<double> >(FACETFLUXES0).swap(FACETFLUXES0);
        SAT0.resize( gref_.Nodes() );
        vector<double>(SAT0).swap(SAT0);
        LTDSATS0.resize( gref_.Cells() );
        vector<vector<double> >(LTDSATS0).swap(LTDSATS0);

        typename vector<vector<double> >::iterator  fit0 = FACETFLUXES0.begin(),
                lit0 = LTDSATS0.begin();

        cout <<"\nNodeCenteredFiniteVolumeTransport(constructor): Initializing data arrays."<< endl;
        for ( typename vector<Element<dim>*>::const_iterator
              eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++, fit0++, lit0++ )
        {
            // resizing sector flux arrays for each element
            (*fit0).resize( (*eit)->FV()->Facets() );
            vector<double>((*fit0)).swap((*fit0));
            // resizing slope-limited saturation array for each element
            (*lit0).resize( (*eit)->FV()->Facets() );
            vector<double>((*lit0)).swap((*lit0));
        }

        // the halo stencils are not dealt with at this point
    }

    return true; // create option for when memory cannot be allocated

} // end InitializeArraysForSecondOrderMethod



/** Sums the various current allocations in bytes which are returned.

@return The approximate current allocation in bytes used.
*/
template<uint32_t dim>
size_t  NodeCenteredFiniteVolumeTransport<dim>::MeasureAllocatedMemory() const
{
    size_t  allocated_memory(0);

    // simple arrays
    allocated_memory += sizeof(double) * FVPOREVOL.capacity();
    allocated_memory += sizeof(double) * FLUX_BALANCE.capacity();

    for ( vector<FV_Parameter>::const_iterator fvt = STENCIL_DATA.begin(); fvt!=STENCIL_DATA.end(); fvt++ )
        allocated_memory += (*fvt).Bytes();

    if ( !SAT0.empty() )  allocated_memory += sizeof(double) * SAT0.capacity();
    if ( !SMINMAX.empty() ) allocated_memory += sizeof(pair<double,double>) * SMINMAX.capacity();

    // composite arrays
    if ( !FACETFLUXES0.empty() )
        for ( typename vector<vector<double> >::const_iterator i=FACETFLUXES0.begin(); i<FACETFLUXES0.end(); i++ )
            allocated_memory += sizeof(double) * (*i).capacity();

    if ( !LTDSATS0.empty() )
        for ( typename vector<vector<double> >::const_iterator i=LTDSATS0.begin(); i<LTDSATS0.end(); i++ )
            allocated_memory += sizeof(double) * (*i).capacity();

    return allocated_memory;

} // end MeasureAllocatedMemory()



/** Checks the input variables for correct placement and type compatibility.
 */
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables() const
{
    if ( (phi_key_.place != ELEMENT && phi_key_.place != SECTOR_INTEGRATION_POINT && phi_key_.place != NODE) 
        || phi_key_.type != SCALAR ) {
        string  msg ="The 'porosity' variable '";
        msg       += pref_.Name( phi_key_ );
        msg       +="' must be a scalar-type element variable";
        throw csmp::Exception( FATAL_ERROR, "NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables", msg.c_str() );
    }

    if ( (diff_key_.place != ELEMENT || diff_key_.type != SCALAR) and diff_key_.index != ULONG_MAX ) {
        string  msg ="The 'diffusivity' variable '";
        msg       += pref_.Name( diff_key_ );
        msg       +="' must be a scalar-type element variable";
        throw csmp::Exception( FATAL_ERROR, "NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables", msg.c_str() );
    }

    if ( adv1_key_.place != NODE ){
        string  msg ="The advected variable '";
        msg       += pref_.Name( adv1_key_ );
        msg       +="' must be a scalar-type node variable";
        throw csmp::Exception( FATAL_ERROR, "NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables", msg.c_str() );
    }

    if ( vel_key_.place != ELEMENT || vel_key_.type != VECTOR ) {
        string  msg ="The transport velocity variable '";
        msg       += pref_.Name( vel_key_ );
        msg       +="' must be a vector-type element variable";
        throw csmp::Exception( FATAL_ERROR, "NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables", msg.c_str() );
    }

    if ( src_key_.place != NODE ){
        string  msg ="The source variable '";
        msg       += pref_.Name( src_key_ );
        msg       +="' must be a scalar-type node variable";
        throw csmp::Exception( FATAL_ERROR, "NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables", msg.c_str() );
    }

} // end CheckTransportVariables



/** Finds initial values of advected property at the FV centers = nodes and
at facet-midpoints

@section arguments Input Arguments

Output Arguments&amp; Return Value

The current saturations are returned into the vector<double> SAT0 and their
minimum and maximum values are stored in the last two method arguments.
 */
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::InitialAdvectedPropertyValues( const csmp::Index& adv_key )
{
    MinMaxAdvectedProperty();

    // SAT0
    for ( size_t i{0U}; i<gref_.Nodes(); i++ )
        // the saturation at the node is assigned to SAT0 vector
        SAT0[i] = gref_.N(i)->Read( adv_key );

    // LTDSATS0
    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
    {
        stencil_.eidx_ = (*eit)->Idx();

        // computing upstream and interpolating advected property values to FV segments
        stencil_.InitializeSecondOrder( STENCIL_DATA[stencil_.eidx_], *(*eit) );

        // 2.3  update and compute values of advected variable at faces (spatial limiting)
        stencil_.IsotropicallyLimitTransportProperties( *(*eit), SMINMAX );

        LTDSATS0[stencil_.eidx_] = stencil_.ipsi1_;
    }

} // end InitialAdvectedPropertyValues



/** Initializes the SAT0 vector<double> and determines min/max values for the
current saturation.

@section arguments Input Arguments

Region and csmp::Index for the transport variable.

SAT0 is initialized and the range of the transport property upstream is returned
into the last two method arguments.

*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::InitialAdvectedPropertyValues( const csmp::Index& adv_key,
                                                                            double& smin,
                                                                            double& smax )
{
    smin = smax = gref_.N(0)->Read( adv_key );

    for ( size_t i{0U}; i<gref_.Nodes(); i++ ) {
        // the saturation at the node is assigned to SAT0 vector
        SAT0[i] = gref_.N(i)->Read( adv_key );
        smin = std::min( smin, SAT0[i] );
        smax = std::max( smax, SAT0[i] );
    }

} // end InitialAdvectedPropertyValues



/**

Loops over the nodes recording min/max values of the transported property
for of each node and its neighbouring nodes. These will be used to
place bounds on upstream values of the transported property in order to
perform the slope (spatial) limiting required to make  higher-order
transport schemes non-oscillatory.

@note min/max values must be updated in each non-linear iteration to allow
the scheme to converge.

@section arguments Input Arguments

To read the values of the transported property, the method needs access
to the current model (argument 1) and a graph of the connectivity among
nodes (argument 2).

Min/max value pairs for each node are returned into the last method
argument, a vector<double> of pairs of floats.

@section implementation Implementation

The method uses a node connectivity graph to find and access the neighbors
of each node. It is ascertained that the determined smin/smax values do
not exceed the provided bounds in the course of the non-linear
iteration loop.

@section application Application

The determined min/max values are required as bounds for the spatial
limiting in higher-order transport schemes.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::MinMaxAdvectedProperty()
{
    typename vector<pair<double,double> >::iterator  sit(SMINMAX.begin());

    for ( typename vector<Node<dim>*>::const_iterator
          nit=gref_.NodesBegin(); nit!=gref_.NodesEnd(); nit++, sit++ ) {
        // 1. the advected property value at the current node is assigned to min-max pair
        (*sit).first = (*sit).second = (*nit)->Read( adv1_key_ );
        for ( auto i{0U}; i<(*nit)->Neighbors(); i++ ) {
            const double adv_var((*nit)->Neighbor(i)->Read( adv1_key_ ));
            // if element value is smaller the current minimum is assigned etc.
            (*sit).first  = std::min( (*sit).first,  adv_var );
            (*sit).second = std::max( (*sit).second, adv_var );
        }
    } // end for all nodes

} // end MinMaxAdvectedProperty




/**

The facet fluxes from the first time-level are backuped into the vector<double>
FACETFLUXES0 so that they can used in the next calculation to obtain
second-order accuracy in time for the transport scheme.

*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::BackupFacetFluxes()
{
    // all stencils including those for the halo elements
    for ( size_t i{0U}; i<gref_.Cells(); i++ )
        for ( size_t j{0U}; j<FACETFLUXES0[i].size(); j++ )
            FACETFLUXES0[i][j] = STENCIL_DATA[i].FacetNormalVelocity(j) * STENCIL_DATA[i].FacetArea(j);

} // end BackupFacetFluxes



/**

Loops over the finite elements and computes the (costly) projections of
the transport velocities onto the facet normals. Then it uses these to
calculate flux balances for each FV cell which are later used to compensate
for divergence terms that may originate in transient flow fields.

Flux balances are calculated only where entire finite volumes can be
formed around nodes. This is the case for the cells inside a model and
based on the concept of halo stencils when the transport region
is located inside a mesh (i.e. the algorithm is restricted to a group).


The exception are finite volumes at no-flow boundaries which are
truncated but for which a flux balance can be obtained because there
should not be flux across their boundaries.

For finite volume cells which are cut by the boundary of the Region
no balance is computed. Here the value of the flux-balance is set to
zero.

@section arguments Input Arguments

The method needs access to the Region to read the velocity values.

The projected velocities are stored in the private variable
'FV_STENCIL_DATA'. The method initializes the FLUX_BALANCE vector<double> that
is supplied as second argument.

@section application Application

The method is called within those methods that compute the Courant
number for the current grid.

@todo fix! - flux balances are not correct for nodes that lie in the inside of the model
but on the perimeter of the target region.

 */
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::UpdateProjectedVelocitiesAndFluxBalances()
{
    gref_.UpdateMemberIndexes();

    VectorVariable<dim>  velo;
    uint32_t             inside_node, outside_node;
    double               flux;

    fill( FLUX_BALANCE.begin(), FLUX_BALANCE.end(), static_cast<double>(0.) );

    typename vector<Element<dim>*>::const_iterator  eit(gref_.CellsBegin());
    typename vector<FV_Parameter>::iterator         stit(STENCIL_DATA.begin());
    assert( gref_.Cells() == STENCIL_DATA.size() );

#if !defined(_OPENMP)
    // for all inside stencils
    while ( eit != gref_.CellsEnd() )
    {
        (*eit)->Read( vel_key_, velo );
        for ( auto i{0U}; i<(*stit).Facets(); i++ )
        {
            (*eit)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
            // inside node
            // computing facet normal velocity
            flux = (*stit).FacetNormalProjection( i, velo );
            // storing it in FV_Parameter object
            (*stit).FacetNormalVelocity( i, flux );
            // updating the flux balance
            FLUX_BALANCE[ (*eit)->N(inside_node)->Idx() ] += flux * (*stit).FacetArea(i);
            // outside node
            // updating facet velocities
            FLUX_BALANCE[ (*eit)->N(outside_node)->Idx() ] -= flux * (*stit).FacetArea(i);
        }
        eit++;
        stit++;
    }
#else
#pragma omp parallel
    {
        size_t tid = omp_get_thread_num();
        Element<dim>* ep;
        double             thread_flux;
        VectorVariable<dim> thread_velo;
        size_t               thread_inside_node, thread_outside_node;
        // for all inside stencils
#pragma omp for
        for ( size_t e = 0U ; e < gref_.Cells();e++)
        {
            ep = gref_.E(e);
            ep->Read( vel_key_, thread_velo );

            //----------------------------------------------------
            //change the element stencil to one for this thread, temporarily.
            FiniteElement* fe_tmp=ep->FE();
            // change pointer here
            ep->Assign(this->femgrs_[tid].E(ep->FE_Type()));
            //----------------------------------------------------

            //----------------------------------------------------
            //change the Finite Volume Stencil to one for this thread, temporarily.
            const FiniteVolumeStencil<dim>* tmp_fvstencil= ep->FV();
            ep->Assign(this->fvmgrs_[tid].Stencil( ep->FE_Type()));
            //----------------------------------------------------

            for ( auto i{0U}; i<STENCIL_DATA[e].Facets(); i++ )
            {
                ep->FV()->FacetEdgeNodes( i, thread_inside_node, thread_outside_node );
                // inside node
                // computing facet normal velocity
                thread_flux = STENCIL_DATA[e].FacetNormalProjection( i, thread_velo );
                // storing it in FV_Parameter object
                STENCIL_DATA[e].FacetNormalVelocity( i, thread_flux );
                // updating the flux balance
                FLUX_BALANCE[ ep->N(thread_inside_node)->Idx() ] += thread_flux * STENCIL_DATA[e].FacetArea(i);
                // outside node
                // updating facet velocities
                FLUX_BALANCE[ ep->N(thread_outside_node)->Idx() ] -= thread_flux * STENCIL_DATA[e].FacetArea(i);
            }
            // put the FEM back here
            ep->Assign(fe_tmp);

            //put the FVM stencil back
            ep->Assign(tmp_fvstencil);
        }
    }
#endif

    for ( typename vector<Node<dim>*>::const_iterator
          nit=gref_.PerimeterNodesBegin(); nit!=gref_.NodesEnd(); nit++ )
        FLUX_BALANCE[ (*nit)->Idx() ] = 0.;

} // end UpdateProjectedVelocitiesAndFluxBalances






/** Method loops over FLUX_BALANCE vector<double> determining min/max values.

Returns Min and max flux balances for the finite volumes are returned in the
first and second method arguments, respectively. Usually, the highest
imbalances will occur where the finite volumes are truncated across
model boundaries.

*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::FluxBalance( double& fmin, double& fmax ) const
{
    typename vector<double>::const_iterator it=FLUX_BALANCE.begin();
    fmax = fmin = *it++;

    while ( it != FLUX_BALANCE.end() ) {
        fmin = std::min( fmin, *it );
        fmax = std::max( fmax, *it );
        it++;
    }
}



/**

Method computes the minimum CFL criterion for the current model on the
basis of the highest velocity in the cells contributing to the current
finite volume cell.

@section arguments Input Arguments

The method needs a reference to the current model in order to determine
the flow velocity distribution.

If diffusion is important, a separate velocity computation is necessary.
In this case, an estimate usind N-gradients is used for the diffusion
velocity.

@return The method returns the minimum global timestep for the current mesh
and velocity distribution for which the the transport distance in each
element is smaller or equal to its diameter.

@section implementation Implementation

This is a slightly expensive version because it loops over the
finite volumes rather than the finite elements, then checking each
finite element connected to FV for its velocity etc.
*/

/*
 *
 *
 *
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::CourantIncrement()
{
    if ( dim != 1U )
        throw logic_error("NodeCenteredFiniteVolumeTransport<dim>::CourantIncrement: \
                          this method should only be applied in 1D");

    UpdateProjectedVelocitiesAndFluxBalances();

    static DenseMatrix<DM_MIN>  DN;
    vector<double>            grad(dim);
    VectorVariable<dim>         vc;
    const double              hundred_days(8640000.);
    double                    courant_increment(hundred_days);

    // 0. diffusion is taken into account if the diffusion key is initialized
    const bool with_diffusion( (diff_key_ == csmp::Index()) ? false : true );

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
    {
        // limit imposed by advection
        // --------------------------
        (*eit)->Read( vel_key_, vc );
        double velocity = vc.Length();

        // limit due to diffusion-driven flow
        // ----------------------------------
        if ( with_diffusion ) {
            fill( grad.begin(), grad.end(), 0. );
            (*(*eit)).dN_AtBaryCenter( DN );
            // grad transport variable-based formulation
            for ( size_t j{0U}; j<(*eit)->Nodes(); j++ ) {
                double adv = (*eit)->N(j)->Read( adv1_key_ );
                for ( size_t k{0U}; k<dim; k++ ) grad[k] += DN(k,j) * adv;
            }
            // getting maximum gradient value leading to the maximum diffusive flux
            double  magnitude_grad(grad[0]);

            // using data from barycenter: O.K. as long as grad does not increase during iterations
            velocity += fabs( magnitude_grad * 2. * (*eit)->Read( diff_key_ ) );
        }

        // calculating the CFL criterion from the cell diameter of each sector (pore) volume
        for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
            double cell_diameter = STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
            courant_increment = std::min( courant_increment, cell_diameter / velocity );
        }
    }

    if ( courant_increment >= hundred_days ) {
        cout <<"\nNodeCenteredFiniteVolumeTransport::CourantIncrement: could not be determined."<< endl;
        cout <<"\nthere appears to be no flow in the model domain, returning CFL=3600."<< endl;
        return hundred_days;
    }

    cout <<"\nNodeCenteredFiniteVolumeTransport::CourantIncrement: "<< courant_increment;
    if ( courant_increment <= 1.0e-3 ) cout <<" (CFL constraint is very tight).";
    cout << endl;

    return courant_increment;

} // end CourantIncrement
*/


/**

CFL condition based on the distance from one side of the element to the
other along a line through the element's barycenter and parallel to
the current flow velocity vector.

@section arguments Input Arguments

The method needs a reference to the current model in order to determine
the flow velocity distribution.

@attention Thickness should not be multiplied if velocity was assumed to be de-scaled.
for legacy purposes, the transport class assumes that the supplied velocity has already been scaled
to adjust for facet areas coming from the thickness attribute.  This is perhaps a little unintuitive,
since any new user needs to know beforehand that this has to happen.  The user should always provide the
"real" velocity (not scaled) and the transport class should deal with thickness through internal scaling of its
facet areas. -- Julian M. 23-09-2013

@return The method returns the smallest global timestep for which the transport
distance in each element is less than or equal to its intersection.
*/
template<uint32_t dim>
double  NodeCenteredFiniteVolumeTransport<dim>::AnisotropicCourantIncrement()
{
    //if ( dim == 1U ) return CourantIncrement();

    UpdateProjectedVelocitiesAndFluxBalances();

    DenseMatrix<DM_MIN>  DN;
    vector<double>            grad(dim);
    const double              zero(0.);
    double                    courant_increment(8640000.); // 100 days
    size_t                      counter(0U);
    const bool                  multiply_with_cell_thickess = (thi_key_ == csmp::Index()) ? false : true;
    const bool                  velocity_multiplier         = (velo_mult_key_ == csmp::Index()) ? false : true;

    // diffusion is taken into account if the diffusion key is initialized
    const bool with_diffusion( (diff_key_ == csmp::Index()) ? false : true );

    // loop over the elements finding their transsect length in the direction of flow
#if !defined(_OPENMP)
    VectorVariable<dim>         vc;
    for ( typename vector<Element<dim>*>::const_iterator
          it=gref_.CellsBegin(); it!=gref_.CellsEnd(); it++ ){
        Element<dim> * eit = *it;
#else
    vector<double> thread_courant_increments(omp_get_max_threads());
    vector<uint32_t> thread_counters(omp_get_max_threads());
#pragma omp parallel
    {
        VectorVariable<dim>         vc;
        size_t tid = omp_get_thread_num();
        thread_courant_increments[tid]=courant_increment;
#pragma omp for
        for ( int32_t e = 0 ; e  < gref_.Cells(); e++ ){
            Element<dim>* eit = gref_.E(e);
            //----------------------------------------------------
            //change the element stencil to one for this thread, temporarily.
            FiniteElement* fe_tmp=eit->FE();
            // change pointer here
            eit->Assign(this->femgrs_[tid].E(eit->FE_Type()));
            //----------------------------------------------------
#endif

            (eit)->Read( vel_key_, vc );
            double       velocity(vc.Length());
            vc       /= velocity; // normalize vc to avoid round-off error during geometrical projection
            ScalarVariable poro, velo_mult;
            (eit)->PropertyValueAtBaryCenter( phi_key_, poro );
            double ediameter((eit)->LengthInDirection( vc ) );

            if (velocity_multiplier)
            {
                (eit)->PropertyValueAtBaryCenter( velo_mult_key_, velo_mult );
                velocity *= velo_mult();
            }
            /// Thickness should not be multiplied if velocity was assumed to be de-scaled.
            /// for legacy purposes, the transport class assumes that the supplied velocity has already been scaled
            /// to adjust for facet areas coming from the thickness attribute.  This is perhaps a little unintuitive,
            /// since any new user needs to know beforehand that this has to happen.  The user should always provide the
            /// "real" velocity (not scaled) and the transport class should deal with thickness through internal scaling of its
            /// facet areas.
            /// -- Julian M. 23-09-2013
            if ( multiply_with_cell_thickess ) ediameter *= (eit)->Read( this->thi_key_ );

            // limit due to diffusion-driven flow
            // ----------------------------------
            //if ( with_diffusion && adv1_key_.type == SCALAR) {
            //    fill( grad.begin(), grad.end(), zero );
            //    (*(eit)).dN_AtBaryCenter( DN );
            //    // grad transport variable-based formulation
            //    for ( size_t j{0U}; j<(eit)->Nodes(); j++ ) {
            //        double adv = (eit)->N(j)->Read( adv1_key_ );
            //        for ( size_t k{0U}; k<dim; k++ ) grad[k] += DN(k,j) * adv;
            //    }
            //    // getting maximum gradient value leading to the maximum diffusive flux
            //    double  magnitude_grad(grad[0]);
            //    if      ( dim == 2U ) magnitude_grad = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
            //    else if ( dim == 3U ) magnitude_grad = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);

            //    // using data from barycenter: O.K. as long as grad does not increase during iterations
            //    velocity += fabs( magnitude_grad * 2. * (eit)->Read( diff_key_ ) );
            //}

            // limit due to diffusion-driven flow
            // ----------------------------------
            if ( with_diffusion )
            {
                velocity +=  2. * (eit)->Read( diff_key_ ) / ediameter;
            }

            // guarding against degenerate cases
            if ( velocity > zero and ediameter > zero ) {
#if !defined(_OPENMP)
                courant_increment = std::min( courant_increment, ediameter*poro() / velocity );
                counter++;
#else
                thread_courant_increments[tid] = std::min( thread_courant_increments[tid], ediameter*poro() / velocity );
                thread_counters[tid]++;
#endif
            }

#if defined(_OPENMP )
            // put the FEM back here in case of openmp
            eit->Assign(fe_tmp);
#endif

        }
#if defined(_OPENMP )
    } // close the parallel pragma.
    for (auto i = 0 ; i < thread_counters.size();i++)
        counter+=thread_counters[i];

    courant_increment = *(std::min_element(thread_courant_increments.begin(),thread_courant_increments.end()));
#endif

    if ( counter == 0U )
        throw csmp::Exception( ERROR, "NodeCenteredFiniteVolumeTransport<dim>::AnisotropicCourantIncrement: ",
                               "anisotropic CFL could not be determined." );

    if (this->Verbose()) cout <<"\nNodeCenteredFiniteVolumeTransport<dim>::AnisotropicCourantIncrement: "<< courant_increment << endl;
    if ( courant_increment <= 1.0e-3 ) cout <<" \n (CFL constraint is very tight).";

    return courant_increment;

} // end AnisotropicCourantIncrement





/**

Method computes the minimum CFL criterion for the current model on the
basis of the volumetric troughput through each finite volume cell. This
is done taking into account that the velocity of the non-wetting phase
differs from the total velocity and that the inflow into a cell is
not necessarily equivalent to the outflow (as in the case of transient
flows or nodalsources or sinks).

The calculations are performed for the supplied relative permeability
model.

@section arguments Input Arguments

The method needs a reference to the current model in order to determine
the flow velocity distribution.

The specific relative permeability model that shall be used in the
computations is returned as second method argument.

@return The method returns the minimum global timestep for the current mesh
and velocity distribution for which the the transport distance in each
element is smaller or equal to its diameter.

@section implementation Implementation

The criteria are evaluated on an element by element basis using the
FV sector volumes to find CFL.

@section application Application

THIS IS RECOMMENDED FOR USE (1D VERSION) for an explicit transport
calculation. It checks everything: advection, diffusion and gravity-driven flow.

Overload this method if you use a less restrictive transport scheme than IMPES is.
*/
/*
 *
 *
 *
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::CourantIncrement( TwoPhaseModel<dim>& relperm )
{
    if ( dim != 1U )
        throw logic_error("NodeCenteredFiniteVolumeTransport<dim>::CourantIncrement: \
                          this method should only be applied in 1D");

    UpdateProjectedVelocitiesAndFluxBalances();

    static DenseMatrix<DM_MIN>  DN;
    double                    grad_psi;
    VectorVariable<dim>         vc;
    const double              hundred_days(8640000.);
    double                    velocity,
    courant_increment(hundred_days);

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
    {
        relperm.Initialize( *(*eit) );
        relperm.InitializeForBaryCenter( *(*eit) );
        relperm.EffectiveSaturation();

        // 1. limit imposed by advection
        // -----------------------------
        (*eit)->Read( vel_key_, vc );
        velocity = vc.Length() * relperm.MaxFractionalFlowDerivative();

        // 2. limit due to buoyancy-related flow
        // -------------------------------------
        velocity += fabs( relperm.GravityMultiplier_G() );

        // 3. limit due to capillary-driven flow
        // -------------------------------------
        const double k_lambda_overbar(relperm.Permeability() * relperm.G());
        if ( k_lambda_overbar > numeric_limits<double>::epsilon() ) {
            (*(*eit)).dN_AtBaryCenter( DN );
            // computing the saturation gradient
            grad_psi = 0.;
            for ( size_t j{0U}; j<(*eit)->Nodes(); j++ )
                grad_psi += DN(0U,j) * (*eit)->N(j)->Read( adv1_key_ );

            // limit is related to element volume and saturation gradient
            // this stability limity F*t < dT/(dx)^2 is explained in Leveque, p. 60
            velocity += k_lambda_overbar * relperm.dpcds_Phase(2U) * (grad_psi * grad_psi);
        }

        // 4. calculating the CFL criterion from the cell diameter for each sector (pore) volume
        // -------------------------------------------------------------------------------------
        for ( auto i{0U}; i<(*eit)->Nodes(); i++ )     // porevolume is taken as proxy of sector length
            courant_increment = std::min( courant_increment, STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i) / velocity );
    }

    if ( courant_increment >= hundred_days ) {
        cout <<"\nNodeCenteredFiniteVolumeTransport::CourantIncrement (2-phase flow): could not be determined."<< endl;
        cout <<"\nthere appears to be no flow in the model domain, returning CFL=3600."<< endl;
        return hundred_days;
    }

    cout <<"\nNodeCenteredFiniteVolumeTransport::CourantIncrement (2-phase flow): "<< courant_increment << endl;
    if ( courant_increment <= 1.0e-3 ) cout <<"\n\tCFL constraint is very tight: "<< courant_increment << endl;

    return courant_increment;

} // end CourantIncrement (2-phase flow - all)

*/


/**

Computes the Courant time increment (CFL citerion) taking into account
viscous, gravitational and capillary fluid displacements using the
contraints from the provided relative permeability model. The CFL
criterion is calculated using the element diameter in the direction
of the flow but not account for the deviation from this vector<double> of
the flow of the considered phase.

For the viscous flow the shock speed is used as a multiplier for the
total velocity. This may lead to a too tight constraint but is safe
for the case where CFL is computed only once at the onset of series
of advection steps with an explicit scheme.

@section arguments Input Arguments

First, a const reference to the current model, a reference to the
TwoPhaseModel used (this will be a subclass), and a maximum
time-increment imposed by external constraints.

@return  The CFL criterion.

@section implementation Implementation

The criterion is computed sector by sector for all the finite
elements in the entire model.

@section application Application

Method has been designed primarily to give stability to an explicit
transport scheme but it also serves as a guide for the timestepping
using an implicit approach.

@section messages Messages

The method alerts the user to the case where CFL is larger than the
input max_time_increment which may be the case if there is no flow at all
in the domain. Equally, the user is informed if the CFL increment is
less than a millisecond (usually a prohibitively small increment).
*/
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::AnisotropicCourantIncrement( TwoPhaseModel<dim>& relperm,
                                                                              double max_time_increment )
{

    //if ( dim == 1U ) return CourantIncrement( relperm);

    UpdateProjectedVelocitiesAndFluxBalances();

    DenseMatrix<DM_MIN>  DN;
    vector<double>     gradPc(dim);
    VectorVariable<dim>  vc;
    double             velocity, courant_increment(max_time_increment);
    const double       millisecond(1.0e-3);
    const bool           multiply_with_cell_thickess = (thi_key_ == csmp::Index()) ? false : true;

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
    {
        relperm.Initialize( *(*eit) );
        relperm.InitializeForBaryCenter( *(*eit) );
        relperm.EffectiveSaturation();

        // 1. limit imposed by advection of the fluid mixture
        // --------------------------------------------------
        (*eit)->Read( vel_key_, vc );
        velocity = vc.Length();
        // NB: This may be a too conservative estimate for the implicit scheme but is necessary for the explicit one
        velocity = std::max( velocity, velocity * relperm.MaxFractionalFlowDerivative() );
        double cell_diameter = (*eit)->LengthInDirection(vc) * (*eit)->Read( phi_key_ );
        if ( multiply_with_cell_thickess ) cell_diameter *= (*eit)->Read( this->thi_key_ );


        // 2. limit due to buoyancy-driven flow in any zone of intermediate saturations
        // ----------------------------------------------------------------------------
        // alternatively use: k * max(lambda_overbar) * delta_rho g  (max(lambda_overbar)=Gmax is achieved at around sn=0.5)
        velocity += fabs( relperm.GravityMultiplier_G() );


        // 3. limit due to potential capillary spreading
        // ---------------------------------------------
        const double k_lambda_overbar(relperm.Permeability() * relperm.G());
        // computing the capillary pressure gradient
        if ( k_lambda_overbar > numeric_limits<double>::epsilon() )
        {
            fill( gradPc.begin(), gradPc.end(), 0. );
            (*(*eit)).dN_AtBaryCenter( DN );
            for ( auto j{0U}; j<(*eit)->Nodes(); j++ ) {
                double sn = (*eit)->N(j)->Read( adv1_key_ );
                relperm.SaturationWettingPhase( 1. - sn );
                relperm.EffectiveSaturation();
                double pc = relperm.pc_Phase( );
                for ( auto k{0U}; k<dim; k++ ) gradPc[k] += DN(k,j) * pc;
            }
            // getting the maximum capillary flux (G= lambda overbar)
            double  magnitude_grad_pc(gradPc[0]); // 1D
            if      ( dim == 3U ) magnitude_grad_pc = sqrt(gradPc[0]*gradPc[0]+gradPc[1]*gradPc[1]+gradPc[2]*gradPc[2]);
            else if ( dim == 2U ) magnitude_grad_pc = hypot(gradPc[0],gradPc[1]);

            // use data from barycenter: O.K. as long as grad_pc does not increase during iterations
            velocity += fabs( magnitude_grad_pc * k_lambda_overbar );
        }

        // 4. calculating the CFL criterion from the cell diameter
        // -------------------------------------------------------
        // guarding against degenerate cases
        if ( velocity > numeric_limits<double>::epsilon() and cell_diameter > numeric_limits<double>::epsilon() )
            courant_increment = std::min( courant_increment, cell_diameter / velocity );
    }

    if ( courant_increment >= max_time_increment ) {
        throw csmp::Exception( INFO, "NodeCenteredFiniteVolumeTransport::AnisotropicCourantIncrement (2-phase flow)",
                               "calculated courant increment is larger than maximum permitted increment, there may be no flow in the model domain" );
        cout <<"\nCFL is set to "<< max_time_increment << endl;
        return max_time_increment;
    }

    if ( courant_increment <= millisecond ) {
        throw csmp::Exception( WARNING, "NodeCenteredFiniteVolumeTransport::AnisotropicCourantIncrement (2-phase flow)",
                               "calculated courant increment is smaller than a millisecond" );
    }

    cout <<"\nNodeCenteredFiniteVolumeTransport::AnisotropicCourantIncrement (2-phase flow): "<< courant_increment << endl;
    return courant_increment;

} // end AnisotropicCourantIncrement (2-phase flow - all)




template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::CFL_Multiplier( double desired_value )
{
    cfl_multiplier_ = desired_value;
}



template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::CFL_Multiplier() const
{
    return cfl_multiplier_;
}

template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::WithLsmGradientLimiter(Model<dim>& sg)
{
    with_lsmgrad_limiter_=true;

    grad_advprop_limiter_=new GenericNodePropertyGradientLimiter<dim>( sg, gref_.Name().c_str(), advected_variable_.c_str() );
    mass_center_key_=sg.Database().StorageKey("mass center");
    grad_advprop_key_=sg.Database().StorageKey((advected_variable_+std::string(" gradient")).c_str());
    grad_advprop_limiter_key_=sg.Database().StorageKey((advected_variable_+std::string(" limiter")).c_str());
}


/**

For the finite volume corresponding to the node supplied via a pointer.
The flow across the model boundary is calculated together with its
flux balance so that this can be split of as a source term if
necessary.

The flux balance is calculated for the stencils that the node is connected to
and which give respective sector contributions for the finite volume of
interest. Distinctions are made with regard to whether stencils are inside
the group for which advection is calculated or outside.

To get model in and outfluxes the interior stencils are used.

@return The method returns the inflow into the model domain and whether the flux
balance could be evaluated as indicated by a bool.

The convention adopted here is that inflow is positive while outflow
is negative since the normals point out of the finite volume.

The flux_balance is returned provided that it can be evaluated.
This is only possible if the finite volume lies inside the model domain so
that it is not truncated. In this case the method returns true, else
false.

The flux balance is equal to the divergence of the flux for the current
finite volume (if there is more outflow, the flux balance and divergence
are greater than zero).

@section application Application

Use this method only on finite volumes that are located on the boundary
of the computational domain.

 */
template<uint32_t dim>
bool NodeCenteredFiniteVolumeTransport<dim>::FluxThroughBoundaryFiniteVolume(
        const Node<dim>* nd_ptr,
        double& inflow, double& flux_balance ) const
{
    inflow = flux_balance = static_cast<double>(0.);
    VectorVariable<dim>  vel;

    // for all those sectors of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
    for ( auto t=0U; t<nd_ptr->Parents(); t++ ) {
        //if(nd_ptr->Parent(t)->Idx() > STENCIL_DATA.size())
        //  cout<<"\n **";
        assert( nd_ptr->Parent(t) != nullptr );
        const auto nid(nd_ptr->ParentNodeNumber(t));
        double  flux(0.);
        // for all facets surrounding the finite volume at the boundary
        for ( auto i{0U}; i<nd_ptr->Parent(t)->FV()->FacetsPerSector(nid); i++ )
        {
            auto iFacet( nd_ptr->Parent(t)->FV()->FacetSurroundingSector(nid,i) );
            // fluxes are determined for the sectors inside and outside of the advection region
            nd_ptr->Parent(t)->Read( vel_key_, vel );
// SKM FIX - ignoring Face objects
if ( dim == 2U and nd_ptr->Parent(t)->IsLine() ) continue;
          
            if ( nid == nd_ptr->Parent(t)->FV()->InsideNode( iFacet ) )
                flux += STENCIL_DATA[ nd_ptr->Parent(t)->Idx() ].FacetArea(iFacet) *
                        STENCIL_DATA[ nd_ptr->Parent(t)->Idx() ].FacetNormalProjection( iFacet, vel );
            else
                flux -= STENCIL_DATA[ nd_ptr->Parent(t)->Idx() ].FacetArea(iFacet) * // facet normal velocity
                        STENCIL_DATA[ nd_ptr->Parent(t)->Idx() ].FacetNormalProjection( iFacet, vel );
        }
        // inflow and outflow are measured using the stencils in the interior of the computational region
        // thus inflows originate as positive and outflows as negative
        // -----------------------------------------------------------
        // if this is an interior stencil (must be if node is on the model boundary)
        if ( gref_.IsPerimeterNode( nd_ptr ) ) inflow += flux; // +to satisfy convention above (that inflow is positive)
        // the balance can  be evaluated if there is a halo stencil
        else flux_balance -= flux;
    }

    // give an indication whether the flux balance was evaluated
    if ( gref_.IsPerimeterNode( nd_ptr ) ) return false;

    // the flux balance can be evaluated
    return true;
}


/* slow code cut out
              if ( nid == nd_ptr->Parent(t)->FV()->InsideNode( iFacet ) )
                flux += ( *nd_ptr->Parent(t) ).FacetArea(iFacet) * // facet normal velocity
                        ( *nd_ptr->Parent(t) ).ProjectionOnFacetNormal( iFacet, vel_key_ );
              else
                flux -= ( *nd_ptr->Parent(t) ).FacetArea(iFacet) * // facet normal velocity
                        ( *nd_ptr->Parent(t) ).ProjectionOnFacetNormal( iFacet, vel_key_ );
*/




/**

Compensates for inflow and outflow at model boundaries where the finite
volumes are truncated.

NB: This method does not compensate for diffusive fluxes which require
the explicit specification of boundary conditions by the user.

@section arguments Input Arguments

Apart from references to the instances of the current model and
transport algorithm, this method allows the user to specify whether
the saturation from the first time level shall be used (this is required
if second-order accuracy in time is desired) and whether it is permitted
to treat in and outfluxes implicitly. While this is the preferred option,
it may lead to a negative matrix diagonal.

@section implementation Implementation

Method uses FluxThroughBoundaryFiniteVolume() to evaluate in and
outfluxes of the domain. This implies that only 1st-order fluxes will
be considered since because the transport variable is not interpolated
to any of the facets of the halo stencil.

A rationale for this comes from the frequent situation that the
transported variable is not expected to be continuous across the
domain boundary.

@section application Application

This method is used in implicit advection calculations.
 */
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AssignFluxBoundaryConditions( NodeCenteredFiniteVolumeAlgorithm<dim>& advection_algorithm,
                                                                           bool use_t0_saturation ) const
{
    double        inflow, flux_balance;
    const double  zero(0.);

    // loop over the boundary cells and adjust fluxes
    for ( size_t i=gref_.InteriorNodes(); i<gref_.Nodes(); i++ )
    {
        // if we have a flux balance from the part of the boundary finite volume which extends beyond the group boundary
        // the balance has the opposite sign than the inflow. Thus inflow is associated with a negative balance.

        // if flux balance cannot be evaluated because we are at the model boundary
        if ( !FluxThroughBoundaryFiniteVolume( gref_.N(i), inflow, flux_balance ) ) {
            // inflow compensation
            if ( inflow > zero ) {
                if ( use_t0_saturation )
                    advection_algorithm.AddToRHS( i, SAT0[i] * inflow );
                else
                    advection_algorithm.AddToRHS( i, gref_.N(i)->Read( adv1_key_ ) * inflow );
            }
            // outflow compensation
            else if ( inflow < zero )
                // no minus sign, since this contribution is added to the LHS
                advection_algorithm.AddToLHS( i, i, -inflow );
        }
        // if the flux balance can be evaluated
        else {
            if ( flux_balance > zero ) {
                if ( use_t0_saturation )
                    advection_algorithm.AddToRHS( i, SAT0[i] * flux_balance );
                else
                    advection_algorithm.AddToRHS( i, gref_.N(i)->Read( adv1_key_ ) * flux_balance );
            }
            else if ( flux_balance < zero )
                advection_algorithm.AddToLHS( i, i, -flux_balance );
        }
    }

} // end AssignBoundaryConditions




/**

Checks volume conservation by FV scheme. Boundary flux mismatches are
used elsewhere to impose zero order Neumann boundary conditions at inflow and
outflow boundaries.

@section arguments Input Arguments
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::RecordFluxBalances( const char* balance_var,
                                                                 double& total_surplus,
                                                                 double& total_deficit ) const
{
    cout <<"\nNodeCenteredFiniteVolumeTransport<"<< dim;
    cout <<">::SolveAdvectionEquationImplicitly: Testing flux balance..."<< endl;
    const double            zero(0.);
    ScalarVariable  sc;
    csmp::Index         var_key_ = pref_.StorageKey(balance_var);
    string              b;
    total_surplus = total_deficit = zero;

    for ( auto n=0U; n<gref_.Nodes(); n++ ) {
        if ( FLUX_BALANCE[n] < zero ) total_deficit += fabs(FLUX_BALANCE[n]);
        else                          total_surplus += FLUX_BALANCE[n];

        if ( gref_.N(n)->AtBoundary() != NOT ) {
            cout <<"\nNode<"<< dim <<"> "<< gref_.N(n)->Idx();
            cout <<" at boundary "<< (b=parseBoundary(gref_.N(n)->AtBoundary()));
            cout <<":  flux balance: "<< FLUX_BALANCE[n];
        }
        // recording the flux mismatch. The model boundaries were not calculated previously
        gref_.N(n)->Store( var_key_, sc=FLUX_BALANCE[n] );
    }

} // end RecordFluxBalances




/** Computes volume fluxes into and out of the model.

The current total inflow and outflow are computed and output to the first
and second method argument, respectively.

@return returns difference between inflow and outflow: flux-balance (|inflow| - |outflow| for entire model)

@attention when the thickness attribute is used (i.e., when lower-dimensional elements are incorporated
into the simulation, this method expects the thickness-weighted velocity as input;
else it will not return the correct flux.

@attention SKM 27/5/14: modified to get correct fluxes using 2 approaches: 1) box-boundary flagging, and 2)
based on perimeter nodes. The former is to be deprecated once the perimeter flagging for the boundary
works.

@test SKM 29/5/14 refactored, and tested for box-boundary case.
      Roman, 14/10/14 corrected local facet number
*/
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::BoundaryFluxes( double& in_flow,
                                                                 double& out_flow,
                                                                 bool box_shaped,                           ///< special case (to be deprecated)
                                                                 bool use_advected_variable,                ///< use advected variable or not?
                                                                 const char* advected_variable,             ///< advected variable
                                                                 const char* variable_to_determine_no_flow, ///< variable to determine no-flow boundaries
                                                                 VARIABLE_FLAG exclude_flag                 ///< variable flag for node/FV to be excluded
                                                               ) const
{
    double inflow(0.), outflow(0.);
    VectorVariable<dim>  vc;

    csmp::Index var_key             = mref_.Database().StorageKey(variable_to_determine_no_flow);

    csmp::Index     advected_var_key;
    ScalarVariable  advected_var;
    if( use_advected_variable )
    {
        advected_var_key = mref_.Database().StorageKey(advected_variable);

        if ( advected_var_key.place != NODE && advected_var_key.type != SCALAR ){
            string  msg ="The advected variable '";
            msg       += pref_.Name( advected_var_key );
            msg       +="' must be a scalar-type node variable";
            throw csmp::Exception( FATAL_ERROR,
                                   "NodeCenteredFiniteVolumeTransport<dim>::BoundaryFluxes",
                                   msg.c_str() );
        }
    }

    // 1. case box shaped model: for all nodes at model boundary velocity is projected on the facets of parent elements
    // ----------------------------------------------------------------------------------------------------------------
    // (influxes are accrued as positive contributions because they are aligned with the facets of boundary cells,
    //  outgoing ones as negative)
    if ( box_shaped ) {
        cout <<"\nNodeCenteredFiniteVolumeTransport<dim>::BoundaryFluxes: carrying out calculation assuming a box-shaped model.\n";
        cout.flush();
        for ( typename std::vector<Node<dim>*>::const_iterator
              nit=gref_.NodesBegin(); nit!=gref_.NodesEnd(); ++nit )
          if ( (*nit)->AtBoundary() != NOT  and (*nit)->Status(var_key) != exclude_flag )
            {  // for the adjacent sectors of the node's parent elements
               double finite_volume_influx(0.);
               if( use_advected_variable )
                   (*nit)->Read( advected_var_key, advected_var );
               for ( auto i=0U; i<(*nit)->Parents(); ++i ) {
                    auto nd = (*nit)->ParentNodeNumber(i);
                    // the velocity is expected to already take potential thickness attributes into account
                    (*nit)->Parent(i)->Read( vel_key_, vc );
                    // loop over all facets j, integrating the velocity over their area
                    for ( auto j{0U}; j<(*nit)->Parent(i)->FV()->FacetsPerSector(nd); ++j ) {
                         const auto facet((*nit)->Parent(i)->FV()->FacetSurroundingSector(nd,j));
                         double facet_flux = (*nit)->Parent(i)->ProjectionOnFacetNormal( facet, vc );
                         facet_flux *= (*nit)->Parent(i)->FacetArea(facet);
                         // distinguishing 2 cases:
                         auto inside_node = (*nit)->Parent(i)->FV()->InsideNode(facet);
                         //   1. outward point normal (in this case a positive flux indicates outflow)
                         if ( inside_node == nd )
                             finite_volume_influx += facet_flux;
                         //   2. inward pointing normal
                         else
                             finite_volume_influx -= facet_flux;
                      }
                 }
               if( use_advected_variable )
                   finite_volume_influx *= advected_var();
               // if we are dealing with a FV marking an inflow boundary
               if ( finite_volume_influx > 0. ) inflow  += finite_volume_influx;
               else                             outflow += fabs(finite_volume_influx);
            }
      
         // output assignment
         // (NB: given a divergence free velocity field, FVs that record outflow must lie an inflow
         //      model boundary as such FVs have no boundary facets. The opposite is true for inflow-dominated FVs).
         in_flow  = fabs(outflow);
         out_flow = fabs(inflow);
      
         return in_flow - out_flow;
      }
  
    // 2. For an arbitrarily shaped model with a correctly flagged perimeter
    // ---------------------------------------------------------------------
    for ( typename std::vector<Node<dim>*>::const_iterator
          nit=gref_.PerimeterNodesBegin(); nit!=gref_.NodesEnd(); ++nit )
      if ( (*nit)->AtBoundary() != NOT and (*nit)->Status(var_key) != exclude_flag )
        {  // for the adjacent sectors of the node's parent elements
           double finite_volume_influx(0.);
           if( use_advected_variable )
               (*nit)->Read( advected_var_key, advected_var );
           for ( auto i=0U; i<(*nit)->Parents(); ++i ) {
                auto nd = (*nit)->ParentNodeNumber(i);
                (*nit)->Parent(i)->Read( vel_key_, vc );
                // loop over all facets j, integrating the velocity over their area
                for ( auto j{0U}; j<(*nit)->Parent(i)->FV()->FacetsPerSector(nd); ++j ) {
                     const auto facet((*nit)->Parent(i)->FV()->FacetSurroundingSector(nd,j));
                     double facet_flux = (*nit)->Parent(i)->ProjectionOnFacetNormal( facet, vc );
                     facet_flux *= (*nit)->Parent(i)->FacetArea(facet);
                     // distinguishing 2 cases:
                     size_t inside_node = (*nit)->Parent(i)->FV()->InsideNode(facet);
                     //   1. outward point normal (in this case a positive flux indicates outflow)
                     if ( inside_node == nd )
                         finite_volume_influx += facet_flux;
                     //   2. inward pointing normal
                     else
                         finite_volume_influx -= facet_flux;
                  }
             }
           if( use_advected_variable )
               finite_volume_influx *= advected_var();
           // if we are dealing with a FV marking an inflow boundary
           if ( finite_volume_influx > 0. ) inflow  += finite_volume_influx;
           else                             outflow += fabs(finite_volume_influx);
        }
  
     // output assignment
     // (NB: given a divergence free velocity field, FVs that record outflow must lie an inflow
     //      model boundary as such FVs have no boundary facets. The opposite is true for inflow-dominated FVs).
     in_flow  = fabs(outflow);
     out_flow = fabs(inflow);
  
     return in_flow - out_flow;

} // end BoundaryFluxes




/** Computes and returns inflow into model or subregion of it.
*/
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::ModelInflow() const
{
    double in_flow, out_flow;
    BoundaryFluxes( in_flow, out_flow );
    return in_flow;
}


/** Computes and returns outflow into model or subregion of it.
 */
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::ModelOutflow() const
{
    double in_flow, out_flow;
    BoundaryFluxes( in_flow, out_flow );
    return out_flow;
}



/** Second-order single-phase transport.

@section arguments Input Arguments

If the model is transient so that there are poro-elastic sources or
sinks, these can be dealt with via a correction that uses the
FLUX_BALANCE for each finite volume.

@todo (1) The SAMG solver should be created one time or the SAMG_Setting should be modified,
otherwise it doesn't converge.
*/
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable( double time_interval,
                                                                 double cfl_multiplication_factor,
                                                                 bool apply_flux_balance_correction,
                                                                 bool update_pore_volumes)
{
    gref_.RenumberNodes();

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( !apply_flux_balance_correction )
        csmp_error.notice( WARNING, "NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable",
                           "without flux balance correction boundary condition assignment may not be suitable; watch for erratic concentrations");

    // 0. backup of fluxes from last time level
    // ----------------------------------------
    if ( SecondOrderInTime() ) BackupFacetFluxes();

    // 1. update the sector volumes and FVPOREVOL if so required (as in heat transport for instance)
    // ---------------------------------------------------------------------------------------------
    if ( update_pore_volumes ) {
        InitializeFiniteVolumeData();
        InitializeArraysForFirstOrderMethod();
    }

    // 2. compute the CFL condition to identify value for overstepping
    // ---------------------------------------------------------------
    const double  courant_increment(AnisotropicCourantIncrement());
    double        time(0.), time_increment = cfl_multiplication_factor * courant_increment;

    // 3. compute solution
    // -------------------

    if (firstCall_){
        baseAdvector_= new NodeCenteredFiniteVolumeAlgorithm<dim>( gref_ );
        firstCall_ = false;
    }

    if (Verbose()) cout <<"\n\n\nNodeCenteredFiniteVolumeTransport<"<< dim;
    if (Verbose()) cout <<">::AdvectVariable: Advecting transport variable";
    while ( time < time_interval ) {
        if ( (time_interval - time) < 2.0*time_increment ) time_increment = time_interval - time;

        if ( SecondOrderInSpace() ) {
            if ( SecondOrderInTime() )
                AdvectVariable2ndOrderInSpaceAndTime( (*baseAdvector_), time_increment, apply_flux_balance_correction );
            else
                AdvectVariable2ndOrder( (*baseAdvector_), time_increment, apply_flux_balance_correction );
        }
        else AdvectVariable1stOrder( (*baseAdvector_), time_increment, apply_flux_balance_correction );

        time += time_increment;
    }

    // 4. check results
    // ----------------
#ifdef DEBUG_NodeCenteredFiniteVolumeTransport
    double surplus, deficit;
    RecordFluxBalances( "nodal flux mismatch", surplus, deficit );
    double prop_min, prop_max;
    gref_.MinMaxOf( advected_variable_.c_str(), prop_min, prop_max );
    cout <<"\nNodeCenteredFiniteVolumeTransport<"<<  dim;
    cout <<">::AdvectVariable: The range of the advected ";
    cout <<"property after advection is: "<< prop_min <<" to "<< prop_max << endl;
    if ( surplus > 1.0e-17 || deficit > 1.0e-17 )
        cout <<"\nNodeCenteredFiniteVolumeTransport<"<<  dim
            <<">::AdvectVariable: total flux-surplus and defizit: "<< surplus <<" "<< deficit << endl;
#endif

    return courant_increment;

} // end AdvectVariable

/** Single phase transport, single step
@author Julian E. Mindel

@section arguments Input Arguments

This method calculates a single timestep of transport, assumming nothing about the time interval
it is "fed".  It will advect using that time interval, assuming the user has externally determined
that this is the correct time interval length.  In contrast with AdvectVariable(, this method does
not break the time interval into sub-parts and guarantee stability. This also means that the courant
increment is not checked for in this method (and should be, externally, of course).

If the model is transient so that there are poro-elastic sources or
sinks, these can be dealt with via a correction that uses the
FLUX_BALANCE for each finite volume.

@todo (1) The SAMG solver should be created one time or the SAMG_Setting should be modified,
otherwise it doesn't converge.

*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AdvectVariableSingleStep( double time_increment,
                                                                       bool apply_flux_balance_correction,
                                                                       bool update_pore_volumes)
{
    gref_.RenumberNodes();

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( !apply_flux_balance_correction )
        csmp_error.notice( WARNING, "NodeCenteredFiniteVolumeTransport<dim>::AdvectVariableSingleStep",
                           "without flux balance correction boundary condition assignment may not be suitable; watch for erratic concentrations");

    // 0. backup of fluxes from last time level
    // ----------------------------------------
    if ( SecondOrderInTime() ) BackupFacetFluxes();

    // 1. update the sector volumes and FVPOREVOL if so required (as in heat transport for instance)
    // ---------------------------------------------------------------------------------------------
    if ( update_pore_volumes ) {
        InitializeFiniteVolumeData();
        InitializeArraysForFirstOrderMethod();
    }

    // 2. compute solution
    // -------------------

    if (firstCall_){
        baseAdvector_= new NodeCenteredFiniteVolumeAlgorithm<dim>( gref_ );
        firstCall_ = false;
    }

    if (Verbose()) cout <<"\n\n\nNodeCenteredFiniteVolumeTransport<"<< dim;
    if (Verbose()) cout <<">::AdvectVariableSingleStep: Advecting transport variable.";


    if ( SecondOrderInSpace() ) {
        if ( SecondOrderInTime() )
            AdvectVariable2ndOrderInSpaceAndTime( (*baseAdvector_), time_increment, apply_flux_balance_correction );
        else
            AdvectVariable2ndOrder( (*baseAdvector_), time_increment, apply_flux_balance_correction );
    }
    else AdvectVariable1stOrder( (*baseAdvector_), time_increment, apply_flux_balance_correction );

    // 3. check results
    // ----------------
#ifdef DEBUG_NodeCenteredFiniteVolumeTransport
    double surplus, deficit;
    RecordFluxBalances( "nodal flux mismatch", surplus, deficit );
    double prop_min, prop_max;
    gref_.MinMaxOf( advected_variable_.c_str(), prop_min, prop_max );
    cout <<"\nNodeCenteredFiniteVolumeTransport<"<<  dim;
    cout <<">::AdvectVariable: The range of the advected ";
    cout <<"property after advection is: "<< prop_min <<" to "<< prop_max << endl;
    if ( surplus > 1.0e-17 || deficit > 1.0e-17 )
        cout <<"\nNodeCenteredFiniteVolumeTransport<"<<  dim
            <<">::AdvectVariable: total flux-surplus and defizit: "<< surplus <<" "<< deficit << endl;
#endif

} // end AdvectVariableSingleStep


/** Stub for method which needs to be implemented in subclass.
 */
// function stub for single-phase flow base class
template<uint32_t dim>
double  NodeCenteredFiniteVolumeTransport<dim>::TransportPhase( TwoPhaseModel<dim>& ff, double )
{
    cout <<"\nNodeCenteredFiniteVolumeTransport<"<<  dim;
    cout <<">::TransportPhase: you called virtual method stub which cannot ";
    cout <<" be used for two-phase flow. You need to call specific subclasses in stead."<< endl;
    ff.Out();
    gref_.Nodes();
    return std::numeric_limits<double>::signaling_NaN();
}



/** First order implicit transport scheme for the advection-only of passive
tracers.

@section arguments Input Arguments

A flux-balance correction can be applied to correct for poro-elastic
fluid sources and sinks.

@section application Application

Used by the method AdvectVariable().
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable1stOrder( NodeCenteredFiniteVolumeAlgorithm<dim>& advector,
                                                                     double time_increment,
                                                                     bool with_flux_balance_correction )
{
    const bool  with_diffusion = (diff_key_ == csmp::Index())  ? false : true;
    if (var_ncomponents_==0)
        throw csmp::Exception (FATAL_ERROR,"NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable1stOrder",
                              "Type of the advected variable is not known or number of variable components not initialized properly.  Will not advect!\
                              Check the initialization process of your NCFVT class.");

    for (size_t ncom = 0; ncom<var_ncomponents_; ncom++ ){
        vector<FV_Parameter>::const_iterator  fvt(STENCIL_DATA.begin());

        for ( typename vector<Element<dim>*>::const_iterator
              eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++, fvt++ )
        {
            // 2.1 getting all the necessary information from each element (second order=false)
            stencil_.InitializeFirstOrder( (*fvt), *(*eit) ,adv1_key_.type, ncom );

            // 2.2 first-order concentrations are accumulated into lefthand side
            advector.AccumulateLHS( stencil_, time_increment );

            // 2.3 source terms due a divergence of vt
            if ( with_flux_balance_correction ) 
            {
                stencil_.ComputeBoundaryFluxMismatch((*fvt), *(*eit));
                advector.AccumulateSectorSourceTermsInLHS( stencil_ );
            }
            // 2.4 conductance matrix for diffusion
            if ( with_diffusion ) advector.AccumulateIntegral_DNT_op_DN_dV_LHS( stencil_ );

            // 2.5 previous solution multiplied by storage and divided by time increment is
            //     accumulated into righthandside for Backward-Euler time stepping
            advector.AccumulateRHS( stencil_, time_increment );

        } // end of accumulation

        // 2.6 algebraic multigrid solver is applied to compute FV saturations
        advector.SolveMatrixEquation();

        // 2.7 Saving the computed new saturations at the FV centers
        advector.OutputResults( pref_, adv1_key_, true, ncom );

        advector.ResetLHS( gref_.Nodes() );
        advector.ResetRHS( gref_.Nodes() );
    }

} // end AdvectVariable1stOrder



/**

See Matthai et al. 2009, TIPM for description of 2nd-order in space
algorithm implemented here.
@todo (1) There're two flux balance corrections in item 2.3 and item 2.8. Should be only one correction. By the moment the item 2.8 is commented.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable2ndOrder(
        NodeCenteredFiniteVolumeAlgorithm<dim>& advector,
        double time_increment,
        bool with_flux_balance_correction )
{
    // diffusion is taken into account if the diffusion key is initialized
    const bool with_diffusion( (diff_key_ == csmp::Index()) ? false : true );

    // 1.0 finding smin/smax of trial solution in the neighborhood of each element

    if (adv1_key_.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR,"NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable2ndOrder"," Second order advection does not support another variable type other than SCALAR.");

#ifdef CSMP_WITH_SAMG_SOLVER
    baseAdvector_->GetSolverSettings().Set_iswit(4);
#else
    /// add extra functionality for alternative solver if needed
#endif

    if(!with_lsmgrad_limiter_){

        size_t iter(1U);
        double res(numeric_limits<double>::max());

        //Initialization section
        double smin(0.0), smax(0.0);
        this->SAT0.resize( this->gref_.Nodes() );
        vector<double>(this->SAT0).swap(this->SAT0);
        this->InitialAdvectedPropertyValues(this->adv1_key_,smin,smax);
        //this->InitialAdvectedPropertyValues(this->ad1_key_);

        // First iteration:
        // 1. empty the sparse matrix and righthand vector
        advector.ResetLHS( this->gref_.Nodes() );
        advector.ResetRHS( this->gref_.Nodes() );

        MinMaxAdvectedProperty();

        for ( typename vector<Element<dim>*>::const_iterator eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
        {
            // 2.0 global to local ID conversion
            stencil_.eidx_ = (*eit)->Idx();

            // 2.1 during FIRST ITERATION: computing facet fluxes and interpolating saturations in FV sectors
            stencil_.InitializeSecondOrder( STENCIL_DATA[stencil_.eidx_], *(*eit) ); // O.K.

            // 2.2 FV volumes/dt and first-order fluxes are accumulated into lefthand side (matrix)
            advector.AccumulateLHS( stencil_, time_increment ); // O.K.

            // 2.3 source terms due a divergence of vt
            if ( with_flux_balance_correction ) 
            {
                stencil_.ComputeBoundaryFluxMismatch( STENCIL_DATA[stencil_.eidx_], *(*eit) );
                advector.AccumulateSectorSourceTermsInLHS( stencil_ );
            }

            // 2.4 conductance matrix for diffusion
            if ( with_diffusion ) advector.AccumulateIntegral_DNT_op_DN_dV_LHS( stencil_ );
            // alternative formulation for diffusion
            //if ( with_diffusion ) advector.AccumulateIntegral_DN_op_dS_LHS( stencil_ );

            // 2.5  update and compute values of advected variable at faces (spatial limiting)
            stencil_.IsotropicallyLimitTransportProperties(*(*eit), SMINMAX);

            // 2.6 account for the storage
            advector.AccumulateRHS( stencil_, SAT0, time_increment );

            // 2.7 add first-order Backward Euler solution to RHS
            //     the latest first-order solution is used here
            advector.AddToRHS( stencil_);

            // 2.8 previous solution multiplied by storage and divided by time increment is
            //     accumulated into righthandside for Backward-Euler time stepping
            advector.AccumulateHigherOrderRHS( stencil_ );

        }

        // 2.9 algebraic multigrid solver is applied to compute FV saturations
        advector.SolveMatrixEquation();

        // 2.10 Saving the computed new saturations at the FV centers
        if ( iter == this->max_nonlinear_limiting_case_iterations_ )
            res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, true );  // O.K.
        else
            res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, false ); // O.K.

        // measuring convergence
        //cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< res << endl;

        if ( res > this->target_nonlinear_limiting_case_residual_){

            //baseAdvector_->GetSolverSettings().Set_iswit(1);

            iter+=1;

            while( (res > this->target_nonlinear_limiting_case_residual_) && (iter < this->max_nonlinear_limiting_case_iterations_)){

                // 3. empty the righthand vector
                advector.ResetRHS( this->gref_.Nodes() );

                MinMaxAdvectedProperty();

                for ( typename vector<Element<dim>*>::const_iterator
                      eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ ){
                    // 3.0 global to local ID conversion
                    stencil_.eidx_ = (*eit)->Idx();

                    // 3.1 computing upstream and interpolating advected property values to FV segments
                    stencil_.InitializeAdvectedVariableValues( *(*eit) ); // O.K.

                    // 3.2  update and compute values of advected variable at faces (spatial limiting)
                    stencil_.IsotropicallyLimitTransportProperties(*(*eit), SMINMAX);

                    // 3.3 account for the storage
                    advector.AccumulateRHS( stencil_, SAT0, time_increment );

                    // 3.4 add first-order Backward Euler solution to RHS
                    //     the latest first-order solution is used here
                    advector.AddToRHS( stencil_);

                    // 3.5 previous solution multiplied by storage and divided by time increment is
                    //     accumulated into righthandside for Backward-Euler time stepping
                    advector.AccumulateHigherOrderRHS( stencil_ );

                } // end of accumulation


                // 3.6 algebraic multigrid solver is applied to compute FV saturations
                advector.SolveMatrixEquation();

                // 3.7 Saving the computed new saturations at the FV centers
                if ( iter == this->max_nonlinear_limiting_case_iterations_ )
                    res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, true );  // O.K.
                else
                    res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, false ); // O.K.

                // measuring convergence
                //cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< res << endl;

                if ( res < this->target_nonlinear_limiting_case_residual_) break;
            }

        }

        // measuring convergence
        if (Verbose())
        {
            cout <<"\n\n\nIterateAdvectionEquation: Number of iterations: "<< iter << endl;
            cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< res << endl;
        }
        
        // 4.0 resetting the LHS, XVEC and RHS for the next iteration
        advector.ResetRHS( gref_.Nodes() );
        advector.ResetLHS( gref_.Nodes() );


    }else{

        size_t iter(1U);
        double res(numeric_limits<double>::max());

        //Initialization section
        double smin(0.0), smax(0.0);
        this->SAT0.resize( this->gref_.Nodes() );
        vector<double>(this->SAT0).swap(this->SAT0);
        this->InitialAdvectedPropertyValues(this->adv1_key_,smin,smax);
        //this->InitialAdvectedPropertyValues(this->ad1_key_);

        // First iteration:
        // 1. empty the sparse matrix and righthand vector
        advector.ResetLHS( this->gref_.Nodes() );
        advector.ResetRHS( this->gref_.Nodes() );

        MinMaxAdvectedProperty();

        this->grad_advprop_limiter_->CalculateGenericNodalGradient();
        this->grad_advprop_limiter_->CalculateSlopeLimiter( this->gref_, this->SMINMAX );

        for ( typename vector<Element<dim>*>::const_iterator eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
        {
            // 2.0 global to local ID conversion
            stencil_.eidx_ = (*eit)->Idx();

            // 2.1 during FIRST ITERATION: computing facet fluxes and interpolating saturations in FV sectors
            stencil_.InitializeSecondOrder( STENCIL_DATA[stencil_.eidx_], *(*eit) ); // O.K.

            // 2.2 FV volumes/dt and first-order fluxes are accumulated into lefthand side (matrix)
            advector.AccumulateLHS( stencil_, time_increment ); // O.K.

            // 2.3 source terms due a divergence of vt
            if ( with_flux_balance_correction ) 
            {
                stencil_.ComputeBoundaryFluxMismatch( STENCIL_DATA[stencil_.eidx_], *(*eit) );
                advector.AccumulateSectorSourceTermsInLHS( stencil_ );
            }

            // 2.4 conductance matrix for diffusion
            if ( with_diffusion ) advector.AccumulateIntegral_DNT_op_DN_dV_LHS( stencil_ );
            // alternative formulation for diffusion
            //if ( with_diffusion ) advector.AccumulateIntegral_DN_op_dS_LHS( stencil_ );

            // 2.5  update and compute values of advected variable at faces (spatial limiting)
            stencil_.ApplyLeastSquareMethodToLimitTransportProperties(*(*eit), SMINMAX, this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);

            // 2.6 account for the storage
            advector.AccumulateRHS( stencil_, SAT0, time_increment );

            // 2.7 add first-order Backward Euler solution to RHS
            //     the latest first-order solution is used here
            advector.AddToRHS( stencil_);

            // 2.8 previous solution multiplied by storage and divided by time increment is
            //     accumulated into righthandside for Backward-Euler time stepping
            advector.AccumulateHigherOrderRHS( stencil_ );

        }

        // 2.9 algebraic multigrid solver is applied to compute FV saturations
        advector.SolveMatrixEquation();

        // 2.10 Saving the computed new saturations at the FV centers
        if ( iter == this->max_nonlinear_limiting_case_iterations_ )
            res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, true );  // O.K.
        else
            res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, false ); // O.K.

        // measuring convergence
        //cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< res << endl;


        if ( res > this->target_nonlinear_limiting_case_residual_){

            //baseAdvector_->GetSolverSettings().Set_iswit(1);

            iter+=1;

            while( (res > this->target_nonlinear_limiting_case_residual_) && (iter < this->max_nonlinear_limiting_case_iterations_)){

                // 3. empty the righthand vector
                advector.ResetRHS( this->gref_.Nodes() );

                MinMaxAdvectedProperty();

                this->grad_advprop_limiter_->CalculateGenericNodalGradient();
                this->grad_advprop_limiter_->CalculateSlopeLimiter( this->gref_, this->SMINMAX );

                for ( typename vector<Element<dim>*>::const_iterator
                      eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ ){
                    // 3.0 global to local ID conversion
                    stencil_.eidx_ = (*eit)->Idx();

                    // 3.1 computing upstream and interpolating advected property values to FV segments
                    stencil_.InitializeAdvectedVariableValues( *(*eit) ); // O.K.

                    // 3.2  update and compute values of advected variable at faces (spatial limiting)
                    stencil_.ApplyLeastSquareMethodToLimitTransportProperties(*(*eit), SMINMAX, this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);

                    // 3.4 add first-order Backward Euler solution to RHS
                    //     the latest first-order solution is used here
                    advector.AddToRHS( stencil_);

                    // 3.5 previous solution multiplied by storage and divided by time increment is
                    //     accumulated into righthandside for Backward-Euler time stepping
                    advector.AccumulateHigherOrderRHS( stencil_ );

                } // end of accumulation


                // 3.6 algebraic multigrid solver is applied to compute FV saturations
                advector.SolveMatrixEquation();

                // 3.7 Saving the computed new saturations at the FV centers
                if ( iter == this->max_nonlinear_limiting_case_iterations_ )
                    res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, true );  // O.K.
                else
                    res = advector.OutputResultsWithL2NormRes( pref_, adv1_key_, false ); // O.K.

                // measuring convergence
                //cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< res << endl;

                if ( res < this->target_nonlinear_limiting_case_residual_) break;
            }

        }

        // measuring convergence
        if (Verbose())
        {
            cout <<"\n\n\nIterateAdvectionEquation: Number of iterations: "<< iter << endl;
            cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< res << endl;
        }
        // 4.0 resetting the LHS, XVEC and RHS for the next iteration
        advector.ResetRHS( gref_.Nodes() );
        advector.ResetLHS( gref_.Nodes() );
    }

} // end AdvectVariable2ndOrder (in space)



/**

See Matthai et al. 2009, TIPM for description of theta-limited
algorithm implemented here.

@todo (1) Does not create sharp profile yet and fails in the presence of sources and sinks */
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable2ndOrderInSpaceAndTime(
        NodeCenteredFiniteVolumeAlgorithm<dim>& advector,
        double time_increment,
        bool with_flux_balance_correction )
{
    static bool first_call(true);

    if (adv1_key_.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR,"NodeCenteredFiniteVolumeTransport<dim>::AdvectVariable2ndOrderInSpaceAndTime"," Second order advection does not support another variable type other than SCALAR.");
    // 1. store initial advected variable values in SAT0 and initialize
    //    the array LTDSATS0
    // ----------------------------------------------------------------
    InitialAdvectedPropertyValues( adv1_key_ ); // also gets min/max values // O.K.
    // FLUX_BALANCE is already updated by UpdateProjectedVelocitiesAndFluxBalances()

    // 2. Backward Euler time-stepping in the frame of non-linear iteration loop:
    //  ([vol]/dt + Vf){S}t+dt = {S}t {vol}/dt + Vf{S_upstr} + ...
    // --------------------------------------------------------------------------------------
    double             change_after_iteration;
    const size_t   iterations(100U);

    // diffusion is taken into account if the diffusion key is initialized
    const bool with_diffusion( (diff_key_ == csmp::Index()) ? false : true );


    if(!with_lsmgrad_limiter_){


        for ( size_t iteration=1; iteration<=iterations; iteration++ )
        {
            cout <<"\n\n\n\n\nNodeCenteredFiniteVolumeTransport<"<<  dim;
            cout <<">::AdvectVariable2ndOrderInSpaceAndTime: non-linear iteration "<< iteration <<" of "<< iterations;
            for ( typename vector<Element<dim>*>::const_iterator eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
            {
                // 2.0 global to local ID conversion
                stencil_.eidx_ = (*eit)->Idx();

                // 2.1 during FIRST ITERATION: computing facet fluxes and interpolating saturations in FV sectors
                if ( iteration == 1U ) {
                    stencil_.InitializeSecondOrder( STENCIL_DATA[stencil_.eidx_], *(*eit) ); // O.K.
                    // initializing t0 value arrays for higher-order solution procedure
                    if ( first_call ) FACETFLUXES0[stencil_.eidx_] = stencil_.facet_flux_;

                    // 2.2 FV volumes/dt and first-order fluxes are accumulated into lefthand side (matrix)
                    advector.AccumulateLHS( stencil_, time_increment ); // O.K.
                }
                else
                    // 2.3 computing upstream and interpolating advected property values to FV segments
                    stencil_.InitializeAdvectedVariableValues( *(*eit) ); // O.K.

                if ( with_diffusion ) advector.AccumulateIntegral_DNT_op_DN_dV_LHS( stencil_ );
                // alternative formulation for diffusion
                //if ( with_diffusion ) advector.AccumulateIntegral_DN_op_dS_LHS( stencil_ );

                // 2.4  update and compute values of advected variable at faces (spatial limiting)
                stencil_.IsotropicallyLimitTransportProperties( *(*eit), SMINMAX ); // O.K.

                // 2.5 compute theta limiter for each face
                if ( iteration < 50U )
                    stencil_.EvaluateThetaValues( *(*eit), FVPOREVOL, SAT0,
                                                  FACETFLUXES0, LTDSATS0, time_increment, false );
                else
                    stencil_.EvaluateThetaValues( *(*eit), FVPOREVOL, SAT0,
                                                  FACETFLUXES0, LTDSATS0, time_increment, true ); // true = max theta is chosen

                // 2.6 previous solution multiplied by storage and divided by time increment is
                //     accumulated into righthandside for Backward-Euler time stepping
                advector.AccumulateRHS( stencil_, SAT0, time_increment ); // must be SAT0 // O.K.

                // 2.7 higher-order saturations (and sources and sinks) are accumulated into righthand side
                advector.AccumulateHigherOrderRHS( stencil_, FACETFLUXES0, LTDSATS0 );

                // 2.8 add first-order Backward Euler solution to RHS
                //     the latest first-order solution is used here
                advector.AddToRHS( stencil_ );

            } // end of accumulation

            // 2.8 adjust for fluxbalances and correct LHS and RHS for out- & inflow contributions
            if ( iteration == 1U ) {
                // compensate for flow divergence related source or sink terms (all is done implicitly)
                if ( with_flux_balance_correction )
                    for ( size_t nidx=0U; nidx<gref_.Nodes(); nidx++ ) advector.AddToLHS( nidx, nidx, -FLUX_BALANCE[nidx] ); // O.K. for 1st-order
                // model and region boundaries (to left and righthand side)
                AssignFluxBoundaryConditions( advector, SecondOrderInTime() );
            }
            // model and region boundaries (to righthand side only)
            else AssignFluxBoundaryConditions( advector, SecondOrderInTime() );

            // 2.9 algebraic multigrid solver is applied to compute FV saturations
            advector.SolveMatrixEquation(); // O.K.

            // 2.10 resetting the LHS, XVEC and RHS for the next iteration
            advector.ResetRHS( gref_.Nodes() ); // O.K.

            // 2.11 Saving the computed new saturations at the FV centers
            if ( iteration == iterations ) change_after_iteration = advector.OutputResults( pref_, adv1_key_, true );  // O.K.
            else                           change_after_iteration = advector.OutputResults( pref_, adv1_key_, false ); // O.K.

            // measuring convergence
            cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< change_after_iteration << endl;
            //      	 if ( change_after_iteration < 1.0e-3 && iteration >= 3U ) break;
            if ( change_after_iteration < 1.0e-3 ) break;

            // 3.0 finding smin/smax of trial solution in the neighborhoos each element for next iteration
            // -------------------------------------------------------------------------------------------
            MinMaxAdvectedProperty(); // O.K.
        }

    }else{

        this->grad_advprop_limiter_->CalculateGenericNodalGradient();
        this->grad_advprop_limiter_->CalculateSlopeLimiter( this->gref_, this->SMINMAX );

        for ( size_t iteration=1; iteration<=iterations; iteration++ )
        {
            cout <<"\n\n\n\n\nNodeCenteredFiniteVolumeTransport<"<<  dim;
            cout <<">::AdvectVariable2ndOrderInSpaceAndTime: non-linear iteration "<< iteration <<" of "<< iterations;
            for ( typename vector<Element<dim>*>::const_iterator eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
            {
                // 2.0 global to local ID conversion
                stencil_.eidx_ = (*eit)->Idx();

                // 2.1 during FIRST ITERATION: computing facet fluxes and interpolating saturations in FV sectors
                if ( iteration == 1U ) {
                    stencil_.InitializeSecondOrder( STENCIL_DATA[stencil_.eidx_], *(*eit) ); // O.K.
                    // initializing t0 value arrays for higher-order solution procedure
                    if ( first_call ) FACETFLUXES0[stencil_.eidx_] = stencil_.facet_flux_;

                    // 2.2 FV volumes/dt and first-order fluxes are accumulated into lefthand side (matrix)
                    advector.AccumulateLHS( stencil_, time_increment ); // O.K.
                }
                else
                    // 2.3 computing upstream and interpolating advected property values to FV segments
                    stencil_.InitializeAdvectedVariableValues( *(*eit) ); // O.K.

                if ( with_diffusion ) advector.AccumulateIntegral_DNT_op_DN_dV_LHS( stencil_ );
                // alternative formulation for diffusion
                //if ( with_diffusion ) advector.AccumulateIntegral_DN_op_dS_LHS( stencil_ );

                // 2.4  update and compute values of advected variable at faces (spatial limiting)
                stencil_.ApplyLeastSquareMethodToLimitTransportProperties(*(*eit), SMINMAX, this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);

                // 2.5 compute theta limiter for each face
                if ( iteration < 50U )
                    stencil_.EvaluateThetaValues( *(*eit), FVPOREVOL, SAT0,
                                                  FACETFLUXES0, LTDSATS0, time_increment, false );
                else
                    stencil_.EvaluateThetaValues( *(*eit), FVPOREVOL, SAT0,
                                                  FACETFLUXES0, LTDSATS0, time_increment, true ); // true = max theta is chosen

                // 2.6 previous solution multiplied by storage and divided by time increment is
                //     accumulated into righthandside for Backward-Euler time stepping
                advector.AccumulateRHS( stencil_, SAT0, time_increment ); // must be SAT0 // O.K.

                // 2.7 higher-order saturations (and sources and sinks) are accumulated into righthand side
                advector.AccumulateHigherOrderRHS( stencil_, FACETFLUXES0, LTDSATS0 );

            } // end of accumulation

            // 2.8 adjust for fluxbalances and correct LHS and RHS for out- & inflow contributions
            if ( iteration == 1U ) {
                // compensate for flow divergence related source or sink terms (all is done implicitly)
                if ( with_flux_balance_correction )
                    for ( size_t nidx=0U; nidx<gref_.Nodes(); nidx++ ) advector.AddToLHS( nidx, nidx, -FLUX_BALANCE[nidx] ); // O.K. for 1st-order
                // model and region boundaries (to left and righthand side)
                AssignFluxBoundaryConditions( advector, SecondOrderInTime() );
            }
            // model and region boundaries (to righthand side only)
            else AssignFluxBoundaryConditions( advector, SecondOrderInTime() );

            // 2.9 algebraic multigrid solver is applied to compute FV saturations
            advector.SolveMatrixEquation(); // O.K.

            // 2.10 resetting the LHS, XVEC and RHS for the next iteration
            advector.ResetRHS( gref_.Nodes() ); // O.K.

            // 2.11 Saving the computed new saturations at the FV centers
            if ( iteration == iterations ) change_after_iteration = advector.OutputResults( pref_, adv1_key_, true );  // O.K.
            else                           change_after_iteration = advector.OutputResults( pref_, adv1_key_, false ); // O.K.

            // measuring convergence
            cout <<"\n\n\nIterateAdvectionEquation: Change since last iteration: "<< change_after_iteration << endl;
            //      	 if ( change_after_iteration < 1.0e-3 && iteration >= 3U ) break;
            if ( change_after_iteration < 1.0e-3 ) break;

            // 3.0 finding smin/smax of trial solution in the neighborhoos each element for next iteration
            // -------------------------------------------------------------------------------------------
            MinMaxAdvectedProperty(); // O.K.
        }

    }

    advector.ResetLHS( gref_.Nodes() );

    first_call = false;

} // end AdvectVariable2ndOrderInSpaceAndTime



/**

Reads the scalar node variable in the region identified as group and
multiplies it by the finite volume or area surrounding each node. The
result is returned into the variable.

@section implementation Implementation

Uses a map to store the temporary variable values.

@section application Application

For instance, to use node variables to assign Neumann boundary
conditions to a finite element model.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::MultiplyScalarNodePropertyByFiniteVolume( const char* property )
{
    csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type != SCALAR || prop_key.place != NODE ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::MultiplyScalarNodePropertyByFiniteVolume",
                               "This method only handles scalar node properties / node-centered finite volume variables" );
        return;
      }

    vector<double>  node_data( gref_.Nodes(), 0. );

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
        for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
            // read variable
            double prop_val = (*eit)->N(i)->Read( prop_key );
            // multiply with sector volume
            prop_val *= ( *(*eit) ).SectorVolume(i);
            node_data[ (*eit)->N(i)->Idx() ] += prop_val;
        }

    // write variable
    ScalarVariable  sc;

    for ( auto i{0U}; i<node_data.size(); i++ ) {
        sc.Flag( ) = gref_.N(i)->Status( prop_key );
        sc        = node_data[i];
        gref_.N(i)->Store( prop_key, sc );
    }

    cout <<"\nNodeCenteredFiniteVolumeTransport<"<< dim <<">::MultiplyScalarNodePropertyByFiniteVolume: ";
    cout <<"multiplied '"<< property <<"' by node-centered finite volume / area."<< endl;

} // end MultiplyScalarNodePropertyByFiniteVolume





/**  VolumeIntegrateScalarFiniteVolumeVariable

To integrate node variables over the entire simulation model using the finite volume framework.

@section arguments Input Arguments

The current model (Region), and the name of the variable of interest
as specified in the property database.

@return The method returns the volume integral of the variable in the model.

@section application Application

This method will not integrate element variables. The obvious choice
for their volume integration is the finite element method.
*/
template<uint32_t dim>
double  NodeCenteredFiniteVolumeTransport<dim>::VolumeIntegrateScalarFiniteVolumeVariable( const char* property,
                                                                                             bool take_porosity_into_account ) const
{
    csmp::Index  prop_key = pref_.StorageKey(property);
    double     result(0.);

    if ( prop_key.type != SCALAR || prop_key.place == ELEMENT ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrateScalarFiniteVolumeVariable",
                               "This method only handles scalar node properties / node-centered finite volume variables" );
        return result;
    }

    if ( take_porosity_into_account ) {
        csmp::Index  phi_key = pref_.StorageKey("porosity");

        for ( typename vector<Element<dim>*>::const_iterator
              eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
        {
            assert( (*eit)->FV() != NULL );
            double phi = (*eit)->Read( phi_key );
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
                result += (*eit)->N(i)->Read( prop_key ) * phi * ( *(*eit) ).SectorVolume(i);
        }
    }
    else {
        for ( typename vector<Element<dim>*>::const_iterator
              eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ ) {
            assert( (*eit)->FV() != NULL );
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
                result += (*eit)->N(i)->Read( prop_key ) * ( *(*eit) ).SectorVolume(i);
        }
    }

    return result;

} // end VolumeIntegrateScalarFiniteVolumeVariable






/**  VolumeIntegrateScalarFiniteVolumeVariable

This method integrates the variable of interest over the region , i.e.
sub-region of the current model as specified by the user.
Upon request, the result variables are normalized by the variable 'porosity'

@section arguments Input Arguments

A reference to the current model, the name of the subregion and the
nodal property of interest, and a boolean variable that determines whether
the porosity of the elements shall be taken into account by the
integration.

@return The method returns the volume integral of the property of interest in
the model subregion.

@section application Application

For monitoring of variables in model subregions.
*/
template<uint32_t dim>
double  NodeCenteredFiniteVolumeTransport<dim>::VolumeIntegrateScalarFiniteVolumeVariable( const char* group,
                                                                                             const Model<dim>& sg,
                                                                                             const char* property,
                                                                                             bool take_porosity_into_account ) const
{
    assert( gref_.Cells() == sg.Region("Model").Cells() );

    csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type != SCALAR || prop_key.place != NODE )
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrateScalarFiniteVolumeVariable",
                               "This method only handles scalar node properties, i.e. node-centered finite volume variables" );

    double  interim_result, result(0.);
    const Region<dim>&  gref(sg.Region(group));

    if ( take_porosity_into_account ) {
        csmp::Index  phi_key = pref_.StorageKey("porosity");
        assert( phi_key.type  == SCALAR );
        assert( phi_key.place == ELEMENT );

        for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ )
        {
            assert( (*eit)->FV() != NULL );
            interim_result = 0.;
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
                interim_result += (*eit)->N(i)->Read( prop_key ) * ( *(*eit) ).SectorVolume(i);
            result += interim_result * (*eit)->Read( phi_key );
        }
    }
    else {
        for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ ) {
            assert( (*eit)->FV() != NULL );
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
                result += (*eit)->N(i)->Read( prop_key ) * ( *(*eit) ).SectorVolume(i);
        }
    }

    return result;

} // end VolumeIntegrateScalarFiniteVolumeVariable (group)







/** This method integrates the nodal finite-element variable after
interpolation to the finite-volume sector integration points.
23-09-2013 - Julian M added the possibility to integrate over a particular region,
which is part of the model.

@section arguments Input Arguments

A reference to the current model and a char string pointer to the
finite element node variable of interest.

@return The method returns the (finite) volume integral of the finite element
variable over the current model.

@section application Application

There has to be a good reason for the application of this more costly
method as opposed to the direct integration of nodal variables.
*/
template<uint32_t dim>
double NodeCenteredFiniteVolumeTransport<dim>::VolumeIntegrateScalarFiniteElementVariable(
        const char* property,
        bool take_porosity_into_account,const char* region ) const
{
    csmp::Index  prop_key = pref_.StorageKey(property);
    string region_name;
    if (region==NULL)
        region_name=gref_.Name();
    else{
        region_name=string(region);
        if (!mref_.ContainsRegion(region_name.c_str()))
            throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrateScalarFiniteElementVariable",
                                   "Region requested for integral is not contained in the model." );
    }

    Region<dim>& rref=mref_.Region(region_name.c_str());

    if ( prop_key.type != SCALAR )
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrateScalarFiniteElementVariable",
                               "This method only handles scalar properties" );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    double result(0.);

    if ( take_porosity_into_account ) {
        csmp::Index  phi_key = pref_.StorageKey("porosity");
        assert( phi_key.type  == SCALAR );
        assert( phi_key.place == ELEMENT );

        if ( phi_key != prop_key )
            csmp_error.notice( WARNING, "NodeCenteredFiniteVolumeTransport<dim>::VolumeIntegrateScalarFiniteElementVariable",
                               "you are trying to scale the integral over the porosity with porosity");

        double interim_result;

        for ( typename vector<Element<dim>*>::const_iterator
              eit=rref.CellsBegin(); eit!=rref.CellsEnd(); eit++ )
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
                // interpolate property to finite volume sector integration points
                interim_result  = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, prop_key );
                interim_result *= ( *(*eit) ).SectorVolume(i) * (*eit)->Read( phi_key );
                result         += interim_result;
            }
    }
    else {
        for ( typename vector<Element<dim>*>::const_iterator
              eit=rref.CellsBegin(); eit!=rref.CellsEnd(); eit++ )
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
                // interpolate property to finite volume sector integration points
                double interim_result  = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, prop_key );
                interim_result *= ( *(*eit) ).SectorVolume(i);
                result         += interim_result;
            }
    }

    return result;

} // end VolumeIntegrateScalarFiniteElementVariable




/** The volumes of the finite volumes is reported into the result variable.

Also returns the total volume or porevolume of the model as calculated
from the finite volume discretization.
*/
template<uint32_t dim>
double  NodeCenteredFiniteVolumeTransport<dim>::FiniteVolume( const char* volume_property ) const
{
    csmp::Index  prop_key = pref_.StorageKey(volume_property);

    if ( prop_key.type != SCALAR || prop_key.place != NODE )
        throw csmp::Exception( ERROR, "NodeCenteredFiniteVolumeTransport<dim>::FiniteVolume",
                               volume_property, "to which the finite volume is assigned, must be a scalar variable placed on the nodes" );

    double          result(0.);
    vector<double>  volumes( gref_.Nodes(), 0. );

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
        for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
            volumes[ (*eit)->N(i)->Idx() ] += ( *(*eit) ).SectorVolume(i);

    auto nit=gref_.NodesBegin();
    for ( typename vector<double>::const_iterator it=volumes.begin(); it!=volumes.end(); it++, nit++ ) {
        (*nit)->Store( prop_key, makeScalar( (*nit)->Status( prop_key ), (*it)) );
        result += (*it);
    }

    return result;

} // end FiniteVolumes



/**

Computes integrals for finite-element computations on the entire Region.
The first parameter specifies the property to be integrated and the
second one the property into which the result will be stored.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::VolumeIntegrate( const char* integrand_property,
                                                              const char* result_property ) const
{
    csmp::Index  iprop_key = pref_.StorageKey(integrand_property);
    csmp::Index  rprop_key = pref_.StorageKey(result_property);

    if ( iprop_key.type != SCALAR ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "This method only integrates scalar properties" );
    }
    if ( rprop_key.type != SCALAR || rprop_key.place != NODE ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "Result property must be scalar placed on the nodes / node-centered finite volumes" );
    }

    vector<double>  nresult( gref_.Nodes(), 0. );

    if ( iprop_key.place == NODE ) {
        for ( typename vector<Element<dim>*>::const_iterator
              eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
                nresult[ (*eit)->N(i)->Idx() ] += (*eit)->N(i)->Read( iprop_key ) * STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
    }
    else if ( iprop_key.place == ELEMENT ) {
        for ( typename vector<Element<dim>*>::const_iterator
              eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
            for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
                size_t nidx = (*eit)->N(i)->Idx();
                // interpolate property to finite volume sector integration points
                double res = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, iprop_key );

                nresult[nidx] += res * STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
            }
    }
    else {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "Method does not handle Face or IntegrationPoint variables" );
        return;
    }

    for ( auto i{0U}; i<gref_.Nodes(); i++ )
        gref_.N(i)->Store( rprop_key, ScalarVariable( gref_.N(i)->Status( rprop_key ), nresult[i]) );

} // end VolumeIntegrate (Region, single property)



template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::VolumeIntegrate( const char* integrand_property,
                                                              const char* integrand_multiplier,
                                                              const char* result_property ) const
{
    const csmp::Index  iprop_key = pref_.StorageKey(integrand_property);
    const csmp::Index  mprop_key = pref_.StorageKey(integrand_multiplier);
    const csmp::Index  rprop_key = pref_.StorageKey(result_property);

    if ( iprop_key.type != SCALAR ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "This method only integrates scalar properties" );
    }
    if ( rprop_key.type != SCALAR || rprop_key.place != NODE ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "Result property must be scalar placed on the nodes / node-centered finite volumes" );
    }

    vector<double>  nresult( gref_.Nodes(), 0. );

    if ( iprop_key.place == NODE ) {
        if ( mprop_key.place == NODE ) {
            for ( typename vector<Element<dim>*>::const_iterator
                  eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
                for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
                    nresult[ (*eit)->N(i)->Idx() ] +=
                            (*eit)->N(i)->Read( iprop_key ) *
                            (*eit)->N(i)->Read( mprop_key ) *
                            STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
        }
        else if ( mprop_key.place == ELEMENT ) {
            for ( typename vector<Element<dim>*>::const_iterator
                  eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
                for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
                    double int_mult = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, mprop_key );
                    nresult[ (*eit)->N(i)->Idx() ] +=
                            (*eit)->N(i)->Read( iprop_key ) * int_mult *
                            STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
                }
        }
        else csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                              "Method does not handle Face or IntegrationPoint integrand multipliers" );
    }
    else if ( iprop_key.place == ELEMENT ) {
        if ( mprop_key.place == ELEMENT ) {
            for ( typename vector<Element<dim>*>::const_iterator
                  eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
                for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
                    // interpolate element properties to finite volume sector integration points
                    double iprop = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, iprop_key );
                    double mprop = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, mprop_key );

                    nresult[ (*eit)->N(i)->Idx() ] += iprop * mprop * STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
                }
        }
        else if ( mprop_key.place == NODE ) {
            for ( typename vector<Element<dim>*>::const_iterator
                  eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ )
                for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
                    // interpolate element properties to finite volume sector integration points
                    double iprop = ( *(*eit) ).PropertyValueAtSectorIntegrationPoint( i, 0U, iprop_key );

                    double mprop = (*eit)->N(i)->Read( mprop_key );

                    nresult[ (*eit)->N(i)->Idx() ] += iprop * mprop * STENCIL_DATA[ (*eit)->Idx() ].SectorVolume(i);
                }
        }
        else csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                              "Method does not handle Face or IntegrationPoint integrand multipliers" );
    }
    else {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "Method does not handle Face or IntegrationPoint variables" );
        return;
    }

    for ( auto i{0U}; i<gref_.Nodes(); i++ )
        gref_.N(i)->Store( rprop_key, makeScalar( gref_.N(i)->Status( rprop_key ), nresult[i] ) );

} // end VolumeIntegrate( entire model, with integral multiplier )




/**

Element vector<double> property is projected on and integrated over finite volume
facets; the measured divergence is stored in the variable 'result_property'.
 */
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::Divergence( const char* div_property,
                                                         const char* result_property ) const
{
    const csmp::Index  dprop_key = pref_.StorageKey(div_property);
    const csmp::Index  rprop_key = pref_.StorageKey(result_property);

    if ( dprop_key.type != VECTOR  or  dprop_key.place != ELEMENT  ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               div_property, "must be a vector<double> variable placed on the element" );
    }
    if ( rprop_key.type != SCALAR || rprop_key.place != NODE ) {
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::VolumeIntegrate",
                               "Result property must be scalar placed on the node / node-centered finite volume" );
    }

    // perform surface integration of projected velocities
    vector<double>       nresult( gref_.Nodes(), 0. );
    uint32_t             inside_node, outside_node;
    VectorVariable<dim>  velo;

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.CellsBegin(); eit!=gref_.CellsEnd(); eit++ ) {
        (*eit)->Read( dprop_key, velo );
        for ( auto i{0U}; i<(*eit)->FV()->Facets(); i++ )
        {
            // projecting velocity onto facet normal i
            (*eit)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
            double proj = STENCIL_DATA[ (*eit)->Idx() ].FacetNormalProjection( i, velo );
            proj   *= STENCIL_DATA[ (*eit)->Idx() ].FacetArea( i );

            nresult[ (*eit)->N(inside_node)->Idx() ]  += proj;
            nresult[ (*eit)->N(outside_node)->Idx() ] -= proj;
        }
    }

    for ( auto i{0U}; i<nresult.size(); i++ )
        gref_.N(i)->Store( rprop_key,
                           ScalarVariable( gref_.N(i)->Status( rprop_key ), nresult[i]) );

} // end Divergence (Region)




/**
Computes gradients normal to the outer model boundary from influxes (sources)
assigned and flagged as NEUMANN elsewhere. To do this it needs material parameters
from the corresponding flow law, for instance the hydraulic conductivity
in the case of Darcy's law.

The convention used is that incoming fluxes have to be positive and outgoing ones negative.

@attention The method only works for box-shaped models and will only operate
on its boundaries. Nodal Neumann flagged terms elsewhere will remain untouched.

@param boundary is the boundary flag of the box-shaped model, IRREGULAR
is valid as well.

@param node_property stores the boundary flux which shall be applied.

@param propertionality_constant is the element variable the method needs in
order to compute the boundary normal gradient which should ensue from the applied
flux term, for instance 'hydraulic conductivity' in the case of Darcy's law.

@param distribute_total_amount tells the method how to interpret the boundary
flux term it reads. If true, the term will be treated as integrated total for
the boundary, else it will be assigned as a flux per unit area.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::TransformScalarBoundaryValuesIntoNeumannConditions(
        BOX_BOUNDARY boundary,
        const char* node_property,
        const char* propertionality_constant,
        bool distribute_total_amount )
{
    const csmp::Index  prop_key = pref_.StorageKey(node_property);
    const csmp::Index  cond_key = pref_.StorageKey(propertionality_constant);

    // here a check whether the model is box-shaped needs to be performed

    if ( prop_key.type != SCALAR || prop_key.place != NODE )
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::TransformScalarBoundaryValuesIntoNeumannConditions",
                               "This method only handles scalar node properties or node-centered finite volume variables" );

    if ( cond_key.type != SCALAR || cond_key.place != ELEMENT )
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::TransformScalarBoundaryValuesIntoNeumannConditions",
                               "The proportionality constant used in the Neumann gradient calculation must be a scalar element property" );

    // setting the normal up for the projection taking into account the model boundary
    vector<double>     bnormal(3U);
    Box                  box;  box.UnitNormalTo( boundary, 3U, bnormal );
    VectorVariable<dim>  vc;   for ( auto i{0U}; i<bnormal.size(); i++ ) vc(i) = -bnormal[i];

    ScalarVariable  sc( NEUMANN, 0.);
    double          total_surface_area(0.), projection;
    uint32_t        inside_node, outside_node;

    // 1. if the specified quantity shall be distributed over the entire surface area of the boundary,
    //    this area must be computed first
    // -----------------------------------
    if ( distribute_total_amount and dim != 1U ) {
        // add the total volumes up
        for ( auto n=gref_.InteriorNodes(); n<gref_.Nodes(); n++ )
            if ( gref_.N(n)->AtBoundary() == boundary )
                // loop over the parent elements of the node
                for ( auto t=0U; t<gref_.N(n)->Parents(); t++ ) {
                    auto nid = gref_.N(n)->ParentNodeNumber(t);
                    // integrate surface area of finite volume facets as projected onto unit normal to model boundary
                    for ( auto i{0U}; i<gref_.N(n)->Parent(t)->FV()->Facets(); i++ )
                    {
                        gref_.N(n)->Parent(t)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                        // project only if necessary
                        if ( nid == inside_node || nid == outside_node ) {
                            projection = ( *gref_.N(n)->Parent(t) ).ProjectionOnFacetNormal( i, vc );
                            // the area is added to the total surface area of the model boundary
                            // only those facets that delimited sector of target node are considered
                            if      ( nid == inside_node )
                                total_surface_area += ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                            else if ( nid == outside_node )
                                total_surface_area -= ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                        }
                    }
                }
        cout <<"\nFNodeCenteredFiniteVolumeTransport<"<< dim;
        cout <<">::TransformScalarBoundaryValuesIntoNeumannConditions: ";
        cout <<"Surface area of target boundary is: "<< total_surface_area;
    }
    else total_surface_area = 1.;


    // 2. assigment of finite-volume cross-sectional area dependent nodal property terms
    // ---------------------------------------------------------------------------------
    for ( auto n=gref_.InteriorNodes(); n<gref_.Nodes(); n++ )
        if ( gref_.N(n)->AtBoundary() == boundary )
        {
            // reading the nodal source term counting it if it is flagged Neumann
            double  val = (gref_.N(n)->Status( prop_key ) == NEUMANN ) ? gref_.N(n)->Read( prop_key ) : 0.;

            // computing the harmonic mean of the material parameters of the elements surrounding the node
            double  material_param(0.);
            for ( auto t=0U; t<gref_.N(n)->Parents(); t++ )
                material_param += 1. / gref_.N(n)->Parent(t)->Read( cond_key );
            material_param = static_cast<double>(gref_.N(n)->Parents()) / material_param;

            // computing the Neumann gradient terms
            sc = (dim == 1U) ? 1. : 0.;
            if ( dim != 1U )
                for ( auto t=0U; t<gref_.N(n)->Parents(); t++ ) {
                    auto nid = gref_.N(n)->ParentNodeNumber(t);
                    for ( auto i{0U}; i<gref_.N(n)->Parent(t)->FV()->Facets(); i++ )
                    {
                        gref_.N(n)->Parent(t)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                        // project only if necessary
                        if ( nid == inside_node || nid == outside_node ) {
                            projection = ( *gref_.N(n)->Parent(t) ).ProjectionOnFacetNormal( i, vc );
                            // the area is added to the total surface area of the model boundary
                            // only those facets that delimited sector of target node are considered
                            if      ( nid == inside_node )
                                sc += ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                            else if ( nid == outside_node )
                                sc -= ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                        }
                    }
                }

            // integrating the influx term over the FVM cross-sectional area
            sc *= val / total_surface_area;
            // dividing it by the harmonic mean of the material property to obtain the target gradient
            sc /= material_param;
            // writing the integrated value to the target variable for later assignment by PointSource_rhsop.
            gref_.N(n)->Store( prop_key, sc );
        }

} // end TransformScalarBoundaryValuesIntoNeumannConditions




/**

For box-shaped models, this method allows to assign boundary normal
surface integrated variables using the finite volume discretization.
The method works only for box-shaped models because otherwise it
would need to have calculations of the surface normal which
it can't do without access to face variables.

This could for instance be a surface normal flux that shall be integrated
over the surface area of the finite elements.

@section application Application

Use method to integrate variables like "nodal volume source". Resulting
quantities can be accumulated directly (without integration) into the
finite element equations. For this purpose, PDE operators like the
PointSource_rhsop can be used.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::AssignScalarBoundaryValues(
        BOX_BOUNDARY boundary,
        const char* property,
        VARIABLE_FLAG bcond, double val,
        bool distribute_total_amount )
{
    const csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type != SCALAR || prop_key.place != ELEMENT )
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::AssignScalarBoundaryValues",
                               "This method only handles scalar node properties / node-centered finite volume variables" );

    // setting the normal up for the projection taking into account the model boundary
    vector<double>     bnormal(3U);
    Box                  box;      box.UnitNormalTo( boundary, 3U, bnormal );
    VectorVariable<dim>  vc;       for ( auto i{0U}; i<bnormal.size(); i++ ) vc(i) = -bnormal[i];

    ScalarVariable  sc( bcond, 0.);
    double          total_surface_area(0.), projection;
    uint32_t        inside_node, outside_node;

    // 1. if the quantity shall be distributed over the surface area of the boundary,
    //    this area must be computed first
    // -----------------------------------
    if ( distribute_total_amount and dim != 1U ) {
        // add the total volumes up
        for ( auto n=gref_.InteriorNodes(); n<gref_.Nodes(); n++ )
            if ( gref_.N(n)->AtBoundary() == boundary )
                // loop over the parent elements of the node
                for ( auto t=0U; t<gref_.N(n)->Parents(); t++ ) {
                    auto nid = gref_.N(n)->ParentNodeNumber(t);
                    // integrate surface area of finite volume facets as projected onto unit normal to model boundary
                    for ( auto i{0U}; i<gref_.N(n)->Parent(t)->FV()->Facets(); i++ )
                      {
                          gref_.N(n)->Parent(t)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                          // project only if necessary
                          if ( nid == inside_node || nid == outside_node ) {
                              projection = ( *gref_.N(n)->Parent(t) ).ProjectionOnFacetNormal( i, vc );
                              // the area is added to the total surface area of the model boundary
                              // only those facets that delimited sector of target node are considered
                              if      ( nid == inside_node )
                                  total_surface_area += ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                              else if ( nid == outside_node )
                                  total_surface_area -= ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                          }
                      }
                  }
        cout <<"\nNodeCenteredFiniteVolumeTransport<"<< dim;
        cout <<">::AssignScalarBoundaryValues: Surface area of target boundary is: "<< total_surface_area;
    }
    else total_surface_area = 1.;

    // 2. assigment of finite-volume cross-sectional area dependent nodal property terms
    // ---------------------------------------------------------------------------------
    for ( auto n=gref_.InteriorNodes(); n<gref_.Nodes(); n++ )
        if ( gref_.N(n)->AtBoundary() == boundary )
        {
            sc = (dim == 1U) ? 1. : 0.;
            if ( dim != 1U )
                for ( auto t=0U; t<gref_.N(n)->Parents(); t++ ) {
                    auto nid = gref_.N(n)->ParentNodeNumber(t);
                    for ( auto i{0U}; i<gref_.N(n)->Parent(t)->FV()->Facets(); i++ )
                    {
                        gref_.N(n)->Parent(t)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                        // project only if necessary
                        if ( nid == inside_node || nid == outside_node ) {
                            projection = ( *gref_.N(n)->Parent(t) ).ProjectionOnFacetNormal( i, vc );
                            // the area is added to the total surface area of the model boundary
                            // only those facets that delimited sector of target node are considered
                            if      ( nid == inside_node )
                                sc += ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                            else if ( nid == outside_node )
                                sc -= ( *gref_.N(n)->Parent(t) ).FacetArea( i ) * projection;
                        }
                    }
                }

            sc *= val / total_surface_area;
            // write the integrated value to the target variable
            gref_.N(n)->Store( prop_key, sc );
        }

} // end AssignScalarBoundaryValues




/** As above but reads values from the nodes before converting them into
equivalent sources.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::MultiplyScalarBoundaryValuesByFiniteVolumeCrossSectionalArea(
        Model<dim>& model,BOX_BOUNDARY boundary, const char* property )
{
    csmp::Index  prop_key = pref_.StorageKey(property);
    Region<dim>  rref(model.Region(parseBoundary(boundary).c_str()));

    if ( prop_key.type != SCALAR || prop_key.place != NODE )
        throw csmp::Exception( ERROR, "FiniteVolumeTransport::ConvertScalarBoundaryValues",
                               "This method only handles scalar node properties / node-centered finite volume variables" );

    // setting the normal up for the projection taking into account the model boundary
    vector<double>       bnormal(3U);
    Box                  box;  box.UnitNormalTo( boundary, 3U, bnormal );
    VectorVariable<dim>  vc;   for ( auto i{0U}; i<bnormal.size(); i++ ) vc(i) = -bnormal[i];
    ScalarVariable       sc;
    double               surface_area, projection;
    uint32_t             inside_node, outside_node;

    // 1. assigment of finite-volume cross-sectional area dependent nodal property terms
    // ---------------------------------------------------------------------------------
    cout <<"\nINFO: MultiplyScalarBoundaryValuesByFiniteVolumeCrossSectionalArea::ConvertScalarBoundaryValues";
    cout <<" Converting the original nodal values to ones weighted by the surface area. ";
    cout <<"Concerned variable '"<< property <<"'"<< endl;

    for ( auto n=0U; n<rref.Nodes(); n++ )
    {
        surface_area = 0.;
        for ( auto t=0U; t<rref.N(n)->Parents(); t++ ) {
            auto nid = rref.N(n)->ParentNodeNumber(t);
            for ( auto i{0U}; i<rref.N(n)->Parent(t)->FV()->Facets(); i++ )
            {
                rref.N(n)->Parent(t)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                // project only if necessary
                if ( nid == inside_node || nid == outside_node ) {
                    projection = ( *rref.N(n)->Parent(t) ).ProjectionOnFacetNormal( i, vc );
                    // the area is added to the total surface area of the model boundary
                    // only those facets that delimited sector of target node are considered
                    if ( nid == inside_node )
                        surface_area += ( *rref.N(n)->Parent(t) ).FacetArea( i ) * projection;
                    else if ( nid == outside_node )
                        surface_area -= ( *rref.N(n)->Parent(t) ).FacetArea( i ) * projection;
                }
            }
        }

        // converting the assigned value with the  surface area of the FV intersection surrounding the node
        rref.N(n)->Read( prop_key, sc );
        sc *= surface_area;
        // write the integrated value to the target variable
        rref.N(n)->Store( prop_key, sc );
    }

} // end ConvertScalarBoundaryValues



/**

Writes internal (state) variables of the NCFVT to the screen. Note that
these may vary dependent on the choice of constructor which was used
when the object was built.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeTransport<dim>::Out() const
{
    cout <<"\nNodeCenteredFiniteVolumeTransport<"<<  dim <<">::Out: ";
    cout <<"Variable to advect '"<< advected_variable_ <<"'"<< endl;

    cout <<"\nFV_STENCIL_DATA - sectors, facets, facet-areas and normals for each finite element: "<< endl;
    for ( size_t i{0U}; i<STENCIL_DATA.size(); i++ ) {
        cout <<"\nStencil "<< i+1;
        STENCIL_DATA[i].Out();
    }
    cout << endl;

    cout <<"\nFVPOREVOL, FLUX_BALANCE pore volumes and flux balances for each finite volume: "<< endl;
    for ( size_t i{0U}; i<FVPOREVOL.size(); i++ )
        cout << i+1 <<": "<< FVPOREVOL[i] <<" "<< FLUX_BALANCE[i] << endl;

    if ( SecondOrderInSpace() ) {
        cout <<"\nThe transport algorithm was constructed for a higher order method; It also contains the following arrays."<< endl;

        cout <<"\nSMINMAX - min,max of transport variable in the neighborhood of finite volume:"<< endl;
        for ( size_t i{0U}; i<SMINMAX.size(); i++ )
            cout << i+1 <<": "<< SMINMAX[i].first <<" "<< SMINMAX[i].second << endl;
        cout << endl;
    }

    if ( SecondOrderInTime() ) {
        cout <<"\nSAT0 - initial transport variable values for each finite volume: "<< endl;
        for ( size_t i{0U}; i<FVPOREVOL.size(); i++ )
            cout << i+1 <<": "<< SAT0[i] <<" ";
        cout << endl;

        cout <<"\nFACETFLUXES0 - fluxes across finite volume facets in each finite element at time levels 0:"<< endl;
        for ( size_t i{0U}; i<FACETFLUXES0.size(); i++ ) {
            cout <<" element "<< i+1;
            for ( size_t j{0U}; j<FACETFLUXES0[i].size(); j++ )
                cout <<"\n\tfacet "<< j <<": "<< FACETFLUXES0[i][j] <<" ";
            cout << endl;
        }
        cout << endl;

        cout <<"\nLTDSATS0 - limited facet values of transported variable at time levels 0:"<< endl;
        for ( size_t i{0U}; i<LTDSATS0.size(); i++ ) {
            cout <<" element "<< i+1;
            for ( size_t j{0U}; j<LTDSATS0[i].size(); j++ )
                cout <<"\n\tfacet "<< j <<": "<< LTDSATS0[i][j] <<" ";
            cout << endl;
        }
        cout << endl;
    }

} // end Out



/**

Method performs a number of tests on the generic finite volume
scheme implementation that underpins the SinglePhaseExplicitNodeCenteredFVTransport
algorithm. The results of the tests are reported on the screen.

@section arguments Input Arguments

The method needs access to the current model and the user should specify
the precision with which the results shall be output. While a smaller
precision makes the results more readable, remember that the plotted
numbers will have been rounded up or down accordingly.

Function needs the nodal variable fluid pressure to be defined.

@section messages Messages

The method reports properties like FV cell and FV sector volumes,
projected velocities and FV facet areas.
*/
template<uint32_t dim>
bool testFiniteVolumeStencil( const PropertyDatabase<dim>& p, const Region<dim>& gref,
                              NodeCenteredFiniteVolumeTransport<dim>& ncvft, long data_precision )
{
    // renumbering nodes and elements
    gref.UpdateMemberIndexes();
    //                           in- and outflow for each FV
    vector<pair<double,double> >  THROUGHPUT( gref.Nodes(), make_pair(0.,0.) );
    VectorVariable<dim>               velo;
    double                            evolume, flux, area, sum, val;
    const double                      zero(0.);
    uint32_t                          inside_node, outside_node;
    vector<ScalarVariable >           node_prop;
    csmp::Index                       prop_key = p.StorageKey("fluid pressure");

    long prec = cout.precision(data_precision);

    // 1. testing the flux integration over the FV surfaces with a prescribed velocity
    // -------------------------------------------------------------------------------------
    cout <<"\ntestFiniteVolumeStencil(NodeCenteredFiniteVolumeTransport): Test 1, velocity integration in X-direction, velo: "<< endl;
    velo    = 0.;
    velo(0) = 1.;
    velo.Out();
    vector<double>  rst(3), IPOL;

    for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ )
    {
        cout <<"\nElement: "<< (*eit)->Idx() <<", vol: "<< (evolume=(*eit)->Volume()) <<", type: ";
        cout << parseFiniteElementType( (*eit)->FE_Type() );
        cout <<"\nXYZ of nodes: "<< endl <<"\t";
        for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
            cout <<"("<< i+1 <<") "<< (*eit)->N(i)->Coordinate() <<" ";
        cout << endl;

        // 2. testing the surface areas and velocity projections
        // -----------------------------------------------------
        cout <<"\nFacet fluxes, areas and inside/outside nodes: "<< endl <<"\t";
        for ( auto i{0U}; i<(*eit)->FV()->Facets(); i++ )
        {
            // identifying the finite volumes to which the flux will be distributed
            (*eit)->FV()->FacetEdgeNodes( i, inside_node, outside_node );

            // projecting velocities onto normals to segments and integrating over area       (ip)
            flux = ( *(*eit) ).ProjectionOnFacetNormal( i, velo );
            area = ( *(*eit) ).FacetArea( i );

            cout <<"("<< i+1 <<") "<< flux <<" "<< area <<" ["<< inside_node <<","<< outside_node <<"] ";

            flux *= area;

            // incoming & outgoing fluxes are added to first and second member of pair, respectively
            // (resulting values are all stored as positive)
            if ( flux < zero ) {
                THROUGHPUT[ (*eit)->N(inside_node)->Idx() ].first   += fabs(flux);
                THROUGHPUT[ (*eit)->N(outside_node)->Idx() ].second += fabs(flux);
            }
            else  {
                THROUGHPUT[ (*eit)->N(outside_node)->Idx() ].first  += fabs(flux);
                THROUGHPUT[ (*eit)->N(inside_node)->Idx() ].second  += fabs(flux);
            }
        }
        cout << endl;


        // 3. testing the interpolation to facet integration points
        // --------------------------------------------------------
        cout <<"\nProperty values interpolated to facet integration points:  "<< endl <<"\t";
        for ( auto i{0U}; i<(*eit)->FV()->Facets(); i++ ) {
            // testing sum of interpolation functions at integration point
            ( *(*eit) ).N_At( (*eit)->FV()->FacetIntegrationPoint( i, 0U ), IPOL );
            if ( fabs(accumulate( IPOL.begin(), IPOL.end(), 0. ) - 1.0) > 1.0e-15 )
                cout <<"\nERROR: at ip of facet "<< i <<", interpolation functions did not sum to one."<< endl;

            // trying property interpolation
            (*eit)->NodePropertyVector( prop_key, node_prop );
            sum = zero;
            for ( auto j{0U}; j<(*eit)->Nodes(); j++ ) sum += IPOL[j] * node_prop[j]();
            cout <<"("<< i+1 <<") "<< sum <<" ";

        }
        cout << endl;

        // 4. testing the sector volume computations
        // -----------------------------------------
        cout <<"\nSector volumes: "<< endl <<"\t";
        sum = zero;
        for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) {
            cout <<"("<< i+1 <<") "<< (val=( *(*eit) ).SectorVolume(i)) <<" ";
            sum += val;
        }
        if ( fabs(sum-evolume) > numeric_limits<double>::epsilon()*10. )
            cout <<"\nERROR: Sector volumes sum up to a different value than the element volume (evol-sum): "<<  evolume-sum << endl;

        cout <<"\n\n\n";
    }


    // 2. calculating the normalized flux balance for each finite volume cell
    // -------------------------------------------------------------------------------------
    double   overall_flux_balance(zero), flux_balance;
    size_t   cells_counted(0);

    cout <<"\ntestFiniteVolumeStencil(NodeCenteredFiniteVolumeTransport): Testing flux-balance in FV cells away from model boundary; ";
    cout <<" cell idx and corresponding balance normalized by FV cell volume: "<< endl;
    for ( size_t i{0U}; i<THROUGHPUT.size(); i++ ) {
        // only for cells that are not located at the model boundary
        if ( gref.N(i)->AtBoundary() == NOT )
        {
            flux_balance = fabs(THROUGHPUT[i].first - THROUGHPUT[i].second) / ncvft.PoreVolume(i);
            overall_flux_balance = std::max( flux_balance, overall_flux_balance );
            cout << i+1 <<": "<< flux_balance <<" ";
            cells_counted++;
        }
        // zeroing vector<double> for next call
        THROUGHPUT[i].first = THROUGHPUT[i].second = zero;
    }

    if ( cells_counted > 0U )
        cout <<"\n\nThe overall normalized balance was: "<< overall_flux_balance << endl << endl;
    else
        cout <<"\n\ntest could not be performed because all cells were located at model boundaries."<< endl << endl;

    cout.precision(prec);

    return true;

} // end testFiniteVolumeStencil


template
bool testFiniteVolumeStencil( const PropertyDatabase<1>&, const Region<1>&,
NodeCenteredFiniteVolumeTransport<1U>&, long );
template
bool testFiniteVolumeStencil( const PropertyDatabase<2>&, const Region<2>&,
NodeCenteredFiniteVolumeTransport<2U>&, long );
template
bool testFiniteVolumeStencil( const PropertyDatabase<3>&, const Region<3>&,
NodeCenteredFiniteVolumeTransport<3U>&, long );


template class NodeCenteredFiniteVolumeTransport<1U>;
template class NodeCenteredFiniteVolumeTransport<2U>;
template class NodeCenteredFiniteVolumeTransport<3U>;

} // end namespace csmp

