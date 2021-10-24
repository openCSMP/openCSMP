#include "DES2PhaseTransport.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#include "FlowFunctionsModule.h"
#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
DES2PhaseTransport<dim,FLOW_FUNCTIONS>::DES2PhaseTransport(  Model<dim>& m,
                                                             const char* target_region,
                                                             bool with_gravity_forces,
                                                             bool with_capillary_spreading,
                                                             double64 cfl_multiplier,
                                                             double64 PEP_multiplier, 
                                                             double64 relaxing_factor,                                             
                                                             bool tensor_k,
                                                             FLOW_FUNCTIONS<dim>& ff )
  : sg_(m), 
    gref_(m.Region(target_region)),
    db_(m.Database()),
    with_capillary_spreading_(with_capillary_spreading),
    with_gravity_forces_(with_gravity_forces),    
    CFL_multiplier_(cfl_multiplier), 
    PEP_multiplier_(PEP_multiplier),
    relaxing_factor_(relaxing_factor),
    tensor_k_(tensor_k),
    flowfunctions_(ff),
    rate_count_(0U), update_count_(0U),
    T_RateOfChange_(0.), T_Schedule_(0.), T_InsertToHeap_(0.), T_Update_(0.), T_Synchronize_(0.), T_RemoveFromHeap_(0.), T_AdvectVariable_(0.)
{
    InitializeBasicVariablsAndKeys();
    InitializeFiniteVolumeProperties();
    cout<<"DES2PhaseTransport constructed"<<endl;
} // end constructor 




template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseTransport<dim,FLOW_FUNCTIONS>::InitializeBasicVariablsAndKeys()
{
    //creating new variables if not defined yet from input file   
    if(!db_.IsDefined("event index")) sg_.CreateProperty( "event index", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!db_.IsDefined("update count")) sg_.CreateProperty( "update count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);  
    if(!db_.IsDefined("rate count")) sg_.CreateProperty( "rate count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!db_.IsDefined("schedule count")) sg_.CreateProperty( "schedule count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!db_.IsDefined("synchronize count")) sg_.CreateProperty( "synchronize count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!db_.IsDefined("DES array")) sg_.CreateProperty( "DES array", "none", ARRAY, NODE, 7, -1.00E+10 ,1.00E+10);
    if(!db_.IsDefined("porosity")) sg_.CreateProperty( "porosity", "none", SCALAR, ELEMENT, 1, 1.00E-05, 1.00E+01); 
    if(!db_.IsDefined("thickness")) sg_.CreateProperty( "thickness", "m", SCALAR, ELEMENT, 1, 0.0E+0, 1.00E+10); 
    if(!db_.IsDefined("facet area")) sg_.CreateProperty( "facet area", "m2", SCALAR, FACET_INTEGRATION_POINT, 1, -1.00E+10, 1.00E+10); 
    if(!db_.IsDefined("facet normal")) sg_.CreateProperty( "facet normal", "m2 s-1", VECTOR, FACET_INTEGRATION_POINT, 3, 0., 1.);
    if(!db_.IsDefined("FV pore volume")) sg_.CreateProperty( "FV pore volume", "m3", SCALAR, NODE, 1, 0.00E+00 ,1.00E+8); 
    if(!db_.IsDefined("saturation carbonic phase")) sg_.CreateProperty( "saturation carbonic phase", "m3/m3", SCALAR, NODE, 1, 0 ,1);
    if(!db_.IsDefined("saturation aqueous phase")) sg_.CreateProperty( "saturation aqueous phase", "m3/m3", SCALAR, NODE, 1, 0 ,1);
    if(!db_.IsDefined("shock saturation aqueous phase")) sg_.CreateProperty( "shock saturation aqueous phase", "none", SCALAR, ELEMENT, 1, -5.00E-02 ,1.05E+00); 
    if(!db_.IsDefined("truncated FV")) sg_.CreateProperty( "truncated FV", "none", SCALAR, NODE, 1, 0 ,1);
    if(!db_.IsDefined("cfl multiplier")) sg_.CreateProperty( "cfl multiplier", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10); 
    if(!db_.IsDefined("saturation gradient")) sg_.CreateProperty( "saturation gradient", "none", VECTOR, ELEMENT, 3, -1.00E+08 ,1.00E+08);
    if(!db_.IsDefined("pressure gradient")) sg_.CreateProperty( "pressure gradient", "none", VECTOR, ELEMENT, 3, -1.00E+10 ,1.00E+10);
    if(!db_.IsDefined("fluid pressure")) sg_.CreateProperty( "fluid pressure", "Pa", SCALAR, NODE, 1, 0 ,1.00E+10);
    if(!db_.IsDefined("permeability")) sg_.CreateProperty( "permeability", "m2", SCALAR, ELEMENT, 1, 1E-21, 1.0e-7);
    
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
    key_fvPV = INDEX<SCALAR,NODE>( db_.StorageKey("FV pore volume") );
    key_sCO2 = INDEX<SCALAR,NODE>( db_.StorageKey("saturation carbonic phase") );
    key_sH2O = INDEX<SCALAR,NODE>( db_.StorageKey("saturation aqueous phase") );
    key_ssH2O = INDEX<SCALAR,ELEMENT>( db_.StorageKey("shock saturation aqueous phase") );
    key_cut = INDEX<SCALAR,NODE> ( db_.StorageKey("truncated FV") );
    key_CFL = INDEX<SCALAR,NODE>( db_.StorageKey("cfl multiplier") );
    key_gradSn = INDEX<VECTOR,ELEMENT>( db_.StorageKey("saturation gradient") );
    key_gradP = INDEX<VECTOR,ELEMENT>( db_.StorageKey("pressure gradient") );
    key_pf = INDEX<SCALAR,NODE>( db_.StorageKey("fluid pressure") );
    key_k = INDEX<SCALAR,ELEMENT>( db_.StorageKey("permeability") );
    
    //checking keys
    if ( key_EventIndex.place != NODE || key_EventIndex.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'event index' variable must be SCALAR and placed on NODE"  );       
    if ( key_update.place != NODE || key_update.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( key_rate.place != NODE || key_rate.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( key_schedule.place != NODE || key_schedule.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( key_synchronize.place != NODE || key_synchronize.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  );    
    if ( key_time.place != NODE || key_time.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'DES array' variable must be ARRAY and placed on NODE"  );   
    if ( key_phi.place != ELEMENT || key_phi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );    
    if ( key_thi.place != ELEMENT || key_thi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  ); 
    if ( key_fA.place != FACET_INTEGRATION_POINT || key_fA.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );       
    if ( key_fn.place != FACET_INTEGRATION_POINT || key_fn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  ); 
    if ( key_fvPV.place != NODE || key_fvPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );    
    if ( key_sCO2.place != NODE || key_sCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'saturation carbonic phase' variable must be SCALAR and placed on NODE"  );  
    if ( key_sH2O.place != NODE || key_sH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'saturation aqueous phase' variable must be SCALAR and placed on NODE"  );                                       
    if ( key_ssH2O.place != ELEMENT || key_ssH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'shock saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  ); 
    if ( key_cut.place != NODE || key_cut.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'truncated FV' variable must be SCALAR and placed on NODE"  );  
    if ( key_CFL.place != NODE || key_CFL.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'cfl multiplier' variable must be SCALAR and placed on NODE"  );
    if ( key_gradSn.place != ELEMENT || key_gradSn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'saturation gradient' variable must be VECTOR and placed on ELEMENT"  ); 
    if ( key_gradP.place != ELEMENT || key_gradP.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'pressure gradient' variable must be VECTOR and placed on ELEMENT"  ); 
    if ( key_pf.place != NODE || key_pf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'fluid pressure' variable must be SCALAR and placed on NODE"  );         
    if ( key_k.place != ELEMENT || key_k.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "DES2PhaseTransport::initializeKeys:",
        "The 'permeability' variable must be SCALAR and placed on ELEMENT"  );

    // model-wide initialisation
    Region<dim> region(sg_.Region("Model"));
    region.InputPropertyValue( "update count", makeScalar(PLAIN,0.), COMPLETE );
    region.InputPropertyValue( "rate count", makeScalar(PLAIN,0.), COMPLETE );
    region.InputPropertyValue( "schedule count", makeScalar(PLAIN,0.), COMPLETE );
    region.InputPropertyValue( "synchronize count", makeScalar(PLAIN,0.), COMPLETE );
    region.InputPropertyValue( "DES array", ArrayVariable(7,0.,PLAIN), COMPLETE);
    region.InputPropertyValue( "truncated FV", makeScalar(PLAIN,0), COMPLETE);
    db_.RangeOf( db_.Name(key_sCO2), lower_limit_, upper_limit_ );    
}


template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseTransport<dim,FLOW_FUNCTIONS>::InitializeFiniteVolumeProperties()
 {
    //zeroing FV pore volumes for accumulation in element loop
    gref_.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE ); 
    
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref_.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref_.ElementsBegin(); it!=it_end; ++it )
    {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // computing sector pore volumes
         double64 phi = (*it)->Read( key_phi);
         const double64 thickness = (*it)->Read( key_thi );
         if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
         
         for ( size_t i=0U; i<sectors; ++i ) {
              const double64 sector_volume = (*it)->SectorVolume(i);  
              double64 pore_volume   = (*it)->N(i)->Read( key_fvPV );
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( key_fvPV, makeScalar(PLAIN,pore_volume) );
         }

         // computing facet normals and areas
         for ( size_t j=0U; j<facets; ++j ) {
              const double64 facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, key_fA, makeScalar( PLAIN, facet_area ) );
              Point<dim> nrml = (*it)->FacetNormal(j);
              VectorVariable<dim>  fnrml;
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, key_fn, fnrml );
         }     
         
         //compute wetting phase saturation at shock
         double64 sw_shock = flowfunctions_.ShockHeight(*it);
         (*it)->Store( key_ssH2O, makeScalar( (*it)->Status( key_ssH2O), sw_shock ) );                  
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
                  eptr->Store( j, 0U, key_fA, makeScalar( PLAIN, facet_area ) );
                  Point<dim> nrml = eptr->FacetNormal(j);
                  VectorVariable<dim>  fnrml;
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, key_fn, fnrml );
             }    
             
             //compute wetting phase saturation at shock
             double64 sw_shock = flowfunctions_.ShockHeight(eptr);
             eptr->Store( key_ssH2O, makeScalar( eptr->Status( key_ssH2O), sw_shock ) );  
               
             //determine whether FV node is truncated
             if ((*nit)->AtBoundary() == NOT) { 
                 if(!gref_.Contains(eptr)) { //parent elment located outside domain
                     halo_stencils_.insert(eptr);
                     truncated_node = true;
                 }  
             }                        
        }
        if (truncated_node) (*nit)->Store( key_cut, makeScalar( (*nit)->Status( key_cut), 1 ) );
   }
 } // end initializeFiniteVolumeProperties




//resest cfl multipliers to default value = CFL_multiplier_*relaxing_factor_ for all nodes
template<size_t dim, template<size_t> class FLOW_FUNCTIONS>
void DES2PhaseTransport<dim,FLOW_FUNCTIONS>::ResetCFLMultiplier()
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


template class DES2PhaseTransport<1U,FlowFunctionsModule1>;
template class DES2PhaseTransport<2U,FlowFunctionsModule1>;
template class DES2PhaseTransport<3U,FlowFunctionsModule1>;

template class DES2PhaseTransport<1U,FlowFunctionsModule2>;
template class DES2PhaseTransport<2U,FlowFunctionsModule2>;
template class DES2PhaseTransport<3U,FlowFunctionsModule2>;

template class DES2PhaseTransport<1U,FlowFunctionsModule3>;
template class DES2PhaseTransport<2U,FlowFunctionsModule3>;
template class DES2PhaseTransport<3U,FlowFunctionsModule3>;

template class DES2PhaseTransport<1U,FlowFunctionsModule4>;
template class DES2PhaseTransport<2U,FlowFunctionsModule4>;
template class DES2PhaseTransport<3U,FlowFunctionsModule4>;

template class DES2PhaseTransport<1U,FlowFunctionsModule5>;
template class DES2PhaseTransport<2U,FlowFunctionsModule5>;
template class DES2PhaseTransport<3U,FlowFunctionsModule5>;

template class DES2PhaseTransport<1U,FlowFunctionsModule6>;
template class DES2PhaseTransport<2U,FlowFunctionsModule6>;
template class DES2PhaseTransport<3U,FlowFunctionsModule6>;

template class DES2PhaseTransport<1U,FlowFunctionsModule7>;
template class DES2PhaseTransport<2U,FlowFunctionsModule7>;
template class DES2PhaseTransport<3U,FlowFunctionsModule7>;

} // end csmp  
    
    
