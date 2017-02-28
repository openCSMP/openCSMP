#include "ExplicitFiniteVolumeTransportPHX.h"
//#include "ElementFiniteVolumeTraits.h"
//#include "CSP_ErrorHandler.h"
//#include "CSP_STL_utilities.h"
//#include <numeric>

using namespace std;

namespace csmp {


/** custom constructor */
template<size_t dim>
ExplicitFiniteVolumeTransportPHX<dim>::ExplicitFiniteVolumeTransportPHX( Model<dim>& model,
                                                  std::vector<std::string>& lhs_property,
                                                  std::vector<std::string>& rhs_property,
                                                  const char* velocity,
                                                  const char* source,
                                                  const char* density)
 : 
   NodeCenteredFiniteVolumeTransport<dim>( "Model", model, "porosity",
                                              lhs_property[0].c_str(), velocity, 
                                              source, false, false ),
   model_ref( model ),
   region_ref( model.Region("Model") ),
   flux_in_vectors(lhs_property.size()),
   flux_out_vectors(lhs_property.size()),
   property_vectors(lhs_property.size()),
   mass_balance_vectors(lhs_property.size()),
   facet_flux_vectors(lhs_property.size()),
   max_time_step(8640000.),
   with_gravity(false)
  {
   
   cout << "\n Constructing ExplicitFiniteVolumeTransportPHX " << endl;

   // resize flux and property vectors
   for ( size_t i=0; i<lhs_property.size(); i++ )
       {
        flux_in_vectors[i].resize(region_ref.Nodes());
        flux_out_vectors[i].resize(region_ref.Nodes());
        property_vectors[i].resize(region_ref.Nodes());
        mass_balance_vectors[i].resize(region_ref.Nodes());
        facet_flux_vectors[i].resize( region_ref.Elements() );
        for ( typename vector<Element<dim>*>::const_iterator 
               eit=region_ref.ElementsBegin(); 
               eit!=region_ref.ElementsEnd(); eit++)
              facet_flux_vectors[i][(*(*eit)).Idx()].resize( (*(*eit)).FV()->Facets() );
       }
    
   // set and check keys
    k_key  =  model.Database().StorageKey( "permeability" );
    pv_key  =  model.Database().StorageKey( "pore volume" );
    rho_key  =  model.Database().StorageKey( density );

    pl_key.resize( lhs_property.size() );
    pr_key.resize( rhs_property.size() );

    for ( size_t i=0; i<lhs_property.size(); i++ ) {
	    pl_key[i] =  model.Database().StorageKey( lhs_property[i].c_str() );
	    pr_key[i] =  model.Database().StorageKey( rhs_property[i].c_str() );
	    
	    if ( pl_key[i].type != SCALAR || pl_key[i].place != NODE )
           throw csmp::Exception( ERROR, "ExplicitFiniteVolumeTransportPHX.<double64, dim>::ExplicitFiniteVolumeTransportPHX\n", 
                             lhs_property[i].c_str(), " must be a nodal scalar property." );

	    if ( pr_key[i].type != SCALAR || (pr_key[i].place != NODE && pr_key[i].place != ELEMENT) )
           throw csmp::Exception( ERROR, "ExplicitFiniteVolumeTransportPHX.<double64, dim>::ExplicitFiniteVolumeTransportPHX\n", 
                             rhs_property[i].c_str(), " must be a nodal or element scalar property." );

   }     

    cout <<"\n\ExplicitFiniteVolumeTransportPHX<"<< typeid(double64).name() <<","<< dim;
    cout <<">: Constructed successfully."<< endl;
    
 } // end constructor 


/** default destructor */
template<size_t dim>
ExplicitFiniteVolumeTransportPHX<dim>::~ExplicitFiniteVolumeTransportPHX()
 {
 }



/** change lhs variable by adding and subtracting flux in and out of control volume */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::ComposeAdvection( )
 {

     for ( size_t nidx=0U; nidx<region_ref.Nodes(); nidx++ )
       for ( size_t i=0; i<flux_in_vectors.size(); i++ )
          {
            property_vectors[i][nidx] -= flux_in_vectors[i][nidx];
            property_vectors[i][nidx] -= flux_out_vectors[i][nidx];
          }
       
 } // end ComposeAdvection                                               


/** store lhs variables */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::WriteResults( )// const
 {
    double64                   pmin(0), pmax(0), value(0);
    ScalarVariable   sc;
    
    for ( size_t i=0; i<property_vectors.size(); i++ )
      {
        model_ref.Database().RangeOf( model_ref.Database().Name(pl_key[i]), pmin, pmax );
	    for ( size_t idx=0; idx<property_vectors[i].size(); idx++ )
	       {
	        value = property_vectors[i][idx];

	        // store property
		    if ( region_ref.N(idx)->Status( pl_key[i] ) != DIRICH )
		      {
  		       // reading the pre-existing value (use it to keep the flag if results are overwritten)
  		       region_ref.N(idx)->Read( pl_key[i], sc );
  		       if ( value <= pmax && value >= pmin )
  		         {
  		           sc()=value;
  		           region_ref.N(idx)->Store( pl_key[i], sc );
  		         }
  		       else
  		         {
  		           int temp;
  		           cerr <<"\nproperty value: "<< value <<" at FV " << idx+1 <<" versus range from PropertyDatabase: "<< pmin <<"-"<< pmax << endl;
  		           if (value < pmin) sc()= pmin;
  		           else if (value > pmax) sc()= pmax;
  		           region_ref.N(idx)->Store( pl_key[i], sc );
                       throw csmp::Exception( ERROR, "ExplicitFiniteVolumeTransportPHX::WriteResults",
                                     "Output property was out of range, legal (min/max) was stored instead");
  		           cin >> temp;
  		         }
  	          }
	      } // end finite volumes
  
    } // end property vectors
 } // end WriteResults

/** change size of maximum time step */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::SetMaximumTimeStep( const double64& new_max_time_step )
   {
     max_time_step = new_max_time_step;
   }

/** access to flux out of control volume - main property */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFluxOut(unsigned int index)
   {
    return flux_out_vectors[0][index];
   }
   
/** access to flux into control volume - main property */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFluxIn(unsigned int index)
   {
    return flux_in_vectors[0][index];
   }

/** access to flux out of control volume - specified property */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFluxOut(unsigned int index, unsigned int property_idx)
   {
    return flux_out_vectors[property_idx][index];
   }

/** access to flux into control volume - specified property */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFluxIn(unsigned int index, unsigned int property_idx)
   {
    return flux_in_vectors[property_idx][index];
   }

/** get value of main property at specified control volume */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetPropertyValue(unsigned int index)
   {
    return region_ref.N(index)->Read( pl_key[0] );
   }

/** get facte flux of specified propoerty at indicated facet */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFacetFlux( Element<dim>& e, unsigned int facet_idx, unsigned int property_idx )
   {
    size_t eidx_,inside_node_,outside_node_;
    double64 flux_;
    eidx_ = e.Idx();
    e.FV()->FacetEdgeNodes( facet_idx, inside_node_, outside_node_ );
    flux_ = facet_flux_vectors[property_idx][eidx_][facet_idx];
    if (flux_<0)
      flux_ *= mass_balance_vectors[property_idx][e.N(outside_node_)->Idx()];
    else
      flux_ *= mass_balance_vectors[property_idx][e.N( inside_node_)->Idx()];
    return flux_;
   }

/** get facte flux of specified propoerty at indicated facet */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetMainPropertyLHS(unsigned int index)
   {
    return property_vectors[0][index];
   }

/** calculate and store facet fluxes with suggested time step */
template<size_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::DetermineFacetFlux( const double64& time_increment, std::vector<DenseMatrix<DM_MIN> >& upwind_visitor )
  {
  
    internal_time_step = time_increment;
   		    
    vector<FV_Parameter>::const_iterator      fvt;
    typename vector<Element<dim>*>::const_iterator  eit;

    for ( size_t i=0; i<flux_in_vectors.size(); i++ )
      {
       fill(    flux_in_vectors[i].begin(),    flux_in_vectors[i].end(), static_cast<double64>(0.) );
       fill(   flux_out_vectors[i].begin(),   flux_out_vectors[i].end(), static_cast<double64>(0.) );
       fill(   property_vectors[i].begin(),   property_vectors[i].end(), static_cast<double64>(0.) );
       fill(mass_balance_vectors[i].begin(),mass_balance_vectors[i].end(), static_cast<double64>(1.) );
       for ( size_t j=0; j<facet_flux_vectors[i].size(); j++ )
           fill( facet_flux_vectors[i][j].begin(), facet_flux_vectors[i][j].end(), static_cast<double64>(0.) );       
      }
    
    for ( size_t i=0; i<flux_out_vectors.size(); i++ )
      {
        fvt = NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA.begin();
        for ( eit=region_ref.ElementsBegin(); eit!=region_ref.ElementsEnd(); eit++, fvt++ )
           {
            GetUpwindMatrix( *(*eit), upwind_visitor );
            if (with_gravity)
              stencilPHX.DetermineFluxOutWithGravity( (*fvt), *(*eit), upwind,
                                                       facet_flux_vectors[i], flux_out_vectors[i],
                                                       pr_key[i], rho_key, k_key);
            else
              stencilPHX.DetermineFluxOutWithoutGravity( (*fvt), *(*eit), upwind,
                                                        facet_flux_vectors[i], flux_out_vectors[i],
                                                        pr_key[i]);
           }
       } 

    CalculateOutflowPerPoreVolume();

 } // end DetermineFacetFlux

/** calculate total flux out of control volumes */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::CalculateOutflowPerPoreVolume()
   {

    double64 pore_vol;

    for ( typename vector<Node<dim>*>::const_iterator 
          fvit=region_ref.NodesBegin();
          fvit!=region_ref.NodesEnd(); fvit++ )
          { 
            size_t idx = (*fvit)->Idx();
   		      pore_vol = (*fvit)->Read( pv_key);
            for( size_t i=0; i<flux_out_vectors.size(); i++ ) {
                property_vectors[i][idx] = (*fvit)->Read( pl_key[i] );
                flux_out_vectors[i][idx] *= internal_time_step;
                flux_out_vectors[i][idx] /= pore_vol; 
   		     } // end flux vectors
          } // end finite volumes

   } // end CalculateOutflowPerPoreVolume

/** adjust fluxes, perform advection calculations and store results*/
template<size_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::AdjustAndPerformFacetFlux( const double64& time_factor ){
   	
   	internal_time_step *= time_factor;
   	
    AdjustFluxOut( time_factor );
    
    CalculateFluxIn();
    
    CalculateInflowPerPoreVolume();
    
    ComposeAdvection();
   
    WriteResults();
     
     } // end AdjustAndPerformFacetFlux

/** adjust fluxes */
template<size_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::AdjustFluxOut( const double64& time_factor ){
   	
   	
   	double64 mass_balance_factor(1.0);
   	
    for ( typename vector<Node<dim>*>::const_iterator 
          fvit=region_ref.NodesBegin(); 
          fvit!=region_ref.NodesEnd(); fvit++ )
          { 
            unsigned int idx = (*fvit)->Idx();
            for ( size_t i=0; i<flux_out_vectors.size(); i++ )
               {
                 flux_out_vectors[i][idx] *= time_factor;
                 if (flux_out_vectors[i][idx] > 0.)
                     mass_balance_factor = property_vectors[i][idx] / flux_out_vectors[i][idx];
                 else mass_balance_factor = 1.0;
                     
                 if (mass_balance_factor< 1.0)
                    {
                      mass_balance_vectors[i][idx] = mass_balance_factor;
                      flux_out_vectors[i][idx] = property_vectors[i][idx];                      
                    }
                 }
          }

     } // end AdjustFluxOut

/** calculate flux into control volumes */
template<size_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::CalculateFluxIn( )
    {
   	
    vector<FV_Parameter>::const_iterator      fvt;
    typename vector<Element<dim>*>::const_iterator  eit;
    
    for ( size_t i=0; i<flux_in_vectors.size(); i++ )
      {
        fvt = NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA.begin();
        for ( eit=region_ref.ElementsBegin(); eit!=region_ref.ElementsEnd(); eit++, fvt++ )
            stencilPHX.DetermineFluxIn( *(*eit),
                                        facet_flux_vectors[i], flux_in_vectors[i],
                                        mass_balance_vectors[i]);
       } 
    } // end CalculateFluxIn

/** calculate total flux into control volumes */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::CalculateInflowPerPoreVolume()
   {

    double64 pore_vol;

    for ( typename vector<Node<dim>*>::const_iterator 
          fvit=region_ref.NodesBegin(); 
          fvit!=region_ref.NodesEnd(); fvit++ )
          { 
            unsigned int idx = (*fvit)->Idx();
   		    pore_vol = (*fvit)->Read( pv_key );
            for ( size_t i=0; i<flux_in_vectors.size(); i++ )
               {
                flux_in_vectors[i][idx] *= internal_time_step;
                flux_in_vectors[i][idx] /= pore_vol; 
   		       } // end flux vectors
          } // end finite volumes

   } // end CalculateInflowPerPoreVolume

/** get projected velocity for specified facet */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetProjectedVelocities( Element<dim>& e, size_t i)
   {
     return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[e.Idx()].FacetNormalVelocity(i)*NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[e.Idx()].FacetArea(i);
   }

/** activate gravity component */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::WithGravityComponent( )
   {
       with_gravity = true;      
   }

/** get matrix for upwind nodes */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::GetUpwindMatrix( const Element<dim>& e, std::vector<DenseMatrix<DM_MIN> >& upwind_visitor )
{

  upwind = upwind_visitor[e.Idx()];    

} // GetUpwindMatrix

/** get index of phase density */
template<size_t dim>
csmp::Index ExplicitFiniteVolumeTransportPHX<dim>::GetDensityKey(  )
{

  return rho_key;    

} // GetDensityKey

/** update projection of velocity onto the facet */
template<size_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::UpdateProjection(  )
{
    NodeCenteredFiniteVolumeTransport<dim>::UpdateProjectedVelocitiesAndFluxBalances( );

} // UpdateProjection

/** update projection of velocity onto the facet */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFacetNormalVelocity( size_t element, size_t facet )
{
   return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[element].FacetNormalVelocity(facet);
} // GetFacetNormalVelocity

/** get normal component of indicated facet and direction */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFacetNormalComponent( size_t element, size_t facet, size_t x_or_y_or_z )
{
   return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[element].FacetNormalComponent(facet, x_or_y_or_z);
} // GetFacetNormalComponent

/** get area of indicated facet */
template<size_t dim>
double64 ExplicitFiniteVolumeTransportPHX<dim>::GetFacetArea( size_t element, size_t facet )
{
   return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[element].FacetArea(facet);
} // GetFacetArea

/** add new advection variables for finite volume calculations */
template<size_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::AddAdvectionVariable( const char* new_lhs, const char* new_rhs)
{

   pl_key.push_back(model_ref.Database().StorageKey( new_lhs ));
   pr_key.push_back(model_ref.Database().StorageKey( new_rhs ));

   flux_in_vectors.resize(pl_key.size());
   flux_out_vectors.resize(pl_key.size());
   property_vectors.resize(pl_key.size());
   mass_balance_vectors.resize(pl_key.size());
   facet_flux_vectors.resize(pl_key.size());

   size_t i = pl_key.size()-1;
   flux_in_vectors[i].resize(region_ref.Nodes());
   flux_out_vectors[i].resize(region_ref.Nodes());
   property_vectors[i].resize(region_ref.Nodes());
   mass_balance_vectors[i].resize(region_ref.Nodes());
   facet_flux_vectors[i].resize( region_ref.Elements() );
   for ( typename vector<Element<dim>*>::const_iterator
         eit=region_ref.ElementsBegin();
         eit!=region_ref.ElementsEnd(); eit++)
       facet_flux_vectors[i][(*(*eit)).Idx()].resize( (*(*eit)).FV()->Facets() );

}


/** single-phase version of mass balanced finite volume calculations */
template<size_t dim>
double64  ExplicitFiniteVolumeTransportPHX<dim>::AdvectMassConservedSinglePhase( const double64& time_increment )
 {

   internal_time_step = time_increment;

   // project velocities to face normals
   UpdateProjection();

   // determine facet fluxes and flux out
   // get suggestions for time step adjustments

   DetermineFacetFluxSinglePhase( );

   CheckDryFiniteVolumesAndAdjustTimestepForSinglePhase( );

   // determine flux in and perform and write fluxes
   AdjustAndPerformFacetFlux( single_phase_time_step_factor );

   // return the time step used for the advection
   internal_time_step *= single_phase_time_step_factor;

   return internal_time_step;

 } // end AdvectMassConserved

/** mass-balance check, reducing time step if criterion is not met */
template<size_t dim>
void  ExplicitFiniteVolumeTransportPHX<dim>::CheckDryFiniteVolumesAndAdjustTimestepForSinglePhase( )
 {

    // check!

    single_phase_time_step_factor = 1.0;

    unsigned int idx;
    double64 LHS, Outflow, temp_factor;

    // time step is only cut if the primary variable (liquid + vapor mass) is running dry
    for ( typename vector<Node<dim>*>::const_iterator fvit = region_ref.NodesBegin();
         fvit != region_ref.NodesEnd(); fvit++ )
         {
           idx = (*fvit)->Idx();

           LHS      = GetMainPropertyLHS(idx);
           Outflow  = GetFluxOut(idx);
           if (Outflow != 0.) temp_factor = LHS / Outflow;
           else temp_factor = 1.;

           if ( temp_factor < single_phase_time_step_factor)
             {
               single_phase_time_step_factor = temp_factor;
             }
         }

 } // end CheckDryFiniteVolumesAndAdjustTimestepForSinglePhase

/** calculate facet flux for single-phase version */
template<size_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::DetermineFacetFluxSinglePhase( )
  {

    vector<FV_Parameter>::const_iterator      fvt;
    typename vector<Element<dim>*>::const_iterator  eit;

    for ( size_t i=0; i<flux_in_vectors.size(); i++ )
      {
       fill(    flux_in_vectors[i].begin(),    flux_in_vectors[i].end(), static_cast<double64>(0.) );
       fill(   flux_out_vectors[i].begin(),   flux_out_vectors[i].end(), static_cast<double64>(0.) );
       fill(   property_vectors[i].begin(),   property_vectors[i].end(), static_cast<double64>(0.) );
       fill(mass_balance_vectors[i].begin(),mass_balance_vectors[i].end(), static_cast<double64>(1.) );
       for ( size_t j=0; j<facet_flux_vectors[i].size(); j++ )
           fill( facet_flux_vectors[i][j].begin(), facet_flux_vectors[i][j].end(), static_cast<double64>(0.) );
      }

    for ( size_t i=0; i<flux_out_vectors.size(); i++ )
      {
        fvt = NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA.begin();
        for ( eit=region_ref.ElementsBegin(); eit!=region_ref.ElementsEnd(); eit++, fvt++ )
              stencilPHX.DetermineFluxOut((*fvt), *(*eit),
                                          facet_flux_vectors[i], flux_out_vectors[i],
                                          pr_key[i]);
      }

    CalculateOutflowPerPoreVolume();

 } // end DetermineFacetFluxSinglePhase

template class ExplicitFiniteVolumeTransportPHX<1U>;
template class ExplicitFiniteVolumeTransportPHX<2U>;
template class ExplicitFiniteVolumeTransportPHX<3U>;
 
} // end namespace csmp
 







   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
