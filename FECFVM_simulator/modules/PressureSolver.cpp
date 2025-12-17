#include "PressureSolver.h"
#include "LinearSolver.h"
#include "Model.h"
#include "Exception.h"
#include "meshManagementUtilities.h"

using namespace std;

namespace csmp {

#ifdef CSMP_WITH_SAMG_SOLVER
template<uint32_t dim, template<uint32_t> class CELL>
PressureSolver<dim, CELL>::PressureSolver( Model<dim>& sg, bool with_gravity, bool with_capillary )
:  solver_(SAMG_Solver(&settings_)),
   PDE_Integrator<dim, CELL>(solver_),
   total_mobility_( sg.Database(), "total mobility", "fluid pressure", "fluid pressure" ),
   nodal_source_( sg.Database(), "nodal fluid volume source", "fluid pressure" ),
   gravity_term_( sg.Database(), "gravity term", "fluid pressure" ),
   capillary_term_( sg.Database(), "capillary term", "fluid pressure" )
{
   	if ( !isoparametricElementMesh( sg ) )
  		throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):", 
  		                                    "elements are not isoparametric; this works only for isoparametric elements." );

    const PropertyDatabase<dim>&   p_ref = sg.Database();
    
    csmp::Index  mob_key = p_ref.StorageKey("total mobility");
    if ( mob_key.place != ELEMENT || mob_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):", 
                      "total mobility", "variable must be  a scalar placed on the element." ); 

    csmp::Index  pressure_key = p_ref.StorageKey("fluid pressure");
    if ( pressure_key.place != NODE || pressure_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):", 
                      "fluid pressure", "variable must be a scalar placed on the node." ); 

    csmp::Index  nsource_key = p_ref.StorageKey("nodal fluid volume source");
    if ( nsource_key.place != NODE  and  nsource_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):", 
                      "nodal fluid volume source", "variable must be a scalar placed on the node." ); 

    csmp::Index  gterm_key = p_ref.StorageKey("gravity term");
    if ( gterm_key.place != ELEMENT  and  gterm_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):", 
                      "gravity term", "variable must be a vector placed on the element." ); 
          
    csmp::Index  cterm_key = p_ref.StorageKey("capillary term");
    if ( cterm_key.place != ELEMENT  and  cterm_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):", 
                      "capillary term", "variable must be a vector placed on the element." );
    
    settings_.Set_idmp(0);
    settings_.Set_iout1(0);

    this->Add( &total_mobility_ );
    this->Add( &nodal_source_ );
    if ( with_gravity ) this->Add( &gravity_term_ );
    if ( with_capillary )this->Add( &capillary_term_ );
}

template<uint32_t dim, template<uint32_t> class CELL>
SAMG_Settings& PressureSolver<dim, CELL>::PressureSolverSettings()
{
    return settings_;
}

#else /* !CSMP_WITH_SAMG_SOLVER */


template<uint32_t dim, template<uint32_t> class CELL>
PressureSolver<dim, CELL>::PressureSolver( Model<dim>& sg, bool with_gravity, bool with_capillary )
:  PDE_Integrator<dim, CELL>(solver_),
   total_mobility_( sg.Database(), "total mobility", "fluid pressure", "fluid pressure" ),
   nodal_source_( sg.Database(), "nodal fluid volume source", "fluid pressure" ),
   gravity_term_( sg.Database(), "gravity term", "fluid pressure" ),
   capillary_term_( sg.Database(), "capillary term", "fluid pressure" )
{

   	if ( !isoparametricElementMesh( sg ) )
  		throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):",
  		                                    "elements are not isoparametric; this works only for isoparametric elements." );

    const PropertyDatabase<dim>&   p_ref = sg.Database();
    
    csmp::Index  mob_key = p_ref.StorageKey("total mobility");
    if ( mob_key.place != ELEMENT || mob_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):",
                      "total mobility", "variable must be  a scalar placed on the element." );

    csmp::Index  pressure_key = p_ref.StorageKey("fluid pressure");
    if ( pressure_key.place != NODE || pressure_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):",
                      "fluid pressure", "variable must be a scalar placed on the node." );

    csmp::Index  nsource_key = p_ref.StorageKey("nodal fluid volume source");
    if ( nsource_key.place != NODE  and  nsource_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):",
                      "nodal fluid volume source", "variable must be a scalar placed on the node." );

    csmp::Index  gterm_key = p_ref.StorageKey("gravity term");
    if ( gterm_key.place != ELEMENT  and  gterm_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):",
                      "gravity term", "variable must be a vector placed on the element." );
          
    csmp::Index  cterm_key = p_ref.StorageKey("capillary term");
    if ( cterm_key.place != ELEMENT  and  cterm_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "PressureSolver<dim>::(constructor):",
                      "capillary term", "variable must be a vector placed on the element." );

    this->Add( &total_mobility_ );
    this->Add( &nodal_source_ );
    if ( with_gravity ) this->Add( &gravity_term_ );
    if ( with_capillary )this->Add( &capillary_term_ );
}

#endif /* !CSMP_WITH_SAMG_SOLVER */




template<uint32_t dim, template<uint32_t> class CELL>
PressureSolver<dim, CELL>::~PressureSolver()
{}


template class PressureSolver<1U,Element>;
template class PressureSolver<2U,Element>;
template class PressureSolver<3U,Element>;
template class PressureSolver<1U,Face>;
template class PressureSolver<2U,Face>;
template class PressureSolver<3U,Face>;

} // end csmp
