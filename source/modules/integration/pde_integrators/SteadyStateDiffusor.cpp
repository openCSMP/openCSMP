#include "SteadyStateDiffusor.h"
#include "Exception.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "MeshManagementUtilities.h"
#include "PL_Utilities.h"

using namespace std;

namespace csmp {

	template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
	SteadyStateDiffusor<dim, COMPUTATION_DOMAIN>::SteadyStateDiffusor( Model<dim>& sg,
                                                                     const char* diffusivity,
                                                                     const char* diffusing_variable,
                                                                     const char* spatial_source_variable,
                                                                     bool LumpedRHS )
		:
#ifdef CSMP_WITH_SAMG_SOLVER
		settings_(),
		PDE_Integrator<dim, COMPUTATION_DOMAIN>(new SAMG_Solver(&settings_)),
#else
		/// add extra functionality for alternative solver if needed
		PDE_Integrator<dim, COMPUTATION_DOMAIN>(new CSMP_DEFAULT_LINEAR_SOLVER()),
#endif
		conductance_(sg.Database(), diffusivity, diffusing_variable, diffusing_variable),
		source_(new NumIntegral_NT_op_N_dV<dim, ComputationCell>(sg.Database(), spatial_source_variable, diffusing_variable)),
		nodal_source_(0),
		gravity_(0),
		grad_multiplier_(1.),
		dep_var_name_(diffusing_variable),
		firstCall_(true)
	{
		if (!isoparametricElementMesh(sg))
			throw csmp::Exception(FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
				"elements are not isoparametric; use other Algorithm");

		const PropertyDatabase<dim>&  p_ref = sg.Database();

		csmp::Index  conductivity_key = p_ref.StorageKey(diffusivity);
		if (conductivity_key.place != ELEMENT and conductivity_key.place != ELEMENT_INTEGRATION_POINT
			and conductivity_key.place != FACE)
			throw csmp::Exception(FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
				diffusivity, "variable must be placed on element, face or element integration point");

		csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
		if (diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR)
			throw csmp::Exception(FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
				diffusing_variable, "variable must be a scalar placed on the node");

		csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
		if (source_variable_key.type != SCALAR)
			throw csmp::Exception(FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
				spatial_source_variable, "variable must be a scalar");

		if (source_variable_key.place != ELEMENT and source_variable_key.place != ELEMENT_INTEGRATION_POINT
			and source_variable_key.place != FACE)
			throw csmp::Exception(FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
				spatial_source_variable, "variable must be placed on element, face or their integration points");

		// initializing the algorithm  
		this->Add(&conductance_);
		if (LumpedRHS) source_->LumpedFormulation();
		this->Add(source_);

	} // end constructor

template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::SteadyStateDiffusor( Model<dim>& sg,
                                                                  const char* diffusivity,
                                                                  const char* diffusing_variable,
                                                                  const char* spatial_source_variable )
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
  settings_(),
  PDE_Integrator<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
   /// add extra functionality for alternative solver if needed
   PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_( new NumIntegral_NT_op_N_dV<dim,ComputationCell>( sg.Database(), spatial_source_variable, diffusing_variable ) ),
   nodal_source_(0),
   gravity_(0),
   grad_multiplier_(1.),
   dep_var_name_(diffusing_variable),
   firstCall_(true)
 {
    if ( !isoparametricElementMesh( sg ) )
        throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
                                            "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&  p_ref = sg.Database();
    
    csmp::Index  conductivity_key = p_ref.StorageKey(diffusivity);
    if ( conductivity_key.place != ELEMENT and conductivity_key.place != ELEMENT_INTEGRATION_POINT 
         and conductivity_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusivity, "variable must be placed on element, face or element integration point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
    if ( source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be a scalar" ); 

    if ( source_variable_key.place != ELEMENT and source_variable_key.place != ELEMENT_INTEGRATION_POINT
         and source_variable_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be placed on element, face or their integration points" ); 

    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
     
} // end constructor






template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::SteadyStateDiffusor( Model<dim>& sg, 
                                                                                const char* diffusivity,
                                                                                const char* diffusing_variable,
                                                                                const char* spatial_source_variable,  
                                                                                const char* point_source_variable ) 
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
     PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), spatial_source_variable, diffusing_variable) ),
   nodal_source_(new PointSource_rhsop<dim,ComputationCell>(sg.Database(), point_source_variable, diffusing_variable) ),
   gravity_(0),
   grad_multiplier_(1.),
   dep_var_name_(diffusing_variable),
   firstCall_(true)
 {
    if ( !isoparametricElementMesh( sg ) )
        throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
                                            "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&  p_ref = sg.Database();
    
    csmp::Index  conductivity_key = p_ref.StorageKey(diffusivity);
    if ( conductivity_key.place != ELEMENT  and  conductivity_key.place != ELEMENT_INTEGRATION_POINT 
         and conductivity_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                      diffusivity, "variable must be placed on element, face or element integration point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                      diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
    if ( source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be a scalar" ); 

    if ( source_variable_key.place != ELEMENT and source_variable_key.place != ELEMENT_INTEGRATION_POINT
         and source_variable_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be placed on element, face or their integration points" ); 

    csmp::Index  nsource_variable_key = p_ref.StorageKey(point_source_variable);
    if ( nsource_variable_key.place != NODE  and  nsource_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                      point_source_variable, "variable must be a scalar placed on the node" ); 
    
    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
    this->Add(  nodal_source_ );
     
} // end constructor





template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::SteadyStateDiffusor( Model<dim>& sg,
                                                                                const char* diffusivity,
                                                                                const char* diffusing_variable,
                                                                                const char* gradient_variable, 
                                                                                double gradient_multiplier )
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_(0),
   nodal_source_(0),
   gravity_(new NumIntegral_dNT_op_dV<dim,ComputationCell>(sg.Database(), gradient_variable, diffusing_variable) ),
   grad_multiplier_(gradient_multiplier),
   dep_var_name_(diffusing_variable),
   firstCall_(true)
 {
    if ( !isoparametricElementMesh( sg ) )
        throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
                                            "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&  p_ref = sg.Database();
    
    csmp::Index  diffusivity_key = p_ref.StorageKey(diffusivity);
    if ( diffusivity_key.place != ELEMENT  and diffusivity_key.place != ELEMENT_INTEGRATION_POINT 
         and diffusivity_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusivity, "variable must be placed on element, face or element integration point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  gradient_variable_key = p_ref.StorageKey(gradient_variable);
    if ( gradient_variable_key.place != ELEMENT and gradient_variable_key.place != ELEMENT_INTEGRATION_POINT
         and gradient_variable_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             gradient_variable, 
                            "variable must be a scalar placed on element, face or element integration point" ); 

    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  gravity_ );
     
} // end constructor (gravity without sources)




template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::SteadyStateDiffusor( Model<dim>& sg,
                                                                                const char* diffusivity,
                                                                                const char* diffusing_variable,
                                                                                const char* spatial_source_variable,                    
                                                                                const char* gradient_variable, 
                                                                                double gradient_multiplier )
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), spatial_source_variable, diffusing_variable) ),
   nodal_source_(0),
   gravity_(new NumIntegral_dNT_op_dV<dim,ComputationCell>(sg.Database(), gradient_variable, diffusing_variable) ),
   grad_multiplier_(gradient_multiplier),
   dep_var_name_(diffusing_variable),
   firstCall_(true)
 {
    if ( !isoparametricElementMesh( sg ) )
        throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
                                            "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&  p_ref = sg.Database();
    
    csmp::Index  diffusivity_key = p_ref.StorageKey(diffusivity);
    if ( diffusivity_key.place != ELEMENT  and  diffusivity_key.place != ELEMENT_INTEGRATION_POINT 
         and diffusivity_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusivity, "variable must be placed on element, face or element integration point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
    if ( source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be a scalar" ); 

    if ( source_variable_key.place != ELEMENT and source_variable_key.place != ELEMENT_INTEGRATION_POINT
         and source_variable_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be placed on element, face or their integration points" ); 

    csmp::Index  gradient_variable_key = p_ref.StorageKey(gradient_variable);
    if ( (gradient_variable_key.place != ELEMENT and gradient_variable_key.place != FACE and gradient_variable_key.place != ELEMENT_INTEGRATION_POINT)
          && gradient_variable_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             gradient_variable, 
                            "variable must be a vector placed on element, face or integration point" ); 

    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
    this->Add(  gravity_ );
     
} // end constructor (gravity with sources)






// discerning lhs and rhs diffusivity terms
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::SteadyStateDiffusor( Model<dim>& sg,
                                                                                const char* lhs_diffusivity,
                                                                                const char* rhs_diffusivity,
                                                                                const char* diffusing_variable,
                                                                                const char* spatial_source_variable,                    
                                                                                const char* gradient_variable, 
                                                                                double gradient_multiplier )
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), lhs_diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), spatial_source_variable, diffusing_variable) ),
   nodal_source_(0),
   gravity_(new NumIntegral_dNT_op_dV<dim,ComputationCell>(sg.Database(), gradient_variable, diffusing_variable) ),
   grad_multiplier_(gradient_multiplier),
   dep_var_name_(diffusing_variable),
   firstCall_(true)
 {
    if ( !isoparametricElementMesh( sg ) )
        throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):",
                                            "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&  p_ref = sg.Database();
    
    csmp::Index  lhs_diffusivity_key = p_ref.StorageKey(lhs_diffusivity);
    if ( lhs_diffusivity_key.place != ELEMENT  and  lhs_diffusivity_key.place != ELEMENT_INTEGRATION_POINT
         and lhs_diffusivity_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             lhs_diffusivity, "variable must be placed on element, face or element integration point" ); 

    csmp::Index  rhs_diffusivity_key = p_ref.StorageKey(rhs_diffusivity);
    if ( rhs_diffusivity_key.place != ELEMENT  and  rhs_diffusivity_key.place != ELEMENT_INTEGRATION_POINT
         and rhs_diffusivity_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             rhs_diffusivity, "variable must be placed on element, face or element integration point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
    if ( source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be a scalar" ); 

    if ( source_variable_key.place != ELEMENT and source_variable_key.place != ELEMENT_INTEGRATION_POINT
         and source_variable_key.place != FACE )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             spatial_source_variable, "variable must be placed on element, face or their integration points" ); 

    csmp::Index  gradient_variable_key = p_ref.StorageKey(gradient_variable);
    if ( (gradient_variable_key.place != ELEMENT and gradient_variable_key.place != FACE and
          gradient_variable_key.place != ELEMENT_INTEGRATION_POINT) && gradient_variable_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "SteadyStateDiffusor<dim>::(constructor):", 
                             gradient_variable, 
                            "variable must be a vector placed on element, face or element integration point" ); 

    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
    this->Add(  gravity_ );
     
} // end constructor (gravity with sources)






template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::~SteadyStateDiffusor()
 {
    delete source_;
    delete nodal_source_;
    delete gravity_;
 } 





/**
*/




#ifdef CSMP_WITH_SAMG_SOLVER
/// From the SAMG solver profile brought here because many users have no clue that profile exists
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::Adjust_SAMG_ForSubsequentSolves()
 {
    assert ( firstCall_ );

    #ifdef SAMG_MULTIPLE_INSTANCES
        if ( firstCall_ ) {
           /// Use solution of previous timestep as an initial guess
          settings_.Set_itypu(0);
          cout << "\n\n*** SteadyStateDiffusor::Adjust_SAMG_ForSubsequentSolves: first call: Set SAMG input parameter itypu = 0";
          cout.flush();

          /// Reuse SAMG solver setup, this step is only necessary if initial setting was iswit(4)
          #ifdef NO_PRIMARY_SOLVER_CONTROL
              settings_.Set_iswit(3); // re-use solver setup from previous timestep
              cout << ", iswit = 3,";
              cout.flush();
          #endif

          #ifdef RELATIVE_CONVERGENCE
              /// Relative convergence is used as stopping criterion "res <= eps.res0"
              /// (res0 = starting residual) - for IMPES without SAMG Multiple Instances only
              settings_.Get_rel_eps() );
              cout << " and relative solution criterion eps = " << settings_.Get_eps() << " after first solver call' ***\n\n";
              cout.flush();
          #else
              /// Absolute convergence is used as stopping criterion "res <= eps" (res0 = starting residual)
              settings_.Set_eps( 1.E-14 );
              cout << " and absolute solution criterion eps = " << -settings_.Get_eps() << " after first solver call' ***\n\n";
              cout.flush();
          #endif
        }
    #endif

    #ifndef SAMG_MULTIPLE_INSTANCES
         if ( firstCall_ ) {
         #ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES

              /// Use solution of previous timestep as an initial guess - for IMPES without SAMG Multiple Instances only
              settings_.Set_itypu(0);
              cout << "\n\n*** SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::Adjust_SAMG_ForSubsequentSolves: ";
              cout <<" first call, setting SAMG input parameter itypu = 0' ***\n";
              cout.flush();

              /// Reuse SAMG solver setup - for IMPES without SAMG Multiple Instances only - this step is only necessary if initial setting was iswit(4)
              #ifdef NO_PRIMARY_SOLVER_CONTROL
                  settings_.Set_iswit(3); // re-use solver setup from previous timestep
                  cout << "\n\n*** SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::Adjust_SAMG_ForSubsequentSolves: first call: ";
                  cout <<" setting SAMG input parameter iswit = 3' ***\n";
                  cout.flush();
              #endif

          #endif

          #ifdef RELATIVE_CONVERGENCE
              /// Relative convergence is used as stopping criterion "res <= eps.res0" (res0 = starting residual) - for IMPES without SAMG Multiple Instances only
              settings_.Set_eps( this->Solver().GetSolverSettings().Get_rel_eps() );
              cout << "\n\n*** SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::Adjust_SAMG_ForSubsequentSolves: first call: ";
              cout <<" setting SAMG relative solution criterion eps = " << settings_.Get_eps() << "' ***\n\n";
              cout.flush();
          #else
              /// Absolute convergence is used as stopping criterion "res <= eps" (res0 = starting residual)
              settings_.Set_eps( 1.E-14 );
              cout << "\n\n*** SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::Adjust_SAMG_ForSubsequentSolves first call: ";
              cout <<" setting SAMG absolute solution criterion eps = " << -settings_.Get_eps() << "' ***\n\n";
              cout.flush();
          #endif
          }
    #endif

 } // end Adjust_SAMG_ForSubsequentSolves
#else
   /// add extra functionality for alternative solver if needed
#endif


/// Apply PDE integrator to entire model
/* does not work with the explicit template instantiations
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::ComputeSteadyState( Model<dim>& sg, bool verbose )
 {
#ifdef CSMP_WITH_SAMG_SOLVER
    if ( verbose ) {
          cout <<" \n\n\nSteadyStateDiffusor<"<< dim <<">: ComputeSteadyState: ";
          cout <<" Computing '"<< dep_var_name_ <<"'"<< endl;
          std::cout << "\n*** SteadyStateDiffusor::ComputeSteadyState 'Received iswit = " << settings_.Get_iswit() << " and";
          std::cout << " itypu = " << settings_.Get_ifirst() << "' ***\n";
      }
   
    if ( firstCall_ ) settings_.Set_itypu(1); // initial guess is produced by SAMG
    this->IntegrateOver( sg.Region("Model"), verbose );

    if ( firstCall_ ) {
         Adjust_SAMG_ForSubsequentSolves();
         firstCall_=false;
      }
#else
   /// add extra functionality for alternative solver if needed
    sg.Apply( *this, verbose );
#endif

 } // ComputeSteadyState
*/





/// Apply PDE integrator to Region / Boundary / SplitBoundary
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::ComputeSteadyState( COMPUTATION_DOMAIN<dim>& sd, bool verbose )
 {
#ifdef CSMP_WITH_SAMG_SOLVER
    if ( verbose ) {
          cout <<" \n\n\nSteadyStateDiffusor<"<< dim <<">: ComputeSteadyState: ";
          cout <<" Computing '"<< dep_var_name_ <<"'"<< endl;
          cout << "\n*** SteadyStateDiffusor::ComputeSteadyState 'Received iswit = " << settings_.Get_iswit() << " and";
          cout << " itypu = " << settings_.Get_ifirst() << "' ***\n";
          cout.flush();
      }

    if ( firstCall_ ) settings_.Set_itypu(1); // initial guess is produced by SAMG
    this->IntegrateOver(sd);

/*
    if ( firstCall_ ) {
         Adjust_SAMG_ForSubsequentSolves();
         firstCall_=false;
      }
*/
#else
   /// add extra functionality for alternative solver if needed
   this->IntegrateOver(sd);
#endif

 } // ComputeSteadyState




/// empirically established best settings for fluid pressure diffusion
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::AdjustSolverSettings()
 {
#ifdef CSMP_WITH_SAMG_SOLVER
    cout <<"\n\n*** SteadyStateDiffusor::AdjustSolverSettings() to best settings for 'fluid pressure diffusion' ***\n\n";
    cout.flush();

   /// Basic solver settings for hybrid element meshes
    settings_.Set_itypu(1);     // create initial guess yourself
    settings_.Set_ncgtyp(5);
    settings_.Set_nxtyp(1);     // 1=ILU, 0=Gauss Seidel
    settings_.Set_ndefault(40);
    settings_.Set_iout1(0);
    settings_.Set_iout2(0);

    /// SAMG output to file
    //settings_.Set_idmp( 8 );          // define SAMG command and file output
    //settings_.Set_ioform( "f" );      // define SAMG file output format for reduced file size, idmp > 1 is required
    //settings_.Set_filnam_dump( "SAMG_Diffusion" ); // set filename for SAMG file output other than default "level", idmp > 1 is required
#else
   /// add extra functionality for alternative solver if needed
#endif

 } // end AdjustSolverSettings



#ifdef CSMP_WITH_SAMG_SOLVER
/// returns a reference to the current settings of the SAMG Solver used by the SteadyStateDiffusor
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
SAMG_Settings& SteadyStateDiffusor<dim,COMPUTATION_DOMAIN>::GetSolverSettings()
 { return settings_; }
#else
   /// add extra functionality for alternative solver if needed
#endif

template class SteadyStateDiffusor<1U,Region>;
template class SteadyStateDiffusor<2U,Region>;
template class SteadyStateDiffusor<3U,Region>;

template class SteadyStateDiffusor<1U,Boundary>;
template class SteadyStateDiffusor<2U,Boundary>;
template class SteadyStateDiffusor<3U,Boundary>;

} // end csmp
