//
//  Integrator.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 2/8/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "Integrator.h"
#include "LinearAlgebraicSystem.h"
#include "Region.h"
#include "ImplicitTransport.h"

using namespace std;

namespace csmp {

template<size_t dim,template<size_t> class USER>
Integrator<dim,USER>::Integrator( size_t m_x_n, double lower_limit, double upper_limit )
 : solver_(&settings_),
   lower_limit_(lower_limit), upper_limit_(upper_limit)
 {
    // configuring SAMG for first use
    settings_.SetSolverInstance(2);
 }
 



template<size_t dim, template<size_t> class USER>
void Integrator<dim,USER>::ReconfigureSolverForRepeatedUse()
 {
    // minimize output
    settings_.Set_iout1( 0 );
    settings_.Set_iout2( 0 );

    // use relative error solution criterion
    settings_.Set_eps(0.);
    settings_.Set_rel_eps(1.E-10);

    // trigger file dump
    //settings_.Set_idmp( 8 );        // define SAMG command and file output
    //settings_.Set_ioform( "f" );    // define SAMG file output format for reduced file size, idmp > 1 is required
    //settings_.Set_filnam_dump( "SAMG_Transport" ); // set filename for SAMG file output other than default "level", idmp > 1 is required
 }
 
 
 
/**
   assuming that the IntegralEquation and Accumulator have already been applied  to build system and Solver to solve it
*/
template<size_t dim, template<size_t> class USER>
void Integrator<dim,USER>::IntegrateOver( double time_increment )
 {
//    if ( time_increment > 0. ) AssignInitialConditions( subdomain );
//    AssignEssentialConditions( subdomain );
    SolveLinearAlgebraicSystem();
    const bool show_range(false), do_range_check(true); 
    VerifyAndAssignResults( show_range, do_range_check );
    
 } // end IntegrateOver
 
 

/*
template<size_t dim, template<size_t> class USER>
void Integrator<dim,USER>::AssignInitialConditions( const COMPUTATION_DOMAIN<dim>& subdomain )
 {
 } // end AssignInitialConditions




template<size_t dim, template<size_t> class USER>
void Integrator<dim,USER>::AssignEssentialConditions( const COMPUTATION_DOMAIN<dim>& subdomain )
 {
 } // end AssignEssentialConditions
*/



template<size_t dim, template<size_t> class USER>
void Integrator<dim,USER>::SolveLinearAlgebraicSystem()
 {
    if ( User()->Verbose() ) {
         cout <<"\n\nIntegrator::SolveLinearAlgebraicSystem: calling solver instance: ";
         cout << settings_.GetSolverInstance() << endl;
         cout <<"\tSAMG settings:";
         cout <<"\n\t\tiswit  = " << settings_.Get_iswit();
         cout <<"\n\t\titypu  = " << settings_.Get_ifirst();
         cout <<"\n\t\tlevelx = " << settings_.Get_levelx();
         cout << endl;
      }
     solver_.Solve( User()->LinearSystem().LHS, User()->LinearSystem().RHS, User()->LinearSystem().X );
     
 } // end SolveLinearAlgebraicSystem
 
 
 
 
/** 
    Checks the range of 'new saturation oil' against that specified in the property database.
    Where the new values comply, they are used to replace the previous ones stored as 'saturation oil'.
    If not, the deviations are reported and the nearest maximum or minimum permitted values
    of saturation are stored.
*/
template<size_t dim, template<size_t> class USER>
double Integrator<dim,USER>::VerifyAndAssignResults( bool show_range, bool do_range_check )
  {
    double amin(+std::numeric_limits<double>::max());
    double amax(-std::numeric_limits<double>::max());
    double difference_to_last_output(0.);
    size_t   error_counter(0);
    
    const typename vector<Node<dim>*>::iterator  nodes_end(User()->ComputationDomain().NodesEnd());
    for ( typename vector<Node<dim>*>::iterator  nit = User()->ComputationDomain().NodesBegin(); nit != nodes_end; ++nit )
    {
      const VARIABLE_FLAG status((*nit)->Status( User()->key_C0 ));
      if ( status != DIRICH )
        {
          // reading the newly computed saturation values
          const double c1 = User()->LinearSystem().X[ (*nit)->Idx() ];
          amin = std::min( amin, c1 );
          amax = std::max( amax, c1 );
          
          // reading the previous values and calculating the maximum change per node
          const double C0 = (*nit)->Read( User()->key_C0 );

          difference_to_last_output = std::max( difference_to_last_output, fabs(c1 - C0) );
          
          // result checking and assignment
          if ( c1 <= upper_limit_ && c1 >= lower_limit_ ) (*nit)->Store( User()->key_C0, makeScalar(status, c1) );
          else {
            cerr <<"\nExplicitTransport<dim>::VerifyAndAssignResults: ";
            cerr <<"value: "<< c1 <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
            if ( c1 > upper_limit_ ) (*nit)->Store( User()->key_C0, makeScalar( status, upper_limit_ ) );
            else if ( c1 < lower_limit_ ) (*nit)->Store( User()->key_C0, makeScalar( status, lower_limit_ ) );
            error_counter++;
          }
        }
      }
    
    if ( do_range_check ) {
      if ( error_counter > User()->ComputationDomain().Nodes() )
        throw out_of_range("Integrator<dim>::VerifyAndAssignResults: Advected variable out of range");
    }
    
    if ( show_range ) {
      cout <<"\n\nExplicitTransport<"<< dim;
      cout <<">::VerifyAndAssignResults: Variable range after advection: ";
      cout << amin <<" to "<< amax << endl << endl;
    }
    
    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );
    
 } // end VerifyAndAssignResults

 


template class Integrator<3U,ImplicitTransport>;

} // end csmp
