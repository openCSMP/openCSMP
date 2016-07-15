#include "StressesAndStrains2.h"
#include "mechanics.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#include "Element.h"

using namespace std;

namespace csmp {

/**
 
Constructs the post-processing math operator for 2D calculations. Thus,
the user must specify whether he/she is interested in plane stress or
plane strain computations (this must be exactly the same setting as for
the actual computation), and they have the option of computing the 
additional variables (character strings):  

@code
name 				type 		placement
----				----		---------
strain				tensor		cpoint
'strain1' 			vector		cpoint	(principal strain axes)
'strain2'			vector		cpoint
('strain3')			vector		cpoint
'dilatation'		scalar		cpoint
'principal strain' 	vector		cpoint	(magnitudes of elongation in the principal directions)
'stress'			tensor		cpoint
'sigma1'  			vector		cpoint	(principal stress axes)
'sigma2'			vector		cpoint
('sigma3')			vector		cpoint
'mean stress'  		scalar 		cpoint	(= (sigma1 + sigma2 +sigma3 / 3)
'principal stress'  vector		cpoint	(magnitudes of the principal stresses)
 @endcode

Clearly, these variables must be present in the variable database file
for the operator to work.  

The stiffness matrices needed for the stress computation are calculated
from the user-specified variables for 'Young's modulus' and
'Poisson's ratio'. 

@section arguments Input Arguments 

Three character strings for the input variables of the calculation
 'Young's modulus', 'Poisson's ratio', and 'displacement', a boolean
 variable to toggle between plane stress and plane strain and another
 boolean variable specifying whether principal strains and stresses
 shall be computed.  

@section messages Messages 

The names of the output variables are hard-wired in this operator. Thus,
checks are performed whether these variables actually exist in the
variable database and errors are reported if they don't or if they
have the wrong placement or type.  
*/

StressesAndStrains<2U>::StressesAndStrains( const Model<2U> & sg, 
                                               const char* oper,  // Young's modulus
                                               const char* basic, // Poisson's ratio
                                               const char* test,  // displacement
                                               bool  planestrain,
                                               bool  principal_vectors )
  : MathOperatorLHS<2U>(sg.Database(),
                        oper,basic,test),    
    // getting the necessary csmp::Index keys
    strain_key_(sg.Database().StorageKey("strain")),
    stress_key_(sg.Database().StorageKey("stress")),
    DISPL_(3*2,1), 
    STIFF_(1,DenseMatrix<DM_MIN>(3,3)), EGP_(3,1), SGP_(3,1),
    STRAIN_(3,3), 
    STRESS_(3,3),
    IPSTRAIN_(3*3), // 3 strain components
    IPSTRESS_(3*3), // 3 strain components
    NSTRAIN_(3*3),
    NSTRESS_(3*3),
    PR_(2,2), EIG_(2,2), EVAL_(2),
    eps_(3), sigma_(3), sum_(3),
    components_(3),
    plane_strain_(planestrain),
    principal_e_and_sigma_(principal_vectors),
    verbose_(false)
    
 {
    MathOperatorLHS<2U>::Name("StressesAndStrains", oper, basic, test );

    // testing the Operands 
    if ( (strain_key_.place != NODE and strain_key_.place != ELEMENT_INTEGRATION_POINT) || strain_key_.type != TENSOR ) 
      throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                     "'strain' must be a tensor property placed on the constraint points or element." );

    if ( (stress_key_.place != NODE and stress_key_.place != ELEMENT_INTEGRATION_POINT) || stress_key_.type != TENSOR ) 
      throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                     "'stress' must be a tensor property placed on the constraint points or element." );

    // displacement
    if ( MathOperatorLHS<2U>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<2U>::TestOperandType() != VECTOR )
    throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                    test, "Operand 'displacement' must be a vector property placed on the nodes." );

    // Young's modulus & Poisson's ratio
    if ( (MathOperatorLHS<2U>::MaterialOperandPlacement() != ELEMENT &&
          MathOperatorLHS<2U>::MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT) ||
         MathOperatorLHS<2U>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                    oper, "Operand 'Young's modulus' must be a scalar property placed on the constraint points or element." );

    // Poisson's ratio
    if ( (MathOperatorLHS<2U>::BasicOperandPlacement() != ELEMENT &&
          MathOperatorLHS<2U>::BasicOperandPlacement() != ELEMENT_INTEGRATION_POINT) ||
         MathOperatorLHS<2U>::BasicOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                    basic, "Basic Operand 'Poisson's ratio' must be a scalar property placed on the constraint points or element." );
                   
    if ( strain_key_.place != stress_key_.place or  
         MathOperatorLHS<2U>::MaterialOperandPlacement() != MathOperatorLHS<2U>::BasicOperandPlacement() )
      throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                     "'Young's modulus' and 'Poisson's ratio', and 'strain' and 'stress' must have the same placement." );

    // principal stresses and strains
    if ( principal_vectors ) 
      {
         strain1_key_ = sg.Database().StorageKey("strain1");
         strain2_key_ = sg.Database().StorageKey("strain2");
         sigma1_key_  = sg.Database().StorageKey("sigma1");
         sigma2_key_  = sg.Database().StorageKey("sigma2");
         means_key_   = sg.Database().StorageKey("mean stress");
         dilat_key_   = sg.Database().StorageKey("dilatation");
         pstrain_key_ = sg.Database().StorageKey("principal strain");
         pstress_key_ = sg.Database().StorageKey("principal stress");
        
         if ( (strain1_key_.place != NODE and strain1_key_.place != ELEMENT_INTEGRATION_POINT) || strain1_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'strain1' must be a vector property placed on the constraint points or nodes." );

         if ( (strain2_key_.place != NODE and strain2_key_.place != ELEMENT_INTEGRATION_POINT) || strain2_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'strain2' must be a vector property placed on the constraint points or nodes." );

         if ( (sigma1_key_.place != NODE and sigma1_key_.place != ELEMENT_INTEGRATION_POINT) || sigma1_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'sigma1' principal stress must be a vector property placed on the constraint points or nodes." );

         if ( (sigma2_key_.place != NODE and sigma2_key_.place != ELEMENT_INTEGRATION_POINT) || sigma2_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'sigma2' principal stress must be a vector property placed on the constraint points or nodes." );
         
         if ( (pstrain_key_.place != NODE and pstrain_key_.place != ELEMENT_INTEGRATION_POINT) || pstrain_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'principal strain' must be a vector property placed on the constraint points or nodes." );
                   
         if ( (pstress_key_.place != NODE and pstress_key_.place != ELEMENT_INTEGRATION_POINT) || pstress_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'principal stress' must be a vector property placed on the constraint points or nodes." );

         if ( means_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'mean stress' must be a scalar property placed on the constraint points or nodes." );

         if ( dilat_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
                   "Operand 'dilatation' must be a scalar property placed on the constraint points or nodes." );

	     if ( strain1_key_.place != strain2_key_.place or 
	          strain2_key_.place != strain_key_.place ) 
	       throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
	                      "strains and principal strains must have the same placement." );

	     if ( sigma1_key_.place != sigma2_key_.place or 
	          sigma2_key_.place != stress_key_.place ) 
	       throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
	                      "stresses and principal stresses must have the same placement." );

	     if ( stress_key_.place  != means_key_.place or 
	          means_key_.place   != dilat_key_.place or 
	          dilat_key_.place   != strain1_key_.place or
	          strain1_key_.place != sigma1_key_.place ) 
	       throw csmp::Exception( ERROR, "StressesAndStrains<2D>::(constructor)", 
	                      "'dilation', mean stress', 'strain', and 'sigma' must have the same placement." );
      }               
                   
    // setting how many times the PDE operator shall be applied during post-processing
    // 2 cycles are required if nodal averaging of the results is wanted
    if ( stress_key_.place == ELEMENT_INTEGRATION_POINT ) MathOperatorLHS<2U>::ApplicationCycles(1);
    else {
         MathOperatorLHS<2U>::ApplicationCycles(2);
         node_output_.resize(sg.Region("Model").Nodes(),false);
         temp_strains_.resize(sg.Region("Model").Nodes());
         temp_stresses_.resize(sg.Region("Model").Nodes());
      }
    
 } // end 2D constructor












/**
 
Allows to toggle between plane stress and plane strain computations in
two-dimensional models. Clearly, this method only has an effect in 2D
calculations and the chosen computation must match the preceding
Algorithm.  

@section arguments Input Arguments 

A boolean variable specifying whether the material property matrix
used shall be initialized for plane strain or plane stress.   

@section messages Messages 

If the method is called in a 3D calculation, the user is warned that
it will have no effect.  

*/

void StressesAndStrains<2U>::PlaneStress( bool yes_no ) 
 { 
    plane_strain_ = yes_no; 
 }











/**
 
Allows to toggle the addional computation of the variables:

@code
name 				type 		placement
----				----		---------
strain				tensor		node
'strain1' 			vector		node	(principal strain axes)
'strain2'			vector		node
('strain3')			vector		node
'dilatation'		scalar		node
'principal strain' 	vector		node	(magnitudes of elongation in the principal directions)
'stress'			tensor		node
'sigma1'  			vector		node	(principal stress axes)
'sigma2'			vector		node
('sigma3')			vector		node
'mean stress'  		scalar 		node	(= (sigma1 + sigma2 +sigma3 / 3)
'principal stress'  vector		node	(magnitudes of the principal stresses)
 @endcode

on or off. Clearly, these variables must be present in the variable 
database file for the method to work.  

@section arguments Input Arguments 

The boolean variable which specifies whether one wants the additional
computations or not. The default is 'true'.  
*/

void StressesAndStrains<2U>::PrincipalStrainsAndStresses( bool yes_no )
 {
    principal_e_and_sigma_ = yes_no;
 }



/**

Reads the displacements (specified as test-function Operand) of each node
for further processing by the ComputeContribution() method.  

Then Youngs modulus and the Poisson's ratio are read and the stiffness
matrix is build   

@section arguments Input Arguments 

The element from which accumulation into the global matrix takes place
and the time-increment over which the force will be applied.  

@section implementation Implementation

Since the forces are nodal properties the GetOperands method of the 
base class MathOperatorRHS must be overloaded, such that the nodal Operand
can be retrieved. 
*/

void StressesAndStrains<2U>::GetOperands( Element<2U>& e )
{
    vector<VectorVariable<2U> >  NVAR(e.Nodes());
    size_t                       k;
    VectorVariable<2U>           avg;
                    
    if ( MathOperatorLHS<2U>::ApplicationCycle() == 1 ) {
         // 1. For strain computation at integration points from nodal displacements
         // ------------------------------------------------------------------------
         // nodal displacements
         e.NodePropertyVector( MathOperatorLHS<2U>::TestOperandKey(), NVAR );
    
         // sticking displacements into a matrix of size nodes*dof x 1
         // and recording the elements spatial translation as the average of
         // its displacements
         DISPL_.Resize( e.Nodes()*2U, 1 );
         k = 0, avg = 0.;
         for ( size_t i=0; i<e.Nodes(); i++ ) {
              for ( size_t j=0; j<2U; j++ ) DISPL_(k++,0) = NVAR[i].Value(j);
              avg         += NVAR[i];
           }
         avg /= static_cast<double64>(e.Nodes());

         if ( verbose_ ) {     
              cout <<"\nStressesAndStrains<"<< 2U;
              cout <<">::GetOperands: Nodal displacements, element: "<< e.Idx() << endl;
              DISPL_.Out();  
           }      
        
        // 2. For stress computation from strains at integration points
        // ------------------------------------------------------------
        // getting Young's modulus and Poisson's ratio 
        if ( MathOperatorLHS<2U>::MaterialOperandPlacement() == ELEMENT ) {
             youngs_.resize(1U);
             e.Read( MathOperatorLHS<2U>::MaterialOperandKey(), youngs_[0] );
             pratio_.resize(1U);
             e.Read( MathOperatorLHS<2U>::BasicOperandKey(), pratio_[0] );
          }
        // it must be a constraint-point property  
        else {
             e.IntegrationPointPropertyVector( MathOperatorLHS<2U>::MaterialOperandKey(), youngs_ );
             e.IntegrationPointPropertyVector( MathOperatorLHS<2U>::BasicOperandKey(), pratio_ );
          }

        // computing the stiffness matrix for the respective element
        if ( plane_strain_ ) planeStrainMatrix( youngs_, pratio_, STIFF_ );
        else                 planeStressMatrix( youngs_, pratio_, STIFF_ );
     } 
    
} // end GetOperands







/**
 
1. From the nodal displacements, compute the strains at the integration
   points (=optimal locations, Barlow, 1977), using the interpolation
   function derivative matrices at these points.  

2. Linearly extrapolate the strains back from the integration points 
   to the nodes. In this step, bilinear extrapolation is carried out
   by the specific finite element used. Only linear variations in strain
   over the element are supported at this stage (C1 elements). 
   The nodal strains are averaged among the adjacent elements since 
   they diverge at these locations from one-another. This is done
   in a second application cycle of this postprocessing 
   mathoperator.  

@section arguments Input Arguments 

A reference to the Element for which the post-processing is done.  
*/

void StressesAndStrains<2U>::ComputeContribution( Element<2U>& e )
{
   // 1. Computing STRAINS and STRESSES at the integration points from the nodal 
   //    displacements and the interpolation function derivative matrices
   //    at the integration points
   // -------------------------------------------------------------------
   if ( MathOperatorLHS<2U>::ApplicationCycle() == 1 ) {
        if ( verbose_ )
          cout <<"\n\nStressesAndStrains<2U>::ComputeContribution: Element: "<< e.Idx() << endl; 

        // strain - stress components of symmetric matrix at integration points 
        IPSTRAIN_.resize( e.FE()->IntegrationPoints() * components_ );
        IPSTRESS_.resize( e.FE()->IntegrationPoints() * components_ );
        
        // 1.1 Computing the strains at the integration points
        // ---------------------------------------------------
        for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) {
             //  getting DN matrices at the node points
            (e).dN_AtIntegrationPoint( EGP_, i, 2U );
             
             // compute strain e_i = [B]_i{d}_i 3x12 * 12x1 -> 3x1 (in 2D)
             EGP_ *= DISPL_;

             // compute {sigma} = [E]([B]{d})
             if ( MathOperatorLHS<2U>::MaterialOperandPlacement() == ELEMENT )
               SGP_ = STIFF_[0];
             else // if a stiffness matrix is defined at each constraint point
               SGP_ = STIFF_[i];
             SGP_ *= EGP_;

             if ( verbose_ ) { 
                  cout <<"\nstrain at integration point: "<< i << endl; 
                  EGP_.Out();
                  cout <<"\nstress at integration point: "<< i << endl; 
                  SGP_.Out();
               }       

             // inserting strains and stresses sequentially into temporary 
             // STL vectors
             for ( size_t k=0; k<components_; k++ ) {
                  IPSTRAIN_[ i * components_ + k ] = EGP_(k,0);
                  IPSTRESS_[ i * components_ + k ] = SGP_(k,0);
               }
          }

        // 1.2 Extrapolating stresses and strains to the nodes
        //     storing them in temporary vectors if later node
        //     averaging is requested.
        // --------------------------------------------------------
        if ( MathOperatorLHS<2U>::ApplicationCycles() == 2 ) {
             // nodal components
             NSTRAIN_.resize( e.Nodes() * components_ );
             NSTRESS_.resize( e.Nodes() * components_ );
             e.ExtrapolateIntegrationPointVariableToNodes( components_, IPSTRAIN_, NSTRAIN_ );
             e.ExtrapolateIntegrationPointVariableToNodes( components_, IPSTRESS_, NSTRESS_ );
     
             // extracting the nodal strain and stress components
             for ( size_t i=0; i<e.Nodes(); i++ ) {
                  for ( size_t k=0; k<components_; k++ ) {
                       eps_[k]   = NSTRAIN_[ i * components_ + k ];
                       sigma_[k] = NSTRESS_[ i * components_ + k ];
                    }
                  // storing strain and stress components in temporary vectors
                  temp_strains_[ e.N(i)->Idx() ].push_back( eps_ );
                  temp_stresses_[ e.N(i)->Idx() ].push_back( sigma_ );
              }
          }
         
     }  // end of first application cycle  


   // 2. During the second visitation the nodal strains and stresses computed by each element
   //    are stored in lists of vectors with their components. The strain and stress
   //    values of the different elements stored in the nodal lists are now
   //    averaged and stored in a vector for output.
   //    -------------------------------------------
   if ( MathOperatorLHS<2U>::ApplicationCycle() == 2 ) {
        STRAIN_.Resize(components_,e.Nodes());
        STRESS_.Resize(components_,e.Nodes());

        for ( size_t i=0; i<e.Nodes(); i++ ) {
          // duplicate calculations are avoided via the boolean vector
          if ( !node_output_[ e.N(i)->Idx() ] )
            {
               for ( size_t j=0U; j<components_; j++ ) 
                 {  
                    // averaging strain components
                    double64 sum(0.);
                    for ( deque<vector<double64> >::const_iterator 
                          lit=temp_strains_[ e.N(i)->Idx() ].begin();
                          lit!=temp_strains_[ e.N(i)->Idx() ].end(); lit++ ) sum += (*lit)[j];
                    sum /= static_cast<double64>(temp_strains_[ e.N(i)->Idx() ].size());
                    STRAIN_(j,i) = sum;

                    // averaging stress components
                    sum = 0.;
                    for ( deque<vector<double64> >::const_iterator
                          lit=temp_stresses_[ e.N(i)->Idx() ].begin();
                          lit!=temp_stresses_[ e.N(i)->Idx() ].end(); lit++ ) sum += (*lit)[j];
                    sum /= static_cast<double64>(temp_stresses_[ e.N(i)->Idx() ].size());
                    STRESS_(j,i) = sum;
                 }
             }

          // now the vector list is no longer needed and therefore erased for the 
          // next application of the post-processing operator
          temp_strains_[ e.N(i)->Idx() ].erase( temp_strains_[ e.N(i)->Idx() ].begin(),
                                                temp_strains_[ e.N(i)->Idx() ].end() );

          temp_stresses_[ e.N(i)->Idx() ].erase( temp_stresses_[ e.N(i)->Idx() ].begin(),
                                                 temp_stresses_[ e.N(i)->Idx() ].end() );
       }

   } // end application cycle 2

} // end ComputeContribution








/** Outputs the computed strains to the Model<2U>  storage checking the range
of each variable using the ranges associated with each of the Operands. x

If principal stresses and strains were requested this method also does 
the necessary computations before outputting these variables.  

@section arguments Input Arguments 

A reference to the property memory manager and the Element for which the 
variables are output.  
*/

void StressesAndStrains<2U>::WriteOperands( Element<2U>& e )
 {
    // if a single stage computation is desired, excluding extrapolations to the nodes
    if ( stress_key_.place == ELEMENT_INTEGRATION_POINT )
      for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ )
        {
            // 1. assigning the strain & stress values
            // ---------------------------------------
            convertTo( IPSTRAIN_, i, ts_ );
            e.Store( i, strain_key_, ts_ ); // O.K.
            convertTo( IPSTRESS_, i, ts_ );
            e.Store( i, stress_key_, ts_ ); // O.K.

            // 2. Computing principal strain and stress axis if requested
            // ----------------------------------------------------------
            if ( principal_e_and_sigma_ )
              {
                 // principal stresses
                 // ------------------
                 ts_.Eigen( evals_, evecs_, true );
                 
                 e.Store( i, sigma1_key_, evecs_.Row(0) );
                 e.Store( i, sigma2_key_, evecs_.Row(1) );
                 
                 // principal stresses = stress magnitudes in the principal stress 
                 e.Store( i, pstress_key_, evals_ );

                 // mean stress = arithmetic average of stress Eigenvalues
                 e.Store( i, means_key_, ScalarVariable(PLAIN,evals_.Average()) );
                 
                 // principal strains
                 // -----------------
                 convertTo( IPSTRAIN_, i, ts_ );
                 ts_.Eigen( evals_, evecs_, true );
                
                 e.Store( i, strain1_key_, evecs_.Row(0) );
                 e.Store( i, strain2_key_, evecs_.Row(1) ); // O.K.
                 
                 // principal strains = stretches in the principal elongation 
                 // directions = Eigenvalues of strain tensor   
                 e.Store( i, pstrain_key_, evals_ ); // O.K.
                                    
                 // dilatation = sum of absolute principal strains i.e. Eigenvalues
                 e.Store( i, dilat_key_, ScalarVariable(PLAIN,evals_[0]+evals_[1]) );
             }
         }

    
    // only once the strains and stresses have been computed, these can be output to Model<2U> 
    if ( MathOperatorLHS<2U>::ApplicationCycle() == 2U ) 
      {
         for ( size_t i=0; i<e.Nodes(); i++ )
           // doing this operation only once per node
           if ( !node_output_[ e.N(i)->Idx()] )
             {
                // 1. assigning the strain & stress values
                // ---------------------------------------
                convertColumnTo( STRAIN_, i, ts_ );
                e.N(i)->Store( strain_key_, ts_ );
                convertColumnTo( STRESS_, i, ts_ );
                e.N(i)->Store( stress_key_, ts_ );

                // 2. Computing principal strain and stress axis if requested
                // ----------------------------------------------------------
                if ( principal_e_and_sigma_ )
                  {
                     // principal stresses
                     // ------------------
                     ts_.Eigen( evals_, evecs_, true );
                     
                     e.N(i)->Store( sigma1_key_, evecs_.Row(0) );
                     e.N(i)->Store( sigma2_key_, evecs_.Row(1) );
                     
                     // principal stresses = stress magnitudes in the principal stress 
                     e.N(i)->Store( pstress_key_, evals_ );

                     // mean stress = arithmetic average of stress Eigenvalues
                     e.N(i)->Store( means_key_, ScalarVariable(PLAIN,evals_.Average()) );

                     // principal strains
                     // -----------------
                     // getting the strain back
                     convertColumnTo( STRAIN_, i, ts_ );
                     ts_.Eigen( evals_, evecs_, true );
                     
                     e.N(i)->Store( strain1_key_, evecs_.Row(0) );
                     e.N(i)->Store( strain2_key_, evecs_.Row(1) ); // O.K.
                     
                     // principal strains = stretches in the principal elongation 
                     // directions = Eigenvalues of strain tensor   
                     e.N(i)->Store( pstrain_key_, evals_ ); // O.K.
                                        
                     // dilatation = sum of absolute principal strains i.e. Eigenvalues
                     e.N(i)->Store( dilat_key_, ScalarVariable(PLAIN,evals_[0]+evals_[1]) );
                 }
                  
                // 3. flagging the node to prevent further computations
                // ----------------------------------------------------
                node_output_[ e.N(i)->Idx() ] = true;
             }
        }

 } // end WriteOperands

// USEFUL OUTPUT
//    cout <<"\nstrain (3-components x ips): ";
//    for ( typename vector::const_iterator 
//          it=IPSTRAIN_.begin(); it!=IPSTRAIN_.end(); it++ ) cout << (*it) <<" ";



 
#ifndef _MSC_VER
template class StressesAndStrains<2U>;
#endif

} // csmp

