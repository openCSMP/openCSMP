#include "StressesAndStrains.h"
#include "mechanics.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#include "Element.h"

#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

// 3D VERSION

/**
 
This constructor for the post-processing operator is different from
its 2D counterpart as it does not require a specification whether
plane stress or plane strain is desired. All other aspects are 
the same.

@attention if this operator gets applied to model subregion,
the arrays it creates for node acccess are too large.

@attention post-processing operator relies on valid node idx values.

*/
StressesAndStrains<3U>::StressesAndStrains( const Model<3U>& sg, 
                                            const char* oper, 
                                            const char* basic, 
                                            const char* test,
                                            bool  principal_vectors,
                                            bool extrapolate_results_to_nodes,
                                            bool geomechanics )
  : MathOperatorLHS<3U>(sg.Database(),oper,basic,test),
    strain_key_(sg.Database().StorageKey("strain")),
    stress_key_(sg.Database().StorageKey("stress")),
    DISPL_(3*3,1), 
    STIFF_(1,DenseMatrix<DM_MIN>(6,6)), EGP_(6,1), SGP_(6,1),
    STRAIN_(6,3), 
    STRESS_(6,3),
    IPSTRAIN_(3*6), // 6 strain components
    IPSTRESS_(3*6), // 6 stress components
    NSTRAIN_(3*6),
    NSTRESS_(3*6),
    NVAR_(6),
    eps_(6), sigma_(6), sum_(6),
    components_(6),
    node_output_(sg.Region("Model").Nodes(),false),
    temp_strains_(sg.Region("Model").Nodes()),
    temp_stresses_(sg.Region("Model").Nodes()),
    principal_e_and_sigma_(principal_vectors),
    geomechanics_conventions_(geomechanics),
    verbose_(false)
 {
    MathOperatorLHS<3U>::Name("StressesAndStrains", oper, basic, test );

    // testing the Operands 
    if ( strain_key_.type != TENSOR ) 
      throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                     "'strain' must be a tensor property." );

    if ( stress_key_.type != TENSOR ) 
      throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                     "'stress' must be a tensor property." );

    // displacement
    if ( MathOperatorLHS<3U>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<3U>::TestOperandType() != VECTOR )
    throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                    test, "Operand 'displacement' must be a vector property placed on the nodes." );

    // Young's modulus
    if ( (MathOperatorLHS<3U>::MaterialOperandPlacement() != ELEMENT &&
          MathOperatorLHS<3U>::MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT) ||
         MathOperatorLHS<3U>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                    oper, "Operand 'Young's modulus' must be a scalar property placed on the constraint points or element." );

    // Poisson's ratio
    if ( (MathOperatorLHS<3U>::BasicOperandPlacement() != ELEMENT &&
          MathOperatorLHS<3U>::BasicOperandPlacement() != ELEMENT_INTEGRATION_POINT) ||
         MathOperatorLHS<3U>::BasicOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                    basic, "Basic Operand 'Poisson's ratio' must be a scalar property placed on the constraint points or element." );

   // principal stresses and strains
    if ( principal_vectors ) 
      {
         strain1_key_ = sg.Database().StorageKey("strain1");
         strain2_key_ = sg.Database().StorageKey("strain2");
         strain3_key_ = sg.Database().StorageKey("strain3");
         sigma1_key_  = sg.Database().StorageKey("sigma1");
         sigma2_key_  = sg.Database().StorageKey("sigma2");
         sigma3_key_  = sg.Database().StorageKey("sigma3");
         means_key_   = sg.Database().StorageKey("mean stress");
         dilat_key_   = sg.Database().StorageKey("dilatation");
         shear_key_   = sg.Database().StorageKey("max shear stress");
        
         if ( strain1_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'strain1' must be a vector property." );

         if ( strain2_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'strain2' must be a vector property." );

         if ( strain3_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'strain3' must be a vector property." );

         if ( sigma1_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'sigma1' principal stress must be a vector property." );

         if ( sigma2_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'sigma2' principal stress must be a vector property." );

         if ( sigma3_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'sigma3' principal stress must be a vector property." );

         if ( means_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'mean stress' must be a scalar property." );

         if ( dilat_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'dilatation' must be a scalar property." );

         if ( shear_key_.type != SCALAR or (shear_key_.place != ELEMENT and shear_key_.place != ELEMENT_INTEGRATION_POINT) )
           throw csmp::Exception( ERROR, "StressesAndStrains<3D>::(constructor)", 
                   "Operand 'max shear stress' must be a scalar property placed on the element." );
     }

    // setting how many times the PDE operator shall be applied during post-processing
    if ( extrapolate_results_to_nodes ) MathOperatorLHS<3U>::ApplicationCycles(2);
    else MathOperatorLHS<3U>::ApplicationCycles(1);
    
 } // end 3D constructor








/**
 
Allows to toggle the addional computation of the variables:

@code
name          type 		placement
----          ----		---------
strain				tensor		node
'strain1' 		vector		node	(principal strain axes)
'strain2'			vector		node
('strain3')		vector		node
'dilatation'	scalar		node
'principal strain' 	vector		node	(magnitudes of elongation in the principal directions)
'stress'			tensor		node
'sigma1'  		vector		node	(principal stress axes)
'sigma2'			vector		node
('sigma3')		vector		node
'mean stress'  		scalar 		node	(= (sigma1 + sigma2 +sigma3 / 3)
'principal stress'  vector		node	(magnitudes of the principal stresses)
 
on or off. Clearly, these variables must be present in the variable 
database file for the method to work.  

@section arguments Input Arguments 

The boolean variable which specifies whether one wants the additional
computations or not. The default is 'true'.  
*/

void StressesAndStrains<3U>::PrincipalStrainsAndStresses( bool yes_no )
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
void StressesAndStrains<3U>::GetOperands( Element<3U>& e )
{
    if ( MathOperatorLHS<3U>::ApplicationCycle() == 1 )
      {
         // 1. For strain computation at integration points from nodal displacements
         // ------------------------------------------------------------------------
         // read nodal displacements
         e.NodePropertyVector( MathOperatorLHS<3U>::TestOperandKey(), NVAR_ );
    
         // 2. Subtracting mean displacement to remove translation & rotation effects
         // ------------------------------------------------------------------------
         const auto  nodes(e.Nodes());
//         vc_ = NVAR_[0];
//         for ( size_t i=1U; i<nodes; i++ ) vc_ += NVAR_[i];
//         for ( auto i{0}; i<nodes; i++ ) NVAR_[i] - vc_;

         // 3. collecting nodal displacements into a matrix of size nodes * dof x 1
         // and recording the elements spatial translation as the average of
         // its displacements
         // ----------------------------
         DISPL_.Resize( nodes * 3U, 1 );
         uint32_t k(0U);
         for ( auto i=0; i<nodes; i++ )
           for ( auto j=0; j<3; j++ ) DISPL_(k++,0) = NVAR_[i](j);

         if ( verbose_ ) {     
              cout <<"\nStressesAndStrains<"<< 3U;
              cout <<">::GetOperands: Nodal displacements, element: "<< e.Idx() << endl;
              DISPL_.Out();  
           }      
        
        // 4. For stress computation from strains at integration points
        // ------------------------------------------------------------
        // getting Young's modulus and Poisson's ratio 
        if ( MathOperatorLHS<3U>::MaterialOperandPlacement() == ELEMENT ) {
             youngs_.resize(1U);
             e.Read( MathOperatorLHS<3U>::MaterialOperandKey(),  youngs_[0] );
          }
        // it must be a constraint-point property  
        else e.IntegrationPointPropertyVector( MathOperatorLHS<3U>::MaterialOperandKey(), youngs_ );
        
        if ( MathOperatorLHS<3U>::BasicOperandPlacement() == ELEMENT ) {
             pratio_.resize(1U);
             e.Read( MathOperatorLHS<3U>::BasicOperandKey(),  pratio_[0] );
          }
        else e.IntegrationPointPropertyVector( MathOperatorLHS<3U>::BasicOperandKey(), pratio_ );

        // computing the stiffness matrix for the respective element
        stiffnessMatrix( youngs_, pratio_, STIFF_ );
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

@todo SKM: when output is element, do averaging of integration point 
stresses first and then calculate the eigenvalues, vectors, etc. 
only once.
 
*/
void StressesAndStrains<3U>::ComputeContribution( Element<3U>& e )
{
   // 1. Computing STRAINS and STRESSES at the integration points from the nodal 
   //    displacements and the interpolation function derivative matrices
   //    at the integration points
   // -------------------------------------------------------------------
   if ( MathOperatorLHS<3U>::ApplicationCycle() == 1 ) {
        if ( verbose_ )
          cout <<"\n\nStressesAndStrains<3U>::ComputeContribution: Element: "<< e.Idx() << endl; 

        // strain - stress components of symmetric matrix at integration points
        const auto integration_points(e.IntegrationPoints());
        IPSTRAIN_.resize( integration_points * components_ );
        IPSTRESS_.resize( integration_points * components_ );
        
        // 1.1 Computing the strains at the integration points
        // ---------------------------------------------------
        for ( auto i=0; i<integration_points; i++ )
          {
             //  getting DN matrices at the node points
             e.dN_AtIntegrationPoint( EGP_, i, 3U );
             
             // compute strain e_i = [B]_i{d}_i 3x12 * 12x1 -> 3x1 (in 2D)
             EGP_ *= DISPL_;

             // compute {sigma} = [E]([B]{d})
             if ( MathOperatorLHS<3U>::MaterialOperandPlacement() == ELEMENT ) 
               SGP_ = STIFF_[0];
             else // a separate stiffness matrix is defined at each integration point
               SGP_ = STIFF_[i];
             // multiplying the shear strains by a factor of 2 before computing the stress
             // via multiplication with the stiffness matrix
             EGP_(3,0) *= 2.;
             EGP_(4,0) *= 2.;
             EGP_(5,0) *= 2.;
             SGP_ *= EGP_;
             // returning them to previous value
             EGP_(3,0) /= 2.;
             EGP_(4,0) /= 2.;
             EGP_(5,0) /= 2.;

             if ( verbose_ ) { 
                  cout <<"\nstrain at integration point: "<< i << endl;
                  EGP_.Out();
                  cout <<"\nstress at integration point: "<< i << endl; 
                  SGP_.Out();
               }       

             // inserting strains and stresses sequentially into temporary 
             // STL vectors
             for ( auto k=0; k<components_; k++ ) {
                  IPSTRAIN_[ i*components_ + k ] = EGP_(k,0);
                  IPSTRESS_[ i*components_ + k ] = SGP_(k,0);
               }
          }

        // 1.2 Extrapolating stresses and strains to the nodes
        //     storing them in temporary vectors for later node
        //     averaging.
        // --------------------------------------------------------
        if ( MathOperatorLHS<3U>::ApplicationCycles() == 2 ) {
             // nodal components
             NSTRAIN_.resize( e.Nodes() * components_ );
             NSTRESS_.resize( e.Nodes() * components_ );
             e.ExtrapolateIntegrationPointVariableToNodes( components_, IPSTRAIN_, NSTRAIN_ );
             e.ExtrapolateIntegrationPointVariableToNodes( components_, IPSTRESS_, NSTRESS_ );
         
             for ( auto i=0; i<e.Nodes(); i++ )
               {
                  // extracting the nodal strain and stress components
                  for ( auto k=0; k<components_; k++ ) {
                       eps_[k]   = NSTRAIN_[ i*components_ + k ];
                       sigma_[k] = NSTRESS_[ i*components_ + k ];
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
   if ( MathOperatorLHS<3U>::ApplicationCycle() == 2 ) {
        STRAIN_.Resize(components_,e.Nodes());
        STRESS_.Resize(components_,e.Nodes());

        for ( auto i=0; i<e.Nodes(); i++ ) {
          // duplicate calculations are avoided via the boolean vector
          if ( !node_output_[ e.N(i)->Idx() ] )
            {
               for ( auto j=0; j<components_; j++ )
                 {  
                    // averaging strain components
                    double sum(0.);
                    for ( deque<vector<double> >::const_iterator 
                          lit=temp_strains_[ e.N(i)->Idx() ].begin();
                          lit!=temp_strains_[ e.N(i)->Idx() ].end(); lit++ ) sum += (*lit)[j];
                    sum /= static_cast<double>(temp_strains_[ e.N(i)->Idx() ].size());
                    STRAIN_(j,i) = sum;

                    // averaging stress components
                    sum = 0.;
                    for ( deque<vector<double> >::const_iterator
                          lit=temp_stresses_[ e.N(i)->Idx() ].begin();
                          lit!=temp_stresses_[ e.N(i)->Idx() ].end(); lit++ ) sum += (*lit)[j];
                    sum /= static_cast<double>(temp_stresses_[ e.N(i)->Idx() ].size());
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








/**
 
Outputs the computed strains to the Model<3U>  storage checking the range
of each variable using the ranges associated with each of the Operands. 
 

If principal stresses and strains were requested this method also does 
the necessary computations before outputting these variables.  

@section arguments Input Arguments 

A reference to the property memory manager and the Element for which the 
variables are output.  

@todo (1) SKM: fix so that this also does the same thing for integration 
point stresses.

*/
void StressesAndStrains<3U>::WriteOperands( Element<3U>& e )
 {
    // if stress and strain are element properties, their
    // integration point values are averaged on the element
    if ( stress_key_.place == ELEMENT )
      {
        IP_STRAIN_TENSOR_ = 0.;
        IP_STRESS_TENSOR_ = 0.;

        // 1. accumulating and averaging the strain & stress values
        // --------------------------------------------------------
        // accumulation
        for ( auto i=0; i<e.IntegrationPoints(); i++ ) {
              convertTo( IPSTRAIN_, i, ts_ );
              IP_STRAIN_TENSOR_ += ts_;
              convertTo( IPSTRESS_, i, ts_ );
              IP_STRESS_TENSOR_ += ts_;
           }
        // averaging
        const double ips(static_cast<double>(e.IntegrationPoints()));
        IP_STRAIN_TENSOR_ /= ips;
        IP_STRESS_TENSOR_ /= ips;
       
        // 2. assigning averaged values to element
        // ---------------------------------------
        e.Store( strain_key_, IP_STRAIN_TENSOR_ );
        e.Store( stress_key_, IP_STRESS_TENSOR_ );
        
        if ( principal_e_and_sigma_ ) {
             // principal strains
             if ( !IP_STRAIN_TENSOR_.Eigen( evals_, evecs_, false ) ) {
                    e.Out();
                    IP_STRAIN_TENSOR_.Out();
                    throw csmp::Exception( ERROR, "StressesAndStrains<3U>::WriteOperands:",
                                          "Eigen decompostion of strain tensor failed.");
                }
             sortEigenVectorsAndValues( evals_, evecs_ );
             e.Store( strain1_key_, evecs_.Row(0) );
             e.Store( strain2_key_, evecs_.Row(1) );
             e.Store( strain3_key_, evecs_.Row(2) );
             e.Store( dilat_key_, makeScalar(PLAIN,evals_[0]+evals_[1]+evals_[2]) );
             // principal stresses
             if ( !IP_STRESS_TENSOR_.Eigen( evals_, evecs_, false ) ) {
                  e.Out();
                  IP_STRESS_TENSOR_.Out();
                  throw csmp::Exception( ERROR, "StressesAndStrains<3U>::WriteOperands:",
                                        "Eigen decompostion of stress tensor failed.");
               }
             sortEigenVectorsAndValues( evals_, evecs_ );
             // geomechanics convention: changing the sign of the stresses
             if ( geomechanics_conventions_ ) {
                  evals_ *= -1.;
                  //evecs_ *= -1.;
               }
             // assuming that each row contains one Eigenvector
             e.Store( sigma1_key_, evecs_.Row(0) );
             e.Store( sigma2_key_, evecs_.Row(1) );
             e.Store( sigma3_key_, evecs_.Row(2) );
             e.Store( means_key_, makeScalar(PLAIN,(evals_[0]+evals_[1]+evals_[2])/3.) );
             // max shear stress
             e.Store( shear_key_, makeScalar(PLAIN,(evals_[0]-evals_[2])/2.) );
          }

      } // end element stresses and strains


    // ELEMENT_INTEGRATION_POINT OUTPUT
    // --------------------------------
    // if a single stage computation is desired, excluding extrapolations to the nodes
    if ( stress_key_.place == ELEMENT_INTEGRATION_POINT )
      for ( auto i=0; i<e.IntegrationPoints(); i++ )
        {
            // 1. assigning the strain & stress values
            // ---------------------------------------
            convertTo( IPSTRAIN_, i, ts_ );
            e.Store( i, strain_key_, ts_ );
            // principal strains
            // -----------------
            if ( principal_e_and_sigma_ ) {
                 // principal strains = stretches in the principal elongation 
                 // directions = Eigenvalues of strain tensor   
                 ts_.Eigen( evals_, evecs_, true ); // do not normalize but scale by eigenvalues
                 sortEigenVectorsAndValues( evals_, evecs_ );
                  vc_ = evecs_.Row(0);
                  e.Store( i, strain1_key_, vc_ ); 
                  vc_ = evecs_.Row(1);
                  e.Store( i, strain2_key_, vc_ ); 
                  vc_ = evecs_.Row(2);
                  e.Store( i, strain3_key_, vc_ ); 

                 // dilatation = sum of principal strains 
                 // -------------------------------------
                 sc_() = evals_[0] + evals_[1] + evals_[2];
                 e.Store( i, dilat_key_, sc_ );
              }
            
            convertTo( IPSTRESS_, i, ts_ );
            e.Store( i, stress_key_, ts_ ); 
            // principal stresses
            // ------------------
            if ( principal_e_and_sigma_ ) {
                  ts_.Eigen( evals_, evecs_, false ); // do not normalize but scale by eigenvalues
                  sortEigenVectorsAndValues( evals_, evecs_ );
                  if ( geomechanics_conventions_ ) {
                       evals_ *= -1.;
                       //evecs_ *= -1.;
                    }
                  vc_ = evecs_.Row(0);
                  e.Store( i, sigma1_key_, vc_ ); 
                  vc_ = evecs_.Row(1);
                  e.Store( i, sigma2_key_, vc_ ); 
                  vc_ = evecs_.Row(2);
                  e.Store( i, sigma3_key_, vc_ ); 

                 // mean stress = average of principal stresses
                 // -------------------------------------------
                 sc_() = (evals_[0] + evals_[1] + evals_[2]) / 3.;
                 e.Store( i, means_key_, sc_ );
                 // max shear stress
                 // ----------------
                 e.Store( i, shear_key_, makeScalar(PLAIN,(evals_[0]-evals_[2])/2.) );
              }
         }


    // NODE OUTPUT
    // -----------
    // only once the strains and stresses have been computed, these can be output to Model<3U> 
    if ( MathOperatorLHS<3U>::ApplicationCycle() == 2 ) 
      {
         for ( auto i{0}; i<e.Nodes(); i++ )
           // doing this operation only once per node
           if ( !node_output_[ e.N(i)->Idx() ] )
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
                     // principal strains
                     // -----------------
                     convertColumnTo( STRAIN_, i, ts_ );
                     ts_.Eigen( evals_, evecs_, false ); // do not normalize but scale by eigenvalues
                     sortEigenVectorsAndValues( evals_, evecs_ );
                     vc_ = evecs_.Row(0);
                     e.N(i)->Store( strain1_key_, vc_ ); 
                     vc_ = evecs_.Row(1);
                     e.N(i)->Store( strain2_key_, vc_ ); 
                     vc_ = evecs_.Row(2);
                     e.N(i)->Store( strain3_key_, vc_ ); 

                     // dilatation = sum of principal strains 
                     // -------------------------------------
                     sc_() = evals_[0] + evals_[1] + evals_[2];
                     e.N(i)->Store( dilat_key_, sc_ );
                       
                     // principal stresses
                     // ------------------
                     convertColumnTo( STRESS_, i, ts_ );
                     ts_.Eigen( evals_, evecs_, false ); // do not normalize but scale by eigenvalues
                     sortEigenVectorsAndValues( evals_, evecs_ );
                     if ( geomechanics_conventions_ ) {
                          evals_ *= -1.;
                          //evecs_ *= -1.;
                       }
                     vc_ = evecs_.Row(0);
                     e.N(i)->Store( sigma1_key_, vc_ ); 
                     vc_ = evecs_.Row(1);
                     e.N(i)->Store( sigma2_key_, vc_ ); 
                     vc_ = evecs_.Row(2);
                     e.N(i)->Store( sigma3_key_, vc_ ); 

                     // mean stress = average of stress Eigenvalues
                     // -----------
                     sc_() = (evals_[0] + evals_[1] + evals_[2]) / 3.;
                     e.N(i)->Store( means_key_, sc_ );
                     // max shear stress
                     // ----------------
                     e.N(i)->Store( shear_key_, makeScalar(PLAIN,(evals_[0]-evals_[2])/2.) );
                 }
                  
                // 3. flagging the node to prevent further computations
                // ----------------------------------------------------
                node_output_[ e.N(i)->Idx() ] = true;
             }
        }

 } // end WriteOperands




// TESTING TENSOR CALCULATIONS
/*
cerr <<"\nCartesian stress tensor:";
IP_STRESS_TENSOR_.Out();
cerr <<"\nEigenvalues:";
evals_.Out();
cerr <<"\nEigenvectors:";
evecs_.Out();
cerr <<"\n";
*/

} // csmp
