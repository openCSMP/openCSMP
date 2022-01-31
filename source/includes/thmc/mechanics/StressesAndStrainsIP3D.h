#ifndef STRESSES_AND_STRAINS_IP3D_H
#define STRESSES_AND_STRAINS_IP3D_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "DenseMatrix.h"
#include "TensorVariable.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Model;

// TODO: add Biot constraints as a post-processing step, Biot-Willis coefficient alpha or drained and undrained moduli are needed

/// 3D stress & strain tensor variables placed at integration points
class StressesAndStrainsIP3D : public MathOperatorLHS<3U> {
  public:
    /// 3D constructor: note that - in contrast with civil engineering - in geomechanics compressive stresses are positive
    StressesAndStrainsIP3D( const Model<3U>&, 
                            const char* oper,   // Young's modulus
                            const char* basic,  // Poisson's ratio
                            const char* test,   // displacement
                            bool  principal_vectors,
                            bool geomechanics_conventions=true  );
  
    /// get displacement values to compute {d} vector
    virtual void GetOperands( Element<3U>& );

    /// {e} = [B]{d}, {s} = [D]{e} at integration points
    virtual void ComputeContribution( Element<3U>& );

    /// write {e} to each elements
    virtual void WriteOperands( Element<3U>& );
    
    void PrincipalStrainsAndStresses( bool yes_no );
    
  private:
    const size_t  components_;  // stress strain components
    
    csmp::Index strain_key_, stress_key_,                  ///< tensor variables
                strain1_key_, strain2_key_, strain3_key_,  ///< vector<double> variables (Eigenvectors)
                sigma1_key_, sigma2_key_, sigma3_key_,
                means_key_, dilat_key_,                    ///< mean stress (scalar), dilatation (scalar)
                shear_key_;                                ///< maximum shear stress (scalar)

    DenseMatrix<DM_MIN>                   DISPL_,
                                          STRESS_, STRAIN_,
                                          EGP_, SGP_, PEGP_, PSGP_;
    std::vector<DenseMatrix<DM_MIN> >     STIFF_;
    std::vector<double>                 IPSTRAIN_, IPSTRESS_, ///< stresses and strains at the element integration points
                                          NSTRAIN_,  NSTRESS_, 
                                          eps_, sigma_, sum_;
    TensorVariable<3U>                    ts_, evecs_;
    TensorVariable<3U>                    IP_STRAIN_TENSOR_;   ///< for averaging of integration point variables to barycentre
    TensorVariable<3U>                    IP_STRESS_TENSOR_;   ///< for averaging of integration point variables to barycentre
    VectorVariable<3U>                    vc_, evals_;         ///< some universally useful vector variable,  eigen values
    VectorVariable<3U>                    IP_e1_, IP_e2_, IP_e3_,
                                          IP_s1_, IP_s2_, IP_s3_;
    ScalarVariable                        sc_;
    std::vector<ScalarVariable>           youngs_, pratio_;
    std::vector<VectorVariable<3U> >      NVAR_;              ///< temporary container for nodal displacements
    bool                                  verbose_,
                                          principal_e_and_sigma_,
                                          geomechanics_conventions_; ///< as opposed to civil engineering where tensile stresses are positive 
};

} // csmp

/**
 
@class StressesAndStrainsIP3D<3U> StressesAndStrainsIP3D<3U> "mechanics/StressesAndStrainsIP3D.h"
@author S.K. Matthaei
@date 2000
 
@section motivation Motivation

Post-Processing: Use this operator to compute strains and stresses from nodal 
displacements using the matrix equation:

 
{e} = [B]{d}
 

In full form for two-dimensional analysis this is:
 

e_x            | d       |
               | -    0  |
               | dx      |
               |         |
               |      d  |    u
{ e_y }      = | 0    -  |  {   }
               |      dy |    v
               |         |
               | d    d  | 
gamma_xy       | -    -  |
               | dy   dx |
          
 
or for the quadratic triangular element and the shape function derivatives at
each node point:
 

    | dN0dx   0   dN1dx   0   dN2dx   0   dN3dx   0   dN4dx   0   dN5dx   0   |
B = |   0   dN0dy   0   dN1dy   0   dN2dy   0   dN3dy   0   dN4dy   0   dN5dy |
    | dN0dx dN0dy dN1dx dN1dy dN2dx dN2dy dN3dx dN3dy dN4dx dN4dy dN5dx dN5dy |            

             
which must be computed at each node for higher order elements.       

The strains and stresses are computed at the elements integration points and
then back extrapolated to the nodes. For this extrapolation linear functions
are used which use the integration points as constraint points.  

In addition to the computation of the tensor variables 'strain' and 'stress'
one has the option to compute the principal stresses and strains as well.
The method then finds the vectors which correspond to the sigma1, 2, and 3
directions which correspond to the Eigenvectors of the stress tensors.
It also calculates the mean stress and the dilation as nodal properties.
The Eigenvalues which correspond to the magnitudes of the principal
strains and stresses are stored in the variables 'principal strain' and
'principal stress'.  

The strains and stresses are averaged at the nodes among the elements 
which share each node since they diverge from element to element but
VTK needs nodal properties. This has the disadvantage that stresses
on material boundaries are somewhat unrealistic mixes.  

 
 
@section implementation Implementation

Procedure: In several post-processing application cycles of this class 
object the following is performed: 

1. From the nodal displacements and the shape function derivatives, at the
   integration points the strain and stress tensors at these locations are 
   computed (=optimal locations, see Barlow, 1977).  

2. Linearly extrapolate the strains back from the integration points 
   to the nodes. In this step, bilinear extrapolation is carried out
   by the specific finite element used. Only linear variations in strain
   over the element are supported at this stage. The nodal strains are 
   averaged among the adjacent elements since they diverge at these locations
   from one-another.  
   
3. From the nodally averaged strains and stresses all other variables are
   calculated.  
   
Steps 1 and 2 are accomplished in the first application cycle of the 
post-processing operator. The nodal averaging and following operations
are done in a second cycle.  
*/

#endif
















