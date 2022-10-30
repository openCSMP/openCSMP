#ifndef STRESSES_AND_STRAINS2_H
#define STRESSES_AND_STRAINS2_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "TensorVariable.h"
#include "DenseMatrix.h"

namespace csmp {

template<uint32_t> class StressesAndStrains;
template<uint32_t> class Element;
template<uint32_t> class Model;

/**
    @brief 2D stress & strain are output as tensor variables, if they
    are nodal variables they are extrapolated to the nodes in two
    application cycles.
    
    @note Compressive strains and stresses are positive.
    
    @note the strain and stress component vectors are normalised to 1.
*/
template<>
class StressesAndStrains<2U> : public MathOperatorLHS<2U,Element> {
  public:
    StressesAndStrains( const Model<2U> & sg, 
                        const char* oper,   ///< Young's modulus
                        const char* basic,  ///< Poisson's ratio
                        const char* test,   ///< displacement
                        bool  plane_strain, ///< stresses will be computed suppressing any displacement perpendicular to the 2D model
                        bool  principal_vectors ); ///< compute principal stresses and strains
    
    /// get displacement values to compute {d} vector
    virtual void GetOperands( const Element<2U>& );

    /// {e} = [B]{d}, {s} = [D]{e} at integration points
    virtual void ComputeContribution( const Element<2U>& );

    /// write {e} to each elements
    virtual void WriteOperands( Element<2U>& );
    
    void PlaneStress( bool yes_no=true );
    
    inline bool IsPlaneStress() const { return !plane_strain_; }
    
    void PrincipalStrainsAndStresses( bool yes_no );
    virtual StressesAndStrains<2U>* clone() const { return new StressesAndStrains<2U> (*this); }
  protected:
    std::vector<ScalarVariable >      youngs_, pratio_;
    std::vector<DenseMatrix<DM_MIN> >  STIFF_;

  private:
    const uint32_t  components_;  /// < stress strain components
    
    csmp::Index strain_key_, stress_key_,                 ///< tensor variables
               strain1_key_, strain2_key_, strain3_key_,  ///< vector<double> variables (Eigenvectors)
               sigma1_key_, sigma2_key_, sigma3_key_,     ///< eigenvectors scaled to unit length
               pstrain_key_, pstress_key_,                ///< principal strains / stresses (Eigenvalues)
               means_key_, dilat_key_;                    ///< scalar variables

    DenseMatrix<DM_MIN>     DISPL_,
                            STRESS_, STRAIN_,
                            EGP_, SGP_, PEGP_, PSGP_,
                            PR_, EIG_;
    
    std::vector<double>   EVAL_;
    std::vector<double>   IPSTRAIN_, IPSTRESS_,
                            NSTRAIN_,  NSTRESS_,
                            eps_, sigma_, sum_;
    TensorVariable<2U>      ts_;
    VectorVariable<2U>      vc_;
    ScalarVariable          sc_;
    VectorVariable<2U>      evals_;
    TensorVariable<2U>      evecs_;
    
    bool  plane_strain_, verbose_, principal_e_and_sigma_;
  
    std::vector<std::deque<std::vector<double> > >  temp_strains_,
                                                      temp_stresses_;
    // to prevent duplicate node output
    std::vector<bool>  node_output_;
};

/**

class StressesAndStrains<2U> StressesAndStrains<2U> "mechanics/StressesAndStrains2.h"
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

}

#endif

















