#ifndef ALTERATION_VISITOR_H
#define ALTERATION_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"
#include "ODE_StiffSolver.h"
#include "Reactant.h"

// copyright (c) 8-15-1999 by Stephan Matthai and Sebastian Geiger (OSU) 
// last modified 8-23-99 SG
// last modified 9-2-99 SKM
// completely revised and commented 6-10-99 SKM

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/**

@brief Permforms node-by-node computations of aqueous - solid phase equilibria
using the ODE_StiffSolver in a reaction based aproach.
 
@class AlterationVisitor AlterationVisitor "visitors/AlterationVisitor.h"
@author S.K. Matthaei
@author S. Geiger
@date 1999

@section motivation Motivation
 

In order to carry out rate-dependent reaction calculations, one needs 
to visit the nodes or elements of the finite-element mesh, retrieve 
aqueous species concentrations (moles m-3) and mineral fractions (X), 
respectively, carry out speciation calculation for the alloted time 
increment, and map the results back to the Model. 

Since full equilibrium calculations are meaningful only for 
fracture-channelled flow and its interaction with the matrix but not for 
matrix flow where disequilibrium may persist between pore fluid and matrix 
minerals, the scheme should also allow rate-dependent calculations.   

The scheme should be simple. This is achieved by the SNIA approach. 
 
 
@section design Design Intent

Encapsulation appears critical in order to deal accurately with the 
complexity of the aqueous systems. It is achieved by delegating 
the numerics of the equilibration calculation to a solver object
for stiff systems of ordinary differential equations. Critical properties
of aqueous solutes, minerals, and gases are encapsulated into
Reactant objects.  

Zero-dimensional reaction calculations can be carried out outside of 
a more complex CSP-reactive-transport model. This is the method of choice
for designing and testing chemical systems before looking at their spatial
implications by means of a CSP model. Design and testing is done
by iteratively creating and modifying a reaction file which is used
to initialize the solver for ordinary differential equations, both inside
and outside of CSP. The design-proofed text file is then read by the 
AlterationVisitor to initialize the reaction calculations. In this 
approach, the dependent variables defined in the input file must 
also be known to the PropertyDatabase. Hence, they must be incorporated
into the 'CSP_variables.txt' file. In addition to the reactions file
the constructor of the AlterationVisitor also expects a datafile 
which lists the reactants again, specifying their phase state (solid, aqueous, 
gas) charge, density, and molar weight, as applicable. This file is 
read after the ODE solver is initialized and it is tested for a 
mismatch of species names. 

Robustness is quint-essential since an inability of the A.V. to treat
the chemical system at just a single node would bring the whole reactive
transport model to a halt. This goal is achieved only partly. 
Importantly, as is documented in more detail for the 'ODE_StiffSolver' 
object, rate coefficients in a single system may not vary by more than
about 12 orders of magnitude among the different dependent variables. 
Otherwise a solution cannot be assured, while a calculation may still 
perform correctly in most cases. Usually, a large number of bad
steps taken by the ODE solver is an indication that a chemical 
system is problematic.  

Flexibility is dictated by the need to iteratively refine the aqueous 
system calculation and to test the functionality of such refinements for 
the modeled geological system of interest. To achieve this goal, 
hardcoding of reactant names is avoided except for the variable H2O_aq
which must be present and which is also used to decide whether a 
chemical calculation is carried out or not: If the variable H2O_aq
is flagged DIRICH, chemical calculations are prohibited at the respective
Node.  

Apart from the variable "H2O_aq" some physical properties are also
hardwired since they are likely to be needed for every conceivable reactive
transport calculation.Their names are  

@code
"nodal porosity"
"fluid density"
"rock density"
"nodal fluid source"
"H2O_aq"
@endcode

They must be specified in the "CSP_variables.txt" file if the 
AlterationVisitor object shall be used. 

Efficiency of the AlterationVisitor is very important since calls
to its Visit() method usually prove to be the main time sink of a 
transient CSP calculation (worse than advection requirements). The 
main culprit in terms of runtime is the ODE_StiffSolver::Solve()
method. For this reason, the invocation of the solver is limited 
to scenarios where, either, the absolute concentration / volume 
fraction of a mineral has changed by a user-specified factor.
This factor can be set with the method 'ThresholdReactantChangeFactor()',
An additional criterion is whether equilibration was obtained in the 
previous timestep. The latter is indicated by the rates of change 
of the reactants in the final finite-difference timestep taken by the 
solver during the last visitation of the node. Again, a threshold rate 
above which the solver will be invoked can be set using the method 
'ThresholdReactantChangeRate()'. 
 
 
@section applicability Applicability
 
This specific implementation of the AlterationVisitor is suitable for 
arbitrary chemical calculations in which the reactants are specified as 
CSP node variables.
 
 
@section structure Structure

See general description of visitors in the User's Guide.
 
 
@section participants Participants

Participating objects are the 'ODE_StiffSolver' and the 'Reactant'
classes.
 
 
@section collaboration Collaboration

Alike other visitors, the AlterationVisitor achieves its functionality 
through collaborations with the Model, Region and Node objects.
 
 
@section implementation Implementation

Each visitation involves reading nodal properties from the Model
storage (1), a conversion of these mineral fractions and concentrations 
per unit volume (m3) into molalities (2), an invocation of the ODE solver if
reactant concentrations / fractions changed above a tolerance value 
since the last visitation or if the previous rate of change was still high
(4). The computed molalities are converted back into concentrations and 
fractions and stored back into the Model (5). From these results, 
a new porosity and rock density is calculated (6), and a potential 
"nodal fluid source" term is computed and stored (7). Finally,
the absolute reactant quantities and reaction rates for each dependent 
variable are stored as solver invocation criteria for the next timestep (8).

The results of each visitation can be output if the visitor is set 
to verbose mode. This is done using the second argument of the 
Equilibrate() method and is the default.
 
 
@section examples Application Examples

An AlterationVisitor is instantiated after a Model was constructed:  

@code
     AlterationVisitor  reactor( sg_model, "rea-file1", time_increment );
@endcode

and is applied, for instance in a time-stepping loop by calling:  

@code
  	reactor.Equilibrate( model, false );
@endcode

In this example the verbose output of each speciation calculation 
is disabled by setting the second 'boolean' argument.
 
*/
template<uint32_t dim>
class AlterationVisitor : public Visitor<dim> {
  public:
    AlterationVisitor( Model<dim>& sg, const char* rea_file, double dt );
    virtual ~AlterationVisitor();

    void Equilibrate( Model<dim>& sg, bool show_results=true );
    
    void AdjustTimeIncrement( double new_dt );

    void Precision( double tol );

    void ThresholdReactantChangeFactor( double fac );

    void ThresholdReactantChangeRate( double rate );

    virtual void Visit( Node<dim>* ); 
    
    void Verbose( bool print_all_results ); // true the default

  private:
    std::string              group_name;  // if vis is restricted to Region
    const PropertyDatabase<dim>&  pref;
    double                 val_tolerance, 
                             rate_tolerance, 
                             fvol, fmass, 
                             rho_r_new, sum;
    ScalarVariable           sc, phi, rho_f, fsource;
    std::vector<bool>        dirichlet;
    bool                     verbose, visited;
    csmp::Index              node_porosity_key,
                             rhof_key,
                             rhor_key,
                             source_key,
                             h2o_key;  
    size_t                   h2o_index;

    // ODE solver related variables 
    ODE_StiffSolver                               solver;
    double                                      time_increment, tolerance;
    std::vector<std::pair<int32_t,std::string> >    coefs, dependent_comps, independent_comps;
     // concentrations, and rates of change from the last timestep
    typename std::vector<std::vector<double> >  amounts, rates;

    typedef enum { SOLID=1, AQUEOUS=2, GAS=3 }  PHASE_STATE;
    // index (ODE solver)  charge or density
    std::map<uint32_t,std::pair<Index,Reactant> >  phase_properties;   

    double TestChargeBalance();
    size_t   ReadDependentVariablePhaseState( const char* fname );
};

} // csmp

#endif

// Potential Improvements
// ----------------------------------------------------------------------------------

// why don't read in nodal conc. of H2O and treat it as a "normal" aq. species as long as H2O 
// is not used as an element property to calculate stress or pressure diffusion in the matrix?
// this should allow to map the change in H2O due to mineral reactions...
// NB: H2O is a quantity conserved by the fluid flow equation. Changes due to reactions need
//     to be incorporated in these as fluid source or sink terms.
        
// CAN WE ADD A CHECK HERE THAT CHECKS HOW MUCH AND IF THE AQ. CONCENTRATIONS HAVE CHANGED SO SPECIATION
// CALCULATIONS ARE ONLY CARRIED OUT IF THERE IS ANY CHANGE IN CONCENTRATION DUE TO DIFFUSION???










