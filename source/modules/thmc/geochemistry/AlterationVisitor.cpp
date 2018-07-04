#include "AlterationVisitor.h"
#include "Region.h"
#include "Model.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {


/**
 
The constructor of the AlterationVisitor reads the property indices 
of the required variables from the property database and initializes, 
from two corresponding text files(".txt", "-prop.txt"), the solver for 
the system of ordinary differential equations which describes the reactions.
A check is performed whether the reactants in the 2 input files 
correspond to another or whether reactants are missing from the 
"-prop.txt" file.  

The constructor also initializes the vectors which will hold the reactant
quantities and reaction rates used to test whether the system
is out of equilibrium such that the solver must be invoked. In
order to always invoke the solver during first visitations, the absolute
quantities are set to 0.0 and the rates are set to 1.0.  

@section arguments Input Arguments 

The first constructor argument is a reference to the Model which the
AlterationVisitor uses to query and initialize a range of variables. 

The second input argument specifies the name of a text file.
This name specifies from where the chemical reactions
are input to initialize the ODE solver. It is specified without the 
extension ".txt" which is however expected to follow the filename.
A second textfile must exist containing critical properties of the 
reactants (state: solid, aqueous, gas; charge or density, and molar weight).
The name of this file must be the same as that of the reaction file, but
with the suffix "-props" attached before the extension ".txt". 

The third argument fixes the time increment over which reactions shall
be carried out at each node that is visited. If equilibrium calculations
shall be carried out the time increment must be set to a large value.
1.0e+10 will do. Because of a sophisticated algorithm used this larger
increment will rarely noticeably slow the calculations down. 

@section implementation Implementation

The input text files are read by the 'ODE_StiffSolver' and by the private
method 'ReadDependentVariablePhaseState()'. The variable indices are 
retrieved from the PropertyDatabase. 

@section messages Messages 

The 'ODE_StiffSolver' and the method 'ReadDependentVariablePhaseState()'
will report reading problems and/or the absence of the necessary text 
files. In this constructor a test is performed whether the number of 
depedendent variables in the reaction file is the same as in the 
properties input file (it should be). If not, an error is reported. 
*/
template<size_t dim>
AlterationVisitor<dim>::AlterationVisitor( Model<dim>& sg, const char* rea_file, double64 dt )
    : pref(sg.Database()), 
      group_name("Model"),
      time_increment(dt),
      solver( vector<SumOfProductsWithExponents>(1) ),
      tolerance(1.0e-12),      // tolerance for the ODE Solver
      val_tolerance(1.5),      // tolerance factor for absolute concentration values
      rate_tolerance(1.0e-13), // maximal rates of change before solver is invoked
      verbose(true),
      node_porosity_key(pref.StorageKey("nodal porosity")),
      rhof_key(pref.StorageKey("fluid density")),
      rhor_key(pref.StorageKey("rock density")),
      source_key(pref.StorageKey("nodal fluid source")),
      h2o_key(pref.StorageKey("H2O_aq")),
      amounts(sg.Region("Model").Nodes()),
      rates(sg.Region("Model").Nodes()),
      dirichlet(sg.Region("Model").Nodes(),false),
      visited(false),
      h2o_index(UINT_MAX)
  { 
     Visitor<dim>::ApplicationLevel(MODEL);    // This visitor does not visit elements but the nodes
     Visitor<dim>::ApplicationTarget(NODE);

     // setting up the ODE solver -> must species be in same order as ODE ins ODEsolver?
     // --------------------------------------------------------------------------------
     char reaction_file_name[NAME_STRING], prop_file_name[NAME_STRING];
     strcpy( reaction_file_name, rea_file );
     strcpy( prop_file_name,     rea_file );
     strcat( prop_file_name,     "-prop" );
     strcat( prop_file_name,     ".txt" );
     strcat( reaction_file_name, ".txt" );

     solver.ReadODEsFrom( reaction_file_name, coefs, dependent_comps, independent_comps );
     solver.Out();
 
     // finding the ODE solver index for H2O
     // ------------------------------------ 
     typename vector<pair<int32,string> >::iterator  dit;
     
     for ( dit=dependent_comps.begin(); dit!=dependent_comps.end(); dit++ ) 
       if ( (*dit).second == "H2O_aq" ) h2o_index = static_cast<size_t>((*dit).first + 1);
     if ( h2o_index == UINT_MAX )
       cout <<"\nAlterationVisitor::(ctor): Error: variable H2O_aq not found in input file."<< endl;
     
     // reading the Index and properties for each of the reactants
     // --------------------------------------------------------------
     if ( ReadDependentVariablePhaseState( prop_file_name ) != dependent_comps.size() )
       {
          cout <<"\nAlterationVisitor::(ctor): Error: wrong number of species read from file:"<< endl;
          cout <<"Read: "<< phase_properties.size() <<" vs. "<< dependent_comps.size() << endl;
       }
     
     // setting up the vectors which will hold the absolute amounts of species and the 
     // rates of change at the last visitation of the element
     // -----------------------------------------------------
     typename vector<vector<double64> >::iterator  it;
     for ( it=amounts.begin(); it!=amounts.end(); it++ )
       {
          // reserve storage in single foul swoops then initialize it
          (*it).reserve( dependent_comps.size() );
          for ( dit=dependent_comps.begin(); dit!=dependent_comps.end(); dit++ ) 
            (*it).push_back(0.0);
       }
     for ( it=rates.begin(); it!=rates.end(); it++ )
       {
          (*it).reserve( dependent_comps.size() );
          // rates are initialized to some finite value
          for ( dit=dependent_comps.begin(); dit!=dependent_comps.end(); dit++ ) 
            (*it).push_back(1.0);
       }
       
  } // end constructor
         
         
         
   
template<size_t dim>
AlterationVisitor<dim>::~AlterationVisitor() 
 {  
 }
 


// switch on and off
template<size_t dim>
void AlterationVisitor<dim>::Verbose( bool print_output ) 
 {
    verbose = print_output;  
 }


 
/**
 
Sets the integration time for the ODE solver.  If the realistic reaction
rates are used this time corresponds to the time which the chemical 
system has to evolve.  

@section arguments Input Arguments 

The new integration time for the ODE solver. 
*/
template<size_t dim>
void AlterationVisitor<dim>::AdjustTimeIncrement( double64 new_dt ) 
 { time_increment = new_dt; }



/**
 
Sets the accuracy to which a reactant or product concentration
shall be calculated at each reaction step of the solver. The default 
value is 1.0e-12. This is a fairly high precision which is used to 
ascertain that the solver progresses to a final solution. At lower 
accuracy the risk that the solver stalls before integrating over the 
desired timespan is increased. The maximum accuracy that can be achieved
theoretically is about 1.0e-16. This is due to the fact that the 
depedendent variables are represented in double precision
(18 significant digits). 

@section arguments Input Arguments 

The accuracy to which all reactant and product concentrations shall be
calculated. 
 */
template<size_t dim>
void AlterationVisitor<dim>::Precision( double64 tol )
 { tolerance = tol; }




/**
 
Sets the factor which controls the invocation of the solver during a
visitation. Thus, if a reactant concentration changed by more than
the the total quantity times the factor, the solver will be called.
The default factor is 1.5. 

@section arguments Input Arguments 

The multiplication factor for absolute reactant concentrations. 
*/
template<size_t dim>
void AlterationVisitor<dim>::ThresholdReactantChangeFactor( double64 fac )
 {
    val_tolerance = fac;  // tolerance factor for absolute concentration values
 }
 
 
 
 
/**
 
Sets the threshold rate of change of any reactant concentration. If this
rate was exceeded during the final step of the ODE solver in the last 
invocation of the solver, the solver will be called again in the 
current node visitation. 

@section arguments Input Arguments 

The new rate above which the solver shall be applied to the chemical
system of interest.  
*/
template<size_t dim>
void AlterationVisitor<dim>::ThresholdReactantChangeRate( double64 rate )
 {
    rate_tolerance = rate; // maximal rates of change before solver is invoked
 }
 

 
 

/**
 
Submits the AlterationVisitor to the Model initiating the visitation
of all nodes during which the reaction calculations will be carried out.

@section arguments Input Arguments 

A reference to the model (Model) to which the AlterationVisitor is
applied and a boolean flag which determines whether the visitation is
done quietly or whether the results of the speciation calculation are
output as the nodes get visited. The default is verbose. Depending on
the 'stdout' interface of the computer system in use, calling the 
visitor in verbose mode is about 5-25 % slower than if the verbose 
mode is disabled.  

@section messages Messages 

The Visit() method as well as the ODE solver report errors and other
information during their execution. 
*/
template<size_t dim>
void AlterationVisitor<dim>::Equilibrate( Model<dim>& sg, bool show_results )
 {
     cout <<"\nAlterationVisitor::Equilibrate: Carrying out speciation calculations, ";
     cout <<"visiting each of the "<< sg.Region("Model").Nodes() << " nodes..."<< endl;
     verbose = show_results;
     sg.Accept( *this );
     // acknowledging that first visit was completed
     visited = true;
 }





/**
 
The method Visit() provides the core functionality of the AlterationVisitor.
It queries the visited nodes for the aqueous species concentrations and
mineral fractions, converting these into molalities. It then initializes
the 'ODE_StiffSolver' object which carries out a speciation calculation 
evolving the chemical system toward a new equilibrium. If the timespan 
allotted for this evolution (= ODE integration time) is set to a 
large value, equilibrium will be obtained provided that the system
of ordinary differential equations which are being integrated and which 
represent the chemical reactions is well posed.  

The computed molalities are converted back into
moles per cubic metre of fluid and mineral volume fractions. From these
the new porosity and composite density of the fluid saturated rock are
computed. The amount of water which was liberated or consumed by the
reactions is translated into the variable "nodal fluid source". If
this quantity exceeds the tolerance specified for the solver, it
is added to this CSMP variable stored in the Model.  

In summary, apart from the speciation calculation, Visit() calculates the 
variable "nodal porosity", "rock density", and "nodal fluid source".
The rock density is the composite density of the minerals
and the fluid stored in the pore space. The effects of the aqueous 
solutes are ignored at this stage. 

@section arguments Input Arguments 

The method takes a pointer to the visited node as argument. It permits
access to the public interface of the Node. 

@section implementation Implementation

1. Dirichlet nodes are tested by checking the variable "H2O_aq" and 
ignoring these nodes ever after if "H2O_aq" is flagged DIRICH. 

The initial composition of the fluid rock system at a node is 
specified in terms of per-unit-volume properties. The composition of the 
pore fluid is specified in [moles m-3] and that of the rock in terms of
mineral volume fractions. The latter is equivalent to the mode of
a rock and is convienient to visualize using transparency.  

It appears as a useful convention adhered to in chemical reaction calculations 
to consider what is happening in 1 kg of fluid with dissolved species and
in contact with elements /minerals. To achieve this, one needs molalities 
of the aqueous species dissolved in the pore fluid of the rock and molar 
quantities for the minerals making up the rock matrix. This requires 
conversion of the fluid-flow / geologically oriented per-unit-volume 
quantities.  

In CSP reactive-transport simulations we carry out these conversion
before and after the speciation calculation, partly already because the 
scale of a model and the porosity of the represented rock may vary greatly 
such that a molality-based treatment is more illustrative. However, 
for a rock with a specific pore volume this requires a mapping 
of the solid phases to the per-kg reaction calculation. This mapping is
done by first converting the mineral volume fractions into moles and 
then dividing the latter by the fluid mass which resides in the pore
space of this unit volume of rock. In detail, the conversion steps
are:  

1. Calculation of fluid mass M'f from porosity, phi, and fluid
density, rho'f  

M'f = phi rho'f                                               


2. Calculation of how many moles H2O reside in one kg of fluid  

N'h2o = 1000 / molar weight'h2o


3. Conversion of dissolved aqueous species a'i in moles per cubic
metre of fluid into molalities am'i (moles per kg fluid)  

am'i = a'i / rho'f                                            


4. Conversion of mineral fractions X'i to molalities N'i   

N'i = (X'i rho'i / molar-weight'i) / M'f                       


5. Re-calculation of fluid mass after reaction occurred   

M'f(t+1) = N'f(t+1) molar-weigth'f phi                                 


6. Calculation of new porosity as 1 - sum of the newly calculated mineral 
fractions  

phi(t+1) = 1 - sum'1-n over (N'i(t+1) molar-weigth'i fmass(t+1)) / rho'i    

It follows that the molar weights of H2O and the minerals must be known
in order to calculate the porosity changes due to the reactions. These 
are queried from the Reactant objects which store the properties of 
the reactants and reside in the reactant STL map 'phase_properties'.  

Summing up the mineral-fraction mineral-density products and adding
the porosity-fluid density product gives the new composite density of
the fluid-saturated rock.  

From the deviation of the H2O molality from 55.2 upon completion
of the solver one finally calculates a "nodal fluid source", src, for the
next fluid-pressure computation step:    

src = ((am'h2o(t+1) - 55.2)/55.2) 18.02 1.0e-3 fmass) / rho'f  


@section application Application

Reaction calculations where the reactants are specified as node variables.
 

@section messages Messages 

Visit() always reports to 'stdout' how many successful as well as useless
steps the ODE solver performed to obtain a solution. Many bad steps
usually indicate a problematic reaction system.   

If the 'verbose' flag evaluates to true, Visit() will report initial and 
final molalities, final reaction rates, as well as new nodal porosities, 
fluid-rock system densities (kg m-3), and fluid sources or sinks caused 
by the reactions.  

If the ODE solver cannot find a solution it will eventually terminate 
reporting errors. It can fail either because the matrix coupling the
ODEs becomes singular or because the rate of change of some equations
becomes so fast that the minimal timestep size limit is reached.  
*/
template<size_t dim>
void AlterationVisitor<dim>::Visit( Node<dim>* n ) 
  { 
     bool    invoke_solver(false);
     double64      val;
     
     if ( verbose ) cout <<"\n\nAlterationVisitor::Visit: Visiting Node: "<< n->Idx() << endl;
     // ------------------------------------------------------------------------
     // 0. The variable H2O is used to test node whether its flagged Dirichlet.
     //    If so, the Node is ignored in future loops.
     // ------------------------------------------------------------------------
     if ( !visited ) {
          n->Read( h2o_key, sc );
          if ( sc.Flag() == DIRICH ) dirichlet[ n->Idx() ] = true;
       }
     if ( dirichlet[ n->Idx() ] ) {
          if ( verbose ) cout <<"    Dirichlet node; nothing is done."<< endl;
          return;
       }

     // ----------------------------------------------------------------------------
     // 1. assuming fluid saturation, specific fluid volume is equal to the porosity
     // ----------------------------------------------------------------------------
     n->Read( node_porosity_key, phi );
     n->Read( rhof_key, rho_f );
     fmass = phi() * rho_f(); // fluid mass in unit volume of rock (kg m-3) 

	
     // ----------------------------------------------------------------------------
     // 2. read all the moles (per unit volume [m-3] of fluid) from the Model 
     //    and convert these to moles per kg. 
     // ----------------------------------------------------------------------------
     for ( typename map<uint32,pair<Index,Reactant>  >::iterator
           it=phase_properties.begin(); it!=phase_properties.end(); it++ )
       {
          // reading concentration/amount of aqueous species / mineral
          n->Read( (*it).second.first, sc );
               
          // comparing newly read fractions / concentrations with previous values. 
          // If the new values deviate by a factor of more than val_tolerance, the 
          // ODE solver will be invoked.
          val = amounts[ n->Idx() ][ (*it).first ];
          if ( (val + val*val_tolerance) <= sc() || (val - val*val_tolerance) >= sc() ) 
            invoke_solver = true;
               
          // checking whether one of the rates of change in the previous timestep
          // was above the threshold value. If so, the solver will be invoked
          if ( fabs(rates[ n->Idx() ][ (*it).first ]) >= rate_tolerance ) 
            invoke_solver = true;
           
         // assigning concentration as initial value to the solver after
          // converting it into a molality
          // (remember that indices in species lists are 0...n-1 vs. 1...n in solver)
          if ( (*it).second.second.State() == AQUEOUS ) 
            solver( (*it).first+1 ) = (*it).second.second.MolesPerVolumeToMolality( sc(), rho_f() );
          else                                  
            solver( (*it).first+1 ) = (*it).second.second.VolumeFractionToMolality( sc(), fmass );
       }
     // if the fractions / concentrations did not change by more than the threshold 
     // factor and the rates of change in the last solution step of the previous visitation
     // indicate that equilibrium was obtained, the solver is not called and the visitation
     // is terminated.
     if ( !invoke_solver ) return;  

     // if the solver has to be called, the new initial fractions / concentrations are 
     // memorized (else, the old quantities are kept for continued monitoring of slow 
     // cumulative change in the reactant quantities).      
     for ( typename map<uint32,pair<Index,Reactant>  >::iterator
           it=phase_properties.begin(); it!=phase_properties.end(); it++ ) {
          n->Read( (*it).second.first, sc );
          amounts[ n->Idx() ][ (*it).first ] = sc();
       }
       
     if ( verbose ) {
          cout <<"\nAlterationVisitor::Visit(Node "<< n->Idx() <<"): Initial system: ";
          cout << solver.Equations() <<" equations: "<< endl << endl;
          for ( typename map<uint32,pair<Index,Reactant>  >::iterator
                it=phase_properties.begin(); it!=phase_properties.end(); it++ )
            {
               cout.width(30);
               cout << (*it).second.second.Name() <<":  "<< solver( (*it).first+1 ) << endl;
            }
          cout << endl;
          cout <<"\nInitial charge balance: "<< TestChargeBalance() << endl; 
          cout.flush();
      }   


     // ----------------------------------------------------------------------------------------------
     // 3. speciation calculation with the ODE solver object.
     // ----------------------------------------------------------------------------------------------
     // maxtime is just the time needed for full equilibration (no other meaning)
     //            maxtime        init.t-step  save-steps  desired tolerance
     solver.Solve( time_increment,  1.0e-17,       0,      tolerance, verbose );

     if ( verbose ) {
          cout <<"\nAlterationVisitor::Visit(Node "<< n->Idx() <<"): Equilibrium solution: ";
          cout << solver.Equations() <<" equations: "<< endl << endl;
          for ( typename map<uint32,pair<Index,Reactant>  >::iterator
                it=phase_properties.begin(); it!=phase_properties.end(); it++ )
            {
               cout.width(30);
               cout << (*it).second.second.Name() <<":  "<< solver( (*it).first+1 ) << endl;
            }
          cout << endl;
          cout <<"\nFinal charge balance: "<< TestChargeBalance() << endl; 
          cout.flush();
      }   
        

     // ----------------------------------------------------------------------------------------------
     // 4. converting concentrations back to moles per m3, assuming that fluid density is not affected
	 //    by the change in solute concentrations. Storing new concentrations back on the nodes.
     // ----------------------------------------------------------------------------------------------
     if ( verbose ) cout <<"\nDissolved volume fraction (X) of: "<< endl << endl;

     rho_r_new=sum=0.;
     for ( typename map<uint32,pair<Index,Reactant>  >::iterator
           it=phase_properties.begin(); it!=phase_properties.end(); it++ )
       {
          if ( (*it).second.second.State() == AQUEOUS ) 
            {
               sc() = (*it).second.second.MolalityToMolesPerVolume( solver((*it).first+1), rho_f() );
               n->Store( (*it).second.first, sc );
            }
          else if ( (*it).second.second.State() == SOLID )
            {
               sc() = (*it).second.second.MolalityToVolumeFraction( solver((*it).first+1), fmass );
               n->Store( (*it).second.first, sc );
 
               // summing mineral fraction for later calculation of new porosity
               sum       += sc();
               // computing new rock density
               rho_r_new += sc() * (*it).second.second.Density();
               
               // calculate specific solid changes (dissolution = positive change)
               if ( verbose )
                 {
                    cout.width(30);
                    cout << (*it).second.second.Name() <<": ";
                    // adding minus-sign space before positive numbers 
                    if ( (amounts[ n->Idx() ][ (*it).first ] - sc()) >= 0. ) cout <<" ";
                    cout << (amounts[ n->Idx() ][ (*it).first ] - sc()) << endl;
                 }
            }
          // memorizing fractions/molalities of reaction products for next visitation
          amounts[ n->Idx() ][ (*it).first ] = sc();

          // memorizing rates of change (dX/dt) for last step of obtained solution from ODE solver
          rates[ n->Idx() ][ (*it).first ] = solver[ (*it).first+1 ];
       }
     if ( verbose ) {
          cout <<"\nRates of reactant change in last ODE finite-difference step:"<< endl << endl;
          for ( typename map<uint32,pair<Index,Reactant>  >::iterator
                it=phase_properties.begin(); it!=phase_properties.end(); it++ )
            {  
               cout.width(27); // since /dt is appended   
               cout << (*it).second.second.Name() <<"/dt: ";
               // adding minus-sign space before positive numbers 
               if ( solver[ (*it).first+1 ] >= 0.0 ) cout <<" "; 
               cout << solver[ (*it).first+1 ] << endl;
            }
          cout << endl;
          cout.flush();
       }
      

     // ----------------------------------------------------------------------------------------------
     // 5. calculating change of nodal porosity and new rock density
     // ----------------------------------------------------------------------------------------------
     sc() = 1.0 - sum;
     // adding the pore fluid contribution to the rock density
     rho_r_new += rho_f() * phi();
     if ( verbose ) {
          cout <<"\nNew fluid-saturated rock density (kg m-3):    "<< rho_r_new;
          cout <<"\nNew nodal porosity (X):                       "<< sc();
          cout <<"\nPorosity increase due to reaction (per mil): "<< (sc()-phi()) * 1.0e+3;
       }
     // avoiding propagation of erratic porosity due to slight inaccuracies 
     sc() < 0. ? 0. : sc();
     sc() > 1. ? 1. : sc();
     // storing new rock porosity and density
     n->Store( node_porosity_key, sc );
     n->Store( rhor_key, makeScalar( n->Status(rhor_key),rho_r_new) );

        
     // ----------------------------------------------------------------------------------------------
     // 6. calculating how much H2O was consumed, converting it to volume and storing it
     // ----------------------------------------------------------------------------------------------
     //         change in mole per kg fluid    ->    kg kg-1    -> kg m-3 -> m3 m-3
     sc() = (((solver(h2o_index) - 55.2)/55.2) * (18.02*1.0e-3) *  fmass) /  rho_f();
     if ( verbose ) 
       cout <<"\nFluid volume source sink (m3 m-3):            "<< sc() << endl << endl;
          
     // only if fluid source term is greater than the accuracy of ODE solver, the nodal source term is modified
     if ( fabs(solver(h2o_index) - 55.2) >= tolerance )
       {
          n->Read( source_key, fsource );
          fsource -= sc;
          n->Store( source_key, fsource );
       }
       
} // end Visit
    
    
    
    
    
    

/**
 
Reads a text file which specifies for each of the reactants in the current
system defined in the ODE solver, the phase state (solid, aqueous, gas),
either charge (aqueous) or density (solids), and the molar weight.  

The text file must have the format: 

@code
headline
species-name-as-in-reaction-file  
                   SOLID-or-AQUEOUS-or-GAS   
                                 charge-or-density-(double kg m-3)  
                                                              molar-weigth
@endcode

Since solids cannot be charged (we ignore for the moment piezo effects)
and the density of aqueous species is difficult to estimate, the third 
column in the file contains either the charge if the species is a solute
or the density if the reactant is a mineral.  

The reactant properties are stored in Reactant objects.  

The PropertyDatabase is queried with the reactant name to retrieve the 
Index of the reactant.  

For each reactant, the information obtained from the file is stored 
in a map of Index - Reactant pairs. This map will be used by the 
Visit() method to query reactant data.  

@section arguments Input Arguments 

The name (with extension) of the file which contains the information
about the reactants. 

@return The method returns for how many reactants it has retrieved information
from the text file.  

@section application Application

Currently, the method is used in the constructor of the A.V. 

@section messages Messages 

Errors are reported if the file cannot be opened, if the reactant 
names do not correspond to those in the file from which the ODE solver
is initialized, or if the reactant is unknown to the PropertyDatabase. 
Also if the phase state is mispelled an error will be reported.  
*/
template<size_t dim>
size_t  AlterationVisitor<dim>::ReadDependentVariablePhaseState( const char* fname )
 {
    phase_properties.erase( phase_properties.begin(), phase_properties.end() );
 
    typename vector<pair<int32,string> >::const_iterator  it;
    char                      intext[INFO_STRING], name[NAME_STRING];
    const char* const         white_delims =" ,\t", *token;
    int32                     index;
    int32                     charge(0);
    double64                  density(0.0), molar_weight(0.0);
    PHASE_STATE               pstate;
    pair<Index,Reactant>      entry;
    ifstream                  ifs(fname);
    
    if ( !ifs )
      {
         cout <<"\nAlterationVisitor::ReadDependentVariablePhaseState: file: "<< fname;
         cout <<" could not be openend. Nothing was done."<< endl;
         ifs.close();
         return 0;
      }
    
    // echoing the first line of the input file  
    ifs.getline( intext, INFO_STRING );
    cout <<"\nAlterationVisitor::ReadDependentVariablePhaseState: Reading file: "<< fname << endl;
    cout <<"Header: "<< intext << endl;
      
    while ( !ifs.eof() )
      {
         // 0. reading line:  "spec-name  state\n"
         ifs.getline( intext, INFO_STRING );
         // exit condition
         if ( strlen( intext ) <= 1 ) break;
         
         // 1. getting reactant name and csp_index
         token       = strtok( intext, white_delims );
         entry.first = pref.StorageKey( token ); // Index
         strcpy( name, token ); 
         
         // 2. finding the index number of the species
         for ( it=dependent_comps.begin(); it!=dependent_comps.end(); it++ )
           if ( (*it).second == name )
             {
                index = (*it).first;
                break;
             }
         if ( it == dependent_comps.end() )
           {
              cout <<"\nAlterationVisitor::ReadDependentVariablePhaseState: ";
              cout <<" Species name: "<< token <<" cannot be parsed."<< endl;
              cout <<" Nothing could be done."<< endl;
              return 0;
           }
           
         // 3. reading and parsing the state permitting capitalization
         token = strtok( NULL, white_delims );
         if      ( !strcmp( token, "SOLID" )   || !strcmp( token, "solid" ) )   pstate = SOLID;
         else if ( !strcmp( token, "AQUEOUS" ) || !strcmp( token, "aqueous" ) ) pstate = AQUEOUS;
         else if ( !strcmp( token, "GAS" )     || !strcmp( token, "gas" ) )     pstate = GAS;
         else
           {
              cout <<"\nAlterationVisitor::ReadDependentVariablePhaseState: ";
              cout <<" Species state: "<< token <<" could not be interpreted."<< endl;
              cout <<" Nothing could be done."<< endl;
              ifs.close();
              return 0;
           }

         // 4. reading the charge or density
         token = strtok( NULL, white_delims );
         if ( pstate == AQUEOUS )
           {
              // the charge is read 
              charge = atoi( token );
           }
         else // if state is solid or gaseous
           {
              // the density is read
              density = atof( token );
           }
           
         // 5. reading the molar weight of the reactant
          token = strtok( NULL, white_delims );
          molar_weight = atof( token );         
         
         // 6. assigning found values to Reactant object
         entry.second.Name( name );
         entry.second.State( static_cast<short>(pstate) );
         entry.second.Charge( charge );
         entry.second.Density( density );
         entry.second.MolarWeight( molar_weight );
           
         // 7. assigning index, state and charge or density to map
         phase_properties[ index ] = entry;
      }
    return phase_properties.size();
    
 } // end ReadDependentVariablePhaseState
    


/**
 
Sums up the products aqueous-species concentrations with their charge 
for the current reaction system and returns these. 

@return a double which holds the sum of the concentration-charge products
for the current reaction system. 

@section implementation Implementation

The species concentrations are queried from the ODE solver object. 

@section application Application

To test whether charge balance is maintained in the reacting system. 
*/
template<size_t dim>
double64 AlterationVisitor<dim>::TestChargeBalance()
 {
    double64 balance(0.);
 
    for ( typename map<uint32,pair<Index,Reactant> >::iterator
          it=phase_properties.begin(); it!= phase_properties.end(); it++ )
      // only aqueous species can be charged
      if ( (*it).second.second.State() == AQUEOUS )
        balance += solver( (*it).first+1 ) * (*it).second.second.Charge();     

    return balance;
 }
   
    
template class AlterationVisitor<1U>;
template class AlterationVisitor<2U>;
template class AlterationVisitor<3U>;

} // csmp
