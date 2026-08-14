//
//  FaultFlowPropertyCalculator.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai & Caroline Milliotte on 12/18/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//
#include "Model.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "mechanics.h"
#include "MechanicalProperties.h"
#include "FaultFlowPropertyCalculator.h"
#include "SteadyStateDiffusor.h"
#include "CSMP_physical_constants.h"
#include "DeformationInducedPermeability.h"
#include "OpeningModeFracture.h"

using namespace std;

namespace csmp {


/**
    Background permeability of tectonically-active portions of the continental crust  
   (Manning and Ingebritsen, 1999, Rev. Geophysics 37, 127-150) 
   
    returns permeability (m2) as a function of subsurface depth (m).
    
*/
inline double  background_k_ManningIngebritsen99( double subsurface_depth )
{
   assert( subsurface_depth > 1.0e-5 );
   double temp = std::log10(0.001 * subsurface_depth);
   return std::pow( 10., -14. - 3.22 * temp);
}






/**   DistanceToPerimeter()

For each element of the fault region, the shortest distance
of its barycenter to the model boundary is computed. 
This is an analytical calculation
that does not involve any PDE integration.

@param model reference to the current model
@param fault lower-dimensional region to which the calculation is applied
@param distance_to_tipline_variable of the variable into which the resulting scalar is returned.

Shortest point-to-boundary distances are returned into the 'radius'
variable which must be placed on the element.

The function returns the maximum distance of any element barycenter from the tipline.

@section Application Application

The tipline distance is an important parameter in semi-analytical fracture aperture calculations.

@attention works only for stand-alone contiguous regions; else other nearby
objects will be considered as well.

@attention if the region is a curved surface, the tipline distance 
approximation is not accurate.


@section Messages Messages

Method reports if the output property does not have the right
placement and type.

*/
template<uint32_t dim>
double  FaultFlowPropertyCalculator<dim>::DistanceToPerimeter( Model<dim>& model,
                                                               const char* fault,
                                                               const char* distance_to_tipline_variable )
 {
     if ( model.Database().Type(distance_to_tipline_variable) != SCALAR or 
          model.Database().Placement(distance_to_tipline_variable) != ELEMENT )
       throw csmp::Exception( ERROR, "FaultFlowPropertyCalculator<dim>::DistanceToPerimeter: property '",
                              distance_to_tipline_variable, "' must be a scalar property placed on the element." );
   
    Region<dim>& gr = model.Region( fault );
    const csmp::Index  fi_key = model.Database().StorageKey(distance_to_tipline_variable);
   
    if ( (gr.Nodes() - gr.InteriorNodes()) < 1 )
      throw csmp::Exception( ERROR, "FaultFlowPropertyCalculator<dim>::DistanceToPerimeter: fault '",
                             fault, "' has no perimeter nodes?!" );
   
    // for each element barycenter inside the region, find the minimum distance to the region boundary
    double  rmax(0.);
    for (  auto it=gr.CellsBegin(); it!=gr.CellsEnd(); ++it )
      {
         // getting the element's barycentre
         Point<dim> bctr = (*it)->BaryCenter();
         double  rmin(DBL_MAX);
         for ( auto nit=gr.PerimeterNodesBegin(); nit!=gr.NodesEnd(); ++nit ) {
              double distance = ((*nit)->Coordinate() - bctr).Length();
              rmin = std::min( rmin, distance );
           }
         // storing the result
         (*it)->Store( fi_key, makeScalar(FIELD_DATA,rmin) );
         rmax = std::max( rmax, rmin );
      }

    printRangeOfVariable( model, fault, distance_to_tipline_variable, true );
    return rmax;

 } // end DistanceToPerimeter






/**  DistanceFromCenter

Using the 'shortest distance to fault tip' variable,
this method identifies the centerline / plane through the current
region = fault. In the case of a circle this will just be a point.

The distance of each node point to the farthest margin is
computed and output to the target variable.

@param model csmp model
@param fault_region name of the region of surface elements to which the computation shall be applied.
@param distance_to_tipline_variable the variable that tracks the distance from the fault center
@param distance_from_center_variable the 'shortest distance to fault tip' variable computed by ShortestDistanceToTipLine().

@return Method returns 1/2 of the maximum dimension of the region. If the
computation cannot be performed NaN is returned.

@section Implementation Implementation

Relies on the distance field computed with ShortestDistanceToTipLine().

@section application Application

Implemented for semi-analytical fracture aperture calculations.

@section messages Messages

If a region only contains a single element, it is valid,
yet the tipline distance is going to be very small.
In this case a warning message is issued.

The method tests for appropriate placement of the input and 
output variables. 

@test O.K. SKM 27/8/2014

*/
template<uint32_t dim>
double  FaultFlowPropertyCalculator<dim>::DistanceFromCenter( Model<dim>& model,
                                                              const char* fault_region,
                                                              const char* distance_to_tipline_variable,
                                                              const char* distance_from_center_variable )
 {
    if ( model.Database().Type(distance_to_tipline_variable) != SCALAR or 
         model.Database().Placement(distance_to_tipline_variable) != ELEMENT )
      throw csmp::Exception( ERROR, "FaultFlowPropertyCalculator<dim>::DistanceFromCenter: property '",
                             distance_to_tipline_variable, "' must be a scalar element property." );

    if ( model.Database().Type(distance_from_center_variable) != SCALAR or
         model.Database().Placement(distance_from_center_variable) != ELEMENT )
      throw csmp::Exception( ERROR, "FaultFlowPropertyCalculator<dim>::DistanceFromCenter: property '",
                             distance_from_center_variable, "' must be a scalar element property." );

    Region<dim>& gr = model.Region( fault_region );

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    // has the input variable meaningful values?
    if ( gr.Cells() == 1U )
      csmp_error.Note( WARNING, "FaultFlowPropertyCalculator<dim>::DistanceFromCenter:", fault_region,
                        "region only contains a single element; result variable will be zero." );
    else {
        double  rmin, rmax;
        gr.MinMaxOf( distance_to_tipline_variable, rmin, rmax );
        // the distance can never be zero and should not be larger than 1000-km
        const double largest_distance(1.0e6);
        if ( rmin < 0. or rmax > largest_distance ) {
             cerr<< "\n "<< distance_to_tipline_variable<< " max = "<< rmax<<endl;
             cerr<< "\n "<< distance_to_tipline_variable<< " min = "<< rmin<<endl;
             cerr<< "\n "<< "versus rmax_limit = "<< largest_distance <<endl;
             csmp_error.Note( ERROR, "FaultFlowPropertyCalculator<dim>::DistanceFromCenter: input property '", distance_to_tipline_variable, "' has erratic values." );
             return strtod("NAN",NULL);
          }
      }

    const csmp::Index  ra_key = model.Database().StorageKey(distance_to_tipline_variable);
    const csmp::Index  gm_key = model.Database().StorageKey(distance_from_center_variable);
    double           quot, distance, max_quot, target_distance;
    
    // for each element barycenter point inside the region, find minimum distance to center line
    // -----------------------------------------------------------------------------------------
    // TODO: SKM: this is a typical N-squared operation that gets very slow! - use other algorithm
    // distance is found by maximizing the quotient value / distance
    const double tolerance(1.0e-6);
    double  min_distance(DBL_MAX);
    for ( auto it=gr.CellsBegin(); it!=gr.CellsEnd(); it++ )
      {
         max_quot = target_distance = 0.;
         for ( auto eit=gr.CellsBegin(); eit!=gr.CellsEnd(); eit++ ) {
              distance = ((*eit)->BaryCenter() - (*it)->BaryCenter()).Length();
              quot = distance / std::max( tolerance, (*it)->Read( ra_key ) );  
              if ( quot > max_quot ) {
                   max_quot        = quot;
                   target_distance = distance;
                }
           }
         // storing the result  
         (*it)->Store( gm_key, makeScalar(FIELD_DATA,target_distance) );
         min_distance = std::min( min_distance, target_distance );
      }

    // now the distance is reduced by its minimum value
    for ( auto it=gr.CellsBegin(); it!=gr.CellsEnd(); it++ ) {
         double  distance_from_center((*it)->Read( gm_key ));
         distance_from_center -= min_distance;
         (*it)->Store( gm_key, makeScalar(FIELD_DATA,distance_from_center) );
      } 

    return printRangeOfVariable( model, fault_region, distance_from_center_variable, true );

 } // end DistanceFromCenter







/**  CalculateTipAndCenterLineDistances

     Contours faults with distance from tipline and centerline attributes
     and returms the maximum radius of the fault.
     
     @attention centerline distance is interpreted correctly only for contiguous (stand-alone) fault patches
*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::CalculateTipAndCenterLineDistances( Model<dim>& model, const set<string>& faults )
 {
    const string  distance_to_tipline_variable("tipline distance");
    const string  distance_from_center_variable("centerline distance");
   
    const csmp::Index  fzs_key = model.Database().StorageKey("fault size");
 
    for ( set<string>::const_iterator ft=faults.begin(); ft!=faults.end(); ++ft ) 
      {
         assert( (*ft) != "Model" );
         // calculates largest distance of any element barycenter to object perimeter and returns it
         double max_distance = DistanceToPerimeter( model, (*ft).c_str(), distance_to_tipline_variable.c_str() );

         // computes shortest distance of any element barycenter to center of region
         double distance_from_barycentre = DistanceFromCenter( model, (*ft).c_str(),
                                                                 distance_to_tipline_variable.c_str(),
                                                                 distance_from_center_variable.c_str() );
        
         // maximum center distance is stored in variable 'fault size'
         model.Region( (*ft).c_str() ).Store( fzs_key, makeScalar(PLAIN, max_distance * 2.) );
         cout <<"\nFaultFlowPropertyCalculator<dim>::CalculateTipAndCenterLineDistances: fault '"<< (*ft); 
         cout <<"' \n\tmax distance to tipline (m): "<< max_distance;
         cout <<", \n\tmax distance from fault barycenter (m): "<< distance_from_barycentre << endl;
      }

 } // end CalculateTipAndCenterLineDistances








/** 

Evaluates shear-, normal- and effective stress acting on fractures, then evaluate failure criteria.

Evaluates failure potential of fracture- or fault planes from the in situ stress
using the method 'IsCriticallyStressed()'.

Failure criteria: the variable 'failure' is used that can assume the following values:

                  0 -> stable state
                  1 -> opening of pre-existing fracture
                  2 -> tensile failure
                  3 -> frictional sliding
                  4 -> shear failure
                 -1 -> compaction failure
                 
Input variables used are supplied via MechanicalProperties object:

*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::EvaluateStressesAndFailurePotential( Model<dim>& model,
                                                                            const set<string>& fractures,
                                                                            const InSituStress& stress )
 {
    Standard_IO_Handler  stdio;
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    // assessing compatibility of stress measurement location with model position
    Point<dim>    model_barycenter = centerOfGravity( model );
    Point<3>      stress_measurement_location = stress.SampleLocation();
    InSituStress  stress_state(stress);
    // if the location of the models barycenter is offset by more than 100 meter relative to stress measurement
    // the user is given the choice to offset the location
    if constexpr ( dim == 3 ) {
          double offset = model_barycenter.DistanceTo(stress_measurement_location);
          if (  offset > 100. ) {
               csmp_error.Note( WARNING, "FaultFlowPropertyCalculator<dim>::EvaluateStressesAndFailurePotential:",
                                 "location of in situ stress measurement is offset by more than 100-m from model barycenter");
            
               cout <<"\n\tThe provide stress measurement location is vertically offset by "<< offset <<" m from model barycenter.";
               if ( stdio.YesNo("\tDo you want to reset stress measurement location to model barycenter for this run") ) {
                    stress_state.Depth( model_barycenter[1] );
                 }
            }
      }
   
    // input variables
    const csmp::Index  E_key(model.Database().StorageKey("Youngs modulus"));
    const csmp::Index  sxy_key(model.Database().StorageKey("shear modulus"));
    const csmp::Index  Poi_key(model.Database().StorageKey("Poissons ratio"));
    const csmp::Index  Ts_key(model.Database().StorageKey("tensile strength"));
    const csmp::Index  UCS_key(model.Database().StorageKey("unconfined compressive strength"));
    const csmp::Index  fric_key(model.Database().StorageKey("friction coefficient"));
    const csmp::Index  coh_key(model.Database().StorageKey("cohesion"));
   
    // resetting all failure modes
    model.InputPropertyValue( "failure", makeScalar(INIT_COND,0.) );
   
    for ( set<string>::const_iterator ft=fractures.begin(); ft!=fractures.end(); ++ft )
      {
         csmp::Region<dim>&    gref = model.Region((*ft).c_str());
         pair<int32_t,int32_t> region_dim = gref.ElementSpatialDimensions();
         if constexpr ( dim == 3 ) {
             // highest element dim          number of spatial element dims in region
             if ( region_dim.second != 2U or region_dim.first != 1U )
               csmp_error.Note( ERROR, "FaultFlowPropertyCalculator<dim>::EvaluateStressesAndFailurePotential:",
                                (*ft).c_str(), "failure criteria can only be evaluated on surface elements; nothing was done for this region");
           }
         else if constexpr ( dim == 2 ) {
             if ( region_dim.second != 1U or region_dim.first != 1U )
               csmp_error.Note( ERROR, "FaultFlowPropertyCalculator<dim>::EvaluateStressesAndFailurePotential:",
                                (*ft).c_str(), "failure criteria can only be evaluated on line elements; nothing was done for this region");
           }

         // read mechanical properties from regions
         // ---------------------------------------
         double  G = gref.Read( sxy_key );            // G = shear modulus = modulus of rigidity
         double  nu = gref.Read( Poi_key );           // nu = Poisson's ratio
         // double  E(2. * G*(1.-nu));                // E = Young's modulus derived from shear modulus (Mavko & Dvorkin)
         double  E = gref.Read( E_key );              // E = Young's modulus
         double  K(G*((2.*(1.+nu))/(3.*(1.-2.*nu)))); // K = bulk modulus = compressibility
         double  TS = gref.Read( Ts_key );            // TS = tensile strength
         double  mu  = gref.Read( fric_key );         // mu = friction coefficient
         MechanicalProperties  rprops( E, K, nu, TS, mu );
         rprops.C_ = gref.Read( coh_key );            // C = cohesive strength
         rprops.UCS_ = gref.Read( UCS_key );          // UCS = unconfined compressive strength

         // evaluate stresses and failure criteria
         // --------------------------------------
         // failure potential called 'failure' is assessed element by element considering
         // "overburden stress", "shear stress", "normal stress", "effective stress", "failure", "slip patch area"
         IsCriticallyStressed( model, rprops, (*ft).c_str(), stress_state );
        
      } // end region loop
      
 } // end EvaluateStressesAndFailurePotential








/**  DilatationFromFarFieldStress

Evaluates failure potential of fracture- or fault planes from the in situ stress
using the method 'IsCriticallyStressed'.

Failure criteria: the variable 'failure' is used that can assume the following values:

                  0 -> stable state
                  1 -> opening of pre-existing fracture
                  2 -> tensile failure
                  3 -> frictional sliding
                  4 -> shear failure
                 -1 -> compaction failure


Calculation of fault dilatation

Estimates of fault dilation are limited by the user-defined maximum of fault porosity.

Based on the assumption that fault width is determined by Riedel shears 
(second or third order fractures that branch at the friction angle into the rock 
adjacent to the fault trace), their length is calculated from the fault zone 
thickness and their aperture, using the extended Mode I criterion from
Cruikshank et al. (1991).

Results are returned into the variable 'dilatation'.

@attention for the algorithms to work properly each fracture or fault region
must only include discontiguous structures, ie., not conjugate crossing or 
multi-orientation groups of structures.

*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::EvaluateFailureAndDilatation( Model<dim>& model,
                                                                     const set<string>& faults,
                                                                     const MechanicalProperties& rprops,
                                                                     const InSituStress& stress,
                                                                     bool recompute_hydrostatic_fluid_pressure )
 {
    // input variables
    const csmp::Index  fai_key(model.Database().StorageKey("failure"));
    const csmp::Index  phi_key(model.Database().StorageKey("porosity"));
    const csmp::Index  thi_key(model.Database().StorageKey("thickness"));
    const csmp::Index  int_key(model.Database().StorageKey("fracture intensity"));
    const csmp::Index  sxy_key(model.Database().StorageKey("shear stress"));
    // output variables
    const csmp::Index  sn_key(model.Database().StorageKey("normal stress"));
    const csmp::Index  pf_key(model.Database().StorageKey("fluid pressure"));
    const csmp::Index  dil_key(model.Database().StorageKey("dilatation"));

    const double     phi_residual(0.00001);
    double           dilatation(0.);
   
    // 1. constructing a Mode I fracture, that is later used to evaluate the dilatation
    //    of Riedel shears and fractures traversing the fault at the friction angle
    // ----------------------------------------------------------------------------
    //                            length   Poisson's ratio, Young's modulus
    OpeningModeFracture  fracture( 100.,   rprops.nu_,       rprops.E_ );

    // initializing the variables 'dilatation' and 'failure' to zero
    model.InputPropertyValue( "dilatation", makeScalar(INIT_COND,0.) );
    model.InputPropertyValue( "failure", makeScalar(INIT_COND,0.) );
    if ( recompute_hydrostatic_fluid_pressure ) {
          cout <<"\nFaultFlowPropertyCalculator<dim>::EvaluateFailureAndDilatation: ";
          cout <<"To compute hydrostatic pf gradient; enter fluid pressure at model top (Pa): ";
          double pf_model_top;
          cin >> pf_model_top;
          CalculateLinearFluidPressure( model, pf_model_top, 1000. );
          cout <<"\n\n\tRecomputed 'fluid pressure' assuming it is hydrostatic (rho_f=1000kg/m3).\n";
          printRangeOfVariable( model,"fluid pressure", true );
      }
   
    // 2. apply failure criteria on all fault surface elements and calculate
    //    the dilatation of critically stressed fault segments
    // -------------------------------------------------------
    for ( set<string>::const_iterator ft=faults.begin(); ft!=faults.end(); ++ft ) 
      {
         csmp::Region<dim>&  gref = model.Region((*ft).c_str());

         // diagnosing failure potential of current element
         // computes: "overburden stress", "shear stress", "normal stress", "effective stress", "failure", "slip patch area"
         IsCriticallyStressed( model, rprops, (*ft).c_str(), stress );

         const auto it_end(gref.CellsEnd());
         for ( auto it=gref.CellsBegin(); it!=it_end; ++it )
           {
              // reading variable 'failure' to guide dilatation calculation
              long failure_type = static_cast<long>((*it)->Read( fai_key ));
              double syy = (*it)->Read( sn_key );
              double sxy = (*it)->Read( sxy_key );
              ScalarVariable pf;
              (*it)->PropertyValueAtBaryCenter( pf_key, pf );
              double phi = (*it)->Read( phi_key );
             
              // calculating length of fault-traversing fractures that are used as proxy for fault
              //  length = sin(friction angle) / fault half-width
              const double fault_half_width = (*it)->Read( thi_key ) / 2.;
              const double length = sin(rprops.alpha_) / fault_half_width;
              // setting length of fracture object
              fracture.Length( length );
              // reading fracture intensity that is used later
              const double frac_intensity   = (*it)->Read( int_key );

              // case 1: if fault patch is stable there is no dilatation
              if ( failure_type == 0 ) (*it)->Store( dil_key, makeScalar(PLAIN,0.) );
             
              // case 2: compressive failure is assumed to destroy current porosity down to residual value
              //         this is used to calculate negative dilatation
              else if ( failure_type == -1 ) {
                   dilatation = -(phi-phi_residual);
                   (*it)->Store( dil_key, makeScalar(PLAIN,dilatation) );
                }
              // case 3: dilatation by fracture opening or tensile failure
              //         Idea: dilatation is due to aperture changes in network
              //               of mode I fractures with given intensity
              //              (to calculate this, a block of rock with a size of the largest fracture is assumed)
              else if ( (failure_type == 1 or failure_type == 2) ) {
                   // calculating maximum dilatation of fault, ignoring shear stress
                   // (using syy + pf = sse,  sxy = 0)
                   dilatation = fracture.Dilatation( frac_intensity, pf(), syy, sxy );
                   (*it)->Store( dil_key, makeScalar(PLAIN,dilatation) );
                }
              // case 4: dilatation during slip or shear-failure
              //         Idea: effective stress control, but only 1/2 of max phi can be realized
              else if ( (failure_type == 3 or failure_type == 4) ) {
                   // calculating maximum dilatation of fault in shear, ignoring opening
                   const double shear_stress = fabs((*it)->Read( sxy_key ));
                   dilatation = fracture.Dilatation( frac_intensity, 0., 0., shear_stress );
                   (*it)->Store( dil_key, makeScalar(PLAIN,dilatation) );
                }
              else
              throw csmp::Exception( ERROR, "FaultFlowPropertyCalculator<dim>::EvaluateFailureAndDilatation:",
                                    "'failure' diagnostics were not recognized (valid range: -1-4).");
           }
        
      } // end region loop
      
 } // end EvaluateFailureAndDilatation






/**
    Semi-analytic aperture model that uses maximum fracture dimension to compute constant aperture value
    that gets stored to variable 'thickness'.
    
    Input parameters read from the model: 'tipline distance', 'failure', 'Youngs modulus', 'Poissons ratio', 'normal stress', 'fluid pressure'

    Output parameters written to model: aperture called 'thickness' and 'permeability'
    
    @attention minimum aperture is limited to original value.
 
*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::ConstantFractureApertureFromCruikshankModel( Model<dim>& model,
                                                                                    const set<string>& fractures )
 {
    // input variables
    const csmp::Index  fai_key(model.Database().StorageKey("failure"));
    const csmp::Index  syy_key(model.Database().StorageKey("normal stress"));
    const csmp::Index  pf_key(model.Database().StorageKey("fluid pressure"));
    ScalarVariable     pf;
    // regional properties
    const csmp::Index  Poi_key(model.Database().StorageKey("Poissons ratio"));
    const csmp::Index  E_key(model.Database().StorageKey("Youngs modulus"));
    // output variables
    const csmp::Index  thi_key(model.Database().StorageKey("thickness")); // = aperture
    const csmp::Index  kef_key(model.Database().StorageKey("permeability"));
    double           radius_min, radius_max;
   
    // 2. apply failure criteria on all fault surface elements and calculate
    //    the dilatation of critically stressed fault segments
    // -------------------------------------------------------
    for ( set<string>::const_iterator ft=fractures.begin(); ft!=fractures.end(); ++ft )
      {
         csmp::Region<dim>&  gref = model.Region((*ft).c_str());
         // getting the size of the fractures
         gref.MinMaxOf( "tipline distance", radius_min, radius_max );
         // Young's modulus E is also equivalent to:  2 G (1 + nu), see Mavko et al. p. 23
         const double E_modulus = gref.Read( E_key );

         // Mode I fracture model, that is later used to evaluate the dilatation fractures
         // ------------------------------------------------------------------------------
         //                             length           Poisson's ratio,  Young's modulus
         OpeningModeFracture  fracture( 2. * radius_max, gref.Read(Poi_key),  E_modulus );

         const auto it_end(gref.CellsEnd());
         for ( auto it=gref.CellsBegin(); it!=it_end; ++it )
           {
              // reading variable 'failure' to guide dilatation calculation
              double failure_type = (*it)->Read( fai_key );

              // computing aperture in center: since effective stress is already known pf is ignored
              if ( failure_type > 0. ) {
                   double syy = (*it)->Read( syy_key );
                   double original_aperture = (*it)->Read( thi_key );
                   (*it)->PropertyValueAtBaryCenter( pf_key, pf );
                   // fracture aperture stored as thickness
                   double apt = fracture.CenterAperture( pf(), syy );
                   (*it)->Store( thi_key, makeScalar( (*it)->Status(thi_key), apt=std::max( apt, original_aperture ) ) );
                   // fracture permeability from parallel-plate law
                   (*it)->Store( kef_key, makeScalar( (*it)->Status(kef_key), (apt*apt)/12. ) );
                }
           }
        
      } // end region loop

 } // end ConstantFractureApertureFromCruikshankModel









/**
     Semi-analytic aperture model computing values from distance to fracture tip and centerline

     fracture opening, W = aperture: Cruikshank et al., 91, p. 873, eqn. 4a, but compressive stress is positive

     dx = distance from the center of the fracture
     r  = frac_radius_max
 
        W = 2. * (-syy + pf) * ((1. - nu) / mu) * sqrt(r * r - dx * dx)
     
     fracture shear displacement, U: Cruikshank et al., 91, p. 873, eqn. 4b
        
        U = 2. * sxy * ((1. - nu) / mu) * sqrt(r * r - dx * dx)
     
     applying Moab rule, U/W = -0.4: Cruikshank et al., 91, p. 873, eqn. 4c
     to get minimum dilatation due to shear:
     
     W = std::max( W, U * UW_ratio );
 
     @attention the minimum aperture is limited to the original aperture value.
*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::VariableFractureApertureFromCruikshankModel( Model<dim>& model, const set<string>& fractures )
 {
    // input variables
    const csmp::Index  tdx_key(model.Database().StorageKey("tipline distance"));
    const csmp::Index  fai_key(model.Database().StorageKey("failure"));
    const csmp::Index  syy_key(model.Database().StorageKey("normal stress"));
    const csmp::Index  sxy_key(model.Database().StorageKey("shear stress"));
    const csmp::Index  pf_key(model.Database().StorageKey("fluid pressure"));
    ScalarVariable     pf;
    // regional properties
    const csmp::Index  Poi_key(model.Database().StorageKey("Poissons ratio"));
    const csmp::Index  G_key(model.Database().StorageKey("shear modulus"));
    // output variables
    const csmp::Index  thi_key(model.Database().StorageKey("thickness")); // = aperture
    const csmp::Index  kef_key(model.Database().StorageKey("permeability"));
    ScalarVariable     x;
    double           radius_min, radius_max;
   
    // apply failure criteria on all fault surface elements and calculate
    // the dilatation of critically stressed fault segments
    const double  UW_ratio(0.4);  // Moab rule, Cruikshank et al., 1991, p.873 used to compute minimum dilatation by shear

    for ( set<string>::const_iterator ft=fractures.begin(); ft!=fractures.end(); ++ft )
      {
         csmp::Region<dim>&  gref = model.Region((*ft).c_str());
         // getting the size of the fractures
         gref.MinMaxOf( "tipline distance", radius_min, radius_max );
        
         // getting the size of the elliptical fracture
         // and the material properties required by semi-analytical crack model
         const double nu = gref.Read( Poi_key );
         const double mu = gref.Read( G_key ); // shear modulus

         const auto it_end(gref.CellsEnd());
         for ( auto it=gref.CellsBegin(); it!=it_end; ++it )
           {
              // reading variable 'failure' to guide dilatation calculation
              double failure_type = (*it)->Read( fai_key );
             
              // computing aperture in center: fluid pressure is treated as negative as it counteracts normal stress
              // Note: if there is no dilatation, the original aperture is kept
              if ( failure_type > 0. ) {
                   const double syy = (*it)->Read( syy_key );
                   const double sxy = (*it)->Read( sxy_key );

                   // fracture opening, W: Cruikshank et al., 91, p. 873, eqn. 4a, but compressive stress is positive
                   double original_aperture((*it)->Read(thi_key));

                   // interpolating fluid pressure to fracture center
                   (*it)->PropertyValueAtBaryCenter( pf_key, pf );

                   // distance from the tipline of the fracture
                   double dx = std::max( 0., radius_max - (*it)->Read(tdx_key) );
                   double W = 2. * (-syy + pf()) * ((1. - nu) / mu) * sqrt(radius_max * radius_max - dx * dx);
                   
                   // fracture shear displacement, U: Cruikshank et al., 91, p. 873, eqn. 4b
                   double U = 2. * sxy * ((1. - nu) / mu) * sqrt(radius_max * radius_max - dx * dx);
                   
                   // exact Moab rule, U/W = -0.4: Cruikshank et al., 91, p. 873, eqn. 4c
                   // used to get minimum dilatation due to shear
                   W = std::max( W, U * UW_ratio );
                
                   // limiting minimum aperture by original value, TODO: is there a better way?
                   (*it)->Store( thi_key, makeScalar( (*it)->Status(thi_key), W=std::max( W, original_aperture ) ) );

                   // fracture permeability from parallel-plate law
                   (*it)->Store( kef_key, makeScalar( (*it)->Status(kef_key), (W*W)/12. ) );
                }
           }
        
      } // end region loop

 } // end VariableFractureApertureFromCruikshankModel

  
  
  



/**
    Similar to VariableFractureApertureFromCruikshankModel, but in stead of the tipline distance,
    the centerline distance is used. 
    This has the effect that an aperture of zero is arrived at only at the most remove tips of the elongated
    fractures, elswhere the fracture is open.
    
    This option makes sense for bed-confined rectangular fractures that completely fragment the layer.
    It addresses the case documented by mineral veins, where the fractures are open along the bed boundaries.
*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::VariableApertureForElongatedFractures( Model<dim>& model, const set<string>& fractures )
 {
    // input variables
    const csmp::Index  cdx_key(model.Database().StorageKey("centerline distance"));
    const csmp::Index  fai_key(model.Database().StorageKey("failure"));
    const csmp::Index  syy_key(model.Database().StorageKey("normal stress"));
    const csmp::Index  sxy_key(model.Database().StorageKey("shear stress"));
    const csmp::Index  pf_key(model.Database().StorageKey("fluid pressure"));
    ScalarVariable     pf;
    // regional properties
    const csmp::Index  Poi_key(model.Database().StorageKey("Poissons ratio"));
    const csmp::Index  G_key(model.Database().StorageKey("shear modulus"));
    // output variables
    const csmp::Index  thi_key(model.Database().StorageKey("thickness")); // = aperture
    const csmp::Index  kef_key(model.Database().StorageKey("permeability"));
    ScalarVariable     tx, cx;
    double           radius_min, radius_max;
   
    // apply failure criteria on all fault surface elements and calculate
    // the dilatation of critically stressed fault segments
    const double  UW_ratio(0.4);  // Moab rule, Cruikshank et al., 1991, p.873 used to compute minimum dilatation by shear

    for ( set<string>::const_iterator ft=fractures.begin(); ft!=fractures.end(); ++ft )
      {
         csmp::Region<dim>&  gref = model.Region((*ft).c_str());
         // getting the size of the fractures
         gref.MinMaxOf( "centerline distance", radius_min, radius_max );
        
         // getting the size of the elliptical fracture
         // and the material properties required by semi-analytical crack model
         const double nu = gref.Read( Poi_key );
         const double mu = gref.Read( G_key ); // shear modulus

         const auto it_end(gref.CellsEnd());
         for ( auto it=gref.CellsBegin(); it!=it_end; ++it )
           {
              // reading variable 'failure' to guide dilatation calculation
              double failure_type = (*it)->Read( fai_key );
             
              // computing aperture in center: fluid pressure is treated as negative as it counteracts normal stress
              // Note: if there is no dilatation, the original aperture is kept
              if ( failure_type > 0. ) {
                   const double syy = (*it)->Read( syy_key );
                   const double sxy = (*it)->Read( sxy_key );

                   // fracture opening, W: Cruikshank et al., 91, p. 873, eqn. 4a, but compressive stress is positive
                   double original_aperture((*it)->Read(thi_key));
                
                   // interpolating fluid pressure to fracture center
                   (*it)->PropertyValueAtBaryCenter( pf_key, pf );

                   // distance from the center of the fracture
                   double dx = std::max( 0., (*it)->Read(cdx_key) );
                   double W = 2. * (-syy + pf()) * ((1. - nu) / mu) * sqrt(radius_max * radius_max - dx * dx);
                   
                   // fracture shear displacement, U: Cruikshank et al., 91, p. 873, eqn. 4b
                   double U = 2. * sxy * ((1. - nu) / mu) * sqrt(radius_max * radius_max - dx * dx);
                   
                   // exact Moab rule, U/W = -0.4: Cruikshank et al., 91, p. 873, eqn. 4c
                   // used to get minimum dilatation due to shear
                   W = std::max( W, U * UW_ratio );
                
                   // limiting minimum aperture by original value, TODO: is there a better way?
                   (*it)->Store( thi_key, makeScalar( (*it)->Status(thi_key), W=std::max( W, original_aperture ) ) );

                   // fracture permeability from parallel-plate law
                   (*it)->Store( kef_key, makeScalar( (*it)->Status(kef_key), (W*W)/12. ) );
                }
           }
        
      } // end region loop

 } // end VariableApertureForElongatedFractures





  
  
  
  
  
/**
    Compacts aperture of fractures where s_eff >= UCS, to user-specified mimimum value.
    Fixed values for porosity and permeability are assigned.
*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::CompressiveFailureAperture( Model<dim>& model, const set<string>& fractures, double closure_aperture )
 {
    const double     phi_residual(0.001); // 0.1% for the trapped trail of fluid inclusions
    const double     k_residual(1.0e-21); // fault-gouge or a sealed fracture
    // input variables
    const csmp::Index  fai_key(model.Database().StorageKey("failure"));
    // output variables
    const csmp::Index  thi_key(model.Database().StorageKey("thickness")); // aperture
    const csmp::Index  phi_key(model.Database().StorageKey("porosity"));  // residual disconnected porosity of sealed fracture
    const csmp::Index  kef_key(model.Database().StorageKey("permeability"));
    size_t elements_failed_under_compression(0U);
    bool first_call(true);
   
    // compacts fracture segments that undergo compressive failure to user-defined minimum value
    // -----------------------------------------------------------------------------------------
    for ( set<string>::const_iterator ft=fractures.begin(); ft!=fractures.end(); ++ft )
      {
         csmp::Region<dim>&  gref = model.Region((*ft).c_str());
        
         const auto it_end(gref.CellsEnd());
         for ( auto it=gref.CellsBegin(); it!=it_end; ++it )
           {
              // reading variable 'failure' to guide dilatation calculation
              long failure_type = static_cast<long>((*it)->Read( fai_key ));
             
              // case 1: compressive failure is assumed to destroy current porosity down to residual value
              //         this is used to calculate negative dilatation
              if ( failure_type == -1 ) {
                   if ( first_call ) {
                        cout <<"\nFaultFlowPropertyCalculator<dim>::ClosureAperture: detected compressive failure.\n";
                        cout <<"will assign closure aperture, porosity and permeability of: ";
                        cout << closure_aperture <<" m, "<< phi_residual <<" X, "<< k_residual <<" m2, respectively.\n";
                        first_call=false;
                     }
                   (*it)->Store( thi_key, makeScalar((*it)->Status(kef_key),closure_aperture) );
                   (*it)->Store( phi_key, makeScalar((*it)->Status(kef_key),phi_residual) );
                   (*it)->Store( kef_key, makeScalar((*it)->Status(kef_key),k_residual) );
                   elements_failed_under_compression++;
                }
           }
      } // end region loop
   
      if ( elements_failed_under_compression > 0U )
        cout <<"\nFaultFlowPropertyCalculator<dim>::ClosureAperture: compacted "<< elements_failed_under_compression <<" failed fracture elements.\n";

 } // end CompressiveFailureAperture

  








  /**  PermeabilityFromStressAndDilatation
  
       Calculates fault permeability from distance to tip, fault size, effective stress, porosity and dilatation.
       The permeability of tectonically active continental crust (Manning & Ingebritsen, 1999, Rev. Geophys.)
       is used as an initial guess of the fault permeability. 
       
       @attention the method assumes that the y-coordinate of the model = TVD
       
       @param max_aperture  is used to put a cap on fracture permeability when integrating the aperture-frequency
                            histogram for permeability.
      
       @param flow_tortuosity  is used as a parameter in the Kozeny-Carman calculation of the breccia permeability.
  */
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::PermeabilityFromStressAndDilatation( Model<dim>& model,
                                                                            const set<string>& faults,
                                                                            double max_aperture,
                                                                            double flow_tortuosity,
                                                                            bool compute_halo_permeability )
   {
       // input variables
      const csmp::Index  fin_key(model.Database().StorageKey("fracture intensity"));
      const csmp::Index  phi_key(model.Database().StorageKey("porosity"));
      const csmp::Index  FAILURE_key(model.Database().StorageKey("failure"));
      const csmp::Index  fsz_key(model.Database().StorageKey("fault size"));
      const csmp::Index  tld_key(model.Database().StorageKey("tipline distance"));
      const csmp::Index  sss_key(model.Database().StorageKey("shear stress"));
      const csmp::Index  dil_key(model.Database().StorageKey("dilatation"));
// USE DISPLACEMENT FROM SKUA in future, to replace 'tipline distance'
//      csmp::Index  dis_key(model.Database().StorageKey("cumulative fault displacement"));
      // output variables
      const csmp::Index   kfr_key(model.Database().StorageKey("fracture permeability"));
      const csmp::Index   kbr_key(model.Database().StorageKey("fault breccia permeability"));
      const csmp::Index   kef_key(model.Database().StorageKey("permeability"));
      const csmp::Index   fth_key(model.Database().StorageKey("thickness"));
      // array variable with 6 entries to represent potential k variation across fault zone
      const csmp::Index   kfz_key = (compute_halo_permeability==true) ? model.Database().StorageKey("fault zone permeability") : csmp::Index();

      // 0. setting model up for computation, initializing fluid pressure
      model.InputPropertyValue( "fracture permeability", makeScalar(INIT_COND,0.) );
      model.InputPropertyValue( "fault breccia permeability", makeScalar(INIT_COND,0.) );
     
      DeformationInducedPermeability  permCalculator(max_aperture);
      const double modeOfFragmentSize(0.03); // (m), = mode of mineral fragments making up the fault breccia / proto-cataclasite
     
      // for all faults
      for ( set<string>::const_iterator ft=faults.begin(); ft!=faults.end(); ++ft ) 
        {
           // for each fault
           csmp::Region<dim>&  gref = model.Region((*ft).c_str());
           const double  fault_size = gref.Read( fsz_key );
           assert( fault_size > 0. );
           ScalarVariable  tip_distance;
          
           // 0. establishing the range of shear stress on elements that did not fail
           double shear_stress_min(1e9), shear_stress_max(0.);
           for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
             // if the fault element is not critically stressed failure=0
             if ( fabs((*it)->Read( FAILURE_key )) < numeric_limits<double>::epsilon() ) {
                   const double  shear_stress = (*it)->Read( sss_key );
                   shear_stress_min = std::min( shear_stress_min, shear_stress );
                   shear_stress_max = std::max( shear_stress_max, shear_stress );
                }
           const double delta_shear_stress( shear_stress_max - shear_stress_min );
          
           // computing fault permeability dependent on failure status of element
           // -------------------------------------------------------------------
           // compaction failure            ->  -1. k=0.1 x (Manning & Ingebritsen (1999)
           // no failure value              ->   0 k= (Manning & Ingebritsen (1999)
           // frictional sliding,           ->   1,
           // shear failure,                ->   2,
           // opening pre-existing fracture ->   3,
           // tensile failure               ->   4,
           for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
             {
                 if constexpr ( dim == 3 ) assert( (*it)->IsLine() == false );
                 // 1. reading input properties
                 Point<dim> bxyx = (*it)->BaryCenter(); // get TVD = fabs(y-coord)
                 double  background_k      = background_k_ManningIngebritsen99( fabs(bxyx[1]) );
                 double  fractureIntensity = (*it)->Read( fin_key );
                 double  FAILURE           = (*it)->Read( FAILURE_key );
                 double  porosity          = (*it)->Read( phi_key );
                 double  dilation          = (*it)->Read( dil_key );

                 // ========================================================
                 // 2. computing deformation induced k added to background k
                 // ========================================================
                 assert( porosity > 0. );
                 const double min_porosity(1.0e-5);
                 double  keff_frac = permCalculator.keff_UniformAperture( std::max( min_porosity, porosity + dilation ), background_k, fractureIntensity );
                 double  keff_brec = permCalculator.kKC_BrecciaFromFragSize( std::max( min_porosity, porosity + dilation ), background_k, modeOfFragmentSize, flow_tortuosity );
                 assert( keff_frac > 0. );
                 assert( keff_brec > 0. );
                 // 2.1 if no failure occurred, background permeability is scaled by the shear stress
                 if ( fabs(FAILURE) < numeric_limits<double>::epsilon() ) {
                      const double scale_fac = ((*it)->Read(sss_key) - shear_stress_min) / delta_shear_stress;
                      // harmonic mean
                      keff_frac = 1. / (scale_fac / keff_frac + (1. - scale_fac) / background_k);
                      keff_brec = 1. / (scale_fac / keff_brec + (1. - scale_fac) / background_k);
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }

                 // k values dependent on failure mode (only where the fault is critically stressed, deformation k is added)
                 // --------------------------------------------------------------------------------------------------------
                 // 2.2: compressive failure: void-space collapse as unconfined compressive strength is exceeded
                 //     (this is treated in terms of Mannings and Ingebritsen (1999)
                 if ( FAILURE < 0. ) {
                      keff_frac = keff_brec = background_k;
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }
                 // 2.3: dilatation or tensile failure
                 else if ( FAILURE > 2. ) { 
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }
                 // 2.4: failure due to fault slip or shear failure (background k can be enhanced by a factor of 100)
                 else if ( (FAILURE > 0. && FAILURE <= 2.) ) {
                      const double scale_fac = ((*it)->Read(sss_key) - shear_stress_min) / delta_shear_stress;
                      // harmonic mean
                      keff_frac = 1. / (scale_fac / keff_frac + (1. - scale_fac) / (background_k*100.));
                      keff_brec = 1. / (scale_fac / keff_brec + (1. - scale_fac) / (background_k*100.));
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }
                 // 3. 'permeability' on fault plane is treated as a mix of breccia and discrete-flow aligned fracture permeability
                 /*  
                     assumptions: 
                     - faults with a size of less than 1-km size do not contain any breccia
                     - in larger faults up to 500-m distance from fault tip there also is no breccia
                     - above this distance, a weighted average of 10% breccia k and 90% fracture k is applied
                       except for where dilation is negative and k is just set 10 x the background value
                 */
                 if ( fault_size < 3000. ) {
                      (*it)->Store( kef_key, makeScalar(PLAIN,keff_frac) );
                   }
                 else {
                         (*it)->PropertyValueAtBaryCenter( tld_key, tip_distance );
                         assert( tip_distance() >= 0. );
                         if ( tip_distance() < 500. ) (*it)->Store( kef_key, makeScalar(PLAIN,keff_frac) );
                         else {
                             double fault_k = 0.9 * keff_frac + 0.1 * keff_brec;
                             (*it)->Store( kef_key, makeScalar(PLAIN,fault_k) );
                          }
                   }
                // 5. if the fault zone permeability is supposed to vary laterally these data are captured
                // ---------------------------------------------------------------------------------------
                //    At 6 points on either side stored in an array variable.
                //    These are later integrated across the FZ to obtain an average value
                //    NB: - this array contains 11 values to store the k profile
                //          the first point represents the center of the fault zone
                //        - since the distribution is symmetric about the fault, only one side is stored
                if ( compute_halo_permeability )
                  {
                     ArrayVariable  k_lateral( "fault zone permeability", model.Database() );
                     const double halfThickness((*it)->Read( fth_key ) / 2.);
                     const double platFac(0.25);
                     const double skewFac(5.);
                     // initializing the fault core data value at a[0]
                     // applying a skewed Gaussian distribution from combined sine and exponential functions
                     k_lateral(0) = keff_brec + permCalculator.k_DamageZoneFracture( halfThickness, 0., skewFac, keff_frac );
                    
                     const double dx = halfThickness / 5.;
                     double       x(dx), fault_k(k_lateral[0]);
                    
                     for ( uint32_t i=1U; i<=5U; i++ ) {
                          // applying negative exponential distribution with a flat center part the width of which can be controlled with the pf 
                          k_lateral(i)  = permCalculator.k_FaultCoreBreccia( halfThickness, x, platFac, keff_brec );
                          k_lateral(i) += permCalculator.k_DamageZoneFracture( halfThickness, x, skewFac, keff_frac );
                          fault_k      += k_lateral[i];
                          x += dx;
                       }
                     // storing the permeability array, and its thickness averaged value
                     (*it)->Store( kfz_key, k_lateral );
                     (*it)->Store( kef_key, makeScalar(PLAIN,fault_k/6.) );
                  }
             }
        }

     printRangeOfVariable( model, "fracture permeability", true );
     printRangeOfVariable( model, "fault breccia permeability", true );
     printRangeOfVariable( model, "permeability", true );


  } // end PermeabilityFromStressAndDilatation




/** ORIGINAL FUNCTION (prior to March 2013 visit)

  void FaultFlowPropertyCalculator<dim>::PermeabilityFromStressAndDilatation( Model<dim>& model,
                                                                         const set<string>& faults,
                                                                         double max_aperture,
                                                                         double flow_tortuosity,
                                                                         bool compute_halo_permeability )
   {
       // input variables
      const csmp::Index  fin_key(model.Database().StorageKey("fracture intensity"));
      const csmp::Index  phi_key(model.Database().StorageKey("porosity"));
      const csmp::Index  FAILURE_key(model.Database().StorageKey("failure"));
      const csmp::Index  fsz_key(model.Database().StorageKey("fault size"));
      const csmp::Index  tld_key(model.Database().StorageKey("tipline distance"));
      const csmp::Index  sss_key(model.Database().StorageKey("shear stress"));
// USE DISPLACEMENT FROM SKUA in future, to replace 'tipline distance'
//      csmp::Index  dis_key(model.Database().StorageKey("cumulative fault displacement"));
      // output variables
      const csmp::Index   kfr_key(model.Database().StorageKey("fracture permeability"));
      const csmp::Index   kbr_key(model.Database().StorageKey("fault breccia permeability"));
      const csmp::Index   kef_key(model.Database().StorageKey("permeability"));
      const csmp::Index   kfz_key(model.Database().StorageKey("fault zone permeability"));
      const csmp::Index   fth_key(model.Database().StorageKey("thickness"));

      // 0. setting model up for computation, initializing fluid pressure
      model.InputPropertyValue( "fracture permeability", makeScalar(INIT_COND,0.) );
      model.InputPropertyValue( "fault breccia permeability", makeScalar(INIT_COND,0.) );
     
      DeformationInducedPermeability  permCalculator(max_aperture);
      const double modeOfFragmentSize(0.03); // (m), = mode of mineral fragments making up the fault breccia / proto-cataclasite
      ArrayVariable  k_lateral( "fault zone permeability", model.Database() );
     
      // for all faults
      for ( set<string>::const_iterator ft=faults.begin(); ft!=faults.end(); ++ft ) 
        {
           // for each fault
           csmp::Region<dim>&  gref = model.Region((*ft).c_str());
           const double  fault_size = gref.Read( fsz_key );
           assert( fault_size > 0. );
           ScalarVariable  tip_distance;
          
           // 0. establishing the range of shear stress on elements that did not fail
           double shear_stress_min(1e9), shear_stress_max(0.);
           for ( vector<Element<dim>*>::const_iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
             // if the fault element is not critically stressed failure=0
             if ( fabs((*it)->Read( FAILURE_key )) < numeric_limits<double>::epsilon() ) {
                   const double  shear_stress = (*it)->Read( sss_key );
                   shear_stress_min = std::min( shear_stress_min, shear_stress );
                   shear_stress_max = std::max( shear_stress_max, shear_stress );
                }
           const double delta_shear_stress( shear_stress_max - shear_stress_min );
          
           // computing fault permeability
           for ( vector<Element<dim>*>::iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
             {
                 // 1. reading input properties
                 Point<dim> bxyx = (*it)->BaryCenter(); // get TVD = fabs(y-coord)
                 double  matrix_k          = background_k_ManningIngebritsen99( fabs(bxyx[1]) );
                 double  fractureIntensity = (*it)->Read( fin_key );
                 double  FAILURE           = (*it)->Read( FAILURE_key );
                 double  porosity          = (*it)->Read( phi_key );

                 // 2. computing deformation induced k added to background k
assert( porosity > 0. );
                 double  keff_frac = permCalculator.keff_UniformAperture( porosity, matrix_k, fractureIntensity );
                 double  keff_brec = permCalculator.kKC_BrecciaFromFragSize( porosity, matrix_k, modeOfFragmentSize, flow_tortuosity );
assert( keff_frac > 0. );
assert( keff_brec > 0. );
                 // 2.1 if no failure occurred the permeability is scaled with the shear stress acting on the fault
                 //     shear-stress variations are considered only for elements which are not critically stressed 
                 if ( fabs(FAILURE) < numeric_limits<double>::epsilon() ) {
                      const double scale_fac = ((*it)->Read(sss_key) - shear_stress_min) / delta_shear_stress;
                      // harmonic mean
                      keff_frac = 1. / (scale_fac / keff_frac + (1. - scale_fac) / matrix_k);
                      keff_brec = 1. / (scale_fac / keff_brec + (1. - scale_fac) / matrix_k);
assert( keff_frac > 0. );
assert( keff_brec > 0. );
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }

                 // 2.2 assigning values dependent on failure mode (only where the fault is critically stressed, all deformation k is added)
                 //     compressive failure: pore-space collapse when unconfined compressive strength is exceeded
                 if ( FAILURE < 0. ) {
                      keff_frac = keff_brec = matrix_k;
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }
                 // dilatation or tensile failure
                 else if ( FAILURE > 2. ) { 
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec) );
                   }
                 // 2.3 failure due to fault slip or shear failure
                 else if ( (FAILURE > 0. && FAILURE <= 2.) ) {
                      (*it)->Store( kfr_key, makeScalar(PLAIN,keff_frac/100.) );
                      (*it)->Store( kbr_key, makeScalar(PLAIN,keff_brec/100.) );
                   }
                 // 4. 'permeability' on fault plane is treated as a mix of breccia and discrete-flow aligned fracture permeability
                 //
                 //    assumptions:
                 //    - faults with a size of less than 1-km size do not contain any breccia
                 //    - in larger faults up to 500-m distance from fault tip there also is no breccia
                 //    - above this distance, a weighted average of 10% breccia k and 90% fracture k is applied
                 //      except for where dilation is negative and k is just set 10 x the background value
                 //
                 if ( fault_size < 1000. ) {
                      (*it)->Store( kef_key, makeScalar(PLAIN,keff_frac) );
                   }
                 else {
                         (*it)->PropertyValueAtBaryCenter( tld_key, tip_distance );
                         assert( tip_distance >= 0. );
                         if ( tip_distance < 500. ) (*it)->Store( kef_key, makeScalar(PLAIN,keff_frac) );
                         else {
                             double fault_k = 0.9 * keff_frac + 0.1 * keff_brec;
                             (*it)->Store( kef_key, makeScalar(PLAIN,fault_k) );
                          }
                   }
                // 5. if the fault zone permeability is supposed to vary laterally these data are captured
                // ---------------------------------------------------------------------------------------
                //    At 6 points on either side stored in an array variable.
                //    These are later integrated across the FZ to obtain an average value
                //    NB: - this array contains 11 values to store the k profile
                //          the first point represents the center of the fault zone
                //        - since the distribution is symmetric about the fault, only one side is stored
                if ( compute_halo_permeability )
                  {
                     const double halfThickness((*it)->Read( fth_key ) / 2.);
                     const double platFac(0.25);
                     const double skewFac(5.);
                     // initializing the fault core data value at a[0]
                     // applying a skewed Gaussian distribution from combined sine and exponential functions
                     k_lateral(0) = keff_brec + permCalculator.k_DamageZoneFracture( halfThickness, 0., skewFac, keff_frac );
                    
                     const double dx = halfThickness / 5.;
                     double       fault_k(k_lateral[0]);
                    
                     for ( uint32_t i=1U, x=dx; i<=5U; i++ ) {
                          // applying negative exponential distribution with a flat center part the width of which can be controlled with the pf 
                          k_lateral(i)  = permCalculator.k_FaultCoreBreccia( halfThickness, x, platFac, keff_brec );
                          k_lateral(i) += permCalculator.k_DamageZoneFracture( halfThickness, x, skewFac, keff_frac );
                          fault_k      += k_lateral[i];
                          x += dx;
                       }
                     // storing the permeability array, and its thickness averaged value
                     (*it)->Store( kfz_key, k_lateral );
                     (*it)->Store( kef_key, makeScalar(PLAIN,fault_k/6.) );
                  }
             }
        }

     printRangeOfVariable( model, "fracture permeability", true );
     printRangeOfVariable( model, "fault breccia permeability", true );
     printRangeOfVariable( model, "permeability", true );


  } // end PermeabilityFromStressAndDilatation

*/







/**  IsCriticallyStressed

Taking into account the variation of fluid and confining pressure with depth, this method evaluates the likelihood of: 

(1) opening existing fracture, 3
(2) tensile failure,           4
(3) frictional sliding,        1
(4) shear failure,             2
(5) compaction failure        -1

The integer values given in the list above are the failure criteria that are reported.

Failure is evaluated taking into account the mechanical properties: 'tensile strength', 'shear modulus',
and 'unconfined compressive strength = UCS'.

@attention For fluid pressure, a hydrostatic gradient is assumed and a vertical offset relative to the stress
measurement given as subsurface depth in the InSituStress input object is calculated.
To avoid this, you can set the depth to the center of your model.

In InSituStress object

- subsurface depth refers to the position along the y-coordinate where the stress measurement was taken
- overburden pressure is provided in Pa and would include a potential water column offshore
- rock density at the subsurface depth from where the stress measurement stems and deeper
- far-field stress = tensor with absolute magnitudes of stress

The variable 'failure' is set accordingly to a corresponding value of 1 to 4,
and -1 for compaction failure.
 
This method writes shear- normal- and overburden stress to corresponding element variables.
Element variable 'failure' records the deformation state of the surface element.

@return method returns area of the fault plane over which tensile of shear failure occured
and outputs this number to the region variable 'slip patch area' 
 
*/ 
template<uint32_t dim>
double FaultFlowPropertyCalculator<dim>::IsCriticallyStressed( Model<dim>& model,
                                                               const MechanicalProperties& rprops,
                                                               const char* model_subregion,
                                                               const InSituStress& far_field_stress )
 {
    // input / (computed) output variables
    const csmp::Index  pf_key(model.Database().StorageKey("fluid pressure"));
    // output
    const csmp::Index  sso_key(model.Database().StorageKey("overburden stress"));
    const csmp::Index  sss_key(model.Database().StorageKey("shear stress"));
    const csmp::Index  ssn_key(model.Database().StorageKey("normal stress"));
    const csmp::Index  sse_key(model.Database().StorageKey("effective stress"));
    const csmp::Index  fai_key(model.Database().StorageKey("failure"));
    const csmp::Index  spa_key(model.Database().StorageKey("slip patch area"));

    Point<dim>         nrml;
    double             sigma_n, sigma_s;
    TensorVariable<3>  insitu_stress; // field measurement that is always 3D
    ScalarVariable     pf; // fluid pressure
    double             rupture_area(0.);
   
    csmp::Region<dim>&  gref = model.Region(model_subregion);

    gref.InputPropertyValue( "failure", makeScalar(INIT_COND,0.) );
   
    cout <<"\n\n\nFaultFlowPropertyCalculator<dim>::IsCriticallyStressed: Evaluating '";
    cout << model_subregion <<"', printing '.' for each failed element: ";
    const auto it_end(gref.CellsEnd());
    for ( auto it=gref.CellsBegin(); it!=it_end; ++it )
     {
        if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL )
          nrml = normalAtFacetCenter( (*it)->N(0)->Coordinate(), (*it)->N(1)->Coordinate(),
                                      (*it)->N(2)->Coordinate(), (*it)->N(3)->Coordinate() );
       else
       if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE )
          nrml = normalOfTriangle( (*it)->N(0)->Coordinate(), (*it)->N(1)->Coordinate(), (*it)->N(2)->Coordinate() );

       else
       if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_BAR )
          nrml = (*it)->UnitNormal();
       else
       throw Exception( ERROR, "FaultFlowPropertyCalculator<dim>::IsCriticallyStressed:",
                       "attempt to compute normal on volume element rather than surface element.");
       
       // adjusting the stress tensor for the vertical position of the element
       // --------------------------------------------------------------------
      if constexpr ( dim == 3 ) {
            // computing the element-barycenter burial depth
            Point<dim> bc = (*it)->BaryCenter();
            // calculating offset in isostatic stress due to vertical distance from point of stress measurement
            // (assumption: there is only one rock density)
            double dsigma = (bc[1] - far_field_stress.Depth()) * far_field_stress.Density() * -ACC_GRAVITY;
            // storing the overburden stress
            (*it)->Store( sso_key, makeScalar(PLAIN, far_field_stress.P_conf() * far_field_stress.Sv() + dsigma) );
            // computation of depth adjusted full stress tensor
            far_field_stress.CartesianStress( insitu_stress, dsigma );
            normalAndShearStressOnPlane( insitu_stress, nrml, sigma_n, sigma_s );
         }
       // assuming that model lies in the horizontal plane, the stress measurement depth is taken as equivalent to depth
       else if constexpr ( dim == 2 ) {
            // storing the overburden stress
            (*it)->Store( sso_key, makeScalar(PLAIN, far_field_stress.P_conf() * far_field_stress.Sv()) );
            // computation of depth adjusted full stress tensor (no further rotation)
            far_field_stress.StressState().CartesianStressTensor( insitu_stress, 0. );
// TODO: provide choice between 2D in XY versus XZ plane
            // 2D model in horizontal plane (in which the normal lies)
            Point<3>  nrmlToXZ_Plane( nrml[0], 0., nrml[1] );
            normalAndShearStressOnPlane( insitu_stress, nrmlToXZ_Plane, sigma_n, sigma_s ); // plane stress
         }
         
       // following the conventions for the subsurface
       assert( sigma_n > 0. );
       (*it)->Store( ssn_key, makeScalar(PLAIN,sigma_n) );
       (*it)->Store( sss_key, makeScalar(PLAIN,sigma_s) );
       
       // evaluating the likelihood of dilation / tensile failure or slip / shear failure:
       // --------------------------------------------------------------------------------
       // (assumes that fracture-normal compressive stress is positive as is dextral displacement)
       // no failure value -> 0
       // (1) opening pre-existing fracture ->   1,
       // (2) tensile failure               ->   2,
       // (3) frictional sliding,           ->   3,
       // (4) shear failure,                ->   4,
       // (5) compaction failure            ->  -1.
       // the variable 'failure' is set to a corresponding value of 1 to 4
       (*it)->PropertyValueAtBaryCenter( pf_key, pf );
       // recording the effective stress
       (*it)->Store( sse_key, makeScalar(PLAIN,sigma_n-pf()) );
       // tensile opening of a pre-existing fracture
       if ( (sigma_n - pf()) < 0. ) (*it)->Store( fai_key, makeScalar(PLAIN,1.) );
       if ( pf() - sigma_n > rprops.TS_ ) (*it)->Store( fai_key, makeScalar(PLAIN,2.) );
       // compressive failure of fracture if UCS of rock is exceeded
       if ( sigma_n > 0. and sigma_n - pf() > rprops.UCS_ ) (*it)->Store( fai_key, makeScalar(PLAIN,-1.) );
       // Amonton's law for frictional sliding, eqn. 4.39, p. 123, Zoback 2007
       if ( (sigma_n - pf()) > 0. and sigma_s / (sigma_n - pf()) > rprops.mu_ ) (*it)->Store( fai_key, makeScalar(PLAIN,3.) );
       // Coulomb's shear failure equation
       if  ( (sigma_n - pf()) > 0. and sigma_s - rprops.mu_ * (sigma_n - pf()) > 0. ) (*it)->Store( fai_key, makeScalar(PLAIN,4.) );
       // recording the area on the plane of the fault where tensile or shear failure would occur
       if ( (*it)->Read( fai_key ) > 0. ) {
             rupture_area += (*it)->Volume();
             cout <<".";
             cout.flush();
          }
     }
   
   // writing out the area over which slip occurred in the current fault
   if ( rupture_area > 0. )
     cout <<"\n\nFaultFlowPropertyCalculator<dim>::IsCriticallyStressed: total area of failed segment(s): "<< rupture_area <<" m2.\n";
   gref.Store( spa_key, makeScalar(PLAIN,rupture_area) );

   return rupture_area;

 } // end IsCriticallyStressed




 /**
     Linear variation of fluid pressure with depth for the mechanical calculations
 */
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::CalculateLinearFluidPressure( Model<dim>& model, double pf_model_top, double fluid_density )
  {
     const csmp::Index  pf_key(model.Database().StorageKey("fluid pressure"));
    
     // finding the max Y-coordinate to set minimum depth
     Region<dim>& mref(model.Region("Model"));
     auto nit=mref.NodesBegin();
     double min_depth((*nit)->y());
     while( nit != mref.NodesEnd() ) {
          min_depth = std::max( min_depth, (*nit)->y() );
          nit++;
       }

     // assigning fluid pressure values from pf = rho_f * g * z + p_top
     for ( nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ ) {
           const double pf = pf_model_top + (min_depth - (*nit)->y()) * fluid_density * fabs(ACC_GRAVITY);
           (*nit)->Store( pf_key, makeScalar(PLAIN,pf) );
       }
    
 } // end CalculateLinearFluidPressure
  



/**  DamageZonePermeability

   Extrapolates 'fracture permeability' and 'breccia zone permeability' along normal into damage zone for given fault width, taking into account cumulative displacement.
   
   The extrapolation is done on 5 evenly spaced points using a negative exponential function.
   The result is stored in an array variable with 11 entries, one for center.
    
*/
template<uint32_t dim>
void FaultFlowPropertyCalculator<dim>::CalculateLateralPermeabilityVariation( Element<dim>* surf_e,
                                                                              double /* breccia_width */,
                                                                              double /* damage_zone_width */,
                                                                              const csmp::Index& )
 {
    throw csmp::Exception( ERROR, "FaultFlowPropertyCalculator<dim>::CalculateLateralPermeabilityVariation",
                                  "method not fully implemented yet");
                                  
    assert( surf_e->FE()->IsSurface() );
   
    // 1. getting unit normal to surface element
    Point<dim>  nrml;
    if ( surf_e->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL )
      nrml = normalAtFacetCenter( surf_e->N(0)->Coordinate(), surf_e->N(1)->Coordinate(),
                                  surf_e->N(2)->Coordinate(), surf_e->N(3)->Coordinate() );
    else
    if ( surf_e->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE )
       nrml = normalOfTriangle( surf_e->N(0)->Coordinate(), surf_e->N(1)->Coordinate(), surf_e->N(2)->Coordinate() );
    else
    throw Exception( ERROR, "FaultFlowPropertyCalculator<dim>::CalculateLateralPermeabilityVariation",
                    "attempt to compute normal on volume element rather than fault surface element");
                    
    // TODO: actual code is missing!
   
 } // end DamageZonePermeability


template class FaultFlowPropertyCalculator<2>;
template class FaultFlowPropertyCalculator<3>;


} // end csmp
