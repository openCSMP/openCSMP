#ifndef IAPWS_H20_PROPERTIES_VISITOR_H
#define IAPWS_H20_PROPERTIES_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#include "DynamicViscosity.h"
// watch out: this file defined global constants ONE, TWO etc.
#include "steam4.h"

namespace csmp {

template<size_t> class Model;

template<size_t dim>
class IAPWS_H2OPropertiesVisitor : public Visitor<dim> {
  public:
    /// custom constructor using internally predefined variable names for, H, cp, alpha, and beta
    IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                const char* fluidPressure,
                                const char* fluidDensity,
                                const char* fluidViscosity,
                                bool with_steam=false,
                                bool verbose=false );

    /// custom constructor allowing the user to choose names for, H, cp, alpha, and beta
    IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                Index& fluidPressure,
                                Index& fluidTemperature,
                                Index& fluidDensity,
                                Index& fluidViscosity,
                                Index& fluidEnthalpy,
                                Index& fluidHeatCapacity,
                                Index& fluidExpansivity,
                                Index& fluidCompressibility,
                                bool with_steam=false,
                                bool verbose=false );

    IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                const char* fluidPressure,
                                const char* fluidTemperature,
                                const char* fluidDensity,
                                const char* fluidViscosity,
                                const char* fluidEnthalpy,
                                const char* fluidHeatCapacity,
                                const char* fluidExpansivity,
                                const char* fluidCompressibility,
                                bool with_steam=false,
                                bool verbose=false );

    /// constructor when only the fluid density and viscosity are required
    IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                const char* temperature,
                                const char* fluidPressure,
                                const char* fluidDensity,
                                const char* fluidViscosity,
                                bool with_steam=false );

    virtual ~IAPWS_H2OPropertiesVisitor();

    /// checks whether the user specified variables have a placement that is compatible with this visitor
    void PlacementChecks();

    /// when all properties are calculated as Node properties
    virtual void Visit( Node<dim>* );
  
    /// when density and viscosity are calculated at the element integration points
    virtual void Visit( Element<dim>* );
  
    /// to ascertain that nothing gets done at the level of the model
    virtual void Visit( Model<dim>* ) {}

    /// optional screen output of calculation results
    void Verbose();

  private:
	  DynamicViscosity   viscosity;
    bool               verbose, with_steam, failed;
    Prop               *props, *lprops, *sprops;   ///< PROST output structures
    ScalarVariable      Tf, Pf, rho, cp, h, x, mu, alpha, beta;
    csmp::Index         P_key,     ///< fluid pressure
                        T_key,     ///< temperature
                        cp_key,    ///< fluid heat capacity
                        h_key,     ///< fluid enthalpy
                        s_key,     ///< saturation steam
                        mu_key,    ///< fluid viscosity
                        alpha_key, ///< fluid expansivity
                        beta_key,  ///< fluid compressibility
                        rho_key;   ///< fluid density
    double            temk, dens;
    const double      dp, dt, kelvin, bar; ///< convergence for PROST         
};


/**
 
@class   IAPWS_H2OPropertiesVisitor IAPWS_H2OPropertiesVisitor "eos\IAPWS_H2OPropertiesVisitor.h"

@author Stephan K. Matthaei
@date 1999

@section motivation Motivation
 
Calculates the PVT properties of pure H2O at the Nodes of a mesh, because this
is also where the dependent variables P,T are located.
 
 
@section applicability Applicability

The visitor calculates the physical variables "nodal fluid density" [kg m-3], 
"nodal relative fluid density" [kg m-3] for a reference pressure, 
"nodal heat capacity liquid" [J kg K-1], "nodal enthalpy liquid" [J kg],
"nodal fluid viscosity" [Pa s-1] and, optionally, also the same properties 
of the gas phase with the names "nodal heat capacity steam", 
"nodal enthalpy steam", as well as the "nodal saturation".
*/

} // end csmp

#endif
