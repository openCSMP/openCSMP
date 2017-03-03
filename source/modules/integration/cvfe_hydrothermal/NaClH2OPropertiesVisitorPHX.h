#ifndef NACLH2OPROPERTIESVISITORPHX_H
#define NACLH2OPROPERTIESVISITORPHX_H

#include "CSMP_definitions.h"
#include "Model.h"
#include "Visitor.h"
#include "ErrorHandler.h"

#include "H2ONaClFluidProperties.h"
#include "H2ONaClThermalEquilibrator.h"

#include "Rock.h"

namespace csmp
{

  template<size_t dim> 
    class NaClH2OPropertiesVisitorPHX : public Visitor<dim> 
    {
    public:

      NaClH2OPropertiesVisitorPHX( Model<dim>& );
      ~NaClH2OPropertiesVisitorPHX();
      
      virtual void Visit(Node<dim>* n);  
	  virtual void Visit(Region<dim>* n);

      void SetTimeIncrement( double64 time_increment );
      void InitialPropertiesFromPTX();
      void TopBoundaryHandling( bool TB );
      
      void TemperatureDependentHeatCapacityRock( double64 cpr_min_ext, double64 t_min_ext,
                                                 double64 cpr_max_ext, double64 t_max_ext );
      void WithOpenBoundaries( double64 reference_specific_enthalpy, double64 reference_salinity );
      void WithOpenBoundaries( double64 reference_salinity );

      void SetAdjustCompressibilityAfterPhaseChangeBoolTo( bool adjust ); // default true, set to false if you want to switch off
      void SetMaximumCompressibilityCyclesAfterPhaseChange( double64  max );    // number of timesteps after phase change where this should be applied
      // potentially problematic if yet another phase change within this number of timesteps!

    private:
      // Accessing CSMP memory
      PropertyDatabase<dim>&     pref;
      /* csmp::MemoryManager<dim>&  pmem; */
      csmp::MeshManager<dim>&    pmesh;

      bool                         top_boundary;
      bool                         bogus_variables;
      bool                         first;
      bool                         pure_halite;
      bool                         fixed_temperature;
      bool                         T_dependent_cpr;
      bool                         open_boundaries;
      bool                         const_top_pressure;
      bool                         verbose_eq;

      bool                         adjust_compressibility_after_phasechange;
      double64                     max_after_phasechange_counter;

      // Variables for the "FluidRockEquilibrator" object
      double64                     m_rock_; 
      double64                     cp_rock_; 
      double64                     phi_; 
      double64                     m_fluid_; 
      double64                     tp_; 
      double64                     p_current_; 
      double64                     H_current_; 
      double64                     MS_current_;
      double64                     H_previous_; 
      double64                     t_fixed; 
      double64                     rho_rock_; 
      double64                     t_diffusion_;
      double64                     dT_diff_; 
      double64                     dh_rock_diff_; 
      double64                     dh_fluid_diff_; 
      double64                     dhCl_; 
      double64                     dhCv_;
      double64                     dxCl_; 
      double64                     dxCv_; 
      double64                     dx_diff_; 
      double64                     x_;
      double64                     dmt_src; 
      double64                     dh_src; 
      double64                     ds_src; 
      double64                     dt_; 
      double64                     dml_; 
      double64                     dmv_; 
      double64                     h_fluid_; 
      double64                     wt_;
      double64                     t_; 
      double64                     p_bar_; 
      double64                     previous_state;
      double64                     mlv_;
      double64                     volume_factor; 
      double64                     volume_factor_LHS; 
      double64                     volume_factor_RHS;
      double64                     pore_volume;
      double64                     rock_volume;
      double64                     kelvin;
      double64                     res_sl;
      double64                     res_sv;
      double64                     krl;
      double64                     krv;
      double64                     ref_spec_h; // reference specific enthalpy (only for constant p at top)
      double64                     ref_sal; // reference salinity at open top
      double64                     time_factor_p;
      // next two unused for the moment
      double64                     time_factor_h;
      double64                     time_factor_s;
      double64                     expected_dp;
      double64                     hrock_prev;
      double64                     hrock_curr;

      int                          old_state;
      int                          current_state;

      // primary variables
      ScalarVariable     t;
      ScalarVariable     tp;
      ScalarVariable     p;

      // Saturation variables
      ScalarVariable     sl;
      ScalarVariable     sv;
      ScalarVariable     slp;
      ScalarVariable     svp;
      ScalarVariable     sh;
      
      // densities
      ScalarVariable     mt;
      ScalarVariable     mtp;
      ScalarVariable     ml;
      ScalarVariable     mv;
      ScalarVariable     mlp;
      ScalarVariable     mvp;
      ScalarVariable     rl;
      ScalarVariable     rv;
      ScalarVariable     rlp;
      ScalarVariable     rvp;
      ScalarVariable     mf;
      ScalarVariable     mh;
      ScalarVariable     rh;

      // viscosities
      ScalarVariable     mul;
      ScalarVariable     muv;
      ScalarVariable     mulp;
      ScalarVariable     muvp;

      // enthalpy variables
      ScalarVariable     hCl;
      ScalarVariable     hCv;
      ScalarVariable     hClp;
      ScalarVariable     hCvp;
      ScalarVariable     Htp;
      ScalarVariable     hf;
      ScalarVariable     hl;
      ScalarVariable     hv;
      ScalarVariable     cpf;
      ScalarVariable     hVl;
      ScalarVariable     hVv;
      ScalarVariable     hh;
      ScalarVariable     hVh;
      ScalarVariable     hCh;

      // compressibility variables
      ScalarVariable     beta;
      ScalarVariable     beta_p;
      ScalarVariable     beta_ref;
      ScalarVariable     nQ;
      ScalarVariable     CT;

      // rock variables
      ScalarVariable     cpr;
      ScalarVariable     phi;
      ScalarVariable     rr;
      ScalarVariable     cpv;
      ScalarVariable     beta_rock;
      ScalarVariable     ncp;

      // fluid source variables
      ScalarVariable     src_h;
      ScalarVariable     src_rate;
      ScalarVariable     src_wt;
      
      ScalarVariable     rl_transport;
      ScalarVariable     rv_transport;
      ScalarVariable     vol_fac;
      ScalarVariable     rho_bulk;
      
      // salinity variables
      ScalarVariable     wt;
      ScalarVariable     xf;
      ScalarVariable     xl;
      ScalarVariable     xVl;
      ScalarVariable     xCl;
      ScalarVariable     xClp;
      ScalarVariable     xv;
      ScalarVariable     xVv;
      ScalarVariable     xCv;
      ScalarVariable     xCvp;
      ScalarVariable     xh;
      ScalarVariable     xVh;
      ScalarVariable     xCh;
      ScalarVariable     xCf;
      ScalarVariable     xCfp;
      ScalarVariable     ms;
      ScalarVariable     msp;
      ScalarVariable     xlp;
      ScalarVariable     xvp;
      ScalarVariable     hlp;
      ScalarVariable     hvp;
      ScalarVariable     hhp;
      
      // to store fluid state
      ScalarVariable     state;
      ScalarVariable     state_p;
      ScalarVariable     visited; 

      // to monitor phase changes that may lead to p < atm; TD May 2011
      ScalarVariable     dangerous_phase_change;
      ScalarVariable     after_phasechange_counter;

      // new JPW
      ScalarVariable     mml;
      ScalarVariable     mmv;
      ScalarVariable     mmld;
      ScalarVariable     mmvd;
      ScalarVariable     eml;
      ScalarVariable     emv;
      ScalarVariable     emld;
      ScalarVariable     emvd;
      ScalarVariable     xml;
      ScalarVariable     xmv;
      ScalarVariable     pv;
      ScalarVariable     rvl;
      ScalarVariable     rvv; 
      ScalarVariable     bfm;
      ScalarVariable     bfe;
      ScalarVariable     bfs; 

      ScalarVariable     time_factor;
      ScalarVariable     reference_enthalpy_top;
      
      // Index to variables
      // primary variables
      csmp::Index                  p_key;     // fluid pressure
      csmp::Index                  t_key;     // temperature
      csmp::Index                  tp_key;    // temperature
	
      // Saturation variables
      csmp::Index                  sl_key;    // saturation liquid
      csmp::Index                  sv_key;    // saturation vapor
      csmp::Index                  sh_key;    // saturation halite
      // fluid state
      csmp::Index                  state_key;
      //*** new TD May 2011
      csmp::Index                  state_p_key;

      // Phase change tracking
      csmp::Index                  apc_key;
      csmp::Index                  dpc_key;
      //*** end new
      // visitation
      csmp::Index                  visit_key; //added JPW 27.10.2010
	
      // densities
      csmp::Index                  mt_key;    // density (mixture)
      csmp::Index                  mtp_key;   // previous density (mixture)
      csmp::Index                  rl_key;    // density liquid
      csmp::Index                  rv_key;    // density vapor
      csmp::Index                  rh_key;    // density halite
      csmp::Index                  mf_key;    // fluid mass
      csmp::Index                  ml_key;    // fluid mass liquid
      csmp::Index                  mv_key;    // fluid mass vapor
      csmp::Index                  mlp_key;   // previous fluid mass liquid
      csmp::Index                  mvp_key;   // previous fluid mass vapor
      csmp::Index                  mh_key;    // solid mass halite
      //extras
      csmp::Index                  rl_transport_key;    // density liquid transport
      csmp::Index                  rv_transport_key;    // density vapor transport
      csmp::Index                  vol_fac_key;
      csmp::Index                  rho_bulk_key;
      csmp::Index                  nQ_key;    // nodal fluid volume source
      // viscosities
      csmp::Index                  mul_key;   // viscosity halite
      csmp::Index                  muv_key;   // viscosity vapor
      // enthalpy variales
      csmp::Index                  hf_key;    // enthalpy (mixture)
      csmp::Index                  hl_key;    // enthalpy liquid
      csmp::Index                  hv_key;    // enthalpy vapor
      csmp::Index                  hh_key;    // enthalpy halite
      csmp::Index                  hCl_key;   // enthalpy content liquid (S * rho * h)
      csmp::Index                  hClp_key;  // previous enthalpy content liquid (S * rho * h)
      csmp::Index                  hCv_key;   // enthalpy content vapor (S * rho * h)
      csmp::Index                  hCvp_key;  // previous enthalpy content vapor (S * rho * h)
      csmp::Index                  hCh_key;   // enthalpy content halite (S * rho * h)
      csmp::Index                  hVl_key;   // volumetric fluid enthalpy liquid (rho * h)
      csmp::Index                  hVv_key;   // volumetric fluid enthalpy vapor (rho * h)
      csmp::Index                  hVh_key;   // volumetric fluid enthalpy halite (rho * h)
      csmp::Index                  cpf_key;   // heat capacity of bulk fluid
      csmp::Index                  wt_key;    // salinity
      csmp::Index                  xf_key;    // mass fraction salt fluid = salinity
      csmp::Index                  xl_key;    // mass fraction salt liquid
      csmp::Index                  xv_key;    // mass fraction salt vapor
      csmp::Index                  xh_key;    // mass fraction salt halite
      csmp::Index                  xVl_key;   // volumetric salinityliquid (rho * x)
      csmp::Index                  xVv_key;   // volumetric salinity vapor (rho * x)
      csmp::Index                  xVh_key;   // volumetric salinity halite (rho * x)
      csmp::Index                  xCl_key;   // salt content liquid (S * rho * x)
      csmp::Index                  xCv_key;   // salt content vapor (S * rho * x)
      csmp::Index                  xCh_key;   // salt content halite (S * rho * x)
      csmp::Index                  xClp_key;  // previous salt content liquid (S * rho * x)
      csmp::Index                  xCvp_key;  // previous salt content vapor (S * rho * x)
      csmp::Index                  xCf_key;   // salt content fluid (rhof * xf)
      //xCfp_key;  // previous salt content fluid (rhof * xf)
      csmp::Index                  msp_key;   // previous mass salt
	
      //                       Ht_key;    // total energy in fluid + rock in J
      csmp::Index                  Htp_key;   // previous total energy fluid + rock in J
	
      // expansivities
      csmp::Index                  beta_key;  // fluid compressibility
      //*** new TD May 2011
      csmp::Index                  beta_p_key;  // previous fluid compressibility
      csmp::Index                  beta_ref_key; // reference compressibility when crossing phase boundaries
      //*** end new
      csmp::Index                  CT_key;  // nodal total compressibility
      // rock variables
      csmp::Index                  phi_key;   // nodal porosity
      csmp::Index                  cpv_key;   // control pore volume
      csmp::Index                  cpr_key;   // heat capacity rock
      csmp::Index                  ncp_key;   // nodal heat capacity (added JPW)
      csmp::Index                  rr_key;    // density rock
      csmp::Index                  betar_key;    // nodal compressibility rock
      // velocities
      // fluid source     
      csmp::Index                  src_h_key;
      csmp::Index                  src_rate_key;
      csmp::Index                  src_wt_key;
	
      // new JPW
      csmp::Index                  mml_key; // liquid mass mobility
      csmp::Index                  mmv_key; // vapor mass mobility
      csmp::Index                  mmld_key; // liquid mass mobility density
      csmp::Index                  mmvd_key; // vapor mass mobility density
      csmp::Index                  eml_key; // liquid enthalpy mobility
      csmp::Index                  emv_key; // vapor enthalpy mobility
      csmp::Index                  emld_key; // liquid enthalpy mobility density
      csmp::Index                  emvd_key; // vapor enthalpy mobility density
      csmp::Index                  xml_key; // liquid salt mobility
      csmp::Index                  xmv_key; // vapor salt mobility
      csmp::Index                  rvl_key; // relperm viscosity liquid
      csmp::Index                  rvv_key; // relperm viscosity vapor

      // JPW Nov 2010
      csmp::Index                  bfm_key; // boundary flow mass
      csmp::Index                  bfe_key; // boundary flow enthalpy
      csmp::Index                  bfs_key; // boundary flow salinity
      csmp::Index                  time_factor_key;
      csmp::Index                  ref_h_top_key;
      
      // Sowat objects
      Rock                         rock;
      H2ONaClThermalEquilibrator   equilibrator;
      Fluidproperties              Liquid;
      Fluidproperties              Vapor;
      Fluidproperties              Bulk;
      Fluidproperties              Salt;

      // Private Memberfunctions
      double64                     EffectiveLiquidSaturationHalitePresent( double64 sat_liquid, double64 sat_vapor ) const;
      double64                     RelativePermeabilityLiquid( double64 sat_liquid, double64 sat_vapor ) const;
      void                         ReadAllVariables( Node<dim>* n );

      void                         CheckForOutOfRange( Node<dim>* n );
      void                         CheckBoundaryFlags( Node<dim>* n );

      void                         CalculateAbsoluteVariables( );
      void                         UpdateSowatVariables();
      void                         Equilibrate( Node<dim>* n );
      void                         ScreenOutputSowatVariables();
      void                         UpdateCSMPVariables( Node<dim>& n );

      double64                           TwoPhasePureWaterCompressibility(double64 cpl, double64 cpv);

      void                         Output_PTXState();

      void                         StorePropertiesAndFlags( Node<dim>& n );  
      void                         StoreInitialPropertiesAndFlags( Node<dim>& n );  

      void                         BoundaryHandling( Node<dim>* n );
      double64                           BoundaryFlow( );
      void                         BoundaryIteration( );

      void                         CheckVolumeMismatchCompensation();
      void                         CheckPhaseChange();

      void                         TimeStepAdjustment();
      void                         VolumeFactorComputations( Node<dim>* n );
      void                         PrepareVariablesForStorage();

      ErrorHandler&                csmp_error;
      //      void BoundaryIteration( Node<dim>* n );
      //      void BoundaryIteration( Node<dim>& n );
      //	                         void UpdatePorosity();
      //	double64  CalculateTimeMultiplicationFactor( Node<dim>* n );
    };
  
  
  /*!
    <!-- ================================================================== -->
    <H1>NaClH2OPropertiesVisitorPHX</H1> <!-- ======================================= -->
    <!-- ================================================================== -->

    <!-- ================================================================== -->
    <H3>MOTIVATION</H3>
    <!-- ================================================================== -->
    In simulations of multiphase fluid flow involving heat and mass transport,
    fluid properties are typically stored at nodes, and heat and mass transport
    via diffusion and advection is done in a decoupled manner. This ultimately
    means that after the diffusion and advection steps, the fluid phase state
    may have changed and, even worse, very likely means that fluid and rock
    are not in thermal equilibrium. NaClH2OPropertiesVisitorPHX will
    perform a thermal equilibration between fluid and rock (using class
    H2ONaClThermalEquilibrator) at a given pressure, total
    enthalpy and salinity and update fluid properties for the equilibrated 
    state.

    <p>
    <!-- ================================================================== -->
    <H3>DESIGN INTENT</H3>
    <!-- ================================================================== -->
    NaClH2OPropertiesVisitorPHX is an FEM_Visitor within the CSMP framework,
    see also documentation there. It applies to nodes only.

    <p>
    <!-- ================================================================== -->
    <H3>APPLICABILITY</H3>
    <!-- ================================================================== -->
    NaClH2OPropertiesVisitorPHX can be applied in simulations with pure H2O
    fluids in combination with heat and mass transport. Simpler simulations 
    should use other tools in the CSMP++ package. A fundamental assumption is
    that of thermal equilibrium between fluid and rock is requested, otherwise
    using NaClH2OPropertiesVisitorPHX makes no sense at all.

    NaClH2OPropertiesVisitorPHX also assumes that fluid properties are stored
    on the nodes (i.e., it relates to a node-centered finite volume schemes for
    solving the advection equations) and that rock properties discretized on the
    element are uniquely defined such that their interpolated nodal values can 
    meaningfully be used (*** is that fully correct? ***). 

    <p>
    <!-- ================================================================== -->
    <H3>STRUCTURE</H3>
    <!-- ================================================================== -->
    <p>
    <!-- ================================================================== -->
    <H3>PARTICIPANTS</H3>
    <!-- ================================================================== -->

    NaClH2OPropertiesVisitorPHX uses H2ONaClThermalEquilibrator,
    which itself calls a number of classes related to fluid properties, in
    particular FluidLookup and the specialized classes for phases in the H2O
    system (see documentation for FluidLookup etc.). 

    <p>
    <!-- ================================================================== -->
    <H3>COLLABORATIONS</H3>
    <!-- ================================================================== -->

    <p>
    <!-- ================================================================== -->
    <H3>CONSEQUENCES</H3>
    <!-- ================================================================== -->

    During 
    <p>
    <!-- ================================================================== -->
    <H3>IMPLEMENTATION</H3>
    <!-- ================================================================== -->
    <p>
    <!-- ================================================================== -->
    <H3>APPLICATION EXAMPLES</H3>
    <!-- ================================================================== -->
    <p>
    <!-- ================================================================== -->
    copyright (c) 2004-2009 by Drs. Thomas Driesner, Sebastian Geiger, Dim Coumou 
  */

  // NaClH2OPropertiesVisitorPHX Methods ====================

  template<size_t dim>
    inline double64  NaClH2OPropertiesVisitorPHX<dim>::EffectiveLiquidSaturationHalitePresent( double64 sat_liquid, double64 sat_vapor ) const
    {
      double64 seff;
      double64 my_res_sl(res_sl*(sat_liquid+sat_vapor)),
        my_res_sv(res_sv*(sat_liquid+sat_vapor));

      // Changed JPW 26.5.2011
      //      seff = (sat_liquid - res_sl) / (sat_liquid + sat_vapor - res_sl - res_sv);
      //      seff = (sat_liquid - my_res_sl) / (sat_liquid + sat_vapor - my_res_sl - my_res_sv);
      seff = (sat_liquid - my_res_sl) / (sat_liquid + sat_vapor - my_res_sl - my_res_sv);
      if      ( seff > static_cast<double64>(1.0) ) return 1.0;
      else if ( seff < static_cast<double64>(0.0) ) return 0.0;
      else                                    return seff;

    } // end EffectiveSaturation

  template<size_t dim>
    inline double64  NaClH2OPropertiesVisitorPHX<dim>::RelativePermeabilityLiquid( double64 sat_liquid, double64 sat_vapor ) const
    {
      double64 seff, sat_total(sat_liquid+sat_vapor);
      double64 my_res_sl(res_sl*sat_total),
        my_res_sv(res_sv*sat_total);

      // Changed JPW 26.5.2011
      seff = sat_total*(sat_liquid - my_res_sl) / (sat_liquid + sat_vapor - my_res_sl - my_res_sv);
      if      ( seff > static_cast<double64>(sat_total) ) return sat_total;
      else if ( seff < static_cast<double64>(0.0) ) return 0.0;
      else                                    return seff;

    } // end RelativePermeabilityLiquid

  //*** new TD May 2011
  template<size_t dim>
    inline void NaClH2OPropertiesVisitorPHX<dim>::SetMaximumCompressibilityCyclesAfterPhaseChange( double64 max )
    {
      max_after_phasechange_counter = max;
    }

  template<size_t dim>
    inline void NaClH2OPropertiesVisitorPHX<dim>::SetAdjustCompressibilityAfterPhaseChangeBoolTo( bool adjust )
    {
      adjust_compressibility_after_phasechange = adjust;
    }
  //*** end new
}

#endif

