#ifndef SAMG_SETTINGS_H
#define SAMG_SETTINGS_H

#include "CSMP_definitions.h"
#include "SolverSettings.h"

namespace csmp {
  
  class SAMG_Settings : public SolverSettings {
    public:
      SAMG_Settings();
      
      int32_t Get_matrix() const;
      int32_t Get_ifirst() const;
      double Get_eps() const;
      double Get_rel_eps() const;
      int32_t Get_nsolve() const;
      int32_t Get_ncyc() const;
      int32_t Get_iswit() const;
      int32_t Get_iswtch() const;
      double Get_chktol() const;
      int32_t Get_idump() const;
      int32_t Get_idmp() const;
      int32_t Get_iout() const;
      double Get_a_cmplx() const;
      double Get_g_cmplx() const;
      double Get_p_cmplx() const;
      double Get_w_avrge() const;
      int32_t Get_ncg() const;
      int32_t Get_ncyc_done() const;
      int32_t Get_ncyc_best() const;
      int32_t Get_levelx() const;
      int32_t Get_ioform() const;
      int32_t Get_ioform_length() const;
      int32_t Get_iter_pre() const;
      int* Get_filnam_dump();
      int32_t Get_filnam_dump_length() const;
      bool ExplicitSecondary() const;
      bool UsePointBasedApproach() const;
      int32_t GetSolverInstance() const;
      int32_t Get_mode_mess() const;
      int32_t Get_nxtyp() const;
      int32_t Get_nrd() const;
      int32_t Get_nru() const;
      
      void Set_isym( int32_t isym );
      void Set_irow( int32_t irow );
      void Set_itypu( int32_t itypu );
      void Set_eps( double eps );
      void Set_rel_eps( double rel_eps );
      void Set_napproach( int32_t napproach );
      void Set_nxtyp( int32_t nxtyp );
      void Set_nrd( int32_t nrd );
      void Set_nru( int32_t nrd );
      void Set_internal( int32_t internal );
      void Set_nprim( int32_t nprim );
      void Set_npr_is_dummy( int32_t npr_is_dummy );
      void Set_nint_weights( int32_t nint_weights );
      void Set_nint_pat( int32_t nint_pat );
      void Set_igam( int32_t igam );
      void Set_ncgrad( int32_t ncgrad );
      void Set_nkdim( int32_t nkdim );
      void Set_ncycle( int32_t ncycle );
      void Set_ncyc_done( int32_t ncyc_done );
      void Set_ncyc_best( int32_t ncyc_best );
      void Set_iswit( int32_t iswit );
      void Set_iextent( int32_t iextent );
      void Set_ndefault( int32_t ndefault );
      void Set_norm_typ( int32_t norm_typ );
      void Set_ioscratch( int32_t ioscratch );
      void Set_chktol( double chktol );
      void Set_idmp( int32_t idmp );
      void Set_igdp( int32_t igdp );
      void Set_iadp( int32_t iadp );
      void Set_iwdp( int32_t iwdp );
      void Set_iout1( int32_t iout1 );
      void Set_iout2( int32_t iout2 );
      void Set_a_cmplx( double a_cmplx );
      void Set_g_cmplx( double g_cmplx );
      void Set_p_cmplx( double p_cmplx );
      void Set_w_avrge( double w_avrge );
      void Set_ncgtyp( int32_t ncgtyp );
      void Set_nred( int32_t nred );
      void Set_nredlev( int32_t nredlev );
      void Set_nxf_clean( int32_t nxf_clean );
      void Set_neg_diag( int neg_diag );
      void Set_npcol( int32_t npcol );
      void Set_levelx( int32_t levelx );
      void Set_ioform( std::string ioform );
      void Set_filnam_dump(  const std::string& filnam_dump );
      void Set_iter_pre(int32_t iter_pre);
      void Set_mode_mess(int32_t);
    
      void SetNegative_nsolve( bool negative_nsolve );
      void SetNegative_ncyc( bool negative_ncyc );
      void SetNegative_idump( bool negative_idump );
      void SetNegative_iout( bool negative_iout );
      /// if this is not set 'true' all attempts to set secondary parameters will fail (default=false)
      void ExplicitSecondary( bool explicit_secondary );
      void SetSolverInstance(  int32_t instance  );

    private:
      int32_t       isym_;
      int32_t       irow_;
      int32_t       itypu_;
      double    eps_;
      double    rel_eps_;
      int32_t       napproach_;
      int32_t       nxtyp_;
      int32_t       nrd_;
      int32_t       nru_;
      int32_t       internal_;
      int32_t       nprim_;
      int32_t       npr_is_dummy_;
      int32_t       nint_weights_;
      int32_t       nint_pat_;
      int32_t       igam_;
      int32_t       ncgrad_;
      int32_t       nkdim_;
      int32_t       ncycle_;
      int32_t       ncyc_done_;    /// < total number of cycles (iterations) performed
      int32_t       ncyc_best_;    /// < stores lowest number of cycles achieved
      int32_t       iswit_;
      int32_t       iextent_;
      int32_t       ndefault_;
      int32_t       norm_typ_;
      int32_t       ioscratch_;
      double    chktol_;
      int32_t       idmp_;
      int32_t       igdp_;
      int32_t       iadp_;
      int32_t       iwdp_;
      int32_t       iout1_;
      int32_t       iout2_;
      double    a_cmplx_;
      double    g_cmplx_;
      double    p_cmplx_;
      double    w_avrge_;
      int32_t       ncgtyp_;
      int32_t       ioform_;
      int32_t       ioform_length_;
      int32_t       levelx_;
      std::string filnam_dump_;
      int32_t       filnam_dump_length_;
      int         filnam_dump_Array_[50];
      int         neg_diag_;
      int32_t       nred_;
      int32_t       nredlev_;
      int32_t       nxf_clean_;
      int32_t       npcol_;
      int32_t       solver_instance_;
      int32_t       iter_pre_;
      int32_t       mode_mess_;

      bool negative_nsolve_;
      bool negative_ncyc_;
      bool negative_idump_;
      bool negative_iout_;
      bool explicit_secondary_;
  };

/**
 
@class SAMG_Settings SAMG_Settings "solver/SAMG_Settings.h"
@author S.K. Matthaei
@author S. Geiger
@author G. Roberts
@date 2001
 
@section motivation Motivation
 
Factoring out SAMG tuning parameter in a separate class so that they
can be changed indepently of the solver object.
 
 
@section design Design Intent
 
SAMG_Settings objects store SAMG tuning parameter. These tuning parameter
are divided into switches and their subswitches.

Each switch can be accessed by its respective Get-method and its subswitches
can be set using the respective Set-methods.
 
Every Set-method checks if the input data conforms to the restrictions
imposed by the SAMG solver. If the input value is invald, a domain_error
exception is thrown.

@warning Only the rules for the subswitch itself are checked! Interdependencies
between subswitches are NOT checked! The user is responsible to guarantee that
all subswitches together form a meaniful parametrization for the SAMG solver.
*/
  
} // end namespace csp

#endif
