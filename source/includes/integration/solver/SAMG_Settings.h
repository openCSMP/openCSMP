#ifndef SAMG_SETTINGS_H
#define SAMG_SETTINGS_H

#include "CSMP_definitions.h"
#include "SolverSettings.h"

namespace csmp {
  
  class SAMG_Settings : public SolverSettings {
    public:
      SAMG_Settings();
      
      int32 Get_matrix() const;
      int32 Get_ifirst() const;
      double64 Get_eps() const;
      double64 Get_rel_eps() const;
      int32 Get_nsolve() const;
      int32 Get_ncyc() const;
      int32 Get_iswit() const;
      int32 Get_iswtch() const;
      double64 Get_chktol() const;
      int32 Get_idump() const;
      int32 Get_idmp() const;
      int32 Get_iout() const;
      double64 Get_a_cmplx() const;
      double64 Get_g_cmplx() const;
      double64 Get_p_cmplx() const;
      double64 Get_w_avrge() const;
      int32 Get_ncg() const;
      int32 Get_ncyc_done() const;
      int32 Get_ncyc_best() const;
      int32 Get_levelx() const;
      int32 Get_ioform() const;
      int32 Get_ioform_length() const;
      int32 Get_iter_pre() const;
      int* Get_filnam_dump();
      int32 Get_filnam_dump_length() const;
      bool ExplicitSecondary() const;
      bool UsePointBasedApproach() const;
      int32 GetSolverInstance() const;
      int32 Get_mode_mess() const;
      int32 Get_nxtyp() const;
      int32 Get_nrd() const;
      int32 Get_nru() const;
      
      void Set_isym( int32 isym );
      void Set_irow( int32 irow );
      void Set_itypu( int32 itypu );
      void Set_eps( double64 eps );
      void Set_rel_eps( double64 rel_eps );
      void Set_napproach( int32 napproach );
      void Set_nxtyp( int32 nxtyp );
      void Set_nrd( int32 nrd );
      void Set_nru( int32 nrd );
      void Set_internal( int32 internal );
      void Set_nprim( int32 nprim );
      void Set_npr_is_dummy( int32 npr_is_dummy );
      void Set_nint_weights( int32 nint_weights );
      void Set_nint_pat( int32 nint_pat );
      void Set_igam( int32 igam );
      void Set_ncgrad( int32 ncgrad );
      void Set_nkdim( int32 nkdim );
      void Set_ncycle( int32 ncycle );
      void Set_ncyc_done( int32 ncyc_done );
      void Set_ncyc_best( int32 ncyc_best );
      void Set_iswit( int32 iswit );
      void Set_iextent( int32 iextent );
      void Set_ndefault( int32 ndefault );
      void Set_norm_typ( int32 norm_typ );
      void Set_ioscratch( int32 ioscratch );
      void Set_chktol( double64 chktol );
      void Set_idmp( int32 idmp );
      void Set_igdp( int32 igdp );
      void Set_iadp( int32 iadp );
      void Set_iwdp( int32 iwdp );
      void Set_iout1( int32 iout1 );
      void Set_iout2( int32 iout2 );
      void Set_a_cmplx( double64 a_cmplx );
      void Set_g_cmplx( double64 g_cmplx );
      void Set_p_cmplx( double64 p_cmplx );
      void Set_w_avrge( double64 w_avrge );
      void Set_ncgtyp( int32 ncgtyp );
      void Set_nred( int32 nred );
      void Set_nredlev( int32 nredlev );
      void Set_nxf_clean( int32 nxf_clean );
      void Set_neg_diag( int neg_diag );
      void Set_npcol( int32 npcol );
      void Set_levelx( int32 levelx );
      void Set_ioform( std::string ioform );
      void Set_filnam_dump(  const std::string& filnam_dump );
      void Set_iter_pre(int32 iter_pre);
      void Set_mode_mess(int32);
    
      void SetNegative_nsolve( bool negative_nsolve );
      void SetNegative_ncyc( bool negative_ncyc );
      void SetNegative_idump( bool negative_idump );
      void SetNegative_iout( bool negative_iout );
      void ExplicitSecondary( bool explicit_secondary );
      void SetSolverInstance(  int32 instance  );

    private:
      int32       isym_;
      int32       irow_;
      int32       itypu_;
      double64    eps_;
      double64    rel_eps_;
      int32       napproach_;
      int32       nxtyp_;
      int32       nrd_;
      int32       nru_;
      int32       internal_;
      int32       nprim_;
      int32       npr_is_dummy_;
      int32       nint_weights_;
      int32       nint_pat_;
      int32       igam_;
      int32       ncgrad_;
      int32       nkdim_;
      int32       ncycle_;
      int32       ncyc_done_;    /// < total number of cycles (iterations) performed
      int32       ncyc_best_;    /// < stores lowest number of cycles achieved
      int32       iswit_;
      int32       iextent_;
      int32       ndefault_;
      int32       norm_typ_;
      int32       ioscratch_;
      double64    chktol_;
      int32       idmp_;
      int32       igdp_;
      int32       iadp_;
      int32       iwdp_;
      int32       iout1_;
      int32       iout2_;
      double64    a_cmplx_;
      double64    g_cmplx_;
      double64    p_cmplx_;
      double64    w_avrge_;
      int32       ncgtyp_;
      int32       ioform_;
      int32       ioform_length_;
      int32       levelx_;
      std::string filnam_dump_;
      int32       filnam_dump_length_;
      int         filnam_dump_Array_[50];
      int         neg_diag_;
      int32       nred_;
      int32       nredlev_;
      int32       nxf_clean_;
      int32       npcol_;
      int32       solver_instance_;
      int32       iter_pre_;
      int32       mode_mess_;

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
