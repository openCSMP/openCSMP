// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

// NaClH2OPropertiesVisitorPHX_initializer_list.h
//
// Constructor initializer list.  Most members use in-class defaults from the
// header; only reference-bound members and objects with constructor arguments
// appear here.  Members are initialized in declaration order regardless of
// the order listed below.

: pref( model.Database() ),
  pmesh( model.Mesh() ),
  csmp_error( ErrorHandler::Instance() ),

  // ── Equilibrator: takes references to working variables ─────────────────
  // These doubles must be declared (and default-initialized) before
  // equilibrator in the header — guaranteed by declaration order.
  equilibrator(
      m_rock_,
      cp_rock_,
      rho_rock_,
      m_air_,
      cp_air_,
      phi_,
      m_fluid_,
      wt_,
      tp_,
      p_current_,
      H_current_,
      H_previous_,
      fixed_temperature,
      t_fixed,
      t_diffusion_,
      verbose_equilibration ),

  // ── Fluid EOS object: takes references for T, P, X, H, cp_rock, rho_rock, phi
  fluid( t_, p_, x_, h_fluid_, cp_rock_, rho_rock_, phi_, false )
