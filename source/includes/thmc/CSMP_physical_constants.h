// Copyright (c) 2012 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_PHYSICAL_CONSTANTS_H
#define CSMP_PHYSICAL_CONSTANTS_H

#include "CSMP_definitions.h"

namespace csmp {

constexpr double CSMP_PI(3.14159265358979);      ///<  pi (m) perimeter of circle with a diameter of 1

/** common constants in SI units, alphabetical order, notation (units in brackets)
*/

constexpr double CSMP_NO_DATA_VALUE(-1.0e+30);   ///<  any undefined value of a material property or physical variable

constexpr double ACC_GRAVITY(9.80655);           ///<  g (m/s2), acceleration due to gravity (32.2 ft/s2), extra decimal places require regional information
constexpr double ATMOSPHERIC_PRESSURE(101325.);  ///<  standard p_atm (Pa)
constexpr double AVOGADRO(6.0221415e23);         ///<  N_A (mol-1)
constexpr double BOHR_RADIUS(0.529177249e-10);   ///<  a0 (m)
constexpr double BOLTZMANN(1.3806505e-23);       ///<  R (J/K)
constexpr double COULOMB(8.987552e9);            ///<  K (1/4 pi e0), N m2 C-2
constexpr double FARADAY(96485.3383);            ///<  F (C/mol)
constexpr double G_GRAVITY(6.67259e-11);         ///<  G (m3 kg-1 s-2)
constexpr double HEAT_FLOW_UNIT(0.04184);        ///<  HFU, W m-2
constexpr double IDEAL_GAS_VOLUME(22.413996e-3); ///<  V_m=RT/p (m3/mol)
constexpr double LOSCHMIDT(2.6867773e25);        ///<  m-3 (0oC, 1 bar)
constexpr double LIGHT_SPEED(299792458.);        ///<  m/s (in vacuum)
constexpr double MASS_OF_ATOM(1.66053886e-27);   ///<  m_U (kg)
constexpr double MOLAR_GAS_CONSTANT(8.31451);    ///<  R (J/ mol.K) = m2·kg/s2·K·mol
constexpr double MOLAR_PLANCK(3.990312716e-10);  ///<  N_A h (J.s / mol)
constexpr double RYDBERG(10973731.534);          ///<  R_inf (m-1)
constexpr double STEFAN_BOLTZMANN(5.6704e-8);    ///<  sigma (W m-2 K-4)
constexpr double WIEN_DISPLACEMENT(2.897756e-3); ///<  b (m . K)

} // end csmp

#endif
