// Copyright © 2019 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LAYERED_COMPOSITE_PROCESSOR1_H
#define LAYERED_COMPOSITE_PROCESSOR1_H

#include "CSMP_definitions.h"

namespace csmp {

/**
   Implementation of relative permeability model of Maartje Boon
   for layered rocktypes = composities.
   The values are supplied by parameter table that initialises model for rocktype.
   Documentation in Maple worksheet 'Maartje_relperms'.
 
   Calculation procedure
 
   0. for given rocktype, table is initialised from file
   1. Average saturation values, sw are calculated for each row for the given capillary number
   2. The relative permeability is calculated by linear interpolation between the rocks of interest, for the wetting and the nonwetting phase
 
   This model is meant to be extended to the case where horizontal relperms differ from vertical ones.
 
   The input table is supposed to have the following format (# comment lines to be ignored).
 
   @code
   HARP_FSst_Slt_new_names.txt - Input data table to find Heterogeneity Aware Relative Permeability (HARP ?) as a function of sw and Nc = k grad P / sigma, sigma = interfacial tension H2O-CO2
    # comment line
    2 # total number of rocktypes in file
    1 # rocktype identifier 1: FSst-Slt (Fine Sandstone - Silt) Planar Bedding
    30 6 # table size(rows vs. columns)
    # F-flow    Sw_VL (cell average)  Sw_CL(cell average)  Sw_low_CL  Sw_high_CL    Sfactor
    0.990000000000000  0.958327612404074  0.978851642918739  1  0.957703285837478  137051613.363363
    0.956206896551724  0.935003844772692  0.964807268982292  1  0.929614537964584  28739845.1229135
    0.922413793103448  0.922536293954733  0.959418405274487  1  0.918836810548975  43986865.1983833
    @endcode
 
    // saturations
    double Sw_VL_, Sw_CL_;          ///< cell saturations at viscous and capillary limit
    double Sw_low_CL_, Sw_high_CL_; ///< average saturations at CL in the different layers a
    double Sw_at_Nc_;               ///< saturation at the capillary number of interest
    double Sfactor_;                ///< scaling factor
*/
class LayeredCompositeProcessor1 {
  public:
    explicit LayeredCompositeProcessor1( const char* RRT_data );
    void ReadRockTypeData( const char* datafile );
  
    // KEY METHODS FOR THE USER - non-constant as they modify the RRT table
  
    double krw( double sw, double Nc, int rocktype );
    double krn( double sw, double Nc, int rocktype );
  
    /// computes wetting- and non-wetting phase relative permeabilities in one operation
    std::pair<double,double> RelativePermeability( double sw, double Nc, int rocktype );
  
    double pc( double sw, int rocktype ) const;
  
    /// writes textfile with sw, krw(sw,Nc), krn(sw,Nc) values computed for rocktype in 0.05 saturation increments
    void WriteRelativePermeabilityTable( const char* filename, int rocktype, double Nc );
  
    /// reports the stored rocktype table(s) etc.
    void Out() const;

  private:
  
    // permeability averages
    /// vertical thickness weighted (harmonic) mean of the vertical layer permeabilities
    double k_AverageY() const;
    /// horizontal thickness weighted average of the horizontal layer permeabilities
    double k_AverageX() const;
    /// water saturation in cell at given capillary number; calculated from entries in table
    double SwAtNc( double Nc, double Sw_VL, double Sw_CL, double Sfactor ) const;
    /// corresponding sw in low-k laminations
    double SwAtNc_Low_k_Layer( double Nc, double Sw_VL, double Sw_low_CL, double Sfactor ) const;
    /// corresponding sw in high-k laminations
    double SwAtNc_High_k_Layer( double Nc, double Sw_VL, double Sw_high_CL, double Sfactor ) const;
    /// effective sw in low-k laminations
    double SwStarLow( double Sw_low_at_Nc ) const;
    /// effective sw in low-k laminations
    double SwStarHigh( double Sw_high_at_Nc ) const;
  
    /// krw relperm of low-k layer at given saturation and capillary number
    double KrwLow( double SwStar_low ) const;
    /// krw relperm of high-k layer at given saturation and capillary number
    double KrwHigh( double SwStar_high ) const;

    /// krn relperm of low-k layer at given saturation and capillary number
    double KrnLow( double SwStar_low ) const;
    /// krn relperm of high-k layer at given saturation and capillary number
    double KrnHigh( double SwStar_high ) const;
  
    /// relative permeability at the given water saturation and capillary number
    double KrwComposite( double swAtNc_Low_k_Layer, double swAtNc_High_k_Layer ) const;
    double KrnComposite( double swAtNc_Low_k_Layer, double swAtNc_High_k_Layer ) const;
  
  private: // parameters & data
  
    enum { f, Sw_VL, Sw_CL, Sw_CL_low_k, Sw_CL_high_k, Sfactor };

    // material properties - should be part of rocktype
    double k_low_, k_high_;         ///< layer permeabilities
    double phi_low_, phi_high_;     ///< porosities of low and high k layers
    double LY_low_, LY_high_;       ///< cumulative layer thickness in the vertical direction (Y), normalized summing to 1
    double Swi_low_, Swi_high_;     ///< irreducible saturations of the 2 different layers
    double m_low_, m_high_;         ///< van Genuchten exponents for the 2 different layers
    double pd_low_, pd_high_;       ///< capillary (drainage) entry pressure of low and high k layers
    double bcp_low_, bcp_high_;     ///< Brooks-Corey 64' exponents for low and high k layers
  
    std::vector<std::vector<std::vector<double> > > RRT_;  ///< data table for n-rocktypes with extra 3 columns for sw, krw, krn
};

} // end csmp

#endif /* LAYERED_COMPOSITE_PROCESSOR1_H */
