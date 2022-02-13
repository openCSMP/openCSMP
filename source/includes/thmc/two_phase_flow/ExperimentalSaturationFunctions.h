#ifndef EXPERIMENTAL_SATURATION_FUNCTIONS_H
#define EXPERIMENTAL_SATURATION_FUNCTIONS_H

#include "CSMP_definitions.h"
#include "CubicSpline.h"
#include "Element.h"

namespace csmp {
  
/**
   Saturation functions parameterised with experimental data, including
   relperm- and pc-s relations read from user-specified data (ASCII-text) 
   files are interpolated with cubic splines.
   
   Each curve is indentified by a unique integer code that represents the rocktype
   of interest.
   
   Input file format:
   
   @code
      To use experimental data for rock types you need to define
      RELPERM_MODEL in Numerical.txt file as EXPERIMENTAL
      Input file for rock types should be put in the execution folder
      and the file name that contains relative permeability data and
      capillary values should be 'rock type.txt'
      File structure is as follow:
     2 ->number of tables
     21 ->number of entries for table 0 (first table) excluding the first row
     -1.99513	-0.00488	0.000125	3.709875	-10557.3	-519.567 ->derivatives at start and end for kro, krw and pc
     0	1	0	5000 -> data (sw  kro krw pc)
     0.05	0.90024375	0.00000625	4472.135955
     0.1	0.8019	0.0001	3162.27766
     0.15	0.70624375	0.00050625	2581.988897
     0.2	0.6144	0.0016	2236.067977
     0.25	0.52734375	0.00390625	2000
     0.3	0.4459	0.0081	1825.741858
     0.35	0.37074375	0.01500625	1690.308509
     0.4	0.3024	0.0256	1581.13883
     0.45	0.24124375	0.04100625	1490.711985
     0.5	0.1875	0.0625	1414.213562
     0.55	0.14124375	0.09150625	1348.399725
     0.6	0.1024	0.1296	1290.994449
     0.65	0.07074375	0.17850625	1240.347346
     0.7	0.0459	0.2401	1195.228609
     0.75	0.02734375	0.31640625	1154.700538
     0.8	0.0144	0.4096	1118.033989
     0.85	0.00624375	0.52200625	1084.652289
     0.9	0.0019	0.6561	1054.092553
     0.95	0.00024375	0.81450625	1025.978352
     1	0	1	1000
     6 ->number of entries for table 1 (second table) excluding the first row
     -1	-1	1	1	0	0 ->derivatives
     0	1	0	0 ->data
     0.2	0.8	0.2	0
     0.4	0.6	0.4	0
     0.6	0.4	0.6	0
     0.8	0.2	0.8	0
     1	0	1	0
  @endcode
 
  @attention DO NOT SET ANY SATURATION ENDPOINTS (swr or snr) when using experimentally determined saturation functions
  as these will naturally emerge from the curves.
 
  @attention consequently, effective saturation also is already accounted for in terms of the spline values
  negative values are avoided by returning 0 or 1 at saturations below and above the end-point saturations.

   @author SKM 
   @date 14/1/2007
   @update 9/2/2019
 
*/
template<uint32_t dim, template<uint32_t> class USER>
class ExperimentalSaturationFunctions {
  public:
    /// reads capillary preessure and saturation functions from "rocktype" file
    explicit ExperimentalSaturationFunctions( const char* filename );
  
    /// TODO: revise to extract this from the curves in a consistent fashion: (sw-swr) / (1-swr-snr); use only in 2-phase flow simulations
    double EffectiveSaturation( Element<dim>* const ) const;
    
    /// TODO: revise to extract this from the curves in a consistent fashion: (sw-swr) / (1-swr-snr); use only in 2-phase flow simulations
    double EffectiveSaturation_at( Element<dim>* const, double sw ) const;

    // pc, kri, and derivative methods that use a user supplied saturation value
  
    double pc_at( Element<dim>* const, double sw ) const;
 
    double dpcds_at( Element<dim>* const, double sw ) const;

    double krw_at( Element<dim>* const, double sw ) const;
  
    double krn_at( Element<dim>* const, double sw ) const;

    double dkrwds_at( Element<dim>* const, double sw ) const;
 
    double dkrnds_at( Element<dim>* const, double sw ) const;

    double dpcdsw_at_Numerical( Element<dim>* const p, double sw, double h ) const;

    double dkrwds_at_Numerical( Element<dim>* const p, double sw, double h ) const ;

    double dkrnds_at_Numerical( Element<dim>* const p, double sw, double h ) const;
  
  
    // pc, kri, and derivative methods that use saturation values at the element barycentre (for FE mobility calculations etc.)
  
    double pc( Element<dim>* const ) const;
  
    double dpcds( Element<dim>* const ) const;
  
    double krw( Element<dim>* const ) const;
  
    double krn( Element<dim>* const ) const;
  
    double dkrwds( Element<dim>* const ) const;
  
    double dkrnds( Element<dim>* const ) const;

    double dpcdsw_Numerical( Element<dim>* const p, double h ) const;

    double dkrwds_Numerical( Element<dim>* const p, double h ) const;

    double dkrnds_Numerical( Element<dim>* const p, double h ) const;
  
    /// @return number of RRTs (reservoir rock types) for which saturation function values are stored
    size_t RockTypes() const;
  
    void Out() const;
    
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    /// reads rock types from file, establishing how many there are and returning this value
    size_t InitialiseReservoirRockTypes( const char* rt_file_name );
  
    /// converts the floating-point rocktype identifier into a positive integer 0..254
    size_t RockType( Element<dim>* const ) const;
  
    std::vector<csmp::CubicSpline> kr1_, kr2_, pc_;
    const double max_derivative_         = 5.0e+8; ///< the absolute value of any derivative calculated herein must be less than this value
    const double max_capillary_pressure_ = 2.0e+7; ///< 20 MPa ~ tensile strength of the rock
};
  
} // end namespace csmp

#endif /* CSMP_EXPERIMENTAL_SATURATION_FUNCTIONS_H */
