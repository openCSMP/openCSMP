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

   @author SKM 
   @date 14/1/2007
 
*/
template<size_t dim, template<size_t> class USER>
class ExperimentalSaturationFunctions {
  public:
    /// TODO: reads capillary preessure and saturation functions from "rocktype" file
    explicit ExperimentalSaturationFunctions( const char* filename );
    
    /// Effective Saturation function
    double64 EffectiveSaturation( const Element<dim>* const ) const ;
  
    /// Effective Saturation function fro saturation S
    double64 EffectiveSaturation_at( const Element<dim>* const , double64 ) const ;
  
    /// The Capillary pressure Eq. (2) from Skjaeveland et al. 2000
    double64 pc( const Element<dim>* const );
  
    /// The first derivative of Capillary pressure.
    double64 dpcds( const Element<dim>* const )  const ;
  
    /// The water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 krw( const Element<dim>* const ) const ;
    
    double64 krw_at( const Element<dim>* const, double64 sw ) const ;
  
    /// The CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 krn( const Element<dim>* const ) const ;
    
    double64 krn_at( const Element<dim>* const, double64 sw ) const ;
  
    /// The first relative of water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 dkrwds( const Element<dim>* const ) const ;
    
    double64 dkrwds_at( const Element<dim>* const, double64 sw ) const ;
  
    /// The first relative of CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 dkrnds( const Element<dim>* const ) const ;
    
    double64 dkrnds_at( const Element<dim>* const, double64 sw ) const ;

    /// Numerical derivatives of first derivatives of relative permeability of water and CO2
    double64 dkrwds_Numerical( const Element<dim>* const p, double64 h ) const ;

    double64 dkrwds_at_Numerical( const Element<dim>* const p, double64 sw, double64 h ) const ;

    double64 dkrnds_Numerical( const Element<dim>* const p, double64 h ) const ;
  
    double64 dkrnds_at_Numerical( const Element<dim>* const p, double64 sw, double64 h ) const;
  
    /// returns number of RRTs (reservoir rock types) for which saturation function values are stored
    size_t RockTypes() const;
  
    void Out() const;
  
    
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    /// reads rock types from file, establishing how many there are and returning this value
    size_t InitialiseReservoirRockTypes( const char* rt_file_name );
  
    /// converts the floating-point rocktype identifier into a positive integer 0..254
    size_t RockType( const Element<dim>* const ) const;
  
    std::vector<csmp::CubicSpline> kr1_, kr2_, pc_;
    const double64 max_derivative_ = 1.0e+8; ///< the absolute value of any derivative should be less than this
};
  
} // end namespace csmp

#endif /* CSMP_EXPERIMENTAL_SATURATION_FUNCTIONS_H */
