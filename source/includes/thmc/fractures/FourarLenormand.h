#ifndef FOURAR_LENORMAND_H
#define FOURAR_LENORMAND_H

#include "TwoPhaseModel.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;

template<uint32_t dim>
class FourarLenormand : public TwoPhaseModel<dim> {
  public:
    FourarLenormand( const PropertyDatabase<dim>& database, const char* fractureAperture );
    FourarLenormand( const PropertyDatabase<dim>& database,
                     const char* fractureAperture,
                     const char* permeability,
                     const char* my_non_wetting,
                     const char* my_wetting,
                     const char* rho_non_wetting,
                     const char* rho_wetting,
                     const char* sat_wetting,
                     const char* res_sat_non_wetting,
                     const char* res_sat_wetting );

    void Initialize( const Element<dim>& e ) override final;

    // relative permeabilities
    double krn_Phase() const override final;
    double krw_Phase() const override final;

    // derivatives of relative permeabilities
    double dkrnds_Phase() const override final;
    double dkrwds_Phase() const override final;

    // capillary pressure
    double pc_Phase( ) const override final;

    // capillary pressure derivatives
    double dpcds_Phase( ) const override final;

    // inverse capillary pressure function
    double Sw_Phase( double pc_Phase ) const override final;

    // inverse capillary pressure derivative
    double dsdpc_Phase( double pc_Phase ) const override final;

  private:

    FourarLenormand();

    double  fractureAperture_;
    Index   fractureApertureKey_;
};




/**
@class FourarLenormand FourarLenormand "two_phase_flow/FourarLenormand.h"
@author P. Lang
@date 2010

@warning Needs debugging, simulation not stable!!!

Relative Permeability model for NFR after Fourar and Lenormand, 1998

*/

} // end namespace csmp

#endif // FOURAR_LENORMAND_H
