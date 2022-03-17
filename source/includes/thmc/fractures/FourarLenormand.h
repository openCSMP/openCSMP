#ifndef FOURAR_LENORMAND_H
#define FOURAR_LENORMAND_H


#include <cmath>

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

    virtual ~FourarLenormand();

    virtual void Initialize( const Element<dim>& e );

    // relative permeabilities
    virtual double krn_Phase() const;
    virtual double krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double dkrnds_Phase() const;
    virtual double dkrwds_Phase() const;

    // capillary pressure
    virtual double pc_Phase( ) const;

    // capillary pressure derivatives
    virtual double dpcds_Phase( ) const;

    // inverse capillary pressure function
    virtual double Sw_Phase( double pc_Phase ) const;

    // inverse capillary pressure derivative
    virtual double dsdpc_Phase( double pc_Phase ) const;

  private:

    FourarLenormand();

    double  fractureAperture_;
    Index     fractureApertureKey_;

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
