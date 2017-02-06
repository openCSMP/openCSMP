#ifndef FRACTURE_MATRIX_UPSCALED_H
#define FRACTURE_MATRIX_UPSCALED_H

#include "TwoPhaseModel.h"
#include "GenericTransferFunction.h"

namespace csmp {

/**

ensemble relative permeability model for fracture matrix systems
parameterized with, qf/qm ratio, Af(sw), and GTF
the BrooksCorey model is retained to calculate lambda overbar and pc
for the transfer function*/
template<size_t dim>
class FractureMatrixUpscaled : public TwoPhaseModel<dim> {
  public:
    FractureMatrixUpscaled( const PropertyDatabase<dim>& database,
                            const char* permeability,
                            const char* visc_nw, const char* visc_w,
                            const char* rho_nw, const char* rho_w,
                            // Brooks-Corey params
                            const char* lamda, const char* pc_entry,
                            // fracture matrix terms
                            const char* specific_fracture_matrix_interface_area, // new
                            const char* phi_fracture, // new element variable (per unit volume of rock)
                            const char* phi_matrix,
                            const char* qfqm_ratio,  // new element variable
                            const char* volume_flux,
                            // for the GTF
                            const char* sw_initial,  // new element variable
                            const char* block_radius, // new element variable
                            const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    FractureMatrixUpscaled( const PropertyDatabase<dim>& database,
                            const char* permeability,
                            const char* visc_nw, const char* visc_w,
                            const char* rho_nw, const char* rho_w,
                            // Brooks-Corey params
                            const char* lamda, const char* pc_entry,
                            // fracture matrix terms
                            const char* specific_fracture_matrix_interface_area, // new
                            const char* phi_fracture, // new element variable (per unit volume of rock)
                            const char* phi_matrix,
                            const char* qfqm_ratio,  // new element variable
                            const char* volume_flux,
                            // for the GTF
                            const char* sw_initial,  // new element variable
                            const char* block_radius, // new element variable
                            const char* rsnw, const char* rsw,
                            const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    FractureMatrixUpscaled( double permeability,
                            double visc_nw, double visc_w,
                            double rho_nw, double rho_w,
                            double lambda, double pc_entry,
                            double specific_fracture_matrix_interface_area, 
                            double phi_fracture, 
                            double phi_matrix,
                            double qfqm_ratio,  
                            double volume_flux,
                            double sw_initial,  
                            double block_radius,
                            double swr_, double snr_ ); 

    virtual ~FractureMatrixUpscaled();

    virtual void Initialize( const Element<dim>& e );

    virtual void InitializeForNode( const Element<dim>& e,
                                    size_t node );

    void InitializeForSw( double sw );

    // relative permeabilities
    virtual double64 krn_Phase() const;
    virtual double64 krw_Phase() const;


    // derivative of fractional flow (advection multipliers)
    virtual double64 dfds() const;

    // maximum absolute value returned by dfdS
    virtual double64 MaxFractionalFlowDerivative() const;

    // derivatives of gravitational flow (advection multipliers)
    virtual double64 dGds( ) const;

    virtual double64 pc_Phase( ) const;

    // capillary pressure derivatives (treat seff as for previous function)
    virtual double64 dpcds_Phase( ) const { std::cout <<"\ndpcds_Phase: not needed\n"; return 0.; };

    virtual void Out( std::ostream& os, size_t phase ) const;


  private:

    FractureMatrixUpscaled();

    double64  Af_sw() const;

    csmp::Index  pd_key_, lambda_key_, Af_key_, phim_key_, phif_key_,
                 qfqm_key_, flux_key_, r_key_, swi_key_;
    double64     Af_, pd_, lambda_, phim_, phif_, qfqm_,
                 qv_, radius_;

    ScalarVariable  swi_;
};


/// exponential function from 0.1 to 1 at sw=1, should be fitted to simulation results
template<size_t dim>
inline double64 FractureMatrixUpscaled<dim>::Af_sw() const
 {
    // .1 * exp(2.3*x^3)
    return Af_ * 0.1 * std::exp( 2.3 * std::pow( TwoPhaseModel<dim>::sat_, 3. ) );
 }

} // end namespace csmp

#endif
