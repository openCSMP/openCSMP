#ifndef CSMP_VELOCITY_AND_VOLUME_FLUX_VISITOR_H
#define CSMP_VELOCITY_AND_VOLUME_FLUX_VISITOR_H

#include "Visitor.h"
#include "FiniteElementManager.h"

namespace csmp {

template<uint32_t> class VectorVariable;
template<uint32_t> class Model;

/// IMPORTANT: for legacy purposes, the transport class assumes that the supplied velocity has already been scaled
/// to adjust for facet areas coming from the thickness attribute (LDE's).  This is perhaps a little unintuitive,
/// since any new user needs to know beforehand that this has to happen.  The user should always provide the
/// "real" velocity (not scaled) and the transport class should deal with thickness through internal scaling of its
/// facet areas.
/// This class calculates single phase velocity in its "de-scaled" form.  No thickness factor is used.  Thus, using this
/// calculation with a problem containing LDE's and using legacy transport classes (e.g. NCFVT) requires
/// an extra step of multiplication by a "thickness" attribute.
///
/// The class supports and will take advantage of Open MP by dividing the evaluation of the contributions equally among available
/// threads.
///
/// -- Julian M. 26-09-2015

template<uint32_t dim>
class SinglePhaseVelocityVisitor : public Visitor<dim> {
  public:

    SinglePhaseVelocityVisitor(Model<dim>&,
                               const char* porosity,
                               const char* conductivity,
                               const char* fluid_density,
                               const char* fluid_pressure,
                               const char* velocity,               // results: Darcy Velocity
                               const char* model_gravity_vector = NULL,
                               const char* pore_velocity = NULL,   // results: True velocity (v/phi)
                               const char* volume_flux = NULL,     // results:
                               const char* nodal_velocity = NULL,
                               const char* nodal_pore_velocity = NULL,
                               const char* nodal_volume_flux = NULL);

    ~SinglePhaseVelocityVisitor();

    virtual void Visit(Element<dim>* e);
    void ComputeContribution(Element<dim>* e);
    virtual void Visit(Region<dim>* region);
    virtual void Visit(Model<dim>* model);
    
  private:
    const size_t components_;

    csmp::Index phi_key_,  rho_key_,   conductivity_key_, fpres_key_,
                velo_key_, ivelo_key_, nvelo_key_, nivelo_key_,
                flux_key_, nflux_key_;

    VectorVariable<dim>           gravity_unit_vector_;
    bool                          with_gravity_;
    double                      gravitational_acceleration_;
    uint32_t                        application_cycle_;
    std::pair<double,double>  minmaxV_, minmaxF_;

#if defined(_OPENMP )
    std::vector<FiniteElementManager> femgrs_; // one manager per thread
#endif
};

} // csmp

#endif // CSMP_VELOCITY_AND_VOLUME_FLUX_VISITOR_H
