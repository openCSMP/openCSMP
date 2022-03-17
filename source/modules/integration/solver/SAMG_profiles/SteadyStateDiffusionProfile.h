#ifndef STEADYSTATEDIFFUSIONPROFILE_H
#define STEADYSTATEDIFFUSIONPROFILE_H

#include "SAMG_Profile.h"
#include "SteadyStateDiffusor.h"
#include "Region.h"

namespace csmp{

template<uint32_t> class Model;
template<uint32_t dim>
class SteadyStateDiffusionProfile : public SAMG_Profile {
  public:
    SteadyStateDiffusionProfile( Model<dim>& model,
                                 const char* diffusivity,
                                 const char* diffusingVariable,
                                 const char* spatialSourceVariable );

    SteadyStateDiffusionProfile( Model<dim>& model,
                                 const char* diffusivity,
                                 const char* diffusingVariable,
                                 const char* spatialSourceVariable,
                                 const char* gradientVariable,
                                 double gradient_multiplier );


    /// Solve for steadystate pressure diffusion
    virtual bool Solve( double modelTime );

    /// bridging to solver
    SteadyStateDiffusor<dim,Region>& Solver() { return ssds_; }

  private:
    SteadyStateDiffusionProfile();

    bool AdjustSolverSettings();

  private:
    Model<dim>& model_;

    SteadyStateDiffusor<dim,Region> ssds_;

    bool firstCall_;
    std::string   dumpFileName_;

};

} // csmp


#endif // STEADYSTATEDIFFUSIONPROFILE_H
