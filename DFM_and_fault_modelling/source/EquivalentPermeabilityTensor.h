#ifndef EQUIVALENTPERMEABILITYTENSOR_H
#define EQUIVALENTPERMEABILITYTENSOR_H

#include "CSMP_definitions.h"
#include "Matrix.h"

// Model
#include "RegionInterface.h"

// PDE solution framework
#include "SAMG_Settings.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"

namespace csmp {

template<uint32_t> class Model;

template<uint32_t dim>
class EquivalentPermeabilityTensor
{
public:
    EquivalentPermeabilityTensor( Model<dim>& );
    ~EquivalentPermeabilityTensor();
  
    void InitializeViscosity( double vis_val );

    void ConstantViscosity( bool vis_flag );

    /// initialize random sampling parameters
    void InitializeRandomSampling( size_t no_bins = 50, size_t max_no_sample = 20, size_t min_no_elements = 100, double sample_size_reducing_factor = -1. );

    void ComputeModelVelocityAndPressureGradient( Model<dim>& );

    /// computes the full equivalent permeability tensor of the model.
    void ComputeEquivalentPermeabilityTensor( Model<dim>& model );

    /// computes the full equivalent permeability tensor of the region of interest.
    void ComputeEquivalentPermeabilityTensor( Model<dim>& model, const std::string regionName );

    /// computes the full equivalent permeability tensor for the random region sampling.
    void ComputeEquivalentPermeabilityTensorByRandomSampling( Model<dim>& model );

    /// return components of permeability tensor
    double operator()( size_t, size_t ) const;
    
    /// in the calculation of off-diagonal terms: maximum difference between diagonal and off-diagonal terms
    void PermeabilityTolerance( double );

    /// sets output files name
    void SpecifyOutputFileNames( std::string tensorFileName = "tensor.txt", std::string statisticsFileName = "statistics.txt" );

    /// outputs computed tensor statistics by random sampling to file
    void OutputStatistics() const;
  
    /// output tensor
    void Out( TensorVariable<dim>& ) const;

    /// outputs computed tensor to the file
    void Out() const;

 private:
    void Left2RightFlow( Model<dim>& );
    void Bottom2TopFlow( Model<dim>& );
    void Back2FrontFlow( Model<dim>& );
    void PressureField ( Model<dim>&, const std::string, PDE_Integrator<dim,Element>& );
    void VelocityAndPressureGradientField( Model<dim>&, const Index&, VectorVariable<dim>&, const Index&, VectorVariable<dim>&, const Index& );
    void SolveOverdeterminedSystem( const Matrix&, const std::vector<double>&, std::vector<double>& );
    void VelocityAndPressureGradientVolumeAverage( const Model<dim>&, const std::string regionName);
    /// if viscosity varies across the domain of interest
    void ViscosityVolumeAverage( const Model<dim>&, const std::string regionName );
    void ComputeTensor();
    void TensorPrincipals();

    void InitializeSolver();

    void InitializeIndexKeys( Model<dim>& );

	void InitializeTransmissivity(Model<dim>&, const std::string regionName);
  
 private:
    VectorVariable<dim> velocity_left2right_;
    VectorVariable<dim> velocity_bottom2top_;
    VectorVariable<dim> velocity_back2front_;
    VectorVariable<dim> gradp_left2right_;
    VectorVariable<dim> gradp_bottom2top_;
    VectorVariable<dim> gradp_back2front_;
    TensorVariable<dim> keq_; // the computed equivalent full permeability tensor
private:
    Index pressure_left2right_key_;
    Index pressure_bottom2top_key_;
    Index pressure_back2front_key_;
    Index velocity_left2right_key_;
    Index velocity_back2front_key_;
    Index velocity_bottom2top_key_;
    Index gradp_left2right_key_;
    Index gradp_bottom2top_key_;
    Index gradp_back2front_key_;
    Index permeability_key_;
    Index viscosity_key_;
    Index trans_key_;
    Index thic_key_;
    Index unitnormal_key_;

    double theta_;
    double kmin_, kmax_;
    double plunge_min_, trend_min_;
    double plunge_max_, trend_max_;
    double vx1_avg_, vx2_avg_, vx3_avg_, vy1_avg_, vy2_avg_, vy3_avg_, vz1_avg_, vz2_avg_, vz3_avg_;
    double dpx1_avg_, dpx2_avg_, dpx3_avg_, dpy1_avg_, dpy2_avg_, dpy3_avg_, dpz1_avg_, dpz2_avg_, dpz3_avg_;
    double viscosity_avg_;
    double eps_;
    double viscosity_;
    double matrix_perm_;

    // sampling settings    
    size_t no_bins_;
    size_t max_no_sample_;
    size_t min_no_elements_;
    size_t sample_elements_;
    double sample_size_;
    double sampleSizeReducingFactor_;
    std::vector<double> ld_;
    double modelCharLength_;
	  std::vector<double> model_sample_size_;
	  std::vector<std::vector<double>> model_bins_;
	  std::vector<std::vector<size_t>> model_freqs_;
    std::string tensor_file_name_;
	  std::string statistics_file_name_;
  
 // SKM FIX: use encapsulation! - this is C++
 private:
    std::set<int32_t>  fracture_elmt_indices_, matrix_elmt_indices_;
    SAMG_Settings      settings_;
    bool               const_viscosity_;
};

} // end csmp

#endif // EQUIVALENTPERMEABILITYTENSOR_H
