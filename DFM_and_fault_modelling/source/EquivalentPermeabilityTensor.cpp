#include "EquivalentPermeabilityTensor.h"

#include "CSMP_highLevelUtilities.h"
#include "PropertyHandle.h"
#include "Matrix.h"

// Model
#include "Model.h"
#include "Region.h"

// Solver
#include "LinearSolver.h"
#ifdef CSMP_WITH_SAMG
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#endif
#include "LUdcmp_Solver.h"

// PDE solution framework
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"

// Interface
#include "VTK_Interface.h"

#include "StatisticalDistributionGenerator.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
EquivalentPermeabilityTensor<dim>::EquivalentPermeabilityTensor( Model<dim>& model )
{
    model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );
    model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );

    InitializeSolver();

    matrix_perm_ = printRangeOfVariable( model, "permeability", false );
    PermeabilityTolerance( 1e-02 * matrix_perm_ );

    InitializeIndexKeys( model );

    InitializeViscosity( 1.0e-3 );

    ConstantViscosity( true );

    InitializeTransmissivity( model, "Model" );
	
    ComputeModelVelocityAndPressureGradient( model );
}






template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::InitializeSolver()
{
#ifdef CSMP_WITH_SAMG
    // targeting SAMG DLL 1 for this pressure solver
    // iout
    settings_.Set_iout1( 0 );
    settings_.Set_iout2( 0 );
    settings_.Set_nxtyp(1);
    // if ( !verbose_ ) settings_.Set_idmp( -1 );
    // one-time solver set-up: nothing is remembered for next try
    settings_.Set_iswit(5);
    // SAMG solution criteria
    settings_.Set_eps(0.); // absolute criterion
    // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
    settings_.Set_a_cmplx(2);
    // Pre-adjust SAMG mesh complexity, based on solver output
    settings_.Set_g_cmplx(1.5);
    settings_.Set_w_avrge(2);

    //settings_.Set_eps( 1.E-20 );
  #ifdef SAMG_OUTPUT_TO_FILE
    // trigger output to file
    fluid_pressure.GetSolverSettings().Set_idmp( 8 );
    // // define SAMG file output format for reduced file size, idmp > 1 is required
    fluid_pressure.GetSolverSettings().Set_ioform( "f" );
    // set filename for SAMG file output other than default "level", idmp > 1 is required
    fluid_pressure.GetSolverSettings().Set_filnam_dump( "ReservoirSimulator_steady_state_p" );
  #endif
  #endif
}



template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::InitializeIndexKeys( Model<dim>& model )
{
    pressure_left2right_key_ = model.Database().StorageKey("fluid pressure1");
    velocity_left2right_key_ = model.Database().StorageKey("velocity1");
    gradp_left2right_key_    = model.Database().StorageKey("fluid pressure gradient1");
    pressure_bottom2top_key_ = model.Database().StorageKey("fluid pressure2");
    velocity_bottom2top_key_ = model.Database().StorageKey("velocity2");
    gradp_bottom2top_key_    = model.Database().StorageKey("fluid pressure gradient2");
    permeability_key_        = model.Database().StorageKey("permeability");
    viscosity_key_           = model.Database().StorageKey("fluid viscosity");
    trans_key_               = model.Database().StorageKey("transmissivity");
    thic_key_                = model.Database().StorageKey("thickness");
    unitnormal_key_          = model.Database().StorageKey("unit normal");
    if ( dim == 3) {
        pressure_back2front_key_  = model.Database().StorageKey("fluid pressure3");
        velocity_back2front_key_  = model.Database().StorageKey("velocity3");
        gradp_back2front_key_     = model.Database().StorageKey("fluid pressure gradient3");
    }
}




/// returns components of permeability tensor
template<uint32_t dim>
double EquivalentPermeabilityTensor<dim>::operator()( size_t i, size_t j ) const
{
   return keq_( static_cast<uint32_t>(i), static_cast<uint32_t>(j) );
}
  
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::InitializeViscosity(double vis_val)
{
    viscosity_ = vis_val;
}

template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ConstantViscosity( bool vis_flag )
{
    const_viscosity_ = vis_flag;
}




template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::InitializeRandomSampling( size_t no_bins, size_t max_no_sample, size_t min_no_elements, double sample_size_reducing_factor )
{
    no_bins_                  = no_bins;
    max_no_sample_            = max_no_sample;
    min_no_elements_          = min_no_elements;
    sampleSizeReducingFactor_ = sample_size_reducing_factor;
    // dimensionless sample size
    if ( sampleSizeReducingFactor_ < 0 ) {
		double init_list[] = { 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.08, 0.06, 0.04, 0.02 };
		ld_.assign(init_list, init_list + 13);
    }
    else {
        double modelToSampleRatio = 1;
        while ( 1. / modelToSampleRatio > 0.01 ) {
            modelToSampleRatio += sampleSizeReducingFactor_;
            ld_.push_back( 1. / modelToSampleRatio );
        }
    }
}





template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::InitializeTransmissivity( Model<dim>& model, const string regionName )
{
    double transmissivity;
    Region<dim>& r_ref(model.Region(regionName.c_str()) );
    if ( const_viscosity_ ) {
         for ( auto it_e = r_ref.CellsBegin(); it_e != r_ref.CellsEnd(); ++it_e ) {
             transmissivity = ((*it_e)->Read(permeability_key_) / viscosity_) * (*it_e)->Read(thic_key_);
             (*it_e)->Store(trans_key_, makeScalar(PLAIN, transmissivity));
           }
    }
    else {
         for ( auto it_e = r_ref.CellsBegin(); it_e != r_ref.CellsEnd(); ++it_e ) {
             transmissivity = ((*it_e)->Read(permeability_key_) / (*it_e)->Read(viscosity_key_)) * (*it_e)->Read(thic_key_);
             (*it_e)->Store(trans_key_, makeScalar(PLAIN, transmissivity));
         }
    }
}

/// Computes the full equivalent permeability tensor for the model.
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ComputeEquivalentPermeabilityTensor( Model<dim>& model )
{   
    Region<dim>& r_ref(model.Region("Model") );
    sample_elements_ = r_ref.Cells();
    if ( dim == 2 ) sample_size_ = pow( r_ref.Volume(), 1. / 2. );
    if ( dim == 3 ) sample_size_ = pow( r_ref.Volume(), 1. / 3. );
    ComputeEquivalentPermeabilityTensor( model, "Model");
    Out();
}





/// Computes the full equivalent permeability tensor in the region of interest.
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ComputeEquivalentPermeabilityTensor( Model<dim>& model, const string regionName )
{
    VelocityAndPressureGradientVolumeAverage( model, regionName );
    if ( const_viscosity_ ) viscosity_avg_ = viscosity_;
    else ViscosityVolumeAverage( model, regionName );
    ComputeTensor();
    TensorPrincipals();
}



/// Computes the full equivalent permeability tensor in the region of interest.
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ComputeEquivalentPermeabilityTensorByRandomSampling( Model<dim>& model )
{
    Point<dim> xyz_min, xyz_max;
    model.MinMaxCoordinates(xyz_min, xyz_max);
    double modelxLength = xyz_max[0] - xyz_min[0];
    double modelyLength = xyz_max[1] - xyz_min[1];
    double modelzLength;
    if ( dim == 3) modelzLength = xyz_max[2] - xyz_min[2];
    modelCharLength_ = modelxLength;
    if (modelyLength < modelCharLength_) modelCharLength_  = modelyLength;
    if ( dim == 3) {
        if (modelzLength < modelCharLength_) modelCharLength_  = modelzLength;
    }
    bool minSmapleElements = false;
    vector<double> kmax_to_km;
    Point<dim> sampleCenter, sampleMinPoint, sampleMaxPoint;
    size_t sampleCounter = 0;
    while ( !minSmapleElements && sampleCounter < ld_.size() ) {
        cout << "----------------------" << endl;
        cout << " sample " << sampleCounter + 1 << " of " << ld_.size() << endl;
        cout << "----------------------";
        sample_size_ = modelCharLength_ * ld_[sampleCounter];
        size_t sameSizeSampleCounter = 0;
        while ( sameSizeSampleCounter < max_no_sample_ ) {
            bool pickedSample = false;
            while ( !pickedSample ) {
                // randomly taking a point between the model's lower left and upper right corners
                sampleCenter = xyz_min + double(rand()) / double(RAND_MAX) * (xyz_max - xyz_min);
                sampleMinPoint = sampleCenter - (sample_size_ / 2.);
                sampleMaxPoint = sampleCenter + (sample_size_ / 2.);
                if ( sampleMinPoint.IsBetween(xyz_min, xyz_max) && sampleMaxPoint.IsBetween(xyz_min, xyz_max) ) pickedSample = true;
            }
            const string sampleRegionName = "SamplingRegion";
            // creats a new region within the model
            model.FormRectangularRegion( sampleRegionName.c_str(), sampleMinPoint, sampleMaxPoint );
            // calculates equivalent permeability tensor for the new region
            Region<dim>& r_ref( model.Region(sampleRegionName.c_str()) );
            sample_elements_ = r_ref.Cells();
            if ( sample_elements_ < min_no_elements_ ) {
                minSmapleElements = true;
                if ( sameSizeSampleCounter != 0 ) model_sample_size_.push_back(sample_size_);
                break;
            }
            else {
                ComputeEquivalentPermeabilityTensor( model, sampleRegionName );
                if ( kmin_ > 0. && kmax_ > 0. ) {
                    Out();
                    sameSizeSampleCounter += 1;
                    if ( max_no_sample_ != 1 ) {
                        kmax_to_km.push_back( kmax_ / matrix_perm_ );
                    }
                }
            }
            // removing current sampling region
            model.RemoveRegion( sampleRegionName.c_str(), true );
        }
        if ( !minSmapleElements ) model_sample_size_.push_back(sample_size_);
        // creats statistics of the sampling method
        if ( max_no_sample_ != 1 ) {
            StatisticalDistributionGenerator same_size_stat;
            same_size_stat.Histogram(kmax_to_km, no_bins_);
            model_bins_.push_back(same_size_stat.Bins());
            model_freqs_.push_back(same_size_stat.Frequencies());
            kmax_to_km.clear();
        }
        sampleCounter += 1;
    }
    if ( max_no_sample_ != 1 ) OutputStatistics();
}




template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ComputeModelVelocityAndPressureGradient( Model<dim>& model )
{
// flow calculations
    Left2RightFlow( model );
    Bottom2TopFlow( model );
    if ( dim == 3 ) Back2FrontFlow( model );
}




// left to right velocity and pressure gradient calculation
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::Left2RightFlow( Model<dim>& model )
{
// Dirichlet boundary condition
	model.InputBoundaryValue( LEFT, "fluid pressure1", makeScalar(DIRICH, 1.0e8) );
    model.InputBoundaryValue( RIGHT, "fluid pressure1", makeScalar(DIRICH, 1.0e5) );
/*
// linear boundary condition
	vector<ScalarVariable> bvalues = { makeScalar(DIRICH, 1.0e8), makeScalar(DIRICH, 1.0e5) };
	if (dim == 3) {
		bvalues.push_back(makeScalar(DIRICH, 1.0e8));
		bvalues.push_back(makeScalar(DIRICH, 1.0e5));
	}
	model.LinearBoundaryCondition( BOTTOM, "fluid pressure1", bvalues );
	model.LinearBoundaryCondition( TOP, "fluid pressure1", bvalues );
*/
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                  samg_solver( &settings );
    PDE_Integrator<dim,Element>  pressure_left2right(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER   linear_solver;
    PDE_Integrator<dim,Element>  pressure_left2right(linear_solver);
#endif
  PressureField( model, "fluid pressure1", pressure_left2right );
  VelocityAndPressureGradientField( model, pressure_left2right_key_, velocity_left2right_, 
                                    velocity_left2right_key_, gradp_left2right_, gradp_left2right_key_ );

// output
//    VTK_Interface<dim> vtk_output;
//    vtk_output.OutputDataToVTK(model, "pressure_left2right", "fluid pressure1", 0);
//    vtk_output.OutputDataToVTK(model, "velocity_left2right", "velocity1", 0);
//    vtk_output.OutputDataToVTK(model, "gradp_left2right", "fluid pressure gradient1", 0);

}




// bottom to top velocity and pressure gradient calculation
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::Bottom2TopFlow( Model<dim>& model )
{
// Dirichlet boundary condition
    model.InputBoundaryValue( BOTTOM, "fluid pressure2", makeScalar(DIRICH, 1.0e8) );
    model.InputBoundaryValue( TOP, "fluid pressure2", makeScalar(DIRICH, 1.0e5) );
/*
// linear boundary condition
	vector<ScalarVariable> bvalues = { makeScalar(DIRICH, 1.0e8), makeScalar(DIRICH, 1.0e5) };
	if (dim == 3) {
		bvalues.push_back(makeScalar(DIRICH, 1.0e8));
		bvalues.push_back(makeScalar(DIRICH, 1.0e5));
	}
	model.LinearBoundaryCondition( LEFT, "fluid pressure2", bvalues );
	model.LinearBoundaryCondition( RIGHT, "fluid pressure2", bvalues );
*/
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                  samg_solver( &settings );
    PDE_Integrator<dim,Element>  pressure_bottom2top(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER   linear_solver;
    PDE_Integrator<dim,Element>  pressure_bottom2top(linear_solver);
#endif
    // releasing any pre-assigned boundary conditions
  PressureField( model, "fluid pressure2", pressure_bottom2top );
  VelocityAndPressureGradientField( model, pressure_bottom2top_key_, velocity_bottom2top_, 
                                    velocity_bottom2top_key_, gradp_bottom2top_, gradp_bottom2top_key_ );

// output
//    VTK_Interface<dim> vtk_output;
//    vtk_output.OutputDataToVTK(model, "pressure_bottom2top", "fluid pressure2", 0);
//    vtk_output.OutputDataToVTK(model, "velocity_bottom2top", "velocity2", 0);
//    vtk_output.OutputDataToVTK(model, "gradp_bottom2top", "fluid pressure gradient2", 0);

}




// back to front velocity and pressure gradient calculation
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::Back2FrontFlow( Model<dim>& model )
{
    model.InputBoundaryValue( BACK, "fluid pressure3", makeScalar(DIRICH, 1.0e8) );
    model.InputBoundaryValue( FRONT, "fluid pressure3", makeScalar(DIRICH, 1.0e5) );
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                  samg_solver( &settings );
    PDE_Integrator<dim,Element>  pressure_back2front(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER   linear_solver;
    PDE_Integrator<dim,Element>  pressure_back2front(linear_solver);
#endif
    PressureField( model, "fluid pressure3", pressure_back2front );
    VelocityAndPressureGradientField( model, pressure_back2front_key_, velocity_back2front_, 
                                      velocity_back2front_key_, gradp_back2front_, gradp_back2front_key_ );

//output
//    VTK_Interface<dim> vtk_output;
//    vtk_output.OutputDataToVTK(model, "pressure_back2front", "fluid pressure3", 0);
//    vtk_output.OutputDataToVTK(model, "velocity_back2front", "velocity3", 0);
//    vtk_output.OutputDataToVTK(model, "gradp_back2front", "fluid pressure gradient3", 0);

}




template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::PressureField( Model<dim>& model,
                                                       const string& pressureVariableName,
                                                       PDE_Integrator<dim,Element>& fluid_pressure )
{
    // conductance matrix [K] on the left-hand side
    NumIntegral_dNT_lhsop_dN_dV<dim> conductance(model.Database(), "transmissivity", pressureVariableName.c_str(), pressureVariableName.c_str());
    // source vector {Q} on the right-hand side
    NumIntegral_NT_rhsop_N_dV<dim> source( model.Database(), "fluid volume source", pressureVariableName.c_str() );
    // add PDE_Operators to the FE Algorithm
    fluid_pressure.Add( &conductance );
    fluid_pressure.Add( &source );
    // solve [K]{p} = {Q} by passing the algorithm to the model
    model.Apply( fluid_pressure );
}




template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::VelocityAndPressureGradientField( Model<dim>& model, const Index& pressure_key_,
                                                     VectorVariable<dim>& velocity, const Index& velocity_key_,
                                                     VectorVariable<dim>& gradp, const Index& gradp_key_ )
{
    DenseMatrix<DM_MIN> dNT;
    // creating a reference to the model region
    Region<dim>& r_ref(model.Region("Model") );
    double element_mobility;
    for ( auto it_e = r_ref.CellsBegin(); it_e != r_ref.CellsEnd(); ++it_e ) {
        // total mobility at element barycenter
        element_mobility = (*it_e)->Read(trans_key_) / (*it_e)->Read(thic_key_);
        // interpolation function derivatives at element barycenter
        (*(*it_e)).dN_AtBaryCenter(dNT, 1U);
        // velocity at element barycenter
        velocity = 0.;
        gradp = 0.;
        double dpdx =0.;
        double dpdy =0.;
        double dpdz =0.;
        for (uint32_t i = 0; i<(*it_e)->Nodes(); ++i)
        {
            // pressure at element nodes
            double pf = (*it_e)->N(i)->Read(pressure_key_);

            // pressure gradient at element barycenter
            dpdx += pf * dNT(0,i);
            dpdy += pf * dNT(1,i);
            if ( dim == 3) dpdz += pf * dNT(2,i);
        }
        gradp(0) = dpdx;
        gradp(1) = dpdy;
        if ( dim == 3) gradp(2) = dpdz;
        velocity(0) = - element_mobility * dpdx;
        velocity(1) = - element_mobility * dpdy;
        if ( dim == 3) velocity(2) = - element_mobility * dpdz;
        (*it_e)->Store(velocity_key_, velocity);
        (*it_e)->Store(gradp_key_, gradp);
    }
}




/// volume averaging of the velocity and pressure gradient fields.
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::VelocityAndPressureGradientVolumeAverage( const Model<dim>& model, const string regionName )
{
    double tot_vol(0.);
    double tot_vx1_vol   = 0.; double tot_vx2_vol   = 0.; double tot_vx3_vol   = 0.;
    double tot_vy1_vol   = 0.; double tot_vy2_vol   = 0.; double tot_vy3_vol   = 0.;
    double tot_vz1_vol   = 0.; double tot_vz2_vol   = 0.; double tot_vz3_vol   = 0.;
    double tot_dpdx1_vol = 0.; double tot_dpdx2_vol = 0.; double tot_dpdx3_vol = 0.;
    double tot_dpdy1_vol = 0.; double tot_dpdy2_vol = 0.; double tot_dpdy3_vol = 0.;
    double tot_dpdz1_vol = 0.; double tot_dpdz2_vol = 0.; double tot_dpdz3_vol = 0.;
	double element_vol;

    const Region<dim>& r_ref(model.Region(regionName.c_str()));
    for ( typename vector<Element<dim>*>::const_iterator it_e = r_ref.CellsBegin(); it_e != r_ref.CellsEnd(); ++it_e )
      {
         element_vol = (*it_e)->Volume() * (*it_e)->Read(thic_key_);
        
        // if the element belongs to any of the fracture regions
        // if ( (*fracture_elmt_indices_.find(rint((*it_e)->Read(idx_key_))) != fracture_elmt_indices_.end() )

          tot_vol += element_vol;
          (*it_e)->Read(velocity_left2right_key_,velocity_left2right_);
          double element_vx1 = velocity_left2right_[0];
          double element_vy1 = velocity_left2right_[1];
          double element_vz1 = velocity_left2right_[2];
          tot_vx1_vol += element_vx1 * element_vol;
          tot_vy1_vol += element_vy1 * element_vol;
          tot_vz1_vol += element_vz1 * element_vol;

          (*it_e)->Read(velocity_bottom2top_key_,velocity_bottom2top_);
          double element_vx2 = velocity_bottom2top_[0];
          double element_vy2 = velocity_bottom2top_[1];
          double element_vz2 = velocity_bottom2top_[2];
          tot_vx2_vol += element_vx2 * element_vol;
          tot_vy2_vol += element_vy2 * element_vol;
          tot_vz2_vol += element_vz2 * element_vol;

          if ( dim == 3) {
              (*it_e)->Read(velocity_back2front_key_,velocity_back2front_);
              double element_vx3 = velocity_back2front_[0];
              double element_vy3 = velocity_back2front_[1];
              double element_vz3 = velocity_back2front_[2];
              tot_vx3_vol += element_vx3 * element_vol;
              tot_vy3_vol += element_vy3 * element_vol;
              tot_vz3_vol += element_vz3 * element_vol;
          }

          (*it_e)->Read(gradp_left2right_key_,gradp_left2right_);
          double element_dpdx1 = gradp_left2right_[0];
          double element_dpdy1 = gradp_left2right_[1];
          double element_dpdz1 = gradp_left2right_[2];
          tot_dpdx1_vol += element_dpdx1 * element_vol;
          tot_dpdy1_vol += element_dpdy1 * element_vol;
          tot_dpdz1_vol += element_dpdz1 * element_vol;

          (*it_e)->Read(gradp_bottom2top_key_,gradp_bottom2top_);
          double element_dpdx2 = gradp_bottom2top_[0];
          double element_dpdy2 = gradp_bottom2top_[1];
          double element_dpdz2 = gradp_bottom2top_[2];
          tot_dpdx2_vol += element_dpdx2 * element_vol;
          tot_dpdy2_vol += element_dpdy2 * element_vol;
          tot_dpdz2_vol += element_dpdz2 * element_vol;

          if ( dim == 3) {
              (*it_e)->Read(gradp_back2front_key_,gradp_back2front_);
              double element_dpdx3 = gradp_back2front_[0];
              double element_dpdy3 = gradp_back2front_[1];
              double element_dpdz3 = gradp_back2front_[2];
              tot_dpdx3_vol += element_dpdx3 * element_vol;
              tot_dpdy3_vol += element_dpdy3 * element_vol;
              tot_dpdz3_vol += element_dpdz3 * element_vol;
          }
    }

    vx1_avg_  = tot_vx1_vol / tot_vol;
    vy1_avg_  = tot_vy1_vol / tot_vol;
    vz1_avg_  = tot_vz1_vol / tot_vol;
    vx2_avg_  = tot_vx2_vol / tot_vol;
    vy2_avg_  = tot_vy2_vol / tot_vol;
    vz2_avg_  = tot_vz2_vol / tot_vol;

    dpx1_avg_ = tot_dpdx1_vol / tot_vol;
    dpy1_avg_ = tot_dpdy1_vol / tot_vol;
    dpz1_avg_ = tot_dpdz1_vol / tot_vol;
    dpx2_avg_ = tot_dpdx2_vol / tot_vol;
    dpy2_avg_ = tot_dpdy2_vol / tot_vol;
    dpz2_avg_ = tot_dpdz2_vol / tot_vol;

    if ( dim == 3) {
        vx3_avg_  = tot_vx3_vol   / tot_vol;
        vy3_avg_  = tot_vy3_vol   / tot_vol;
        vz3_avg_  = tot_vz3_vol   / tot_vol;
        dpx3_avg_ = tot_dpdx3_vol / tot_vol;
        dpy3_avg_ = tot_dpdy3_vol / tot_vol;
        dpz3_avg_ = tot_dpdz3_vol / tot_vol;
    }
}




/// sets the private variable 'viscosity_avg_'
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ViscosityVolumeAverage( const Model<dim>& model, const string regionName )
{
    double tot_vol     = 0.;
    double tot_vis_vol = 0.;
    const Region<dim>& r_ref(model.Region(regionName.c_str()));
    for ( typename vector<Element<dim>*>::const_iterator it_e = r_ref.CellsBegin(); it_e != r_ref.CellsEnd(); ++it_e ) {
        double element_vol = (*it_e)->Volume() * (*it_e)->Read(thic_key_);
        tot_vol += element_vol;
        double element_vis = (*it_e)->Read(viscosity_key_);
        tot_vis_vol += element_vis * element_vol;
    }
  
    viscosity_avg_ = tot_vis_vol / tot_vol;
}




/// computes the equivalent permeability tensor and stores it in the class
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::ComputeTensor()
{
    // symmetric effective permeability tensor
    size_t m, n;
    Matrix A;
    vector<double> b;

    if ( dim == 2) {
        m = 5; n = 4;
        A.Resize(m,n);
        A.Zero();

        A(0,0) = dpx1_avg_ / viscosity_avg_;
        A(0,1) = dpy1_avg_ / viscosity_avg_;

        A(1,2) = dpx1_avg_ / viscosity_avg_;
        A(1,3) = dpy1_avg_ / viscosity_avg_;

        A(2,0) = dpx2_avg_ / viscosity_avg_;
        A(2,1) = dpy2_avg_ / viscosity_avg_;

        A(3,2) = dpx2_avg_ / viscosity_avg_;
        A(3,3) = dpy2_avg_ / viscosity_avg_;

        A(4,0) = 0.;
        A(4,1) = 1.;
        A(4,2) = -1.;
        A(4,3) = 0.;

        b.resize(m);
        b[0] = - vx1_avg_;
        b[1] = - vy1_avg_;
        b[2] = - vx2_avg_;
        b[3] = - vy2_avg_;
        b[4] = 0.;
    }

    if ( dim == 3) {
        m = 12; n = 9;
        A.Resize(m,n);
        A.Zero();

        A(0,0) = dpx1_avg_ / viscosity_avg_;
        A(0,1) = dpy1_avg_ / viscosity_avg_;
        A(0,2) = dpz1_avg_ / viscosity_avg_;

        A(1,3) = dpx1_avg_ / viscosity_avg_;
        A(1,4) = dpy1_avg_ / viscosity_avg_;
        A(1,5) = dpz1_avg_ / viscosity_avg_;

        A(2,6) = dpx1_avg_ / viscosity_avg_;
        A(2,7) = dpy1_avg_ / viscosity_avg_;
        A(2,8) = dpz1_avg_ / viscosity_avg_;

        A(3,0) = dpx2_avg_ / viscosity_avg_;
        A(3,1) = dpy2_avg_ / viscosity_avg_;
        A(3,2) = dpz2_avg_ / viscosity_avg_;

        A(4,3) = dpx2_avg_ / viscosity_avg_;
        A(4,4) = dpy2_avg_ / viscosity_avg_;
        A(4,5) = dpz2_avg_ / viscosity_avg_;

        A(5,6) = dpx2_avg_ / viscosity_avg_;
        A(5,7) = dpy2_avg_ / viscosity_avg_;
        A(5,8) = dpz2_avg_ / viscosity_avg_;

        A(6,0) = dpx3_avg_ / viscosity_avg_;
        A(6,1) = dpy3_avg_ / viscosity_avg_;
        A(6,2) = dpz3_avg_ / viscosity_avg_;

        A(7,3) = dpx3_avg_ / viscosity_avg_;
        A(7,4) = dpy3_avg_ / viscosity_avg_;
        A(7,5) = dpz3_avg_ / viscosity_avg_;

        A(8,6) = dpx3_avg_ / viscosity_avg_;
        A(8,7) = dpy3_avg_ / viscosity_avg_;
        A(8,8) = dpz3_avg_ / viscosity_avg_;

        A(9,1) = 1.;
        A(9,3) = -1.;

        A(10,2) = 1.;
        A(10,6) = -1.;

        A(11,5) = 1.;
        A(11,7) = -1.;

        b.resize(m);
        b[0] = - vx1_avg_;
        b[1] = - vy1_avg_;
        b[2] = - vz1_avg_;
        b[3] = - vx2_avg_;
        b[4] = - vy2_avg_;
        b[5] = - vz2_avg_;
        b[6] = - vx3_avg_;
        b[7] = - vy3_avg_;
        b[8] = - vz3_avg_;
        b[9] =  0.;
        b[10] = 0.;
        b[11] = 0.;
    }

    vector<double> x(n);
    SolveOverdeterminedSystem(A, b, x);

    if ( dim == 2) {
        keq_(0,0) = x[0];
        keq_(0,1) = (x[1] + x[2]) / 2.;
        keq_(1,0) = (x[1] + x[2]) / 2.;
        keq_(1,1) = x[3];
    }
    if ( dim == 3) {
        keq_(0,0) = x[0];
        keq_(1,1) = x[4];
        keq_(2,2) = x[8];

        keq_(0,1) = (x[1] + x[3]) / 2.;
        keq_(1,0) = (x[1] + x[3]) / 2.;

        keq_(0,2) = (x[2] + x[6]) / 2.;
        keq_(2,0) = (x[2] + x[6]) / 2.;

        keq_(1,2) = (x[5] + x[7]) / 2.;
        keq_(2,1) = (x[5] + x[7]) / 2.;
    }
} // end ComputeTensor







template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::SolveOverdeterminedSystem( const Matrix& A, const vector<double>& b, vector<double>& x)
{
    size_t m = A.Rows();
    size_t n = A.Cols();    
    Matrix W(m,m);
    W.Zero();

    if ( dim == 2) {
        for ( size_t i = 0; i < 4; ++i)
            W(i,i) = 1.;
        W(4,4) = 1.e3;
    }
    if ( dim == 3) {
        for ( size_t i = 0; i < 9; ++i)
            W(i,i) = 1.;
        for ( size_t i = 9; i < 12; ++i)
            W(i,i) = 1.e3;
    }

    Matrix AT(n,m);
    A.Transposed(AT);

    Matrix ATW(n,m);
    ATW = AT * W;

    vector<double> ATWb(n);
    ATWb = ATW * b;

    Matrix ATWA(n,n);
    ATWA = ATW * A;

    SparseMatrix K(n);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            K.Assign(i,j,ATWA(i,j));
        }
    }

    vector<double> f(n);
    f = ATW * b;

    LUdcmp_Solver mySystem;
    mySystem.Solve(K, f, x, n);
}




/**
    computes eigenvalues, eigenvectors and orientation of permeability tensor
*/
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::TensorPrincipals()
{
    Matrix A(dim,dim);
    for ( uint32_t i{0}; i<dim; ++i )
      for ( uint32_t j{0}; j<dim; ++j )
        A(i,j) = keq_(i,j);
    // w: eigenvalues in ascending order
    vector<double> w(dim, 0.);
    // Q: corresponding eigenvectors (in columns)
    Matrix Q(dim,dim);
    Q.Zero();

    StatisticalDistributionGenerator::EigenDecomposition( A, Q, w );

    if ( dim == 2) {
        // orientation of the major axes of the ellipse
        double theta_rad = acos(Q(1,0));
        theta_ = theta_rad * 180 / PI;
        // principal components of the tensor (diagonal tensor) (Bear 1972)
        kmax_ = (keq_(0,0) + keq_(1,1)) / 2 + (keq_(0,0) - keq_(1,1)) / 2 * cos(2 * theta_rad) + keq_(0,1) * sin(2 * theta_rad);
        kmin_ = (keq_(0,0) + keq_(1,1)) / 2 - (keq_(0,0) - keq_(1,1)) / 2 * cos(2 * theta_rad) - keq_(0,1) * sin(2 * theta_rad);
        kmin_ = w[0];
        kmax_ = w[1];
        if (abs(theta_) < 0.1) theta_ = 0.;
    }
    if ( dim == 3) {
        // trend & plunge
        plunge_min_ = asin(Q(2,0)) * 180 / PI;
        if (Q(0,0) == 0.) {
            trend_min_ = atan(Q(1,0) / (Q(0,0) + 1.0e-6)) * 180 / PI ;
        }
        else {
            trend_min_ = atan(Q(1,0) / Q(0,0)) * 180 / PI ;
        }
        if (abs(plunge_min_) < 0.1) plunge_min_ = 0.;
        if (abs(trend_min_)  < 0.1) trend_min_  = 0.;

        plunge_max_ = asin(Q(2,2)) * 180 / PI;
        if (Q(0,2) == 0.) {
            trend_max_ = atan(Q(1,2) / (Q(0,2) + 1.0e-6)) * 180 / PI ;
        }
        else {
            trend_max_ = atan(Q(1,2) / Q(0,2)) * 180 / PI ;
        }
        if (abs(plunge_max_) < 0.1) plunge_max_ = 0.;
        if (abs(trend_max_)  < 0.1) trend_max_  = 0.;
        // principal components of the tensor (diagonal tensor)
        kmin_ = w[0];
        kmax_ = w[2];
    }
    for (uint32_t i = 0; i < dim; ++i) {
        for (uint32_t j = 0; j < dim; ++j) {
            if (fabs(keq_(i,j)) < eps_) keq_(i,j) = 0.;
        }
    }
}




template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::PermeabilityTolerance( const double val )
{
    eps_ = val;
}




template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::SpecifyOutputFileNames( string tensorFileName, string statisticsFileName )
{    
    tensor_file_name_ = tensorFileName;
    statistics_file_name_ = statisticsFileName;
}



/// output tensor
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::Out( TensorVariable<dim>& k_equiv_tensor ) const
 {
    k_equiv_tensor = keq_;
 }
 



/// output tensor components and its orientation to the specified file
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::Out() const
{
    ofstream file( tensor_file_name_, ios_base::app );
    if ( dim == 2 ) {
        file << setw(3)  << setprecision(3) << sample_size_;
        file << setw(10) << setprecision(3) << sample_elements_;
        file << setw(15) << setprecision(3) << keq_(0,0);
        file << setw(15) << setprecision(3) << keq_(0,1);
        file << setw(15) << setprecision(3) << keq_(1,0);
        file << setw(15) << setprecision(3) << keq_(1,1);
        file << setw(10) << fixed << setprecision(1) << theta_;
		file << setw(20) << scientific << setprecision(2) << kmin_;
		file << setw(15) << scientific << setprecision(2) << kmax_ << endl;
	}
    if ( dim == 3 ) {
        // 1st line
        file << setw(28) << setprecision(3) << keq_(0,0);
        file << setw(20) << setprecision(3) << keq_(0,1);
        file << setw(20) << setprecision(3) << keq_(0,2) << endl;
        // 2nd line
        file << setw(3)  << setprecision(3) << sample_size_;
        file << setw(10) << setprecision(3) << sample_elements_;
        file << setw(15) << setprecision(3) << keq_(1,0);
        file << setw(20) << setprecision(3) << keq_(1,1);
        file << setw(20) << setprecision(3) << keq_(1,2);
        file << setw(22) << setprecision(3) << kmin_;
        file << setw(22) << setprecision(3) << kmax_ << endl;
        // 3rd line
        file << setw(28) << setprecision(3) << keq_(2,0);
        file << setw(20) << setprecision(3) << keq_(2,1);
        file << setw(20) << setprecision(3) << keq_(2,2);
        file << setw(17) << fixed << setprecision(1) << trend_min_;
        file << setw(8)  << fixed << setprecision(1) << plunge_min_;
        file << setw(12) << fixed << setprecision(1) << trend_max_;
        file << setw(9)  << fixed << setprecision(1) << plunge_max_ << endl;
        file << "----------------------------------------------------------------------------------------------------------------------" << endl;
    }
}





/// output tensor statistics by random sampling to the specified file
template<uint32_t dim>
void EquivalentPermeabilityTensor<dim>::OutputStatistics() const
{
    ofstream file( statistics_file_name_ );
    size_t counter;
    for ( size_t j = 0; j < model_sample_size_.size(); ++j) {
        counter = j * max_no_sample_;
        file << "ld = " << fixed << setprecision(3) << ld_[j] << "    sample size = " << model_sample_size_[j] << endl;
        for (size_t i = 0; i < no_bins_; ++i) {
            file << setw(20) << fixed << setprecision(3) << model_bins_  [j][i];
            file << setw(20) << fixed << setprecision(3) << model_freqs_ [j][i] << endl;
        }
        file << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    }
   file.close();
}

template class EquivalentPermeabilityTensor<2U>;
template class EquivalentPermeabilityTensor<3U>;

} // end csmp

