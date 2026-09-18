// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "InputDataManager.h"
#include "CSMP_highLevelUtilities.h"

// Model
#include "Region.h"
#include "Model.h"

// Interface
#include "VTK_Interface.h"

#include "EquivalentPermeabilityTensor.h"
#include "fracturePropertyModelling.h"
#include "StatisticalDistributionGenerator.h"

using namespace std;

namespace csmp {

/**
    Calculates the equivalent full tensor permeability of the model
    using the velocity-volume averaging method described in Durlofsky (2005).

    @author Siroos Azizmohammadi (2013)
*/
template<uint32_t dim>
void fullTensorPermeability_Durlofsky( const char* model_name )
{
// ============================================
// setup csmp binary model
// ============================================
    Model<dim>  model(string{model_name});
    string      modelName(model_name);

// ============================================
// merge all fracture regions to "fractures"
// and matrix region to "matrix"
// ============================================
    string regionName;
    set<string> fracureRegionsList, matrixRegionsList;
    for ( auto it_r = model.UniqueRegionsBegin(); it_r != model.UniqueRegionsEnd(); ++it_r)
      {
         regionName = (*it_r).first;
         if ( regionName.compare(0,3,"FRA") == 0 ) fracureRegionsList.insert(regionName);
         if ( regionName.compare(0,3,"SET") == 0 ) fracureRegionsList.insert(regionName);
         if ( regionName.compare(0,3,"MAT") == 0 ) matrixRegionsList.insert(regionName);
      }
    if ( !model.ContainsRegion("fractures") && !fracureRegionsList.empty() ) model.MergeRegions(fracureRegionsList,"fractures");
    if ( !model.ContainsRegion("matrix") && !matrixRegionsList.empty() ) model.MergeRegions(matrixRegionsList,"matrix");

    VTK_Interface<dim> vtk_output;
    vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0);

// ============================================
// equivalent permeability tensor calculations
// ============================================
    double CPU_time;
    EquivalentPermeabilityTensor<dim> keq_tensor( model );
    string fileName = string(model.Name()) + "-k_eq_tensor.txt";
    keq_tensor.SpecifyOutputFileNames( fileName.c_str() );
    ofstream file( fileName.c_str() );
    if ( dim == 2 ) {
        file << "==========================================================================================================================" << endl;
        file << setw(62)  << "Equivalent Pemeability Tensor" << setw(25) << "Orientation" << setw(32) << "Tensor Principal Components" << endl;
        file << setw(4) << "Size" << setw(11) << "Elements" << setw(10) << "kxx" << setw(15) << "kxy" << setw(15) << "kyx" << setw(15) << "kyy" << setw(14) << "(Deg)" << setw(15) << "min" << setw(16) << "max" << endl;
        file << "--------------------------------------------------------------------------------------------------------------------------" << endl;
    }
    if ( dim == 3 ) {
        file << "======================================================================================================================" << endl;
        file << setw(57) << "Equivalent Pemeability Tensor" << setw(50) << "Principal Components" << endl;
        file << setw(25) << "kxx" << setw(21) << "kxy" << setw(19) << "kxz" << endl;
        file << setw(4) << "Size" << setw(10) << "Elements" << setw(11) << "kyx" << setw(21) << "kyy" << setw(19) << "kyz" << setw(24) << "kmin" << setw(20) << "kmax" << endl;
        file << setw(25) << "kzx" << setw(21) << "kzy" << setw(19) << "kzz" << setw(20) << "Trend" << setw(10) << "Plunge" << setw(10) << "Trend" << setw(10) << "Plunge" << endl;
        file << "----------------------------------------------------------------------------------------------------------------------" << endl;
    }
    bool loop_flag = true;
    size_t respond;
    while (loop_flag) {
        cout << "\n\n";
        cout << "---------------------------------------------" << endl;
        cout << " Equivalent Permeability Tensor Calculations " << endl;
        cout << "---------------------------------------------" << endl;
        cout << "(1) calculate for whole model" << endl;
        cout << "(2) calculate with random sampling" << endl;
        cout << "(3) quit" << endl;
        cout << "---------------------------------------------" << endl;
        cin >> respond;
        switch (respond) {
        case 1: {
            // assuming that the permeability of the model is single-valued
            auto time_start = clock();
            keq_tensor.ComputeEquivalentPermeabilityTensor( model );
            auto time_end = clock();
            CPU_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
            cout << endl << " cpu time = " << CPU_time << " sec" << endl;
            loop_flag = false;
        }
            break;
        case 2: {
            auto time_start = clock();
            // 4th argument, "sample size reducing factor", -1 is default value.
            // sample size = 1 / ( 1 + sample size reducing factor ) * model characteristic length
            keq_tensor.InitializeRandomSampling( 10,1,100,0.05 );
            keq_tensor.ComputeEquivalentPermeabilityTensorByRandomSampling( model );
            auto time_end = clock();
            CPU_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
            cout << endl << " cpu time = " << CPU_time << " sec";
            loop_flag = false;
        }
            break;
        case 3: { // SKM COMMENT: fix this illegal exit from switch
            loop_flag = false;
			      file.close();
            return;
        }
        default: {
            loop_flag = true;
        }
            break;
        }
    }
	file.close();
  cout << endl << "\nfullTensorPermeability_Durlofsky: Full Equivallent Permeability Tensor Analysis: DONE SUCCESFULLY\n"<< endl;
  cout << "see output file(s)! \n\n"<< endl;
    
} // end fullTensorPermeability_Durlofsky

template void fullTensorPermeability_Durlofsky<2U>( const char* );
template void fullTensorPermeability_Durlofsky<3U>( const char* );

} // end csmp
